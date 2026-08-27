#!/bin/bash
set -uo pipefail

## =============================================================================
## Per-BED Region Concordance Analysis
## Analyzes concordance for specific genomic regions defined by BED files
## Excludes samples 12, 21, 22, 23
## Outputs: TSV files with concordance metrics (no per-sample plots)
## =============================================================================

## Directories
BASE_DIR="/mnt/genomics/pilot_PPR/concordance_perbed"
VCF_SOURCE_DIR="/mnt/genomics/pilot_PPR/concordance_IlluminaAvitiGATKPPR"
BED_BASE_DIR="/mnt/genomics/GIAB_stratifications/v3.6/with_chr"

mkdir -p "${BASE_DIR}"

## BED files to process (explicit list, excluding 25GC65)
declare -a BED_FILES=(
    "${BED_BASE_DIR}/OMIM.bed"
    "${BED_BASE_DIR}/ACMG.bed"
)

## Reference genome
REF_GENOME="/mnt/genomics/illumina/ref_chr/Homo_sapiens.GRCh38.dna.primary_assembly.fasta"

## Sample mapping (ALL samples, but we'll skip 12, 21, 22, 23)
declare -A SAMPLE_MAP=(
    ["ICDG1"]="1ICDG_S1"
    ["ICDG2"]="2ICDG_S4"
    ["ICDG3"]="3ICDG_S7"
    ["ICDG4"]="4ICDG_S10"
    ["ICDG5"]="5ICDG_S13"
    ["ICDG6"]="6ICDG_S16"
    ["ICDG7"]="7ICDG_S19"
    ["ICDG8"]="8ICDG_S22"
    ["ICDG9"]="9ICDG_S2"
    ["ICDG10"]="10ICDG_S5"
    ["ICDG11"]="11ICDG_S8"
    ["ICDG12"]="12ICDG_S11"
    ["ICDG13"]="13ICDG_S14"
    ["ICDG14"]="14ICDG_S17"
    ["ICDG15"]="15ICDG_S20"
    ["ICDG16"]="16ICDG_S23"
    ["ICDG17"]="17ICDG_S3"
    ["ICDG18"]="18ICDG_S6"
    ["ICDG19"]="19ICDG_S9"
    ["ICDG20"]="20ICDG_S12"
    ["ICDG21"]="21ICDG_S15"
    ["ICDG22"]="22ICDG_S18"
    ["ICDG23"]="23ICDG_S21"
    ["ICDG24"]="24ICDG_S24"
)

## Samples to EXCLUDE
EXCLUDE_SAMPLES=("ICDG12" "ICDG21" "ICDG22" "ICDG23")

echo "==================================================================="
echo "Per-BED Region Concordance Analysis"
echo "==================================================================="
echo "Base directory: ${BASE_DIR}"
echo "BED files to process: ${#BED_FILES[@]}"
for bed in "${BED_FILES[@]}"; do
    echo "  - $(basename ${bed})"
done
echo "Excluded samples: ${EXCLUDE_SAMPLES[@]}"
echo "==================================================================="

## Create R analysis script
cat > "${BASE_DIR}/calculate_concordance_metrics.R" << 'RSCRIPT_EOF'
#!/usr/bin/env Rscript

## Calculate concordance metrics between two TSV files
## Args: file1 file2 method1_name method2_name variant_type output_file sample_id bed_name

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 8) {
    cat("Usage: Rscript calculate_concordance_metrics.R file1 file2 method1 method2 variant_type output_file sample_id bed_name\n")
    quit(status = 1)
}

file1 <- args[1]
file2 <- args[2]
method1_name <- args[3]
method2_name <- args[4]
variant_type <- args[5]
output_file <- args[6]
sample_id <- args[7]
bed_name <- args[8]

## Read data
if (!file.exists(file1) || !file.exists(file2)) {
    cat("WARNING: Missing files\n")
    result <- data.frame(
        Sample = sample_id,
        BED_Region = bed_name,
        Variant_Type = variant_type,
        Comparison = paste0(method1_name, "_vs_", method2_name),
        Pos_Conc_Pct = NA,
        GT_Conc_Pct = NA,
        N_Method1_Unique = NA,
        N_Method2_Unique = NA,
        N_Common_Pos = NA,
        N_Common_GT = NA,
        N_Total_Union = NA
    )
    
    if (file.exists(output_file)) {
        write.table(result, output_file, sep = "\t", quote = FALSE, row.names = FALSE, append = TRUE, col.names = FALSE)
    } else {
        write.table(result, output_file, sep = "\t", quote = FALSE, row.names = FALSE, append = FALSE, col.names = TRUE)
    }
    quit(status = 0)
}

df1 <- read.table(file1, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                  col.names = c("CHROM", "POS", "REF", "ALT", "GT"))
df2 <- read.table(file2, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                  col.names = c("CHROM", "POS", "REF", "ALT", "GT"))

## Position keys
df1$pos_key <- paste(df1$CHROM, df1$POS, df1$REF, df1$ALT, sep = ":")
df2$pos_key <- paste(df2$CHROM, df2$POS, df2$REF, df2$ALT, sep = ":")

## Genotype keys
df1$gt_key <- paste(df1$pos_key, df1$GT, sep = ":")
df2$gt_key <- paste(df2$pos_key, df2$GT, sep = ":")

## Unique sets
pos1 <- unique(df1$pos_key)
pos2 <- unique(df2$pos_key)
gt1 <- unique(df1$gt_key)
gt2 <- unique(df2$gt_key)

## Calculate metrics
total_pos <- length(union(pos1, pos2))
common_pos <- length(intersect(pos1, pos2))
common_gt <- length(intersect(gt1, gt2))

unique_1 <- length(setdiff(pos1, pos2))
unique_2 <- length(setdiff(pos2, pos1))

pos_conc_pct <- if(total_pos > 0) round(100 * common_pos / total_pos, 2) else NA
gt_conc_pct <- if(total_pos > 0) round(100 * common_gt / total_pos, 2) else NA

## Create result
result <- data.frame(
    Sample = sample_id,
    BED_Region = bed_name,
    Variant_Type = variant_type,
    Comparison = paste0(method1_name, "_vs_", method2_name),
    Pos_Conc_Pct = pos_conc_pct,
    GT_Conc_Pct = gt_conc_pct,
    N_Method1_Unique = unique_1,
    N_Method2_Unique = unique_2,
    N_Common_Pos = common_pos,
    N_Common_GT = common_gt,
    N_Total_Union = total_pos
)

## Append to file (or create with header)
if (file.exists(output_file)) {
    write.table(result, output_file, sep = "\t", quote = FALSE, row.names = FALSE, append = TRUE, col.names = FALSE)
} else {
    write.table(result, output_file, sep = "\t", quote = FALSE, row.names = FALSE, append = FALSE, col.names = TRUE)
}

cat(sprintf("[%s][%s] %s vs %s: Pos=%.2f%% GT=%.2f%% (Total=%d)\n", 
            bed_name, variant_type, method1_name, method2_name, pos_conc_pct, gt_conc_pct, total_pos))
RSCRIPT_EOF
chmod +x "${BASE_DIR}/calculate_concordance_metrics.R"

## =============================================================================
## Main processing loop
## =============================================================================

## Get list of samples to process (excluding 12, 21, 22, 23)
SAMPLES_TO_PROCESS=()
for AVITI_ID in $(printf '%s\n' "${!SAMPLE_MAP[@]}" | sort -V); do
    SKIP=0
    for EXCLUDE in "${EXCLUDE_SAMPLES[@]}"; do
        if [ "$AVITI_ID" == "$EXCLUDE" ]; then
            SKIP=1
            break
        fi
    done
    if [ $SKIP -eq 0 ]; then
        SAMPLES_TO_PROCESS+=("$AVITI_ID")
    fi
done

echo ""
echo "Processing ${#SAMPLES_TO_PROCESS[@]} samples (excluded: ${EXCLUDE_SAMPLES[@]})"
echo ""

## Process each BED file
for BED_FILE in "${BED_FILES[@]}"; do
    BED_NAME=$(basename "${BED_FILE}" .bed)
    
    echo "==================================================================="
    echo "Processing BED region: ${BED_NAME}"
    echo "==================================================================="
    
    ## Create output directory for this BED
    BED_OUT_DIR="${BASE_DIR}/${BED_NAME}"
    mkdir -p "${BED_OUT_DIR}"
    
    ## Check if BED file exists
    if [ ! -f "${BED_FILE}" ]; then
        echo "⚠ BED file not found: ${BED_FILE} - SKIPPING"
        continue
    fi
    
    ## Process each sample
    for AVITI_ID in "${SAMPLES_TO_PROCESS[@]}"; do
        echo ""
        echo "-------------------------------------------------------------------"
        echo "Sample: ${AVITI_ID} | BED: ${BED_NAME}"
        echo "-------------------------------------------------------------------"
        
        SAMPLE_SOURCE_DIR="${VCF_SOURCE_DIR}/sample_${AVITI_ID}"
        SAMPLE_OUT_DIR="${BED_OUT_DIR}/sample_${AVITI_ID}"
        mkdir -p "${SAMPLE_OUT_DIR}"
        
        ## Check if source VCFs exist
        if [ ! -d "${SAMPLE_SOURCE_DIR}" ]; then
            echo "⚠ Source directory not found: ${SAMPLE_SOURCE_DIR} - SKIPPING"
            continue
        fi
        
        ## Define VCF files (from whole-genome analysis)
        ILLUMINA_PPR_VCF="${SAMPLE_SOURCE_DIR}/illumina.ppr.vcf.gz"
        ILLUMINA_GATK_VCF="${SAMPLE_SOURCE_DIR}/illumina.gatk.vcf.gz"
        AVITI_PPR_VCF="${SAMPLE_SOURCE_DIR}/aviti.ppr.vcf.gz"
        AVITI_GATK_VCF="${SAMPLE_SOURCE_DIR}/aviti.gatk.vcf.gz"
        
        ## Check if all VCFs exist
        ALL_EXIST=1
        for VCF in "${ILLUMINA_PPR_VCF}" "${ILLUMINA_GATK_VCF}" "${AVITI_PPR_VCF}" "${AVITI_GATK_VCF}"; do
            if [ ! -f "${VCF}" ]; then
                echo "⚠ Missing VCF: $(basename ${VCF})"
                ALL_EXIST=0
            fi
        done
        
        if [ $ALL_EXIST -eq 0 ]; then
            echo "⚠ Skipping ${AVITI_ID} - incomplete VCF set"
            continue
        fi
        
        ## Intersect VCFs with BED region
        echo "Intersecting VCFs with BED region..."
        
        bcftools view -R "${BED_FILE}" "${ILLUMINA_PPR_VCF}" -Oz -o "${SAMPLE_OUT_DIR}/illumina.ppr.bed.vcf.gz" 2>/dev/null
        tabix -f "${SAMPLE_OUT_DIR}/illumina.ppr.bed.vcf.gz" 2>/dev/null
        
        bcftools view -R "${BED_FILE}" "${ILLUMINA_GATK_VCF}" -Oz -o "${SAMPLE_OUT_DIR}/illumina.gatk.bed.vcf.gz" 2>/dev/null
        tabix -f "${SAMPLE_OUT_DIR}/illumina.gatk.bed.vcf.gz" 2>/dev/null
        
        bcftools view -R "${BED_FILE}" "${AVITI_PPR_VCF}" -Oz -o "${SAMPLE_OUT_DIR}/aviti.ppr.bed.vcf.gz" 2>/dev/null
        tabix -f "${SAMPLE_OUT_DIR}/aviti.ppr.bed.vcf.gz" 2>/dev/null
        
        bcftools view -R "${BED_FILE}" "${AVITI_GATK_VCF}" -Oz -o "${SAMPLE_OUT_DIR}/aviti.gatk.bed.vcf.gz" 2>/dev/null
        tabix -f "${SAMPLE_OUT_DIR}/aviti.gatk.bed.vcf.gz" 2>/dev/null
        
        ## Extract to TSV for ALL, SNPs, INDELs
        echo "Extracting variants (ALL, SNPs, INDELs)..."
        
        ## Illumina PPR
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n' "${SAMPLE_OUT_DIR}/illumina.ppr.bed.vcf.gz" > "${SAMPLE_OUT_DIR}/illumina.ppr.all.tsv"
        bcftools view -v snps "${SAMPLE_OUT_DIR}/illumina.ppr.bed.vcf.gz" | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n' > "${SAMPLE_OUT_DIR}/illumina.ppr.snps.tsv"
        bcftools view -v indels "${SAMPLE_OUT_DIR}/illumina.ppr.bed.vcf.gz" | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n' > "${SAMPLE_OUT_DIR}/illumina.ppr.indels.tsv"
        
        ## Illumina GATK
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n' "${SAMPLE_OUT_DIR}/illumina.gatk.bed.vcf.gz" > "${SAMPLE_OUT_DIR}/illumina.gatk.all.tsv"
        bcftools view -v snps "${SAMPLE_OUT_DIR}/illumina.gatk.bed.vcf.gz" | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n' > "${SAMPLE_OUT_DIR}/illumina.gatk.snps.tsv"
        bcftools view -v indels "${SAMPLE_OUT_DIR}/illumina.gatk.bed.vcf.gz" | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n' > "${SAMPLE_OUT_DIR}/illumina.gatk.indels.tsv"
        
        ## Aviti PPR
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n' "${SAMPLE_OUT_DIR}/aviti.ppr.bed.vcf.gz" > "${SAMPLE_OUT_DIR}/aviti.ppr.all.tsv"
        bcftools view -v snps "${SAMPLE_OUT_DIR}/aviti.ppr.bed.vcf.gz" | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n' > "${SAMPLE_OUT_DIR}/aviti.ppr.snps.tsv"
        bcftools view -v indels "${SAMPLE_OUT_DIR}/aviti.ppr.bed.vcf.gz" | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n' > "${SAMPLE_OUT_DIR}/aviti.ppr.indels.tsv"
        
        ## Aviti GATK
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n' "${SAMPLE_OUT_DIR}/aviti.gatk.bed.vcf.gz" > "${SAMPLE_OUT_DIR}/aviti.gatk.all.tsv"
        bcftools view -v snps "${SAMPLE_OUT_DIR}/aviti.gatk.bed.vcf.gz" | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n' > "${SAMPLE_OUT_DIR}/aviti.gatk.snps.tsv"
        bcftools view -v indels "${SAMPLE_OUT_DIR}/aviti.gatk.bed.vcf.gz" | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n' > "${SAMPLE_OUT_DIR}/aviti.gatk.indels.tsv"
        
        ## Output file for this sample's results
        SAMPLE_RESULTS="${BED_OUT_DIR}/sample_${AVITI_ID}_concordance_results.tsv"
        
        ## Calculate concordance for all comparisons and variant types
        echo "Calculating concordance metrics..."
        
        ## Variant types
        VARIANT_TYPES=("all" "snps" "indels")
        
        for VT in "${VARIANT_TYPES[@]}"; do
            ## Comparison 1: Illumina PPR vs Illumina GATK
            Rscript "${BASE_DIR}/calculate_concordance_metrics.R" \
                "${SAMPLE_OUT_DIR}/illumina.ppr.${VT}.tsv" \
                "${SAMPLE_OUT_DIR}/illumina.gatk.${VT}.tsv" \
                "IlluminaPPR" "IlluminaGATK" "${VT}" \
                "${SAMPLE_RESULTS}" "${AVITI_ID}" "${BED_NAME}"
            
            ## Comparison 2: Illumina PPR vs Aviti PPR
            Rscript "${BASE_DIR}/calculate_concordance_metrics.R" \
                "${SAMPLE_OUT_DIR}/illumina.ppr.${VT}.tsv" \
                "${SAMPLE_OUT_DIR}/aviti.ppr.${VT}.tsv" \
                "IlluminaPPR" "AvitiPPR" "${VT}" \
                "${SAMPLE_RESULTS}" "${AVITI_ID}" "${BED_NAME}"
            
            ## Comparison 3: Aviti PPR vs Aviti GATK
            Rscript "${BASE_DIR}/calculate_concordance_metrics.R" \
                "${SAMPLE_OUT_DIR}/aviti.ppr.${VT}.tsv" \
                "${SAMPLE_OUT_DIR}/aviti.gatk.${VT}.tsv" \
                "AvitiPPR" "AvitiGATK" "${VT}" \
                "${SAMPLE_RESULTS}" "${AVITI_ID}" "${BED_NAME}"
            
            ## Comparison 4: Illumina GATK vs Aviti GATK
            Rscript "${BASE_DIR}/calculate_concordance_metrics.R" \
                "${SAMPLE_OUT_DIR}/illumina.gatk.${VT}.tsv" \
                "${SAMPLE_OUT_DIR}/aviti.gatk.${VT}.tsv" \
                "IlluminaGATK" "AvitiGATK" "${VT}" \
                "${SAMPLE_RESULTS}" "${AVITI_ID}" "${BED_NAME}"
        done
        
        echo "✓ Sample ${AVITI_ID} complete - results: ${SAMPLE_RESULTS}"
        
        ## Clean up intermediate VCF files to save space
        rm -f "${SAMPLE_OUT_DIR}"/*.bed.vcf.gz* "${SAMPLE_OUT_DIR}"/*.tsv
    done
    
    echo ""
    echo "✓ BED region ${BED_NAME} complete!"
done

echo ""
echo "==================================================================="
echo "✓ Per-BED concordance analysis complete!"
echo "==================================================================="
echo "Results directory: ${BASE_DIR}"
echo ""
echo "Structure:"
for BED_FILE in "${BED_FILES[@]}"; do
    BED_NAME=$(basename "${BED_FILE}" .bed)
    echo "  ${BED_NAME}/"
    echo "    └── sample_*_concordance_results.tsv"
done
echo ""
echo "Next step: Run plot_concordance_lineplot.R to generate summary plots"
echo "==================================================================="
