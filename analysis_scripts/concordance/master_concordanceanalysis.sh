#!/bin/bash
set -uo pipefail

## =============================================================================
## Master Concordance Analysis: Illumina vs Aviti, PPR vs GATK
## FIXED VERSION - Resolves MT chromosome errors
## ENHANCEMENTS:
## - SNPs/INDELs/ALL separation
## - Combined Pos/GT heatmaps
## - S3 download for PPR VCFs (direct from DeepVariant output)
## - GATK PASS filtering
## - Non-ref filtering for PPR (like DeepVariant: exclude 0/0 and ./.)
## - Samples 21-24 added
## - MT exclusion BEFORE chromosome renaming (CRITICAL FIX!)
## - Suppressed stderr for normalization warnings
## =============================================================================

BASE_DIR="/mnt/genomics/pilot_PPR/concordance_IlluminaAvitiGATKPPR"
mkdir -p "${BASE_DIR}"

## Reference genomes
REF_GENOME="/mnt/genomics/illumina/ref_chr/Homo_sapiens.GRCh38.dna.primary_assembly.fasta"

## GATK VCF paths
ILLUMINA_GATK_VCF="/mnt/genomics/pilot_PPR/vcfs/pilotillumina.GATK.chr.vcf.gz"
AVITI_GATK_VCF="/mnt/genomics/juls/concordance/pilotaviti.GATK.vcf.gz"

## S3 paths for PPR VCFs
ILLUMINA_PPR_S3="s3://icdg/pilot1/illumina/processed2/PPR"
AVITI_PPR_S3="s3://icdg/pilot1/aviti/processed2/PPR"

## Sample mapping: Aviti ID -> Illumina ID (NOW WITH 21-24!)
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

echo "==================================================================="
echo "Starting FIXED concordance analysis"
echo "- Samples: 1-24 (24 total)"
echo "- GATK: PASS variants only"
echo "- PPR: Downloaded from S3, non-ref filtered"
echo "- Variant types: ALL, SNPs, INDELs"
echo "- MT exclusion: BEFORE chromosome renaming (FIXED!)"
echo "Base directory: ${BASE_DIR}"
echo "==================================================================="

## Create analysis scripts
echo "Creating analysis scripts..."

## 1. VCF preparation script for GATK files (WITH PASS FILTERING + MT FIX!)
cat > "${BASE_DIR}/prepare_gatk_vcfs.sh" << 'EOF'
#!/bin/bash
set -uo pipefail

AVITI_ID="$1"
ILLUMINA_ID="$2"
BASE_DIR="$3"

OUTDIR="${BASE_DIR}/sample_${AVITI_ID}"
mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

REF_GENOME="/mnt/genomics/illumina/ref_chr/Homo_sapiens.GRCh38.dna.primary_assembly.fasta"
ILLUMINA_GATK_VCF="/mnt/genomics/pilot_PPR/vcfs/pilotillumina.GATK.chr.vcf.gz"
AVITI_GATK_VCF="/mnt/genomics/juls/concordance/pilotaviti.GATK.vcf.gz"

echo "Processing GATK VCFs for ${AVITI_ID} (Aviti) / ${ILLUMINA_ID} (Illumina)"

## Process Illumina GATK (already has chr prefix) - PASS ONLY!
echo "  - Illumina GATK (PASS only)..."
bcftools view -s "${ILLUMINA_ID}" "${ILLUMINA_GATK_VCF}" | \
bcftools view -t "^chrM,chrMT" | \
bcftools view -f PASS | \
bcftools view -i 'GT!="0/0" && GT!="0|0" && GT!="./." && GT!=".|."' | \
bcftools norm -m -both -f "${REF_GENOME}" 2>/dev/null | \
awk 'BEGIN{OFS="\t"} /^#/{print; next} {gsub(/\|/, "/", $10); print}' | \
bcftools view -t chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY | \
bgzip -c > illumina.gatk.vcf.gz
tabix -f illumina.gatk.vcf.gz

## Process Aviti GATK - CRITICAL FIX: Exclude MT BEFORE renaming chromosomes!
echo "  - Aviti GATK (PASS only)..."
bcftools view -s "${AVITI_ID}" "${AVITI_GATK_VCF}" | \
bcftools view -t "^M,MT,chrM,chrMT" | \
bcftools view -f PASS | \
bcftools view -i 'GT!="0/0" && GT!="0|0" && GT!="./." && GT!=".|."' | \
bcftools annotate --rename-chrs <(printf "1\tchr1\n2\tchr2\n3\tchr3\n4\tchr4\n5\tchr5\n6\tchr6\n7\tchr7\n8\tchr8\n9\tchr9\n10\tchr10\n11\tchr11\n12\tchr12\n13\tchr13\n14\tchr14\n15\tchr15\n16\tchr16\n17\tchr17\n18\tchr18\n19\tchr19\n20\tchr20\n21\tchr21\n22\tchr22\nX\tchrX\nY\tchrY\n") 2>/dev/null | \
bcftools norm -m -both -f "${REF_GENOME}" 2>/dev/null | \
awk 'BEGIN{OFS="\t"} /^#/{print; next} {gsub(/\|/, "/", $10); print}' | \
bcftools view -t chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY 2>/dev/null | \
bgzip -c > aviti.gatk.vcf.gz
tabix -f aviti.gatk.vcf.gz

## Extract to TSV - ALL, SNPs, and INDELs
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n' illumina.gatk.vcf.gz > illumina.gatk.all.tsv
bcftools view -v snps illumina.gatk.vcf.gz | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n' > illumina.gatk.snps.tsv
bcftools view -v indels illumina.gatk.vcf.gz | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n' > illumina.gatk.indels.tsv

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n' aviti.gatk.vcf.gz > aviti.gatk.all.tsv
bcftools view -v snps aviti.gatk.vcf.gz | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n' > aviti.gatk.snps.tsv
bcftools view -v indels aviti.gatk.vcf.gz | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n' > aviti.gatk.indels.tsv

echo "  ✓ GATK VCFs prepared (Illumina - ALL: $(wc -l < illumina.gatk.all.tsv), SNPs: $(wc -l < illumina.gatk.snps.tsv), INDELs: $(wc -l < illumina.gatk.indels.tsv))"
echo "                        (Aviti - ALL: $(wc -l < aviti.gatk.all.tsv), SNPs: $(wc -l < aviti.gatk.snps.tsv), INDELs: $(wc -l < aviti.gatk.indels.tsv))"
EOF
chmod +x "${BASE_DIR}/prepare_gatk_vcfs.sh"

## 2. Download and process PPR VCFs from S3 (DIRECT FROM DEEPVARIANT!)
cat > "${BASE_DIR}/download_ppr_vcfs.sh" << 'EOF'
#!/bin/bash
set -uo pipefail

AVITI_ID="$1"
ILLUMINA_ID="$2"
BASE_DIR="$3"

OUTDIR="${BASE_DIR}/sample_${AVITI_ID}"
mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

ILLUMINA_PPR_S3="s3://icdg/pilot1/illumina/processed2/PPR"
AVITI_PPR_S3="s3://icdg/pilot1/aviti/processed2/PPR"
REF_GENOME="/mnt/genomics/illumina/ref_chr/Homo_sapiens.GRCh38.dna.primary_assembly.fasta"

echo "Downloading and processing PPR VCFs for ${AVITI_ID}/${ILLUMINA_ID}..."

## Download Illumina PPR from S3
echo "  - Downloading Illumina PPR from S3..."
ILLUMINA_PPR_S3_PATH="${ILLUMINA_PPR_S3}/${ILLUMINA_ID}/small_variants/${ILLUMINA_ID}.deepvariant.vcf.gz"
if aws s3 cp "${ILLUMINA_PPR_S3_PATH}" illumina.ppr.raw.vcf.gz 2>/dev/null; then
    echo "    ✓ Downloaded Illumina PPR"
    
    ## Process: filter non-ref, exclude MT, normalize, chr prefix
    bcftools view -i 'GT!="0/0" && GT!="./."' illumina.ppr.raw.vcf.gz | \
    bcftools view -t ^M,MT,chrM,chrMT | \
    bcftools norm -m -both -f "${REF_GENOME}" 2>/dev/null | \
    awk 'BEGIN{OFS="\t"} /^#/{print; next} $1 !~ /^chr/ {$1="chr"$1} {print}' | \
    bcftools view -t chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY | \
    bgzip -c > illumina.ppr.vcf.gz
    tabix -f illumina.ppr.vcf.gz
    
    ## Extract ALL, SNPs, and INDELs
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n' illumina.ppr.vcf.gz > illumina.ppr.all.tsv
    bcftools view -v snps illumina.ppr.vcf.gz | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n' > illumina.ppr.snps.tsv
    bcftools view -v indels illumina.ppr.vcf.gz | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n' > illumina.ppr.indels.tsv
    
    rm illumina.ppr.raw.vcf.gz
    echo "    ✓ Illumina PPR processed (ALL: $(wc -l < illumina.ppr.all.tsv), SNPs: $(wc -l < illumina.ppr.snps.tsv), INDELs: $(wc -l < illumina.ppr.indels.tsv))"
else
    echo "    ⚠ Illumina PPR not found in S3"
    return 1
fi

## Download Aviti PPR from S3
echo "  - Downloading Aviti PPR from S3..."
AVITI_PPR_S3_PATH="${AVITI_PPR_S3}/${AVITI_ID}/small_variants/${AVITI_ID}.deepvariant.vcf.gz"
if aws s3 cp "${AVITI_PPR_S3_PATH}" aviti.ppr.raw.vcf.gz 2>/dev/null; then
    echo "    ✓ Downloaded Aviti PPR"
    
    ## Process: filter non-ref, exclude MT, normalize, chr prefix
    bcftools view -i 'GT!="0/0" && GT!="./."' aviti.ppr.raw.vcf.gz | \
    bcftools view -t ^M,MT,chrM,chrMT | \
    bcftools norm -m -both -f "${REF_GENOME}" 2>/dev/null | \
    awk 'BEGIN{OFS="\t"} /^#/{print; next} $1 !~ /^chr/ {$1="chr"$1} {print}' | \
    bcftools view -t chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY | \
    bgzip -c > aviti.ppr.vcf.gz
    tabix -f aviti.ppr.vcf.gz
    
    ## Extract ALL, SNPs, and INDELs
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n' aviti.ppr.vcf.gz > aviti.ppr.all.tsv
    bcftools view -v snps aviti.ppr.vcf.gz | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n' > aviti.ppr.snps.tsv
    bcftools view -v indels aviti.ppr.vcf.gz | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n' > aviti.ppr.indels.tsv
    
    rm aviti.ppr.raw.vcf.gz
    echo "    ✓ Aviti PPR processed (ALL: $(wc -l < aviti.ppr.all.tsv), SNPs: $(wc -l < aviti.ppr.snps.tsv), INDELs: $(wc -l < aviti.ppr.indels.tsv))"
else
    echo "    ⚠ Aviti PPR not found in S3"
    return 1
fi
EOF
chmod +x "${BASE_DIR}/download_ppr_vcfs.sh"

## 3. R analysis script - ENHANCED with combined Pos/GT heatmaps
cat > "${BASE_DIR}/analyze_comprehensive_concordance.R" << 'EOF'
#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
    cat("Usage: Rscript analyze_comprehensive_concordance.R SAMPLE_DIR AVITI_ID\n")
    quit(status = 1)
}

sample_dir <- args[1]
aviti_id <- args[2]
setwd(sample_dir)

suppressPackageStartupMessages({
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
})

cat("\n=================================================================\n")
cat("Analyzing comprehensive concordance for", aviti_id, "\n")
cat("=================================================================\n")

## Function to calculate concordance between two methods
calc_concordance <- function(file1, file2, method1_name, method2_name, var_type) {
    if (!file.exists(file1) || !file.exists(file2)) {
        return(list(pos_conc = NA, gt_conc = NA))
    }
    
    df1 <- read.table(file1, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                      col.names = c("CHROM", "POS", "REF", "ALT", "GT"))
    df2 <- read.table(file2, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                      col.names = c("CHROM", "POS", "REF", "ALT", "GT"))
    
    ## Position key
    df1$pos_key <- paste(df1$CHROM, df1$POS, df1$REF, df1$ALT, sep = ":")
    df2$pos_key <- paste(df2$CHROM, df2$POS, df2$REF, df2$ALT, sep = ":")
    
    ## Genotype key
    df1$gt_key <- paste(df1$pos_key, df1$GT, sep = ":")
    df2$gt_key <- paste(df2$pos_key, df2$GT, sep = ":")
    
    pos1 <- unique(df1$pos_key)
    pos2 <- unique(df2$pos_key)
    gt1 <- unique(df1$gt_key)
    gt2 <- unique(df2$gt_key)
    
    total_pos <- length(union(pos1, pos2))
    pos_conc <- length(intersect(pos1, pos2))
    gt_conc <- length(intersect(gt1, gt2))
    
    pos_conc_pct <- if(total_pos > 0) round(100 * pos_conc / total_pos, 2) else NA
    gt_conc_pct <- if(total_pos > 0) round(100 * gt_conc / total_pos, 2) else NA
    
    cat(sprintf("  [%s] %s vs %s: Pos=%.2f%% GT=%.2f%%\n", 
                var_type, method1_name, method2_name, pos_conc_pct, gt_conc_pct))
    
    list(pos_conc = pos_conc_pct, gt_conc = gt_conc_pct)
}

## Process each variant type: ALL, SNPs, INDELs
variant_types <- c("all", "snps", "indels")
variant_labels <- c("ALL", "SNPs", "INDELs")

for (vt_idx in 1:length(variant_types)) {
    vt <- variant_types[vt_idx]
    vt_label <- variant_labels[vt_idx]
    
    cat("\n", strrep("=", 60), "\n", sep="")
    cat("Processing", vt_label, "variants\n")
    cat(strrep("=", 60), "\n")
    
    methods <- c("illumina.ppr", "illumina.gatk", "aviti.ppr", "aviti.gatk")
    method_labels <- c("Illumina\nPPR", "Illumina\nGATK", "Aviti\nPPR", "Aviti\nGATK")
    
    n <- length(methods)
    pos_matrix <- matrix(NA, nrow = n, ncol = n)
    gt_matrix <- matrix(NA, nrow = n, ncol = n)
    rownames(pos_matrix) <- colnames(pos_matrix) <- method_labels
    rownames(gt_matrix) <- colnames(gt_matrix) <- method_labels
    
    ## Fill diagonal with 100%
    diag(pos_matrix) <- 100
    diag(gt_matrix) <- 100
    
    ## Calculate all pairwise concordances
    for (i in 1:(n-1)) {
        for (j in (i+1):n) {
            file1 <- paste0(methods[i], ".", vt, ".tsv")
            file2 <- paste0(methods[j], ".", vt, ".tsv")
            
            result <- calc_concordance(file1, file2, method_labels[i], method_labels[j], vt_label)
            
            pos_matrix[i, j] <- result$pos_conc
            pos_matrix[j, i] <- result$pos_conc
            gt_matrix[i, j] <- result$gt_conc
            gt_matrix[j, i] <- result$gt_conc
        }
    }
    
    ## Create COMBINED Pos/GT matrix
    combined_matrix <- matrix(NA, nrow = n, ncol = n)
    rownames(combined_matrix) <- colnames(combined_matrix) <- method_labels
    
    for (i in 1:n) {
        for (j in 1:n) {
            if (i == j) {
                combined_matrix[i, j] <- "100 / 100"
            } else {
                combined_matrix[i, j] <- sprintf("%.1f / %.1f", pos_matrix[i, j], gt_matrix[i, j])
            }
        }
    }
    
    ## Save matrices
    write.table(pos_matrix, paste0(aviti_id, "_", vt, "_position_concordance.tsv"), 
                sep = "\t", quote = FALSE)
    write.table(gt_matrix, paste0(aviti_id, "_", vt, "_genotype_concordance.tsv"), 
                sep = "\t", quote = FALSE)
    write.table(combined_matrix, paste0(aviti_id, "_", vt, "_combined_concordance.tsv"), 
                sep = "\t", quote = FALSE)
    
    ## Create combined heatmap with Pos/GT in same cell
    create_combined_heatmap <- function(pos_mat, gt_mat, title, filename) {
        ## Convert to dataframes
        pos_df <- melt(pos_mat)
        gt_df <- melt(gt_mat)
        colnames(pos_df) <- c("Method1", "Method2", "Pos_Concordance")
        colnames(gt_df) <- c("Method1", "Method2", "GT_Concordance")
        
        ## Merge
        df <- merge(pos_df, gt_df, by = c("Method1", "Method2"))
        
        ## Create combined label
        df$label <- ifelse(df$Method1 == df$Method2, 
                           "100 / 100",
                           sprintf("%.1f / %.1f", df$Pos_Concordance, df$GT_Concordance))
        
        ## Use position concordance for color scale
        df$color_val <- df$Pos_Concordance
        
        p <- ggplot(df, aes(x = Method1, y = Method2, fill = color_val)) +
            geom_tile(color = "white", size = 1.5) +
            geom_text(aes(label = label), 
                      color = "black", fontface = "bold", size = 6) +
            scale_fill_gradient2(low = "#d73027", mid = "#fee08b", high = "#1a9850",
                                midpoint = 90, limits = c(70, 100),
                                na.value = "grey90") +
            labs(title = title,
                 subtitle = paste("Sample:", aviti_id, "| Format: Position% / Genotype%"),
                 x = "", y = "", fill = "Position\nConcordance (%)") +
            theme_minimal(base_size = 16) +
            theme(
                plot.title = element_text(hjust = 0.5, face = "bold", size = 22, color = "black"),
                plot.subtitle = element_text(hjust = 0.5, size = 14, color = "black"),
                axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, 
                                          face = "bold", size = 14, color = "black"),
                axis.text.y = element_text(face = "bold", size = 14, color = "black"),
                legend.title = element_text(face = "bold", size = 12, color = "black"),
                legend.text = element_text(size = 11, color = "black"),
                panel.grid = element_blank(),
                plot.background = element_rect(fill = "white", color = NA),
                panel.background = element_rect(fill = "white", color = NA)
            ) +
            coord_fixed()
        
        ggsave(filename, plot = p, width = 11, height = 10, dpi = 300, bg = "white")
        cat("  Saved:", filename, "\n")
    }
    
    ## Generate combined heatmap
    create_combined_heatmap(pos_matrix, gt_matrix, 
                           paste("Concordance Analysis -", vt_label),
                           paste0(aviti_id, "_", vt, "_concordance_heatmap.png"))
}

cat("\n✓ Analysis complete for", aviti_id, "\n")
EOF
chmod +x "${BASE_DIR}/analyze_comprehensive_concordance.R"

## 4. Combined plotting script for averages across all samples
cat > "${BASE_DIR}/create_combined_plots.R" << 'EOF'
#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
    cat("Usage: Rscript create_combined_plots.R BASE_DIR\n")
    quit(status = 1)
}

base_dir <- args[1]
setwd(base_dir)

suppressPackageStartupMessages({
    library(ggplot2)
    library(reshape2)
})

cat("\n=================================================================\n")
cat("Creating combined concordance summary plots\n")
cat("=================================================================\n")

## Collect all concordance matrices
sample_dirs <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)
sample_dirs <- sample_dirs[grepl("sample_ICDG", sample_dirs)]

variant_types <- c("all", "snps", "indels")
variant_labels <- c("ALL", "SNPs", "INDELs")

for (vt_idx in 1:length(variant_types)) {
    vt <- variant_types[vt_idx]
    vt_label <- variant_labels[vt_idx]
    
    cat("\nProcessing", vt_label, "variants...\n")
    
    all_pos_data <- list()
    all_gt_data <- list()
    
    for (dir in sample_dirs) {
        sample_id <- basename(dir)
        sample_id <- sub("sample_", "", sample_id)
        
        pos_file <- file.path(dir, paste0(sample_id, "_", vt, "_position_concordance.tsv"))
        gt_file <- file.path(dir, paste0(sample_id, "_", vt, "_genotype_concordance.tsv"))
        
        if (file.exists(pos_file) && file.exists(gt_file)) {
            pos_mat <- read.table(pos_file, header = TRUE, sep = "\t", row.names = 1)
            gt_mat <- read.table(gt_file, header = TRUE, sep = "\t", row.names = 1)
            
            all_pos_data[[sample_id]] <- pos_mat
            all_gt_data[[sample_id]] <- gt_mat
        }
    }
    
    if (length(all_pos_data) == 0) {
        cat("  WARNING: No data found for", vt_label, "\n")
        next
    }
    
    cat("  Found", length(all_pos_data), "samples\n")
    
    ## Calculate averages
    avg_pos <- Reduce("+", all_pos_data) / length(all_pos_data)
    avg_gt <- Reduce("+", all_gt_data) / length(all_gt_data)
    
    ## Create combined heatmap
    create_avg_combined_heatmap <- function(pos_mat, gt_mat, title, filename) {
        pos_df <- melt(as.matrix(pos_mat))
        gt_df <- melt(as.matrix(gt_mat))
        colnames(pos_df) <- c("Method1", "Method2", "Pos_Concordance")
        colnames(gt_df) <- c("Method1", "Method2", "GT_Concordance")
        
        df <- merge(pos_df, gt_df, by = c("Method1", "Method2"))
        df$label <- sprintf("%.1f / %.1f", df$Pos_Concordance, df$GT_Concordance)
        df$color_val <- df$Pos_Concordance
        
        p <- ggplot(df, aes(x = Method1, y = Method2, fill = color_val)) +
            geom_tile(color = "white", size = 2) +
            geom_text(aes(label = label), 
                      color = "black", fontface = "bold", size = 7) +
            scale_fill_gradient2(low = "#d73027", mid = "#fee08b", high = "#1a9850",
                                midpoint = 90, limits = c(70, 100)) +
            labs(title = title,
                 subtitle = sprintf("Average across %d samples | Format: Position%% / Genotype%%", 
                                   length(all_pos_data)),
                 x = "", y = "", fill = "Position\nConcordance (%)") +
            theme_minimal(base_size = 18) +
            theme(
                plot.title = element_text(hjust = 0.5, face = "bold", size = 24, color = "black"),
                plot.subtitle = element_text(hjust = 0.5, size = 15, color = "black"),
                axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, 
                                          face = "bold", size = 16, color = "black"),
                axis.text.y = element_text(face = "bold", size = 16, color = "black"),
                legend.title = element_text(face = "bold", size = 14, color = "black"),
                legend.text = element_text(size = 13, color = "black"),
                panel.grid = element_blank(),
                plot.background = element_rect(fill = "white", color = NA),
                panel.background = element_rect(fill = "white", color = NA)
            ) +
            coord_fixed()
        
        ggsave(filename, plot = p, width = 12, height = 11, dpi = 300, bg = "white")
        cat("  Saved:", filename, "\n")
    }
    
    create_avg_combined_heatmap(avg_pos, avg_gt,
                               paste("Average Concordance -", vt_label),
                               paste0("AVERAGE_", vt, "_concordance_heatmap.png"))
    
    ## Save average matrices
    write.table(avg_pos, paste0("AVERAGE_", vt, "_position_concordance.tsv"), 
                sep = "\t", quote = FALSE)
    write.table(avg_gt, paste0("AVERAGE_", vt, "_genotype_concordance.tsv"), 
                sep = "\t", quote = FALSE)
}

cat("\n✓ Combined plots created successfully\n")
EOF
chmod +x "${BASE_DIR}/create_combined_plots.R"

## =============================================================================
## Main processing loop
## =============================================================================

SUCCESS_COUNT=0
FAIL_COUNT=0
ERROR_LOG="${BASE_DIR}/processing_errors.log"
> "${ERROR_LOG}"

AVITI_SAMPLES=($(printf '%s\n' "${!SAMPLE_MAP[@]}" | sort -V))

for AVITI_ID in "${AVITI_SAMPLES[@]}"; do
    ILLUMINA_ID="${SAMPLE_MAP[$AVITI_ID]}"
    
    echo ""
    echo "==================================================================="
    echo "Processing: ${AVITI_ID} (Aviti) / ${ILLUMINA_ID} (Illumina)"
    echo "==================================================================="
    
    SAMPLE_DIR="${BASE_DIR}/sample_${AVITI_ID}"
    
    ## Step 1: Download PPR VCFs from S3
    echo "[1/3] Downloading PPR VCFs from S3..."
    if bash "${BASE_DIR}/download_ppr_vcfs.sh" "${AVITI_ID}" "${ILLUMINA_ID}" "${BASE_DIR}" 2>&1 | tee -a "${ERROR_LOG}"; then
        echo "  ✓ PPR VCFs downloaded and processed"
    else
        echo "  ⚠ PPR VCFs not available - skipping ${AVITI_ID}" | tee -a "${ERROR_LOG}"
        ((FAIL_COUNT++))
        continue
    fi
    
    ## Step 2: Prepare GATK VCFs
    echo "[2/3] Preparing GATK VCFs (PASS only)..."
    if bash "${BASE_DIR}/prepare_gatk_vcfs.sh" "${AVITI_ID}" "${ILLUMINA_ID}" "${BASE_DIR}" 2>&1 | tee -a "${ERROR_LOG}"; then
        echo "  ✓ GATK VCFs prepared"
    else
        echo "  ⚠ GATK preparation failed - skipping ${AVITI_ID}" | tee -a "${ERROR_LOG}"
        ((FAIL_COUNT++))
        continue
    fi
    
    ## Step 3: Analyze concordance
    echo "[3/3] Analyzing concordance..."
    if Rscript "${BASE_DIR}/analyze_comprehensive_concordance.R" "${SAMPLE_DIR}" "${AVITI_ID}" 2>&1 | tee -a "${ERROR_LOG}"; then
        echo "  ✓ Analysis complete"
        ((SUCCESS_COUNT++))
    else
        echo "  ⚠ Analysis failed - ${AVITI_ID}" | tee -a "${ERROR_LOG}"
        ((FAIL_COUNT++))
    fi
done

echo ""
echo "==================================================================="
echo "Individual sample processing complete"
echo "Successful: ${SUCCESS_COUNT} / Failed: ${FAIL_COUNT}"
echo "==================================================================="

## Create combined summary plots
if [ ${SUCCESS_COUNT} -gt 0 ]; then
    echo ""
    echo "Creating combined summary plots..."
    Rscript "${BASE_DIR}/create_combined_plots.R" "${BASE_DIR}"
fi

echo ""
echo "==================================================================="
echo "✓ FIXED Analysis complete!"
echo "==================================================================="
echo "Results: ${BASE_DIR}"
echo ""
echo "Configuration:"
echo "  - Samples: 1-24 (24 total)"
echo "  - GATK: PASS variants only"
echo "  - PPR: Downloaded from S3, non-ref filtered"
echo "  - Variant types: ALL, SNPs, INDELs"
echo "  - MT exclusion: BEFORE chromosome renaming (FIXED!)"
echo ""
echo "Individual samples:"
echo "  - sample_ICDG*/*_all_concordance_heatmap.png (ALL variants)"
echo "  - sample_ICDG*/*_snps_concordance_heatmap.png (SNPs only)"
echo "  - sample_ICDG*/*_indels_concordance_heatmap.png (INDELs only)"
echo ""
echo "Average across all samples:"
echo "  - AVERAGE_all_concordance_heatmap.png"
echo "  - AVERAGE_snps_concordance_heatmap.png"
echo "  - AVERAGE_indels_concordance_heatmap.png"
echo ""
echo "Format: Each cell shows 'Position% / Genotype%'"
if [ ${FAIL_COUNT} -gt 0 ]; then
    echo ""
    echo "⚠ Error log: ${ERROR_LOG}"
fi
echo "==================================================================="
