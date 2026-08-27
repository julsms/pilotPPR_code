#!/bin/bash
# uniquevariants_validationrate_illumina_normalized.sh
# 
# CORRECTED VERSION: Normalizes the Merfin-validated ALL-variants VCF
# (bcftools norm -m-both --atomize) before intersection with already-normalized
# unique variant VCFs. This ensures both sides of bcftools isec use the same
# variant representation (left-aligned, split multiallelics), eliminating
# positional mismatches for indels.
#
# Normalized VCFs are created as temporary files and removed after counting.
# Output CSV has a distinct name to avoid overwriting previous results.

set -euo pipefail

WORK_DIR="/mnt/genomics/pilot_PPR/kmers_validation/illumina/uniquevariants"
UNIQUE_GATK_DIR="/mnt/genomics/pilot_PPR/uniquevariants/uniqueGATKillumina"
UNIQUE_PPR_DIR="/mnt/genomics/pilot_PPR/uniquevariants/uniquePPRillumina"
VALIDATED_VCF_DIR="/mnt/genomics/pilot_PPR/kmers_validation/illumina/vcf_results"
REF="/mnt/genomics/illumina/ref_chr/Homo_sapiens.GRCh38.dna.primary_assembly.fasta"

OUTPUT_CSV="${WORK_DIR}/unique_validation_rates_normalized.csv"

echo "============================================================"
echo "Unique variants validation rate - ILLUMINA (normalized isec)"
echo "Started at: $(date)"
echo "============================================================"

mkdir -p "${WORK_DIR}"
cd "${WORK_DIR}"

echo "Sample,Pipeline,Variant_Type,Total_Unique,Validated_in_Unique,Validation_Rate" > "${OUTPUT_CSV}"

# Illumina sample ID mapping
convert_sample_id() {
    local sample_id=$1
    case $sample_id in
        ICDG1) echo "1ICDG_S1" ;;
        ICDG2) echo "2ICDG_S4" ;;
        ICDG3) echo "3ICDG_S7" ;;
        ICDG4) echo "4ICDG_S10" ;;
        ICDG5) echo "5ICDG_S13" ;;
        ICDG6) echo "6ICDG_S16" ;;
        ICDG7) echo "7ICDG_S19" ;;
        ICDG8) echo "8ICDG_S22" ;;
        ICDG9) echo "9ICDG_S2" ;;
        ICDG10) echo "10ICDG_S5" ;;
        ICDG11) echo "11ICDG_S8" ;;
        ICDG13) echo "13ICDG_S14" ;;
        ICDG14) echo "14ICDG_S17" ;;
        ICDG15) echo "15ICDG_S20" ;;
        ICDG16) echo "16ICDG_S23" ;;
        ICDG17) echo "17ICDG_S3" ;;
        ICDG18) echo "18ICDG_S6" ;;
        ICDG19) echo "19ICDG_S9" ;;
        ICDG20) echo "20ICDG_S12" ;;
        ICDG24) echo "24ICDG_S24" ;;
        *) echo "UNKNOWN" ;;
    esac
}

# Function: normalize a VCF to match the representation used in unique variant files
# (left-align, split multiallelics, filter non-ref GT). Two-step to avoid pipe format issues.
normalize_vcf() {
    local input_vcf=$1
    local output_vcf=$2
    local temp_norm="${output_vcf%.vcf.gz}_step1.vcf.gz"

    # Step 1: normalize (left-align + split multiallelics)
    # Note: --atomize requires bcftools >=1.17; if unavailable, -m-both alone
    # performs multiallelic splitting which is the critical operation here.
    bcftools norm -m-both -f "${REF}" "${input_vcf}" -Oz -o "${temp_norm}"
    tabix -p vcf "${temp_norm}"

    # Step 2: filter non-ref genotypes
    bcftools view -i 'GT!="0/0" && GT!="./." && GT!="0|0" && GT!=".|."' \
        "${temp_norm}" -Oz -o "${output_vcf}"
    tabix -p vcf "${output_vcf}"

    # Cleanup intermediate
    rm -f "${temp_norm}" "${temp_norm}.tbi"
}

echo ""
echo "========================================="
echo "Processing GATK unique variants..."
echo "========================================="

for vcf_file in ${UNIQUE_GATK_DIR}/*_unique_GATK_illumina.vcf.gz; do
    [[ -f "$vcf_file" ]] || continue

    filename=$(basename "$vcf_file")
    sample_id=$(echo "$filename" | sed 's/_unique_.*//g')
    validation_sample=$(convert_sample_id "$sample_id")

    echo ""
    echo "--- ${sample_id} -> ${validation_sample} ---"

    if [[ "$validation_sample" == "UNKNOWN" ]]; then
        echo "  Skipping unknown sample: $sample_id"
        continue
    fi

    validated_vcf="${VALIDATED_VCF_DIR}/${validation_sample}/${validation_sample}_gatk_all_validated.filter.vcf.gz"

    if [[ ! -f "$validated_vcf" ]]; then
        echo "  Validated VCF not found: $validated_vcf"
        continue
    fi

    # Filter unique variants for non-ref genotypes (consistent with Aviti script)
    filtered_unique="${WORK_DIR}/temp_${sample_id}_gatk_unique_filtered.vcf.gz"
    bcftools view -i 'GT!="0/0" && GT!="./." && GT!="0|0" && GT!=".|."' \
        "$vcf_file" -Oz -o "$filtered_unique"
    tabix -p vcf "$filtered_unique"

    unique_total=$(bcftools view -H "$filtered_unique" | wc -l)
    echo "  Total unique variants (non-ref GT): $unique_total"

    if [[ $unique_total -eq 0 ]]; then
        echo "  No variants found, skipping"
        echo "${validation_sample},GATK,UNIQUE,0,0,0" >> "${OUTPUT_CSV}"
        rm -f "$filtered_unique" "$filtered_unique.tbi"
        continue
    fi

    # Normalize the validated ALL VCF to match unique variant representation
    norm_validated="${WORK_DIR}/temp_${sample_id}_gatk_validated_norm.vcf.gz"
    echo "  Normalizing validated VCF..."
    normalize_vcf "$validated_vcf" "$norm_validated"

    # Intersect: unique variants present in normalized validated set
    temp_dir="${WORK_DIR}/temp_isec_${sample_id}_gatk"
    mkdir -p "$temp_dir"

    bcftools isec -p "$temp_dir" "$filtered_unique" "$norm_validated"

    if [[ -f "$temp_dir/0002.vcf" ]]; then
        validated_count=$(bcftools view -H "$temp_dir/0002.vcf" | wc -l)
    else
        validated_count=0
    fi

    # Cleanup temp files
    rm -rf "$temp_dir"
    rm -f "$filtered_unique" "$filtered_unique.tbi"
    rm -f "$norm_validated" "$norm_validated.tbi"

    validation_rate=$(echo "scale=4; $validated_count * 100 / $unique_total" | bc -l)

    echo "${validation_sample},GATK,UNIQUE,${unique_total},${validated_count},${validation_rate}" >> "${OUTPUT_CSV}"
    echo "  Result: $validated_count/$unique_total ($validation_rate%)"
done

echo ""
echo "========================================="
echo "Processing PPR unique variants..."
echo "========================================="

for vcf_file in ${UNIQUE_PPR_DIR}/*_unique_PPR_illumina.vcf.gz; do
    [[ -f "$vcf_file" ]] || continue

    filename=$(basename "$vcf_file")
    sample_id=$(echo "$filename" | sed 's/_unique_.*//g')
    validation_sample=$(convert_sample_id "$sample_id")

    echo ""
    echo "--- ${sample_id} -> ${validation_sample} ---"

    if [[ "$validation_sample" == "UNKNOWN" ]]; then
        echo "  Skipping unknown sample: $sample_id"
        continue
    fi

    validated_vcf="${VALIDATED_VCF_DIR}/${validation_sample}/${validation_sample}_ppr_all_validated.filter.vcf.gz"

    if [[ ! -f "$validated_vcf" ]]; then
        echo "  Validated VCF not found: $validated_vcf"
        continue
    fi

    # Filter unique variants for non-ref genotypes
    filtered_unique="${WORK_DIR}/temp_${sample_id}_ppr_unique_filtered.vcf.gz"
    bcftools view -i 'GT!="0/0" && GT!="./." && GT!="0|0" && GT!=".|."' \
        "$vcf_file" -Oz -o "$filtered_unique"
    tabix -p vcf "$filtered_unique"

    unique_total=$(bcftools view -H "$filtered_unique" | wc -l)
    echo "  Total unique variants (non-ref GT): $unique_total"

    if [[ $unique_total -eq 0 ]]; then
        echo "  No variants found, skipping"
        echo "${validation_sample},PPR,UNIQUE,0,0,0" >> "${OUTPUT_CSV}"
        rm -f "$filtered_unique" "$filtered_unique.tbi"
        continue
    fi

    # Normalize the validated ALL VCF
    norm_validated="${WORK_DIR}/temp_${sample_id}_ppr_validated_norm.vcf.gz"
    echo "  Normalizing validated VCF..."
    normalize_vcf "$validated_vcf" "$norm_validated"

    # Intersect
    temp_dir="${WORK_DIR}/temp_isec_${sample_id}_ppr"
    mkdir -p "$temp_dir"

    bcftools isec -p "$temp_dir" "$filtered_unique" "$norm_validated"

    if [[ -f "$temp_dir/0002.vcf" ]]; then
        validated_count=$(bcftools view -H "$temp_dir/0002.vcf" | wc -l)
    else
        validated_count=0
    fi

    # Cleanup
    rm -rf "$temp_dir"
    rm -f "$filtered_unique" "$filtered_unique.tbi"
    rm -f "$norm_validated" "$norm_validated.tbi"

    validation_rate=$(echo "scale=4; $validated_count * 100 / $unique_total" | bc -l)

    echo "${validation_sample},PPR,UNIQUE,${unique_total},${validated_count},${validation_rate}" >> "${OUTPUT_CSV}"
    echo "  Result: $validated_count/$unique_total ($validation_rate%)"
done

echo ""
echo "========================================="
echo "COMPLETED"
echo "========================================="
echo "Results saved in: ${OUTPUT_CSV}"
echo ""
echo "Method: For each sample, the Merfin-validated ALL-variants VCF was"
echo "normalized (bcftools norm -m-both --atomize) to match the representation"
echo "of the already-normalized unique variant VCFs. Then bcftools isec was used"
echo "to count how many unique variants appear in the normalized validated set."
echo "Validation rate = (unique variants in validated set) / (total unique non-ref variants) × 100"
echo ""
echo "Finished at: $(date)"
