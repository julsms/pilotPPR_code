#!/bin/bash

set -euo pipefail

BASE_DIR="/mnt/genomics/pilot_PPR/uniquevariants"
OUTPUT_DIR="${BASE_DIR}/snv_indels_plots"

echo "Counting SNVs, MNPs, INDELs, and OTHER variants using bcftools"
echo "Output directory: ${OUTPUT_DIR}"
echo ""

mkdir -p "${OUTPUT_DIR}"

SAMPLE_COUNTS="${OUTPUT_DIR}/sample_snv_indel_counts.tsv"
MULTISAMPLE_COUNTS="${OUTPUT_DIR}/multisample_snv_indel_counts.tsv"

echo -e "Category\tSample\tSNPs\tMNPs\tINDELs\tOTHER\tTotal" > "${SAMPLE_COUNTS}"
echo -e "VCF_Type\tSNPs\tMNPs\tINDELs\tOTHER\tTotal" > "${MULTISAMPLE_COUNTS}"

count_per_sample() {
    local category=$1
    local vcf_dir=$2
    local pattern=$3
    
    echo "Processing: ${category}"
    
    for vcf in ${vcf_dir}/${pattern}; do
        if [[ -f "${vcf}" && ! "${vcf}" =~ multisample ]]; then
            sample=$(basename "${vcf}" | sed 's/_unique.*//')
            
            # Total variants
            total=$(bcftools view -H "${vcf}" | wc -l | tr -d ' \n')
            total=$((total))
            
            # Count SNPs - UPPERCASE pattern
            snps=$(bcftools query -f '%TYPE\n' "${vcf}" 2>/dev/null | grep -c "^SNP$" 2>/dev/null | tr -d ' \n' || echo 0)
            snps=$((snps))
            
            # Count MNPs - UPPERCASE pattern
            mnps=$(bcftools query -f '%TYPE\n' "${vcf}" 2>/dev/null | grep -c "^MNP$" 2>/dev/null | tr -d ' \n' || echo 0)
            mnps=$((mnps))
            
            # Count INDELs - UPPERCASE pattern
            indels=$(bcftools query -f '%TYPE\n' "${vcf}" 2>/dev/null | grep -c "^INDEL$" 2>/dev/null | tr -d ' \n' || echo 0)
            indels=$((indels))
            
            # Count OTHER (OVERLAP, BND, REF, etc.) - everything except SNP, MNP, INDEL
            other=$(bcftools query -f '%TYPE\n' "${vcf}" 2>/dev/null | grep -v "^SNP$\|^MNP$\|^INDEL$" 2>/dev/null | wc -l | tr -d ' \n' || echo 0)
            other=$((other))
            
            calculated=$((snps + mnps + indels + other))
            
            echo -e "${category}\t${sample}\t${snps}\t${mnps}\t${indels}\t${other}\t${total}" >> "${SAMPLE_COUNTS}"
            
            if [[ ${calculated} -eq ${total} ]]; then
                echo "  ${sample}: SNPs=${snps}, MNPs=${mnps}, INDELs=${indels}, OTHER=${other}, Total=${total} ✓"
            else
                echo "  ${sample}: SNPs=${snps}, MNPs=${mnps}, INDELs=${indels}, OTHER=${other}, Total=${total} (calc=${calculated}) ⚠️"
            fi
        fi
    done
    echo ""
}

count_multisample() {
    local vcf_type=$1
    local vcf_file=$2
    
    if [[ ! -f "${vcf_file}" ]]; then
        echo "  WARNING: ${vcf_file} not found, skipping..."
        return
    fi
    
    echo "Processing: ${vcf_type}"
    
    # Total variants
    total=$(bcftools view -H "${vcf_file}" | wc -l | tr -d ' \n')
    total=$((total))
    
    # Count SNPs - UPPERCASE pattern
    snps=$(bcftools query -f '%TYPE\n' "${vcf_file}" 2>/dev/null | grep -c "^SNP$" 2>/dev/null | tr -d ' \n' || echo 0)
    snps=$((snps))
    
    # Count MNPs - UPPERCASE pattern
    mnps=$(bcftools query -f '%TYPE\n' "${vcf_file}" 2>/dev/null | grep -c "^MNP$" 2>/dev/null | tr -d ' \n' || echo 0)
    mnps=$((mnps))
    
    # Count INDELs - UPPERCASE pattern
    indels=$(bcftools query -f '%TYPE\n' "${vcf_file}" 2>/dev/null | grep -c "^INDEL$" 2>/dev/null | tr -d ' \n' || echo 0)
    indels=$((indels))
    
    # Count OTHER (OVERLAP, BND, REF, etc.) - everything except SNP, MNP, INDEL
    other=$(bcftools query -f '%TYPE\n' "${vcf_file}" 2>/dev/null | grep -v "^SNP$\|^MNP$\|^INDEL$" 2>/dev/null | wc -l | tr -d ' \n' || echo 0)
    other=$((other))
    
    calculated=$((snps + mnps + indels + other))
    
    echo -e "${vcf_type}\t${snps}\t${mnps}\t${indels}\t${other}\t${total}" >> "${MULTISAMPLE_COUNTS}"
    
    if [[ ${calculated} -eq ${total} ]]; then
        echo "  ${vcf_type}: SNPs=${snps}, MNPs=${mnps}, INDELs=${indels}, OTHER=${other}, Total=${total} ✓"
    else
        echo "  ${vcf_type}: SNPs=${snps}, MNPs=${mnps}, INDELs=${indels}, OTHER=${other}, Total=${total} (calc=${calculated}) ⚠️"
    fi
    echo ""
}

echo "PART 1: Per-sample counts"
echo "========================================================================"

count_per_sample "GATK_Illumina" "${BASE_DIR}/uniqueGATKillumina" "*_unique_GATK_illumina.vcf.gz"
count_per_sample "PPR_Illumina" "${BASE_DIR}/uniquePPRillumina" "*_unique_PPR_illumina.vcf.gz"
count_per_sample "GATK_Aviti" "${BASE_DIR}/uniqueGATKaviti" "*_unique_GATK_aviti.vcf.gz"
count_per_sample "PPR_Aviti" "${BASE_DIR}/uniquePPRaviti" "*_unique_PPR_aviti.vcf.gz"

echo "PART 2: Multisample counts"
echo "========================================================================"

count_multisample "uniqueGATKillumina" "${BASE_DIR}/uniqueGATKillumina/multisample_uniqueGATKillumina.vcf.gz"
count_multisample "uniquePPRillumina" "${BASE_DIR}/uniquePPRillumina/multisample_uniquePPRillumina.vcf.gz"
count_multisample "uniqueGATKaviti" "${BASE_DIR}/uniqueGATKaviti/multisample_uniqueGATKaviti.vcf.gz"
count_multisample "uniquePPRaviti" "${BASE_DIR}/uniquePPRaviti/multisample_uniquePPRaviti.vcf.gz"
count_multisample "uniquePPR" "${OUTPUT_DIR}/uniquePPR.vcf.gz"
count_multisample "uniqueGATK" "${OUTPUT_DIR}/uniqueGATK.vcf.gz"

echo "========================================================================"
echo "✓ Done! Output files:"
echo "  ${SAMPLE_COUNTS}"
echo "  ${MULTISAMPLE_COUNTS}"
echo ""
echo "Categories: SNPs, MNPs, INDELs, OTHER (+ Total for verification)"
echo "Next: Rscript 03_plot_snvs_indels_UPDATED.R"
