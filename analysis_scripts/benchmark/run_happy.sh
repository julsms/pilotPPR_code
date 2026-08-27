#!/bin/bash

VCF_DIR="/mnt/genomics/HG002Aviti_SRR33551004/benchGATK_v5_checktest"
REFERENCE="/mnt/genomics/illumina/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
OUTPUT_BASE="/mnt/genomics/HG002Aviti_SRR33551004/benchGATK_v5_checktest/happy"

mkdir -p "${OUTPUT_BASE}"

REGIONS=(
    "AD" "15GC85" "25GC65" "BB" "LM"
    "OD" "REFSEQ" "TRHP" "ACMG" "OMIM"
)

for REGION in "${REGIONS[@]}"; do
    echo "--- ${REGION} ---"

    GOLD_VCF="${VCF_DIR}/gold_${REGION}.vcf"
    QUERY_VCF="${VCF_DIR}/query_${REGION}.vcf"
    REGION_OUTPUT="${OUTPUT_BASE}/${REGION}"
    mkdir -p "${REGION_OUTPUT}"

    if [ ! -f "${GOLD_VCF}" ]; then
        echo "gold VCF not found: ${GOLD_VCF}, skipping"
        continue
    fi

    if [ ! -f "${QUERY_VCF}" ]; then
        echo "query VCF not found: ${QUERY_VCF}, skipping"
        continue
    fi

    hap.py \
        "${GOLD_VCF}" \
        "${QUERY_VCF}" \
        -r "${REFERENCE}" \
        -o "${REGION_OUTPUT}/${REGION}" \
        --threads 8 \
        --engine=vcfeval \
        --verbose

    if [ $? -eq 0 ]; then
        echo "done: ${REGION}"
    else
        echo "failed: ${REGION}"
    fi
done

echo "finished. results in ${OUTPUT_BASE}"
ls -lh "${OUTPUT_BASE}"/*/*.summary.csv 2>/dev/null | awk '{print $9}'

