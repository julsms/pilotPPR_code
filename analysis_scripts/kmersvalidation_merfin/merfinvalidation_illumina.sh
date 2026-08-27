#!/bin/bash
#
# Builds the Meryl k-mer database and runs Merfin -filter on the complete
# (ALL variants) GATK and PPR call sets per sample. Shown here for a single
# sample as a worked example of the direct Merfin invocation; run per-sample
# for the full cohort and for Aviti with platform-specific paths substituted.
#
# The Merfin-validated ALL-variants VCFs produced here
# (<sample>_gatk_all_validated.filter.vcf.gz, <sample>_ppr_all_validated.filter.vcf.gz)
# are the input consumed by uniquevariants_validationrate_normalized.sh, which
# computes the unique-variant validation rate by intersection instead of
# re-running Merfin on the unique VCFs directly (see that script's header).

set -e

ILLUMINA_SAMPLES="1ICDG_S1"

PLATFORM="illumina"
OUTDIR="/mnt/genomics/pilot_PPR/kmers_validation/illumina"
GATK_VCF_DIR="/mnt/genomics/pilot_PPR/vcfs/GATK_illumina"
PPR_VCF_DIR="/mnt/genomics/pilot_PPR/vcfs/PPR_illumina"
REFERENCE="/mnt/genomics/illumina/ref_chr/Homo_sapiens.GRCh38.dna.primary_assembly.fasta"
S3_BASE="s3://icdg/pilot1/illumina/processed1"

echo "Illumina batch processing started"
echo "Samples: $(echo $ILLUMINA_SAMPLES | wc -w)"

source /mnt/genomics/tools/miniconda3/etc/profile.d/conda.sh
conda activate meryl_env

mkdir -p ${OUTDIR}/final_results
mkdir -p ${OUTDIR}/vcf_results

for SAMPLE in $ILLUMINA_SAMPLES; do
    echo "--- ${SAMPLE} | $(date) ---"

    WORK_DIR="${OUTDIR}/temp_${SAMPLE}"
    mkdir -p ${WORK_DIR}
    cd ${WORK_DIR}

    # download trimmed FASTQs from S3
    echo "downloading FASTQs for ${SAMPLE}..."
    aws s3 cp ${S3_BASE}/${SAMPLE}/qc/${SAMPLE}_r1_val_1.fq.gz ./${SAMPLE}_R1_val_1.fq.gz
    aws s3 cp ${S3_BASE}/${SAMPLE}/qc/${SAMPLE}_r2_val_2.fq.gz ./${SAMPLE}_R2_val_2.fq.gz

    if [[ ! -f "${SAMPLE}_R1_val_1.fq.gz" ]] || [[ ! -f "${SAMPLE}_R2_val_2.fq.gz" ]]; then
        echo "ERROR: FASTQ download failed for ${SAMPLE}, skipping"
        continue
    fi

    # build meryl k-mer database (k=31), keep k-mers with frequency > 1
    echo "building meryl database for ${SAMPLE}..."
    meryl count k=31 threads=65 memory=145 ${SAMPLE}_R1_val_1.fq.gz ${SAMPLE}_R2_val_2.fq.gz output ${SAMPLE}.meryl
    meryl greater-than 1 ${SAMPLE}.meryl output ${SAMPLE}.gt1.meryl
    rm -rf ${SAMPLE}.meryl

    # extract non-ref genotypes per sample from multisample VCFs
    echo "extracting variants for ${SAMPLE}..."

    bcftools view -s ${SAMPLE} -i 'GT!="0/0" && GT!="./." && GT!="0|0" && GT!=".|."' \
        ${GATK_VCF_DIR}/${SAMPLE}.GATK.filtered.vcf.gz | \
        bcftools view -Oz -o ${SAMPLE}_gatk_all.vcf.gz
    tabix -p vcf ${SAMPLE}_gatk_all.vcf.gz

    bcftools view -s ${SAMPLE} -i 'GT!="0/0" && GT!="./." && GT!="0|0" && GT!=".|."' \
        ${PPR_VCF_DIR}/${SAMPLE}.PPR.filtered.vcf.gz | \
        bcftools view -Oz -o ${SAMPLE}_ppr_all.vcf.gz
    tabix -p vcf ${SAMPLE}_ppr_all.vcf.gz

    # count variants before merfin filtering
    GATK_ALL_TOTAL=$(bcftools view -H ${SAMPLE}_gatk_all.vcf.gz | wc -l)
    PPR_ALL_TOTAL=$(bcftools view -H ${SAMPLE}_ppr_all.vcf.gz | wc -l)

    echo "GATK all: ${GATK_ALL_TOTAL} | PPR all: ${PPR_ALL_TOTAL}"

    # merfin -filter on GATK all variants
    if [[ ${GATK_ALL_TOTAL} -gt 0 ]]; then
        merfin -filter \
            -sequence ${REFERENCE} \
            -readmers ${SAMPLE}.gt1.meryl \
            -vcf ${SAMPLE}_gatk_all.vcf.gz \
            -output ${SAMPLE}_gatk_all_validated \
            -threads 60

        if [[ -f "${SAMPLE}_gatk_all_validated.filter.vcf" ]]; then
            GATK_ALL_VALIDATED=$(bcftools view -H ${SAMPLE}_gatk_all_validated.filter.vcf | wc -l)
            GATK_ALL_RATE=$(echo "scale=4; ${GATK_ALL_VALIDATED} * 100 / ${GATK_ALL_TOTAL}" | bc -l)
            bgzip ${SAMPLE}_gatk_all_validated.filter.vcf
            tabix -p vcf ${SAMPLE}_gatk_all_validated.filter.vcf.gz
        else
            GATK_ALL_VALIDATED=0
            GATK_ALL_RATE=0
        fi
    else
        GATK_ALL_VALIDATED=0
        GATK_ALL_RATE=0
    fi

    # merfin -filter on PPR all variants
    if [[ ${PPR_ALL_TOTAL} -gt 0 ]]; then
        merfin -filter \
            -sequence ${REFERENCE} \
            -readmers ${SAMPLE}.gt1.meryl \
            -vcf ${SAMPLE}_ppr_all.vcf.gz \
            -output ${SAMPLE}_ppr_all_validated \
            -threads 60

        if [[ -f "${SAMPLE}_ppr_all_validated.filter.vcf" ]]; then
            PPR_ALL_VALIDATED=$(bcftools view -H ${SAMPLE}_ppr_all_validated.filter.vcf | wc -l)
            PPR_ALL_RATE=$(echo "scale=4; ${PPR_ALL_VALIDATED} * 100 / ${PPR_ALL_TOTAL}" | bc -l)
            bgzip ${SAMPLE}_ppr_all_validated.filter.vcf
            tabix -p vcf ${SAMPLE}_ppr_all_validated.filter.vcf.gz
        else
            PPR_ALL_VALIDATED=0
            PPR_ALL_RATE=0
        fi
    else
        PPR_ALL_VALIDATED=0
        PPR_ALL_RATE=0
    fi

    # save VCFs
    mkdir -p ${OUTDIR}/vcf_results/${SAMPLE}

    cp ${SAMPLE}_gatk_all.vcf.gz* ${OUTDIR}/vcf_results/${SAMPLE}/
    cp ${SAMPLE}_ppr_all.vcf.gz* ${OUTDIR}/vcf_results/${SAMPLE}/

    [[ -f "${SAMPLE}_gatk_all_validated.filter.vcf.gz" ]] && cp ${SAMPLE}_gatk_all_validated.filter.vcf.gz* ${OUTDIR}/vcf_results/${SAMPLE}/
    [[ -f "${SAMPLE}_ppr_all_validated.filter.vcf.gz" ]] && cp ${SAMPLE}_ppr_all_validated.filter.vcf.gz* ${OUTDIR}/vcf_results/${SAMPLE}/

    # write per-sample CSV
    echo "Pipeline,Variant_Type,Platform,Sample,Total_Variants,Validated_Variants,Validation_Rate" > ${OUTDIR}/final_results/${PLATFORM}_${SAMPLE}_validation_results.csv
    echo "GATK,ALL,${PLATFORM},${SAMPLE},${GATK_ALL_TOTAL},${GATK_ALL_VALIDATED},${GATK_ALL_RATE}" >> ${OUTDIR}/final_results/${PLATFORM}_${SAMPLE}_validation_results.csv
    echo "PPR,ALL,${PLATFORM},${SAMPLE},${PPR_ALL_TOTAL},${PPR_ALL_VALIDATED},${PPR_ALL_RATE}" >> ${OUTDIR}/final_results/${PLATFORM}_${SAMPLE}_validation_results.csv

    # write summary text
    {
        echo "platform: ${PLATFORM}"
        echo "sample: ${SAMPLE}"
        echo "date: $(date)"
        echo ""
        echo "all variants"
        echo "  GATK: ${GATK_ALL_VALIDATED}/${GATK_ALL_TOTAL} (${GATK_ALL_RATE}%)"
        echo "  PPR:  ${PPR_ALL_VALIDATED}/${PPR_ALL_TOTAL} (${PPR_ALL_RATE}%)"
    } > ${OUTDIR}/final_results/summary_${PLATFORM}_${SAMPLE}.txt

    cd ${OUTDIR}
    rm -rf ${WORK_DIR}

    echo "done ${SAMPLE}"
    echo "  GATK all: ${GATK_ALL_VALIDATED}/${GATK_ALL_TOTAL} (${GATK_ALL_RATE}%)"
    echo "  PPR  all: ${PPR_ALL_VALIDATED}/${PPR_ALL_TOTAL} (${PPR_ALL_RATE}%)"
done

echo "finished. results in ${OUTDIR}/final_results/"
