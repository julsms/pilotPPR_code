#!/bin/bash
set -euo pipefail

# ==============================================================================
# GenomeScope Runner — AVITI
# ==============================================================================
# Same logic as run_illumina.sh but for Aviti platform.
# Handles the R1/r1 naming inconsistency in S3.
# Pairs Aviti samples with Illumina counterparts for coverage matching.
#
# Usage:
#   ./run_aviti.sh
#
# Output:
#   ./genomescope_results_aviti/<sample>/genomescope2_output/
#
# Requires: Docker, Nextflow >= 24.04, AWS CLI
# ==============================================================================

# --- Configuration ---
COVERAGE_TSV="/mnt/genomics/pilot_PPR/genomescope/qc_metrics_collected_for_plots.tsv"
S3_BUCKET="s3://icdg/pilot1/aviti/processed1"
OUTPUT_DIR="./genomescope_results_aviti"
TEMP_DIR="./temp_fastq_aviti"
GENOME_SIZE=3100000000  # 3.1 Gb
PIPELINE_DIR="$(cd "$(dirname "$0")" && pwd)"

# Common S3 options
s3_opts=(--endpoint-url https://s3.icdg.auxilio.ai --region us-east-1 --profile s3-genomics)

# Large file transfer settings (30+ GB files need longer timeouts)
export AWS_MAX_ATTEMPTS=5
export AWS_RETRY_MODE=adaptive

mkdir -p "$OUTPUT_DIR" "$TEMP_DIR"

# --- Sample pairing table ---
# Aviti_id -> Illumina_num (for coverage lookup)
# Note: Aviti ICDG21 = Illumina 12, Aviti ICDG12 = Illumina 21
declare -A PAIR_MAP
PAIR_MAP=(
    [ICDG1]=1   [ICDG2]=2   [ICDG3]=3   [ICDG4]=4   [ICDG5]=5
    [ICDG6]=6   [ICDG7]=7   [ICDG8]=8   [ICDG9]=9   [ICDG10]=10
    [ICDG11]=11 [ICDG12]=21  [ICDG13]=13 [ICDG14]=14 [ICDG15]=15
    [ICDG16]=16 [ICDG17]=17 [ICDG18]=18 [ICDG19]=19 [ICDG20]=20
    [ICDG21]=12  [ICDG22]=22 [ICDG24]=24
)

# --- Read coverage from TSV into arrays ---
declare -A ILLUMINA_COV
declare -A AVITI_COV

echo "Reading coverage from: $COVERAGE_TSV"
while IFS=$'\t' read -r platform sample cov cv; do
    [[ "$platform" == "platform" ]] && continue  # skip header
    if [[ "$platform" == "Illumina" ]]; then
        NUM=$(echo "$sample" | grep -oP '^\d+')
        ILLUMINA_COV["$NUM"]="$cov"
    elif [[ "$platform" == "AVITI" ]]; then
        AVITI_COV["$sample"]="$cov"
    fi
done < "$COVERAGE_TSV"

echo "Loaded ${#ILLUMINA_COV[@]} Illumina and ${#AVITI_COV[@]} Aviti coverage values."
echo ""

# --- Process each Aviti sample ---
echo "================================================================"
echo " Processing Aviti samples (with coverage-matched subsampling)"
echo " Output: $OUTPUT_DIR"
echo "================================================================"
echo ""

for SAMPLE in $(echo "${!AVITI_COV[@]}" | tr ' ' '\n' | sort); do
    AVITI_COV_VAL=${AVITI_COV[$SAMPLE]}

    # Get paired Illumina sample number and its coverage
    ILLUM_NUM=${PAIR_MAP[$SAMPLE]:-""}
    if [[ -z "$ILLUM_NUM" ]]; then
        echo "  WARNING: No Illumina pair for Aviti sample $SAMPLE. Skipping."
        continue
    fi
    ILLUM_COV=${ILLUMINA_COV[$ILLUM_NUM]:-""}
    if [[ -z "$ILLUM_COV" ]]; then
        echo "  WARNING: No Illumina coverage for sample $ILLUM_NUM. Skipping."
        continue
    fi

    # Determine target coverage = min(illumina, aviti)
    TARGET_COV=$(python3 -c "print(min(${AVITI_COV_VAL}, ${ILLUM_COV}))")
    TARGET_BASES=$(python3 -c "import math; print(int(math.floor(${TARGET_COV} * ${GENOME_SIZE})))")

    echo "========== ${SAMPLE} =========="
    echo "  Aviti coverage: ${AVITI_COV_VAL}x | Illumina pair (${ILLUM_NUM}ICDG): ${ILLUM_COV}x"
    echo "  Target: ${TARGET_COV}x → ${TARGET_BASES} bases"

    # Skip if already processed
    if [[ -d "${OUTPUT_DIR}/${SAMPLE}/genomescope2_output" ]]; then
        echo "  ✓ Already processed. Skipping."
        echo ""
        continue
    fi

    # --- Detect actual FASTQ filenames (handles both r1/r2 and R1/R2) ---
    S3_QC_PATH="${S3_BUCKET}/${SAMPLE}/qc/"
    echo "  Listing: ${S3_QC_PATH}"
    QC_LISTING=$(aws s3 ls "$S3_QC_PATH" "${s3_opts[@]}") || {
        echo "  ERROR: Could not list ${S3_QC_PATH}. Skipping."
        echo ""
        continue
    }

    # Find the val_1 file (could be r1_val_1, R1_val_1, or just _val_1)
    R1_FILENAME=$(echo "$QC_LISTING" | grep -oP '\S+[rR]1_val_1\.fq\.gz$' | head -1 || true)
    R2_FILENAME=$(echo "$QC_LISTING" | grep -oP '\S+[rR]2_val_2\.fq\.gz$' | head -1 || true)

    # Fallback: try pattern without r1/R1 prefix (e.g., ICDG1_val_1.fq.gz)
    if [[ -z "$R1_FILENAME" ]]; then
        R1_FILENAME=$(echo "$QC_LISTING" | grep -oP '\S+_val_1\.fq\.gz$' | head -1 || true)
    fi
    if [[ -z "$R2_FILENAME" ]]; then
        R2_FILENAME=$(echo "$QC_LISTING" | grep -oP '\S+_val_2\.fq\.gz$' | head -1 || true)
    fi

    if [[ -z "$R1_FILENAME" || -z "$R2_FILENAME" ]]; then
        echo "  ERROR: Could not find val_1/val_2 FASTQ files in ${S3_QC_PATH}"
        echo "  Listing:"
        echo "$QC_LISTING" | sed 's/^/    /'
        echo "  Skipping."
        echo ""
        continue
    fi

    echo "  R1: ${R1_FILENAME}"
    echo "  R2: ${R2_FILENAME}"

    R1="${TEMP_DIR}/${R1_FILENAME}"
    R2="${TEMP_DIR}/${R2_FILENAME}"

    if [[ -f "$R1" && -f "$R2" ]]; then
        echo "  ✓ FASTQ files already local."
    else
        echo "  Downloading: ${R1_FILENAME}, ${R2_FILENAME} (large files, may take a while)..."
        aws s3 cp "${S3_QC_PATH}${R1_FILENAME}" "$R1" "${s3_opts[@]}" \
            --cli-read-timeout 0 --cli-connect-timeout 30 || {
            echo "  ERROR: Failed to download ${R1_FILENAME}. Skipping."
            rm -f "$R1"
            echo ""
            continue
        }
        aws s3 cp "${S3_QC_PATH}${R2_FILENAME}" "$R2" "${s3_opts[@]}" \
            --cli-read-timeout 0 --cli-connect-timeout 30 || {
            echo "  ERROR: Failed to download ${R2_FILENAME}. Skipping."
            rm -f "$R1" "$R2"
            echo ""
            continue
        }
        echo "  ✓ Downloaded"
    fi

    # Determine if subsampling is needed
    SUBSAMPLE_ARGS=()
    if python3 -c "exit(0 if ${AVITI_COV_VAL} > ${TARGET_COV} else 1)"; then
        SUBSAMPLE_ARGS=(--subsample true --target_bases "${TARGET_BASES}")
        echo "  → Will subsample (${AVITI_COV_VAL}x > ${TARGET_COV}x)"
    else
        echo "  → No subsampling needed (already at target)"
    fi

    # Run pipeline
    echo "  Running FastK + GenomeScope..."
    nextflow run "${PIPELINE_DIR}/main.nf" \
        -c "${PIPELINE_DIR}/nextflow.config" \
        --sample_id "$SAMPLE" \
        --fastq_r1 "$R1" \
        --fastq_r2 "$R2" \
        --output_dir "$OUTPUT_DIR" \
        "${SUBSAMPLE_ARGS[@]+"${SUBSAMPLE_ARGS[@]}"}"

    # Cleanup
    rm -rf work/
    rm -f "$R1" "$R2"

    echo "  ✓ ${SAMPLE} done"
    echo ""
done

echo "======================================="
echo " All Aviti samples processed!"
echo " Results: $OUTPUT_DIR"
echo "======================================="
