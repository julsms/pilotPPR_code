#!/bin/bash
set -euo pipefail

# ==============================================================================
# GenomeScope Runner — ILLUMINA
# ==============================================================================
# Reads coverage info from qc_metrics_collected_for_plots.tsv, pairs each
# Illumina sample with its Aviti counterpart, determines the subsampling
# target (min coverage of the pair × genome_size), then runs FastK + GenomeScope.
#
# Subsampling logic:
#   target_cov = min(illumina_cov, aviti_cov)
#   - If illumina_cov > target_cov → subsample Illumina reads to target_bases
#   - If illumina_cov ≤ target_cov → run without subsampling (Illumina already
#     at or below the paired minimum)
#
# Usage:
#   ./run_illumina.sh
#
# Output:
#   ./genomescope_results_illumina/<sample>/genomescope2_output/
#
# Requires: Docker, Nextflow >= 24.04, AWS CLI
# ==============================================================================

# --- Configuration ---
COVERAGE_TSV="/mnt/genomics/pilot_PPR/genomescope/qc_metrics_collected_for_plots.tsv"
S3_BUCKET="s3://icdg/pilot1/illumina/processed1"
OUTPUT_DIR="./genomescope_results_illumina"
TEMP_DIR="./temp_fastq_illumina"
GENOME_SIZE=3100000000  # 3.1 Gb
PIPELINE_DIR="$(cd "$(dirname "$0")" && pwd)"

# Common S3 options
s3_opts=(--endpoint-url https://s3.icdg.auxilio.ai --region us-east-1 --profile s3-genomics)

mkdir -p "$OUTPUT_DIR" "$TEMP_DIR"

# --- Sample pairing table ---
# Illumina_base_id -> Aviti_id
# Note: Illumina 12 = Aviti 21, Illumina 21 = Aviti 12
declare -A PAIR_MAP
PAIR_MAP=(
    [1]=ICDG1   [2]=ICDG2   [3]=ICDG3   [4]=ICDG4   [5]=ICDG5
    [6]=ICDG6   [7]=ICDG7   [8]=ICDG8   [9]=ICDG9   [10]=ICDG10
    [11]=ICDG11 [12]=ICDG21  [13]=ICDG13 [14]=ICDG14 [15]=ICDG15
    [16]=ICDG16 [17]=ICDG17 [18]=ICDG18 [19]=ICDG19 [20]=ICDG20
    [21]=ICDG12  [22]=ICDG22 [24]=ICDG24
)

# --- S3 folder name lookup table (NUM -> full folder name with _Sx suffix) ---
declare -A S3_FOLDER
S3_FOLDER=(
    [1]=1ICDG_S1     [2]=2ICDG_S4     [3]=3ICDG_S7     [4]=4ICDG_S10
    [5]=5ICDG_S13    [6]=6ICDG_S16    [7]=7ICDG_S19    [8]=8ICDG_S22
    [9]=9ICDG_S2     [10]=10ICDG_S5   [11]=11ICDG_S8   [12]=12ICDG_S11
    [13]=13ICDG_S14  [14]=14ICDG_S17  [15]=15ICDG_S20  [16]=16ICDG_S23
    [17]=17ICDG_S3   [18]=18ICDG_S6   [19]=19ICDG_S9   [20]=20ICDG_S12
    [21]=21ICDG_S15  [22]=22ICDG_S18  [23]=23ICDG_S21  [24]=24ICDG_S24
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

# --- Process each Illumina sample ---
echo "================================================================"
echo " Processing Illumina samples (with coverage-matched subsampling)"
echo " Subsampling tool: rasusa (single-pass, multi-threaded)"
echo " Output: $OUTPUT_DIR"
echo "================================================================"
echo ""

for NUM in $(echo "${!ILLUMINA_COV[@]}" | tr ' ' '\n' | sort -n); do
    ILLUM_COV=${ILLUMINA_COV[$NUM]}

    # Get paired Aviti sample and its coverage
    AVITI_ID=${PAIR_MAP[$NUM]:-""}
    if [[ -z "$AVITI_ID" ]]; then
        echo "  WARNING: No Aviti pair for Illumina sample $NUM. Skipping."
        continue
    fi
    AVITI_COV_VAL=${AVITI_COV[$AVITI_ID]:-""}
    if [[ -z "$AVITI_COV_VAL" ]]; then
        echo "  WARNING: No Aviti coverage for $AVITI_ID. Skipping."
        continue
    fi

    # Determine target coverage = min(illumina, aviti)
    TARGET_COV=$(python3 -c "print(min(${ILLUM_COV}, ${AVITI_COV_VAL}))")
    TARGET_BASES=$(python3 -c "import math; print(int(math.floor(${TARGET_COV} * ${GENOME_SIZE})))")

    # Get the S3 folder name (with _Sx suffix)
    S3_SAMPLE_FOLDER=${S3_FOLDER[$NUM]:-""}
    if [[ -z "$S3_SAMPLE_FOLDER" ]]; then
        echo "  WARNING: No S3 folder mapping for sample $NUM. Skipping."
        continue
    fi

    SAMPLE="${NUM}ICDG"

    echo "========== ${SAMPLE} (S3: ${S3_SAMPLE_FOLDER}) =========="
    echo "  Illumina coverage: ${ILLUM_COV}x | Aviti pair ($AVITI_ID): ${AVITI_COV_VAL}x"
    echo "  Target: ${TARGET_COV}x → ${TARGET_BASES} bases"

    # Skip if already processed
    if [[ -d "${OUTPUT_DIR}/${SAMPLE}/genomescope2_output" ]]; then
        echo "  ✓ Already processed. Skipping."
        echo ""
        continue
    fi

    # --- Detect actual FASTQ filenames (handles both r1/r2 and R1/R2) ---
    S3_QC_PATH="${S3_BUCKET}/${S3_SAMPLE_FOLDER}/qc/"
    QC_LISTING=$(aws s3 ls "$S3_QC_PATH" "${s3_opts[@]}")

    # Find the val_1 file (could be r1_val_1 or R1_val_1)
    R1_FILENAME=$(echo "$QC_LISTING" | grep -oP '\S+[rR]1_val_1\.fq\.gz$' | head -1)
    R2_FILENAME=$(echo "$QC_LISTING" | grep -oP '\S+[rR]2_val_2\.fq\.gz$' | head -1)

    if [[ -z "$R1_FILENAME" || -z "$R2_FILENAME" ]]; then
        echo "  ERROR: Could not find val_1/val_2 FASTQ files in ${S3_QC_PATH}"
        echo "  Listing: $QC_LISTING"
        echo "  Skipping."
        echo ""
        continue
    fi

    R1="${TEMP_DIR}/${R1_FILENAME}"
    R2="${TEMP_DIR}/${R2_FILENAME}"

    if [[ -f "$R1" && -f "$R2" ]]; then
        echo "  ✓ FASTQ files already local."
    else
        echo "  Downloading from S3: ${R1_FILENAME}, ${R2_FILENAME}"
        aws s3 cp "${S3_QC_PATH}${R1_FILENAME}" "$R1" --quiet "${s3_opts[@]}"
        aws s3 cp "${S3_QC_PATH}${R2_FILENAME}" "$R2" --quiet "${s3_opts[@]}"
        echo "  ✓ Downloaded"
    fi

    # Determine if subsampling is needed
    SUBSAMPLE_PARAMS=""
    if python3 -c "exit(0 if ${ILLUM_COV} > ${TARGET_COV} else 1)"; then
        SUBSAMPLE_PARAMS="--subsample true --target_bases ${TARGET_BASES}"
        echo "  → Will subsample (${ILLUM_COV}x > ${TARGET_COV}x)"
    else
        echo "  → No subsampling needed (already at or below target)"
    fi

    # Run pipeline
    echo "  Running rasusa (if needed) + FastK + GenomeScope..."
    nextflow run "${PIPELINE_DIR}/main.nf" \
        -c "${PIPELINE_DIR}/nextflow.config" \
        --sample_id "$SAMPLE" \
        --fastq_r1 "$R1" \
        --fastq_r2 "$R2" \
        --output_dir "$OUTPUT_DIR" \
        $SUBSAMPLE_PARAMS 

    # Cleanup
    rm -rf work/
    rm -f "$R1" "$R2"

    echo "  ✓ ${SAMPLE} done"
    echo ""
done

echo "======================================="
echo " All Illumina samples processed!"
echo " Results: $OUTPUT_DIR"
echo "======================================="
