#!/bin/bash
set -euo pipefail

# ==============================================================================
# Stratified Merfin Validation Rates
# ==============================================================================
# Calculates k-mer validation rates per sample, per pipeline, per genomic
# stratification region (BED files).
#
# Logic:
#   ALL VARIANTS:
#     For each sample × pipeline × BED region:
#       total_in_region = bcftools view original.vcf -R region.bed | count variants
#       validated_in_region = bcftools view merfin_output.vcf -R region.bed | count variants
#       validation_rate = validated_in_region / total_in_region * 100
#
#   UNIQUE VARIANTS (corrected logic using bcftools isec):
#     For each sample × pipeline × BED region:
#       total_unique_in_region = variants in unique.vcf overlapping region.bed
#       validated_unique_in_region = unique variants ALSO present in the
#           Merfin-validated complete VCF (all_validated.filter.vcf.gz),
#           determined via bcftools isec intersection within the region
#       validation_rate = validated_unique_in_region / total_unique_in_region * 100
#
# Processes BOTH platforms (Illumina + Aviti) and BOTH variant sets (ALL + UNIQUE).
#
# File structure:
#   ALL variants (original + validated):
#     Illumina: /mnt/genomics/pilot_PPR/kmers_validation/illumina/vcf_results/<SAMPLE_DIR>/
#       <SAMPLE_DIR>_gatk_all.vcf.gz  /  <SAMPLE_DIR>_gatk_all_validated.filter.vcf.gz
#       <SAMPLE_DIR>_ppr_all.vcf.gz   /  <SAMPLE_DIR>_ppr_all_validated.filter.vcf.gz
#     Aviti: /mnt/genomics/pilot_PPR/kmers_validation/aviti/vcf_results/<SAMPLE>/
#       <SAMPLE>_gatk_all_chr.vcf.gz  /  <SAMPLE>_gatk_all_validated.filter.vcf.gz
#       <SAMPLE>_ppr_all.vcf.gz       /  <SAMPLE>_ppr_all_validated.filter.vcf.gz
#
#   UNIQUE variants (original from /uniquevariants/, validated via isec with all_validated):
#     Original unique VCFs:
#       /mnt/genomics/pilot_PPR/uniquevariants/uniqueGATKillumina/<SAMPLE>_unique_GATK_illumina.vcf.gz
#       /mnt/genomics/pilot_PPR/uniquevariants/uniquePPRillumina/<SAMPLE>_unique_PPR_illumina.vcf.gz
#       /mnt/genomics/pilot_PPR/uniquevariants/uniqueGATKaviti/<SAMPLE>_unique_GATK_aviti.vcf.gz
#       /mnt/genomics/pilot_PPR/uniquevariants/uniquePPRaviti/<SAMPLE>_unique_PPR_aviti.vcf.gz
#     Validation is computed by intersecting unique VCFs with the corresponding
#     Merfin-validated complete VCF (all_validated.filter.vcf.gz).
#
# Output: CSV files in /mnt/genomics/pilot_PPR/kmers_validation/validationrates_stratified/
#
# Usage:
#   ./stratified_merfin_validation.sh
#
# Requirements: bcftools (with index support), python3
# ==============================================================================

# ===================== CONFIGURATION =====================

WORK_DIR="/mnt/genomics/pilot_PPR/kmers_validation/validationrates_stratified"
mkdir -p "$WORK_DIR"

# VCF directories
ILLUMINA_VCF_BASE="/mnt/genomics/pilot_PPR/kmers_validation/illumina/vcf_results"
AVITI_VCF_BASE="/mnt/genomics/pilot_PPR/kmers_validation/aviti/vcf_results"

# Unique variants (original, before Merfin)
UNIQUE_BASE="/mnt/genomics/pilot_PPR/uniquevariants"

# Stratification BED files (GIAB v3.6) — same as concordance analysis
STRAT_BED_DIR="/mnt/genomics/GIAB_stratifications/v3.6/with_chr"

# Temporary directory for bcftools isec operations
TMPDIR_BASE="${WORK_DIR}/tmp_isec"
mkdir -p "$TMPDIR_BASE"

# ===================== BED FILES MAPPING =====================

declare -A REGION_BEDS
REGION_BEDS=(
    [BB]="${STRAT_BED_DIR}/BB.bed"
    [LM]="${STRAT_BED_DIR}/LM.bed"
    [TRHP]="${STRAT_BED_DIR}/TRHP.bed"
    [OD]="${STRAT_BED_DIR}/OD.bed"
    [15GC85]="${STRAT_BED_DIR}/15GC85.bed"
    [25GC65]="${STRAT_BED_DIR}/25GC65.bed"
    [AD]="${STRAT_BED_DIR}/AD.bed"
    [ACMG]="${STRAT_BED_DIR}/ACMG.bed"
    [OMIM]="${STRAT_BED_DIR}/OMIM.bed"
    [REFSEQ]="${STRAT_BED_DIR}/REFSEQ.bed"
)

# ===================== SAMPLE LISTS =====================
# Illumina sample directories (as they appear in vcf_results/)
ILLUMINA_SAMPLE_DIRS=(
    1ICDG_S1 2ICDG_S4 3ICDG_S7 4ICDG_S10 5ICDG_S13 6ICDG_S16 7ICDG_S19 8ICDG_S22
    9ICDG_S2 10ICDG_S5 11ICDG_S8 12ICDG_S11 13ICDG_S14 14ICDG_S17 15ICDG_S20
    16ICDG_S23 17ICDG_S3 18ICDG_S6 19ICDG_S9 20ICDG_S12 21ICDG_S15 22ICDG_S18 24ICDG_S24
)

# Aviti sample directories
AVITI_SAMPLES=(
    ICDG1 ICDG2 ICDG3 ICDG4 ICDG5 ICDG6 ICDG7 ICDG8 ICDG9 ICDG10
    ICDG11 ICDG12 ICDG13 ICDG14 ICDG15 ICDG16 ICDG17 ICDG18 ICDG19
    ICDG20 ICDG21 ICDG22 ICDG24
)

# ===================== HELPER FUNCTIONS =====================

count_variants() {
    # Count variants in a VCF restricted to a BED region.
    # Does NOT apply PASS filter (Merfin output is already filtered).
    # Args: $1 = VCF file, $2 = BED file
    local vcf="$1"
    local bed="$2"

    if [[ ! -f "$vcf" ]]; then
        echo "0"
        return
    fi

    local count
    count=$(bcftools view -R "$bed" -H "$vcf" 2>/dev/null | wc -l || echo "0")
    echo "$count"
}

count_validated_unique_in_region() {
    # Count how many variants from the unique VCF are also present in the
    # Merfin-validated complete VCF, restricted to a BED region.
    #
    # Logic: intersect unique_vcf ∩ validated_all_vcf within the BED region.
    # A unique variant is "validated" if Merfin kept it in the complete set.
    #
    # Args: $1 = unique VCF, $2 = all_validated VCF, $3 = BED file, $4 = tmp subdir
    local unique_vcf="$1"
    local validated_all_vcf="$2"
    local bed="$3"
    local tmpdir="$4"

    if [[ ! -f "$unique_vcf" ]] || [[ ! -f "$validated_all_vcf" ]]; then
        echo "0"
        return
    fi

    # Clean tmp
    rm -rf "${tmpdir}"
    mkdir -p "${tmpdir}"

    # Extract unique variants in this region
    bcftools view -R "$bed" "$unique_vcf" -Oz -o "${tmpdir}/unique_region.vcf.gz" 2>/dev/null
    bcftools index "${tmpdir}/unique_region.vcf.gz" 2>/dev/null

    # Extract validated (complete) variants in this region
    bcftools view -R "$bed" "$validated_all_vcf" -Oz -o "${tmpdir}/validated_region.vcf.gz" 2>/dev/null
    bcftools index "${tmpdir}/validated_region.vcf.gz" 2>/dev/null

    # Intersect: find unique variants that are present in the validated set
    # bcftools isec -n=2 outputs variants present in both files
    # -w1 writes output from the first file (unique) that matched
    bcftools isec -n=2 -w1 \
        "${tmpdir}/unique_region.vcf.gz" \
        "${tmpdir}/validated_region.vcf.gz" \
        -Oz -o "${tmpdir}/matched.vcf.gz" 2>/dev/null

    local count
    count=$(bcftools view -H "${tmpdir}/matched.vcf.gz" 2>/dev/null | wc -l || echo "0")

    # Cleanup
    rm -rf "${tmpdir}"

    echo "$count"
}

# Extract the "ICDG" sample name from Illumina directory name
# e.g., "10ICDG_S5" -> "ICDG10"
illumina_dir_to_icdg() {
    local dir_name="$1"
    # Extract number before ICDG
    local num="${dir_name%%ICDG*}"
    echo "ICDG${num}"
}

# ===================== MAIN =====================

echo "================================================================"
echo " Stratified Merfin Validation Rates"
echo " Output: ${WORK_DIR}"
echo "================================================================"
echo ""

# --- Verify BED files ---
echo "Checking BED files..."
MISSING_BEDS=0
for region in "${!REGION_BEDS[@]}"; do
    bed_file="${REGION_BEDS[$region]}"
    if [[ ! -f "$bed_file" ]]; then
        echo "  ✗ ${region}: NOT FOUND — ${bed_file}"
        MISSING_BEDS=$((MISSING_BEDS + 1))
    else
        echo "  ✓ ${region}"
    fi
done
echo ""

if [[ $MISSING_BEDS -gt 0 ]]; then
    echo "WARNING: $MISSING_BEDS BED file(s) not found. Those regions will be skipped."
    echo "Please adjust REGION_BEDS paths if needed."
    echo ""
fi

# ===========================================================================
# PART 1: ALL VARIANTS — ILLUMINA
# ===========================================================================
echo "========== ILLUMINA: ALL VARIANTS =========="

ALL_ILL_OUTPUT="${WORK_DIR}/stratified_validation_ALL_illumina.csv"
echo "Sample,Pipeline,Platform,Region,Total_Variants,Validated_Variants,Validation_Rate" > "$ALL_ILL_OUTPUT"

for SAMPLE_DIR in "${ILLUMINA_SAMPLE_DIRS[@]}"; do
    # Check directory exists
    DIR_PATH="${ILLUMINA_VCF_BASE}/${SAMPLE_DIR}"
    if [[ ! -d "$DIR_PATH" ]]; then
        echo "  SKIP: $SAMPLE_DIR (directory not found)"
        continue
    fi

    SAMPLE_ID=$(illumina_dir_to_icdg "$SAMPLE_DIR")
    echo "  ${SAMPLE_DIR} (${SAMPLE_ID})..."

    for PIPELINE in gatk ppr; do
        PIPELINE_UPPER=$(echo "$PIPELINE" | tr '[:lower:]' '[:upper:]')

        # Original VCF
        ORIGINAL="${DIR_PATH}/${SAMPLE_DIR}_${PIPELINE}_all.vcf.gz"
        # Validated VCF (Merfin output)
        VALIDATED_VCF="${DIR_PATH}/${SAMPLE_DIR}_${PIPELINE}_all_validated.filter.vcf.gz"

        if [[ ! -f "$ORIGINAL" ]]; then
            echo "    WARNING: Missing original: $ORIGINAL"
            continue
        fi
        if [[ ! -f "$VALIDATED_VCF" ]]; then
            echo "    WARNING: Missing validated: $VALIDATED_VCF"
            continue
        fi

        for REGION in "${!REGION_BEDS[@]}"; do
            BED_FILE="${REGION_BEDS[$REGION]}"
            [[ ! -f "$BED_FILE" ]] && continue

            TOTAL=$(count_variants "$ORIGINAL" "$BED_FILE")
            VALIDATED=$(count_variants "$VALIDATED_VCF" "$BED_FILE")

            if [[ "$TOTAL" -gt 0 ]]; then
                RATE=$(python3 -c "print(round(${VALIDATED}/${TOTAL}*100, 4))")
            else
                RATE="NA"
            fi

            echo "${SAMPLE_ID},${PIPELINE_UPPER},Illumina,${REGION},${TOTAL},${VALIDATED},${RATE}" >> "$ALL_ILL_OUTPUT"
        done
    done
done

echo "  ✓ Saved: $ALL_ILL_OUTPUT"
echo ""

# ===========================================================================
# PART 2: ALL VARIANTS — AVITI
# ===========================================================================
echo "========== AVITI: ALL VARIANTS =========="

ALL_AV_OUTPUT="${WORK_DIR}/stratified_validation_ALL_aviti.csv"
echo "Sample,Pipeline,Platform,Region,Total_Variants,Validated_Variants,Validation_Rate" > "$ALL_AV_OUTPUT"

for SAMPLE in "${AVITI_SAMPLES[@]}"; do
    DIR_PATH="${AVITI_VCF_BASE}/${SAMPLE}"
    if [[ ! -d "$DIR_PATH" ]]; then
        echo "  SKIP: $SAMPLE (directory not found)"
        continue
    fi

    echo "  ${SAMPLE}..."

    for PIPELINE in gatk ppr; do
        PIPELINE_UPPER=$(echo "$PIPELINE" | tr '[:lower:]' '[:upper:]')

        # Aviti GATK files have "_chr" suffix, PPR does not
        if [[ "$PIPELINE" == "gatk" ]]; then
            ORIGINAL="${DIR_PATH}/${SAMPLE}_gatk_all_chr.vcf.gz"
        else
            ORIGINAL="${DIR_PATH}/${SAMPLE}_ppr_all.vcf.gz"
        fi
        VALIDATED_VCF="${DIR_PATH}/${SAMPLE}_${PIPELINE}_all_validated.filter.vcf.gz"

        if [[ ! -f "$ORIGINAL" ]]; then
            echo "    WARNING: Missing original: $ORIGINAL"
            continue
        fi
        if [[ ! -f "$VALIDATED_VCF" ]]; then
            echo "    WARNING: Missing validated: $VALIDATED_VCF"
            continue
        fi

        for REGION in "${!REGION_BEDS[@]}"; do
            BED_FILE="${REGION_BEDS[$REGION]}"
            [[ ! -f "$BED_FILE" ]] && continue

            TOTAL=$(count_variants "$ORIGINAL" "$BED_FILE")
            VALIDATED=$(count_variants "$VALIDATED_VCF" "$BED_FILE")

            if [[ "$TOTAL" -gt 0 ]]; then
                RATE=$(python3 -c "print(round(${VALIDATED}/${TOTAL}*100, 4))")
            else
                RATE="NA"
            fi

            echo "${SAMPLE},${PIPELINE_UPPER},Aviti,${REGION},${TOTAL},${VALIDATED},${RATE}" >> "$ALL_AV_OUTPUT"
        done
    done
done

echo "  ✓ Saved: $ALL_AV_OUTPUT"
echo ""

# ===========================================================================
# PART 3: UNIQUE VARIANTS — ILLUMINA (corrected with bcftools isec)
# ===========================================================================
echo "========== ILLUMINA: UNIQUE VARIANTS (bcftools isec) =========="

UNIQUE_ILL_OUTPUT="${WORK_DIR}/stratified_validation_UNIQUE_illumina.csv"
echo "Sample,Pipeline,Platform,Region,Total_Variants,Validated_Variants,Validation_Rate" > "$UNIQUE_ILL_OUTPUT"

for SAMPLE_DIR in "${ILLUMINA_SAMPLE_DIRS[@]}"; do
    DIR_PATH="${ILLUMINA_VCF_BASE}/${SAMPLE_DIR}"
    if [[ ! -d "$DIR_PATH" ]]; then
        continue
    fi

    SAMPLE_ID=$(illumina_dir_to_icdg "$SAMPLE_DIR")
    echo "  ${SAMPLE_DIR} (${SAMPLE_ID})..."

    for PIPELINE in gatk ppr; do
        PIPELINE_UPPER=$(echo "$PIPELINE" | tr '[:lower:]' '[:upper:]')

        # Original unique VCFs are in /uniquevariants/
        ORIGINAL_UNIQUE="${UNIQUE_BASE}/unique${PIPELINE_UPPER}illumina/${SAMPLE_ID}_unique_${PIPELINE_UPPER}_illumina.vcf.gz"

        # Merfin-validated COMPLETE VCF (all variants that passed Merfin filter)
        VALIDATED_ALL="${DIR_PATH}/${SAMPLE_DIR}_${PIPELINE}_all_validated.filter.vcf.gz"

        if [[ ! -f "$ORIGINAL_UNIQUE" ]]; then
            echo "    WARNING: Missing unique original: $ORIGINAL_UNIQUE"
            continue
        fi
        if [[ ! -f "$VALIDATED_ALL" ]]; then
            echo "    WARNING: Missing validated ALL VCF: $VALIDATED_ALL"
            continue
        fi

        for REGION in "${!REGION_BEDS[@]}"; do
            BED_FILE="${REGION_BEDS[$REGION]}"
            [[ ! -f "$BED_FILE" ]] && continue

            TMP_ISEC="${TMPDIR_BASE}/${SAMPLE_ID}_${PIPELINE}_illumina_${REGION}"

            # Total unique variants in region
            TOTAL=$(count_variants "$ORIGINAL_UNIQUE" "$BED_FILE")

            # Validated = unique variants present in the Merfin-validated complete set
            VALIDATED=$(count_validated_unique_in_region "$ORIGINAL_UNIQUE" "$VALIDATED_ALL" "$BED_FILE" "$TMP_ISEC")

            if [[ "$TOTAL" -gt 0 ]]; then
                RATE=$(python3 -c "print(round(${VALIDATED}/${TOTAL}*100, 4))")
            else
                RATE="NA"
            fi

            echo "${SAMPLE_ID},${PIPELINE_UPPER},Illumina,${REGION},${TOTAL},${VALIDATED},${RATE}" >> "$UNIQUE_ILL_OUTPUT"
        done
    done
done

echo "  ✓ Saved: $UNIQUE_ILL_OUTPUT"
echo ""

# ===========================================================================
# PART 4: UNIQUE VARIANTS — AVITI (corrected with bcftools isec)
# ===========================================================================
echo "========== AVITI: UNIQUE VARIANTS (bcftools isec) =========="

UNIQUE_AV_OUTPUT="${WORK_DIR}/stratified_validation_UNIQUE_aviti.csv"
echo "Sample,Pipeline,Platform,Region,Total_Variants,Validated_Variants,Validation_Rate" > "$UNIQUE_AV_OUTPUT"

for SAMPLE in "${AVITI_SAMPLES[@]}"; do
    DIR_PATH="${AVITI_VCF_BASE}/${SAMPLE}"
    if [[ ! -d "$DIR_PATH" ]]; then
        continue
    fi

    echo "  ${SAMPLE}..."

    for PIPELINE in gatk ppr; do
        PIPELINE_UPPER=$(echo "$PIPELINE" | tr '[:lower:]' '[:upper:]')

        # Original unique VCFs
        ORIGINAL_UNIQUE="${UNIQUE_BASE}/unique${PIPELINE_UPPER}aviti/${SAMPLE}_unique_${PIPELINE_UPPER}_aviti.vcf.gz"

        # Merfin-validated COMPLETE VCF
        VALIDATED_ALL="${DIR_PATH}/${SAMPLE}_${PIPELINE}_all_validated.filter.vcf.gz"

        if [[ ! -f "$ORIGINAL_UNIQUE" ]]; then
            echo "    WARNING: Missing unique original: $ORIGINAL_UNIQUE"
            continue
        fi
        if [[ ! -f "$VALIDATED_ALL" ]]; then
            echo "    WARNING: Missing validated ALL VCF: $VALIDATED_ALL"
            continue
        fi

        for REGION in "${!REGION_BEDS[@]}"; do
            BED_FILE="${REGION_BEDS[$REGION]}"
            [[ ! -f "$BED_FILE" ]] && continue

            TMP_ISEC="${TMPDIR_BASE}/${SAMPLE}_${PIPELINE}_aviti_${REGION}"

            TOTAL=$(count_variants "$ORIGINAL_UNIQUE" "$BED_FILE")
            VALIDATED=$(count_validated_unique_in_region "$ORIGINAL_UNIQUE" "$VALIDATED_ALL" "$BED_FILE" "$TMP_ISEC")

            if [[ "$TOTAL" -gt 0 ]]; then
                RATE=$(python3 -c "print(round(${VALIDATED}/${TOTAL}*100, 4))")
            else
                RATE="NA"
            fi

            echo "${SAMPLE},${PIPELINE_UPPER},Aviti,${REGION},${TOTAL},${VALIDATED},${RATE}" >> "$UNIQUE_AV_OUTPUT"
        done
    done
done

echo "  ✓ Saved: $UNIQUE_AV_OUTPUT"
echo ""

# ===========================================================================
# MERGE ALL RESULTS & SUMMARY
# ===========================================================================
echo "========== MERGING RESULTS =========="

COMBINED="${WORK_DIR}/stratified_validation_ALL_combined.csv"
echo "Sample,Pipeline,Platform,Region,Total_Variants,Validated_Variants,Validation_Rate" > "$COMBINED"
tail -n +2 "$ALL_ILL_OUTPUT" >> "$COMBINED"
tail -n +2 "$ALL_AV_OUTPUT" >> "$COMBINED"
echo "  ✓ Combined ALL variants: $COMBINED"

COMBINED_UNIQUE="${WORK_DIR}/stratified_validation_UNIQUE_combined.csv"
echo "Sample,Pipeline,Platform,Region,Total_Variants,Validated_Variants,Validation_Rate" > "$COMBINED_UNIQUE"
tail -n +2 "$UNIQUE_ILL_OUTPUT" >> "$COMBINED_UNIQUE"
tail -n +2 "$UNIQUE_AV_OUTPUT" >> "$COMBINED_UNIQUE"
echo "  ✓ Combined UNIQUE variants: $COMBINED_UNIQUE"

# Cleanup temp directory
rm -rf "$TMPDIR_BASE"

echo ""
echo "========== SUMMARY =========="
echo ""

python3 << 'PYEOF'
import pandas as pd
import sys

work_dir = "/mnt/genomics/pilot_PPR/kmers_validation/validationrates_stratified"

for variant_set, filename in [("ALL", "stratified_validation_ALL_combined.csv"),
                               ("UNIQUE", "stratified_validation_UNIQUE_combined.csv")]:
    csv_path = f"{work_dir}/{filename}"
    try:
        df = pd.read_csv(csv_path)
    except Exception:
        continue

    df = df[df['Validation_Rate'] != 'NA'].copy()
    if df.empty:
        print(f"\n  {variant_set}: No data available yet.")
        continue

    df['Validation_Rate'] = df['Validation_Rate'].astype(float)

    print(f"\n{'='*70}")
    print(f" {variant_set} VARIANTS — Mean Validation Rate (%) by Pipeline × Region")
    print(f"{'='*70}")

    # Per platform
    for platform in sorted(df['Platform'].unique()):
        pdata = df[df['Platform'] == platform]
        pivot = pdata.pivot_table(
            index='Region',
            columns='Pipeline',
            values='Validation_Rate',
            aggfunc='mean'
        ).round(2)

        if not pivot.empty:
            if 'PPR' in pivot.columns and 'GATK' in pivot.columns:
                pivot['PPR-GATK'] = (pivot['PPR'] - pivot['GATK']).round(2)
            pivot = pivot.sort_index()
            print(f"\n  Platform: {platform} (n={len(pdata['Sample'].unique())} samples)")
            print(f"  {'-'*60}")
            print(pivot.to_string())

    # Sanity check for UNIQUE
    if variant_set == "UNIQUE":
        over100 = df[df['Validation_Rate'] > 100.0]
        if len(over100) > 0:
            print(f"\n  ⚠️  WARNING: {len(over100)} entries have validation rate > 100%!")
            print(f"     Check bcftools isec logic or VCF normalization.")
        else:
            print(f"\n  ✓ All UNIQUE validation rates are ≤ 100% — logic correct.")

print()
PYEOF

echo ""
echo "================================================================"
echo " DONE! Results in: ${WORK_DIR}/"
echo "================================================================"
echo "Files:"
echo "  stratified_validation_ALL_illumina.csv"
echo "  stratified_validation_ALL_aviti.csv"
echo "  stratified_validation_ALL_combined.csv"
echo "  stratified_validation_UNIQUE_illumina.csv"
echo "  stratified_validation_UNIQUE_aviti.csv"
echo "  stratified_validation_UNIQUE_combined.csv"
echo ""
echo "Use these CSVs for downstream plotting."
