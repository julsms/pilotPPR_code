#!/usr/bin/env Rscript

## =============================================================================
## GENERATE concordance_summary_for_analysis_corrected.tsv
##
## Removes duplicate/invalid entries from the concordance summary:
##   - ICDG12: ran twice, keep ONLY second occurrence (correct Aviti data)
##   - ICDG21: ran twice, keep ONLY second occurrence (correct Aviti data)
##   - ICDG22: identical duplicates, keep ONLY one copy
##   - ICDG1, ICDG2, ICDG3: OMIM has NA rows then valid rows, keep ONLY valid
##
## Input:  concordance_summary_for_analysis.tsv
## Output: concordance_summary_for_analysis_corrected.tsv
## =============================================================================

library(dplyr)

## --- Configuration ---
## Set this to your actual base directory:
# base_dir <- "/mnt/genomics/pilot_PPR/concordance_perbed"
# setwd(base_dir)

## For testing, use current directory:
input_file <- "concordance_summary_for_analysis.tsv"
output_file <- "concordance_summary_for_analysis_corrected.tsv"

cat("=================================================================\n")
cat("GENERATING CORRECTED concordance_summary_for_analysis.tsv\n")
cat("=================================================================\n\n")

## Read the original data
data <- read.table(input_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

cat("Original data:\n")
cat("  Total records:", nrow(data), "\n")
cat("  Unique samples:", length(unique(data$Sample)), "\n")
cat("  Unique BED regions:", length(unique(data$BED_Region)), "\n\n")

## --- Step 1: Identify duplicates ---
## For each Sample × BED_Region × Variant_Type × Comparison, count occurrences
dup_check <- data %>%
    group_by(Sample, BED_Region, Variant_Type, Comparison) %>%
    mutate(occurrence = row_number(), total_occurrences = n()) %>%
    ungroup()

duplicated_entries <- dup_check %>% filter(total_occurrences > 1)
cat("Duplicated entries found:\n")
cat("  Total duplicate rows:", nrow(duplicated_entries), "\n")
cat("  Affected samples:", paste(sort(unique(duplicated_entries$Sample)), collapse = ", "), "\n")
cat("  Affected regions:", paste(sort(unique(duplicated_entries$BED_Region)), collapse = ", "), "\n\n")

## --- Step 2: Apply corrections ---
##
## A) ICDG12 and ICDG21: These were run twice with different Aviti data.
##    The FIRST occurrence has bad/low-coverage Aviti data, the SECOND is correct.
##    -> Keep ONLY the LAST (second) occurrence for each duplicated combination.
##
## B) ICDG22: Identical duplicates (same values repeated).
##    -> Keep ONLY the FIRST occurrence (or last, doesn't matter - they're identical).
##
## C) ICDG1, ICDG2, ICDG3 for OMIM: First set has NA values (N_Total_Union = 0),
##    second set has valid data.
##    -> Remove rows where ALL count columns are 0 (and concordance is NA).
##    -> Then deduplicate keeping the valid (non-NA) entry.
##
## UNIFIED APPROACH:
## 1. First, remove rows that are clearly invalid (all zeros / NA concordance)
## 2. Then, for remaining duplicates, keep the LAST occurrence
##    (this handles ICDG12/ICDG21 where second run is correct,
##     and ICDG22 where entries are identical so order doesn't matter)

cat("--- Applying corrections ---\n\n")

## Step 2a: Remove completely empty/invalid rows (OMIM NA issue for ICDG1/2/3)
## These are rows where N_Total_Union == 0 and concordance is NA
invalid_rows <- data %>%
    filter(N_Total_Union == 0 & is.na(Pos_Conc_Pct) & is.na(GT_Conc_Pct))

cat("Step 2a: Removing invalid (all-zero) rows:\n")
cat("  Found", nrow(invalid_rows), "invalid rows\n")
if (nrow(invalid_rows) > 0) {
    cat("  Affected samples:", paste(sort(unique(invalid_rows$Sample)), collapse = ", "), "\n")
    cat("  Affected regions:", paste(sort(unique(invalid_rows$BED_Region)), collapse = ", "), "\n")
}
cat("\n")

data_clean <- data %>%
    filter(!(N_Total_Union == 0 & is.na(Pos_Conc_Pct) & is.na(GT_Conc_Pct)))

cat("  Records after removing invalid rows:", nrow(data_clean), "\n\n")

## Step 2b: For remaining duplicates, keep ONLY the LAST occurrence
## This correctly handles:
##   - ICDG12: first run (bad Aviti) removed, second run (correct) kept
##   - ICDG21: first run (bad Aviti) removed, second run (correct) kept
##   - ICDG22: identical duplicates, keeping last (same as first)

data_corrected <- data_clean %>%
    group_by(Sample, BED_Region, Variant_Type, Comparison) %>%
    slice_tail(n = 1) %>%
    ungroup()

n_removed_dedup <- nrow(data_clean) - nrow(data_corrected)
cat("Step 2b: Deduplicating (keeping last occurrence for each combination):\n")
cat("  Removed", n_removed_dedup, "duplicate rows\n")
cat("  Records after deduplication:", nrow(data_corrected), "\n\n")

## --- Step 3: Validation ---
cat("--- Validation ---\n\n")

## Check expected record count: 23 samples × 9 regions × 4 comparisons × 3 variant types = 2484
expected_records <- 23 * 9 * 4 * 3
cat("Expected records: 23 samples × 9 regions × 4 comparisons × 3 variant types =", expected_records, "\n")
cat("Actual records:", nrow(data_corrected), "\n")

if (nrow(data_corrected) == expected_records) {
    cat("Match: data is complete and deduplicated.\n\n")
} else if (nrow(data_corrected) < expected_records) {
    cat("WARNING: fewer records than expected. Some sample-region combinations may be missing.\n\n")
} else {
    cat("WARNING: more records than expected. Deduplication may not be complete.\n\n")
}

## Per-sample validation
sample_counts <- data_corrected %>%
    group_by(Sample) %>%
    summarise(N = n(), .groups = "drop") %>%
    arrange(Sample)

expected_per_sample <- 9 * 4 * 3  # 108 records per sample
cat("Per-sample record counts (expected:", expected_per_sample, "per sample):\n")
for (i in 1:nrow(sample_counts)) {
    flag <- ifelse(sample_counts$N[i] != expected_per_sample, " *** ISSUE", " ok")
    cat(sprintf("  %s: %d%s\n", sample_counts$Sample[i], sample_counts$N[i], flag))
}

## Check no more duplicates exist
remaining_dups <- data_corrected %>%
    group_by(Sample, BED_Region, Variant_Type, Comparison) %>%
    filter(n() > 1) %>%
    ungroup()

cat("\nRemaining duplicates after correction:", nrow(remaining_dups), "\n")
if (nrow(remaining_dups) > 0) {
    cat("WARNING: still have duplicates! Manual inspection needed.\n")
    print(remaining_dups %>% select(Sample, BED_Region, Variant_Type, Comparison) %>% distinct())
}

## --- Step 4: Show what changed for key samples ---
cat("\n--- Impact summary for affected samples ---\n\n")

for (s in c("ICDG12", "ICDG21")) {
    cat(sprintf("  %s (15GC85, all variants, IlluminaPPR_vs_AvitiPPR):\n", s))

    original_vals <- data %>%
        filter(Sample == s, BED_Region == "15GC85",
               Variant_Type == "all", Comparison == "IlluminaPPR_vs_AvitiPPR") %>%
        pull(GT_Conc_Pct)

    corrected_val <- data_corrected %>%
        filter(Sample == s, BED_Region == "15GC85",
               Variant_Type == "all", Comparison == "IlluminaPPR_vs_AvitiPPR") %>%
        pull(GT_Conc_Pct)

    cat(sprintf("    Original values: %s\n", paste(original_vals, collapse = ", ")))
    cat(sprintf("    Corrected value (kept): %s\n\n", corrected_val))
}

## --- Step 5: Write output ---
write.table(data_corrected, output_file, sep = "\t", quote = FALSE, row.names = FALSE)

cat("=================================================================\n")
cat("OUTPUT WRITTEN:", output_file, "\n")
cat(sprintf("Original records: %d -> Corrected records: %d (removed %d)\n",
            nrow(data), nrow(data_corrected), nrow(data) - nrow(data_corrected)))
cat("=================================================================\n")
cat("\nNext step: run wilcoxon_perregion_analysis.R, which reads this\n")
cat("corrected file directly.\n")
