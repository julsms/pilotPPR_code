#!/usr/bin/env Rscript

# ============================================================================
# Script 1: Coverage vs Merfin Validation Rates (CORRECTED)
#
# Investigates whether sequencing coverage (per sample) correlates with
# Merfin k-mer validation rates for:
#   - ALL variants (from summary files in final_results/)
#   - UNIQUE variants (from unique_validation_rates_corrected.csv)
#
# Input:
#   - Coverage: /mnt/genomics/pilot_PPR/coverage/qc_metrics_collected_for_plots.tsv
#   - Merfin ALL: /mnt/genomics/pilot_PPR/kmers_validation/{illumina,aviti}/final_results/summary_*.txt
#   - Merfin UNIQUE: /mnt/genomics/pilot_PPR/kmers_validation/{illumina,aviti}/uniquevariants/unique_validation_rates_corrected.csv
#
# Output:
#   - correlation_coverage_validation_results.csv (all correlation stats)
#   - correlation_coverage_validation_data.csv (merged data for plotting)
#
# Author: Generated for manuscript analysis
# Date: 2026-06-29
# ============================================================================

suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
})

cat("============================================================\n")
cat("Coverage vs Merfin Validation Rates — Correlation Analysis\n")
cat("============================================================\n\n")

# ----------------------------------------------------------------------------
# 1. LOAD COVERAGE DATA
# ----------------------------------------------------------------------------

coverage_file <- "/mnt/genomics/pilot_PPR/coverage_correlations/qc_metrics_collected_for_plots.tsv"
cov_df <- read.delim(coverage_file, header = TRUE, stringsAsFactors = FALSE)

cat("Coverage data loaded:", nrow(cov_df), "rows\n")

# Filter to only Illumina and AVITI (exclude ONT)
cov_df <- cov_df[cov_df$platform %in% c("Illumina", "AVITI"), ]
cat("  After filtering to Illumina/AVITI:", nrow(cov_df), "rows\n")

# Standardize sample IDs to numeric ICDG format
# Illumina: "10ICDG_S5" -> 10, "1ICDG_S1" -> 1, etc.
# Aviti: "ICDG1" -> 1, "ICDG10" -> 10, etc.
cov_df$sample_num <- NA

for (i in 1:nrow(cov_df)) {
    s <- cov_df$sample[i]
    if (cov_df$platform[i] == "AVITI") {
        cov_df$sample_num[i] <- as.integer(sub("ICDG", "", s))
    } else {
        # Illumina: 10ICDG_S5, 1ICDG_S1, etc.
        cov_df$sample_num[i] <- as.integer(sub("ICDG_S.*", "", s))
    }
}

cat("  Illumina samples:", sum(cov_df$platform == "Illumina"), "\n")
cat("  Aviti samples:", sum(cov_df$platform == "AVITI"), "\n")

# ----------------------------------------------------------------------------
# 2. PARSE MERFIN ALL VARIANTS FROM SUMMARY FILES
# ----------------------------------------------------------------------------

parse_merfin_summary <- function(filepath) {
    lines <- readLines(filepath, warn = FALSE)
    result <- list()
    
    for (line in lines) {
        if (grepl("GATK All Variants:", line)) {
            pct <- as.numeric(sub(".*\\((.*)%\\).*", "\\1", line))
            result$gatk_all <- pct
        } else if (grepl("PPR All Variants:", line)) {
            pct <- as.numeric(sub(".*\\((.*)%\\).*", "\\1", line))
            result$ppr_all <- pct
        }
    }
    return(result)
}

# Parse Illumina ALL summaries
illumina_dir <- "/mnt/genomics/pilot_PPR/kmers_validation/illumina/final_results"
aviti_dir <- "/mnt/genomics/pilot_PPR/kmers_validation/aviti/final_results"

merfin_all_data <- data.frame()

# Illumina
illumina_files <- list.files(illumina_dir, pattern = "^summary_illumina_.*\\.txt$", full.names = TRUE)
cat("\nParsing", length(illumina_files), "Illumina Merfin ALL summaries...\n")

for (f in illumina_files) {
    sample_name <- sub("summary_illumina_(.*)\\.txt", "\\1", basename(f))
    sample_num <- as.integer(sub("ICDG_S.*", "", sample_name))
    
    parsed <- parse_merfin_summary(f)
    
    if (!is.null(parsed$gatk_all) && !is.null(parsed$ppr_all)) {
        merfin_all_data <- rbind(merfin_all_data, data.frame(
            platform = "Illumina",
            sample = sample_name,
            sample_num = sample_num,
            gatk_all_validation = parsed$gatk_all,
            ppr_all_validation = parsed$ppr_all,
            stringsAsFactors = FALSE
        ))
    }
}

# Aviti
aviti_files <- list.files(aviti_dir, pattern = "^summary_aviti_.*\\.txt$", full.names = TRUE)
cat("Parsing", length(aviti_files), "Aviti Merfin ALL summaries...\n")

for (f in aviti_files) {
    sample_name <- sub("summary_aviti_(.*)\\.txt", "\\1", basename(f))
    sample_num <- as.integer(sub("ICDG", "", sample_name))
    
    parsed <- parse_merfin_summary(f)
    
    if (!is.null(parsed$gatk_all) && !is.null(parsed$ppr_all)) {
        merfin_all_data <- rbind(merfin_all_data, data.frame(
            platform = "AVITI",
            sample = sample_name,
            sample_num = sample_num,
            gatk_all_validation = parsed$gatk_all,
            ppr_all_validation = parsed$ppr_all,
            stringsAsFactors = FALSE
        ))
    }
}

cat("  Total ALL validation records:", nrow(merfin_all_data), "\n")

# ----------------------------------------------------------------------------
# 3. LOAD UNIQUE VARIANTS VALIDATION FROM CORRECTED CSV FILES
# ----------------------------------------------------------------------------

cat("\nLoading UNIQUE variant validation rates from corrected CSVs...\n")

illumina_unique_file <- "/mnt/genomics/pilot_PPR/kmers_validation/illumina/uniquevariants/unique_validation_rates_corrected.csv"
aviti_unique_file <- "/mnt/genomics/pilot_PPR/kmers_validation/aviti/uniquevariants/unique_validation_rates_corrected.csv"

# Read Illumina unique
illumina_unique <- read.csv(illumina_unique_file, stringsAsFactors = FALSE)
illumina_unique$platform <- "Illumina"
# Extract sample_num from Sample column (e.g., "10ICDG_S5" -> 10)
illumina_unique$sample_num <- as.integer(sub("ICDG_S.*", "", illumina_unique$Sample))

# Read Aviti unique
aviti_unique <- read.csv(aviti_unique_file, stringsAsFactors = FALSE)
aviti_unique$platform <- "AVITI"
# Extract sample_num from Sample column (e.g., "ICDG10" -> 10)
aviti_unique$sample_num <- as.integer(sub("ICDG", "", aviti_unique$Sample))

# Combine
unique_all <- rbind(illumina_unique, aviti_unique)

cat("  Illumina unique records:", nrow(illumina_unique), "\n")
cat("  Aviti unique records:", nrow(aviti_unique), "\n")

# Pivot to get GATK and PPR unique validation per sample
unique_wide <- unique_all %>%
    select(platform, sample_num, Pipeline, Validation_Rate) %>%
    pivot_wider(
        names_from = Pipeline,
        values_from = Validation_Rate,
        names_prefix = "unique_validation_"
    ) %>%
    rename(
        gatk_unique_validation = unique_validation_GATK,
        ppr_unique_validation = unique_validation_PPR
    )

cat("  Unique validation wide:", nrow(unique_wide), "samples\n")

# ----------------------------------------------------------------------------
# 4. MERGE EVERYTHING
# ----------------------------------------------------------------------------

# First merge coverage with ALL validation
merged_df <- merge(
    cov_df[, c("platform", "sample_num", "mean_coverage", "coverage_cv_50kb")],
    merfin_all_data[, c("platform", "sample_num", "gatk_all_validation", "ppr_all_validation")],
    by = c("platform", "sample_num"),
    all.x = FALSE
)

# Then merge with UNIQUE validation
merged_df <- merge(
    merged_df,
    unique_wide[, c("platform", "sample_num", "gatk_unique_validation", "ppr_unique_validation")],
    by = c("platform", "sample_num"),
    all.x = TRUE  # keep ALL rows even if no unique validation
)

# Standardize platform names for output
merged_df$platform[merged_df$platform == "AVITI"] <- "Aviti"

cat("\nMerged data:", nrow(merged_df), "rows\n")
cat("  Illumina:", sum(merged_df$platform == "Illumina"), "\n")
cat("  Aviti:", sum(merged_df$platform == "Aviti"), "\n")
cat("  With unique validation:", sum(!is.na(merged_df$gatk_unique_validation)), "\n")

# Exclude ICDG23 (very low coverage = 1.23x, outlier)
cat("\nExcluding ICDG23 (coverage = 1.23x, outlier)...\n")
merged_df <- merged_df[merged_df$sample_num != 23, ]
cat("  After exclusion:", nrow(merged_df), "rows\n")

# ----------------------------------------------------------------------------
# 5. CORRELATION ANALYSIS
# ----------------------------------------------------------------------------

cat("\n============================================================\n")
cat("CORRELATION RESULTS\n")
cat("============================================================\n\n")

calc_correlation <- function(x, y) {
    complete <- complete.cases(x, y)
    x <- x[complete]
    y <- y[complete]
    n <- length(x)
    
    if (n < 5) return(data.frame(n = n, pearson_r = NA, pearson_p = NA,
                                  spearman_rho = NA, spearman_p = NA))
    
    pearson <- cor.test(x, y, method = "pearson")
    spearman <- cor.test(x, y, method = "spearman", exact = FALSE)
    
    data.frame(
        n = n,
        pearson_r = round(pearson$estimate, 4),
        pearson_p = pearson$p.value,
        spearman_rho = round(spearman$estimate, 4),
        spearman_p = spearman$p.value
    )
}

# Define comparisons
comparisons <- list(
    list(name = "GATK_ALL", col = "gatk_all_validation"),
    list(name = "PPR_ALL", col = "ppr_all_validation"),
    list(name = "GATK_UNIQUE", col = "gatk_unique_validation"),
    list(name = "PPR_UNIQUE", col = "ppr_unique_validation")
)

platforms <- c("Illumina", "Aviti", "Both")
results <- data.frame()

cat("--- Mean Coverage vs Validation ---\n\n")

for (comp in comparisons) {
    for (plat in platforms) {
        if (plat == "Both") {
            subset_df <- merged_df
        } else {
            subset_df <- merged_df[merged_df$platform == plat, ]
        }
        
        corr <- calc_correlation(subset_df$mean_coverage, subset_df[[comp$col]])
        
        result_row <- data.frame(
            comparison = comp$name,
            platform = plat,
            n = corr$n,
            pearson_r = corr$pearson_r,
            pearson_p = corr$pearson_p,
            spearman_rho = corr$spearman_rho,
            spearman_p = corr$spearman_p,
            stringsAsFactors = FALSE
        )
        results <- rbind(results, result_row)
        
        sig_p <- ifelse(is.na(corr$pearson_p), "NA",
                  ifelse(corr$pearson_p < 0.001, "***",
                   ifelse(corr$pearson_p < 0.01, "**",
                    ifelse(corr$pearson_p < 0.05, "*", "ns"))))
        sig_s <- ifelse(is.na(corr$spearman_p), "NA",
                  ifelse(corr$spearman_p < 0.001, "***",
                   ifelse(corr$spearman_p < 0.01, "**",
                    ifelse(corr$spearman_p < 0.05, "*", "ns"))))
        
        cat(sprintf("  %-15s | %-10s | n=%2d | Pearson r=%.3f (%s) | Spearman rho=%.3f (%s)\n",
                    comp$name, plat, corr$n, corr$pearson_r, sig_p, corr$spearman_rho, sig_s))
    }
    cat("\n")
}

# Also test coverage CV vs validation
cat("\n--- Coverage CV (uniformity) vs Validation ---\n\n")

for (comp in comparisons) {
    for (plat in platforms) {
        if (plat == "Both") {
            subset_df <- merged_df
        } else {
            subset_df <- merged_df[merged_df$platform == plat, ]
        }
        
        corr <- calc_correlation(subset_df$coverage_cv_50kb, subset_df[[comp$col]])
        
        result_row <- data.frame(
            comparison = paste0(comp$name, "_vs_CV"),
            platform = plat,
            n = corr$n,
            pearson_r = corr$pearson_r,
            pearson_p = corr$pearson_p,
            spearman_rho = corr$spearman_rho,
            spearman_p = corr$spearman_p,
            stringsAsFactors = FALSE
        )
        results <- rbind(results, result_row)
        
        sig_p <- ifelse(is.na(corr$pearson_p), "NA",
                  ifelse(corr$pearson_p < 0.001, "***",
                   ifelse(corr$pearson_p < 0.01, "**",
                    ifelse(corr$pearson_p < 0.05, "*", "ns"))))
        
        cat(sprintf("  %-20s | %-10s | n=%2d | Pearson r=%.3f (%s) | Spearman rho=%.3f\n",
                    paste0(comp$name, "_vs_CV"), plat, corr$n, corr$pearson_r, sig_p, corr$spearman_rho))
    }
    cat("\n")
}

# ----------------------------------------------------------------------------
# 6. SAVE RESULTS
# ----------------------------------------------------------------------------

output_dir <- "/mnt/genomics/pilot_PPR/coverage_correlations"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(results, file.path(output_dir, "correlation_coverage_validation_results.csv"),
          row.names = FALSE)
write.csv(merged_df, file.path(output_dir, "correlation_coverage_validation_data.csv"),
          row.names = FALSE)

cat("\n============================================================\n")
cat("Results saved to:", output_dir, "\n")
cat("  - correlation_coverage_validation_results.csv\n")
cat("  - correlation_coverage_validation_data.csv\n")
cat("============================================================\n")
