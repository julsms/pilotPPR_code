#!/usr/bin/env Rscript

# ============================================================================
# Script 2: Coverage vs Concordance (UPDATED — uses new concordance directory)
#
# Correlates sequencing coverage with genotype concordance.
# For each concordance comparison, we test SEPARATELY:
#   - Illumina coverage vs concordance
#   - Aviti coverage vs concordance
# (No averaging of coverage between platforms.)
#
# Concordance pairs (6 total from the 4x4 matrix):
#   Cross-platform, same method:
#     1. Illumina PPR vs Aviti PPR
#     2. Illumina GATK vs Aviti GATK
#   Cross-method, same platform:
#     3. Illumina PPR vs Illumina GATK
#     4. Aviti PPR vs Aviti GATK
#   Cross-platform, cross-method:
#     5. Illumina PPR vs Aviti GATK
#     6. Illumina GATK vs Aviti PPR
#
# SWAP: Samples ICDG12 and ICDG21 have a known sample swap between
# Illumina and Aviti. The concordance analysis already applied this swap.
# For ICDG12: concordance used Illumina_12 + Aviti_21
# For ICDG21: concordance used Illumina_21 + Aviti_12
# We swap the Aviti coverage accordingly when correlating.
#
# CHANGE FROM PREVIOUS VERSION:
#   - Concordance path: concordance_IlluminaAvitiGATKPPR_new (was: concordance_IlluminaAvitiGATKPPR)
#   - Output path: coverage_correlations_again (was: coverage_correlations)
#
# Input:
#   - Coverage: /mnt/genomics/pilot_PPR/coverage_correlations/qc_metrics_collected_for_plots.tsv
#   - Concordance: /mnt/genomics/pilot_PPR/concordance_IlluminaAvitiGATKPPR_new/sample_ICDG*/
#
# Output:
#   - correlation_coverage_concordance_results.csv
#   - correlation_coverage_concordance_data.csv
# ============================================================================

suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
})

cat("============================================================\n")
cat("Coverage vs Concordance — Correlation Analysis\n")
cat("Using NEW concordance directory\n")
cat("============================================================\n\n")

# ----------------------------------------------------------------------------
# 1. LOAD COVERAGE DATA
# ----------------------------------------------------------------------------

coverage_file <- "/mnt/genomics/pilot_PPR/coverage_correlations/qc_metrics_collected_for_plots.tsv"
cov_df <- read.delim(coverage_file, header = TRUE, stringsAsFactors = FALSE)

# Filter to only Illumina and AVITI (exclude ONT)
cov_df <- cov_df[cov_df$platform %in% c("Illumina", "AVITI"), ]

# Standardize sample IDs to numeric
cov_df$sample_num <- NA
for (i in 1:nrow(cov_df)) {
    s <- cov_df$sample[i]
    if (cov_df$platform[i] == "AVITI") {
        cov_df$sample_num[i] <- as.integer(sub("ICDG", "", s))
    } else {
        cov_df$sample_num[i] <- as.integer(sub("ICDG_S.*", "", s))
    }
}

# Standardize platform name
cov_df$platform[cov_df$platform == "AVITI"] <- "Aviti"

# Pivot to wide format — one row per sample with Illumina and Aviti coverage
cov_wide <- cov_df %>%
    select(sample_num, platform, mean_coverage) %>%
    group_by(sample_num, platform) %>%
    summarise(mean_coverage = mean(mean_coverage, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = platform, values_from = mean_coverage,
                names_prefix = "cov_")

cat("Coverage data prepared:", nrow(cov_wide), "samples\n")
cat("  Columns:", paste(names(cov_wide), collapse = ", "), "\n")

# Apply swap for Aviti coverage:
# sample_ICDG12 concordance used Aviti ICDG21's data -> use ICDG21 Aviti coverage
# sample_ICDG21 concordance used Aviti ICDG12's data -> use ICDG12 Aviti coverage
cov_wide$cov_Aviti_swapped <- cov_wide$cov_Aviti

idx_12 <- which(cov_wide$sample_num == 12)
idx_21 <- which(cov_wide$sample_num == 21)

if (length(idx_12) == 1 && length(idx_21) == 1) {
    cov_wide$cov_Aviti_swapped[idx_12] <- cov_wide$cov_Aviti[idx_21]
    cov_wide$cov_Aviti_swapped[idx_21] <- cov_wide$cov_Aviti[idx_12]
    cat("Swap applied: ICDG12 <-> ICDG21 Aviti coverage\n")
}

# ----------------------------------------------------------------------------
# 2. PARSE CONCORDANCE MATRICES
# ----------------------------------------------------------------------------

read_concordance_matrix <- function(filepath) {
    lines <- readLines(filepath, warn = FALSE)
    
    # Remove empty lines
    lines <- lines[nchar(trimws(lines)) > 0]
    
    n_methods <- 4
    method_names <- c("Illumina_PPR", "Illumina_GATK", "Aviti_PPR", "Aviti_GATK")
    
    # Extract all numeric values — each data row has exactly 4 numbers
    data_rows <- list()
    
    for (line in lines) {
        parts <- strsplit(line, "\t")[[1]]
        nums <- suppressWarnings(as.numeric(parts))
        valid_nums <- nums[!is.na(nums)]
        
        if (length(valid_nums) == 4) {
            data_rows[[length(data_rows) + 1]] <- valid_nums
        }
    }
    
    # If not found cleanly, try merging consecutive lines
    if (length(data_rows) != n_methods) {
        data_rows <- list()
        i <- 1
        while (i <= length(lines)) {
            parts <- strsplit(lines[i], "\t")[[1]]
            nums <- suppressWarnings(as.numeric(parts))
            valid_nums <- nums[!is.na(nums)]
            
            if (length(valid_nums) == 4) {
                data_rows[[length(data_rows) + 1]] <- valid_nums
                i <- i + 1
            } else if (i < length(lines)) {
                merged <- paste0(lines[i], "\t", lines[i + 1])
                parts <- strsplit(merged, "\t")[[1]]
                nums <- suppressWarnings(as.numeric(parts))
                valid_nums <- nums[!is.na(nums)]
                
                if (length(valid_nums) >= 4) {
                    data_rows[[length(data_rows) + 1]] <- tail(valid_nums, 4)
                    i <- i + 2
                } else {
                    i <- i + 1
                }
            } else {
                i <- i + 1
            }
        }
    }
    
    if (length(data_rows) != n_methods) {
        warning(sprintf("Expected 4 data rows in %s, found %d", filepath, length(data_rows)))
        return(NULL)
    }
    
    mat <- matrix(unlist(data_rows), nrow = n_methods, ncol = n_methods, byrow = TRUE)
    rownames(mat) <- method_names
    colnames(mat) <- method_names
    return(mat)
}

# *** NEW PATH ***
conc_base <- "/mnt/genomics/pilot_PPR/concordance_IlluminaAvitiGATKPPR_new"

variant_types <- c("all", "snps", "indels")
concordance_data <- data.frame()

sample_dirs <- list.dirs(conc_base, recursive = FALSE, full.names = TRUE)
sample_dirs <- sample_dirs[grepl("sample_ICDG", sample_dirs)]

cat("\nReading concordance matrices from", length(sample_dirs), "directories...\n")
cat("  Base path:", conc_base, "\n")

parsed_ok <- 0
parsed_fail <- 0

for (dir in sample_dirs) {
    sample_id <- sub(".*sample_", "", basename(dir))
    sample_num <- as.integer(sub("ICDG", "", sample_id))
    
    for (vt in variant_types) {
        gt_file <- file.path(dir, paste0(sample_id, "_", vt, "_genotype_concordance.tsv"))
        if (!file.exists(gt_file)) next
        
        mat <- read_concordance_matrix(gt_file)
        if (is.null(mat)) {
            parsed_fail <- parsed_fail + 1
            next
        }
        
        parsed_ok <- parsed_ok + 1
        
        # Matrix indices: Illumina_PPR(1), Illumina_GATK(2), Aviti_PPR(3), Aviti_GATK(4)
        pairs <- list(
            list(name = "IlluminaPPR_vs_AvitiPPR", i = 1, j = 3, type = "cross_platform_same_method"),
            list(name = "IlluminaGATK_vs_AvitiGATK", i = 2, j = 4, type = "cross_platform_same_method"),
            list(name = "IlluminaPPR_vs_IlluminaGATK", i = 1, j = 2, type = "same_platform_illumina"),
            list(name = "AvitiPPR_vs_AvitiGATK", i = 3, j = 4, type = "same_platform_aviti"),
            list(name = "IlluminaPPR_vs_AvitiGATK", i = 1, j = 4, type = "cross_both"),
            list(name = "IlluminaGATK_vs_AvitiPPR", i = 2, j = 3, type = "cross_both")
        )
        
        for (pair in pairs) {
            concordance_data <- rbind(concordance_data, data.frame(
                sample_num = sample_num,
                variant_type = vt,
                comparison = pair$name,
                comparison_type = pair$type,
                gt_concordance = mat[pair$i, pair$j],
                stringsAsFactors = FALSE
            ))
        }
    }
}

cat("  Successfully parsed:", parsed_ok, "matrices\n")
cat("  Failed to parse:", parsed_fail, "matrices\n")
cat("  Concordance records:", nrow(concordance_data), "\n")
cat("  Unique samples:", length(unique(concordance_data$sample_num)), "\n")

if (nrow(concordance_data) == 0) {
    stop("ERROR: No concordance data was parsed. Check file format in concordance directories.")
}

# ----------------------------------------------------------------------------
# 3. MERGE & EXCLUDE OUTLIER
# ----------------------------------------------------------------------------

concordance_data <- concordance_data[concordance_data$sample_num != 23, ]
cov_wide <- cov_wide[cov_wide$sample_num != 23, ]

merged_conc <- merge(concordance_data, cov_wide, by = "sample_num")

cat("\nMerged data:", nrow(merged_conc), "rows (excluding ICDG23)\n")

# ----------------------------------------------------------------------------
# 4. CORRELATION ANALYSIS
# ----------------------------------------------------------------------------

cat("\n============================================================\n")
cat("CORRELATION RESULTS: Coverage vs Genotype Concordance\n")
cat("Separate correlations for Illumina and Aviti coverage\n")
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

results <- data.frame()

comparisons_unique <- unique(merged_conc$comparison)
variant_types_unique <- unique(merged_conc$variant_type)

for (vt in variant_types_unique) {
    cat(sprintf("\n--- %s variants ---\n", toupper(vt)))
    
    for (comp in comparisons_unique) {
        subset_df <- merged_conc[merged_conc$variant_type == vt &
                                  merged_conc$comparison == comp, ]
        if (nrow(subset_df) < 5) next
        
        # (a) Illumina coverage vs concordance
        corr_ill <- calc_correlation(subset_df$cov_Illumina, subset_df$gt_concordance)
        
        results <- rbind(results, data.frame(
            variant_type = vt,
            comparison = comp,
            comparison_type = subset_df$comparison_type[1],
            coverage_platform = "Illumina",
            n = corr_ill$n,
            pearson_r = corr_ill$pearson_r,
            pearson_p = corr_ill$pearson_p,
            spearman_rho = corr_ill$spearman_rho,
            spearman_p = corr_ill$spearman_p,
            stringsAsFactors = FALSE
        ))
        
        sig_ill <- ifelse(is.na(corr_ill$pearson_p), "NA",
                   ifelse(corr_ill$pearson_p < 0.001, "***",
                    ifelse(corr_ill$pearson_p < 0.01, "**",
                     ifelse(corr_ill$pearson_p < 0.05, "*", "ns"))))
        
        # (b) Aviti coverage (swapped) vs concordance
        corr_avi <- calc_correlation(subset_df$cov_Aviti_swapped, subset_df$gt_concordance)
        
        results <- rbind(results, data.frame(
            variant_type = vt,
            comparison = comp,
            comparison_type = subset_df$comparison_type[1],
            coverage_platform = "Aviti",
            n = corr_avi$n,
            pearson_r = corr_avi$pearson_r,
            pearson_p = corr_avi$pearson_p,
            spearman_rho = corr_avi$spearman_rho,
            spearman_p = corr_avi$spearman_p,
            stringsAsFactors = FALSE
        ))
        
        sig_avi <- ifelse(is.na(corr_avi$pearson_p), "NA",
                   ifelse(corr_avi$pearson_p < 0.001, "***",
                    ifelse(corr_avi$pearson_p < 0.01, "**",
                     ifelse(corr_avi$pearson_p < 0.05, "*", "ns"))))
        
        cat(sprintf("  %-35s | Illumina cov: r=%.3f (%s) | Aviti cov: r=%.3f (%s)\n",
                    comp, corr_ill$pearson_r, sig_ill, corr_avi$pearson_r, sig_avi))
    }
}

# ----------------------------------------------------------------------------
# 5. SUMMARY
# ----------------------------------------------------------------------------

cat("\n\n============================================================\n")
cat("SUMMARY: Significant correlations (p < 0.05)\n")
cat("============================================================\n\n")

sig_results <- results[!is.na(results$pearson_p) & results$pearson_p < 0.05, ]

if (nrow(sig_results) > 0) {
    sig_results <- sig_results[order(sig_results$pearson_p), ]
    for (i in 1:nrow(sig_results)) {
        cat(sprintf("  %s | %s | cov=%s | r=%.3f (p=%.4f) | rho=%.3f\n",
                    toupper(sig_results$variant_type[i]),
                    sig_results$comparison[i],
                    sig_results$coverage_platform[i],
                    sig_results$pearson_r[i],
                    sig_results$pearson_p[i],
                    sig_results$spearman_rho[i]))
    }
} else {
    cat("  No significant correlations found at p < 0.05\n")
}

# ----------------------------------------------------------------------------
# 6. MANUSCRIPT VALUES CHECK
# The following values appear in the manuscript for the coverage-concordance
# section. After running, compare these CSV values with the claims below:
#
#   Manuscript claim                          | Look for in CSV
#   ------------------------------------------|------------------------------------------
#   "Illumina cov vs PPR-GATK on Illumina:   | all, IlluminaPPR_vs_IlluminaGATK, Illumina
#    r = -0.822, p < 0.001"                   |
#   "Aviti cov vs PPR cross-platform:         | all, IlluminaPPR_vs_AvitiPPR, Aviti
#    r = 0.517, p = 0.012"                    |
#   "Illumina cov vs GATK cross-platform:     | all, IlluminaGATK_vs_AvitiGATK, Illumina
#    r = 0.547, p = 0.007"                    |
#   SNVs within-platform: r = -0.631, p=0.001 | snps, IlluminaPPR_vs_IlluminaGATK, Illumina
#   Indels within-platform: r = -0.523, p=0.011 | indels, IlluminaPPR_vs_IlluminaGATK, Illumina
#   GATK cross-platform indels: r = 0.686    | indels, IlluminaGATK_vs_AvitiGATK, Illumina
#   PPR cross-platform SNVs: r = 0.570       | snps, IlluminaPPR_vs_AvitiPPR, Aviti
#
# If new concordance data changes these values, UPDATE the manuscript.
# ----------------------------------------------------------------------------

cat("\n\n============================================================\n")
cat("MANUSCRIPT REFERENCE VALUES (from old concordance)\n")
cat("Compare with NEW results above. If different, update text.\n")
cat("============================================================\n")
cat("  Illum cov vs PPR-GATK on Illumina:      OLD r = -0.822\n")
cat("  Aviti cov vs PPR cross-platform:         OLD r = 0.517\n")
cat("  Illum cov vs GATK cross-platform:        OLD r = 0.547\n")
cat("  SNVs within-platform (Illum):            OLD r = -0.631\n")
cat("  Indels within-platform (Illum):          OLD r = -0.523\n")
cat("  GATK cross-platform indels (Illum cov):  OLD r = 0.686\n")
cat("  PPR cross-platform SNVs (Aviti cov):     OLD r = 0.570\n")
cat("============================================================\n\n")

# ----------------------------------------------------------------------------
# 7. SAVE
# ----------------------------------------------------------------------------

output_dir <- "/mnt/genomics/pilot_PPR/coverage_correlations_again"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(results, file.path(output_dir, "correlation_coverage_concordance_results.csv"),
          row.names = FALSE)

# Save merged data with both coverage columns for plotting
merged_conc_export <- merged_conc[, c("sample_num", "variant_type", "comparison",
                                       "comparison_type", "gt_concordance",
                                       "cov_Illumina", "cov_Aviti_swapped")]
write.csv(merged_conc_export, file.path(output_dir, "correlation_coverage_concordance_data.csv"),
          row.names = FALSE)

cat("Results saved to:", output_dir, "\n")
cat("  - correlation_coverage_concordance_results.csv\n")
cat("  - correlation_coverage_concordance_data.csv\n")
cat("============================================================\n")
