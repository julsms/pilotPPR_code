#!/usr/bin/env Rscript

# ============================================================================
# STATISTICAL ANALYSIS FOR BENCHMARK COMPARISON
# Friedman test + paired Wilcoxon signed-rank tests with effect sizes
# 
# Rationale: Standard two-way ANOVA assumes independent observations.
# However, the 10 genomic stratification regions are evaluated on the SAME
# sample (HG002) and the same sequencing data — they are not independent
# replicates. Some regions even overlap (AD = union of LM + TRHP + OD).
#
# Friedman test is the non-parametric equivalent of repeated-measures ANOVA:
# - Treats regions as BLOCKS (not independent replicates)
# - Does not assume normality
# - Does not assume independence between observations
# - Tests: "Are the differences between pipelines consistent across blocks?"
#
# Post-hoc: Paired Wilcoxon signed-rank tests with Holm correction
# Effect sizes: Kendall's W (Friedman) and matched-pairs rank-biserial r (Wilcoxon)
# ============================================================================

library(tidyverse)
library(rstatix)     # for friedman_test, wilcox_test, friedman_effsize
library(broom)

# ============================================================================
# INPUT: Adjust this path to your all_benchmark_data.csv
# ============================================================================
# The CSV should have columns: Region, Pipeline, Type, Filter, 
#   METRIC.Precision (or Precision), METRIC.Recall (or Recall), 
#   METRIC.F1_Score (or F1_Score)
#
# Each row = one region × one pipeline combination
# Expected: 10 regions × 4 pipelines = 40 rows per Type × Filter

input_file <- "all_benchmark_data.csv"
output_dir <- "benchmark_stats_results"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

cat("=================================================================\n")
cat("BENCHMARK STATISTICAL ANALYSIS\n")
cat("Friedman test + paired Wilcoxon signed-rank (non-parametric)\n")
cat("=================================================================\n\n")

# ============================================================================
# READ AND PREPARE DATA
# ============================================================================

data <- read.csv(input_file)
cat("Input file:", input_file, "\n")
cat("Dimensions:", nrow(data), "rows x", ncol(data), "columns\n")

# Detect column naming (handle both METRIC.X and X naming)
if ("METRIC.Precision" %in% colnames(data)) {
  data <- data %>%
    rename(Precision = METRIC.Precision, Recall = METRIC.Recall, F1_Score = METRIC.F1_Score)
}

# Filter to PASS variants only
analysis_data <- data %>%
  filter(Filter == "PASS") %>%
  mutate(
    Platform = case_when(
      str_detect(Pipeline, "Illumina") ~ "Illumina",
      str_detect(Pipeline, "Aviti") ~ "Aviti"
    ),
    Method = case_when(
      str_detect(Pipeline, "GATK") ~ "GATK",
      str_detect(Pipeline, "PPR") ~ "PPR"
    ),
    Pipeline = factor(Pipeline),
    Platform = factor(Platform, levels = c("Illumina", "Aviti")),
    Method = factor(Method, levels = c("GATK", "PPR")),
    Type = factor(Type),
    Region = factor(Region)
  ) %>%
  # Pivot to long format for metrics
  pivot_longer(
    cols = c(Precision, Recall, F1_Score),
    names_to = "Metric",
    values_to = "Score"
  ) %>%
  mutate(Metric = factor(Metric, levels = c("Precision", "Recall", "F1_Score")))

cat("Pipelines:", paste(levels(analysis_data$Pipeline), collapse = ", "), "\n")
cat("Regions:", paste(levels(analysis_data$Region), collapse = ", "), "\n")
cat("Types:", paste(levels(analysis_data$Type), collapse = ", "), "\n\n")

# ============================================================================
# 1. FRIEDMAN TEST — overall difference among 4 pipelines
#    (Regions = blocks, Pipelines = treatments)
# ============================================================================

cat("=================================================================\n")
cat("1. FRIEDMAN TEST (non-parametric repeated-measures)\n")
cat("   H0: No difference among the 4 pipelines across regions\n")
cat("   Regions treated as blocks (same sample, different contexts)\n")
cat("=================================================================\n\n")

friedman_results <- data.frame()

for (vtype in c("SNP", "INDEL")) {
  for (metric in c("Precision", "Recall", "F1_Score")) {
    
    df <- analysis_data %>%
      filter(Type == vtype, Metric == metric) %>%
      select(Region, Pipeline, Score)
    
    # Friedman test requires: group (Pipeline), block (Region), value (Score)
    # rstatix::friedman_test wants data in long format with id column
    ft <- df %>%
      friedman_test(Score ~ Pipeline | Region)
    
    # Effect size: Kendall's W
    # W ranges 0-1: 0 = no agreement, 1 = perfect agreement in rankings
    eff <- df %>%
      friedman_effsize(Score ~ Pipeline | Region)
    
    result_row <- data.frame(
      Type = vtype,
      Metric = metric,
      Friedman_chi2 = ft$statistic,
      df = ft$df,
      p_value = ft$p,
      Kendall_W = eff$effsize,
      W_magnitude = eff$magnitude,
      significant = ft$p < 0.05
    )
    
    friedman_results <- rbind(friedman_results, result_row)
    
    cat(sprintf("%-5s %-10s: χ²(%d) = %.2f, p = %.4f, Kendall W = %.3f (%s) %s\n",
                vtype, metric, ft$df, ft$statistic, ft$p,
                eff$effsize, eff$magnitude,
                ifelse(ft$p < 0.05, "***", "ns")))
  }
}

cat("\n")
write.csv(friedman_results, file.path(output_dir, "friedman_test_results.csv"), row.names = FALSE)

# ============================================================================
# 2. POST-HOC: Paired Wilcoxon signed-rank tests
#    Specific comparisons of interest:
#    a) Method effect: GATK vs PPR (on same platform, paired by region)
#    b) Platform effect: Illumina vs Aviti (same method, paired by region)
#    c) All pairwise pipeline comparisons
# ============================================================================

cat("=================================================================\n")
cat("2. PAIRED WILCOXON SIGNED-RANK TESTS (post-hoc)\n")
cat("   Paired by region (same genomic context)\n")
cat("=================================================================\n\n")

wilcox_results <- data.frame()

for (vtype in c("SNP", "INDEL")) {
  for (metric in c("Precision", "Recall", "F1_Score")) {
    
    # Reshape to wide: one column per pipeline, rows = regions
    wide_df <- analysis_data %>%
      filter(Type == vtype, Metric == metric) %>%
      select(Region, Pipeline, Score) %>%
      pivot_wider(names_from = Pipeline, values_from = Score)
    
    # Get pipeline column names
    pipe_cols <- setdiff(colnames(wide_df), "Region")
    
    # Define comparisons of interest
    # Adjust these names to match your actual pipeline names
    comparisons <- list(
      # Method effect (same platform)
      list(name = "Method: PPR vs GATK (Illumina)", 
           col_a = grep("Illumina.*PPR|PPR.*Illumina", pipe_cols, value = TRUE),
           col_b = grep("Illumina.*GATK|GATK.*Illumina", pipe_cols, value = TRUE)),
      list(name = "Method: PPR vs GATK (Aviti)",
           col_a = grep("Aviti.*PPR|PPR.*Aviti", pipe_cols, value = TRUE),
           col_b = grep("Aviti.*GATK|GATK.*Aviti", pipe_cols, value = TRUE)),
      # Platform effect (same method)
      list(name = "Platform: Aviti vs Illumina (GATK)",
           col_a = grep("Aviti.*GATK|GATK.*Aviti", pipe_cols, value = TRUE),
           col_b = grep("Illumina.*GATK|GATK.*Illumina", pipe_cols, value = TRUE)),
      list(name = "Platform: Aviti vs Illumina (PPR)",
           col_a = grep("Aviti.*PPR|PPR.*Aviti", pipe_cols, value = TRUE),
           col_b = grep("Illumina.*PPR|PPR.*Illumina", pipe_cols, value = TRUE))
    )
    
    for (comp in comparisons) {
      if (length(comp$col_a) == 1 && length(comp$col_b) == 1) {
        vals_a <- wide_df[[comp$col_a]]
        vals_b <- wide_df[[comp$col_b]]
        
        # Skip if all differences are zero (no variation)
        diffs <- vals_a - vals_b
        if (all(diffs == 0, na.rm = TRUE)) {
          next
        }
        
        # Paired Wilcoxon signed-rank test
        wt <- wilcox.test(vals_a, vals_b, paired = TRUE, exact = FALSE)
        
        # Effect size: matched-pairs rank-biserial correlation r
        # r = Z / sqrt(N), where N = number of pairs
        # Using rstatix for cleaner calculation
        paired_df <- data.frame(
          score = c(vals_a, vals_b),
          group = factor(rep(c("A", "B"), each = length(vals_a))),
          id = factor(rep(1:length(vals_a), 2))
        )
        
        # Calculate r effect size manually
        # r = 1 - (2*W) / (n*(n+1)/2) where W is the smaller rank sum
        n_pairs <- sum(!is.na(diffs) & diffs != 0)
        if (n_pairs > 0) {
          # Z approximation
          z_val <- qnorm(wt$p.value / 2) * sign(mean(diffs, na.rm = TRUE))
          r_effect <- abs(z_val) / sqrt(n_pairs)
          r_effect <- min(r_effect, 1)  # cap at 1
        } else {
          r_effect <- NA
        }
        
        # Median difference
        med_diff <- median(diffs, na.rm = TRUE)
        mean_diff <- mean(diffs, na.rm = TRUE)
        
        wilcox_row <- data.frame(
          Type = vtype,
          Metric = metric,
          Comparison = comp$name,
          Median_A = median(vals_a, na.rm = TRUE),
          Median_B = median(vals_b, na.rm = TRUE),
          Median_Diff = med_diff,
          Mean_Diff = mean_diff,
          W_statistic = as.numeric(wt$statistic),
          p_value = wt$p.value,
          effect_size_r = r_effect,
          n_pairs = n_pairs,
          significant = wt$p.value < 0.05
        )
        
        wilcox_results <- rbind(wilcox_results, wilcox_row)
      }
    }
  }
}

# Apply multiple comparison corrections
if (nrow(wilcox_results) > 0) {
  wilcox_results <- wilcox_results %>%
    group_by(Type, Metric) %>%
    mutate(
      p_holm = p.adjust(p_value, method = "holm"),
      p_fdr = p.adjust(p_value, method = "fdr"),
      significant_holm = p_holm < 0.05,
      significant_fdr = p_fdr < 0.05
    ) %>%
    ungroup()
}

cat("Post-hoc paired Wilcoxon results:\n\n")
for (vtype in c("SNP", "INDEL")) {
  cat(sprintf("--- %s ---\n", vtype))
  subset_res <- wilcox_results %>% filter(Type == vtype)
  for (i in 1:nrow(subset_res)) {
    row <- subset_res[i, ]
    cat(sprintf("  %-10s | %-35s | W=%5.1f, p=%.4f (Holm: %.4f) | r=%.3f | Δmedian=%+.4f %s\n",
                row$Metric, row$Comparison, row$W_statistic, row$p_value, row$p_holm,
                row$effect_size_r, row$Median_Diff,
                ifelse(row$significant_holm, "***", "ns")))
  }
  cat("\n")
}

write.csv(wilcox_results, file.path(output_dir, "wilcoxon_paired_results.csv"), row.names = FALSE)

# ============================================================================
# 3. TWO-WAY ANOVA (kept for comparison / reviewer expectations)
#    With explicit note that this is supplementary to the non-parametric tests
# ============================================================================

cat("=================================================================\n")
cat("3. TWO-WAY ANOVA (for reference — supplementary to Friedman)\n")
cat("   NOTE: Assumes independence of regions (approximation)\n")
cat("=================================================================\n\n")

anova_results <- data.frame()

for (vtype in c("SNP", "INDEL")) {
  for (metric in c("Precision", "Recall", "F1_Score")) {
    
    df <- analysis_data %>%
      filter(Type == vtype, Metric == metric)
    
    model <- aov(Score ~ Platform * Method, data = df)
    anova_tidy <- tidy(model)
    anova_tidy$Type <- vtype
    anova_tidy$Metric <- metric
    
    # Eta-squared effect sizes
    ss <- anova_tidy$sumsq
    ss_total <- sum(ss)
    anova_tidy$eta_squared <- ss / ss_total
    
    anova_results <- rbind(anova_results, anova_tidy)
    
    # Print compact
    for (j in 1:(nrow(anova_tidy) - 1)) {  # skip Residuals
      row <- anova_tidy[j, ]
      cat(sprintf("  %-5s %-10s | %-18s: F(%d,%d)=%6.2f, p=%.4f, η²=%.3f %s\n",
                  vtype, metric, row$term,
                  row$df, anova_tidy$df[nrow(anova_tidy)],
                  ifelse(is.na(row$statistic), 0, row$statistic),
                  ifelse(is.na(row$p.value), 1, row$p.value),
                  row$eta_squared,
                  ifelse(!is.na(row$p.value) && row$p.value < 0.05, "***", "ns")))
    }
  }
}

cat("\n")
write.csv(anova_results, file.path(output_dir, "anova_results_reference.csv"), row.names = FALSE)

# ============================================================================
# 4. DESCRIPTIVE STATISTICS
# ============================================================================

cat("=================================================================\n")
cat("4. DESCRIPTIVE STATISTICS\n")
cat("=================================================================\n\n")

descriptive <- analysis_data %>%
  group_by(Type, Metric, Pipeline, Platform, Method) %>%
  summarise(
    N = n(),
    Mean = mean(Score, na.rm = TRUE),
    SD = sd(Score, na.rm = TRUE),
    Median = median(Score, na.rm = TRUE),
    Min = min(Score, na.rm = TRUE),
    Max = max(Score, na.rm = TRUE),
    .groups = "drop"
  )

print(as.data.frame(descriptive))
write.csv(descriptive, file.path(output_dir, "descriptive_statistics.csv"), row.names = FALSE)

# ============================================================================
# 5. CONCORDANCE BETWEEN ANOVA AND NON-PARAMETRIC TESTS
# ============================================================================

cat("\n=================================================================\n")
cat("5. CONCORDANCE: ANOVA vs FRIEDMAN/WILCOXON\n")
cat("=================================================================\n\n")

# Compare conclusions
cat("Summary of significant effects:\n\n")
cat(sprintf("%-5s %-10s | %-20s | ANOVA p | Friedman p | Wilcoxon (Method) p\n",
            "Type", "Metric", "Effect"))
cat(paste(rep("-", 85), collapse = ""), "\n")

for (vtype in c("SNP", "INDEL")) {
  for (metric in c("Precision", "Recall", "F1_Score")) {
    
    # ANOVA method effect
    anova_method <- anova_results %>%
      filter(Type == vtype, Metric == metric, term == "Method")
    anova_p <- ifelse(nrow(anova_method) > 0, anova_method$p.value[1], NA)
    
    # Friedman overall
    fried_p <- friedman_results %>%
      filter(Type == vtype, Metric == metric) %>%
      pull(p_value)
    
    # Wilcoxon method effect (take min p across platforms)
    wilcox_method <- wilcox_results %>%
      filter(Type == vtype, Metric == metric, str_detect(Comparison, "Method"))
    wilcox_p <- ifelse(nrow(wilcox_method) > 0, min(wilcox_method$p_holm), NA)
    
    cat(sprintf("%-5s %-10s | %-20s | %.4f  | %.4f     | %.4f\n",
                vtype, metric, "Method (PPR vs GATK)",
                ifelse(is.na(anova_p), 1, anova_p),
                ifelse(length(fried_p) > 0, fried_p, 1),
                ifelse(is.na(wilcox_p), 1, wilcox_p)))
  }
}

# ============================================================================
# 6. GENERATE SUMMARY TABLE FOR PAPER
# ============================================================================

cat("\n=================================================================\n")
cat("6. SUMMARY TABLE FOR PAPER\n")
cat("=================================================================\n\n")

paper_table <- friedman_results %>%
  select(Type, Metric, Friedman_chi2, df, p_value, Kendall_W, W_magnitude, significant) %>%
  mutate(
    across(c(Friedman_chi2, Kendall_W), ~round(., 2)),
    p_value = format.pval(p_value, digits = 3)
  )

cat("Friedman test results (for paper):\n")
print(as.data.frame(paper_table))

cat("\nPaired Wilcoxon key comparisons (for paper):\n")
key_wilcox <- wilcox_results %>%
  select(Type, Metric, Comparison, Median_Diff, W_statistic, p_value, p_holm, effect_size_r) %>%
  mutate(
    across(c(Median_Diff, effect_size_r), ~round(., 4)),
    W_statistic = round(W_statistic, 1),
    p_value = format.pval(p_value, digits = 3),
    p_holm = format.pval(p_holm, digits = 3)
  )
print(as.data.frame(key_wilcox))

write.csv(paper_table, file.path(output_dir, "paper_friedman_table.csv"), row.names = FALSE)
write.csv(key_wilcox, file.path(output_dir, "paper_wilcoxon_table.csv"), row.names = FALSE)

# ============================================================================
# 7. FINAL SUMMARY
# ============================================================================

cat("\n=================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("=================================================================\n")
cat("\nOutput files saved to:", output_dir, "\n")
cat("  - friedman_test_results.csv     (primary non-parametric test)\n")
cat("  - wilcoxon_paired_results.csv   (post-hoc pairwise comparisons)\n")
cat("  - anova_results_reference.csv   (ANOVA for reference/comparison)\n")
cat("  - descriptive_statistics.csv    (summary stats per pipeline)\n")
cat("  - paper_friedman_table.csv      (formatted for paper)\n")
cat("  - paper_wilcoxon_table.csv      (formatted for paper)\n")
cat("\n")
cat("RECOMMENDED REPORTING STRATEGY:\n")
cat("  Primary: Friedman test (overall pipeline difference)\n")
cat("  Post-hoc: Paired Wilcoxon signed-rank with Holm correction\n")
cat("  Supplementary: Two-way ANOVA (for reviewer familiarity)\n")
cat("  Effect sizes: Kendall W (Friedman) + r (Wilcoxon)\n")
cat("=================================================================\n")
