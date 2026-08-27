#!/usr/bin/env Rscript

## =============================================================================
## WILCOXON PER-REGION STATISTICAL ANALYSIS OF CONCORDANCE RESULTS
## Every test runs within a single genomic region -- no cross-region
## averaging or pooling anywhere in this script.
##
## Why: an earlier version of this analysis averaged GT_Conc_Pct across the
## 9 genomic-context BED regions to get one value per sample before testing.
## Those regions overlap substantially (curated functional/difficulty
## categories, not a genome partition), so pooling across them -- weighted
## or not -- does not correspond to a well-defined genome-wide quantity.
## This script never combines information across regions; each of the 9
## regions gets its own, self-contained test.
##
## Independence note: within a single region's test, the 23 samples ARE
## independent of each other (different individuals) -- that is exactly
## what the paired Wilcoxon test requires, and it holds here. What is NOT
## independent is the set of 9 per-region tests taken together, since the
## same 23 samples contribute to every one of them. Holm's step-down
## procedure is proven to control the family-wise error rate under
## ARBITRARY dependence between tests, so the correction remains valid
## here regardless. What changes is only the INTERPRETATION: "9/9 regions
## significant" means the effect replicates consistently across genomic
## contexts within this one 23-sample cohort -- not 9 independent
## confirmations from 9 distinct groups of people.
##
## H1: Cross-Platform (mean of PPR_cross_platform, GATK_cross_platform --
##     both already local to a single region, so this averaging is safe)
##     vs Cross-Method (mean of the two same-platform comparisons, also
##     region-local)
## H2: PPR_cross_platform vs GATK_cross_platform directly (no averaging
##     needed -- both are already single columns for that region)
##
## For each hypothesis: 9 independent-per-region paired Wilcoxon
## signed-rank tests (n = 23 samples each), Holm-corrected across the 9,
## separately for each variant type (all / snps / indels).
##
## Input:  ../concordance_summary_for_analysis_corrected.tsv
##         (produced by generate_concordance_summary_corrected.R)
## Output: posthoc_H1_crossplatform_vs_crossmethod.csv
##         posthoc_H2_PPR_vs_GATK_crossplatform.csv
## =============================================================================

library(tidyverse)

base_dir <- "/mnt/genomics/pilot_PPR/concordance_perbed"
output_dir <- file.path(base_dir, "stats_wilcoxon_perregion")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

data_file <- file.path(base_dir, "concordance_summary_for_analysis_corrected.tsv")
if (!file.exists(data_file)) {
    cat("ERROR: Data file not found:", data_file, "\n")
    quit(status = 1)
}

data <- read.table(data_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

data <- data %>%
    filter(!is.na(GT_Conc_Pct)) %>%
    mutate(
        Specific_Comparison = case_when(
            Comparison == "IlluminaPPR_vs_AvitiPPR" ~ "PPR_cross_platform",
            Comparison == "IlluminaGATK_vs_AvitiGATK" ~ "GATK_cross_platform",
            Comparison == "IlluminaPPR_vs_IlluminaGATK" ~ "PPR_vs_GATK_Illumina",
            Comparison == "AvitiPPR_vs_AvitiGATK" ~ "PPR_vs_GATK_Aviti"
        ),
        Sample = factor(Sample),
        BED_Region = factor(BED_Region)
    )

n_samples <- length(unique(data$Sample))
n_regions <- length(unique(data$BED_Region))

cat("=================================================================\n")
cat("WILCOXON PER-REGION CONCORDANCE ANALYSIS (no cross-region pooling)\n")
cat("=================================================================\n\n")
cat(sprintf("Samples: %d | Regions: %d (%s)\n\n", n_samples, n_regions,
            paste(sort(unique(data$BED_Region)), collapse = ", ")))

if (n_samples != 23) cat(sprintf("WARNING: expected 23 samples, found %d\n\n", n_samples))
if (n_regions != 9) cat(sprintf("WARNING: expected 9 regions, found %d\n\n", n_regions))

posthoc_H1_results <- data.frame()
posthoc_H2_results <- data.frame()

for (vt in c("all", "snps", "indels")) {

    cat(strrep("=", 65), "\n")
    cat(" VARIANT TYPE:", toupper(vt), "\n")
    cat(strrep("=", 65), "\n\n")

    vt_data <- data %>% filter(Variant_Type == vt)

    # -------------------------------------------------------------------
    # H1: Cross-Platform vs Cross-Method, PER REGION
    # -------------------------------------------------------------------
    cat("--- H1: Cross-Platform vs Cross-Method, per region ---\n")

    h1_data <- vt_data %>%
        select(Sample, BED_Region, Specific_Comparison, GT_Conc_Pct) %>%
        pivot_wider(names_from = Specific_Comparison, values_from = GT_Conc_Pct) %>%
        mutate(
            Cross_Platform = (PPR_cross_platform + GATK_cross_platform) / 2,
            Cross_Method = (PPR_vs_GATK_Illumina + PPR_vs_GATK_Aviti) / 2
        ) %>%
        filter(!is.na(Cross_Platform) & !is.na(Cross_Method))

    for (reg in levels(data$BED_Region)) {
        reg_data <- h1_data %>% filter(BED_Region == reg)
        if (nrow(reg_data) < 3) next

        wt <- wilcox.test(reg_data$Cross_Platform, reg_data$Cross_Method,
                          paired = TRUE, exact = FALSE)
        diff_vals <- reg_data$Cross_Platform - reg_data$Cross_Method
        z_val <- qnorm(wt$p.value / 2)
        r_effect <- min(abs(z_val) / sqrt(nrow(reg_data)), 1)

        cat(sprintf("  %-8s N=%d  CrossPlat=%.1f%%  CrossMeth=%.1f%%  Diff=%+.1fpp  V=%.0f  p=%.2e  r=%.2f\n",
                    reg, nrow(reg_data), mean(reg_data$Cross_Platform), mean(reg_data$Cross_Method),
                    mean(diff_vals), wt$statistic, wt$p.value, r_effect))

        posthoc_H1_results <- rbind(posthoc_H1_results, data.frame(
            Variant_Type = vt, BED_Region = reg,
            Mean_CrossPlatform = mean(reg_data$Cross_Platform),
            Mean_CrossMethod = mean(reg_data$Cross_Method),
            Mean_Difference_pp = mean(diff_vals),
            V_statistic = as.numeric(wt$statistic),
            p_value_raw = wt$p.value,
            effect_size_r = r_effect,
            N = nrow(reg_data),
            stringsAsFactors = FALSE
        ))
    }
    cat("\n")

    # -------------------------------------------------------------------
    # H2: PPR cross-platform vs GATK cross-platform, PER REGION
    # -------------------------------------------------------------------
    cat("--- H2: PPR vs GATK cross-platform, per region ---\n")

    h2_data <- vt_data %>%
        select(Sample, BED_Region, Specific_Comparison, GT_Conc_Pct) %>%
        pivot_wider(names_from = Specific_Comparison, values_from = GT_Conc_Pct) %>%
        select(Sample, BED_Region, PPR_cross_platform, GATK_cross_platform) %>%
        filter(!is.na(PPR_cross_platform) & !is.na(GATK_cross_platform))

    for (reg in levels(data$BED_Region)) {
        reg_data <- h2_data %>% filter(BED_Region == reg)
        if (nrow(reg_data) < 3) next

        wt2 <- wilcox.test(reg_data$PPR_cross_platform, reg_data$GATK_cross_platform,
                           paired = TRUE, exact = FALSE)
        diff_vals2 <- reg_data$PPR_cross_platform - reg_data$GATK_cross_platform
        z_val2 <- qnorm(wt2$p.value / 2)
        r_effect2 <- min(abs(z_val2) / sqrt(nrow(reg_data)), 1)

        cat(sprintf("  %-8s N=%d  PPR=%.1f%%  GATK=%.1f%%  Diff=%+.1fpp  V=%.0f  p=%.2e  r=%.2f\n",
                    reg, nrow(reg_data), mean(reg_data$PPR_cross_platform), mean(reg_data$GATK_cross_platform),
                    mean(diff_vals2), wt2$statistic, wt2$p.value, r_effect2))

        posthoc_H2_results <- rbind(posthoc_H2_results, data.frame(
            Variant_Type = vt, BED_Region = reg,
            Mean_PPR = mean(reg_data$PPR_cross_platform),
            Mean_GATK = mean(reg_data$GATK_cross_platform),
            Mean_Difference_pp = mean(diff_vals2),
            V_statistic = as.numeric(wt2$statistic),
            p_value_raw = wt2$p.value,
            effect_size_r = r_effect2,
            N = nrow(reg_data),
            stringsAsFactors = FALSE
        ))
    }
    cat("\n")
}

# =============================================================================
# HOLM CORRECTION across the 9 regions, within each Variant_Type x Hypothesis
# =============================================================================
posthoc_H1_results <- posthoc_H1_results %>%
    group_by(Variant_Type) %>%
    mutate(p_value_holm = p.adjust(p_value_raw, method = "holm"),
           Significant = p_value_holm < 0.05) %>%
    ungroup()

posthoc_H2_results <- posthoc_H2_results %>%
    group_by(Variant_Type) %>%
    mutate(p_value_holm = p.adjust(p_value_raw, method = "holm"),
           Significant = p_value_holm < 0.05) %>%
    ungroup()

# =============================================================================
# SAVE
# =============================================================================
write.csv(posthoc_H1_results, file.path(output_dir, "posthoc_H1_crossplatform_vs_crossmethod.csv"), row.names = FALSE)
write.csv(posthoc_H2_results, file.path(output_dir, "posthoc_H2_PPR_vs_GATK_crossplatform.csv"), row.names = FALSE)

cat(strrep("=", 65), "\n")
cat("DONE. Results saved to:", output_dir, "/\n")
cat("  posthoc_H1_crossplatform_vs_crossmethod.csv   (9 regions x 3 variant types = 27 rows)\n")
cat("  posthoc_H2_PPR_vs_GATK_crossplatform.csv      (9 regions x 3 variant types = 27 rows)\n")
cat(strrep("=", 65), "\n")

# Quick summary
cat("\nSUMMARY (variant type = all):\n")
h1_all <- posthoc_H1_results %>% filter(Variant_Type == "all")
h2_all <- posthoc_H2_results %>% filter(Variant_Type == "all")
cat(sprintf("  H1 (Cross-Platform > Cross-Method): %d/%d regions significant (Holm), all favoring Cross-Platform: %s\n",
            sum(h1_all$Significant), nrow(h1_all),
            all(h1_all$Mean_Difference_pp[h1_all$Significant] > 0)))
cat(sprintf("  H2 (PPR > GATK, cross-platform):    %d/%d regions significant (Holm), all favoring PPR: %s\n",
            sum(h2_all$Significant), nrow(h2_all),
            all(h2_all$Mean_Difference_pp[h2_all$Significant] > 0)))

cat("\nNOTE: significance across the 9 regions is not 9 independent confirmations\n")
cat("(the same 23 samples contribute to every region's test); read it as\n")
cat("'the effect replicates consistently across genomic contexts in this cohort'.\n")
