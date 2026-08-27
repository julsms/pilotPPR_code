# Coverage Correlations

Correlation analysis between sequencing depth and variant quality metrics (Figure 4C–D, Supplementary Figures S9–S10).

## Scripts

| Script | Description |
|--------|-------------|
| `correlation_coverage_validation.R` | Correlates per-sample mean coverage with Merfin k-mer validation rates for all variants and pipeline-specific unique variants. Analyzes each platform independently (n=23 per platform, excluding one low-coverage outlier at 1.23×). |
| `correlation_coverage_concordance.R` | Correlates per-sample mean coverage with genotype concordance across all six pairwise comparisons (cross-platform and cross-method). Tests Illumina and Aviti coverage separately as predictors. Applies sample swap correction for ICDG12/ICDG21. |

## Input Data

- **Coverage metrics**: mosdepth output (`qc_metrics_collected_for_plots.tsv`)
- **Merfin validation (all variants)**: Summary files from `final_results/summary_*.txt`
- **Merfin validation (unique variants)**: `unique_validation_rates_corrected.csv`
- **Concordance matrices**: Per-sample genotype concordance TSVs from the concordance analysis

## Output

- `correlation_coverage_validation_results.csv` — Pearson and Spearman correlation statistics for coverage vs. validation
- `correlation_coverage_validation_data.csv` — Merged per-sample data used for scatter plots
- `correlation_coverage_concordance_results.csv` — Correlation statistics for coverage vs. concordance
- `correlation_coverage_concordance_data.csv` — Merged data for plotting

## Methods Summary

Both Pearson and Spearman correlations are reported. Sample ICDG23 (coverage = 1.23×) is excluded as an outlier. For concordance correlations involving cross-platform comparisons, the Aviti coverage values for samples ICDG12 and ICDG21 are swapped to match the known sample swap applied during concordance analysis. Statistical significance was set at p < 0.05.

## Dependencies

- R ≥ 4.3
- `dplyr`, `tidyr`
