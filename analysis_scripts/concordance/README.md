# Concordance Analysis

Per-sample genotype and position concordance between variant calling pipelines (GATK vs PPR) across sequencing platforms (Illumina vs Element Aviti). Results correspond to Figure 2A–D and Supplementary Figure S4.

## Scripts

| Script | Description |
|--------|-------------|
| `master_concordanceanalysis.sh` | Genome-wide concordance analysis for 23 samples. Downloads PPR VCFs from S3, prepares GATK VCFs (PASS-filtered), computes 4×4 pairwise concordance matrices (Illumina PPR, Illumina GATK, Aviti PPR, Aviti GATK) for ALL/SNP/INDEL variant types. Generates per-sample heatmaps and cohort-average summaries. |
| `concordance_perbed.sh` | Region-stratified concordance using GIAB v3.6 genomic stratification BED files. Intersects whole-genome VCFs with each of the 9 stratification regions and calculates concordance metrics per sample, per region, per variant type. |
| `generate_concordance_summary_corrected.R` | Collects per-sample, per-region concordance results into `concordance_summary_for_analysis.tsv`, then removes duplicate/invalid rows (see "Sample Notes" below) to produce `concordance_summary_for_analysis_corrected.tsv`. Must be run before the statistical analysis. |
| `wilcoxon_perregion_analysis.R` | Statistical analysis. Runs a separate paired Wilcoxon signed-rank test within each of the 9 genomic regions — no averaging or pooling of concordance values across regions (see rationale below). Tests two hypotheses (H1, H2), Holm-corrected across the 9 regions, for ALL/SNPs/INDELs separately. |

## Methodology

1. **VCF preparation**: Per-sample VCFs are extracted from multisample call sets, normalized (`bcftools norm -m -both`), restricted to autosomes + sex chromosomes (chr1–22, chrX, chrY), and filtered to exclude reference/missing genotypes (0/0, ./.).
2. **GATK filtering**: Only PASS variants (post-VQSR) are retained.
3. **PPR filtering**: Non-reference calls from DeepVariant (excludes 0/0 and ./.).
4. **MT exclusion**: Mitochondrial chromosomes are removed BEFORE chromosome renaming to avoid contig mismatch errors.
5. **Concordance calculation**: Position concordance = shared CHROM:POS:REF:ALT / union. Genotype concordance = additionally requires matching GT field.
6. **Variant type separation**: ALL, SNPs (`bcftools view -v snps`), INDELs (`bcftools view -v indels`).
7. **Deduplication**: raw per-region output can contain duplicate or invalid rows for a few samples (see Sample Notes); `generate_concordance_summary_corrected.R` removes these before any statistics are computed.
8. **Statistical testing**: paired Wilcoxon signed-rank tests, run independently within each of the 9 genomic regions (no cross-region averaging), Holm-corrected across regions.

## Statistical Analysis Design

The 9 genomic stratification regions used here are overlapping functional/difficulty categories (for example, AD is the union of LM + TRHP + OD), not a partition of the genome. Averaging concordance across them into one value per sample — the approach used earlier in this project — mixes regions of very different sizes and does not correspond to a well-defined quantity. `wilcoxon_perregion_analysis.R` instead runs one independent test per region, for each of two hypotheses:

### H1: Cross-Platform vs Cross-Method
- **Per region**: Cross-Platform = mean of (PPR Illumina↔Aviti, GATK Illumina↔Aviti); Cross-Method = mean of (PPR↔GATK on Illumina, PPR↔GATK on Aviti). Both quantities are computed within a single region, so this averaging does not mix across regions.
- **Test**: paired Wilcoxon signed-rank per region (n = 23), Holm-corrected across the 9 regions.
- **Interpretation**: higher Cross-Platform than Cross-Method concordance → variant calling method contributes more to detection differences than sequencing platform.

### H2: PPR vs GATK cross-platform reproducibility
- **Per region**: PPR cross-platform concordance vs GATK cross-platform concordance directly, no averaging needed.
- **Test**: paired Wilcoxon signed-rank per region (n = 23), Holm-corrected across the 9 regions.
- **Interpretation**: higher PPR value → PPR more reproducible across sequencing platforms than GATK, in that region.

Within each region's test the 23 samples are independent (different individuals) — the assumption the paired Wilcoxon test requires. What is not independent is the set of 9 per-region tests taken together, since the same 23 samples contribute to all of them; Holm's correction controls the family-wise error rate under this kind of dependence, so it stays valid here. A result of "9/9 regions significant" should be read as "the effect replicates consistently across genomic contexts in this cohort," not as 9 independent confirmations from 9 separate groups of people.

## Sample Notes

- **23 samples** are included in the published concordance analysis (ICDG1–ICDG24, one sample excluded due to low sequencing yield).
- Pairwise relatedness (KING-robust): ICDG6 and ICDG9 identified as second-degree relatives; both retained in concordance analyses (a single related sample is excluded only from singleton counts, in a separate part of the study).
- Sample ICDG23 is excluded from coverage-concordance correlation analyses due to extreme low coverage (outlier by IQR×3 criterion) — see `coverage_correlations/`.
- **Duplicate/invalid rows in the raw per-region output**, corrected by `generate_concordance_summary_corrected.R`:
  - ICDG12 and ICDG21 were each processed twice against Aviti data (the first run used incorrect data); only the second (correct) run is kept.
  - ICDG22 has identical duplicate rows; one copy is kept.
  - ICDG1, ICDG2, ICDG3 have empty/NA rows for the OMIM region from an incomplete first pass; these are dropped, keeping only the valid rows.

## Output

### Genome-wide concordance (`master_concordanceanalysis.sh`)
- Per-sample concordance matrices (`*_position_concordance.tsv`, `*_genotype_concordance.tsv`)
- Per-sample combined heatmaps (`*_concordance_heatmap.png`)
- Cohort-average matrices and heatmaps (`AVERAGE_*`)

### Region-stratified concordance (`concordance_perbed.sh`)
- Per-sample per-region concordance TSVs (`sample_*_concordance_results.tsv`)

### Corrected summary (`generate_concordance_summary_corrected.R`)
- `concordance_summary_for_analysis_corrected.tsv` — deduplicated, 23 samples × 9 regions × 4 comparisons × 3 variant types (2,484 rows)

### Statistical analysis (`wilcoxon_perregion_analysis.R`)
- `stats_wilcoxon_perregion/posthoc_H1_crossplatform_vs_crossmethod.csv` — per-region H1 results (9 regions × 3 variant types)
- `stats_wilcoxon_perregion/posthoc_H2_PPR_vs_GATK_crossplatform.csv` — per-region H2 results (9 regions × 3 variant types)

## Genomic Stratification Regions

Nine regions from GIAB v3.6 stratifications (the benchmark bed, BB, used for Figure 1 HG002 benchmarking, is not part of the concordance stratification):

| Abbreviation | Description |
|-------------|-------------|
| LM | Low-mappability regions |
| TRHP | Tandem repeats and homopolymers |
| AD | All difficult regions (union) |
| OD | Other difficult regions |
| 25GC65 | Moderate GC content (<25% or >65%) |
| 15GC85 | Extreme GC content (<15% or >85%) |
| REFSEQ | RefSeq coding regions |
| ACMG | ACMG Secondary Findings v3.3 genes |
| OMIM | OMIM disease genes |

## Dependencies

- `bcftools` ≥ 1.17, `tabix`, `bgzip`
- AWS CLI (for S3 download of PPR VCFs)
- R ≥ 4.3 with packages: `ggplot2`, `reshape2`, `RColorBrewer`, `tidyverse`, `dplyr`
- Reference genome: GRCh38 primary assembly (`Homo_sapiens.GRCh38.dna.primary_assembly.fasta`)
- GIAB v3.6 stratification BED files (with chr prefix)

## Reproducibility

All shell scripts use `set -uo pipefail` for strict error handling. Random seeds are not applicable (deterministic operations). Chromosome renaming is performed after MT exclusion to prevent contig mismatch errors in bcftools operations.
