# Variant Calling Benchmark — hap.py Pipeline

Scripts for benchmarking GATK and PPR variant calling pipelines against the
GIAB HG002 v5.0 truth set, across ten genomic stratification regions.

## Overview

```
Query VCF + Gold VCF
        │
        ▼
Nextflow: filter to benchmark regions → normalize → extract per BED region
        │
        ▼
hap.py (vcfeval engine) — per region
        │
        ▼
R: plots + statistical analysis
```

## Scripts

### Nextflow pipeline — VCF preparation

| File | Description |
|------|-------------|
| `main.nf` | DSL2 pipeline: filters both VCFs to GIAB high-confidence regions, normalizes, then extracts each of the 10 stratification BED regions in parallel |
| `nextflow.config` | Resource config: local executor, 40 CPUs, 40 GB RAM; process-level overrides for NORMALIZE_VCF and EXTRACT_REGIONS |

**Modules (in `./modules/`):**

- `filter_benchmark.nf` — restricts VCFs to GIAB high-confidence BED
- `normalize_vcf.nf` — bcftools norm against GRCh38
- `extract_regions.nf` — bcftools view per stratification BED

**Parameters:**

- Gold standard: GIAB HG002 GRCh38 v5.0q (`HG002_GRCh38_v5.0q_smvar.vcf.gz`)
- Reference: GRCh38 primary assembly
- Stratification BEDs (GIAB v3.6): TRHP, LM, OD, AD, 25GC65, 15GC85, REFSEQ, BB, ACMG, OMIM

**Output per VCF × region:** `<outdir>/<region>/query_<region>.vcf` and `gold_<region>.vcf`

---

### `run_happy.sh` — hap.py execution

Loops over the 10 stratification regions and runs `hap.py` with the `vcfeval`
engine for each, comparing the query VCF against the gold standard.

**Input:** per-region VCFs produced by the Nextflow pipeline  
**Output per region:** `<happy_outdir>/<region>/<region>.summary.csv` (plus extended CSV and other hap.py outputs)  
**Key hap.py options:** `--engine=vcfeval --threads 8`

---

### `compare_happy_4pipelines.R` — performance plots

Takes four pipeline result directories as command-line arguments and produces
performance comparison plots across all ten genomic stratifications.

**Usage:**

```bash
Rscript compare_happy_4pipelines.R \
    IlluminaGATK /path/to/illumina_gatk/happy \
    IlluminaPPR  /path/to/illumina_ppr/happy  \
    AvitiGATK    /path/to/aviti_gatk/happy    \
    AvitiPPR     /path/to/aviti_ppr/happy
```

**What it does:**

1. Reads all `<region>.summary.csv` files from each pipeline directory
2. Plots Precision, Recall, and F1-score as line plots across regions (solid = ALL variants, dashed = PASS)
3. Produces separate figures for SNPs and INDELs
4. Generates a wide-format comparison table with per-pipeline scores and pairwise differences

**Output files:**

- `HG002_SNP_comparison_<pipelines>.pdf/png`
- `HG002_INDEL_comparison_<pipelines>.pdf/png`
- `comparison_<pipelines>_happy.csv`

**Dependencies:** tidyverse, patchwork

---

### `benchmark_satisticalanalysis_friedman_wilcoxon.R` — statistical analysis

Reads `all_benchmark_data.csv` (aggregated hap.py summary data across all
pipelines and regions) and performs a non-parametric statistical comparison.

**Rationale:** The 10 genomic stratification regions are evaluated on the same
sample (HG002) and the same sequencing data — they are not independent
replicates. Some regions overlap (e.g., AD is the union of LM + TRHP + OD).
Standard ANOVA assumes independence; the Friedman test does not.

**Input:** `all_benchmark_data.csv` in the working directory  
**Output directory:** `benchmark_stats_results/`

**Analysis steps:**

1. **Friedman test** — non-parametric repeated-measures test treating genomic
   regions as blocks and pipelines as treatments. Reports Kendall's W effect
   size. This is the primary statistical test.
2. **Paired Wilcoxon signed-rank tests** — post-hoc pairwise comparisons:
   - Method effect: PPR vs GATK (within each platform)
   - Platform effect: Aviti vs Illumina (within each method)
   - P-values corrected with Holm method
   - Effect size: matched-pairs rank-biserial r
3. **Two-way ANOVA** — Platform × Method interaction model, kept as
   supplementary analysis for comparison with parametric approach. Reports η².
4. **Descriptive statistics** — mean, SD, median, min, max per pipeline ×
   variant type × metric.
5. **Concordance check** — compares ANOVA vs Friedman/Wilcoxon conclusions to
   verify robustness.
6. **Paper-ready tables** — formatted output for direct use in manuscript.

**Output files:**

| File | Content |
|------|---------|
| `friedman_test_results.csv` | Primary non-parametric test (χ², df, p, Kendall W) |
| `wilcoxon_paired_results.csv` | Post-hoc pairwise comparisons (W, p, Holm-corrected p, effect size r, median differences) |
| `anova_results_reference.csv` | Supplementary ANOVA (F, p, η²) |
| `descriptive_statistics.csv` | Summary statistics per pipeline × metric × region |
| `paper_friedman_table.csv` | Friedman results formatted for manuscript |
| `paper_wilcoxon_table.csv` | Wilcoxon results formatted for manuscript |

**Dependencies:** tidyverse, rstatix, broom

---

## Input data format

`all_benchmark_data.csv` should have columns:

```
Pipeline,Region,Type,Filter,F1_Score,Precision,Recall
```

- **Pipeline:** one of `IlluminaGATK`, `IlluminaPPR`, `AvitiGATK`, `AvitiPPR`
- **Region:** one of the 10 stratification codes (see below)
- **Type:** `SNP` or `INDEL`
- **Filter:** `ALL` or `PASS`

Expected: 10 regions × 4 pipelines × 2 types × 2 filters = 160 rows.

---

## Genomic stratification regions

| Code | Description |
|------|-------------|
| BB | Benchmark BED (GIAB v5.0q high-confidence regions) |
| REFSEQ | RefSeq coding regions |
| ACMG | ACMG Secondary Findings v3.3 genes |
| OMIM | OMIM disease genes |
| 25GC65 | Moderate GC content (below 25% or above 65%) |
| 15GC85 | Extreme GC content (below 15% or above 85%) |
| AD | All difficult regions (union of LM + TRHP + OD) |
| TRHP | Tandem repeats and homopolymers |
| OD | Other difficult regions |
| LM | Low-mappability regions |

Stratification BED files from GIAB v3.6 (Dwarshuis et al., 2024).

---

## Reference

Krusche P, et al. Best practices for benchmarking germline small-variant calls
in human genomes. *Nat Biotechnol* 37, 555–560 (2019).
