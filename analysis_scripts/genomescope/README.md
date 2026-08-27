# GenomeScope 2.0 K-mer Spectrum Analysis

K-mer spectrum analysis for platform comparison (Illumina vs Element Aviti) using coverage-matched subsampling.

## Overview

For each sample sequenced on both platforms, reads are subsampled to equal coverage (minimum of the pair) before k-mer counting, ensuring unbiased platform comparison. The pipeline runs FastK for k-mer counting, exports histograms with Histex, and fits GenomeScope 2.0 models.

## Pipeline steps

1. **Subsampling** (conditional) — Rasusa downsamples the higher-coverage platform to match its pair (seed=42)
2. **K-mer counting** — FastK (k=21)
3. **Histogram export** — Histex with `-G` flag (cumulative last bin for accurate genome size estimation)
4. **Model fitting** — GenomeScope 2.0 (diploid model)

## Files

| File | Description |
|------|-------------|
| `main.nf` | Nextflow DSL2 pipeline |
| `nextflow.config` | Pipeline configuration |
| `run_illumina.sh` | Wrapper: runs all Illumina samples with coverage-matched subsampling |
| `run_aviti.sh` | Wrapper: runs all Aviti samples with coverage-matched subsampling |

## Usage

```bash
# Run all Illumina samples
./run_illumina.sh

# Run all Aviti samples
./run_aviti.sh
```

The wrapper scripts read per-sample coverage from a TSV file, determine the subsampling target for each paired sample, and launch the Nextflow pipeline accordingly.

## Tools and versions

- Nextflow ≥24.04
- FastK v1.2
- Histex (from FastK suite)
- Rasusa (Hall 2022)
- GenomeScope 2.0 (Ranallo-Benavidez et al., 2020)

## Output

Per-sample directories containing:
- `*_k21.histo` — k-mer frequency histogram
- `genomescope2_output/` — GenomeScope model fit, summary statistics, and plots
