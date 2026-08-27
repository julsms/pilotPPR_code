# pilotPPRpaper_scripts

Analysis scripts for:

**Improved small variant detection using personalized pangenome references for short-read sequencing data with k-mer-based validation**

## Repository Structure

| Directory | Analysis |
|-----------|----------|
| `benchmark/` | GIAB HG002 benchmarking with hap.py across genomic stratifications |
| `concordance/` | Cross-platform and cross-method genotype concordance |
| `coverage_correlations/` | Coverage vs. validation rate and concordance correlations |
| `genomescope/` | K-mer spectrum analysis with GenomeScope 2.0 |
| `kmersvalidation_merfin/` | K-mer-based variant validation with Merfin |
| `uniquevariants/` | Pipeline-specific variant characterization (annotation, frequencies, novel variants) |

## Overview

Comparison of GATK and personalized pangenome reference (PPR) pipelines across 23 human samples sequenced on Illumina and Element Aviti platforms. Each subdirectory contains its own README with script descriptions and execution order.
