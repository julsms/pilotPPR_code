# pilotPPR — Code and Analysis Scripts

Code accompanying the manuscript comparing GATK and personalized pangenome reference (PPR) small variant calling on Illumina and Element Aviti short-read WGS data.

## Structure

- **`GATK/`** — GATK best-practices small variant calling pipeline (Nextflow DSL2), in two stages: per-sample alignment through GenomicsDB import, then cohort-level joint genotyping and VQSR. See `GATK/README.md`.
- **`PPR/`** — Personalized pangenome reference (PPR) small variant calling pipeline (Nextflow DSL2): k-mer-based haplotype sampling, graph alignment, and DeepVariant calling. See `PPR/README.md`.
- **`analysis_scripts/`** — Downstream analysis and figure-generation scripts used in the manuscript: GIAB benchmarking statistics, genotype concordance, coverage correlation, GenomeScope k-mer spectrum analysis, Merfin k-mer variant validation, and pipeline-specific unique variant characterization. See `analysis_scripts/README.md` and the README in each subfolder.

Each pipeline/script folder has its own README with methodological detail; this file is just the map between them.
