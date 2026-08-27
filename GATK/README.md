# GATK Small Variant Calling Pipeline

Nextflow DSL2 implementation of the GATK best-practices small variant calling pipeline (CPU, no Docker/GPU), run in two stages. Stage 1 is per-sample, parallelized across genomic intervals; stage 2 joins the per-sample GenomicsDB stores into one cohort and jointly genotypes them.

## Stage 1 — `alignment_to_genomicsdb/`

Per-sample: alignment through per-interval GVCF calling and GenomicsDB import.

1. **BWAMEM_SORT** — align reads to GRCh38 (BWA-MEM), sort
2. **MARKDUPLICATES** / **BUILDBAMINDEX** — mark duplicates (Picard), index
3. **BASERECALIBRATOR** / **GATHERBQSRREPORTS** / **APPLYBQSR** — base quality score recalibration, parallelized per genomic interval, gathered back per sample
4. **GATHERBAMFILES** / **INDEXBAMBQSR** — merge per-interval BQSR'd BAMs into the final per-sample BAM, index
5. **HAPLOTYPECALLER** — GVCF-mode variant calling, per interval
6. **INDEXGVCF** — index each per-interval GVCF
7. **GENOMICSDBIMPORT** — import per-interval GVCFs into a GenomicsDB store (one per interval), ready for joint genotyping

Run per sample; in this study, intervals parallelized the workflow across 66 genomic intervals (see `nextflow.config` / `--intervals`).

**Illumina vs. Aviti:** as shown, `main.nf` is written for Illumina data. For Aviti samples, the same workflow is run except the BQSR step is skipped (`BASERECALIBRATOR` / `GATHERBQSRREPORTS` / `APPLYBQSR` omitted, `MARKDUPLICATES` output goes directly into `HAPLOTYPECALLER`), since Aviti's reported base quality scores have been shown to be well-calibrated empirically and BQSR was not applied to Aviti data in this study (see Methods).

## Stage 2 — `joint_genotyping/`

Cohort-level: joins all samples' GenomicsDB stores (`--genomicsdb_dir`, produced by stage 1) and jointly genotypes, recalibrates, and annotates.

1. **GENOTYPEGVCFS_SPLIT** — joint genotyping from the GenomicsDB store, per interval
2. **MERGE_INTERVAL_VCFS** — merge per-interval VCFs into one cohort VCF
3. **VARIANTRECALIBRATOR_SPLIT** — build the VQSR recalibration model (SNPs and indels) using HapMap, Omni, 1000G phase1, Mills, and dbSNP resources
4. **APPLYVQSR_SPLIT** — apply VQSR filtering to the cohort VCF
5. **VCF_ANNOTATION** — functional annotation of the final callset

## Why two directories

These are two separate Nextflow runs, not one workflow: stage 1 runs independently per sample, and stage 2 only starts once every sample's GenomicsDB store from stage 1 is available, since it genotypes the full cohort together. The `jointgenotyping.nf` params (`--genomicsdb_dir`, `--samplename_joint`) are explicitly documented in the script as taking their input from stage 1's output.

## Scope note

This covers the CPU-only pipeline as run for this study and described in the manuscript's Methods (BWA-MEM, Picard duplicate marking, BaseRecalibrator/ApplyBQSR, HaplotypeCaller in GVCF mode across 66 intervals, GenotypeGVCFs, and VQSR using HapMap/1000 Genomes/dbSNP resources). A Docker/GPU-accelerated variant of stage 1 also exists in the development environment (GPU HaplotypeCaller/GenomicsDBImport via NVIDIA Parabricks) but was not used to produce the results reported in the manuscript and is not included here.

## Dependencies

Nextflow (DSL2), GATK4 (BWA-MEM, Picard, BaseRecalibrator/ApplyBQSR, HaplotypeCaller, GenomicsDBImport, GenotypeGVCFs, VariantRecalibrator/ApplyVQSR), a variant annotation tool (VCF_ANNOTATION module)
