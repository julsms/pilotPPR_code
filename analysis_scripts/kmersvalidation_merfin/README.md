# Merfin K-mer Variant Validation

K-mer-based validation of variant calls from GATK and PPR pipelines using [Merfin](https://github.com/arangrhie/merfin), providing alignment-free quality assessment.

## Overview

For each sample, a Meryl k-mer database is built from raw reads (k=31, singletons excluded). Merfin then filters each variant set, retaining only variants supported by k-mer evidence. Validation rate = retained variants / total variants × 100.

Two levels of analysis, computed in two dependent steps:

1. **All variants** — Merfin applied directly to the complete GATK and PPR call sets per sample (`merfinvalidation_illumina.sh`). This is the expensive step (builds the Meryl database, runs `merfin -filter`) and its output VCFs are reused below.
2. **Pipeline-specific unique variants** — rather than re-running Merfin -filter on the unique VCFs (costly, and the unique-variant VCFs were later revised), validation rate is instead computed by intersecting the unique variant VCFs with the already Merfin-validated complete set from step 1, via `bcftools isec` (`uniquevariants_validationrate_normalized.sh`). Both sides are normalized first (`bcftools norm -m-both`) so indel representation matches across the two VCFs.

Additionally, validation rates were stratified across GIAB v3.6 genomic regions. **Only the stratified analysis on all variants was reported in the manuscript**; stratified unique variant rates are included in the script for completeness but were not used due to normalization artifacts (see Limitations below).

## Scripts

| Script | What it does |
|--------|-------------|
| `merfinvalidation_illumina.sh` | Step 1: builds the Meryl DB and runs Merfin `-filter` on the complete (ALL variants) GATK and PPR call sets per sample. Shown for a single Illumina sample as a worked example; run per-sample for the full cohort and for Aviti with platform-specific paths substituted. Its output (`<sample>_gatk_all_validated.filter.vcf.gz`, `<sample>_ppr_all_validated.filter.vcf.gz`) is the input to step 2. |
| `uniquevariants_validationrate_normalized.sh` | Step 2: computes unique-variant validation rates by normalizing (`bcftools norm -m-both`) and intersecting (`bcftools isec`) the unique VCFs against the step-1 Merfin-validated complete set, instead of re-running Merfin on the unique VCFs directly — avoids positional mismatches for indels and the cost of a second Merfin run each time the unique-variant VCFs are revised. |
| `stratified_validation.sh` | Computes per-region validation rates (10 GIAB stratifications) for both all and unique variants, both platforms. |

## Key details

- **Meryl DB**: k=31, frequency >1 (singletons removed to exclude sequencing errors)
- **Merfin mode**: `-filter` (retains k-mer-supported variants)
- **Reference**: GRCh38 primary assembly
- **Unique variant validation**: uses `bcftools isec` to intersect unique VCFs with the validated complete set, rather than running Merfin independently on each subset

## Limitations

The unique variant validation rates computed by the stratified script can exceed 100% due to VCF normalization differences between the unique variant files and the Merfin-validated set (left-alignment and multiallelic splitting cause positional shifts for some indels). This limitation is acknowledged in the manuscript.

## Output

- `final_results/` — per-sample CSVs with validation counts and rates
- `vcf_results/<sample>/` — input and Merfin-filtered VCFs (compressed + indexed)
- `validationrates_stratified/` — stratified validation CSVs (per region × pipeline × platform)

## Dependencies

Meryl, Merfin, bcftools, htslib (bgzip/tabix), AWS CLI, Python 3

## Reference

Formenti G, et al. Merfin: improved variant filtering, assembly evaluation and polishing via k-mer validation. *Nat Methods* 19, 696–704 (2022).
