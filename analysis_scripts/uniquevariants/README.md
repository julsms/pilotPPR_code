# Unique Variants Analysis

Scripts for characterizing pipeline-specific (unique) variants from the PPR vs GATK comparison study. Unique variants are those detected by one variant calling method but not the other on the same sequencing platform.

## Prerequisites

- bcftools, tabix, bgzip
- Ensembl VEP v112 (offline cache, GRCh38)
- Python 3 (pandas, numpy, matplotlib, seaborn, matplotlib-venn, upsetplot)
- R 4.x (VennDiagram, UpSetR, dplyr, readr, ggplot2)
- Nextflow (for annotation pipeline)

## Pipeline Overview

```
1. count_snvs_indels.sh               → Count SNV/INDEL/MNP per VCF
2. annotate_unique_variants.nf        → VEP annotation (main workflow)
   ├── annotate_vcf.nf                → VEP process module
   └── summarize.nf                   → Summary process module
3. extract_consequences.sh            → Extract consequence types (MANE/canonical)
4. extract_dbsnp_status.sh            → Extract dbSNP/ClinVar status
5. extract_gnomad_frequencies.sh      → Extract population allele frequencies
6. extract_mane_canonical_all.sh      → Full CSQ field extraction
7. analyze_unique_variants_regions.sh → Variant distribution across GIAB stratification regions
8. analyse_population_frequencies.py  → Categorize variants by gnomAD/1000G presence
9. analyse_dbsnp_status.py            → Categorize variants by dbSNP presence
10. AF_stats.py                       → Allele frequency distribution statistics
11. analyse_novel_variants.py         → Novel variants analysis (Venn, UpSet; autosomes only)
12. novel_variants_analysis.R         → Novel variants Venn diagrams in R (autosomes only)
```

## Key Outputs

- **Figures**: allele frequency distributions, database presence heatmap, Venn diagrams (novel variants), UpSet plot, consequence type plots
- **Tables**: per-pipeline variant counts, frequency statistics, novel variant lists

## Novel Variant Definition

A variant is classified as **novel** if absent from all three databases:
- No allele frequency > 0 in gnomAD (exomes or genomes)
- No allele frequency > 0 in 1000 Genomes (any population)
- No rsID in dbSNP (via VEP Existing_variation or custom dbSNP annotation)

**Important:** Novel variant analysis is restricted to autosomes (chr1-22) to avoid
false novel calls caused by incomplete sex chromosome annotations in population databases.
Allele frequency and dbSNP/ClinVar analyses include all chromosomes.

## Notes

- Input VCFs are multisample merged files of pipeline-specific unique variants
- VEP annotation uses `--everything --canonical` with custom ClinVar and dbSNP tracks
- Transcript selection priority: MANE Select → Canonical → first transcript
- `analyze_unique_variants_regions.sh` intersects unique variants with GIAB v3.6 stratification BEDs for region-specific counts
