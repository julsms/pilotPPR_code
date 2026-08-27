# Personalized Pangenome References (PPR) — Small Variant Calling Pipeline

Nextflow DSL2 implementation of the personalized pangenome reference (PPR) approach for short-read small variant calling, following [Sirén et al.](https://www.nature.com/articles/s41588-024-01824-y). Constructs a sample-specific subgraph from the HPRC v1.1 pangenome, aligns reads to it, and calls variants with DeepVariant on the surjected alignments.

## Pipeline steps

1. **KMC_COUNT** — k-mer counting from the sample's raw reads (KMC)
2. **HAPLOTYPE_SAMPLING** — imputes which haplotypes from the HPRC v1.1 pangenome (`hprc-v1.1-mc-grch38.gbz`) are present in the sample, based on k-mer counts, and constructs a personalized subgraph
3. **VG_GIRAFFE** — aligns reads to the personalized subgraph
4. **VG_PACK** / **VG_CALL** — graph-based coverage packing and variant calling on the subgraph
5. **VG_VIEW** / **VG_INDEX_XG** / **VG_SURJECT** — converts the personalized graph to linear (XG) form and surjects the graph alignment (GAM) to a linear BAM
6. **BAM_PREPROCESSING** — sorts, fixes header, indexes the surjected BAM
7. **DEEPVARIANT** — calls small variants (SNVs + indels) from the surjected BAM, GPU-accelerated (native Google DeepVariant container, WGS model)
8. **COMPRESS_INDEX_VCF** — bgzip + tabix the output VCF

Steps 1–3 are orchestrated directly in `main.nf`; steps 5–7 are grouped in the `SMALLVARIANTS_DEEPVARIANT` subworkflow.

## Scripts

| File | What it does |
|------|--------------|
| `main.nf` | Top-level workflow: k-mer counting → haplotype sampling → graph alignment, then calls the small-variants subworkflow |
| `nextflow.config` | Reference genome, pangenome graph paths, per-process resource limits (cpus/memory/time), Docker settings |
| `modules/kmc.nf` | K-mer counting (KMC) |
| `modules/haplotype_sampling.nf` | Personalized subgraph construction from HPRC v1.1 |
| `modules/vg_giraffe.nf` | Read alignment to the personalized subgraph |
| `modules/vg_pack.nf`, `modules/vg_call.nf` | Graph-based coverage packing and calling |
| `modules/vg_view.nf`, `modules/vg_index_xg.nf`, `modules/vg_surject.nf` | GBZ→VG conversion, XG indexing, surjection of graph alignment to linear BAM |
| `modules/bam_preprocessing.nf` | BAM sort/header-fix/index |
| `modules/deepvariant.nf` | DeepVariant small variant calling (GPU) |
| `modules/compress_index.nf` | VCF bgzip + tabix |
| `subworkflows/smallvariants_deepvariant.nf` | Groups surjection → BAM preprocessing → DeepVariant into one subworkflow |

## Key parameters (`nextflow.config`)

- `graph_dir` / `graph_gbz` / `graph_hapl` — HPRC v1.1 default pangenome graph
- `reference_dir` / `reference_fasta` — GRCh38 primary assembly (used for DeepVariant's reference-aware calling)
- `max_cpus`, `max_memory`, `max_time` — resource ceilings; per-process overrides under `process { withName: ... }`

## Dependencies

Nextflow (DSL2), Docker, vg toolkit (v1.52.0), KMC, DeepVariant (`google/deepvariant:1.9.0-gpu`), bcftools/htslib (bgzip, tabix)

## Scope note

This repository covers the small-variant calling pipeline as described in the manuscript's Methods. It does not include structural-variant-oriented extensions (vg deconstruct, vcfbub filtering, PanGenie genotyping) or the optional GATK-style VQSR post-processing path that exist in the working development environment — these were not used to produce the results reported in the manuscript (structural variant calling was explicitly out of scope; VQSR is applied to the GATK pipeline only, not PPR).

## Reference

Sirén J, et al. Personalized pangenome references. *Nat Methods* 21, 2017–2023 (2024).
