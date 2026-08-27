#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Parameters
params.query_vcf = null
params.gold_vcf = null
params.reference = null
params.benchmark_bed = null  // BED high-confidence din benchmark (v4.2.1 sau v5.0)
params.bed_dir = null
params.outdir = null

// BED files list
params.bed_files = [
    "TRHP.bed", "LM.bed", "OD.bed", "AD.bed", "25GC65.bed", "15GC85.bed",
    "REFSEQ.bed", "BB.bed", "ACMG.bed", "OMIM.bed"
]

// Include processes
include { FILTER_BENCHMARK_REGIONS } from './modules/filter_benchmark.nf'
include { NORMALIZE_VCF } from './modules/normalize_vcf.nf'
include { EXTRACT_REGIONS } from './modules/extract_regions.nf'

workflow {
    
    println "=== VCF NORMALIZATION AND REGION EXTRACTION PIPELINE ==="
    println "Query VCF: ${params.query_vcf}"
    println "Gold VCF: ${params.gold_vcf}"
    println "Reference: ${params.reference}"
    println "Benchmark BED: ${params.benchmark_bed}"
    println "BED Directory: ${params.bed_dir}"
    println "Output Directory: ${params.outdir}"
    println "Number of BED files: ${params.bed_files.size()}"
    println "============================================================"

    // Create input channels
    query_vcf_ch = Channel.value([file(params.query_vcf), "query"])
    gold_vcf_ch = Channel.value([file(params.gold_vcf), "gold"])
    
    // Combine both VCFs for parallel processing
    vcf_input_ch = query_vcf_ch.mix(gold_vcf_ch)
    
    // Step 0: Filter to benchmark regions (NEW STEP - pre-normalization)
    FILTER_BENCHMARK_REGIONS(
        vcf_input_ch,
        file(params.benchmark_bed)
    )
    
    // Step 1: Normalize VCFs in parallel (uses filtered VCFs)
    NORMALIZE_VCF(
        FILTER_BENCHMARK_REGIONS.out.filtered_vcf,
        file(params.reference)
    )
    
    // Step 2: Create BED files channel
    bed_files_ch = Channel.fromList(params.bed_files)
        .map { bed_file -> 
            [bed_file, file("${params.bed_dir}/${bed_file}")]
        }
    
    // Step 3: Extract regions in parallel for all combinations
    // Combine normalized VCFs with all BED files
    extraction_input_ch = NORMALIZE_VCF.out.normalized_vcf
        .combine(bed_files_ch)  // [vcf_file, vcf_index, sample_type, bed_name, bed_file]
    
    EXTRACT_REGIONS(extraction_input_ch)
    
    // Summary
    EXTRACT_REGIONS.out.region_vcf
        .collect()
        .view { vcf_list ->
            "Pipeline completed successfully!"
            "Generated ${vcf_list.size()} region-specific VCF files"
            "All files saved to: ${params.outdir}"
        }
}
