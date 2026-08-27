#!/usr/bin/env nextflow

nextflow.enable.dsl=2

if (!params.genomicsdb_dir) {
    error "GenomicsDB directory must be specified with --genomicsdb_dir"
}
if (!params.samplename_joint) {
    error "Joint sample name must be specified with --samplename_joint"
}

// GenomicsDB location from pipeline1
params.genomicsdb_dir = null   
params.samplename_joint = null    

// Reference
params.reference = null
params.reference_fai = null
params.reference_dict = null

// VQSR Resource 
params.hapmap = null
params.hapmap_tbi = null
params.omni = null 
params.omni_tbi = null
params.phase1 = null
params.phase1_tbi = null
params.dbsnp = null
params.dbsnp_tbi = null
params.mills = null
params.mills_tbi = null

// Intervals file 
params.intervals = null

// Output directory
params.outdir = "results_joint_calling"

// Annotation BED files dir
params.bed_files_dir = null  

// Module includes 
include { GENOTYPEGVCFS_SPLIT } from './modules/genotypegvcfs_split.nf'
include { MERGE_INTERVAL_VCFS } from './modules/merge_interval_vcfs.nf'
include { VARIANTRECALIBRATOR_SPLIT } from './modules/variantrecalibrator_split.nf'
include { APPLYVQSR_SPLIT } from './modules/applyvqsr_new_split.nf'
include { VCF_ANNOTATION } from './modules/vcf_annotation.nf'

workflow {
    
    println "=== JOINT CALLING PIPELINE CONFIGURATION ==="
    println "GenomicsDB Directory: ${params.genomicsdb_dir}"
    println "Joint Analysis Name: ${params.samplename_joint}"
    println "Output Dir: ${params.outdir}"
    println "Intervals: ${params.intervals}"
    println "Reference: ${params.reference}"
    println "=============================================="

    
    reference_ch = Channel.value([
        file(params.reference),
        file(params.reference_fai),
        file(params.reference_dict)
    ])

    
    intervals_ch = Channel.fromPath(params.intervals)
        .splitText()
        .map { it.trim() }
        .filter { it != "" }

    
    genomicsdb_input_ch = intervals_ch
        .map { interval ->
            def interval_safe = interval.toString().replaceAll('[:\\s]+', '_')
            def genomicsdb_path = "${params.genomicsdb_dir}/${interval_safe}"

            if (!file(genomicsdb_path).exists()) {
                error "GenomicsDB not found for interval ${interval}: ${genomicsdb_path}"
            }     

            [interval, genomicsdb_path]
        }

    println "Found GenomicsDB for ${genomicsdb_input_ch.count()} intervals"

    
    GENOTYPEGVCFS_SPLIT(genomicsdb_input_ch, reference_ch)

    
    interval_vcfs_ch = GENOTYPEGVCFS_SPLIT.out.interval_vcf
        .map { interval, vcf -> vcf }  
        .collect()  

   
    reference_merge_ch = Channel.value([
        file(params.reference),
        file(params.reference_fai), 
        file(params.reference_dict)
    ])

    
    MERGE_INTERVAL_VCFS(interval_vcfs_ch, reference_merge_ch)

    // VQSR resource channels
    hapmap_ch = Channel.value([file(params.hapmap), file(params.hapmap_tbi)])
    omni_ch = Channel.value([file(params.omni), file(params.omni_tbi)])
    phase1_ch = Channel.value([file(params.phase1), file(params.phase1_tbi)])
    dbsnp_vqsr_ch = Channel.value([file(params.dbsnp), file(params.dbsnp_tbi)])
    mills_vqsr_ch = Channel.value([file(params.mills), file(params.mills_tbi)])

    // Reference channel for VQSR
    reference_vqsr_ch = Channel.value([
        file(params.reference),
        file(params.reference_fai),
        file(params.reference_dict)
    ])

    // VariantRecalibrator
    VARIANTRECALIBRATOR_SPLIT(
        MERGE_INTERVAL_VCFS.out.merged_vcf,
        MERGE_INTERVAL_VCFS.out.merged_vcf_index,
        reference_vqsr_ch,
        hapmap_ch,
        omni_ch,
        phase1_ch,
        dbsnp_vqsr_ch,
        mills_vqsr_ch
    )

    // Apply VQSR - Saves recalibrated VCF to results
    APPLYVQSR_SPLIT(
        MERGE_INTERVAL_VCFS.out.merged_vcf,
        MERGE_INTERVAL_VCFS.out.merged_vcf_index,
        reference_vqsr_ch,
        VARIANTRECALIBRATOR_SPLIT.out.snps_recal,
        VARIANTRECALIBRATOR_SPLIT.out.snps_tranches,
        VARIANTRECALIBRATOR_SPLIT.out.snps_recal_idx,
        VARIANTRECALIBRATOR_SPLIT.out.indels_recal,
        VARIANTRECALIBRATOR_SPLIT.out.indels_tranches,
        VARIANTRECALIBRATOR_SPLIT.out.indels_recal_idx,
        VARIANTRECALIBRATOR_SPLIT.out.snps_vcf,
        VARIANTRECALIBRATOR_SPLIT.out.snps_vcf_index,
        VARIANTRECALIBRATOR_SPLIT.out.indels_vcf,
        VARIANTRECALIBRATOR_SPLIT.out.indels_vcf_index
    )

    // Prepare for annotation
    bed_files_ch = Channel.fromPath("${params.bed_files_dir}/*.bed*")

    // Annotate the final VCF
    VCF_ANNOTATION(
        APPLYVQSR_SPLIT.out.merged_recalibrated,
        bed_files_ch.collect()
    )

    // Display final results
    APPLYVQSR_SPLIT.out.merged_recalibrated.view { vcf, idx ->
        "✅ Recalibrated VCF saved: ${vcf}"
    }
    
    VCF_ANNOTATION.out.annotated_vcf.view { vcf ->
        "✅ Joint calling completed! Final annotated VCF: ${vcf}"
    }
}
