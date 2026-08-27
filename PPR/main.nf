#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Include modules
include { KMC_COUNT } from './modules/kmc'
include { HAPLOTYPE_SAMPLING } from './modules/haplotype_sampling'
include { VG_GIRAFFE } from './modules/vg_giraffe'
include { VG_PACK } from './modules/vg_pack'
include { VG_CALL } from './modules/vg_call'
include { COMPRESS_INDEX_VCF } from './modules/compress_index'

// Include subworkflow
include { SMALLVARIANTS_DEEPVARIANT } from './subworkflows/smallvariants_deepvariant'

// Parameters validation
if (!params.reads) {
    error "Please provide reads with --reads"
}
if (!params.sample_name) {
    error "Please provide sample name with --sample_name"
}
if (!params.outdir) {
    error "Please provide output directory with --outdir"
}

workflow {
    // Create input channel for paired reads
    reads_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)
        .map { sample_id, files -> 
            tuple(params.sample_name, files)
        }
    
    // K-mer counting
    KMC_COUNT(reads_ch)
    
    // Haplotype sampling
    HAPLOTYPE_SAMPLING(KMC_COUNT.out.kff)
    
    // Alignment with vg giraffe
    VG_GIRAFFE(HAPLOTYPE_SAMPLING.out.gbz, reads_ch)
    
    // vg pack
    VG_PACK(VG_GIRAFFE.out.gam, HAPLOTYPE_SAMPLING.out.gbz)
    
    // vg call
    VG_CALL(VG_PACK.out.pack, HAPLOTYPE_SAMPLING.out.gbz)
    
    // Compress and index VCF
    COMPRESS_INDEX_VCF(VG_CALL.out.vcf)
    
    // Optional small variants calling with DeepVariant
    // Pornește doar după ce se termină COMPRESS_INDEX_VCF
    if (params.run_deepvariant) {
        // Creează un trigger bazat pe finalizarea COMPRESS_INDEX_VCF
        trigger_ready = COMPRESS_INDEX_VCF.out.vcf_gz.map { sample_name, vcf_gz -> sample_name }
        
        // Combină datele necesare cu trigger-ul
        gam_ready = VG_GIRAFFE.out.gam
            .combine(trigger_ready)
            .filter { gam_sample, gam_file, trigger_sample -> gam_sample == trigger_sample }
            .map { gam_sample, gam_file, trigger_sample -> tuple(gam_sample, gam_file) }
            
        gbz_ready = HAPLOTYPE_SAMPLING.out.gbz
            .combine(trigger_ready)
            .filter { gbz_sample, gbz_file, trigger_sample -> gbz_sample == trigger_sample }
            .map { gbz_sample, gbz_file, trigger_sample -> tuple(gbz_sample, gbz_file) }
        
        SMALLVARIANTS_DEEPVARIANT(
            gam_ready,
            gbz_ready
        )
    }
}

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    if (params.run_deepvariant) {
        println "Small variants calling with DeepVariant: completed"
    } else {
        println "Small variants calling with DeepVariant: skipped (use --run_deepvariant to enable)"
    }
}
