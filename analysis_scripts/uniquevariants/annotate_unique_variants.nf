#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.base_dir = "/mnt/genomics/pilot_PPR/uniquevariants"
params.outdir = "/mnt/genomics/pilot_PPR/uniquevariants/annotated"
params.vep_dir = "/mnt/genomics/tools/ensembl-vep"
params.vep_cache = "/mnt/genomics/tools/vep_data"
params.reference = "/mnt/genomics/illumina/ref_chr/Homo_sapiens.GRCh38.dna.primary_assembly.fasta"
params.clinvar = "/mnt/genomics/tools/vep_data/annotations/clinvar.vcf.gz"
params.dbsnp = "/mnt/genomics/tools/vep_data/annotations/dbSNP_latest_GRCh38.vcf.gz"
params.assembly = "GRCh38"
params.species = "homo_sapiens"
params.fork = 4

include { ANNOTATE_VCF } from './modules/annotate_vcf.nf'
include { GENERATE_CATEGORY_SUMMARY; GENERATE_GLOBAL_SUMMARY } from './modules/summarize.nf'

workflow {
    def categories = [
        [category: 'GATK_illumina', dir: 'uniqueGATKillumina', pattern: '_unique_GATK_illumina.vcf.gz'],
        [category: 'GATK_aviti', dir: 'uniqueGATKaviti', pattern: '_unique_GATK_aviti.vcf.gz'],
        [category: 'PPR_illumina', dir: 'uniquePPRillumina', pattern: '_unique_PPR_illumina.vcf.gz'],
        [category: 'PPR_aviti', dir: 'uniquePPRaviti', pattern: '_unique_PPR_aviti.vcf.gz']
    ]
    
    vcf_ch = Channel.fromList(categories)
        .flatMap { cat ->
            def vcf_dir = file("${params.base_dir}/${cat.dir}")
            
            // FIX: Use file() with glob pattern instead of .matches()
            def glob_pattern = "*${cat.pattern}"
            def vcf_files = file("${vcf_dir}/${glob_pattern}")
            
            // Ensure it's a list
            def file_list = vcf_files instanceof List ? vcf_files : [vcf_files]
            
            // Filter and map
            file_list
                .findAll { it.exists() && it.name.endsWith('.vcf.gz') }
                .collect { vcf_file ->
                    def vcf_idx = file("${vcf_file}.tbi")
                    def sample_name = vcf_file.name.replaceAll(/\.vcf\.gz$/, '')
                    [cat.category, sample_name, vcf_file, vcf_idx]
                }
        }
    
    ANNOTATE_VCF(vcf_ch)
    
    category_summaries = ANNOTATE_VCF.out.annotated_vcf
        .map { category, sample, vcf, idx -> [category, vcf] }
        .groupTuple()
    
    GENERATE_CATEGORY_SUMMARY(category_summaries)
    
    all_vcfs = ANNOTATE_VCF.out.annotated_vcf
        .map { category, sample, vcf, idx -> vcf }
        .collect()
    
    GENERATE_GLOBAL_SUMMARY(all_vcfs, GENERATE_CATEGORY_SUMMARY.out.tsv.collect())
}

workflow.onComplete {
    println "\n${'='*70}"
    println "Completed: ${workflow.complete} | Duration: ${workflow.duration}"
    println "Success: ${workflow.success} | Output: ${params.outdir}"
    println "${'='*70}"
}
