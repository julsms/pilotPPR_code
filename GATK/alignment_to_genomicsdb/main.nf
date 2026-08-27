#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Parameters 
params.samplename = null    // The --samplename parameter, which must be provided from the command line 
params.fastq1 = null        // The --fastq1 parameter, which must be provided from the command line
params.fastq2 = null        // The --fastq2 parameter, which must be provided from the command line

// Reference genome files 
params.reference = null
params.reference_fai = null
params.reference_dict = null
params.reference_ann = null
params.reference_bwt = null
params.reference_pac = null
params.reference_amb = null
params.reference_sa = null

// Tools
params.picard = null

// Resources for Base Recalibrator
params.resources_knownindels = null
params.resources_knownindels_tbi = null
params.resources_mills = null
params.resources_mills_tbi = null

// Intervals file for parallelization
params.intervals = null

// Output directory
params.outdir = "results"

// Module includes - Only what we need up to GenomicsDB
include { BWAMEM_SORT } from './modules/bwamem.nf'
include { MARKDUPLICATES } from './modules/markduplicates.nf'
include { BUILDBAMINDEX } from './modules/buildbamindex.nf'
include { BASERECALIBRATOR } from './modules/baserecalibrator.nf'
include { GATHERBQSRREPORTS } from './modules/gatherbqsrreports.nf'
include { APPLYBQSR } from './modules/applybqsr.nf'
include { GATHERBAMFILES } from './modules/gatherbamfiles.nf'
include { INDEXBAMBQSR } from './modules/indexbambqsr.nf'
include { HAPLOTYPECALLER } from './modules/haplotypecaller.nf'
include { INDEXGVCF } from './modules/indexgvcf.nf'
include { GENOMICSDBIMPORT } from './modules/genomicsdbimport.update.nf'

workflow {
    
    println "=== PIPELINE CONFIGURATION ==="
    println "Sample Name: ${params.samplename}"
    println "FASTQ1: ${params.fastq1}"
    println "FASTQ2: ${params.fastq2}"
    println "Output Dir: ${params.outdir}"
    println "Intervals: ${params.intervals}"
    println "Reference: ${params.reference}"
    println "=============================="

    
    reference_ch = Channel.value([
        file(params.reference),
        file(params.reference_fai),
        file(params.reference_dict),
        file(params.reference_amb),
        file(params.reference_ann),
        file(params.reference_bwt),
        file(params.reference_pac),
        file(params.reference_sa)
    ])

    
    knownindels_ch = Channel.value([
        file(params.resources_knownindels), 
        file(params.resources_knownindels_tbi)
    ])
    
    mills_ch = Channel.value([
        file(params.resources_mills), 
        file(params.resources_mills_tbi)
    ])

   
    input_ch = Channel.of([params.samplename, file(params.fastq1), file(params.fastq2)])

    
    intervals_ch = Channel.fromPath(params.intervals)
        .splitText()
        .map { it.trim() }
        .filter { it != "" }

    
    BWAMEM_SORT(input_ch, reference_ch)
    
    // Mark duplicates
    MARKDUPLICATES(BWAMEM_SORT.out.bam)
    
    // Build BAM index
    BUILDBAMINDEX(MARKDUPLICATES.out.bam)  
    
    
    bam_bai_combined = MARKDUPLICATES.out.bam
        .join(BUILDBAMINDEX.out.bai)  // [sample_id, bam, bai]
    
    
   
    bam_bai_intervals_ch = bam_bai_combined
        .combine(intervals_ch)  // [sample_id, bam, bai, interval]
        .map { sample_id, bam, bai, interval ->
            // Create a unique task ID but keep interval separate
            def task_id = "${sample_id}.${interval.replaceAll(/[:\-]/, '_')}"
            [task_id, sample_id, bam, bai, interval]
        }
    
    
    BASERECALIBRATOR(
        bam_bai_intervals_ch.map { task_id, sample_id, bam, bai, interval -> 
            [task_id, bam, interval] 
        },
        bam_bai_intervals_ch.map { task_id, sample_id, bam, bai, interval -> 
            [task_id, bai] 
        },
        reference_ch,
        knownindels_ch,
        mills_ch
    )
    
   
    bqsr_tables_ch = BASERECALIBRATOR.out.table
        .map { task_id, table ->
            // Extract original sample_id from task_id
            def original_sample_id = task_id.split("\\.")[0]
            [original_sample_id, table]
        }
        .groupTuple()  // Group by sample_id: [sample_id, [table1, table2, ...]]
   
    GATHERBQSRREPORTS(bqsr_tables_ch)
    
    
    bqsr_table_intervals_ch = GATHERBQSRREPORTS.out.table
        .combine(intervals_ch)  // [sample_id, bqsr_table, interval]
        .map { sample_id, bqsr_table, interval ->
            def task_id = "${sample_id}.${interval.replaceAll(/[:\-]/, '_')}"
            [task_id, sample_id, bqsr_table, interval]
        }
    
    
    bam_bai_intervals_applybqsr = MARKDUPLICATES.out.bam
        .join(BUILDBAMINDEX.out.bai)
        .combine(intervals_ch)
        .map { sample_id, bam, bai, interval ->
            def task_id = "${sample_id}.${interval.replaceAll(/[:\-]/, '_')}"
            [task_id, sample_id, bam, bai, interval]
        }
    
    
    applybqsr_input_ch = bam_bai_intervals_applybqsr
        .join(bqsr_table_intervals_ch)  // Join by task_id
        .map { task_id, sample_id1, bam, bai, interval1, sample_id2, bqsr_table, interval2 ->
            [task_id, sample_id1, bam, bai, bqsr_table, interval1]
        }
    
    
    APPLYBQSR(
        applybqsr_input_ch.map { task_id, sample_id, bam, bai, bqsr_table, interval -> 
            [task_id, bam, interval] 
        },
        applybqsr_input_ch.map { task_id, sample_id, bam, bai, bqsr_table, interval -> 
            [task_id, bai] 
        },
        applybqsr_input_ch.map { task_id, sample_id, bam, bai, bqsr_table, interval -> 
            [task_id, bqsr_table] 
        },
        reference_ch
    )
    
    
    bqsr_bams_ch = APPLYBQSR.out.bam
        .map { task_id, bam ->
            def original_sample_id = task_id.split("\\.")[0]
            [original_sample_id, bam]
        }
        .groupTuple()  // Group by sample_id: [sample_id, [bam1, bam2, ...]]
    
    GATHERBAMFILES(bqsr_bams_ch)
    
    
    INDEXBAMBQSR(GATHERBAMFILES.out.bam)
    
    
    bqsr_bam_ch = GATHERBAMFILES.out.bam      // [sample_id, final_bqsr.bam]
    bqsr_bai_ch = INDEXBAMBQSR.out.bai        // [sample_id, final_bqsr.bam.bai]
    
    haplotypecaller_input_ch = bqsr_bam_ch
        .join(bqsr_bai_ch)  // [sample_id, bam, bai]
        .combine(intervals_ch)  // [sample_id, bam, bai, interval]
        .map { sample_id, bam, bai, interval ->
            [ 
                sample_id,
                bam,
                bai,
                file(params.reference),
                file(params.reference_fai),
                file(params.reference_dict),
                file(params.reference_amb),
                file(params.reference_ann),
                file(params.reference_bwt),
                file(params.reference_pac),
                file(params.reference_sa),
                interval.toString()  
            ]
        }
    
    
    HAPLOTYPECALLER(haplotypecaller_input_ch)
    
    
    // OPTIMIZARE: Prepare GVCFs for indexing - fiecare interval individual
    gvcf_for_indexing = HAPLOTYPECALLER.out.gvcf
        .map { interval, gvcf ->
            def sample_id = params.samplename ?: "sample"
            def unique_id = "${sample_id}_${interval.toString().replaceAll('[^A-Za-z0-9]', '_')}"
            [sample_id, unique_id, gvcf]
        }

    INDEXGVCF(gvcf_for_indexing)

    // OPTIMIZARE: Create channel for GenomicsDBImport care rulează imediat pentru fiecare interval gata
    genomicsdb_input_ch = HAPLOTYPECALLER.out.gvcf  // [interval, gvcf]
        .join(
            INDEXGVCF.out.tbi
                .map { unique_id, tbi -> 
                    // Extract original interval from unique_id
                    def sample_id = params.samplename ?: "sample" 
                    def interval_safe = unique_id.replace("${sample_id}_", "")
                    
                    // Convert safe interval back to original format
                    // Example: "22_1_50818468" -> "22:1-50818468"
                    def parts = interval_safe.split('_')
                    def original_interval
                    if (parts.size() >= 3) {
                        // Format: chr_start_end -> chr:start-end
                        original_interval = "${parts[0]}:${parts[1]}-${parts[2..-1].join('_')}"
                    } else if (parts.size() == 2) {
                        // Format: chr_pos -> chr:pos
                        original_interval = "${parts[0]}:${parts[1]}"
                    } else {
                        // Single part, probably just chromosome
                        original_interval = interval_safe.replace('_', ':')
                    }
                    
                    [original_interval, tbi]
                },
                by: 0  // Join by first element (interval)
            )
        // Result: [interval, gvcf, tbi]

    // OPTIMIZARE: Run GenomicsDBImport - acum rulează imediat pentru fiecare interval gata!
    // Nu mai așteaptă toate intervalele să fie gata
    GENOMICSDBIMPORT(genomicsdb_input_ch, reference_ch)

    // Print completion message pentru fiecare interval individual
    GENOMICSDBIMPORT.out.database.view { interval, database ->
        "✅ GenomicsDB created for interval ${interval}: ${database}"
    }
               
} 
