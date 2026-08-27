process GATHERBAMFILES {
    publishDir "${params.outdir}/${sample_id}/bqsr", mode: 'symlink'

    input:
    tuple val(sample_id), path(bam_files)

    output:
    tuple val(sample_id), path("${sample_id}_bqsr.bam"), emit: bam
    path "versions.yml", emit: versions

    script:
    def tmp_dir = "/mnt/genomics/tmp"
    """
    echo "Merging BAM files for sample: ${sample_id}"
    echo "Input BAMs: ${bam_files}"

    
    samtools merge \\
        -@ ${task.cpus} \\
        -o ${sample_id}_bqsr.bam \\
        ${bam_files.join(' ')}

    
    rm ${bam_files.join(' ')}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version 2>&1 | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${sample_id}_bqsr.bam
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: 1.17
    END_VERSIONS
    """
}
