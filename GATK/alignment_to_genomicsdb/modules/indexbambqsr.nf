process INDEXBAMBQSR {
    publishDir "${params.outdir}/${sample_id}/bqsr", mode: 'symlink'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${bam}.bai"), emit: bai
    path "versions.yml", emit: versions

    script:
    """
    echo "Building BAM index for BQSR BAM: ${bam}"
    
    # Verify BAM file is sorted before indexing
    samtools quickcheck ${bam}
    
    # Use samtools index with threading
    samtools index -@ ${task.cpus} ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version 2>&1 | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${bam}.bai
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: 1.17
    END_VERSIONS
    """
}
