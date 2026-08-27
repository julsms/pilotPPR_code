process INDEXGVCF {
    tag "${sample_id}"
    publishDir "${params.outdir}/${sample_id}/gvcf", mode: 'symlink'

    input:
    tuple val(sample_id), val(unique_id), path(gvcf)

    output:
    tuple val(unique_id), path("${gvcf}.tbi"), emit: tbi
    path "versions.yml", emit: versions

    script:
    """
    echo "Indexing GVCF: ${gvcf} for sample: ${sample_id}"
    
    gatk --java-options \\
        "-Xms${task.memory.toGiga()}g \\
        -Xmx${task.memory.toGiga()}g \\
        -Djava.io.tmpdir=/mnt/genomics/tmp" \\
        IndexFeatureFile \\
        -I ${gvcf}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | grep -oP 'Version:\\s*\\K[0-9.]+' || echo "unknown")
    END_VERSIONS
    """

    stub:
    """
    touch ${gvcf}.tbi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: 4.4.0.0
    END_VERSIONS
    """
}
