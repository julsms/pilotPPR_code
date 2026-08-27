process BUILDBAMINDEX {

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}_markdup.bai"), emit: bai
    path "versions.yml", emit: versions

    script:
    def tmp_dir = "/mnt/genomics/tmp"
    """
    java \\
        -XX:+UseParallelGC \\
        -XX:ParallelGCThreads=${task.cpus} \\
        -Xmx${task.memory.toGiga()-2}g \\
        -Djava.io.tmpdir=${tmp_dir} \\
        -jar ${params.picard} BuildBamIndex \\
        I=${bam} \\
        VALIDATION_STRINGENCY=LENIENT \\
        TMP_DIR=${tmp_dir}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(java -jar ${params.picard} BuildBamIndex --version 2>&1 | grep -oP 'Version:\\s*\\K[0-9.]+' || echo "unknown")
    END_VERSIONS
    """

    stub:
    """
    touch ${sample_id}_markdup.bai
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: 3.0.0
    END_VERSIONS
    """
}
