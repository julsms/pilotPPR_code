process MARKDUPLICATES {
    publishDir "${params.outdir}/${sample_id}/markduplicates", mode: 'symlink', pattern: "*_markdup.txt"
    

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}_markdup.bam"), emit: bam
    tuple val(sample_id), path("${sample_id}_markdup.txt"), emit: metrics
    path "versions.yml", emit: versions

    script:
    def tmp_dir = "/mnt/genomics/tmp"
    """
    java \\
        -XX:+UseParallelGC \\
        -XX:ParallelGCThreads=${task.cpus} \\
        -Xms${task.memory.toGiga()-2}g \\
        -Xmx${task.memory.toGiga()-2}g \\
        -Djava.io.tmpdir=${tmp_dir} \\
        -jar ${params.picard} MarkDuplicates \\
        I=${bam} \\
        O=${sample_id}_markdup.bam \\
        M=${sample_id}_markdup.txt \\
        VALIDATION_STRINGENCY=LENIENT \\
        ASSUME_SORT_ORDER=coordinate \\
        CREATE_INDEX=false \\
        TMP_DIR=${tmp_dir}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(java -jar ${params.picard} MarkDuplicates --version 2>&1 | grep -oP 'Version:\\s*\\K[0-9.]+' || echo "unknown")
    END_VERSIONS
    """

    stub:
    """
    touch ${sample_id}_markdup.bam
    touch ${sample_id}_markdup.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: 3.0.0
    END_VERSIONS
    """
}


