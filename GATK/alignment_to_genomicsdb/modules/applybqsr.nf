process APPLYBQSR {
    
    input:
    tuple val(task_id), path(bam), val(interval)
    tuple val(task_id), path(bai)
    tuple val(task_id), path(bqsr_table)
    tuple path(reference), path(reference_fai), path(reference_dict), path(reference_amb), path(reference_ann), path(reference_bwt), path(reference_pac), path(reference_sa)

    output:
    tuple val(task_id), path("${task_id}_bqsr.bam"), emit: bam
    path "versions.yml", emit: versions

    script:
    def tmp_dir = "/mnt/genomics/tmp"

    """
    gatk --java-options \\
        "-Xms${task.memory.toGiga()}g \\
        -Xmx${task.memory.toGiga()}g \\
        -XX:+UseParallelGC \\
        -XX:ParallelGCThreads=${task.cpus} \\
        -Djava.io.tmpdir=${tmp_dir}" \\
        ApplyBQSR \\
        -I ${bam} \\
        -R ${reference} \\
        -L ${interval} \\
        --bqsr-recal-file ${bqsr_table} \\
        -O ${task_id}_bqsr.bam \\
        --tmp-dir ${tmp_dir}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | grep -oP 'Version:\\s*\\K[0-9.]+' || echo "unknown")
    END_VERSIONS
    """

    stub:
    """
    touch ${task_id}_bqsr.bam
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: 4.4.0.0
    END_VERSIONS
    """
}
