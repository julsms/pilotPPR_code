process BASERECALIBRATOR {

    input:
    tuple val(task_id), path(bam), val(interval)
    tuple val(task_id), path(bai)
    tuple path(reference), path(reference_fai), path(reference_dict), path(reference_amb), path(reference_ann), path(reference_bwt), path(reference_pac), path(reference_sa)
    tuple path(resources_knownindels), path(resources_knownindels_tbi)
    tuple path(resources_mills), path(resources_mills_tbi)

    output:
    tuple val(task_id), path("${task_id}.table"), emit: table
    path "versions.yml", emit: versions

    script:
    def tmp_dir = "/mnt/genomics/tmp"
    def output_file = "${task_id}.table"

    """
    echo "Running BaseRecalibrator for task: ${task_id}, region: ${interval}"
    echo "Output file will be: ${output_file}"

    gatk --java-options \\
        "-Xms${task.memory.toGiga()}g \\
        -Xmx${task.memory.toGiga()}g \\
        -XX:+UseParallelGC \\
        -XX:ParallelGCThreads=${task.cpus} \\
        -Djava.io.tmpdir=${tmp_dir}" \\
        BaseRecalibrator \\
        -I ${bam} \\
        -R ${reference} \\
        -L ${interval} \\
        -O ${output_file} \\
        --known-sites ${resources_knownindels} \\
        --known-sites ${resources_mills} \\
        --tmp-dir ${tmp_dir}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | grep -oP 'Version:\\s*\\K[0-9.]+' || echo "unknown")
    END_VERSIONS
    """

    stub:
    def output_file = "${task_id}.table"
    """
    touch ${output_file}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: 4.4.0.0
    END_VERSIONS
    """
}
