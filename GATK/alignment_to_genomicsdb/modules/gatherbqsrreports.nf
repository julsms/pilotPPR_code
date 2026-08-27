process GATHERBQSRREPORTS {
    publishDir "${params.outdir}/${sample_id}/bqsr", mode: 'copy'

    input:
    tuple val(sample_id), path(bqsr_tables)

    output:
    tuple val(sample_id), path("${sample_id}_bqsr.table"), emit: table
    path "versions.yml", emit: versions

    script:
    def tmp_dir = "/mnt/genomics/tmp"
    def input_tables = bqsr_tables.collect{ "-I ${it}" }.join(' ')
    """
    echo "Gathering BQSR reports for sample: ${sample_id}"
    echo "Input tables: ${bqsr_tables}"

    gatk --java-options \\
        "-Xms${task.memory.toGiga()-2}g
        -Xmx${task.memory.toGiga()-2}g \\
        -Djava.io.tmpdir=${tmp_dir}" \\
        GatherBQSRReports \\
        ${input_tables} \\
        -O ${sample_id}_bqsr.table \\
        --tmp-dir ${tmp_dir}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | grep -oP 'Version:\\s*\\K[0-9.]+' || echo "unknown")
    END_VERSIONS
    """

    stub:
    """
    touch ${sample_id}_bqsr.table
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: 4.4.0.0
    END_VERSIONS
    """
}
