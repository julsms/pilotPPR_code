process GENOMICSDBIMPORT {
    tag "${interval}"
    
    input:
    tuple val(interval), path(gvcf), path(tbi)
    tuple path(reference), path(reference_fai), path(reference_dict), path(reference_amb), path(reference_ann), path(reference_bwt), path(reference_pac), path(reference_sa)

    output:
    tuple val(interval), path("${interval.toString().replaceAll('[^A-Za-z0-9]', '_')}_genomicsdb"), emit: database
    path "versions.yml", emit: versions

    script:
    def interval_safe = interval.toString().replaceAll('[^A-Za-z0-9]', '_')
    def tmp_dir = System.getenv('TMP_DIR') ?: "/tmp/genomics"
    def persistent_db_path = "${params.genomicsdb_persistent_path}/${interval_safe}_genomicsdb"
    
    """
    mkdir -p ${tmp_dir}
    mkdir -p ${params.genomicsdb_persistent_path}/
    
    echo "Starting GenomicsDBImport for interval: ${interval}"
    echo "Input GVCF: ${gvcf}"
    echo "Input TBI: ${tbi}"
    echo "Persistent DB path: ${persistent_db_path}"

    # Try to update existing database first, fallback to create if it fails
    if [ -d "${persistent_db_path}" ]; then
        echo "Found existing database, attempting update..."
        gatk --java-options \\
            "-XX:+UseParallelGC \\
            -XX:ParallelGCThreads=${task.cpus} \\
            -Xms${task.memory.toGiga()}g \\
            -Xmx${task.memory.toGiga()}g \\
            -Djava.io.tmpdir=${tmp_dir}" \\
            GenomicsDBImport \\
            -V ${gvcf} \\
            --genomicsdb-update-workspace-path ${persistent_db_path} \\
            -L ${interval} \\
            --tmp-dir ${tmp_dir} \\
        && echo "Database update completed successfully" \\
        || {
            echo "Update failed, this shouldn't happen with existing database"
            exit 1
        }
    else
        echo "No existing database found, creating new one..."
        gatk --java-options \\
            "-XX:+UseParallelGC \\
            -XX:ParallelGCThreads=${task.cpus} \\
            -Xms${task.memory.toGiga()}g \\
            -Xmx${task.memory.toGiga()}g \\
            -Djava.io.tmpdir=${tmp_dir}" \\
            GenomicsDBImport \\
            -V ${gvcf} \\
            --genomicsdb-workspace-path ${persistent_db_path} \\
            -L ${interval} \\
            --tmp-dir ${tmp_dir} \\
        && echo "Database creation completed successfully"
    fi

    # Create symlink for Nextflow output (required for publishDir)
    ln -sf ${persistent_db_path} ${interval_safe}_genomicsdb

    echo "GenomicsDBImport completed successfully for interval: ${interval}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | grep -oP 'Version:\\s*\\K[0-9.]+' || echo "unknown")
    END_VERSIONS
    """

    stub:
    def interval_safe = interval.toString().replaceAll('[^A-Za-z0-9]', '_')
    """
    mkdir -p ${interval_safe}_genomicsdb
    touch ${interval_safe}_genomicsdb/stub_file
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: 4.4.0.0
    END_VERSIONS
    """
}
