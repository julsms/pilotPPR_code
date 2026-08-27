process MERGE_INTERVAL_VCFS {
    publishDir "${params.outdir}/merged_vcf/", mode: 'symlink'

    cpus 60
    memory '60 GB'

    input:
    path(interval_vcfs)
    tuple path(reference), path(reference_fai), path(reference_dict)

    output:
    path "${params.samplename_joint}_merged_raw_variants.vcf.gz", emit: merged_vcf
    path "${params.samplename_joint}_merged_raw_variants.vcf.gz.tbi", emit: merged_vcf_index
    path "versions.yml", emit: versions

    script:
    def tmp_dir = "/mnt/genomics/tmp"
    def vcf_inputs = interval_vcfs.collect{ "--INPUT ${it}" }.join(' ')
    """
    mkdir -p ${tmp_dir}
    echo "Starting merge of interval VCFs..."
    echo "Input interval VCFs: ${interval_vcfs.join(', ')}"
    echo "Number of interval VCFs to merge: ${interval_vcfs.size()}"

    set -euo pipefail

    # List input files for debugging
    ls -lh *.vcf.gz

    # Merge all interval VCFs using Picard MergeVcfs
    echo "Merging all interval VCFs with Picard MergeVcfs..."
    java \\
        -XX:+UseParallelGC \\
        -XX:ParallelGCThreads=${task.cpus} \\
        -Xmx${task.memory.toGiga()-2}g \\
        -Djava.io.tmpdir=${tmp_dir} \\
        -jar ${params.picard} MergeVcfs \\
        ${vcf_inputs} \\
        --OUTPUT ${params.samplename_joint}_merged_raw_variants.vcf.gz \\
        --SEQUENCE_DICTIONARY ${reference_dict} \\
        --TMP_DIR ${tmp_dir}

    echo "VCF merging completed successfully"

    # Index the merged VCF
    echo "Indexing merged VCF..."
    gatk --java-options \\
        "-Xmx4g \\
        -Djava.io.tmpdir=${tmp_dir}" \\
        IndexFeatureFile \\
        -I ${params.samplename_joint}_merged_raw_variants.vcf.gz

    echo "VCF merge and indexing completed successfully!"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(java -jar ${params.picard} MergeVcfs --version 2>&1 | grep -oP 'Version:\\s*\\K[0-9.]+' || echo "unknown")
        gatk4: \$(gatk --version 2>&1 | grep -oP 'Version:\\s*\\K[0-9.]+' || echo "unknown")
    END_VERSIONS
    """

    stub:
    """
    touch ${params.samplename_joint}_merged_raw_variants.vcf.gz
    touch ${params.samplename_joint}_merged_raw_variants.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: 3.0.0
        gatk4: 4.4.0.0
    END_VERSIONS
    """
}
