process HAPLOTYPECALLER {
    tag "${interval}"
    publishDir "${params.outdir}/${sample_id}/gvcf/", mode: 'symlink'
    
    input:
    tuple val(sample_id), path(bam), path(bai), path(reference), path(reference_fai), path(reference_dict), path(reference_amb), path(reference_ann), path(reference_bwt), path(reference_pac), path(reference_sa), val(interval)

    output:
    tuple val(interval), path("*.g.vcf.gz"), emit: gvcf

    script:
    def sample_name = bam.baseName.replaceAll('_bqsr', '')
    def interval_str = interval.toString()
    def interval_name = interval_str.replaceAll('[^A-Za-z0-9]', '_')
    """
    gatk --java-options \\
        "-Xms${task.memory.toGiga()}g \\
        -Xmx${task.memory.toGiga()}g" \\
        HaplotypeCaller \\
        -R ${reference} \\
        -I ${bam} \\
        -L "${interval_str}" \\
        --interval-padding 1000 \\
        --max-reads-per-alignment-start 0 \\
        --native-pair-hmm-threads ${task.cpus} \\
        -ERC GVCF \\
        -O ${sample_name}_${interval_name}.g.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | grep -oP 'Version:\\s*\\K[0-9.]+' || echo "unknown")
    END_VERSIONS
    """

    stub:
    def sample_name = "sample"
    def clean_interval = interval.replaceAll(/[:\-\s]/, "_")
    def output_file = "${sample_name}.${clean_interval}.g.vcf.gz"
    """
    touch ${output_file}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: 4.4.0.0
    END_VERSIONS
    """
}
