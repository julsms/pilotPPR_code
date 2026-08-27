process BWAMEM_SORT {
    tag "$sample_id"
    publishDir "${params.outdir}/aligned_bam/", mode: 'symlink'

    cpus 65
    memory '140 GB'
    time '12.h'

    input:
    tuple val(sample_id), path(read1), path(read2)
    tuple path(reference), path(reference_fai), path(reference_dict), path(reference_amb), path(reference_ann), path(reference_bwt), path(reference_pac), path(reference_sa)

    output:
    tuple val(sample_id), path("${sample_id}_sorted.bam"), emit: bam
    path "versions.yml", emit: versions

    script:
    def read_group = "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA\\tLB:${sample_id}_lib\\tPU:${sample_id}_unit"
    def tmp_dir = "/mnt/genomics/tmp"
    def bwa_cpus = 50
    def sort_cpus = 15
    def memory_per_thread = 8
    
    """
    mkdir -p ${tmp_dir}
    
    bwa mem \\
        -t ${bwa_cpus} \\
        -R "${read_group}" \\
        "${reference}" \\
        "${read1}" \\
        "${read2}" \\
    | samtools sort \\
        -@ ${sort_cpus} \\
        -m ${memory_per_thread}G \\
        -T "${tmp_dir}/${sample_id}_sort_tmp" \\
        -o "${sample_id}_sorted.bam" \\
        -

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(bwa 2>&1 | grep -E '^Version' | sed 's/Version: //')
        samtools: \$(samtools --version 2>&1 | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${sample_id}_sorted.bam
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: 0.7.17
        samtools: 1.17
    END_VERSIONS
    """
}
