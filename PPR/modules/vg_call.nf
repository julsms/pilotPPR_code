process VG_CALL {
    tag "$sample_name"
    publishDir "${params.outdir}/SV_calling", mode: 'copy'
    
    container 'quay.io/vgteam/vg:v1.52.0'
    containerOptions "--user \$(id -u):\$(id -g) -v ${params.graph_dir}:/graph"
    
    input:
    tuple val(sample_name), path(pack_file)
    tuple val(sample_name), path(gbz_file)
    
    output:
    tuple val(sample_name), path("${sample_name}_vgcall.vcf"), emit: vcf
    path "${sample_name}_vgcall.err", emit: log
    
    script:
    """
    vg call /graph/hprc-v1.1-mc-grch38.gbz \\
        -k ${pack_file} \\
        -S GRCh38 \\
        -t ${task.cpus} \\
        -z > ${sample_name}_vgcall.vcf 2> ${sample_name}_vgcall.err
    """
}
