process VG_PACK {
    tag "$sample_name"
    publishDir "${params.outdir}/SV_calling", mode: 'copy'
    
    container 'quay.io/vgteam/vg:v1.52.0'
    containerOptions "--user \$(id -u):\$(id -g) -v ${params.graph_dir}:/graph"
    
    input:
    tuple val(sample_name), path(gam_file)
    tuple val(sample_name), path(gbz_file)
    
    output:
    tuple val(sample_name), path("${sample_name}.pack"), emit: pack
    
    script:
    """
    vg pack -x /graph/hprc-v1.1-mc-grch38.gbz \\
        -g ${gam_file} \\
        -Q 5 \\
        -t ${task.cpus} \\
        -o ${sample_name}.pack
    """
}
