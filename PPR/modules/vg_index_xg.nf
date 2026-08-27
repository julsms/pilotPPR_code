process VG_INDEX_XG {
    tag "$sample_name"
    publishDir "${params.outdir}/small_variants/intermediate", mode: 'copy'
    
    container 'quay.io/vgteam/vg:v1.52.0'
    containerOptions "--user \$(id -u):\$(id -g)"
    
    input:
    tuple val(sample_name), path(vg_file)
    
    output:
    tuple val(sample_name), path("${sample_name}.xg"), emit: xg
    
    script:
    """
    vg index -t ${task.cpus} -x ${sample_name}.xg ${vg_file}
    """
}
