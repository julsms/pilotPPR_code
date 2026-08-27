process VG_VIEW {
    tag "$sample_name"
    publishDir "${params.outdir}/small_variants/intermediate", mode: 'copy'
    
    container 'quay.io/vgteam/vg:v1.52.0'
    containerOptions "--user \$(id -u):\$(id -g)"
    
    input:
    tuple val(sample_name), path(gbz_file)
    
    output:
    tuple val(sample_name), path("${sample_name}.vg"), emit: vg
    
    script:
    """
    vg view -v ${gbz_file} --threads ${task.cpus} > ${sample_name}.vg
    """
}
