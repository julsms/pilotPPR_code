process VG_GIRAFFE {
    tag "$sample_name"
    publishDir "${params.outdir}/alignment", mode: 'copy'
    
    container 'quay.io/vgteam/vg:v1.52.0'
    
    input:
    tuple val(sample_name), path(gbz_file)
    tuple val(sample_name), path(reads)
    
    output:
    tuple val(sample_name), path("${sample_name}.gam"), emit: gam
    path "${sample_name}_align.err", emit: log
    
    script:
    """
    vg giraffe -v 2 -p -t ${task.cpus} -Z ${gbz_file} \\
        -f ${reads[0]} \\
        -f ${reads[1]} \\
        > ${sample_name}.gam 2> ${sample_name}_align.err
    """
}
