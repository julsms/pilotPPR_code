process HAPLOTYPE_SAMPLING {
    tag "$sample_name"
    publishDir "${params.outdir}/haplotype_sampling", mode: 'copy'
    
    container 'quay.io/vgteam/vg:v1.52.0'
    containerOptions "--user \$(id -u):\$(id -g) -v ${params.graph_dir}:/graph"
    
    input:
    tuple val(sample_name), path(kff_file)
    
    output:
    tuple val(sample_name), path("${sample_name}.gbz"), emit: gbz
    
    script:
    """
    vg haplotypes -v 2 -t ${task.cpus} \\
        --include-reference --diploid-sampling \\
        -i /graph/graph.hapl \\
        -k ${kff_file} \\
        -g ${sample_name}.gbz \\
        /graph/hprc-v1.1-mc-grch38.gbz
    """
}
