process VG_SURJECT {
    tag "$sample_name"
    publishDir "${params.outdir}/small_variants/intermediate", mode: 'copy'
    
    container 'quay.io/vgteam/vg:v1.52.0'
    containerOptions "--user \$(id -u):\$(id -g)"
    
    input:
    tuple val(sample_name), path(gam_file)
    tuple val(sample_name), path(xg_file)
    
    output:
    tuple val(sample_name), path("${sample_name}_surject.bam"), emit: bam
    path "${sample_name}_surject.log", emit: log
    
    script:
    """
    vg surject \\
        -x ${xg_file} \\
        -b ${gam_file} \\
        -t ${task.cpus} \\
        --prune-low-cplx \\
        -p GRCh38#0#chr1 -p GRCh38#0#chr2 -p GRCh38#0#chr3 -p GRCh38#0#chr4 -p GRCh38#0#chr5 \\
        -p GRCh38#0#chr6 -p GRCh38#0#chr7 -p GRCh38#0#chr8 -p GRCh38#0#chr9 -p GRCh38#0#chr10 \\
        -p GRCh38#0#chr11 -p GRCh38#0#chr12 -p GRCh38#0#chr13 -p GRCh38#0#chr14 -p GRCh38#0#chr15 \\
        -p GRCh38#0#chr16 -p GRCh38#0#chr17 -p GRCh38#0#chr18 -p GRCh38#0#chr19 -p GRCh38#0#chr20 \\
        -p GRCh38#0#chr21 -p GRCh38#0#chr22 -p GRCh38#0#chrX -p GRCh38#0#chrY \\
        > ${sample_name}_surject.bam 2> ${sample_name}_surject.log
    """
}
