process KMC_COUNT {
    tag "$sample_name"
    publishDir "${params.outdir}/kmc", mode: 'copy'
    
    container 'quay.io/biocontainers/kmc:3.2.4--h5ca1c30_4'
    
    input:
    tuple val(sample_name), path(reads)
    
    output:
    tuple val(sample_name), path("${sample_name}.kff"), emit: kff
    path "kmc_tmp", emit: tmp_dir
    
    script:
    """
    # Create input file list
    echo "${reads[0]}" > kmc_input.txt
    echo "${reads[1]}" >> kmc_input.txt
    
    # Create temporary directory
    mkdir -p kmc_tmp
    
    # Run KMC
    kmc -k29 -m120 -okff -t${task.cpus} -hp \\
        @kmc_input.txt \\
        ${sample_name} \\
        kmc_tmp
    """
}
