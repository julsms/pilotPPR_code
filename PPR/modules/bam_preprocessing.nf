process BAM_PREPROCESSING {
    tag "$sample_name"
    publishDir "${params.outdir}/small_variants", mode: 'copy'
    
    container 'quay.io/biocontainers/samtools:1.17--h00cdaf9_0'
    containerOptions "--user \$(id -u):\$(id -g)"
    
    input:
    tuple val(sample_name), path(bam_file)
    
    output:
    tuple val(sample_name), path("${sample_name}_chr_named_rg.bam"), path("${sample_name}_chr_named_rg.bam.bai"), emit: bam
    
    script:
    // Calculează memoria per thread pentru samtools sort
    def mem_per_thread = Math.max(1, (task.memory.toGiga() / task.cpus).intValue())
    """
    # Sort BAM with optimized memory usage
    samtools sort -@ ${task.cpus} -m ${mem_per_thread}G -o ${sample_name}.sorted.bam ${bam_file}
    
    # Fix header (remove GRCh38#0# prefix)
    samtools view -H ${sample_name}.sorted.bam | sed 's/SN:GRCh38#0#chr/SN:chr/' > new_header.sam
    
    # Apply new header
    samtools reheader new_header.sam ${sample_name}.sorted.bam > ${sample_name}_chr_named.bam
    
    # Add read groups with threads
    samtools addreplacerg \\
        --threads ${task.cpus} \\
        -r '@RG\\tID:${sample_name}\\tSM:${sample_name}\\tPL:ILLUMINA\\tLB:lib1\\tPU:unit1' \\
        -o ${sample_name}_chr_named_rg.bam \\
        ${sample_name}_chr_named.bam
    
    # Index final BAM with threads
    samtools index -@ ${task.cpus} ${sample_name}_chr_named_rg.bam
    
    echo 'BAM processing with read groups completed successfully!'
    """
}
