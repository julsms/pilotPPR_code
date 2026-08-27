process NORMALIZE_VCF {
    tag "${sample_type}"
    publishDir "${params.outdir}", mode: 'copy'
    
    cpus 4
    memory '8 GB'
    
    input:
    tuple path(vcf), val(sample_type)
    path reference
    
    output:
    tuple path("${sample_type}.vcf.gz"), path("${sample_type}.vcf.gz.tbi"), val(sample_type), emit: normalized_vcf
    path "versions.yml", emit: versions
    
    script:
    """
    echo "Normalizing ${sample_type} VCF..."
    echo "Input VCF: ${vcf}"
    echo "Output: ${sample_type}.vcf.gz"
    
    # Normalize VCF
    bcftools norm \\
        -f ${reference} \\
        -m -any \\
        ${vcf} \\
        -Oz \\
        -o ${sample_type}.vcf.gz \\
        --threads ${task.cpus}
    
    # Index the normalized VCF
    tabix -p vcf ${sample_type}.vcf.gz
    
    echo "Normalization completed for ${sample_type}"
    echo "Output size: \$(ls -lh ${sample_type}.vcf.gz | awk '{print \$5}')"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -n1 | sed 's/bcftools //')
        tabix: \$(tabix --version | head -n1 | sed 's/tabix (htslib) //')
    END_VERSIONS
    """
    
    stub:
    """
    touch ${sample_type}.vcf.gz
    touch ${sample_type}.vcf.gz.tbi
    echo '"${task.process}": {"bcftools": "1.17", "tabix": "1.17"}' > versions.yml
    """
}
