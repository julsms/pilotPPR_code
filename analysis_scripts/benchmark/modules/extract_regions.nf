process EXTRACT_REGIONS {
    tag "${sample_type}_${bed_name}"
    publishDir "${params.outdir}", mode: 'copy'
    
    cpus 2
    memory '4 GB'
    
    input:
    tuple path(vcf), path(vcf_index), val(sample_type), val(bed_name), path(bed_file)
    
    output:
    tuple val(sample_type), val(bed_name), path("${sample_type}_${bed_basename}.vcf"), emit: region_vcf
    path "versions.yml", emit: versions
    
    script:
    bed_basename = bed_name.replaceAll('\\.bed$', '')
    """
    echo "Extracting regions for ${sample_type} using ${bed_name}..."
    echo "Input VCF: ${vcf}"
    echo "BED file: ${bed_file}"
    echo "Output: ${sample_type}_${bed_basename}.vcf"
    
    # Extract regions using BED file
    bcftools view \\
        -R ${bed_file} \\
        ${vcf} \\
        -o ${sample_type}_${bed_basename}.vcf \\
        --threads ${task.cpus}
    
    echo "Region extraction completed for ${sample_type}_${bed_basename}"
    echo "Output variants: \$(grep -v '^#' ${sample_type}_${bed_basename}.vcf | wc -l)"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -n1 | sed 's/bcftools //')
    END_VERSIONS
    """
    
    stub:
    """
    touch ${sample_type}_${bed_basename}.vcf
    echo '"${task.process}": {"bcftools": "1.17"}' > versions.yml
    """
}
