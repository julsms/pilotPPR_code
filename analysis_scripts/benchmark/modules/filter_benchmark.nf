process FILTER_BENCHMARK_REGIONS {
    tag "${sample_type}"
    
    cpus 2
    memory '4 GB'
    
    input:
    tuple path(vcf), val(sample_type)
    path benchmark_bed
    
    output:
    tuple path("${sample_type}_filtered.vcf.gz"), val(sample_type), emit: filtered_vcf
    
    script:
    """
    echo "Pre-filtering ${sample_type} VCF to benchmark high-confidence regions..."
    echo "Input VCF: ${vcf}"
    echo "Benchmark BED: ${benchmark_bed}"
    
    # Check if input VCF is compressed
    if [[ ${vcf} == *.gz ]]; then
        input_vcf="${vcf}"
    else
        # Compress if needed
        bgzip -c ${vcf} > ${vcf}.gz
        tabix -p vcf ${vcf}.gz
        input_vcf="${vcf}.gz"
    fi
    
    # Index if not already indexed
    if [ ! -f "\${input_vcf}.tbi" ]; then
        tabix -p vcf \${input_vcf}
    fi
    
    # Filter VCF to benchmark regions
    bcftools view \\
        -R ${benchmark_bed} \\
        \${input_vcf} \\
        -Oz \\
        -o ${sample_type}_filtered.vcf.gz \\
        --threads ${task.cpus}
    
    # Print filtering stats
    original_count=\$(bcftools view -H \${input_vcf} | wc -l)
    filtered_count=\$(bcftools view -H ${sample_type}_filtered.vcf.gz | wc -l)
    
    echo "Original variants: \${original_count}"
    echo "After benchmark filtering: \${filtered_count}"
    echo "Filtered out: \$(( \${original_count} - \${filtered_count} )) variants"
    """
    
    stub:
    """
    touch ${sample_type}_filtered.vcf.gz
    """
}
