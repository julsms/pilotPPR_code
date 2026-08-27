process COMPRESS_INDEX_VCF {
    tag "$sample_name"
    publishDir "${params.outdir}/SV_calling", mode: 'copy'
    
    input:
    tuple val(sample_name), path(vcf_file)
    
    output:
    tuple val(sample_name), path("${sample_name}_vgcall.vcf.gz"), emit: vcf_gz
    tuple val(sample_name), path("${sample_name}_vgcall.vcf.gz.tbi"), emit: vcf_index
    
    script:
    """
    # Compress VCF
    bgzip -c ${vcf_file} > ${sample_name}_vgcall.vcf.gz
    
    # Index compressed VCF
    tabix -p vcf ${sample_name}_vgcall.vcf.gz
    """
}
