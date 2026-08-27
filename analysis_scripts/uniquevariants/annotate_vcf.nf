process ANNOTATE_VCF {
    tag "${category}/${sample_name}"
    
    cpus params.fork
    memory '24 GB'
    time '12h'
    errorStrategy { task.attempt <= 2 ? 'retry' : 'ignore' }
    maxRetries 2

    input:
    tuple val(category), val(sample_name), path(vcf), path(vcf_idx)

    output:
    tuple val(category), val(sample_name), path("${sample_name}.annotated.vcf.gz"), path("${sample_name}.annotated.vcf.gz.tbi"), emit: annotated_vcf
    path "${sample_name}_vep_summary.html", emit: summary
    path "${sample_name}_vep_warnings.txt", optional: true, emit: warnings

    script:
    """
    ${params.vep_dir}/vep \\
        --input_file ${vcf} \\
        --output_file ${sample_name}.annotated.vcf \\
        --format vcf --vcf --force_overwrite \\
        --species ${params.species} --assembly ${params.assembly} \\
        --fork ${task.cpus} --offline --cache \\
        --dir_cache ${params.vep_cache} --fasta ${params.reference} \\
        --everything --canonical --gene_phenotype --variant_class \\
        --check_existing --symbol --biotype \\
        --custom ${params.clinvar},ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN \\
        --custom ${params.dbsnp},dbSNP,vcf,exact,0,RS \\
        --stats_file ${sample_name}_vep_summary.html \\
        --warning_file ${sample_name}_vep_warnings.txt

    bgzip -c ${sample_name}.annotated.vcf > ${sample_name}.annotated.vcf.gz
    tabix -p vcf ${sample_name}.annotated.vcf.gz
    rm ${sample_name}.annotated.vcf
    
    # Copy annotated files back to original directory
    cp ${sample_name}.annotated.vcf.gz ${params.base_dir}/${category}/
    cp ${sample_name}.annotated.vcf.gz.tbi ${params.base_dir}/${category}/
    cp ${sample_name}_vep_summary.html ${params.base_dir}/${category}/
    cp ${sample_name}_vep_warnings.txt ${params.base_dir}/${category}/ 2>/dev/null || true
    """
}
