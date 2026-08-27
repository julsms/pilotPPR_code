process APPLYVQSR_SPLIT {
    tag "ApplyVQSR"
    publishDir "${params.outdir}/applyvqsr/", mode: 'copy'
    publishDir "${params.outdir}/applyvqsr/snps/", mode: 'copy', pattern: "${params.samplename_joint}_snps_recalibrated.vcf.gz*"
    publishDir "${params.outdir}/applyvqsr/indels/", mode: 'copy', pattern: "${params.samplename_joint}_indels_recalibrated.vcf.gz*"
    
    cpus 20
    memory '12 GB'

    input:
    path merged_vcf
    path merged_vcf_index
    tuple path(reference), path(reference_fai), path(reference_dict)
    path snps_recal
    path snps_tranches  
    path snps_recal_idx
    path indels_recal
    path indels_tranches
    path indels_recal_idx
    path snps_vcf
    path snps_vcf_index
    path indels_vcf
    path indels_vcf_index

    output:
    tuple path("${params.samplename_joint}_merged_recalibrated.vcf.gz"), path("${params.samplename_joint}_merged_recalibrated.vcf.gz.tbi"), emit: merged_recalibrated
    tuple path("${params.samplename_joint}_snps_recalibrated.vcf.gz"), path("${params.samplename_joint}_snps_recalibrated.vcf.gz.tbi"), emit: snps_recalibrated
    tuple path("${params.samplename_joint}_indels_recalibrated.vcf.gz"), path("${params.samplename_joint}_indels_recalibrated.vcf.gz.tbi"), emit: indels_recalibrated
    path "versions.yml", emit: versions

    script:
    def tmp_dir = "/mnt/genomics/tmp"
    def memory_gb = task.memory.toGiga() - 6
    """
    mkdir -p ${tmp_dir}
    echo "Starting ApplyVQSR - applying to split SNP and INDEL VCFs"
    echo "SNPs VCF: ${snps_vcf}"
    echo "INDELs VCF: ${indels_vcf}"
    echo "SNPs recalibration: ${snps_recal}"
    echo "INDELs recalibration: ${indels_recal}"
    echo "Reference: ${reference}"

    set -euo pipefail

    # Apply VQSR to SNPs VCF with 99.0% sensitivity
    echo "Applying VQSR to SNPs VCF with 99.0% sensitivity threshold..."
    /usr/bin/java \\
        -XX:+UseParallelGC \\
        -XX:ParallelGCThreads=20 \\
        -Xmx${memory_gb}g \\
        -Djava.io.tmpdir=${tmp_dir} \\
        -jar /mnt/genomics/tools/gatk/gatk-package-4.5.0.0-local.jar ApplyVQSR \\
        -R ${reference} \\
        -V ${snps_vcf} \\
        --recal-file ${snps_recal} \\
        --tranches-file ${snps_tranches} \\
        --truth-sensitivity-filter-level 99.9 \\
        --create-output-variant-index true \\
        -mode SNP \\
        -O ${params.samplename_joint}_snps_recalibrated.vcf.gz \\
        --tmp-dir ${tmp_dir}

    echo "ApplyVQSR for SNPs completed successfully"
    if [[ ! -f "${params.samplename_joint}_snps_recalibrated.vcf.gz" ]]; then
        echo "ERROR: SNPs recalibrated VCF file was not created"
        exit 1
    fi

    echo "SNPs recalibrated VCF file size: \$(ls -lh ${params.samplename_joint}_snps_recalibrated.vcf.gz)"

    # Apply VQSR to INDELs VCF with 99.0% sensitivity  
    echo "Applying VQSR to INDELs VCF with 99.0% sensitivity threshold..."
    /usr/bin/java \\
        -XX:+UseParallelGC \\
        -XX:ParallelGCThreads=20 \\
        -Xmx${memory_gb}g \\
        -Djava.io.tmpdir=${tmp_dir} \\
        -jar /mnt/genomics/tools/gatk/gatk-package-4.5.0.0-local.jar ApplyVQSR \\
        -R ${reference} \\
        -V ${indels_vcf} \\
        --recal-file ${indels_recal} \\
        --tranches-file ${indels_tranches} \\
        --truth-sensitivity-filter-level 99.0 \\
        --create-output-variant-index true \\
        -mode INDEL \\
        -O ${params.samplename_joint}_indels_recalibrated.vcf.gz \\
        --tmp-dir ${tmp_dir}

    echo "ApplyVQSR for INDELs completed successfully"
    if [[ ! -f "${params.samplename_joint}_indels_recalibrated.vcf.gz" ]]; then
        echo "ERROR: INDELs recalibrated VCF file was not created"
        exit 1
    fi

    echo "INDELs recalibrated VCF file size: \$(ls -lh ${params.samplename_joint}_indels_recalibrated.vcf.gz)"

    # Merge the filtered SNPs and INDELs back together
    echo "Merging filtered SNPs and INDELs..."
    java \\
        -XX:+UseParallelGC \\
        -XX:ParallelGCThreads=20 \\
        -Xmx${memory_gb}g \\
        -Djava.io.tmpdir=${tmp_dir} \\
        -jar ${params.picard} MergeVcfs \\
        --INPUT ${params.samplename_joint}_snps_recalibrated.vcf.gz \\
        --INPUT ${params.samplename_joint}_indels_recalibrated.vcf.gz \\
        --OUTPUT ${params.samplename_joint}_merged_recalibrated.vcf.gz \\
        --SEQUENCE_DICTIONARY ${reference_dict} \\
        --TMP_DIR ${tmp_dir}

    echo "Merging completed successfully"

    # Index final merged VCF
    echo "Indexing final merged VCF..."
    /usr/bin/java \\
        -Xmx4g \\
        -Djava.io.tmpdir=${tmp_dir} \\
        -jar /mnt/genomics/tools/gatk/gatk-package-4.5.0.0-local.jar IndexFeatureFile \\
        -I ${params.samplename_joint}_merged_recalibrated.vcf.gz

    echo "Final recalibrated VCF file size: \$(ls -lh ${params.samplename_joint}_merged_recalibrated.vcf.gz)"
    echo "ApplyVQSR completed successfully - SNPs and INDELs filtered and merged!"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(/usr/bin/java -jar /mnt/genomics/tools/gatk/gatk-package-4.5.0.0-local.jar --version 2>&1 | grep -E '^The Genome Analysis Toolkit' | sed 's/The Genome Analysis Toolkit (GATK) v//' || echo "4.5.0.0")
        picard: \$(java -jar ${params.picard} MergeVcfs --version 2>&1 | grep -oP 'Version:\\s*\\K[0-9.]+' || echo "unknown")
        java: \$(/usr/bin/java -version 2>&1 | head -n1 | sed 's/openjdk version "//;s/".*//')
    END_VERSIONS
    """

    stub:
    """
    touch ${params.samplename_joint}_merged_recalibrated.vcf.gz
    touch ${params.samplename_joint}_merged_recalibrated.vcf.gz.tbi
    touch ${params.samplename_joint}_snps_recalibrated.vcf.gz
    touch ${params.samplename_joint}_snps_recalibrated.vcf.gz.tbi
    touch ${params.samplename_joint}_indels_recalibrated.vcf.gz
    touch ${params.samplename_joint}_indels_recalibrated.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: 4.5.0.0
        java: 11.0.27
    END_VERSIONS
    """
}
