process VARIANTRECALIBRATOR_SPLIT {
    publishDir "${params.outdir}/variantrecalibrator/", mode: 'symlink'

    cpus 60
    memory '100 GB'

    input:
    path(merged_vcf)
    path(merged_vcf_index)
    tuple path(reference), path(reference_fai), path(reference_dict)
    path(hapmap_files)
    path(omni_files)
    path(phase1_files)
    path(dbsnp_files)
    path(mills_files)

    output:
    path "${params.samplename_joint}_snps.recal", emit: snps_recal
    path "${params.samplename_joint}_snps.tranches", emit: snps_tranches
    path "${params.samplename_joint}_snps.recal.idx", emit: snps_recal_idx
    path "${params.samplename_joint}_indels.recal", emit: indels_recal
    path "${params.samplename_joint}_indels.tranches", emit: indels_tranches
    path "${params.samplename_joint}_indels.recal.idx", emit: indels_recal_idx
    path "${params.samplename_joint}_snps.vcf.gz", emit: snps_vcf
    path "${params.samplename_joint}_snps.vcf.gz.tbi", emit: snps_vcf_index
    path "${params.samplename_joint}_indels.vcf.gz", emit: indels_vcf
    path "${params.samplename_joint}_indels.vcf.gz.tbi", emit: indels_vcf_index
    path "versions.yml", emit: versions

    script:
    def tmp_dir = "/mnt/genomics/tmp"
    
    // Extract only VCF files from the collections
    def hapmap_vcf = hapmap_files.find { it.toString().endsWith('.vcf.gz') && !it.toString().endsWith('.tbi') }
    def omni_vcf = omni_files.find { it.toString().endsWith('.vcf.gz') && !it.toString().endsWith('.tbi') }
    def phase1_vcf = phase1_files.find { it.toString().endsWith('.vcf.gz') && !it.toString().endsWith('.tbi') }
    def dbsnp_vcf = dbsnp_files.find { it.toString().endsWith('.vcf.gz') && !it.toString().endsWith('.tbi') }
    def mills_vcf = mills_files.find { it.toString().endsWith('.vcf.gz') && !it.toString().endsWith('.tbi') }
    
    """
    mkdir -p ${tmp_dir}
    echo "Starting VariantRecalibrator for whole genome merged VCF"
    echo "Input merged VCF: ${merged_vcf}"
    echo "Reference: ${reference}"
    echo "Resources: HapMap=${hapmap_vcf}, Omni=${omni_vcf}, 1000G=${phase1_vcf}, dbSNP=${dbsnp_vcf}, Mills=${mills_vcf}"

    set -euo pipefail

    # Split SNPs first
    echo "Splitting SNPs from merged VCF..."
    gatk --java-options \\
        "-XX:+UseParallelGC \\
        -XX:ParallelGCThreads=20 \\
        -Xms${task.memory.toGiga()-6}g \\
        -Xmx${task.memory.toGiga()-6}g \\
        -Djava.io.tmpdir=${tmp_dir}" \\
        SelectVariants \\
        -R ${reference} \\
        -V ${merged_vcf} \\
        --select-type-to-include SNP \\
        -O ${params.samplename_joint}_snps.vcf.gz \\
        --tmp-dir ${tmp_dir}

    echo "SNP splitting completed successfully"

    # Index SNPs VCF
    echo "Indexing SNPs VCF..."
    gatk --java-options \\
        "-Xms4g -Xmx4g \\
        -Djava.io.tmpdir=${tmp_dir}" \\
        IndexFeatureFile \\
        -I ${params.samplename_joint}_snps.vcf.gz

    echo "SNPs VCF indexing completed successfully"

    # Split INDELs
    echo "Splitting INDELs from merged VCF..."
    gatk --java-options \\
        "-XX:+UseParallelGC \\
        -XX:ParallelGCThreads=20 \\
        -Xms${task.memory.toGiga()-6}g \\
        -Xmx${task.memory.toGiga()-6}g \\
        -Djava.io.tmpdir=${tmp_dir}" \\
        SelectVariants \\
        -R ${reference} \\
        -V ${merged_vcf} \\
        --select-type-to-include INDEL \\
        -O ${params.samplename_joint}_indels.vcf.gz \\
        --tmp-dir ${tmp_dir}

    echo "INDEL splitting completed successfully"

    # Index INDELs VCF
    echo "Indexing INDELs VCF..."
    gatk --java-options \\
        "-Xms4g -Xmx4g \\
        -Djava.io.tmpdir=${tmp_dir}" \\
        IndexFeatureFile \\
        -I ${params.samplename_joint}_indels.vcf.gz

    echo "INDELs VCF indexing completed successfully"

    # SNP recalibration on whole genome
    echo "Starting SNP recalibration with VQSR on whole genome..."
    gatk --java-options \\
        "-XX:+UseParallelGC \\
        -XX:ParallelGCThreads=${task.cpus} \\
        -Xms${task.memory.toGiga()-6}g \\
        -Xmx${task.memory.toGiga()-6}g \\
        -Djava.io.tmpdir=${tmp_dir}" \\
        VariantRecalibrator \\
        -R ${reference} \\
        -V ${params.samplename_joint}_snps.vcf.gz \\
        --resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${hapmap_vcf} \\
        --resource:omni,known=false,training=true,truth=true,prior=12.0 ${omni_vcf} \\
        --resource:1000G,known=false,training=true,truth=false,prior=10.0 ${phase1_vcf} \\
        --resource:dbsnp,known=true,training=false,truth=false,prior=7.0 ${dbsnp_vcf} \\
        -an QD -an FS -an SOR -an MQ \\
        -mode SNP \\
        --max-gaussians 4 \\
        -O ${params.samplename_joint}_snps.recal \\
        --tranches-file ${params.samplename_joint}_snps.tranches \\
        --tmp-dir ${tmp_dir}

    echo "SNP recalibration completed successfully"

    # INDEL recalibration on whole genome
    echo "Starting INDEL recalibration with VQSR on whole genome..."
    gatk --java-options \\
        "-XX:+UseParallelGC \\
        -XX:ParallelGCThreads=${task.cpus} \\
        -Xms${task.memory.toGiga()-6}g \\
        -Xmx${task.memory.toGiga()-6}g \\
        -Djava.io.tmpdir=${tmp_dir}" \\
        VariantRecalibrator \\
        -R ${reference} \\
        -V ${params.samplename_joint}_indels.vcf.gz \\
        --resource:mills,known=false,training=true,truth=true,prior=12.0 ${mills_vcf} \\
        --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${dbsnp_vcf} \\
        -an QD -an FS -an SOR \\
        -mode INDEL \\
        --max-gaussians 4 \\
        -O ${params.samplename_joint}_indels.recal \\
        --tranches-file ${params.samplename_joint}_indels.tranches \\
        --tmp-dir ${tmp_dir}

    echo "INDEL recalibration completed successfully"
    echo "VariantRecalibrator completed successfully for whole genome!"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | grep -oP 'Version:\\s*\\K[0-9.]+' || echo "unknown")
    END_VERSIONS
    """

    stub:
    """
    touch ${params.samplename_joint}_snps.recal ${params.samplename_joint}_snps.tranches ${params.samplename_joint}_snps.recal.idx
    touch ${params.samplename_joint}_indels.recal ${params.samplename_joint}_indels.tranches ${params.samplename_joint}_indels.recal.idx
    touch ${params.samplename_joint}_snps.vcf.gz ${params.samplename_joint}_snps.vcf.gz.tbi
    touch ${params.samplename_joint}_indels.vcf.gz ${params.samplename_joint}_indels.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: 4.4.0.0
    END_VERSIONS
    """
}
