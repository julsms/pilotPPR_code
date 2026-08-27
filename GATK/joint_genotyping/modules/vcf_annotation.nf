process VCF_ANNOTATION {
    tag "$final_vcf"
    publishDir "${params.outdir}/annotated_vcf/", mode: 'copy'

    cpus 8
    memory '16.GB'

    input:
    tuple path(final_vcf), path(final_vcf_index)
    path(bed_files)

    output:
    path "${params.samplename_joint}_annotated_variants.vcf.gz", emit: annotated_vcf
    path "${params.samplename_joint}_annotated_variants.vcf.gz.tbi", emit: annotated_vcf_index
    path "versions.yml", emit: versions

    script:
    """
    set -euo pipefail

    echo "Starting VCF annotation with vcfanno..."
    echo "Input VCF: ${final_vcf}"

    # List all files in work directory
    echo "All files in work directory:"
    ls -la

    # Index all BED files first
    echo "Indexing BED files..."
    for bed_file in *.bed; do
        [[ ! -f "\$bed_file" ]] && continue
        echo "Indexing \$bed_file..."
        bgzip -c "\$bed_file" > "\${bed_file}.gz"
        tabix -p bed "\${bed_file}.gz"
    done

    # Create vcfanno configuration file dynamically
    > vcfanno_all.conf

    # Process only .bed.gz files (indexed)
    for bed_file in *.bed.gz; do
        # Skip if file doesn't exist (glob didn't match)
        [[ ! -f "\$bed_file" ]] && continue
        
        name=\$(basename "\$bed_file" .bed.gz)
        
        echo "Processing BED file: \$bed_file (name: \$name)"
        
        # For flag operation with BED files, use columns = [1]
        cat >> vcfanno_all.conf << EOFCONF
[[annotation]]
file = "\$bed_file"
columns = [1]
names = ["\$name"]
ops = ["flag"]

EOFCONF
        
        echo "Added \$bed_file to configuration as flag annotation"
    done

    # Debug: show the generated config
    echo "Generated vcfanno config:"
    cat vcfanno_all.conf
    echo "Config file size: \$(wc -l < vcfanno_all.conf) lines"

    # Check if config file is not empty
    if [[ ! -s vcfanno_all.conf ]]; then
        echo "No valid BED files found for annotation. Creating copy of input VCF..."
        cp ${final_vcf} ${params.samplename_joint}_annotated_variants.vcf.gz
        cp ${final_vcf_index} ${params.samplename_joint}_annotated_variants.vcf.gz.tbi
    else
        echo "Running vcfanno annotation with \$(grep -c '\\[\\[annotation\\]\\]' vcfanno_all.conf) annotations..."
        
        # Run vcfanno with error handling
        if vcfanno \\
            -p ${task.cpus} \\
            vcfanno_all.conf \\
            ${final_vcf} \\
            2> vcfanno.log \\
            | bgzip -c > ${params.samplename_joint}_annotated_variants.vcf.gz; then
            
            echo "vcfanno completed successfully"
            
            # Check if output file was created and is not empty
            if [[ -s ${params.samplename_joint}_annotated_variants.vcf.gz ]]; then
                # Index the annotated VCF
                tabix -p vcf ${params.samplename_joint}_annotated_variants.vcf.gz
                echo "VCF indexing completed"
            else
                echo "vcfanno output is empty, using input VCF"
                cp ${final_vcf} ${params.samplename_joint}_annotated_variants.vcf.gz
                cp ${final_vcf_index} ${params.samplename_joint}_annotated_variants.vcf.gz.tbi
            fi
        else
            echo "vcfanno failed. Error log:"
            cat vcfanno.log
            echo "Creating copy of input VCF instead..."
            cp ${final_vcf} ${params.samplename_joint}_annotated_variants.vcf.gz
            cp ${final_vcf_index} ${params.samplename_joint}_annotated_variants.vcf.gz.tbi
        fi
    fi

    echo "VCF annotation process completed!"

    # Create versions file
    cat > versions.yml << EOFVER
"${task.process}":
    vcfanno: \$(vcfanno 2>&1 | head -n 1 | sed 's/vcfanno version //' || echo "unknown")
    bgzip: \$(bgzip --version 2>&1 | head -n 1 | sed 's/bgzip (htslib) //' || echo "unknown")
    tabix: \$(tabix --version 2>&1 | head -n 1 | sed 's/tabix (htslib) //' || echo "unknown")
EOFVER
    """

    stub:
    """
    touch ${params.samplename_joint}_annotated_variants.vcf.gz
    touch ${params.samplename_joint}_annotated_variants.vcf.gz.tbi
    touch versions.yml
    """
}
