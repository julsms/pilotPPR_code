process GENERATE_CATEGORY_SUMMARY {
    tag "${category}"
    publishDir "${params.base_dir}/${category}", mode: 'copy'
    memory '4 GB'
    time '1h'

    input:
    tuple val(category), path(annotated_vcfs)

    output:
    path "${category}_summary.tsv", emit: tsv
    path "${category}_summary.txt", emit: txt

    script:
    """
    #!/bin/bash
    echo -e "Sample\\tTotal\\tAnnotated\\tNovel\\tKnown\\tClinVar\\tPathogenic\\tHighImpact\\tGenes" > ${category}_summary.tsv
    echo "Category: ${category}" > ${category}_summary.txt
    
    for vcf_gz in *.annotated.vcf.gz; do
        [ -f "\${vcf_gz}" ] || continue
        SAMPLE=\$(basename \${vcf_gz} .annotated.vcf.gz)
        TOTAL=\$(bcftools view -H \${vcf_gz} | wc -l)
        ANNOTATED=\$(bcftools view -H \${vcf_gz} | grep -c "CSQ=" || echo 0)
        KNOWN=\$(bcftools query -f '%ID\\n' \${vcf_gz} | awk '\$1 != "."' | wc -l)
        NOVEL=\$((TOTAL - KNOWN))
        CLINVAR=\$(bcftools view -H \${vcf_gz} | grep -c "ClinVar" || echo 0)
        PATHOGENIC=\$(bcftools view -H \${vcf_gz} | grep -E "Pathogenic|Likely_pathogenic" | grep -v "Benign" | wc -l || echo 0)
        HIGH=\$(bcftools view -H \${vcf_gz} | grep -E "HIGH|stop_gained|frameshift|splice" | wc -l || echo 0)
        GENES=\$(bcftools query -f '%INFO/CSQ\\n' \${vcf_gz} 2>/dev/null | grep -oP 'SYMBOL=[^;|]+' | cut -d= -f2 | sort -u | wc -l || echo 0)
        echo -e "\${SAMPLE}\\t\${TOTAL}\\t\${ANNOTATED}\\t\${NOVEL}\\t\${KNOWN}\\t\${CLINVAR}\\t\${PATHOGENIC}\\t\${HIGH}\\t\${GENES}" >> ${category}_summary.tsv
        echo "\${SAMPLE}: \${TOTAL} total, \${PATHOGENIC} pathogenic, \${HIGH} high impact" >> ${category}_summary.txt
    done
    """
}

process GENERATE_GLOBAL_SUMMARY {
    publishDir "${params.base_dir}/summaries", mode: 'copy'
    memory '4 GB'
    time '1h'

    input:
    path all_vcfs
    path category_tsvs

    output:
    path "global_annotation_summary.tsv"
    path "global_annotation_summary.txt"

    script:
    """
    #!/bin/bash
    echo -e "Category\\tSample\\tTotal\\tAnnotated\\tNovel\\tKnown\\tClinVar\\tPathogenic\\tHighImpact\\tGenes" > global_annotation_summary.tsv
    for tsv in *_summary.tsv; do
        CATEGORY=\$(basename \${tsv} _summary.tsv)
        tail -n +2 \${tsv} | awk -v cat="\${CATEGORY}" '{print cat"\\t"\$0}' >> global_annotation_summary.tsv
    done
    
    echo "GLOBAL SUMMARY - 4 Multisample VCF files" > global_annotation_summary.txt
    awk 'NR>1 {cat[\$1]++; total[\$1]+=\$3; path[\$1]+=\$8; high[\$1]+=\$9}
         END {for (c in cat) print c": "cat[c]" files, "total[c]" variants, "path[c]" pathogenic"}' \\
         global_annotation_summary.tsv >> global_annotation_summary.txt
    """
}
