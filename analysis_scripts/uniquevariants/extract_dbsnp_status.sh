#!/bin/bash
# extract_dbsnp_status.sh

BASE_DIR="/mnt/genomics/pilot_PPR/uniquevariants"
OUTPUT_DIR="/mnt/genomics/pilot_PPR/uniquevariants/extract_csq_tsv"
mkdir -p "$OUTPUT_DIR"

VCFS=(
    "$BASE_DIR/uniqueGATKillumina/multisample_uniqueGATKillumina.annotated.vcf.gz"
    "$BASE_DIR/uniqueGATKaviti/multisample_uniqueGATKaviti.annotated.vcf.gz"
    "$BASE_DIR/uniquePPRillumina/multisample_uniquePPRillumina.annotated.vcf.gz"
    "$BASE_DIR/uniquePPRaviti/multisample_uniquePPRaviti.annotated.vcf.gz"
)

NAMES=("uniqueGATKillumina" "uniqueGATKaviti" "uniquePPRillumina" "uniquePPRaviti")

for i in "${!VCFS[@]}"; do
    vcf="${VCFS[$i]}"
    name="${NAMES[$i]}"
    
    echo "Verific status dbSNP pentru $name..."
    
    echo -e "CHROM\tPOS\tREF\tALT\tSYMBOL\tGene\tExisting_variation\tdbSNP_RS\tCLIN_SIG\tIn_dbSNP\tIn_ClinVar" > "$OUTPUT_DIR/${name}_dbsnp_clinvar.tsv"
    
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/CSQ\n' "$vcf" | \
    python3 -c "
import sys

for line in sys.stdin:
    parts = line.strip().split('\t')
    if len(parts) < 5 or parts[4] == '.':
        continue
    
    chrom, pos, ref, alt, csq = parts
    transcripts = csq.split(',')
    
    selected = None
    
    # Selectează MANE_SELECT (poziția 26 în Python)
    for transcript in transcripts:
        fields = transcript.split('|')
        if len(fields) > 26 and fields[26] and fields[26] != '.':
            selected = transcript
            break
    
    # Dacă nu găsește MANE_SELECT, caută CANONICAL (poziția 24 în Python)
    if not selected:
        for transcript in transcripts:
            fields = transcript.split('|')
            if len(fields) > 24 and fields[24] == 'YES':
                selected = transcript
                break
    
    # Dacă nu găsește nici unul, ia primul
    if not selected and transcripts:
        selected = transcripts[0]
    
    if selected:
        fields = selected.split('|')
        symbol = fields[3] if len(fields) > 3 else ''           # poziția 3
        gene = fields[4] if len(fields) > 4 else ''             # poziția 4
        existing_var = fields[17] if len(fields) > 17 else ''   # poziția 17
        clin_sig = fields[72] if len(fields) > 72 else ''       # poziția 72
        dbsnp_rs = fields[86] if len(fields) > 86 else ''       # poziția 86
        
        # Determină dacă e în dbSNP
        in_dbsnp = 'NO'
        if existing_var and existing_var != '.' and 'rs' in existing_var:
            in_dbsnp = 'YES'
        elif dbsnp_rs and dbsnp_rs != '.':
            in_dbsnp = 'YES'
        
        # Determină dacă e în ClinVar
        in_clinvar = 'NO'
        if clin_sig and clin_sig != '.':
            in_clinvar = 'YES'
        
        print(f'{chrom}\t{pos}\t{ref}\t{alt}\t{symbol}\t{gene}\t{existing_var}\t{dbsnp_rs}\t{clin_sig}\t{in_dbsnp}\t{in_clinvar}')
" >> "$OUTPUT_DIR/${name}_dbsnp_clinvar.tsv"
    
    # Generează statistici
    echo "=== Statistici dbSNP și ClinVar pentru $name ===" > "$OUTPUT_DIR/${name}_dbsnp_clinvar_stats.txt"
    echo "Total variante: $(tail -n +2 "$OUTPUT_DIR/${name}_dbsnp_clinvar.tsv" | wc -l)" >> "$OUTPUT_DIR/${name}_dbsnp_clinvar_stats.txt"
    echo "În dbSNP: $(tail -n +2 "$OUTPUT_DIR/${name}_dbsnp_clinvar.tsv" | cut -f10 | grep -c YES)" >> "$OUTPUT_DIR/${name}_dbsnp_clinvar_stats.txt"
    echo "Nu în dbSNP (novel): $(tail -n +2 "$OUTPUT_DIR/${name}_dbsnp_clinvar.tsv" | cut -f10 | grep -c NO)" >> "$OUTPUT_DIR/${name}_dbsnp_clinvar_stats.txt"
    echo "În ClinVar: $(tail -n +2 "$OUTPUT_DIR/${name}_dbsnp_clinvar.tsv" | cut -f11 | grep -c YES)" >> "$OUTPUT_DIR/${name}_dbsnp_clinvar_stats.txt"
    echo "" >> "$OUTPUT_DIR/${name}_dbsnp_clinvar_stats.txt"
    echo "Semnificații clinice:" >> "$OUTPUT_DIR/${name}_dbsnp_clinvar_stats.txt"
    tail -n +2 "$OUTPUT_DIR/${name}_dbsnp_clinvar.tsv" | cut -f9 | grep -v '^$' | grep -v '^\.$' | sort | uniq -c | sort -nr >> "$OUTPUT_DIR/${name}_dbsnp_clinvar_stats.txt"
    
    echo "Salvat: $OUTPUT_DIR/${name}_dbsnp_clinvar.tsv"
done
