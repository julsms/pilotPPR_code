#!/bin/bash
# extract_mane_canonical_all.sh

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
    
    echo "Procesez $name..."
    
    # Extrage header CSQ
    header=$(bcftools view -h "$vcf" | grep "##INFO=<ID=CSQ" | sed 's/.*Format: //' | sed 's/">$//')
    echo "$header" | tr '|' '\t' > "$OUTPUT_DIR/${name}_mane_canonical_all.tsv"
    echo -e "CHROM\tPOS\tREF\tALT\t$(cat "$OUTPUT_DIR/${name}_mane_canonical_all.tsv")" > "$OUTPUT_DIR/${name}_mane_canonical_all.tsv"
    
    # Extrage toate anotările pentru MANE Select sau Canonical
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
    
    # Caută MANE_SELECT (poziția 26 în Python)
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
        print(f'{chrom}\t{pos}\t{ref}\t{alt}\t{selected.replace(\"|\", chr(9))}')
" >> "$OUTPUT_DIR/${name}_mane_canonical_all.tsv"
    
    echo "Salvat: $OUTPUT_DIR/${name}_mane_canonical_all.tsv"
done
