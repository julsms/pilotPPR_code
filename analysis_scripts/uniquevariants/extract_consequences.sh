#!/bin/bash
# extract_consequences.sh

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
    
    echo "Extrag consequence types pentru $name..."
    
    echo -e "CHROM\tPOS\tREF\tALT\tConsequence\tIMPACT\tSYMBOL\tGene\tFeature\tBIOTYPE\tCANONICAL\tMANE_SELECT" > "$OUTPUT_DIR/${name}_consequences.tsv"
    
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
        consequence = fields[1] if len(fields) > 1 else ''      # poziția 1
        impact = fields[2] if len(fields) > 2 else ''           # poziția 2
        symbol = fields[3] if len(fields) > 3 else ''           # poziția 3
        gene = fields[4] if len(fields) > 4 else ''             # poziția 4
        feature = fields[6] if len(fields) > 6 else ''          # poziția 6
        biotype = fields[7] if len(fields) > 7 else ''          # poziția 7
        canonical = fields[24] if len(fields) > 24 else ''      # poziția 24
        mane_select = fields[26] if len(fields) > 26 else ''    # poziția 26
        
        print(f'{chrom}\t{pos}\t{ref}\t{alt}\t{consequence}\t{impact}\t{symbol}\t{gene}\t{feature}\t{biotype}\t{canonical}\t{mane_select}')
" >> "$OUTPUT_DIR/${name}_consequences.tsv"
    
    # Generează statistici consequence types
    echo "=== Statistici consequence types pentru $name ===" > "$OUTPUT_DIR/${name}_consequence_stats.txt"
    echo "Total variante: $(tail -n +2 "$OUTPUT_DIR/${name}_consequences.tsv" | wc -l)" >> "$OUTPUT_DIR/${name}_consequence_stats.txt"
    echo "" >> "$OUTPUT_DIR/${name}_consequence_stats.txt"
    echo "Top consequence types:" >> "$OUTPUT_DIR/${name}_consequence_stats.txt"
    tail -n +2 "$OUTPUT_DIR/${name}_consequences.tsv" | cut -f5 | sort | uniq -c | sort -nr >> "$OUTPUT_DIR/${name}_consequence_stats.txt"
    
    echo "Salvat: $OUTPUT_DIR/${name}_consequences.tsv"
done
