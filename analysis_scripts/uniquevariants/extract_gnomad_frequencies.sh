#!/bin/bash
# extract_gnomad_frequencies.sh

BASE_DIR="/mnt/genomics/pilot_PPR/uniquevariants"
OUTPUT_DIR="/mnt/genomics/pilot_PPR/uniquevariants/extract_csq_tsv"

mkdir -p "$OUTPUT_DIR"

VCFS=(
    "$BASE_DIR/uniqueGATKillumina/multisample_uniqueGATKillumina.annotated.vcf.gz"
    "$BASE_DIR/uniqueGATKaviti/multisample_uniqueGATKaviti.annotated.vcf.gz"
    "$BASE_DIR/uniquePPRillumina/multisample_uniquePPRillumina.annotated.vcf.gz"
    "$BASE_DIR/uniquePPRaviti/multisample_uniquePPRaviti.annotated.vcf.gz"
)

NAMES=(
    "uniqueGATKillumina"
    "uniqueGATKaviti"
    "uniquePPRillumina"
    "uniquePPRaviti"
)

for i in "${!VCFS[@]}"; do
    vcf="${VCFS[$i]}"
    name="${NAMES[$i]}"

    echo "Extrag frecvențe populaționale pentru $name..."

    echo -e "CHROM\tPOS\tREF\tALT\tSYMBOL\tGene\t1000G_AF\t1000G_AFR_AF\t1000G_AMR_AF\t1000G_EAS_AF\t1000G_EUR_AF\t1000G_SAS_AF\tgnomADe_AF\tgnomADe_AFR_AF\tgnomADe_AMR_AF\tgnomADe_ASJ_AF\tgnomADe_EAS_AF\tgnomADe_FIN_AF\tgnomADe_MID_AF\tgnomADe_NFE_AF\tgnomADe_REMAINING_AF\tgnomADe_SAS_AF\tgnomADg_AF\tgnomADg_AFR_AF\tgnomADg_AMI_AF\tgnomADg_AMR_AF\tgnomADg_ASJ_AF\tgnomADg_EAS_AF\tgnomADg_FIN_AF\tgnomADg_MID_AF\tgnomADg_NFE_AF\tgnomADg_REMAINING_AF\tgnomADg_SAS_AF\tMAX_AF\tMAX_AF_POPS\tIn_1000G\tIn_gnomAD" \
    > "$OUTPUT_DIR/${name}_frequencies.tsv"

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

        symbol = fields[3] if len(fields) > 3 else ''
        gene = fields[4] if len(fields) > 4 else ''

        # 1000 Genomes (pozițiile 43-48)
        af_1000g = fields[43] if len(fields) > 43 else ''
        afr_af_1000g = fields[44] if len(fields) > 44 else ''
        amr_af_1000g = fields[45] if len(fields) > 45 else ''
        eas_af_1000g = fields[46] if len(fields) > 46 else ''
        eur_af_1000g = fields[47] if len(fields) > 47 else ''
        sas_af_1000g = fields[48] if len(fields) > 48 else ''

        # gnomAD Exomes (pozițiile 49-58)
        gnomade_af = fields[49] if len(fields) > 49 else ''
        gnomade_afr_af = fields[50] if len(fields) > 50 else ''
        gnomade_amr_af = fields[51] if len(fields) > 51 else ''
        gnomade_asj_af = fields[52] if len(fields) > 52 else ''
        gnomade_eas_af = fields[53] if len(fields) > 53 else ''
        gnomade_fin_af = fields[54] if len(fields) > 54 else ''
        gnomade_mid_af = fields[55] if len(fields) > 55 else ''
        gnomade_nfe_af = fields[56] if len(fields) > 56 else ''
        gnomade_remaining_af = fields[57] if len(fields) > 57 else ''
        gnomade_sas_af = fields[58] if len(fields) > 58 else ''

        # gnomAD Genomes (pozițiile 59-69)
        gnomadg_af = fields[59] if len(fields) > 59 else ''
        gnomadg_afr_af = fields[60] if len(fields) > 60 else ''
        gnomadg_ami_af = fields[61] if len(fields) > 61 else ''
        gnomadg_amr_af = fields[62] if len(fields) > 62 else ''
        gnomadg_asj_af = fields[63] if len(fields) > 63 else ''
        gnomadg_eas_af = fields[64] if len(fields) > 64 else ''
        gnomadg_fin_af = fields[65] if len(fields) > 65 else ''
        gnomadg_mid_af = fields[66] if len(fields) > 66 else ''
        gnomadg_nfe_af = fields[67] if len(fields) > 67 else ''
        gnomadg_remaining_af = fields[68] if len(fields) > 68 else ''
        gnomadg_sas_af = fields[69] if len(fields) > 69 else ''

        # MAX AF (pozițiile 70-71)
        max_af = fields[70] if len(fields) > 70 else ''
        max_af_pops = fields[71] if len(fields) > 71 else ''

        # Determină dacă e în 1000G
        in_1000g = 'NO'
        try:
            freqs_1000g = [
                af_1000g,
                afr_af_1000g,
                amr_af_1000g,
                eas_af_1000g,
                eur_af_1000g,
                sas_af_1000g
            ]

            if any(
                float(freq) > 0
                for freq in freqs_1000g
                if freq and freq != '.'
            ):
                in_1000g = 'YES'
        except:
            pass

        # Determină dacă e în gnomAD
        in_gnomad = 'NO'
        try:
            gnomad_freqs = [gnomade_af, gnomadg_af]

            if any(
                float(freq) > 0
                for freq in gnomad_freqs
                if freq and freq != '.'
            ):
                in_gnomad = 'YES'
        except:
            pass

        print(
            f'{chrom}\t{pos}\t{ref}\t{alt}\t{symbol}\t{gene}\t'
            f'{af_1000g}\t{afr_af_1000g}\t{amr_af_1000g}\t'
            f'{eas_af_1000g}\t{eur_af_1000g}\t{sas_af_1000g}\t'
            f'{gnomade_af}\t{gnomade_afr_af}\t{gnomade_amr_af}\t'
            f'{gnomade_asj_af}\t{gnomade_eas_af}\t{gnomade_fin_af}\t'
            f'{gnomade_mid_af}\t{gnomade_nfe_af}\t'
            f'{gnomade_remaining_af}\t{gnomade_sas_af}\t'
            f'{gnomadg_af}\t{gnomadg_afr_af}\t{gnomadg_ami_af}\t'
            f'{gnomadg_amr_af}\t{gnomadg_asj_af}\t'
            f'{gnomadg_eas_af}\t{gnomadg_fin_af}\t'
            f'{gnomadg_mid_af}\t{gnomadg_nfe_af}\t'
            f'{gnomadg_remaining_af}\t{gnomadg_sas_af}\t'
            f'{max_af}\t{max_af_pops}\t'
            f'{in_1000g}\t{in_gnomad}'
        )
" >> "$OUTPUT_DIR/${name}_frequencies.tsv"

    # Generează statistici
    echo "=== Statistici frecvențe pentru $name ===" \
    > "$OUTPUT_DIR/${name}_frequency_stats.txt"

    echo "Total variante: $(tail -n +2 "$OUTPUT_DIR/${name}_frequencies.tsv" | wc -l)" \
    >> "$OUTPUT_DIR/${name}_frequency_stats.txt"

    echo "În 1000G: $(tail -n +2 "$OUTPUT_DIR/${name}_frequencies.tsv" | cut -f36 | grep -c YES)" \
    >> "$OUTPUT_DIR/${name}_frequency_stats.txt"

    echo "În gnomAD: $(tail -n +2 "$OUTPUT_DIR/${name}_frequencies.tsv" | cut -f37 | grep -c YES)" \
    >> "$OUTPUT_DIR/${name}_frequency_stats.txt"

    echo "Novel (nu în 1000G și nu în gnomAD): $(tail -n +2 "$OUTPUT_DIR/${name}_frequencies.tsv" | awk -F'\t' '\$36==\"NO\" && \$37==\"NO\"' | wc -l)" \
    >> "$OUTPUT_DIR/${name}_frequency_stats.txt"

    echo "Salvat: $OUTPUT_DIR/${name}_frequencies.tsv"
done
