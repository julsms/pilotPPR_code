#!/bin/bash

# Script pentru analiza distribuției variantelor unice pe regiuni genomice
# Autor: Seqera AI
# Data: $(date)

set -euo pipefail

# Directoare
WORK_DIR="/mnt/genomics/pilot_PPR/uniquevariants/genomestratification_distribution"
UNIQUE_DIR="/mnt/genomics/pilot_PPR/uniquevariants"
BED_DIR="/mnt/genomics/GIAB_stratifications/v3.6/with_chr"

# Creează directorul de lucru
mkdir -p "$WORK_DIR"
cd "$WORK_DIR"

# Lista cu fișierele VCF multisample
declare -a VCF_FILES=(
    "uniqueGATKaviti/multisample_uniqueGATKaviti.annotated.vcf.gz"
    "uniqueGATKillumina/multisample_uniqueGATKillumina.annotated.vcf.gz"
    "uniquePPRaviti/multisample_uniquePPRaviti.annotated.vcf.gz"
    "uniquePPRillumina/multisample_uniquePPRillumina.annotated.vcf.gz"
)

# Lista cu regiunile genomice (conform manuscriptului)
declare -a REGIONS=(
    "15GC85"
    "25GC65"
    "ACMG"
    "AD"
    "BB"
    "LM"
    "OD"
    "OMIM"
    "REFSEQ"
    "TRHP"
)

echo "=== ANALIZA DISTRIBUȚIEI VARIANTELOR UNICE PE REGIUNI GENOMICE ==="
echo "Data: $(date)"
echo "Directorul de lucru: $WORK_DIR"
echo ""

# Funcție pentru procesarea unui VCF
process_vcf() {
    local vcf_path="$1"
    local vcf_name="$2"
    
    echo "Procesez: $vcf_name"
    
    # Creează directorul pentru acest VCF
    local output_dir="${vcf_name}_regions"
    mkdir -p "$output_dir"
    
    # Inițializează fișierul de statistici
    local stats_file="${vcf_name}_region_stats.tsv"
    echo -e "Region\tTotal_Variants\tSNVs\tINDELs\tOTHER" > "$stats_file"
    
    # Procesează fiecare regiune
    for region in "${REGIONS[@]}"; do
        local bed_file="$BED_DIR/${region}.bed"
        local output_vcf="${output_dir}/${vcf_name}_${region}.vcf.gz"
        
        if [[ -f "$bed_file" ]]; then
            echo "  - Procesez regiunea: $region"
            
            # Intersectează VCF-ul cu regiunea BED
            bcftools view -R "$bed_file" "$vcf_path" -O z -o "$output_vcf"
            bcftools index "$output_vcf"
            
            # Calculează statistici
            local total_vars=$(bcftools view -H "$output_vcf" | wc -l)
            local snvs=$(bcftools view -H -v snps "$output_vcf" | wc -l)
            local indels=$(bcftools view -H -v indels "$output_vcf" | wc -l)
            local other=$((total_vars - snvs - indels))
            
            # Adaugă în fișierul de statistici
            echo -e "$region\t$total_vars\t$snvs\t$indels\t$other" >> "$stats_file"
            
            echo "    Total variante: $total_vars (SNVs: $snvs, INDELs: $indels, OTHER: $other)"
        else
            echo "  - AVERTISMENT: Fișierul BED nu există: $bed_file"
        fi
    done
    
    echo "  Statistici salvate în: $stats_file"
    echo ""
}

# Procesează fiecare VCF multisample
for vcf_file in "${VCF_FILES[@]}"; do
    vcf_full_path="$UNIQUE_DIR/$vcf_file"
    
    if [[ -f "$vcf_full_path" ]]; then
        # Extrage numele pentru directorul de output
        vcf_basename=$(basename "$vcf_file" .annotated.vcf.gz)
        vcf_basename=$(basename "$vcf_basename" .vcf.gz)
        
        process_vcf "$vcf_full_path" "$vcf_basename"
    else
        echo "EROARE: Fișierul VCF nu există: $vcf_full_path"
    fi
done

# Creează un fișier sumar cu toate statisticile
echo "=== CREEZ FIȘIERUL SUMAR ==="
summary_file="all_regions_summary.tsv"
echo -e "Pipeline\tRegion\tTotal_Variants\tSNVs\tINDELs\tOTHER" > "$summary_file"

for vcf_file in "${VCF_FILES[@]}"; do
    vcf_basename=$(basename "$vcf_file" .annotated.vcf.gz)
    vcf_basename=$(basename "$vcf_basename" .vcf.gz)
    stats_file="${vcf_basename}_region_stats.tsv"
    
    if [[ -f "$stats_file" ]]; then
        # Adaugă numele pipeline-ului la fiecare linie (skip header)
        tail -n +2 "$stats_file" | while IFS=$'\t' read -r region total snvs indels other; do
            echo -e "$vcf_basename\t$region\t$total\t$snvs\t$indels\t$other" >> "$summary_file"
        done
    fi
done

echo "Fișier sumar creat: $summary_file"

# Creează un script Python pentru verificarea rezultatelor
cat > "check_results.py" << 'EOF'
#!/usr/bin/env python3

import pandas as pd
import os

def check_results():
    """Verifică rezultatele analizei"""
    
    # Citește fișierul sumar
    if os.path.exists('all_regions_summary.tsv'):
        df = pd.read_csv('all_regions_summary.tsv', sep='\t')
        
        print("=== VERIFICARE REZULTATE ===")
        print(f"Total înregistrări: {len(df)}")
        print(f"Pipeline-uri: {df['Pipeline'].unique()}")
        print(f"Regiuni: {df['Region'].unique()}")
        print()
        
        # Sumar pe pipeline
        print("=== SUMAR PE PIPELINE ===")
        pipeline_summary = df.groupby('Pipeline')['Total_Variants'].sum().sort_values(ascending=False)
        for pipeline, total in pipeline_summary.items():
            print(f"{pipeline}: {total:,} variante")
        print()
        
        # Sumar pe regiuni
        print("=== SUMAR PE REGIUNI ===")
        region_summary = df.groupby('Region')['Total_Variants'].sum().sort_values(ascending=False)
        for region, total in region_summary.items():
            print(f"{region}: {total:,} variante")
        print()
        
        # Top 5 combinații pipeline-regiune
        print("=== TOP 5 COMBINAȚII PIPELINE-REGIUNE ===")
        top_combinations = df.nlargest(5, 'Total_Variants')[['Pipeline', 'Region', 'Total_Variants']]
        for _, row in top_combinations.iterrows():
            print(f"{row['Pipeline']} - {row['Region']}: {row['Total_Variants']:,} variante")
            
    else:
        print("EROARE: Fișierul all_regions_summary.tsv nu există!")

if __name__ == "__main__":
    check_results()
EOF

chmod +x check_results.py

echo ""
echo "=== ANALIZA COMPLETĂ ==="
echo "Rezultatele sunt salvate în: $WORK_DIR"
echo "Fișiere create:"
echo "  - Directoare cu VCF-uri pe regiuni pentru fiecare pipeline"
echo "  - Fișiere de statistici individuale (*_region_stats.tsv)"
echo "  - Fișier sumar: $summary_file"
echo "  - Script de verificare: check_results.py"
echo ""
echo "Pentru a verifica rezultatele, rulează:"
echo "  cd $WORK_DIR"
echo "  python3 check_results.py"
echo ""
echo "Pentru a vedea structura directoarelor:"
echo "  tree $WORK_DIR"
