#!/usr/bin/env python3
# analyze_allele_frequencies.py

import pandas as pd
import numpy as np
from pathlib import Path

def analyze_allele_frequencies(input_dir, output_dir):
    """Analyze allele frequency distributions for pipeline-specific variants."""
    
    # Create output directory
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    # Define input files
    files = {
        'GATK_Illumina': f'{input_dir}/uniqueGATKillumina_frequencies.tsv',
        'PPR_Illumina': f'{input_dir}/uniquePPRillumina_frequencies.tsv', 
        'GATK_Aviti': f'{input_dir}/uniqueGATKaviti_frequencies.tsv',
        'PPR_Aviti': f'{input_dir}/uniquePPRaviti_frequencies.tsv'
    }
    
    all_results = []
    
    print("=== ALLELE FREQUENCY ANALYSIS ===\n")
    
    for dataset_name, file_path in files.items():
        print(f"Processing {dataset_name}...")
        
        # Read data
        df = pd.read_csv(file_path, sep='\t')
        
        # Clean and convert frequency columns
        gnomad_af = pd.to_numeric(df['gnomADg_AF'].replace(['', '.', 'NA'], np.nan), errors='coerce')
        kg_af = pd.to_numeric(df['1000G_AF'].replace(['', '.', 'NA'], np.nan), errors='coerce')
        
        # Remove NaN and zero values
        gnomad_af_clean = gnomad_af[(gnomad_af > 0) & (~gnomad_af.isna())]
        kg_af_clean = kg_af[(kg_af > 0) & (~kg_af.isna())]
        
        # Analyze gnomAD frequencies
        if len(gnomad_af_clean) > 0:
            gnomad_stats = {
                'Dataset': dataset_name,
                'Database': 'gnomAD_Genomes',
                'Total_variants_with_AF': len(gnomad_af_clean),
                'Median_AF': np.median(gnomad_af_clean),
                'Mean_AF': np.mean(gnomad_af_clean),
                'Q25_AF': np.percentile(gnomad_af_clean, 25),
                'Q75_AF': np.percentile(gnomad_af_clean, 75),
                'Rare_variants_count': np.sum(gnomad_af_clean < 0.01),
                'Rare_variants_pct': (np.sum(gnomad_af_clean < 0.01) / len(gnomad_af_clean)) * 100,
                'Common_variants_count': np.sum(gnomad_af_clean > 0.05),
                'Common_variants_pct': (np.sum(gnomad_af_clean > 0.05) / len(gnomad_af_clean)) * 100,
                'Very_rare_count': np.sum(gnomad_af_clean < 0.001),
                'Very_rare_pct': (np.sum(gnomad_af_clean < 0.001) / len(gnomad_af_clean)) * 100,
                'Min_AF': np.min(gnomad_af_clean),
                'Max_AF': np.max(gnomad_af_clean)
            }
            all_results.append(gnomad_stats)
            
            print(f"  gnomAD Genomes:")
            print(f"    Variants with AF: {len(gnomad_af_clean):,}")
            print(f"    Median AF: {np.median(gnomad_af_clean):.6f}")
            print(f"    Rare variants (AF < 0.01): {np.sum(gnomad_af_clean < 0.01):,} ({(np.sum(gnomad_af_clean < 0.01) / len(gnomad_af_clean)) * 100:.1f}%)")
            print(f"    Common variants (AF > 0.05): {np.sum(gnomad_af_clean > 0.05):,} ({(np.sum(gnomad_af_clean > 0.05) / len(gnomad_af_clean)) * 100:.1f}%)")
        
        # Analyze 1000G frequencies
        if len(kg_af_clean) > 0:
            kg_stats = {
                'Dataset': dataset_name,
                'Database': '1000_Genomes',
                'Total_variants_with_AF': len(kg_af_clean),
                'Median_AF': np.median(kg_af_clean),
                'Mean_AF': np.mean(kg_af_clean),
                'Q25_AF': np.percentile(kg_af_clean, 25),
                'Q75_AF': np.percentile(kg_af_clean, 75),
                'Rare_variants_count': np.sum(kg_af_clean < 0.01),
                'Rare_variants_pct': (np.sum(kg_af_clean < 0.01) / len(kg_af_clean)) * 100,
                'Common_variants_count': np.sum(kg_af_clean > 0.05),
                'Common_variants_pct': (np.sum(kg_af_clean > 0.05) / len(kg_af_clean)) * 100,
                'Very_rare_count': np.sum(kg_af_clean < 0.001),
                'Very_rare_pct': (np.sum(kg_af_clean < 0.001) / len(kg_af_clean)) * 100,
                'Min_AF': np.min(kg_af_clean),
                'Max_AF': np.max(kg_af_clean)
            }
            all_results.append(kg_stats)
            
            print(f"  1000 Genomes:")
            print(f"    Variants with AF: {len(kg_af_clean):,}")
            print(f"    Median AF: {np.median(kg_af_clean):.6f}")
            print(f"    Rare variants (AF < 0.01): {np.sum(kg_af_clean < 0.01):,} ({(np.sum(kg_af_clean < 0.01) / len(kg_af_clean)) * 100:.1f}%)")
            print(f"    Common variants (AF > 0.05): {np.sum(kg_af_clean > 0.05):,} ({(np.sum(kg_af_clean > 0.05) / len(kg_af_clean)) * 100:.1f}%)")
        
        print()
    
    # Convert to DataFrame and save
    results_df = pd.DataFrame(all_results)
    results_df.to_csv(f'{output_dir}/allele_frequency_statistics.csv', index=False)
    
    # Create summary report
    create_summary_report(results_df, output_dir)
    
    print(f"✓ Analysis complete! Results saved in {output_dir}")
    print(f"✓ Created: allele_frequency_statistics.csv")
    print(f"✓ Created: allele_frequency_summary.txt")

def create_summary_report(results_df, output_dir):
    """Create a readable summary report."""
    
    with open(f'{output_dir}/allele_frequency_summary.txt', 'w') as f:
        f.write("ALLELE FREQUENCY ANALYSIS SUMMARY\n")
        f.write("=" * 50 + "\n\n")
        
        # Group by dataset
        datasets = results_df['Dataset'].unique()
        
        for dataset in datasets:
            f.write(f"{dataset.replace('_', ' ')} Specific Variants:\n")
            f.write("-" * 40 + "\n")
            
            dataset_data = results_df[results_df['Dataset'] == dataset]
            
            for _, row in dataset_data.iterrows():
                f.write(f"\n{row['Database']}:\n")
                f.write(f"  Total variants with AF: {row['Total_variants_with_AF']:,}\n")
                f.write(f"  Median AF: {row['Median_AF']:.6f}\n")
                f.write(f"  Mean AF: {row['Mean_AF']:.6f}\n")
                f.write(f"  Q25-Q75 AF: {row['Q25_AF']:.6f} - {row['Q75_AF']:.6f}\n")
                f.write(f"  AF range: {row['Min_AF']:.6f} - {row['Max_AF']:.6f}\n")
                f.write(f"  Rare variants (AF < 0.01): {row['Rare_variants_count']:,} ({row['Rare_variants_pct']:.1f}%)\n")
                f.write(f"  Very rare variants (AF < 0.001): {row['Very_rare_count']:,} ({row['Very_rare_pct']:.1f}%)\n")
                f.write(f"  Common variants (AF > 0.05): {row['Common_variants_count']:,} ({row['Common_variants_pct']:.1f}%)\n")
            
            f.write("\n")
        
        # Comparative analysis
        f.write("COMPARATIVE ANALYSIS:\n")
        f.write("=" * 30 + "\n\n")
        
        # Compare median frequencies
        f.write("Median Allele Frequencies:\n")
        for dataset in datasets:
            dataset_data = results_df[results_df['Dataset'] == dataset]
            f.write(f"\n{dataset.replace('_', ' ')}:\n")
            for _, row in dataset_data.iterrows():
                f.write(f"  {row['Database']}: {row['Median_AF']:.6f}\n")
        
        f.write("\n")
        
        # Compare rare variant proportions
        f.write("Rare Variants (AF < 0.01) Proportions:\n")
        for dataset in datasets:
            dataset_data = results_df[results_df['Dataset'] == dataset]
            f.write(f"\n{dataset.replace('_', ' ')}:\n")
            for _, row in dataset_data.iterrows():
                f.write(f"  {row['Database']}: {row['Rare_variants_pct']:.1f}%\n")
        
        f.write("\n")
        
        # Compare common variant proportions
        f.write("Common Variants (AF > 0.05) Proportions:\n")
        for dataset in datasets:
            dataset_data = results_df[results_df['Dataset'] == dataset]
            f.write(f"\n{dataset.replace('_', ' ')}:\n")
            for _, row in dataset_data.iterrows():
                f.write(f"  {row['Database']}: {row['Common_variants_pct']:.1f}%\n")

def create_comparison_table(results_df, output_dir):
    """Create a comparison table for manuscript."""
    
    # Pivot table for easier comparison
    comparison_data = []
    
    datasets = results_df['Dataset'].unique()
    databases = results_df['Database'].unique()
    
    for dataset in datasets:
        for db in databases:
            row_data = results_df[(results_df['Dataset'] == dataset) & (results_df['Database'] == db)]
            if not row_data.empty:
                row = row_data.iloc[0]
                comparison_data.append({
                    'Pipeline_Platform': dataset.replace('_', ' '),
                    'Database': db.replace('_', ' '),
                    'Median_AF': f"{row['Median_AF']:.6f}",
                    'Rare_Pct': f"{row['Rare_variants_pct']:.1f}%",
                    'Common_Pct': f"{row['Common_variants_pct']:.1f}%",
                    'Total_Variants': f"{row['Total_variants_with_AF']:,}"
                })
    
    comparison_df = pd.DataFrame(comparison_data)
    comparison_df.to_csv(f'{output_dir}/allele_frequency_comparison.csv', index=False)
    print(f"✓ Created: allele_frequency_comparison.csv")

if __name__ == "__main__":
    input_dir = "/mnt/genomics/pilot_PPR/uniquevariants/extract_csq_tsv"
    output_dir = "/mnt/genomics/pilot_PPR/uniquevariants/extract_csq_tsv/plots_frequencies"
    
    analyze_allele_frequencies(input_dir, output_dir)
