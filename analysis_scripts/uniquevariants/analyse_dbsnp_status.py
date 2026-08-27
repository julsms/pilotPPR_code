#!/usr/bin/env python3
# analyze_dbsnp_status.py

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def process_dbsnp_data(input_dir, output_dir):
    """Process dbSNP data and create categorized files."""
    
    # Create output directory
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    # Define input files
    files = {
        'GATK_Illumina': f'{input_dir}/uniqueGATKillumina_dbsnp_clinvar.tsv',
        'PPR_Illumina': f'{input_dir}/uniquePPRillumina_dbsnp_clinvar.tsv', 
        'GATK_Aviti': f'{input_dir}/uniqueGATKaviti_dbsnp_clinvar.tsv',
        'PPR_Aviti': f'{input_dir}/uniquePPRaviti_dbsnp_clinvar.tsv'
    }
    
    all_stats = {}
    
    # Process each file
    for dataset_name, file_path in files.items():
        print(f"Processing dbSNP data for {dataset_name}...")
        
        # Read data
        df = pd.read_csv(file_path, sep='\t')
        
        # Determine dbSNP status based on rsID presence
        # A variant is in dbSNP if it has an rsID (rs followed by numbers)
        df['has_rsid'] = False
        
        # Check Existing_variation column for rsID
        existing_var_mask = df['Existing_variation'].notna() & (df['Existing_variation'] != '.') & (df['Existing_variation'] != '')
        df.loc[existing_var_mask, 'has_rsid'] = df.loc[existing_var_mask, 'Existing_variation'].str.contains(r'rs\d+', na=False)
        
        # Check dbSNP_RS column for rsID (backup check)
        dbsnp_rs_mask = df['dbSNP_RS'].notna() & (df['dbSNP_RS'] != '.') & (df['dbSNP_RS'] != '')
        df.loc[dbsnp_rs_mask & ~df['has_rsid'], 'has_rsid'] = df.loc[dbsnp_rs_mask & ~df['has_rsid'], 'dbSNP_RS'].str.contains(r'rs\d+', na=False)
        
        # Also check the In_dbSNP column from your original script
        in_dbsnp_mask = df['In_dbSNP'] == 'YES'
        df.loc[in_dbsnp_mask, 'has_rsid'] = True
        
        # Base columns to save
        base_cols = ['CHROM', 'POS', 'REF', 'ALT', 'SYMBOL', 'Gene']
        dbsnp_cols = ['Existing_variation', 'dbSNP_RS', 'CLIN_SIG', 'In_dbSNP', 'In_ClinVar']
        
        # Variants in dbSNP (have rsID)
        in_dbsnp_variants = df[df['has_rsid']][base_cols + dbsnp_cols]
        in_dbsnp_variants.to_csv(f'{output_dir}/{dataset_name}_in_dbSNP.csv', index=False)
        
        # Variants NOT in dbSNP (no rsID)
        not_in_dbsnp_variants = df[~df['has_rsid']][base_cols + dbsnp_cols]
        not_in_dbsnp_variants.to_csv(f'{output_dir}/{dataset_name}_NOT_in_dbSNP.csv', index=False)
        
        # Variants in ClinVar
        in_clinvar_variants = df[df['In_ClinVar'] == 'YES'][base_cols + dbsnp_cols]
        in_clinvar_variants.to_csv(f'{output_dir}/{dataset_name}_in_ClinVar.csv', index=False)
        
        # Collect statistics
        stats = {
            'total_variants': len(df),
            'in_dbSNP': df['has_rsid'].sum(),
            'not_in_dbSNP': (~df['has_rsid']).sum(),
            'in_ClinVar': (df['In_ClinVar'] == 'YES').sum(),
            'not_in_ClinVar': (df['In_ClinVar'] != 'YES').sum(),
            'in_both_dbSNP_ClinVar': (df['has_rsid'] & (df['In_ClinVar'] == 'YES')).sum()
        }
        
        all_stats[dataset_name] = stats
        
        print(f"  Total variants: {stats['total_variants']:,}")
        print(f"  In dbSNP (have rsID): {stats['in_dbSNP']:,} ({stats['in_dbSNP']/stats['total_variants']*100:.1f}%)")
        print(f"  NOT in dbSNP (no rsID): {stats['not_in_dbSNP']:,} ({stats['not_in_dbSNP']/stats['total_variants']*100:.1f}%)")
        print(f"  In ClinVar: {stats['in_ClinVar']:,} ({stats['in_ClinVar']/stats['total_variants']*100:.1f}%)")
        print()
    
    # Create summary statistics
    create_dbsnp_summary_stats(all_stats, output_dir)
    
    # Create comparison plots
    create_dbsnp_comparison_plots(all_stats, output_dir)
    
    print(f"✓ dbSNP analysis complete! Results saved in {output_dir}")

def create_dbsnp_summary_stats(all_stats, output_dir):
    """Create summary statistics for dbSNP analysis."""
    
    # Convert to DataFrame
    stats_df = pd.DataFrame(all_stats).T
    
    # Add percentages
    for col in ['in_dbSNP', 'not_in_dbSNP', 'in_ClinVar', 'not_in_ClinVar', 'in_both_dbSNP_ClinVar']:
        stats_df[f'{col}_pct'] = (stats_df[col] / stats_df['total_variants'] * 100).round(1)
    
    # Save detailed stats
    stats_df.to_csv(f'{output_dir}/dbsnp_summary_statistics.csv')
    
    # Create readable summary
    with open(f'{output_dir}/dbsnp_summary_report.txt', 'w') as f:
        f.write("dbSNP AND ClinVar ANALYSIS SUMMARY\n")
        f.write("=" * 50 + "\n\n")
        
        for dataset in stats_df.index:
            f.write(f"{dataset.replace('_', ' ')}:\n")
            f.write(f"  Total variants: {stats_df.loc[dataset, 'total_variants']:,}\n")
            f.write(f"  In dbSNP (have rsID): {stats_df.loc[dataset, 'in_dbSNP']:,} ({stats_df.loc[dataset, 'in_dbSNP_pct']:.1f}%)\n")
            f.write(f"  NOT in dbSNP (no rsID): {stats_df.loc[dataset, 'not_in_dbSNP']:,} ({stats_df.loc[dataset, 'not_in_dbSNP_pct']:.1f}%)\n")
            f.write(f"  In ClinVar: {stats_df.loc[dataset, 'in_ClinVar']:,} ({stats_df.loc[dataset, 'in_ClinVar_pct']:.1f}%)\n")
            f.write(f"  In both dbSNP & ClinVar: {stats_df.loc[dataset, 'in_both_dbSNP_ClinVar']:,} ({stats_df.loc[dataset, 'in_both_dbSNP_ClinVar_pct']:.1f}%)\n")
            f.write("\n")
    
    print("✓ Created: dbsnp_summary_statistics.csv and dbsnp_summary_report.txt")

def create_dbsnp_comparison_plots(all_stats, output_dir):
    """Create comparison bar plots for dbSNP data."""
    
    # Colors for each dataset
    colors = {
        'GATK_Illumina': '#0173B2',
        'PPR_Illumina': '#DE8F05', 
        'GATK_Aviti': '#029E73',
        'PPR_Aviti': '#CC78BC'
    }
    
    # Prepare data for plotting
    datasets = list(all_stats.keys())
    datasets_clean = [d.replace('_', ' ') for d in datasets]
    
    # Calculate percentages
    in_dbsnp_pct = [all_stats[d]['in_dbSNP'] / all_stats[d]['total_variants'] * 100 for d in datasets]
    not_in_dbsnp_pct = [all_stats[d]['not_in_dbSNP'] / all_stats[d]['total_variants'] * 100 for d in datasets]
    in_clinvar_pct = [all_stats[d]['in_ClinVar'] / all_stats[d]['total_variants'] * 100 for d in datasets]
    
    # Create comparison plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 8))
    
    # Plot 1: dbSNP status
    x = np.arange(len(datasets))
    width = 0.35
    
    bars1 = ax1.bar(x - width/2, in_dbsnp_pct, width, label='In dbSNP (have rsID)', color='#2E7D32', alpha=0.8)
    bars2 = ax1.bar(x + width/2, not_in_dbsnp_pct, width, label='NOT in dbSNP (no rsID)', color='#D32F2F', alpha=0.8)
    
    ax1.set_xlabel('Dataset', fontsize=18, fontweight='bold')
    ax1.set_ylabel('Percentage of Variants (%)', fontsize=18, fontweight='bold')
    ax1.set_title('dbSNP Status of Unique Variants', fontsize=20, fontweight='bold', pad=20)
    ax1.set_xticks(x)
    ax1.set_xticklabels(datasets_clean, fontsize=16)
    ax1.legend(fontsize=16)
    ax1.tick_params(axis='y', labelsize=16)
    
    # Add value labels on bars
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            ax1.annotate(f'{height:.1f}%', xy=(bar.get_x() + bar.get_width() / 2, height),
                        xytext=(0, 3), textcoords="offset points", ha='center', va='bottom', fontsize=14)
    
    # Plot 2: ClinVar status
    bars3 = ax2.bar(datasets_clean, in_clinvar_pct, color=[colors[d] for d in datasets], alpha=0.8)
    
    ax2.set_xlabel('Dataset', fontsize=18, fontweight='bold')
    ax2.set_ylabel('Percentage of Variants (%)', fontsize=18, fontweight='bold')
    ax2.set_title('Variants with Clinical Significance (ClinVar)', fontsize=20, fontweight='bold', pad=20)
    ax2.tick_params(axis='both', labelsize=16)
    
    # Add value labels
    for bar in bars3:
        height = bar.get_height()
        ax2.annotate(f'{height:.1f}%', xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, 3), textcoords="offset points", ha='center', va='bottom', fontsize=14)
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/dbsnp_comparison.png", dpi=300, bbox_inches='tight')
    plt.close()
    print("✓ Created: dbsnp_comparison.png")

if __name__ == "__main__":
    input_dir = "/mnt/genomics/pilot_PPR/uniquevariants/extract_csq_tsv"
    output_dir = "/mnt/genomics/pilot_PPR/uniquevariants/extract_csq_tsv/plots_frequencies"
    
    process_dbsnp_data(input_dir, output_dir)
