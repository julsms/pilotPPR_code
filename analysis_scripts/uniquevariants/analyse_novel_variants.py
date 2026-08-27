#!/usr/bin/env python3
# analyze_novel_variants_final.py

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import matplotlib.patches as patches
from matplotlib_venn import venn3, venn3_circles
from upsetplot import UpSet
from itertools import combinations
from matplotlib.patches import Circle
import warnings
warnings.filterwarnings('ignore')

# Autosomes only - exclude sex chromosomes to avoid false novel calls
AUTOSOMES = [f'chr{i}' for i in range(1, 23)]


def filter_autosomes(df):
    """Filter DataFrame to autosomal chromosomes only."""
    return df[df['CHROM'].isin(AUTOSOMES)].copy()


def create_variant_id(df):
    """Create unique variant ID: CHROM:POS:REF:ALT"""
    return df['CHROM'].astype(str) + ':' + df['POS'].astype(str) + ':' + df['REF'] + ':' + df['ALT']


def analyze_novel_variants(base_dir, output_dir):
    """Comprehensive analysis of novel variants across all databases."""
    
    # Create output directory
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    # Define datasets
    datasets = ['GATK_Illumina', 'PPR_Illumina', 'GATK_Aviti', 'PPR_Aviti']
    
    # Load all data
    all_data = {}
    novel_variants = {}
    database_stats = {}
    
    for dataset in datasets:
        print(f"Processing {dataset}...")
        
        # Load data from different sources
        gnomad_file = f"{base_dir}/{dataset}_in_gnomAD.csv"
        not_gnomad_file = f"{base_dir}/{dataset}_NOT_in_gnomAD.csv"
        kg_file = f"{base_dir}/{dataset}_in_1000G.csv"
        not_kg_file = f"{base_dir}/{dataset}_NOT_in_1000G.csv"
        dbsnp_file = f"{base_dir}/{dataset}_in_dbSNP.csv"
        not_dbsnp_file = f"{base_dir}/{dataset}_NOT_in_dbSNP.csv"
        
        # Read all files
        in_gnomad = pd.read_csv(gnomad_file)
        not_in_gnomad = pd.read_csv(not_gnomad_file)
        in_1000g = pd.read_csv(kg_file)
        not_in_1000g = pd.read_csv(not_kg_file)
        in_dbsnp = pd.read_csv(dbsnp_file)
        not_in_dbsnp = pd.read_csv(not_dbsnp_file)
        
        # Filter to autosomes only for novel variant analysis
        in_gnomad = filter_autosomes(in_gnomad)
        not_in_gnomad = filter_autosomes(not_in_gnomad)
        in_1000g = filter_autosomes(in_1000g)
        not_in_1000g = filter_autosomes(not_in_1000g)
        in_dbsnp = filter_autosomes(in_dbsnp)
        not_in_dbsnp = filter_autosomes(not_in_dbsnp)
        
        # Calculate total variants (autosomes only)
        total_variants = len(in_gnomad) + len(not_in_gnomad)
        
        # Create variant IDs for intersection analysis
        not_gnomad_ids = set(create_variant_id(not_in_gnomad))
        not_1000g_ids = set(create_variant_id(not_in_1000g))
        not_dbsnp_ids = set(create_variant_id(not_in_dbsnp))
        
        in_gnomad_ids = set(create_variant_id(in_gnomad))
        in_1000g_ids = set(create_variant_id(in_1000g))
        in_dbsnp_ids = set(create_variant_id(in_dbsnp))
        
        # Find completely novel variants
        # NOVEL = NOT in gnomAD AND NOT in 1000G AND NOT in dbSNP
        novel_ids = not_gnomad_ids & not_1000g_ids & not_dbsnp_ids
        
        # Verify the intersection logic
        print(f"  [Autosomes only] Total variants: {total_variants:,}")
        print(f"  NOT in gnomAD: {len(not_gnomad_ids):,}")
        print(f"  NOT in 1000G: {len(not_1000g_ids):,}")
        print(f"  NOT in dbSNP: {len(not_dbsnp_ids):,}")
        print(f"  Intersection (Novel): {len(novel_ids):,}")
        
        # Create DataFrame with novel variants
        all_variants = pd.concat([in_gnomad, not_in_gnomad], ignore_index=True)
        all_variants['variant_id'] = create_variant_id(all_variants)
        novel_df = all_variants[all_variants['variant_id'].isin(novel_ids)][['CHROM', 'POS', 'REF', 'ALT']].copy()
        
        # Save novel variants
        novel_df.to_csv(f"{output_dir}/{dataset}_NOVEL_variants.csv", index=False)
        novel_variants[dataset] = novel_ids
        
        # Calculate database statistics
        stats = {
            'total_variants': total_variants,
            'in_gnomAD': len(in_gnomad),
            'in_1000G': len(in_1000g),
            'in_dbSNP': len(in_dbsnp),
            'novel_variants': len(novel_ids),
            'in_gnomAD_pct': (len(in_gnomad) / total_variants) * 100 if total_variants > 0 else 0,
            'in_1000G_pct': (len(in_1000g) / total_variants) * 100 if total_variants > 0 else 0,
            'in_dbSNP_pct': (len(in_dbsnp) / total_variants) * 100 if total_variants > 0 else 0,
            'novel_pct': (len(novel_ids) / total_variants) * 100 if total_variants > 0 else 0
        }
        
        database_stats[dataset] = stats
        all_data[dataset] = {
            'not_gnomad': not_gnomad_ids,
            'not_1000g': not_1000g_ids,
            'not_dbsnp': not_dbsnp_ids,
            'novel': novel_ids
        }
        
        print(f"  Novel variants (absent from all DBs): {len(novel_ids):,} ({stats['novel_pct']:.1f}%)")
        print()
    
    # Find common novel variants across all datasets
    find_common_novel_variants(novel_variants, output_dir)
    
    # Create all visualizations
    create_database_heatmap(database_stats, output_dir)
    create_venn_diagrams(all_data, output_dir)
    create_simple_upset_plot(novel_variants, output_dir)
    create_summary_report(database_stats, output_dir)
    
    print(f"✓ Novel variants analysis complete! Results saved in {output_dir}")
    print(f"  NOTE: Analysis restricted to autosomes (chr1-22) to avoid")
    print(f"  false novel calls from incomplete sex chromosome annotations.")


def find_common_novel_variants(novel_variants, output_dir):
    """Find and save common novel variants across all datasets."""
    
    # Find intersection of all novel variants
    all_datasets = list(novel_variants.keys())
    common_variants = novel_variants[all_datasets[0]]
    
    for dataset in all_datasets[1:]:
        common_variants = common_variants & novel_variants[dataset]
    
    print(f"Common novel variants across all pipelines: {len(common_variants):,}")
    
    # Convert variant IDs back to CHROM, POS, REF, ALT
    common_variants_list = []
    for variant_id in common_variants:
        parts = variant_id.split(':')
        if len(parts) == 4:
            common_variants_list.append({
                'CHROM': parts[0],
                'POS': int(parts[1]),
                'REF': parts[2],
                'ALT': parts[3]
            })
    
    # Save common variants
    if common_variants_list:
        common_df = pd.DataFrame(common_variants_list)
        common_df.to_csv(f"{output_dir}/COMMON_novel_variants_all_pipelines.csv", index=False)
        print(f"✓ Saved {len(common_variants_list):,} common novel variants")
    else:
        # Create empty file
        pd.DataFrame(columns=['CHROM', 'POS', 'REF', 'ALT']).to_csv(f"{output_dir}/COMMON_novel_variants_all_pipelines.csv", index=False)
        print("✓ No common novel variants found across all pipelines")


def create_database_heatmap(database_stats, output_dir):
    """Create heatmap showing percentage of variants in each database."""
    
    # Prepare data for heatmap
    datasets = ['GATK Illumina\nSpecific', 'PPR Illumina\nSpecific', 'GATK Aviti\nSpecific', 'PPR Aviti\nSpecific']
    databases = ['Present in\ngnomAD', 'Present in\n1000G', 'Present in\ndbSNP']
    
    # Create matrix
    data_matrix = []
    for dataset in ['GATK_Illumina', 'PPR_Illumina', 'GATK_Aviti', 'PPR_Aviti']:
        row = [
            database_stats[dataset]['in_gnomAD_pct'],
            database_stats[dataset]['in_1000G_pct'],
            database_stats[dataset]['in_dbSNP_pct']
        ]
        data_matrix.append(row)
    
    data_matrix = np.array(data_matrix)
    
    # Create heatmap
    fig, ax = plt.subplots(figsize=(16, 12))
    
    # Create heatmap with custom colormap
    im = ax.imshow(data_matrix, cmap='RdYlBu_r', aspect='auto', vmin=0, vmax=100)
    
    # Set ticks and labels
    ax.set_xticks(np.arange(len(databases)))
    ax.set_yticks(np.arange(len(datasets)))
    ax.set_xticklabels(databases, fontsize=22, fontweight='bold')
    ax.set_yticklabels(datasets, fontsize=22, fontweight='bold')
    
    # Rotate the tick labels with very small angle
    plt.setp(ax.get_xticklabels(), rotation=5, ha="right", rotation_mode="anchor")
    
    # Add text annotations with larger font
    for i in range(len(datasets)):
        for j in range(len(databases)):
            text = ax.text(j, i, f'{data_matrix[i, j]:.1f}%',
                          ha="center", va="center", color="black", fontsize=22, fontweight='bold')
    
    # Add colorbar
    cbar = ax.figure.colorbar(im, ax=ax, shrink=0.8)
    cbar.ax.set_ylabel('Percentage of Variants (%)', rotation=-90, va="bottom", fontsize=20, fontweight='bold')
    cbar.ax.tick_params(labelsize=18)
    
    # Set title with larger font
    ax.set_title('Pipeline-Specific Variants: Presence in Population Databases\n(Autosomes only)', 
                fontsize=26, fontweight='bold', pad=40)
    
    # Adjust layout
    fig.tight_layout()
    plt.savefig(f"{output_dir}/database_presence_heatmap.png", dpi=300, bbox_inches='tight')
    plt.close()
    print("✓ Created: database_presence_heatmap.png")


def create_venn_diagrams(all_data, output_dir):
    """Create Venn diagrams for each dataset showing overlap of variants absent from databases."""
    
    fig, axes = plt.subplots(2, 2, figsize=(28, 24))
    datasets = ['GATK_Illumina', 'PPR_Illumina', 'GATK_Aviti', 'PPR_Aviti']
    titles = ['GATK Illumina Specific', 'PPR Illumina Specific', 'GATK Aviti Specific', 'PPR Aviti Specific']
    
    for idx, (dataset, title) in enumerate(zip(datasets, titles)):
        row = idx // 2
        col = idx % 2
        ax = axes[row, col]
        
        # Get sets for this dataset
        not_gnomad = all_data[dataset]['not_gnomad']
        not_1000g = all_data[dataset]['not_1000g']
        not_dbsnp = all_data[dataset]['not_dbsnp']
        
        print(f"Venn for {dataset}:")
        print(f"  NOT in gnomAD: {len(not_gnomad):,}")
        print(f"  NOT in 1000G: {len(not_1000g):,}")
        print(f"  NOT in dbSNP: {len(not_dbsnp):,}")
        print(f"  Intersection (Novel): {len(not_gnomad & not_1000g & not_dbsnp):,}")
        
        # Create Venn diagram
        venn = venn3([not_gnomad, not_1000g, not_dbsnp], 
                     ('NOT in\ngnomAD', 'NOT in\n1000G', 'NOT in\ndbSNP'), 
                     ax=ax)
        
        # Customize colors
        colors = ['#FF6B6B', '#4ECDC4', '#45B7D1', '#96CEB4', '#FFEAA7', '#DDA0DD', '#FFB347']
        patch_ids = ['100', '010', '001', '110', '101', '011', '111']
        
        for patch_id, color in zip(patch_ids, colors):
            if venn.get_patch_by_id(patch_id):
                venn.get_patch_by_id(patch_id).set_color(color)
                venn.get_patch_by_id(patch_id).set_alpha(0.7)
        
        # Add circles with thicker lines
        venn3_circles([not_gnomad, not_1000g, not_dbsnp], ax=ax, linewidth=4)
        
        # Customize labels with much larger font
        for text in venn.set_labels:
            if text:
                text.set_fontsize(24)
                text.set_fontweight('bold')
        
        # Make subset labels MUCH LARGER and bolder
        for text in venn.subset_labels:
            if text:
                text.set_fontsize(24)
                text.set_fontweight('bold')
                text.set_bbox(dict(boxstyle="round,pad=0.4", facecolor='white', alpha=0.9, edgecolor='black', linewidth=1))
        
        # Set title
        ax.set_title(f'{title}\nVariants Absent from Databases (Autosomes)', 
                    fontsize=24, fontweight='bold', pad=30)
        
        # Add text box with novel variants count
        novel_count = len(all_data[dataset]['novel'])
        ax.text(0.02, 0.98, f'Completely Novel\n(absent from all DBs):\n{novel_count:,} variants', 
                transform=ax.transAxes, fontsize=20, fontweight='bold',
                verticalalignment='top', 
                bbox=dict(boxstyle='round,pad=0.5', facecolor='yellow', alpha=0.9, edgecolor='black', linewidth=2))
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/venn_diagrams_absent_variants.png", dpi=300, bbox_inches='tight')
    plt.close()
    print("✓ Created: venn_diagrams_absent_variants.png")


def create_simple_upset_plot(novel_variants, output_dir):
    """Create custom UpSet plot showing overlap of novel variants between pipelines."""
    
    # Define colors for each dataset
    colors = {
        'GATK_Illumina': '#0173B2',
        'PPR_Illumina': '#DE8F05', 
        'GATK_Aviti': '#029E73',
        'PPR_Aviti': '#CC78BC'
    }
    
    dataset_order = ['GATK_Illumina', 'PPR_Illumina', 'GATK_Aviti', 'PPR_Aviti']
    
    # Compute all possible intersections
    intersection_counts = {}
    
    # 1-way intersections (unique to each dataset)
    for dataset in dataset_order:
        unique_variants = novel_variants[dataset].copy()
        for other_dataset in dataset_order:
            if other_dataset != dataset:
                unique_variants -= novel_variants[other_dataset]
        if len(unique_variants) > 0:
            intersection_counts[(dataset,)] = len(unique_variants)
    
    # 2-way intersections
    for combo in combinations(dataset_order, 2):
        intersection = novel_variants[combo[0]] & novel_variants[combo[1]]
        excluded_datasets = [d for d in dataset_order if d not in combo]
        for dataset in excluded_datasets:
            intersection -= novel_variants[dataset]
        if len(intersection) > 0:
            intersection_counts[combo] = len(intersection)
    
    # 3-way intersections
    for combo in combinations(dataset_order, 3):
        intersection = novel_variants[combo[0]] & novel_variants[combo[1]] & novel_variants[combo[2]]
        excluded_datasets = [d for d in dataset_order if d not in combo]
        for dataset in excluded_datasets:
            intersection -= novel_variants[dataset]
        if len(intersection) > 0:
            intersection_counts[combo] = len(intersection)
    
    # 4-way intersection (all datasets)
    intersection = novel_variants[dataset_order[0]]
    for dataset in dataset_order[1:]:
        intersection &= novel_variants[dataset]
    if len(intersection) > 0:
        intersection_counts[tuple(dataset_order)] = len(intersection)
    
    # Sort by intersection size
    sorted_intersections = sorted(intersection_counts.items(), key=lambda x: x[1], reverse=True)
    
    print(f"Found {len(sorted_intersections)} non-empty intersections")
    for combo, count in sorted_intersections[:10]:
        print(f"  {' & '.join(combo)}: {count:,} variants")
    
    # Create figure
    fig = plt.figure(figsize=(20, 12))
    
    # Top panel: bar chart
    ax_bars = plt.subplot2grid((5, 1), (0, 0), rowspan=3)
    
    # Bottom panel: set membership matrix
    ax_matrix = plt.subplot2grid((5, 1), (3, 0), rowspan=2)
    
    n_intersections = len(sorted_intersections)
    x_positions = np.arange(n_intersections)
    
    # Draw bars with colors based on intersection size
    colors_bar = []
    for combo, count in sorted_intersections:
        if len(combo) == 4:
            colors_bar.append('#e74c3c')
        elif len(combo) == 3:
            colors_bar.append('#f39c12')
        elif len(combo) == 2:
            colors_bar.append('#3498db')
        else:
            colors_bar.append('#95a5a6')
    
    bars = ax_bars.bar(x_positions, [count for _, count in sorted_intersections],
                       color=colors_bar, edgecolor='black', linewidth=2)
    
    # Add count labels on bars
    max_count = max([c for _, c in sorted_intersections]) if sorted_intersections else 1
    for bar, (_, count) in zip(bars, sorted_intersections):
        height = bar.get_height()
        ax_bars.text(bar.get_x() + bar.get_width()/2., height + max_count * 0.01,
                    f'{int(count):,}', ha='center', va='bottom',
                    fontsize=18, fontweight='bold', color='black')
    
    ax_bars.set_ylabel('Number of Novel Variants', fontsize=18, fontweight='bold', color='black')
    ax_bars.set_title('UpSet Plot: Novel Variants Overlap Between Pipeline-Specific Datasets\n(Novel = absent from gnomAD, 1000G, and dbSNP; Autosomes only)', 
                      fontsize=20, fontweight='bold', color='black', pad=20)
    ax_bars.tick_params(axis='y', labelsize=16, colors='black')
    ax_bars.set_xticklabels([])
    ax_bars.grid(axis='y', alpha=0.3, linestyle='--', linewidth=1)
    ax_bars.spines['top'].set_visible(False)
    ax_bars.spines['right'].set_visible(False)
    
    # Set membership matrix
    y_positions = np.arange(len(dataset_order))
    
    for i, (combo, _) in enumerate(sorted_intersections):
        for j, dataset in enumerate(dataset_order):
            if dataset in combo:
                circle = Circle((i, j), 0.25, color=colors[dataset],
                              ec='black', linewidth=2, zorder=3)
                ax_matrix.add_patch(circle)
                
                if len(combo) > 1:
                    combo_indices = [dataset_order.index(d) for d in combo]
                    min_idx, max_idx = min(combo_indices), max(combo_indices)
                    if j >= min_idx and j <= max_idx:
                        ax_matrix.plot([i, i], [min_idx, max_idx],
                                     color='black', linewidth=4, zorder=2)
            else:
                circle = Circle((i, j), 0.12, color='white',
                              ec='gray', linewidth=1, zorder=1)
                ax_matrix.add_patch(circle)
    
    ax_matrix.set_xlim(-0.5, n_intersections - 0.5)
    ax_matrix.set_ylim(-0.5, len(dataset_order) - 0.5)
    ax_matrix.set_yticks(y_positions)
    ax_matrix.set_yticklabels([name.replace('_', ' ') for name in dataset_order],
                             fontsize=16, fontweight='bold', color='black')
    ax_matrix.set_xticks([])
    ax_matrix.set_xlabel('Intersection Groups', fontsize=18, fontweight='bold', color='black')
    ax_matrix.invert_yaxis()
    ax_matrix.grid(False)
    
    for spine in ax_matrix.spines.values():
        spine.set_visible(False)
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/upset_plot_novel_variants.png", dpi=300, bbox_inches='tight')
    plt.close()
    print("✓ Created: upset_plot_novel_variants.png")


def create_summary_report(database_stats, output_dir):
    """Create comprehensive summary report."""
    
    # Convert to DataFrame
    stats_df = pd.DataFrame(database_stats).T
    stats_df.to_csv(f'{output_dir}/novel_variants_summary.csv')
    
    # Create detailed report
    with open(f'{output_dir}/novel_variants_report.txt', 'w') as f:
        f.write("NOVEL VARIANTS ANALYSIS SUMMARY\n")
        f.write("=" * 60 + "\n")
        f.write("NOTE: Analysis restricted to autosomes (chr1-22) to avoid\n")
        f.write("false novel calls from incomplete sex chromosome annotations.\n")
        f.write("=" * 60 + "\n\n")
        
        f.write("Pipeline-Specific Variants Analysis:\n")
        f.write("-" * 40 + "\n")
        
        for dataset in stats_df.index:
            f.write(f"\n{dataset.replace('_', ' ')} Specific Variants:\n")
            f.write(f"  Total variants (autosomes): {stats_df.loc[dataset, 'total_variants']:,}\n")
            f.write(f"  Present in gnomAD: {stats_df.loc[dataset, 'in_gnomAD']:,} ({stats_df.loc[dataset, 'in_gnomAD_pct']:.1f}%)\n")
            f.write(f"  Present in 1000G: {stats_df.loc[dataset, 'in_1000G']:,} ({stats_df.loc[dataset, 'in_1000G_pct']:.1f}%)\n")
            f.write(f"  Present in dbSNP: {stats_df.loc[dataset, 'in_dbSNP']:,} ({stats_df.loc[dataset, 'in_dbSNP_pct']:.1f}%)\n")
            f.write(f"  NOVEL (absent from all DBs): {stats_df.loc[dataset, 'novel_variants']:,} ({stats_df.loc[dataset, 'novel_pct']:.1f}%)\n")
        
        f.write(f"\n\nSummary Statistics:\n")
        f.write("-" * 20 + "\n")
        f.write(f"Total novel variants across all pipelines: {stats_df['novel_variants'].sum():,}\n")
        f.write(f"Average novel variant rate: {stats_df['novel_pct'].mean():.1f}%\n")
        f.write(f"Pipeline with most novel variants: {stats_df['novel_variants'].idxmax().replace('_', ' ')}\n")
        f.write(f"Pipeline with highest novel rate: {stats_df['novel_pct'].idxmax().replace('_', ' ')}\n")
    
    print("✓ Created: novel_variants_summary.csv and novel_variants_report.txt")


if __name__ == "__main__":
    base_dir = "/mnt/genomics/pilot_PPR/uniquevariants/extract_csq_tsv/plots_frequencies"
    output_dir = "/mnt/genomics/pilot_PPR/uniquevariants/extract_csq_tsv/plots_populationdb/novelvariants_autosomesonly_py"
    
    analyze_novel_variants(base_dir, output_dir)
