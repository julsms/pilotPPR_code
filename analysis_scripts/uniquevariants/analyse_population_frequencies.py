#!/usr/bin/env python3
# analyze_population_frequencies.py

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Set style for better plots
plt.style.use('default')
sns.set_palette("husl")


def process_frequency_data(input_dir, output_dir):
    """Process frequency data and create comprehensive analysis."""

    # Create output directory
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # Define input files
    files = {
        'GATK_Illumina': f'{input_dir}/uniqueGATKillumina_frequencies.tsv',
        'PPR_Illumina': f'{input_dir}/uniquePPRillumina_frequencies.tsv',
        'GATK_Aviti': f'{input_dir}/uniqueGATKaviti_frequencies.tsv',
        'PPR_Aviti': f'{input_dir}/uniquePPRaviti_frequencies.tsv'
    }

    # Colors for each dataset
    colors = {
        'GATK_Illumina': '#0173B2',
        'PPR_Illumina': '#DE8F05',
        'GATK_Aviti': '#029E73',
        'PPR_Aviti': '#CC78BC'
    }

    all_stats = {}

    # Process each file
    for dataset_name, file_path in files.items():
        print(f"Processing {dataset_name}...")

        # Read data
        df = pd.read_csv(file_path, sep='\t')

        # Clean and convert frequency columns
        df['gnomADg_AF_clean'] = pd.to_numeric(
            df['gnomADg_AF'].replace(['', '.', 'NA'], np.nan),
            errors='coerce'
        )

        df['1000G_AF_clean'] = pd.to_numeric(
            df['1000G_AF'].replace(['', '.', 'NA'], np.nan),
            errors='coerce'
        )

        # Categorize variants
        df['in_gnomAD'] = (
            (~df['gnomADg_AF_clean'].isna()) &
            (df['gnomADg_AF_clean'] > 0)
        )

        df['in_1000G'] = (
            (~df['1000G_AF_clean'].isna()) &
            (df['1000G_AF_clean'] > 0)
        )

        # Save categorized variants
        base_cols = ['CHROM', 'POS', 'REF', 'ALT', 'SYMBOL', 'Gene']
        freq_cols = [
            'gnomADg_AF',
            '1000G_AF',
            'gnomADg_AF_clean',
            '1000G_AF_clean'
        ]

        # Variants in gnomAD
        gnomad_variants = df[df['in_gnomAD']][base_cols + freq_cols]
        gnomad_variants.to_csv(
            f'{output_dir}/{dataset_name}_in_gnomAD.csv',
            index=False
        )

        # Variants NOT in gnomAD
        not_gnomad_variants = df[~df['in_gnomAD']][base_cols + freq_cols]
        not_gnomad_variants.to_csv(
            f'{output_dir}/{dataset_name}_NOT_in_gnomAD.csv',
            index=False
        )

        # Variants in 1000G
        kg_variants = df[df['in_1000G']][base_cols + freq_cols]
        kg_variants.to_csv(
            f'{output_dir}/{dataset_name}_in_1000G.csv',
            index=False
        )

        # Variants NOT in 1000G
        not_kg_variants = df[~df['in_1000G']][base_cols + freq_cols]
        not_kg_variants.to_csv(
            f'{output_dir}/{dataset_name}_NOT_in_1000G.csv',
            index=False
        )

        # Collect statistics
        stats = {
            'total_variants': len(df),
            'in_gnomAD': df['in_gnomAD'].sum(),
            'not_in_gnomAD': (~df['in_gnomAD']).sum(),
            'in_1000G': df['in_1000G'].sum(),
            'not_in_1000G': (~df['in_1000G']).sum(),
            'in_both_db': (
                df['in_gnomAD'] & df['in_1000G']
            ).sum(),
            'in_gnomAD_only': (
                df['in_gnomAD'] & ~df['in_1000G']
            ).sum(),
            'in_1000G_only': (
                ~df['in_gnomAD'] & df['in_1000G']
            ).sum(),
            'in_neither_db': (
                (~df['in_gnomAD']) & (~df['in_1000G'])
            ).sum()
        }

        all_stats[dataset_name] = stats

        print(f"  Total variants: {stats['total_variants']:,}")
        print(
            f"  In gnomAD: {stats['in_gnomAD']:,} "
            f"({stats['in_gnomAD']/stats['total_variants']*100:.1f}%)"
        )
        print(
            f"  In 1000G: {stats['in_1000G']:,} "
            f"({stats['in_1000G']/stats['total_variants']*100:.1f}%)"
        )
        print(
            f"  In both databases: {stats['in_both_db']:,} "
            f"({stats['in_both_db']/stats['total_variants']*100:.1f}%)"
        )
        print(
            f"  Not in either database: {stats['in_neither_db']:,} "
            f"({stats['in_neither_db']/stats['total_variants']*100:.1f}%)"
        )
        print()

    # Create plots and reports
    create_frequency_plot(files, colors, output_dir)
    create_summary_stats(all_stats, output_dir)
    create_comparison_plots(all_stats, colors, output_dir)

    print(f"✓ Analysis complete! Results saved in {output_dir}")


def create_frequency_plot(files, colors, output_dir):
    """Create allele frequency distribution plots with larger text."""

    fig, axes = plt.subplots(2, 2, figsize=(22, 18))

    for idx, (dataset_name, file_path) in enumerate(files.items()):
        row = idx // 2
        col = idx % 2
        ax = axes[row, col]

        # Read data
        df = pd.read_csv(file_path, sep='\t')

        # Clean frequency data
        gnomad_af = pd.to_numeric(
            df['gnomADg_AF'].replace(['', '.', 'NA'], np.nan),
            errors='coerce'
        )

        kg_af = pd.to_numeric(
            df['1000G_AF'].replace(['', '.', 'NA'], np.nan),
            errors='coerce'
        )

        # Remove NaN and zero values for log transformation
        gnomad_af_clean = gnomad_af[
            (gnomad_af > 0) & (~gnomad_af.isna())
        ]

        kg_af_clean = kg_af[
            (kg_af > 0) & (~kg_af.isna())
        ]

        # Log-transform
        gnomad_log = (
            np.log10(gnomad_af_clean)
            if len(gnomad_af_clean) > 0
            else np.array([])
        )

        kg_log = (
            np.log10(kg_af_clean)
            if len(kg_af_clean) > 0
            else np.array([])
        )

        # Plot histograms
        if len(gnomad_log) > 0:
            ax.hist(
                gnomad_log,
                bins=40,
                alpha=0.75,
                label=f'gnomAD Genomes (n={len(gnomad_log):,})',
                color='#1976D2',
                edgecolor='black',
                linewidth=1.5
            )

        if len(kg_log) > 0:
            ax.hist(
                kg_log,
                bins=40,
                alpha=0.75,
                label=f'1000 Genomes (n={len(kg_log):,})',
                color='#F57C00',
                edgecolor='black',
                linewidth=1.5
            )

        # Styling
        ax.set_title(
            dataset_name.replace("_", " "),
            fontsize=26,
            fontweight='bold',
            pad=20,
            color='black'
        )

        ax.set_xlabel(
            'log₁₀(Allele Frequency)',
            fontsize=22,
            fontweight='bold',
            color='black'
        )

        ax.set_ylabel(
            'Number of Variants',
            fontsize=22,
            fontweight='bold',
            color='black'
        )

        legend = ax.legend(
            loc='upper left',
            fontsize=18,
            frameon=True
        )

        legend.get_frame().set_edgecolor('black')
        legend.get_frame().set_linewidth(2)

        ax.tick_params(
            axis='both',
            labelsize=18,
            colors='black',
            width=2,
            length=6
        )

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_linewidth(2)
        ax.spines['bottom'].set_linewidth(2)

        ax.yaxis.grid(True, alpha=0.3, linewidth=1)
        ax.set_axisbelow(True)

        ax.set_xticks([-6, -5, -4, -3, -2, -1, 0])
        ax.set_xticklabels(
            ['10⁻⁶', '10⁻⁵', '10⁻⁴', '10⁻³', '10⁻²', '10⁻¹', '10⁰'],
            fontsize=18,
            fontweight='bold'
        )

    plt.suptitle(
        'Allele Frequency Distributions in Population Databases',
        fontsize=30,
        fontweight='bold',
        y=0.99,
        color='black'
    )

    plt.tight_layout(rect=[0, 0, 1, 0.96])

    plt.savefig(
        f"{output_dir}/allele_frequency_distributions.png",
        dpi=300,
        bbox_inches='tight'
    )

    plt.close()
    print("✓ Created: allele_frequency_distributions.png")


def create_summary_stats(all_stats, output_dir):
    """Create summary statistics table."""

    stats_df = pd.DataFrame(all_stats).T

    percentage_cols = [
        'in_gnomAD',
        'not_in_gnomAD',
        'in_1000G',
        'not_in_1000G',
        'in_both_db',
        'in_gnomAD_only',
        'in_1000G_only',
        'in_neither_db'
    ]

    for col in percentage_cols:
        stats_df[f'{col}_pct'] = (
            stats_df[col] / stats_df['total_variants'] * 100
        ).round(1)

    # Save detailed stats
    stats_df.to_csv(
        f'{output_dir}/summary_statistics.csv'
    )

    # Create readable summary
    with open(
        f'{output_dir}/summary_report.txt',
        'w'
    ) as f:

        f.write(
            "POPULATION DATABASE ANALYSIS SUMMARY\n"
        )
        f.write("=" * 50 + "\n\n")

        for dataset in stats_df.index:
            f.write(
                f"{dataset.replace('_', ' ')}:\n"
            )
            f.write(
                f"  Total variants: "
                f"{stats_df.loc[dataset, 'total_variants']:,}\n"
            )
            f.write(
                f"  In gnomAD: "
                f"{stats_df.loc[dataset, 'in_gnomAD']:,} "
                f"({stats_df.loc[dataset, 'in_gnomAD_pct']:.1f}%)\n"
            )
            f.write(
                f"  In 1000G: "
                f"{stats_df.loc[dataset, 'in_1000G']:,} "
                f"({stats_df.loc[dataset, 'in_1000G_pct']:.1f}%)\n"
            )
            f.write(
                f"  In both databases: "
                f"{stats_df.loc[dataset, 'in_both_db']:,} "
                f"({stats_df.loc[dataset, 'in_both_db_pct']:.1f}%)\n"
            )
            f.write(
                f"  In gnomAD only: "
                f"{stats_df.loc[dataset, 'in_gnomAD_only']:,} "
                f"({stats_df.loc[dataset, 'in_gnomAD_only_pct']:.1f}%)\n"
            )
            f.write(
                f"  In 1000G only: "
                f"{stats_df.loc[dataset, 'in_1000G_only']:,} "
                f"({stats_df.loc[dataset, 'in_1000G_only_pct']:.1f}%)\n"
            )
            f.write(
                f"  Not in either database: "
                f"{stats_df.loc[dataset, 'in_neither_db']:,} "
                f"({stats_df.loc[dataset, 'in_neither_db_pct']:.1f}%)\n"
            )
            f.write("\n")

    print(
        "✓ Created: summary_statistics.csv "
        "and summary_report.txt"
    )


def create_comparison_plots(
    all_stats,
    colors,
    output_dir
):
    """Create comparison bar plots with larger text."""

    datasets = list(all_stats.keys())
    datasets_clean = [
        d.replace('_', ' ')
        for d in datasets
    ]

    # Calculate percentages
    gnomad_pct = [
        all_stats[d]['in_gnomAD']
        / all_stats[d]['total_variants'] * 100
        for d in datasets
    ]

    kg_pct = [
        all_stats[d]['in_1000G']
        / all_stats[d]['total_variants'] * 100
        for d in datasets
    ]

    both_pct = [
        all_stats[d]['in_both_db']
        / all_stats[d]['total_variants'] * 100
        for d in datasets
    ]

    neither_pct = [
        all_stats[d]['in_neither_db']
        / all_stats[d]['total_variants'] * 100
        for d in datasets
    ]

    fig, (ax1, ax2) = plt.subplots(
        1,
        2,
        figsize=(22, 10)
    )

    # Plot 1
    x = np.arange(len(datasets))
    width = 0.25

    bars1 = ax1.bar(
        x - width,
        gnomad_pct,
        width,
        label='In gnomAD',
        color='#1976D2',
        alpha=0.8,
        edgecolor='black',
        linewidth=1
    )

    bars2 = ax1.bar(
        x,
        kg_pct,
        width,
        label='In 1000G',
        color='#F57C00',
        alpha=0.8,
        edgecolor='black',
        linewidth=1
    )

    bars3 = ax1.bar(
        x + width,
        both_pct,
        width,
        label='In Both',
        color='#4CAF50',
        alpha=0.8,
        edgecolor='black',
        linewidth=1
    )

    ax1.set_xlabel(
        'Dataset',
        fontsize=22,
        fontweight='bold'
    )

    ax1.set_ylabel(
        'Percentage of Variants (%)',
        fontsize=22,
        fontweight='bold'
    )

    ax1.set_title(
        'Variants Present in Population Databases',
        fontsize=24,
        fontweight='bold',
        pad=20
    )

    ax1.set_xticks(x)
    ax1.set_xticklabels(
        datasets_clean,
        fontsize=18,
        fontweight='bold'
    )

    ax1.legend(
        fontsize=18,
        frameon=True,
        edgecolor='black'
    )

    ax1.tick_params(
        axis='y',
        labelsize=18
    )

    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['left'].set_linewidth(2)
    ax1.spines['bottom'].set_linewidth(2)

    for bars in [bars1, bars2, bars3]:
        for bar in bars:
            height = bar.get_height()

            ax1.annotate(
                f'{height:.1f}%',
                xy=(
                    bar.get_x() + bar.get_width() / 2,
                    height
                ),
                xytext=(0, 4),
                textcoords="offset points",
                ha='center',
                va='bottom',
                fontsize=14,
                fontweight='bold'
            )

    # Plot 2
    bars4 = ax2.bar(
        datasets_clean,
        neither_pct,
        color=[colors[d] for d in datasets],
        alpha=0.8,
        edgecolor='black',
        linewidth=1.5
    )

    ax2.set_xlabel(
        'Dataset',
        fontsize=22,
        fontweight='bold'
    )

    ax2.set_ylabel(
        'Percentage of Variants (%)',
        fontsize=22,
        fontweight='bold'
    )

    ax2.set_title(
        'Variants Not in Any Population Database',
        fontsize=24,
        fontweight='bold',
        pad=20
    )

    ax2.tick_params(
        axis='both',
        labelsize=18
    )

    ax2.set_xticklabels(
        datasets_clean,
        fontsize=18,
        fontweight='bold'
    )

    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.spines['left'].set_linewidth(2)
    ax2.spines['bottom'].set_linewidth(2)

    for bar in bars4:
        height = bar.get_height()

        ax2.annotate(
            f'{height:.1f}%',
            xy=(
                bar.get_x() + bar.get_width() / 2,
                height
            ),
            xytext=(0, 4),
            textcoords="offset points",
            ha='center',
            va='bottom',
            fontsize=16,
            fontweight='bold'
        )

    plt.tight_layout()

    plt.savefig(
        f"{output_dir}/database_comparison.png",
        dpi=300,
        bbox_inches='tight'
    )

    plt.close()
    print("✓ Created: database_comparison.png")


if __name__ == "__main__":
    input_dir = "/mnt/genomics/pilot_PPR/uniquevariants/extract_csq_tsv"
    output_dir = (
        "/mnt/genomics/pilot_PPR/uniquevariants/"
        "extract_csq_tsv/plots_frequencies"
    )

    process_frequency_data(
        input_dir,
        output_dir
    )
