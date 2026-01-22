#!/usr/bin/env python3
"""Cross-system comparative analysis for allostery networks

Compares clustering patterns, correlation strengths, and allosteric networks
across different TCR-pMHC complexes.
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import logging

sys.path.insert(0, str(Path(__file__).parent.parent))

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def load_clustering_summary(analysis_dir):
    """Load clustering summary for all systems"""
    analysis_path = Path(analysis_dir)
    data = []

    for task_dir in sorted(analysis_path.iterdir()):
        if not task_dir.is_dir():
            continue

        task_name = task_dir.name
        summary_file = task_dir / "clustering" / "cluster_summary.csv"

        if not summary_file.exists():
            logger.warning(f"Missing clustering data: {task_name}")
            continue

        try:
            # Load cluster summary
            cluster_df = pd.read_csv(summary_file)

            # Load contact labels to get total contacts
            labels_file = task_dir / "contact_labels.csv"
            if labels_file.exists():
                labels_df = pd.read_csv(labels_file)
                n_contacts = len(labels_df)
            else:
                n_contacts = 0

            # Compute statistics
            data.append({
                'System': task_name,
                'N_Clusters': len(cluster_df),
                'N_Contacts': n_contacts,
                'Max_Cluster_Size': cluster_df['Size'].max() if len(cluster_df) > 0 else 0,
                'Mean_Cluster_Size': cluster_df['Size'].mean() if len(cluster_df) > 0 else 0,
                'Total_Clustered': cluster_df['Size'].sum() if len(cluster_df) > 0 else 0,
                'Mean_Intra_Corr': cluster_df['Mean_Intra_Correlation'].mean() if len(cluster_df) > 0 else 0,
                'Std_Intra_Corr': cluster_df['Mean_Intra_Correlation'].std() if len(cluster_df) > 0 else 0
            })

        except Exception as e:
            logger.error(f"Error processing {task_name}: {e}")
            continue

    return pd.DataFrame(data)


def plot_cluster_distribution(df, output_dir):
    """Plot cluster count distribution across systems"""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(2, 2, figsize=(16, 12))

    # 1. Cluster count distribution
    axes[0, 0].hist(df['N_Clusters'], bins=30, edgecolor='black', alpha=0.7)
    axes[0, 0].set_xlabel('Number of Clusters', fontweight='bold')
    axes[0, 0].set_ylabel('Number of Systems', fontweight='bold')
    axes[0, 0].set_title('Distribution of Cluster Counts', fontweight='bold')
    axes[0, 0].axvline(df['N_Clusters'].mean(), color='red', linestyle='--',
                      label=f'Mean: {df["N_Clusters"].mean():.1f}')
    axes[0, 0].legend()

    # 2. Cluster size distribution
    axes[0, 1].scatter(df['N_Clusters'], df['Max_Cluster_Size'], alpha=0.6)
    axes[0, 1].set_xlabel('Number of Clusters', fontweight='bold')
    axes[0, 1].set_ylabel('Max Cluster Size', fontweight='bold')
    axes[0, 1].set_title('Cluster Count vs Max Cluster Size', fontweight='bold')

    # 3. Clustering ratio
    df['Clustering_Ratio'] = df['Total_Clustered'] / df['N_Contacts']
    axes[1, 0].hist(df['Clustering_Ratio'] * 100, bins=30, edgecolor='black', alpha=0.7)
    axes[1, 0].set_xlabel('Clustering Ratio (%)', fontweight='bold')
    axes[1, 0].set_ylabel('Number of Systems', fontweight='bold')
    axes[1, 0].set_title('Distribution of Clustering Ratios', fontweight='bold')
    axes[1, 0].axvline(df['Clustering_Ratio'].mean() * 100, color='red', linestyle='--',
                      label=f'Mean: {df["Clustering_Ratio"].mean()*100:.1f}%')
    axes[1, 0].legend()

    # 4. Mean intra-cluster correlation
    axes[1, 1].hist(df['Mean_Intra_Corr'], bins=30, edgecolor='black', alpha=0.7)
    axes[1, 1].set_xlabel('Mean Intra-Cluster Correlation', fontweight='bold')
    axes[1, 1].set_ylabel('Number of Systems', fontweight='bold')
    axes[1, 1].set_title('Distribution of Intra-Cluster Correlations', fontweight='bold')
    axes[1, 1].axvline(df['Mean_Intra_Corr'].mean(), color='red', linestyle='--',
                      label=f'Mean: {df["Mean_Intra_Corr"].mean():.3f}')
    axes[1, 1].legend()

    plt.tight_layout()
    plt.savefig(output_dir / 'cluster_distribution.png', dpi=300, bbox_inches='tight')
    plt.close()
    logger.info(f"Saved: {output_dir / 'cluster_distribution.png'}")


def plot_top_systems(df, output_dir, n_top=20):
    """Plot top systems by various metrics"""
    output_dir = Path(output_dir)

    fig, axes = plt.subplots(2, 2, figsize=(18, 14))

    # 1. Top systems by cluster count
    top_clusters = df.nlargest(n_top, 'N_Clusters')
    axes[0, 0].barh(range(len(top_clusters)), top_clusters['N_Clusters'])
    axes[0, 0].set_yticks(range(len(top_clusters)))
    axes[0, 0].set_yticklabels(top_clusters['System'], fontsize=9)
    axes[0, 0].set_xlabel('Number of Clusters', fontweight='bold')
    axes[0, 0].set_title(f'Top {n_top} Systems by Cluster Count', fontweight='bold')
    axes[0, 0].invert_yaxis()

    # 2. Top systems by max cluster size
    top_size = df.nlargest(n_top, 'Max_Cluster_Size')
    axes[0, 1].barh(range(len(top_size)), top_size['Max_Cluster_Size'], color='orange')
    axes[0, 1].set_yticks(range(len(top_size)))
    axes[0, 1].set_yticklabels(top_size['System'], fontsize=9)
    axes[0, 1].set_xlabel('Max Cluster Size (residues)', fontweight='bold')
    axes[0, 1].set_title(f'Top {n_top} Systems by Max Cluster Size', fontweight='bold')
    axes[0, 1].invert_yaxis()

    # 3. Top systems by clustering ratio
    df['Clustering_Ratio'] = df['Total_Clustered'] / df['N_Contacts']
    top_ratio = df.nlargest(n_top, 'Clustering_Ratio')
    axes[1, 0].barh(range(len(top_ratio)), top_ratio['Clustering_Ratio'] * 100, color='green')
    axes[1, 0].set_yticks(range(len(top_ratio)))
    axes[1, 0].set_yticklabels(top_ratio['System'], fontsize=9)
    axes[1, 0].set_xlabel('Clustering Ratio (%)', fontweight='bold')
    axes[1, 0].set_title(f'Top {n_top} Systems by Clustering Ratio', fontweight='bold')
    axes[1, 0].invert_yaxis()

    # 4. Top systems by mean correlation
    top_corr = df.nlargest(n_top, 'Mean_Intra_Corr')
    axes[1, 1].barh(range(len(top_corr)), top_corr['Mean_Intra_Corr'], color='red')
    axes[1, 1].set_yticks(range(len(top_corr)))
    axes[1, 1].set_yticklabels(top_corr['System'], fontsize=9)
    axes[1, 1].set_xlabel('Mean Intra-Cluster Correlation', fontweight='bold')
    axes[1, 1].set_title(f'Top {n_top} Systems by Correlation Strength', fontweight='bold')
    axes[1, 1].invert_yaxis()

    plt.tight_layout()
    plt.savefig(output_dir / 'top_systems.png', dpi=300, bbox_inches='tight')
    plt.close()
    logger.info(f"Saved: {output_dir / 'top_systems.png'}")


def generate_summary_report(df, output_dir):
    """Generate text summary report"""
    output_dir = Path(output_dir)
    report_file = output_dir / 'comparative_analysis_report.txt'

    df['Clustering_Ratio'] = df['Total_Clustered'] / df['N_Contacts']

    with open(report_file, 'w') as f:
        f.write("="*80 + "\n")
        f.write("Cross-System Comparative Analysis Report\n")
        f.write("="*80 + "\n\n")

        f.write("Dataset Overview:\n")
        f.write("-"*80 + "\n")
        f.write(f"Total systems analyzed: {len(df)}\n")
        f.write(f"Average contacts per system: {df['N_Contacts'].mean():.1f} +/- {df['N_Contacts'].std():.1f}\n")
        f.write(f"Contact range: {df['N_Contacts'].min():.0f} - {df['N_Contacts'].max():.0f}\n\n")

        f.write("Clustering Statistics:\n")
        f.write("-"*80 + "\n")
        f.write(f"Average clusters per system: {df['N_Clusters'].mean():.1f} +/- {df['N_Clusters'].std():.1f}\n")
        f.write(f"Cluster range: {df['N_Clusters'].min():.0f} - {df['N_Clusters'].max():.0f}\n")
        f.write(f"Average max cluster size: {df['Max_Cluster_Size'].mean():.1f} +/- {df['Max_Cluster_Size'].std():.1f}\n")
        f.write(f"Average clustering ratio: {df['Clustering_Ratio'].mean()*100:.2f}% +/- {df['Clustering_Ratio'].std()*100:.2f}%\n")
        f.write(f"Clustering ratio range: {df['Clustering_Ratio'].min()*100:.2f}% - {df['Clustering_Ratio'].max()*100:.2f}%\n\n")

        f.write("Correlation Statistics:\n")
        f.write("-"*80 + "\n")
        f.write(f"Average intra-cluster correlation: {df['Mean_Intra_Corr'].mean():.3f} +/- {df['Mean_Intra_Corr'].std():.3f}\n")
        f.write(f"Correlation range: {df['Mean_Intra_Corr'].min():.3f} - {df['Mean_Intra_Corr'].max():.3f}\n\n")

        f.write("Top 10 Systems by Cluster Count:\n")
        f.write("-"*80 + "\n")
        top_clusters = df.nlargest(10, 'N_Clusters')
        for idx, row in top_clusters.iterrows():
            f.write(f"{row['System']:20s}: {row['N_Clusters']:.0f} clusters "
                   f"(max size: {row['Max_Cluster_Size']:.0f}, ratio: {row['Clustering_Ratio']*100:.1f}%)\n")

        f.write("\nTop 10 Systems by Clustering Ratio:\n")
        f.write("-"*80 + "\n")
        top_ratio = df.nlargest(10, 'Clustering_Ratio')
        for idx, row in top_ratio.iterrows():
            f.write(f"{row['System']:20s}: {row['Clustering_Ratio']*100:.2f}% "
                   f"({row['Total_Clustered']:.0f}/{row['N_Contacts']:.0f} contacts)\n")

        f.write("\nTop 10 Systems by Max Cluster Size:\n")
        f.write("-"*80 + "\n")
        top_size = df.nlargest(10, 'Max_Cluster_Size')
        for idx, row in top_size.iterrows():
            f.write(f"{row['System']:20s}: {row['Max_Cluster_Size']:.0f} residues "
                   f"({row['N_Clusters']:.0f} total clusters)\n")

        f.write("\n" + "="*80 + "\n")
        f.write("Key Findings:\n")
        f.write("="*80 + "\n")

        high_clustering = df[df['Clustering_Ratio'] > df['Clustering_Ratio'].quantile(0.75)]
        low_clustering = df[df['Clustering_Ratio'] < df['Clustering_Ratio'].quantile(0.25)]

        f.write(f"\n1. High clustering systems (top 25%, n={len(high_clustering)}):\n")
        f.write(f"   - Average clustering ratio: {high_clustering['Clustering_Ratio'].mean()*100:.2f}%\n")
        f.write(f"   - Average cluster count: {high_clustering['N_Clusters'].mean():.1f}\n")
        f.write(f"   - Average max cluster size: {high_clustering['Max_Cluster_Size'].mean():.1f}\n")

        f.write(f"\n2. Low clustering systems (bottom 25%, n={len(low_clustering)}):\n")
        f.write(f"   - Average clustering ratio: {low_clustering['Clustering_Ratio'].mean()*100:.2f}%\n")
        f.write(f"   - Average cluster count: {low_clustering['N_Clusters'].mean():.1f}\n")
        f.write(f"   - Average max cluster size: {low_clustering['Max_Cluster_Size'].mean():.1f}\n")

        f.write(f"\n3. Correlation strength observations:\n")
        f.write(f"   - Systems with strong correlations (>0.5): {len(df[df['Mean_Intra_Corr'] > 0.5])}\n")
        f.write(f"   - Systems with moderate correlations (0.3-0.5): {len(df[(df['Mean_Intra_Corr'] >= 0.3) & (df['Mean_Intra_Corr'] <= 0.5)])}\n")
        f.write(f"   - Systems with weak correlations (<0.3): {len(df[df['Mean_Intra_Corr'] < 0.3])}\n")

    logger.info(f"Saved: {report_file}")


def main():
    analysis_dir = Path("output/allostery_analysis")
    output_dir = Path("output/comparative_analysis")
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("="*80)
    logger.info("Cross-System Comparative Analysis")
    logger.info("="*80)

    # Load data
    logger.info("\nLoading clustering summaries...")
    df = load_clustering_summary(analysis_dir)
    logger.info(f"Loaded data for {len(df)} systems")

    if len(df) == 0:
        logger.error("No data loaded!")
        return 1

    # Save summary table
    summary_file = output_dir / "system_statistics.csv"
    df.to_csv(summary_file, index=False)
    logger.info(f"Saved: {summary_file}")

    # Generate visualizations
    logger.info("\nGenerating visualizations...")
    plot_cluster_distribution(df, output_dir)
    plot_top_systems(df, output_dir)

    # Generate report
    logger.info("\nGenerating summary report...")
    generate_summary_report(df, output_dir)

    logger.info("\n" + "="*80)
    logger.info("Comparative analysis completed!")
    logger.info(f"Output directory: {output_dir}")
    logger.info("="*80)

    return 0


if __name__ == "__main__":
    sys.exit(main())
