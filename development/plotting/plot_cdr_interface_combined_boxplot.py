#!/usr/bin/env python3
"""
Generate combined boxplot for CDR RMSD and contact interface RMSD
Order: Alpha CDR1, CDR2, CDR3, Contact Interface, Beta CDR1, CDR2, CDR3
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path
import argparse
import logging

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class CDRInterfaceBoxplotGenerator:
    """Generate combined boxplot for CDR and contact interface RMSD"""

    def __init__(self, cdr_summary_file: str, interface_rmsd_dir: str, output_dir: str):
        self.cdr_summary_file = Path(cdr_summary_file)
        self.interface_rmsd_dir = Path(interface_rmsd_dir)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def load_cdr_data(self):
        """Load CDR RMSD data"""
        logger.info("Loading CDR RMSD data...")
        if not self.cdr_summary_file.exists():
            raise FileNotFoundError(f"CDR summary file not found: {self.cdr_summary_file}")

        cdr_df = pd.read_csv(self.cdr_summary_file)
        cdr_df = cdr_df[cdr_df['status'] == 'success'].copy()

        logger.info(f"Loaded {len(cdr_df)} CDR RMSD records")
        return cdr_df

    def load_interface_rmsd_data(self):
        """Load contact interface RMSD data from CSV files"""
        logger.info("Loading contact interface RMSD data...")

        csv_files = list(self.interface_rmsd_dir.glob('*_interface_rmsd.csv'))
        logger.info(f"Found {len(csv_files)} interface RMSD CSV files")

        interface_data = []

        for csv_file in csv_files:
            try:
                # Read CSV file
                df = pd.read_csv(csv_file)

                # Filter for post-20ns data (Time_ps >= 20000)
                post_equilibration = df[df['Time_ps'] >= 20000.0]

                if len(post_equilibration) > 0:
                    # Calculate mean RMSD for post-20ns
                    mean_rmsd = post_equilibration['RMSD_nm'].mean()

                    # Extract task name from filename
                    # e.g., "1ao7_run1_interface_rmsd.csv" -> "1ao7_run1"
                    task_name = csv_file.stem.replace('_interface_rmsd', '')

                    interface_data.append({
                        'task_name': task_name,
                        'mean_rmsd_nm': mean_rmsd,
                        'status': 'success'
                    })
                else:
                    logger.warning(f"No post-20ns data in {csv_file.name}")

            except Exception as e:
                logger.error(f"Error processing {csv_file.name}: {e}")
                continue

        interface_df = pd.DataFrame(interface_data)
        logger.info(f"Successfully loaded {len(interface_df)} interface RMSD records")
        logger.info(f"Interface RMSD mean: {interface_df['mean_rmsd_nm'].mean():.4f} nm, "
                   f"std: {interface_df['mean_rmsd_nm'].std():.4f} nm")

        return interface_df

    def prepare_combined_data(self, cdr_df, interface_df):
        """Prepare combined data for boxplot"""
        plot_data = []

        # Order: Alpha CDR1, CDR2, CDR3, Contact Interface, Beta CDR1, CDR2, CDR3
        region_order = [
            ('Alpha', 'CDR1'),
            ('Alpha', 'CDR2'),
            ('Alpha', 'CDR3'),
            ('Contact', 'Interface'),  # Contact Interface
            ('Beta', 'CDR1'),
            ('Beta', 'CDR2'),
            ('Beta', 'CDR3')
        ]

        for chain, region in region_order:
            if chain == 'Contact':
                # Contact Interface data
                subset = interface_df.copy()
                subset['region_label'] = 'Contact Interface'
                subset['rmsd_value'] = subset['mean_rmsd_nm']
            else:
                # CDR data - map display names to CSV column values
                chain_map = {'Alpha': 'TCR_alpha', 'Beta': 'TCR_beta'}
                csv_chain = chain_map[chain]

                subset = cdr_df[
                    (cdr_df['chain'] == csv_chain) &
                    (cdr_df['cdr_region'] == region)
                ].copy()
                subset['region_label'] = f'{chain} {region}'
                subset['rmsd_value'] = subset['mean_rmsd_nm']

            if len(subset) > 0:
                logger.info(f"{subset['region_label'].iloc[0]}: n={len(subset)}, "
                          f"mean={subset['rmsd_value'].mean():.3f} nm, "
                          f"std={subset['rmsd_value'].std():.3f} nm")
                plot_data.append(subset[['region_label', 'rmsd_value']])

        combined_df = pd.concat(plot_data, ignore_index=True)
        return combined_df, region_order

    def generate_boxplot(self, combined_df, region_order):
        """Generate combined boxplot"""
        logger.info("Generating combined boxplot...")

        # Set style
        sns.set_style("whitegrid")
        plt.rcParams['font.size'] = 12
        plt.rcParams['font.family'] = 'Arial'

        # Create figure
        fig, ax = plt.subplots(figsize=(14, 8))

        # Order labels
        ordered_labels = [f'{chain} {region}' if chain != 'Contact' else 'Contact Interface'
                         for chain, region in region_order]

        # Create boxplot
        bp = ax.boxplot(
            [combined_df[combined_df['region_label'] == label]['rmsd_value'].values
             for label in ordered_labels],
            labels=ordered_labels,
            patch_artist=True,
            showmeans=True,
            meanprops=dict(marker='D', markerfacecolor='red', markeredgecolor='red', markersize=6),
            medianprops=dict(color='blue', linewidth=2),
            boxprops=dict(facecolor='lightblue', edgecolor='black', linewidth=1.5),
            whiskerprops=dict(color='black', linewidth=1.5),
            capprops=dict(color='black', linewidth=1.5),
            flierprops=dict(marker='o', markerfacecolor='gray', markersize=4, alpha=0.5)
        )

        # Customize appearance
        ax.set_ylabel('RMSD (nm)', fontsize=14, fontweight='bold')
        ax.set_xlabel('Region', fontsize=14, fontweight='bold')
        ax.set_title('CDR and Contact Interface RMSD Distribution (post-20ns equilibration)',
                    fontsize=16, fontweight='bold', pad=20)

        # Rotate x-axis labels
        plt.xticks(rotation=45, ha='right')

        # Add grid
        ax.yaxis.grid(True, linestyle='--', alpha=0.7)
        ax.set_axisbelow(True)

        # Add legend
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor='lightblue', edgecolor='black', label='IQR'),
            plt.Line2D([0], [0], color='blue', linewidth=2, label='Median'),
            plt.Line2D([0], [0], marker='D', color='w', markerfacecolor='red',
                      markersize=8, label='Mean')
        ]
        ax.legend(handles=legend_elements, loc='upper right', fontsize=11)

        # Add sample size annotations
        for i, label in enumerate(ordered_labels, 1):
            n = len(combined_df[combined_df['region_label'] == label])
            y_pos = ax.get_ylim()[1] * 0.95
            ax.text(i, y_pos, f'n={n}', ha='center', va='top', fontsize=9, color='gray')

        # Adjust layout
        plt.tight_layout()

        # Save figure
        output_file = self.output_dir / 'cdr_interface_combined_boxplot.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        logger.info(f"Boxplot saved to: {output_file}")

        # Save high-res PDF
        output_pdf = self.output_dir / 'cdr_interface_combined_boxplot.pdf'
        plt.savefig(output_pdf, format='pdf', bbox_inches='tight')
        logger.info(f"PDF version saved to: {output_pdf}")

        plt.close()

        return output_file

    def generate_statistics_report(self, combined_df, region_order):
        """Generate statistics report"""
        logger.info("Generating statistics report...")

        ordered_labels = [f'{chain} {region}' if chain != 'Contact' else 'Contact Interface'
                         for chain, region in region_order]

        stats_lines = []
        stats_lines.append("=" * 80)
        stats_lines.append("CDR and Contact Interface RMSD Statistics (post-20ns equilibration)")
        stats_lines.append("=" * 80)
        stats_lines.append("")

        for label in ordered_labels:
            subset = combined_df[combined_df['region_label'] == label]['rmsd_value']
            if len(subset) > 0:
                stats_lines.append(f"{label}:")
                stats_lines.append(f"  N: {len(subset)}")
                stats_lines.append(f"  Mean: {subset.mean():.4f} nm")
                stats_lines.append(f"  Std: {subset.std():.4f} nm")
                stats_lines.append(f"  Median: {subset.median():.4f} nm")
                stats_lines.append(f"  Min: {subset.min():.4f} nm")
                stats_lines.append(f"  Max: {subset.max():.4f} nm")
                stats_lines.append(f"  Q1: {subset.quantile(0.25):.4f} nm")
                stats_lines.append(f"  Q3: {subset.quantile(0.75):.4f} nm")
                stats_lines.append("")

        stats_text = "\n".join(stats_lines)

        # Save to file
        output_file = self.output_dir / 'cdr_interface_combined_statistics.txt'
        with open(output_file, 'w') as f:
            f.write(stats_text)

        logger.info(f"Statistics report saved to: {output_file}")
        return output_file

    def run(self):
        """Run complete boxplot generation"""
        logger.info("Starting combined CDR-Interface boxplot generation...")

        # Load data
        cdr_df = self.load_cdr_data()
        interface_df = self.load_interface_rmsd_data()

        # Prepare combined data
        combined_df, region_order = self.prepare_combined_data(cdr_df, interface_df)

        # Generate boxplot
        plot_file = self.generate_boxplot(combined_df, region_order)

        # Generate statistics
        stats_file = self.generate_statistics_report(combined_df, region_order)

        logger.info("=" * 80)
        logger.info("Combined boxplot generation completed successfully!")
        logger.info(f"Plot: {plot_file}")
        logger.info(f"Stats: {stats_file}")
        logger.info("=" * 80)

        return plot_file, stats_file


def main():
    parser = argparse.ArgumentParser(
        description='Generate combined boxplot for CDR and contact interface RMSD'
    )
    parser.add_argument(
        '--cdr-summary',
        type=str,
        default='output/cdr_rmsd_exact_analysis/batch_cdr_rmsd_summary.csv',
        help='Path to CDR RMSD summary CSV file'
    )
    parser.add_argument(
        '--interface-rmsd-dir',
        type=str,
        default='output/rmsd_contactsurf',
        help='Directory containing interface RMSD CSV files'
    )
    parser.add_argument(
        '--output-dir',
        type=str,
        default='output/cdr_interface_combined_plots',
        help='Output directory for plots'
    )

    args = parser.parse_args()

    generator = CDRInterfaceBoxplotGenerator(
        cdr_summary_file=args.cdr_summary,
        interface_rmsd_dir=args.interface_rmsd_dir,
        output_dir=args.output_dir
    )

    try:
        generator.run()
    except Exception as e:
        logger.error(f"Error: {e}")
        import traceback
        traceback.print_exc()
        return 1

    return 0


if __name__ == '__main__':
    exit(main())
