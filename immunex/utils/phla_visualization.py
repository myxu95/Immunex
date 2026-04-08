#!/usr/bin/env python3
"""
Visualization Module for pHLA-TCR Analysis Results

Generate comprehensive plots for trajectory analysis results.

Author: Immunex Development Team
Date: 2025-12-07
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from typing import Optional, List
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 10


class pHLATCRVisualizer:
    """Visualization toolkit for pHLA-TCR analysis results."""

    def __init__(self, results_dir: str):
        """
        Initialize visualizer.

        Args:
            results_dir: Directory containing analysis results
        """
        self.results_dir = Path(results_dir)
        self.figures_dir = self.results_dir / 'figures'
        self.figures_dir.mkdir(exist_ok=True)

    def plot_rmsd_comparison(self, output_file: Optional[str] = None):
        """
        Plot RMSD comparison for all chains.

        Args:
            output_file: Output file path (default: figures/rmsd_comparison.png)
        """
        logger.info("Plotting RMSD comparison...")

        # Load data
        chain_rmsd_file = self.results_dir / 'stability' / 'chain_rmsd.csv'

        if not chain_rmsd_file.exists():
            logger.error(f"Chain RMSD file not found: {chain_rmsd_file}")
            return

        df = pd.read_csv(chain_rmsd_file)

        # Create figure
        fig, ax = plt.subplots(figsize=(10, 6))

        # Plot each chain
        chains = {
            'Peptide': {'color': '#E74C3C', 'label': 'Peptide'},
            'HLA_alpha': {'color': '#3498DB', 'label': 'HLA-α'},
            'HLA_beta': {'color': '#2ECC71', 'label': 'HLA-β'},
            'TCR_alpha': {'color': '#F39C12', 'label': 'TCR-α'},
            'TCR_beta': {'color': '#9B59B6', 'label': 'TCR-β'}
        }

        for chain, props in chains.items():
            col_name = f'{chain}_rmsd'
            if col_name in df.columns:
                ax.plot(df['time_ns'], df[col_name],
                       label=props['label'], color=props['color'],
                       linewidth=2, alpha=0.8)

        ax.set_xlabel('Time (ns)', fontsize=12, fontweight='bold')
        ax.set_ylabel('RMSD (Å)', fontsize=12, fontweight='bold')
        ax.set_title('Per-Chain RMSD Evolution', fontsize=14, fontweight='bold')
        ax.legend(loc='best', frameon=True, shadow=True)
        ax.grid(True, alpha=0.3)

        plt.tight_layout()

        # Save figure
        if output_file is None:
            output_file = self.figures_dir / 'rmsd_comparison.png'

        plt.savefig(output_file, bbox_inches='tight')
        logger.info(f"RMSD comparison saved to {output_file}")
        plt.close()

    def plot_interface_distances(self, output_file: Optional[str] = None):
        """
        Plot interface distance evolution.

        Args:
            output_file: Output file path
        """
        logger.info("Plotting interface distances...")

        # Load data
        interface_file = self.results_dir / 'interface' / 'interface_distances.csv'

        if not interface_file.exists():
            logger.error(f"Interface distances file not found: {interface_file}")
            return

        df = pd.read_csv(interface_file)

        # Create figure with subplots
        fig, axes = plt.subplots(3, 1, figsize=(10, 10))

        # Plot 1: Peptide-HLA distance
        if 'Peptide_HLA_COM' in df.columns:
            axes[0].plot(df['time_ns'], df['Peptide_HLA_COM'],
                        color='#E74C3C', linewidth=2)
            axes[0].axhline(y=df['Peptide_HLA_COM'].mean(),
                           color='k', linestyle='--', alpha=0.5,
                           label=f"Mean: {df['Peptide_HLA_COM'].mean():.2f} Å")
            axes[0].set_ylabel('Distance (Å)', fontweight='bold')
            axes[0].set_title('Peptide-HLA COM Distance', fontweight='bold')
            axes[0].legend()
            axes[0].grid(True, alpha=0.3)

        # Plot 2: TCR-Peptide distance
        if 'TCR_Peptide_min' in df.columns:
            axes[1].plot(df['time_ns'], df['TCR_Peptide_min'],
                        color='#3498DB', linewidth=2)
            axes[1].axhline(y=df['TCR_Peptide_min'].mean(),
                           color='k', linestyle='--', alpha=0.5,
                           label=f"Mean: {df['TCR_Peptide_min'].mean():.2f} Å")
            axes[1].set_ylabel('Distance (Å)', fontweight='bold')
            axes[1].set_title('TCR-Peptide Minimum Distance', fontweight='bold')
            axes[1].legend()
            axes[1].grid(True, alpha=0.3)

        # Plot 3: TCR-HLA distance
        if 'TCR_HLA_COM' in df.columns:
            axes[2].plot(df['time_ns'], df['TCR_HLA_COM'],
                        color='#2ECC71', linewidth=2)
            axes[2].axhline(y=df['TCR_HLA_COM'].mean(),
                           color='k', linestyle='--', alpha=0.5,
                           label=f"Mean: {df['TCR_HLA_COM'].mean():.2f} Å")
            axes[2].set_xlabel('Time (ns)', fontweight='bold')
            axes[2].set_ylabel('Distance (Å)', fontweight='bold')
            axes[2].set_title('TCR-HLA COM Distance', fontweight='bold')
            axes[2].legend()
            axes[2].grid(True, alpha=0.3)

        plt.tight_layout()

        # Save figure
        if output_file is None:
            output_file = self.figures_dir / 'interface_distances.png'

        plt.savefig(output_file, bbox_inches='tight')
        logger.info(f"Interface distances plot saved to {output_file}")
        plt.close()

    def plot_contact_evolution(self, output_file: Optional[str] = None):
        """
        Plot contact count evolution for key interfaces.

        Args:
            output_file: Output file path
        """
        logger.info("Plotting contact evolution...")

        contacts_dir = self.results_dir / 'contacts'

        if not contacts_dir.exists():
            logger.error(f"Contacts directory not found: {contacts_dir}")
            return

        # Create figure
        fig, axes = plt.subplots(3, 1, figsize=(10, 10))

        # Define interfaces
        interfaces = [
            ('Peptide_HLA_contacts.csv', 'Peptide-HLA Contacts', '#E74C3C'),
            ('TCR_Peptide_contacts.csv', 'TCR-Peptide Contacts', '#3498DB'),
            ('TCR_HLA_contacts.csv', 'TCR-HLA Contacts', '#2ECC71')
        ]

        for idx, (filename, title, color) in enumerate(interfaces):
            file_path = contacts_dir / filename

            if not file_path.exists():
                logger.warning(f"Contact file not found: {file_path}")
                continue

            df = pd.read_csv(file_path)

            axes[idx].plot(df['time_ns'], df['contact_count'],
                          color=color, linewidth=2)
            axes[idx].axhline(y=df['contact_count'].mean(),
                             color='k', linestyle='--', alpha=0.5,
                             label=f"Mean: {df['contact_count'].mean():.1f}")

            if idx == 2:
                axes[idx].set_xlabel('Time (ns)', fontweight='bold')

            axes[idx].set_ylabel('Contact Count', fontweight='bold')
            axes[idx].set_title(title, fontweight='bold')
            axes[idx].legend()
            axes[idx].grid(True, alpha=0.3)

        plt.tight_layout()

        # Save figure
        if output_file is None:
            output_file = self.figures_dir / 'contact_evolution.png'

        plt.savefig(output_file, bbox_inches='tight')
        logger.info(f"Contact evolution plot saved to {output_file}")
        plt.close()

    def plot_radius_gyration(self, output_file: Optional[str] = None):
        """
        Plot radius of gyration evolution.

        Args:
            output_file: Output file path
        """
        logger.info("Plotting radius of gyration...")

        # Load data
        rg_file = self.results_dir / 'stability' / 'radius_gyration.csv'

        if not rg_file.exists():
            logger.error(f"Rg file not found: {rg_file}")
            return

        df = pd.read_csv(rg_file)

        # Create figure
        fig, axes = plt.subplots(2, 1, figsize=(10, 8))

        # Plot 1: Complex Rg
        if 'Rg_complex' in df.columns:
            axes[0].plot(df['time_ns'], df['Rg_complex'],
                        color='#34495E', linewidth=2)
            axes[0].axhline(y=df['Rg_complex'].mean(),
                           color='k', linestyle='--', alpha=0.5,
                           label=f"Mean: {df['Rg_complex'].mean():.2f} Å")
            axes[0].set_ylabel('Rg (Å)', fontweight='bold')
            axes[0].set_title('Complex Radius of Gyration', fontweight='bold')
            axes[0].legend()
            axes[0].grid(True, alpha=0.3)

        # Plot 2: Per-chain Rg
        chains = {
            'Rg_Peptide': {'color': '#E74C3C', 'label': 'Peptide'},
            'Rg_HLA_alpha': {'color': '#3498DB', 'label': 'HLA-α'},
            'Rg_HLA_beta': {'color': '#2ECC71', 'label': 'HLA-β'},
            'Rg_TCR_alpha': {'color': '#F39C12', 'label': 'TCR-α'},
            'Rg_TCR_beta': {'color': '#9B59B6', 'label': 'TCR-β'}
        }

        for col, props in chains.items():
            if col in df.columns:
                axes[1].plot(df['time_ns'], df[col],
                           label=props['label'], color=props['color'],
                           linewidth=2, alpha=0.8)

        axes[1].set_xlabel('Time (ns)', fontweight='bold')
        axes[1].set_ylabel('Rg (Å)', fontweight='bold')
        axes[1].set_title('Per-Chain Radius of Gyration', fontweight='bold')
        axes[1].legend(loc='best', frameon=True, shadow=True)
        axes[1].grid(True, alpha=0.3)

        plt.tight_layout()

        # Save figure
        if output_file is None:
            output_file = self.figures_dir / 'radius_gyration.png'

        plt.savefig(output_file, bbox_inches='tight')
        logger.info(f"Radius of gyration plot saved to {output_file}")
        plt.close()

    def plot_rmsf_comparison(self, output_file: Optional[str] = None):
        """
        Plot RMSF comparison for all chains.

        Generates a dual-panel figure:
        - Top: Per-residue RMSF profiles for all chains (overlaid line plots)
        - Bottom: Mean RMSF per chain (bar plot with error bars)

        Biological Interpretation:
        - Low RMSF peaks: Rigid regions (structural core, binding sites)
        - High RMSF peaks: Flexible regions (loops, CDRs, termini)
        - Peptide RMSF: Indicator of binding stability
        - TCR RMSF: Reflects entropic clamp effect

        Args:
            output_file: Output file path (default: figures/rmsf_comparison.png)
        """
        logger.info("Plotting RMSF comparison...")

        # Load chain RMSF data
        chain_rmsf_file = self.results_dir / 'dynamics' / 'chain_rmsf.csv'

        if not chain_rmsf_file.exists():
            logger.error(f"Chain RMSF file not found: {chain_rmsf_file}")
            return

        df = pd.read_csv(chain_rmsf_file)

        # Create figure with two subplots
        fig, axes = plt.subplots(2, 1, figsize=(12, 8))

        # Color scheme (matching existing RMSD visualization)
        chain_colors = {
            'Peptide': '#E74C3C',
            'HLA_alpha': '#3498DB',
            'HLA_beta': '#2ECC71',
            'TCR_alpha': '#F39C12',
            'TCR_beta': '#9B59B6'
        }

        chain_labels = {
            'Peptide': 'Peptide',
            'HLA_alpha': 'HLA-α',
            'HLA_beta': 'HLA-β',
            'TCR_alpha': 'TCR-α',
            'TCR_beta': 'TCR-β'
        }

        # Plot 1: Per-residue RMSF profiles (overlaid)
        ax1 = axes[0]

        for chain in ['Peptide', 'HLA_alpha', 'HLA_beta', 'TCR_alpha', 'TCR_beta']:
            chain_data = df[df['chain_id'] == chain]
            if len(chain_data) > 0:
                # Use sequential residue index for x-axis
                x_values = np.arange(len(chain_data))
                ax1.plot(x_values, chain_data['rmsf_angstrom'],
                        label=chain_labels[chain],
                        color=chain_colors[chain],
                        linewidth=2, alpha=0.8)

        ax1.set_xlabel('Residue Index (Sequential)', fontsize=12, fontweight='bold')
        ax1.set_ylabel('RMSF (Å)', fontsize=12, fontweight='bold')
        ax1.set_title('Per-Chain RMSF Profile', fontsize=14, fontweight='bold')
        ax1.legend(loc='best', frameon=True, shadow=True)
        ax1.grid(True, alpha=0.3)

        # Add 1.0 Å threshold line as reference
        ax1.axhline(y=1.0, color='gray', linestyle='--', alpha=0.5, linewidth=1.5,
                   label='1.0 Å (rigid threshold)')

        # Plot 2: Mean RMSF per chain (bar plot with error bars)
        ax2 = axes[1]

        chain_means = []
        chain_stds = []
        chain_names_ordered = []

        for chain in ['Peptide', 'HLA_alpha', 'HLA_beta', 'TCR_alpha', 'TCR_beta']:
            chain_data = df[df['chain_id'] == chain]
            if len(chain_data) > 0:
                chain_means.append(chain_data['rmsf_angstrom'].mean())
                chain_stds.append(chain_data['rmsf_angstrom'].std())
                chain_names_ordered.append(chain_labels[chain])

        x_pos = np.arange(len(chain_names_ordered))
        colors = [chain_colors[k] for k in ['Peptide', 'HLA_alpha', 'HLA_beta', 'TCR_alpha', 'TCR_beta']
                  if k in [chain for chain, data in zip(['Peptide', 'HLA_alpha', 'HLA_beta', 'TCR_alpha', 'TCR_beta'],
                                                         [df[df['chain_id'] == c] for c in ['Peptide', 'HLA_alpha', 'HLA_beta', 'TCR_alpha', 'TCR_beta']])
                           if len(data) > 0]]

        bars = ax2.bar(x_pos, chain_means, yerr=chain_stds,
                      color=colors,
                      alpha=0.7, capsize=5, error_kw={'linewidth': 2})

        ax2.set_xlabel('Chain', fontsize=12, fontweight='bold')
        ax2.set_ylabel('Mean RMSF (Å)', fontsize=12, fontweight='bold')
        ax2.set_title('Mean RMSF by Chain', fontsize=14, fontweight='bold')
        ax2.set_xticks(x_pos)
        ax2.set_xticklabels(chain_names_ordered)
        ax2.grid(True, alpha=0.3, axis='y')

        # Add value labels on top of bars
        for i, (mean, std) in enumerate(zip(chain_means, chain_stds)):
            ax2.text(i, mean + std + 0.1, f'{mean:.2f}',
                    ha='center', va='bottom', fontweight='bold', fontsize=9)

        plt.tight_layout()

        # Save figure
        if output_file is None:
            output_file = self.figures_dir / 'rmsf_comparison.png'

        plt.savefig(output_file, bbox_inches='tight')
        logger.info(f"RMSF comparison plot saved to {output_file}")
        plt.close()

    def generate_all_plots(self):
        """Generate all standard plots."""
        logger.info("Generating all plots...")

        self.plot_rmsd_comparison()
        self.plot_rmsf_comparison()
        self.plot_interface_distances()
        self.plot_contact_evolution()
        self.plot_radius_gyration()

        logger.info(f"All plots saved to {self.figures_dir}")


def main():
    """Test visualization."""
    import sys

    if len(sys.argv) < 2:
        print("Usage: python visualization.py <results_dir>")
        sys.exit(1)

    results_dir = sys.argv[1]

    # Create visualizer
    viz = pHLATCRVisualizer(results_dir)

    # Generate all plots
    viz.generate_all_plots()

    print(f"\nPlots saved to {viz.figures_dir}")


if __name__ == "__main__":
    main()
