#!/usr/bin/env python3
"""
Hydrogen Bond Analysis Module for pHLA-TCR Complexes

Analyzes hydrogen bonds at key interfaces:
- Peptide-HLA interface
- TCR-pHLA interface

Author: Immunex Development Team
Date: 2025-12-07
"""

import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA
from typing import Dict, Optional
import logging
from pathlib import Path

logger = logging.getLogger(__name__)


class pHLATCRHydrogenBondAnalyzer:
    """Hydrogen bond analysis for pHLA-TCR interfaces."""

    def __init__(self, universe: mda.Universe):
        """
        Initialize hydrogen bond analyzer.

        Args:
            universe: MDAnalysis Universe object
        """
        self.u = universe

    def analyze_interface_hbonds(self,
                                 donors_selection: str,
                                 acceptors_selection: str,
                                 distance_cutoff: float = 3.5,
                                 angle_cutoff: float = 150.0) -> Dict:
        """
        Analyze hydrogen bonds at a specific interface.

        Args:
            donors_selection: MDAnalysis selection string for donors
            acceptors_selection: MDAnalysis selection string for acceptors
            distance_cutoff: H...A distance cutoff (Angstroms)
            angle_cutoff: D-H...A angle cutoff (degrees)

        Returns:
            Dictionary with hydrogen bond analysis results
        """
        logger.info(f"Analyzing H-bonds: donors='{donors_selection}', acceptors='{acceptors_selection}'")

        try:
            # Create HBA instance
            hbonds = HBA(
                universe=self.u,
                donors_sel=donors_selection,
                acceptors_sel=acceptors_selection,
                d_h_cutoff=1.2,  # D-H distance
                d_a_cutoff=distance_cutoff,  # D-A distance
                d_h_a_angle_cutoff=angle_cutoff
            )

            # Run analysis
            hbonds.run()

            # Get results
            n_frames = len(self.u.trajectory)

            # Count unique H-bonds per frame
            if hasattr(hbonds, 'results') and hasattr(hbonds.results, 'hbonds'):
                hbonds_data = hbonds.results.hbonds

                # Count H-bonds per frame
                frame_counts = []
                times = []

                for frame_idx in range(n_frames):
                    frame_hbonds = hbonds_data[hbonds_data[:, 0] == frame_idx]
                    frame_counts.append(len(frame_hbonds))
                    times.append(self.u.trajectory[frame_idx].time / 1000.0)

                results = {
                    'times': times,
                    'counts': frame_counts,
                    'mean_count': np.mean(frame_counts),
                    'std_count': np.std(frame_counts),
                    'total_unique_hbonds': len(hbonds_data)
                }

            else:
                logger.warning("H-bond results not available")
                results = {
                    'times': [],
                    'counts': [],
                    'mean_count': 0,
                    'std_count': 0,
                    'total_unique_hbonds': 0
                }

            return results

        except Exception as e:
            logger.error(f"H-bond analysis failed: {e}")
            return {
                'times': [],
                'counts': [],
                'mean_count': 0,
                'std_count': 0,
                'total_unique_hbonds': 0,
                'error': str(e)
            }

    def analyze_peptide_hla_hbonds(self) -> Dict:
        """
        Analyze hydrogen bonds at Peptide-HLA interface.

        Returns:
            H-bond analysis results
        """
        logger.info("Analyzing Peptide-HLA hydrogen bonds...")

        return self.analyze_interface_hbonds(
            donors_selection="chainID C",
            acceptors_selection="chainID A or chainID B"
        )

    def analyze_tcr_phla_hbonds(self) -> Dict:
        """
        Analyze hydrogen bonds at TCR-pHLA interface.

        Returns:
            H-bond analysis results
        """
        logger.info("Analyzing TCR-pHLA hydrogen bonds...")

        return self.analyze_interface_hbonds(
            donors_selection="chainID D or chainID E",
            acceptors_selection="chainID A or chainID B or chainID C"
        )

    def save_hbond_timeseries(self, results: Dict, output_file: Path):
        """
        Save hydrogen bond time series to file.

        Args:
            results: H-bond analysis results
            output_file: Output CSV file path
        """
        if not results['times']:
            logger.warning("No H-bond data to save")
            return

        df = pd.DataFrame({
            'time_ns': results['times'],
            'hbond_count': results['counts']
        })

        df.to_csv(output_file, index=False)
        logger.info(f"H-bond timeseries saved to {output_file}")


def main():
    """Test hydrogen bond analyzer."""
    import sys

    if len(sys.argv) < 3:
        print("Usage: python hydrogen_bond_analyzer.py <topology> <trajectory>")
        sys.exit(1)

    topology = sys.argv[1]
    trajectory = sys.argv[2]

    # Load trajectory
    u = mda.Universe(topology, trajectory)

    # Create analyzer
    hb_analyzer = pHLATCRHydrogenBondAnalyzer(u)

    # Analyze Peptide-HLA H-bonds
    results_pep_hla = hb_analyzer.analyze_peptide_hla_hbonds()
    print("\n=== Peptide-HLA H-bonds ===")
    print(f"Mean count: {results_pep_hla['mean_count']:.2f} ± {results_pep_hla['std_count']:.2f}")
    print(f"Total unique: {results_pep_hla['total_unique_hbonds']}")

    # Analyze TCR-pHLA H-bonds
    results_tcr = hb_analyzer.analyze_tcr_phla_hbonds()
    print("\n=== TCR-pHLA H-bonds ===")
    print(f"Mean count: {results_tcr['mean_count']:.2f} ± {results_tcr['std_count']:.2f}")
    print(f"Total unique: {results_tcr['total_unique_hbonds']}")


if __name__ == "__main__":
    main()
