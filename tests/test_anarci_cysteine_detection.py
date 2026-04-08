#!/usr/bin/env python3
"""
Test script for ANARCI-based conserved cysteine detection

This script tests the new integrated ANARCI strategy in DockingAngleAnalyzer.

Usage:
    python development/test_anarci_cysteine_detection.py
"""

import sys
from pathlib import Path

project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from immunex.analysis.angles import DockingAngleAnalyzer
import logging

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)

logger = logging.getLogger(__name__)


def test_single_structure():
    """Test ANARCI cysteine detection on a single structure."""

    # TODO: Update these paths to actual test data
    test_pdb = "path/to/test_structure.pdb"

    if not Path(test_pdb).exists():
        logger.error(f"Test file not found: {test_pdb}")
        logger.info("Please update test_pdb path in script")
        return

    logger.info("="*70)
    logger.info("Test 1: Single Structure Cysteine Detection")
    logger.info("="*70)

    # Initialize analyzer
    analyzer = DockingAngleAnalyzer(test_pdb)

    # Define selections (no need to specify cysteine residues!)
    selections = {
        'mhc': 'chainID A and resid 50:86 140:176 and name CA',
        'tcr_alpha': 'chainID D',
        'tcr_beta': 'chainID E',
        'tcr_v': 'chainID D E and resid 1:115 and name CA',
        'peptide': 'chainID C and name CA'
    }

    # Calculate docking angles (ANARCI will auto-detect cysteines)
    try:
        twist, tilt, swing = analyzer.calculate_docking_angles(
            mhc_selection=selections['mhc'],
            tcr_alpha_selection=selections['tcr_alpha'],
            tcr_beta_selection=selections['tcr_beta'],
            tcr_v_selection=selections['tcr_v'],
            peptide_selection=selections['peptide']
        )

        logger.info(f"\nDocking angles calculated successfully:")
        logger.info(f"  Twist: {twist:.2f} degrees")
        logger.info(f"  Tilt:  {tilt:.2f} degrees")
        logger.info(f"  Swing: {swing:.2f} degrees")

        logger.info("\nTest PASSED: ANARCI cysteine detection working!")

    except Exception as e:
        logger.error(f"Test FAILED: {e}")
        import traceback
        logger.error(traceback.format_exc())


def test_trajectory():
    """Test ANARCI cysteine detection on MD trajectory."""

    # TODO: Update these paths
    test_tpr = "path/to/test.tpr"
    test_xtc = "path/to/test.xtc"

    if not Path(test_tpr).exists() or not Path(test_xtc).exists():
        logger.error("Test files not found")
        logger.info("Please update test file paths in script")
        return

    logger.info("="*70)
    logger.info("Test 2: MD Trajectory Cysteine Detection")
    logger.info("="*70)

    # Initialize analyzer
    analyzer = DockingAngleAnalyzer(test_tpr, test_xtc)

    # Define selections
    selections = {
        'mhc': 'chainID A and resid 50:86 140:176 and name CA',
        'tcr_alpha': 'chainID D',
        'tcr_beta': 'chainID E',
        'tcr_v': 'chainID D E and resid 1:115 and name CA',
        'peptide': 'chainID C and name CA'
    }

    # Calculate for first 10 frames
    try:
        times, twist, tilt, swing = analyzer.calculate_docking_angles_trajectory(
            mhc_selection=selections['mhc'],
            tcr_alpha_selection=selections['tcr_alpha'],
            tcr_beta_selection=selections['tcr_beta'],
            tcr_v_selection=selections['tcr_v'],
            peptide_selection=selections['peptide'],
            stride=10,
            output_file="development/test_docking_angles.csv"
        )

        logger.info(f"\nTrajectory analysis completed:")
        logger.info(f"  Frames: {len(times)}")
        logger.info(f"  Twist: {twist.mean():.2f} ± {twist.std():.2f} degrees")
        logger.info(f"  Tilt:  {tilt.mean():.2f} ± {tilt.std():.2f} degrees")
        logger.info(f"  Swing: {swing.mean():.2f} ± {swing.std():.2f} degrees")

        logger.info("\nTest PASSED: Trajectory analysis working!")

    except Exception as e:
        logger.error(f"Test FAILED: {e}")
        import traceback
        logger.error(traceback.format_exc())


def main():
    """Run all tests."""
    logger.info("Testing ANARCI-based conserved cysteine detection")
    logger.info("This replaces hardcoded residue numbers with automatic IMGT detection\n")

    # Test 1: Single structure
    test_single_structure()

    print("\n")

    # Test 2: Trajectory
    test_trajectory()

    logger.info("\n" + "="*70)
    logger.info("All tests completed!")
    logger.info("="*70)


if __name__ == "__main__":
    main()
