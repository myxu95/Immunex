#!/usr/bin/env python3
"""
Unified RMSD Calculation Script

This script provides a command-line interface for calculating RMSD with flexible
component selection for alignment and calculation.

Supported components:
- HLA: MHC alpha and beta chains
- pHLA: pMHC complex (MHC + peptide)
- peptide: Antigenic peptide
- TCR: T-cell receptor (both chains)
- TCR_alpha: TCR alpha chain
- TCR_beta: TCR beta chain
- CDR3_alpha: CDR3 region of TCR alpha (requires sequence)
- CDR3_beta: CDR3 region of TCR beta (requires sequence)

Examples:
    # Calculate TCR RMSD aligned to pHLA
    python calculate_rmsd.py md.tpr md.xtc --align pHLA --calc TCR

    # Calculate CDR3 RMSD aligned to TCR
    python calculate_rmsd.py md.tpr md.xtc --align TCR --calc CDR3_beta --cdr3-beta CASSLGQAYEQYF

    # Batch calculation with multiple alignments
    python calculate_rmsd.py md.tpr md.xtc --batch-config rmsd_config.json

Author: Immunex Development Team
"""

import sys
import argparse
import logging
import json
from pathlib import Path
from typing import Dict, List, Optional

sys.path.insert(0, str(Path(__file__).parent.parent))

from immunex.analysis.trajectory import RMSDInterface

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Unified RMSD Calculation with Flexible Component Selection",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage: TCR RMSD aligned to pHLA
  python %(prog)s md.tpr md.xtc --align pHLA --calc TCR

  # CDR3 RMSD with sequence
  python %(prog)s md.tpr md.xtc --align TCR --calc CDR3_beta --cdr3-beta CASSLGQAYEQYF

  # Multiple calculations
  python %(prog)s md.tpr md.xtc --align pHLA --calc TCR --calc peptide

  # Batch mode with config file
  python %(prog)s md.tpr md.xtc --batch-config config.json

Available components:
  HLA, pHLA, peptide, TCR, TCR_alpha, TCR_beta, CDR3_alpha, CDR3_beta
        """
    )

    parser.add_argument(
        "topology",
        help="Topology file (.tpr, .gro, .pdb)"
    )

    parser.add_argument(
        "trajectory",
        help="Trajectory file (.xtc, .trr)"
    )

    parser.add_argument(
        "--align",
        default="pHLA",
        help="Component for alignment (default: pHLA)"
    )

    parser.add_argument(
        "--calc",
        nargs="+",
        default=["TCR"],
        help="Component(s) for RMSD calculation (default: TCR)"
    )

    parser.add_argument(
        "--cdr3-alpha",
        metavar="SEQUENCE",
        help="CDR3 alpha sequence (required for CDR3_alpha calculation)"
    )

    parser.add_argument(
        "--cdr3-beta",
        metavar="SEQUENCE",
        help="CDR3 beta sequence (required for CDR3_beta calculation)"
    )

    parser.add_argument(
        "-o", "--output-dir",
        default="rmsd_results",
        help="Output directory for results (default: rmsd_results)"
    )

    parser.add_argument(
        "--no-standardize",
        action="store_true",
        help="Skip automatic chain standardization"
    )

    parser.add_argument(
        "--batch-config",
        metavar="FILE",
        help="JSON config file for batch calculations"
    )

    parser.add_argument(
        "--output-csv",
        metavar="FILE",
        help="CSV file for batch results summary"
    )

    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable verbose output"
    )

    return parser.parse_args()


def load_batch_config(config_file: str) -> List[Dict]:
    """
    Load batch calculation configuration from JSON file.

    Args:
        config_file: Path to JSON config file

    Returns:
        List of calculation configurations

    Example JSON format:
        [
            {
                "align": "pHLA",
                "calc": "TCR",
                "output_file": "tcr_rmsd.xvg"
            },
            {
                "align": "TCR",
                "calc": "CDR3_beta",
                "cdr3_sequences": {"beta": "CASSLGQAYEQYF"}
            }
        ]
    """
    with open(config_file, 'r') as f:
        config = json.load(f)

    if not isinstance(config, list):
        raise ValueError("Config file must contain a JSON array")

    return config


def main():
    """Main function."""
    args = parse_arguments()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    logger.info("="*80)
    logger.info("Unified RMSD Calculation")
    logger.info("="*80)
    logger.info(f"Topology: {args.topology}")
    logger.info(f"Trajectory: {args.trajectory}")
    logger.info(f"Output directory: {args.output_dir}")
    logger.info(f"Chain standardization: {'Disabled' if args.no_standardize else 'Enabled'}")
    logger.info("="*80 + "\n")

    try:
        rmsd_interface = RMSDInterface(
            topology=args.topology,
            trajectory=args.trajectory,
            auto_standardize=not args.no_standardize,
            output_dir=args.output_dir
        )

        cdr3_sequences = {}
        if args.cdr3_alpha:
            cdr3_sequences['alpha'] = args.cdr3_alpha
        if args.cdr3_beta:
            cdr3_sequences['beta'] = args.cdr3_beta

        if args.batch_config:
            logger.info("Running in batch mode\n")
            calculations = load_batch_config(args.batch_config)

            df = rmsd_interface.batch_calculate(
                calculations=calculations,
                output_csv=args.output_csv
            )

            print("\nBatch Results Summary:")
            print(df.to_string(index=False))

        else:
            calculations = []
            for calc_component in args.calc:
                calculations.append({
                    'align': args.align,
                    'calc': calc_component,
                    'cdr3_sequences': cdr3_sequences if calc_component.startswith('CDR3') else None
                })

            if len(calculations) > 1:
                df = rmsd_interface.batch_calculate(
                    calculations=calculations,
                    output_csv=args.output_csv
                )

                print("\nResults Summary:")
                print(df[['align_component', 'calc_component', 'rmsd_mean_nm', 'rmsd_std_nm']].to_string(index=False))

            else:
                result = rmsd_interface.calculate(
                    align=args.align,
                    calc=args.calc[0],
                    cdr3_sequences=cdr3_sequences if args.calc[0].startswith('CDR3') else None
                )

                print(f"\nRMSD Statistics ({args.align} -> {args.calc[0]}):")
                print(f"  Mean: {result['mean']:.4f} nm")
                print(f"  Std:  {result['std']:.4f} nm")
                print(f"  Min:  {result['min']:.4f} nm")
                print(f"  Max:  {result['max']:.4f} nm")
                print(f"  Frames: {result['n_frames']}")
                print(f"\nOutput: {result['output_file']}")

        logger.info("\n" + "="*80)
        logger.info("RMSD Calculation Complete")
        logger.info("="*80)

    except Exception as e:
        logger.error(f"\nFailed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
