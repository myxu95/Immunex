#!/usr/bin/env python3
"""
IMN - Immunex Command Line Interface

Main entry point for all Immunex CLI commands.

Author: Immunex Development Team
Date: 2026-03-15
"""

import argparse
import sys
from importlib import import_module


SUBCOMMANDS = {
    'preprocess': 'Process MD trajectories (PBC correction)',
    'quality': 'Quality assessment and validation',
    'rmsd': 'RMSD calculation and analysis',
    'rmsf': 'RMSF calculation and analysis',
    'contact': 'Contact analysis',
    'identity': 'Biological identity annotation',
    'bsa': 'Buried surface area analysis',
    'nma': 'Normal mode / ENM analysis',
    'inter_cluster': 'Interface-aware state clustering',
    'angle': 'Docking angle analysis',
    'pdb': 'PDB file processing',
    'batch': 'Batch processing workflows',
    'report': 'Generate analysis reports',
}


def print_usage():
    """Print usage information"""
    print("""
IMN - Immunex MD Analysis Toolkit

Usage:
    imn <subcommand> [options]

Available subcommands:
    preprocess    Process MD trajectories (PBC correction)
    quality       Quality assessment and validation
    rmsd          RMSD calculation and analysis
    rmsf          RMSF calculation and analysis
    contact       Contact analysis
    identity      Biological identity annotation
    bsa           Buried surface area analysis
    nma           Normal mode / ENM analysis
    inter_cluster Interface-aware state clustering
    angle         Docking angle analysis
    pdb           PDB file processing
    batch         Batch processing workflows
    report        Generate analysis reports

Examples:
    imn preprocess -f md.xtc -t md.tpr -o processed.xtc
    imn quality -f processed.xtc -t md.tpr
    imn rmsd -f processed.xtc -t md.tpr -s backbone
    imn identity --structure md_processed_converted.pdb -o ./output/identity_demo
    imn bsa -f md_processed.xtc -s md.tpr --structure md_processed_converted.pdb -o ./bsa
    imn nma --structure md_processed_converted.pdb -o ./output/nma_demo
    imn inter_cluster -f md_processed.xtc -s md.tpr --structure md_processed_converted.pdb -o ./output/interface_cluster
    imn report interaction --base-dir output/interaction_case_1OGA --system-id 1OGA_sd_run2
    imn pdb download 1ao7 2ckb
    imn batch preprocess /data/md/ --workers 4

For help on specific subcommand:
    imn preprocess --help
    imn quality --help
    """)


def main():
    """Main entry point for IMN CLI"""

    # If no arguments, print usage
    if len(sys.argv) < 2:
        print_usage()
        return 0

    # Check if subcommand is valid
    subcommand = sys.argv[1]

    # Handle help
    if subcommand in ['-h', '--help', 'help']:
        print_usage()
        return 0

    # Handle version
    if subcommand in ['-v', '--version', 'version']:
        from immunex import __version__
        print(f"IMN (Immunex) version {__version__}")
        return 0

    # Validate subcommand
    if subcommand not in SUBCOMMANDS:
        print(f"Error: Unknown subcommand '{subcommand}'")
        print(f"\nAvailable subcommands: {', '.join(SUBCOMMANDS.keys())}")
        print(f"\nRun 'imn --help' for more information.")
        return 1

    # Import and run subcommand
    try:
        module = import_module(f'immunex.cli.commands.{subcommand}')
        # Pass remaining arguments to subcommand
        return module.main(sys.argv[2:])
    except ImportError as e:
        print(f"Error loading subcommand '{subcommand}': {e}")
        import traceback
        traceback.print_exc()
        return 1
    except Exception as e:
        print(f"Error running subcommand '{subcommand}': {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == '__main__':
    sys.exit(main())
