#!/usr/bin/env python3
"""
Concatenate Multi-Model PDB Files

Command-line tool to merge multi-model PDB files into single PDB files
with continuous atom numbering. No alignment or coordinate changes.

Usage:
    # Single file
    python scripts/concatenate_multimodel_pdbs.py input.pdb -o output.pdb

    # Batch processing
    python scripts/concatenate_multimodel_pdbs.py input/*.pdb -o output/ --batch

    # With report
    python scripts/concatenate_multimodel_pdbs.py input.pdb -o output.pdb --report result.json

Author: Immunex Development Team
Date: 2026-01
"""

import argparse
import sys
import json
import logging
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from immunex.utils.multimodel_concatenator import MultiModelConcatenator


def setup_logging(verbose: bool = False):
    """Configure logging."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(levelname)s: %(message)s'
    )


def main():
    parser = argparse.ArgumentParser(
        description='Concatenate multi-model PDB into single PDB with continuous atom numbering',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Single file processing
  python scripts/concatenate_multimodel_pdbs.py input/6VMA.pdb -o output/6VMA_concatenated.pdb

  # Batch processing
  python scripts/concatenate_multimodel_pdbs.py input/*.pdb -o output/ --batch

  # Generate JSON report
  python scripts/concatenate_multimodel_pdbs.py input/6VMA.pdb -o output/6VMA.pdb --report result.json

  # Without preserving remarks
  python scripts/concatenate_multimodel_pdbs.py input/6VMA.pdb -o output/6VMA.pdb --no-remarks

Notes:
  - All MODEL records are removed from output
  - Atom numbering is continuous (1, 2, 3, ...)
  - Original coordinates are preserved unchanged
  - Residue numbering is preserved (may have duplicates across models)
        """
    )

    parser.add_argument(
        'input',
        nargs='+',
        help='Input multi-model PDB file(s)'
    )
    parser.add_argument(
        '-o', '--output',
        required=True,
        help='Output PDB file or directory (for batch mode)'
    )
    parser.add_argument(
        '--batch',
        action='store_true',
        help='Batch mode: process multiple files to output directory'
    )
    parser.add_argument(
        '--report',
        type=str,
        metavar='FILE',
        help='Save processing report to JSON file'
    )
    parser.add_argument(
        '--no-remarks',
        action='store_true',
        help='Do not preserve HEADER/TITLE/CRYST1/REMARK lines'
    )
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Verbose output'
    )

    args = parser.parse_args()

    # Setup logging
    setup_logging(args.verbose)

    # Initialize concatenator
    concatenator = MultiModelConcatenator()

    if args.batch:
        # Batch processing mode
        output_dir = Path(args.output)
        output_dir.mkdir(parents=True, exist_ok=True)

        print(f"Batch processing {len(args.input)} file(s)...")
        print(f"Output directory: {output_dir}")
        print()

        results = []
        for input_file in args.input:
            input_path = Path(input_file)
            if not input_path.exists():
                print(f"✗ {input_path.name}: File not found")
                results.append({
                    'input_file': str(input_path),
                    'status': 'failed',
                    'error': 'File not found'
                })
                continue

            output_file = output_dir / f"{input_path.stem}_concatenated.pdb"

            print(f"Processing {input_path.name}...", end=' ')
            try:
                result = concatenator.concatenate_models(
                    str(input_path),
                    str(output_file),
                    preserve_remarks=not args.no_remarks
                )
                print(f"✓ {result['n_models']} models -> {result['total_atoms']} atoms")
                results.append(result)
            except Exception as e:
                print(f"✗ Failed: {e}")
                results.append({
                    'input_file': str(input_path),
                    'status': 'failed',
                    'error': str(e)
                })

        # Print summary
        print()
        success_count = sum(1 for r in results if r.get('status') == 'success')
        failed_count = len(results) - success_count
        print(f"Summary: {success_count} succeeded, {failed_count} failed")

        # Save report
        if args.report:
            with open(args.report, 'w') as f:
                json.dump(results, f, indent=2)
            print(f"Report saved to {args.report}")

    else:
        # Single file mode
        if len(args.input) != 1:
            print("Error: Single file mode requires exactly one input file", file=sys.stderr)
            print("       Use --batch for multiple files", file=sys.stderr)
            sys.exit(1)

        input_file = args.input[0]
        output_file = args.output

        # Check input exists
        if not Path(input_file).exists():
            print(f"Error: Input file not found: {input_file}", file=sys.stderr)
            sys.exit(1)

        print(f"Input:  {input_file}")
        print(f"Output: {output_file}")
        print()

        try:
            result = concatenator.concatenate_models(
                input_file,
                output_file,
                preserve_remarks=not args.no_remarks
            )

            print(f"✓ Successfully merged {result['n_models']} model(s)")
            print(f"  Total atoms: {result['total_atoms']}")
            print(f"  Atoms per model: {result['atoms_per_model']}")

            # Save report
            if args.report:
                with open(args.report, 'w') as f:
                    json.dump(result, f, indent=2)
                print(f"\nReport saved to {args.report}")

        except Exception as e:
            print(f"✗ Failed: {e}", file=sys.stderr)
            sys.exit(1)


if __name__ == '__main__':
    main()
