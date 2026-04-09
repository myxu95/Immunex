#!/usr/bin/env python3
"""
Batch trim PDB files by distance from peptide center.

Process multiple PDB files in parallel, removing residues beyond distance threshold.
"""
import sys
import argparse
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import json

# Import the single-file trimming function
sys.path.insert(0, str(Path(__file__).parent))
from trim_pdb_by_distance import trim_pdb_by_distance


def process_single_pdb(args):
    """Process single PDB file (for multiprocessing)."""
    input_pdb, output_dir, distance_cutoff = args

    try:
        input_path = Path(input_pdb)
        output_path = output_dir / f"{input_path.stem}_trimmed.pdb"

        print(f"\nProcessing: {input_path.name}")
        print("-" * 80)

        stats = trim_pdb_by_distance(
            str(input_pdb),
            str(output_path),
            distance_cutoff_nm=distance_cutoff,
            verbose=False  # Less verbose for batch processing
        )

        return {
            'status': 'success',
            'input': str(input_pdb),
            'output': str(output_path),
            'stats': stats
        }

    except Exception as e:
        print(f"ERROR processing {input_pdb}: {e}")
        return {
            'status': 'failed',
            'input': str(input_pdb),
            'error': str(e)
        }


def batch_trim_pdbs(input_dir, output_dir, distance_cutoff=8.0, pattern="*.pdb", max_workers=4):
    """
    Batch process multiple PDB files.

    Args:
        input_dir: Input directory containing PDB files
        output_dir: Output directory for trimmed PDBs
        distance_cutoff: Distance cutoff in nanometers
        pattern: File pattern for PDB files (default: *.pdb)
        max_workers: Number of parallel workers

    Returns:
        List of result dictionaries
    """
    input_path = Path(input_dir)
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Find all PDB files
    pdb_files = sorted(input_path.glob(pattern))

    if not pdb_files:
        print(f"No PDB files found in {input_dir} with pattern {pattern}")
        return []

    print("=" * 80)
    print("Batch PDB Trimming")
    print("=" * 80)
    print(f"Input directory:  {input_dir}")
    print(f"Output directory: {output_dir}")
    print(f"Distance cutoff:  {distance_cutoff} nm")
    print(f"Total files:      {len(pdb_files)}")
    print(f"Workers:          {max_workers}")
    print("=" * 80)

    # Prepare arguments for parallel processing
    args_list = [(str(pdb), output_path, distance_cutoff) for pdb in pdb_files]

    # Process in parallel
    results = []
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        future_to_pdb = {executor.submit(process_single_pdb, args): args[0]
                        for args in args_list}

        for future in as_completed(future_to_pdb):
            result = future.result()
            results.append(result)

    # Print summary
    print("\n" + "=" * 80)
    print("Batch Processing Summary")
    print("=" * 80)

    success_count = sum(1 for r in results if r['status'] == 'success')
    failed_count = sum(1 for r in results if r['status'] == 'failed')

    print(f"Total processed:    {len(results)}")
    print(f"Successful:         {success_count}")
    print(f"Failed:             {failed_count}")

    if success_count > 0:
        # Calculate statistics
        total_original = sum(r['stats']['total_residues'] for r in results if r['status'] == 'success')
        total_removed = sum(r['stats']['removed_residues'] for r in results if r['status'] == 'success')
        avg_retention = sum(r['stats']['kept_residues'] / r['stats']['total_residues'] * 100
                           for r in results if r['status'] == 'success') / success_count

        print(f"\nStatistics:")
        print(f"  Total residues:    {total_original}")
        print(f"  Residues removed:  {total_removed} ({total_removed/total_original*100:.1f}%)")
        print(f"  Avg retention:     {avg_retention:.1f}%")

    if failed_count > 0:
        print(f"\nFailed files:")
        for r in results:
            if r['status'] == 'failed':
                print(f"  {Path(r['input']).name}: {r['error']}")

    print("=" * 80)

    return results


def main():
    parser = argparse.ArgumentParser(
        description="Batch trim PDB files by distance from peptide center"
    )
    parser.add_argument(
        "input_dir",
        help="Input directory containing PDB files"
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output directory for trimmed PDB files"
    )
    parser.add_argument(
        "-d", "--distance",
        type=float,
        default=8.0,
        help="Distance cutoff in nanometers (default: 8.0 nm)"
    )
    parser.add_argument(
        "-p", "--pattern",
        default="*.pdb",
        help="File pattern for PDB files (default: *.pdb)"
    )
    parser.add_argument(
        "-w", "--workers",
        type=int,
        default=4,
        help="Number of parallel workers (default: 4)"
    )
    parser.add_argument(
        "--save-report",
        help="Save detailed report to JSON file"
    )

    args = parser.parse_args()

    # Run batch processing
    results = batch_trim_pdbs(
        input_dir=args.input_dir,
        output_dir=args.output,
        distance_cutoff=args.distance,
        pattern=args.pattern,
        max_workers=args.workers
    )

    # Save report if requested
    if args.save_report:
        report_path = Path(args.save_report)
        with open(report_path, 'w') as f:
            json.dump(results, f, indent=2)
        print(f"\nDetailed report saved to: {report_path}")

    # Exit code based on failures
    failed_count = sum(1 for r in results if r['status'] == 'failed')
    return 1 if failed_count > 0 else 0


if __name__ == "__main__":
    sys.exit(main())
