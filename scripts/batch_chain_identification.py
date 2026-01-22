#!/usr/bin/env python3
"""
Batch Chain Identification for pHLA-TCR Complexes

Applies new identification strategy to multiple PDB files:
1. Peptide: <=20 AA (definitive)
2. Beta2m: 90-110 AA (definitive)
3. Remaining 3 chains -> ANARCI -> TCR identification
4. Remaining chain -> HLA-alpha (by elimination)

Usage:
    python scripts/batch_chain_identification.py <pdb_directory> [output_csv]

Author: AfterMD Development Team
Date: 2026-01-22
"""

import sys
import csv
from pathlib import Path
from typing import List, Dict, Optional
import logging

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from aftermd.utils import IntelligentChainIdentifier, PDBSequenceExtractor

# Setup logging
logging.basicConfig(
    level=logging.WARNING,  # Only show warnings and errors
    format='%(message)s'
)
logger = logging.getLogger(__name__)


def identify_pdb_chains(pdb_file: Path) -> Dict:
    """
    Identify chains in a single PDB file.

    Returns:
        Dictionary with identification results
    """
    try:
        identifier = IntelligentChainIdentifier(use_anarci=True)
        identifications = identifier.identify_chains(str(pdb_file))
        mapping = identifier.create_standardization_mapping(identifications)

        result = {
            'pdb_file': pdb_file.name,
            'pdb_id': pdb_file.stem.split('_')[0],
            'status': 'success',
            'num_chains': len(identifications),
        }

        # Add chain information
        for chain_id, info in identifications.items():
            result[f'chain_{chain_id}_type'] = info.chain_type
            result[f'chain_{chain_id}_length'] = info.length
            result[f'chain_{chain_id}_confidence'] = f"{info.confidence:.2f}"
            result[f'chain_{chain_id}_mapped_to'] = mapping.get(chain_id, chain_id)

        # Validation
        expected_types = {'peptide', 'beta2m', 'TCR_alpha', 'TCR_beta', 'HLA_alpha'}
        identified_types = {info.chain_type for info in identifications.values()}

        if expected_types == identified_types:
            result['validation'] = 'PASS'
        else:
            missing = expected_types - identified_types
            result['validation'] = f'FAIL (missing: {missing})'

        # Low confidence warning
        low_conf = [
            (cid, info.chain_type, info.confidence)
            for cid, info in identifications.items()
            if info.confidence < 0.8
        ]

        if low_conf:
            result['warnings'] = f'Low confidence: {low_conf}'
        else:
            result['warnings'] = ''

        return result

    except Exception as e:
        return {
            'pdb_file': pdb_file.name,
            'pdb_id': pdb_file.stem.split('_')[0],
            'status': 'error',
            'error_message': str(e)
        }


def batch_identify(pdb_directory: str, output_csv: Optional[str] = None) -> List[Dict]:
    """
    Batch identify chains in all PDB files in directory.

    Args:
        pdb_directory: Directory containing PDB files
        output_csv: Optional output CSV file path

    Returns:
        List of identification results
    """
    pdb_dir = Path(pdb_directory)

    if not pdb_dir.exists():
        raise FileNotFoundError(f"Directory not found: {pdb_directory}")

    # Find all PDB files
    pdb_files = sorted(pdb_dir.glob('*.pdb'))

    if not pdb_files:
        raise ValueError(f"No PDB files found in: {pdb_directory}")

    print("="*80)
    print("Batch Chain Identification")
    print("="*80)
    print(f"\nDirectory: {pdb_directory}")
    print(f"Found {len(pdb_files)} PDB files\n")
    print("Strategy:")
    print("  1. Peptide: <=20 AA (definitive)")
    print("  2. Beta2m: 90-110 AA (definitive)")
    print("  3. Remaining 3 chains -> ANARCI")
    print("  4. Remaining -> HLA-alpha (elimination)")
    print("\n" + "="*80 + "\n")

    results = []

    for i, pdb_file in enumerate(pdb_files, 1):
        print(f"[{i}/{len(pdb_files)}] Processing {pdb_file.name}...", end=' ')

        result = identify_pdb_chains(pdb_file)
        results.append(result)

        if result['status'] == 'success':
            validation = result.get('validation', 'UNKNOWN')
            if validation == 'PASS':
                print(f"✓ {validation}")
            else:
                print(f"⚠ {validation}")
        else:
            print(f"✗ ERROR: {result.get('error_message', 'Unknown error')}")

    # Summary
    print("\n" + "="*80)
    print("Summary")
    print("="*80)

    success_count = sum(1 for r in results if r['status'] == 'success')
    pass_count = sum(1 for r in results if r.get('validation') == 'PASS')
    fail_count = success_count - pass_count
    error_count = len(results) - success_count

    print(f"\nTotal: {len(results)} files")
    print(f"  ✓ Success: {success_count}")
    print(f"    - Validation PASS: {pass_count}")
    print(f"    - Validation FAIL: {fail_count}")
    print(f"  ✗ Errors: {error_count}")

    # Write CSV if requested
    if output_csv:
        output_path = Path(output_csv)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # Collect all possible column names
        all_keys = set()
        for result in results:
            all_keys.update(result.keys())

        # Sort columns for better readability
        priority_cols = ['pdb_file', 'pdb_id', 'status', 'validation', 'num_chains']
        other_cols = sorted(all_keys - set(priority_cols))
        fieldnames = priority_cols + other_cols

        with open(output_path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(results)

        print(f"\n✓ Results saved to: {output_csv}")

    print("\n" + "="*80 + "\n")

    return results


def main():
    """Main entry point."""
    if len(sys.argv) < 2:
        print("Usage: python scripts/batch_chain_identification.py <pdb_directory> [output_csv]")
        print("\nExample:")
        print("  python scripts/batch_chain_identification.py input/standardizedpdbs/protein_only")
        print("  python scripts/batch_chain_identification.py input/standardizedpdbs/protein_only results.csv")
        sys.exit(1)

    pdb_directory = sys.argv[1]
    output_csv = sys.argv[2] if len(sys.argv) > 2 else None

    try:
        results = batch_identify(pdb_directory, output_csv)

        # Exit with error code if any failures
        if any(r['status'] == 'error' or r.get('validation') != 'PASS' for r in results):
            sys.exit(1)

    except Exception as e:
        print(f"\n✗ Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
