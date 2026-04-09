"""
Test script for new chain identification strategy (2026-01-22)

New Strategy:
1. Peptide: <=20 AA (definitive)
2. Beta2m: 90-110 AA (definitive)
3. Remaining 3 chains -> ANARCI -> identify TCR
4. After TCR identified -> remaining = HLA-alpha
"""

import sys
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from immunex.utils import IntelligentChainIdentifier, PDBSequenceExtractor

def test_chain_identification(pdb_file: str):
    """Test chain identification with new strategy."""

    print("="*80)
    print("Testing New Chain Identification Strategy")
    print("="*80)
    print(f"\nInput PDB: {pdb_file}\n")

    # Step 1: Extract sequences
    print("-"*80)
    print("Step 1: Extract sequences")
    print("-"*80)
    extractor = PDBSequenceExtractor()
    chains_data = extractor.extract_sequences_from_pdb(pdb_file)

    print(f"\nFound {len(chains_data)} chains:")
    for chain_id, info in sorted(chains_data.items(), key=lambda x: x[1]['length']):
        print(f"  Chain {chain_id}: {info['length']:3d} AA - {info['sequence'][:20]}...")

    # Step 2: Run intelligent identification
    print("\n" + "-"*80)
    print("Step 2: Intelligent Chain Identification (New Strategy)")
    print("-"*80)
    print("\nStrategy:")
    print("  1. Peptide: <=20 AA (definitive)")
    print("  2. Beta2m: 90-110 AA (definitive)")
    print("  3. Remaining 3 chains -> ANARCI")
    print("  4. After TCR identified -> remaining = HLA-alpha\n")

    identifier = IntelligentChainIdentifier(use_anarci=True)
    identifications = identifier.identify_chains(pdb_file)

    # Step 3: Show results
    print("\n" + "-"*80)
    print("Step 3: Identification Results")
    print("-"*80)
    print()

    # Sort by chain type order
    type_order = ['peptide', 'beta2m', 'TCR_alpha', 'TCR_beta', 'HLA_alpha', 'unknown']
    sorted_ids = sorted(
        identifications.items(),
        key=lambda x: (type_order.index(x[1].chain_type), x[0])
    )

    for chain_id, info in sorted_ids:
        confidence_str = f"{info.confidence:.2f}"
        print(f"  Chain {chain_id}: {info.chain_type:12s} ({info.length:3d} AA, confidence={confidence_str})")

    # Step 4: Create standardization mapping
    print("\n" + "-"*80)
    print("Step 4: Standard Chain Mapping")
    print("-"*80)
    print()

    mapping = identifier.create_standardization_mapping(identifications)

    print("Standardization mapping:")
    for old_id, new_id in sorted(mapping.items()):
        chain_type = identifications[old_id].chain_type
        print(f"  {old_id} -> {new_id} ({chain_type})")

    # Step 5: Validation
    print("\n" + "-"*80)
    print("Step 5: Validation")
    print("-"*80)
    print()

    expected_types = {'peptide', 'beta2m', 'TCR_alpha', 'TCR_beta', 'HLA_alpha'}
    identified_types = {info.chain_type for info in identifications.values()}

    if expected_types == identified_types:
        print("✓ All 5 chain types identified correctly!")
    else:
        missing = expected_types - identified_types
        extra = identified_types - expected_types
        if missing:
            print(f"✗ Missing types: {missing}")
        if extra:
            print(f"✗ Unexpected types: {extra}")

    # Check confidence
    low_confidence = [
        (cid, info.chain_type, info.confidence)
        for cid, info in identifications.items()
        if info.confidence < 0.8
    ]

    if low_confidence:
        print("\n⚠ Low confidence identifications:")
        for cid, ctype, conf in low_confidence:
            print(f"  Chain {cid} ({ctype}): confidence={conf:.2f}")
    else:
        print("\n✓ All identifications have high confidence (>=0.8)")

    print("\n" + "="*80)
    print("Test Complete")
    print("="*80)
    print()

    return identifications, mapping


if __name__ == "__main__":
    # Get PDB file from command line or use default
    if len(sys.argv) > 1:
        test_pdb = sys.argv[1]
    else:
        test_pdb = "test_data/test_complex.pdb"

    if not Path(test_pdb).exists():
        print(f"Error: Test PDB not found: {test_pdb}")
        print("\nCreating synthetic test PDB...")

        # Create test data directory
        Path("test_data").mkdir(exist_ok=True)

        # Create synthetic PDB with 5 chains
        with open(test_pdb, 'w') as f:
            # Chain V: 50 residues (should not match any pattern)
            for i in range(1, 51):
                f.write(f"ATOM  {i:5d}  CA  ALA V{i:4d}       0.000   0.000   0.000  1.00  0.00           C\n")

            # Chain W: 45 residues
            for i in range(1, 46):
                f.write(f"ATOM  {50+i:5d}  CA  ILE W{i:4d}       0.000   0.000   0.000  1.00  0.00           C\n")

            # Chain X: 10 residues (peptide)
            for i in range(1, 11):
                f.write(f"ATOM  {95+i:5d}  CA  VAL X{i:4d}       0.000   0.000   0.000  1.00  0.00           C\n")

            # Chain Y: 30 residues
            for i in range(1, 31):
                f.write(f"ATOM  {105+i:5d}  CA  GLY Y{i:4d}       0.000   0.000   0.000  1.00  0.00           C\n")

            # Chain Z: 40 residues
            for i in range(1, 41):
                f.write(f"ATOM  {135+i:5d}  CA  LEU Z{i:4d}       0.000   0.000   0.000  1.00  0.00           C\n")

            f.write("END\n")

        print(f"Created synthetic test PDB: {test_pdb}\n")

    # Run test
    identifications, mapping = test_chain_identification(test_pdb)
