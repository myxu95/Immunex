#!/usr/bin/env python3
"""
End-to-end test of intelligent chain standardization with real PDB file
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))


def test_complete_workflow():
    """Test complete workflow with test PDB."""
    print("="*80)
    print("End-to-End Test: Intelligent Chain Standardization")
    print("="*80)

    from aftermd.utils import (
        PDBSequenceExtractor,
        IntelligentChainIdentifier,
        PDBChainStandardizer
    )

    test_pdb = Path("test_data/test_complex.pdb")

    if not test_pdb.exists():
        print(f"✗ Test PDB not found: {test_pdb}")
        return False

    print(f"\n✓ Test PDB found: {test_pdb}")

    # Step 1: Extract sequences
    print("\n" + "-"*80)
    print("Step 1: Extract Sequences")
    print("-"*80)

    try:
        extractor = PDBSequenceExtractor()
        chains = extractor.extract_sequences_from_pdb(str(test_pdb))

        print(f"\n✓ Extracted {len(chains)} chains:")
        for chain_id in sorted(chains.keys()):
            info = chains[chain_id]
            print(f"  Chain {chain_id}: {info['length']:3d} residues - {info['sequence'][:20]}...")

    except Exception as e:
        print(f"✗ Failed: {e}")
        import traceback
        traceback.print_exc()
        return False

    # Step 2: Length-based identification (legacy)
    print("\n" + "-"*80)
    print("Step 2: Length-based Standardization (Legacy)")
    print("-"*80)

    try:
        std_legacy = PDBChainStandardizer(use_intelligent_identification=False)
        result_legacy = std_legacy.process_single(
            input_pdb=str(test_pdb),
            output_pdb="test_data/test_complex_legacy.pdb",
            task_name="test_legacy",
            skip_if_standard=False
        )

        print(f"\n✓ Status: {result_legacy.status}")
        print(f"  Mapping: {result_legacy.chain_mapping}")

    except Exception as e:
        print(f"✗ Failed: {e}")
        import traceback
        traceback.print_exc()
        return False

    # Step 3: Intelligent identification
    print("\n" + "-"*80)
    print("Step 3: Intelligent Chain Identification")
    print("-"*80)

    try:
        identifier = IntelligentChainIdentifier(use_anarci=True)
        identifications = identifier.identify_chains(str(test_pdb))

        print(f"\n✓ Identified {len(identifications)} chains:")
        for chain_id in sorted(identifications.keys()):
            info = identifications[chain_id]
            print(f"  Chain {chain_id}: {info.chain_type:15s} "
                  f"({info.length:3d} aa, confidence={info.confidence:.2f})")

        # Create mapping
        mapping = identifier.create_standardization_mapping(identifications)
        print(f"\n✓ Proposed mapping:")
        for old_id in sorted(mapping.keys()):
            new_id = mapping[old_id]
            chain_type = identifications[old_id].chain_type
            print(f"  {old_id} → {new_id} ({chain_type})")

    except Exception as e:
        print(f"✗ Failed: {e}")
        import traceback
        traceback.print_exc()
        return False

    # Step 4: Intelligent standardization
    print("\n" + "-"*80)
    print("Step 4: Intelligent Standardization")
    print("-"*80)

    try:
        std_intelligent = PDBChainStandardizer(use_intelligent_identification=True)
        result_intelligent = std_intelligent.process_single(
            input_pdb=str(test_pdb),
            output_pdb="test_data/test_complex_intelligent.pdb",
            task_name="test_intelligent",
            skip_if_standard=False
        )

        print(f"\n✓ Status: {result_intelligent.status}")
        print(f"  Mapping: {result_intelligent.chain_mapping}")
        print(f"  Processing time: {result_intelligent.processing_time:.3f}s")

    except Exception as e:
        print(f"✗ Failed: {e}")
        import traceback
        traceback.print_exc()
        return False

    # Step 5: Compare methods
    print("\n" + "-"*80)
    print("Step 5: Compare Methods")
    print("-"*80)

    print(f"\nLegacy mapping:      {result_legacy.chain_mapping}")
    print(f"Intelligent mapping: {result_intelligent.chain_mapping}")

    if result_legacy.chain_mapping == result_intelligent.chain_mapping:
        print("\n✓ Both methods produced identical mapping")
    else:
        print("\n⚠ Different mappings detected!")
        print("  This is expected if TCR chains were swapped in legacy method")

        # Show differences
        all_chains = set(result_legacy.chain_mapping.keys()) | set(result_intelligent.chain_mapping.keys())
        for chain in sorted(all_chains):
            legacy_target = result_legacy.chain_mapping.get(chain, '?')
            intelligent_target = result_intelligent.chain_mapping.get(chain, '?')

            if legacy_target != intelligent_target:
                print(f"  Chain {chain}: {legacy_target} (legacy) vs {intelligent_target} (intelligent)")

    print("\n" + "="*80)
    print("✓ End-to-End Test Completed Successfully!")
    print("="*80)

    return True


if __name__ == "__main__":
    success = test_complete_workflow()
    sys.exit(0 if success else 1)
