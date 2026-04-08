#!/usr/bin/env python3
"""
Example usage of Intelligent Chain Standardization

Demonstrates the new ANARCI-based chain identification system for pHLA-TCR complexes.
"""

from pathlib import Path
from immunex.utils import (
    PDBSequenceExtractor,
    IntelligentChainIdentifier,
    PDBChainStandardizer
)


def example1_extract_sequences():
    """Example 1: Extract sequences from PDB to CSV."""
    print("="*80)
    print("Example 1: Extract Chain Sequences")
    print("="*80)

    extractor = PDBSequenceExtractor()

    # Extract sequences from single PDB
    chains = extractor.extract_and_save(
        pdb_file="input/1ao7.pdb",
        output_csv="output/1ao7_sequences.csv",
        task_name="1ao7"
    )

    print(f"\nExtracted {len(chains)} chains:")
    for chain_id, info in chains.items():
        print(f"  Chain {chain_id}: {info['length']} residues")
        print(f"    Sequence: {info['sequence'][:50]}...")


def example2_intelligent_identification():
    """Example 2: Identify chains using ANARCI."""
    print("\n" + "="*80)
    print("Example 2: Intelligent Chain Identification")
    print("="*80)

    identifier = IntelligentChainIdentifier(use_anarci=True)

    # Identify chains in PDB file
    identifications = identifier.identify_chains("input/1ao7.pdb")

    print(f"\nIdentified {len(identifications)} chains:")
    for chain_id, info in identifications.items():
        print(f"  Chain {chain_id}: {info.chain_type}")
        print(f"    Length: {info.length} aa")
        print(f"    Confidence: {info.confidence:.2f}")

    # Create standardization mapping
    mapping = identifier.create_standardization_mapping(identifications)
    print(f"\nProposed chain mapping:")
    for old_id, new_id in sorted(mapping.items()):
        chain_type = identifications[old_id].chain_type
        print(f"  {old_id} → {new_id} ({chain_type})")


def example3_intelligent_standardization():
    """Example 3: Standardize PDB with intelligent identification."""
    print("\n" + "="*80)
    print("Example 3: Intelligent Chain Standardization")
    print("="*80)

    # Create standardizer with intelligent mode enabled
    standardizer = PDBChainStandardizer(use_intelligent_identification=True)

    # Process PDB file
    result = standardizer.process_single(
        input_pdb="input/1ao7.pdb",
        output_pdb="output/1ao7_standardized_intelligent.pdb",
        task_name="1ao7"
    )

    print(f"\nStandardization result:")
    print(f"  Status: {result.status}")
    print(f"  Chains: {result.num_chains}")
    print(f"  Mapping: {result.chain_mapping}")
    print(f"  Time: {result.processing_time:.2f}s")

    if result.status == "OK":
        print(f"  Output: {result.output_file}")


def example4_batch_with_sequences():
    """Example 4: Batch extract sequences from multiple PDBs."""
    print("\n" + "="*80)
    print("Example 4: Batch Sequence Extraction")
    print("="*80)

    extractor = PDBSequenceExtractor()

    # Extract from multiple PDB files
    pdb_files = [
        "input/1ao7.pdb",
        "input/2vlk.pdb",
        "input/3kxf.pdb"
    ]

    extractor.batch_extract(
        pdb_files=pdb_files,
        output_csv="output/all_sequences.csv"
    )

    print(f"\nExtracted sequences from {len(pdb_files)} PDB files")
    print("Saved to: output/all_sequences.csv")


def example5_comparison_modes():
    """Example 5: Compare length-based vs intelligent identification."""
    print("\n" + "="*80)
    print("Example 5: Compare Identification Methods")
    print("="*80)

    pdb_file = "input/1ao7.pdb"

    # Method 1: Length-based (legacy)
    print("\nMethod 1: Length-based standardization")
    standardizer_length = PDBChainStandardizer(use_intelligent_identification=False)
    result_length = standardizer_length.process_single(
        input_pdb=pdb_file,
        output_pdb="output/1ao7_length_based.pdb",
        task_name="1ao7_length"
    )
    print(f"  Mapping: {result_length.chain_mapping}")

    # Method 2: Intelligent (ANARCI)
    print("\nMethod 2: Intelligent (ANARCI) standardization")
    standardizer_intelligent = PDBChainStandardizer(use_intelligent_identification=True)
    result_intelligent = standardizer_intelligent.process_single(
        input_pdb=pdb_file,
        output_pdb="output/1ao7_intelligent.pdb",
        task_name="1ao7_intelligent"
    )
    print(f"  Mapping: {result_intelligent.chain_mapping}")

    # Compare
    if result_length.chain_mapping != result_intelligent.chain_mapping:
        print("\n⚠️  Warning: Different mappings detected!")
        print("  This may indicate TCR alpha/beta chain swap in length-based method")
    else:
        print("\n✓ Both methods produced identical mapping")


def example6_handle_edge_cases():
    """Example 6: Handle edge cases gracefully."""
    print("\n" + "="*80)
    print("Example 6: Edge Cases and Error Handling")
    print("="*80)

    identifier = IntelligentChainIdentifier(use_anarci=True)

    # Case 1: File with non-standard chain count
    print("\nCase 1: Non-standard chain count")
    try:
        identifications = identifier.identify_chains("input/non_standard.pdb")
        print(f"  Identified {len(identifications)} chains (expected 5)")
    except FileNotFoundError:
        print("  File not found (expected for this example)")

    # Case 2: ANARCI not available (fallback)
    print("\nCase 2: ANARCI fallback")
    identifier_no_anarci = IntelligentChainIdentifier(use_anarci=False)
    identifications = identifier_no_anarci.identify_chains("input/1ao7.pdb")
    print(f"  Used length-based fallback for {len(identifications)} chains")

    # Case 3: Already standardized PDB
    print("\nCase 3: Already standardized PDB")
    standardizer = PDBChainStandardizer(use_intelligent_identification=True)
    result = standardizer.process_single(
        input_pdb="output/1ao7_standardized_intelligent.pdb",
        output_pdb="output/1ao7_check.pdb",
        task_name="1ao7_check",
        skip_if_standard=True
    )
    print(f"  Status: {result.status}")


if __name__ == "__main__":
    print("Intelligent Chain Standardization - Usage Examples")
    print("="*80)
    print()

    # Create output directory
    Path("output").mkdir(exist_ok=True)

    # Run examples
    try:
        example1_extract_sequences()
    except Exception as e:
        print(f"Example 1 failed: {e}")

    try:
        example2_intelligent_identification()
    except Exception as e:
        print(f"Example 2 failed: {e}")

    try:
        example3_intelligent_standardization()
    except Exception as e:
        print(f"Example 3 failed: {e}")

    try:
        example4_batch_with_sequences()
    except Exception as e:
        print(f"Example 4 failed: {e}")

    try:
        example5_comparison_modes()
    except Exception as e:
        print(f"Example 5 failed: {e}")

    try:
        example6_handle_edge_cases()
    except Exception as e:
        print(f"Example 6 failed: {e}")

    print("\n" + "="*80)
    print("All examples completed!")
    print("="*80)
