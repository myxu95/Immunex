#!/usr/bin/env python3
"""
Test script for chain standardization enhancements

Tests:
1. ChainBasedIndexGenerator with PDB input
2. Chain reading and validation
3. Index file generation

Author: Immunex Development Team
Date: 2026-01-21
"""

import sys
import os
import tempfile
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from immunex.utils import (
    ChainBasedIndexGenerator,
    generate_peptide_index_from_pdb
)


def test_chain_based_index_generator():
    """Test ChainBasedIndexGenerator with standardized PDB"""
    print("=" * 80)
    print("Test: ChainBasedIndexGenerator")
    print("=" * 80)

    # Use a standardized PDB file
    test_pdb = project_root / "input/standardized_pdbs/1ao7_standardized.pdb"

    if not test_pdb.exists():
        print(f"Test PDB not found: {test_pdb}")
        return False

    print(f"\nTest PDB: {test_pdb}")

    # Create temporary output directory
    with tempfile.TemporaryDirectory() as temp_dir:
        print(f"Temp output dir: {temp_dir}")

        try:
            # Test 1: Initialize generator
            print("\n[Test 1] Initialize ChainBasedIndexGenerator")
            generator = ChainBasedIndexGenerator(str(test_pdb))
            print(f"  Input type: {generator.input_type}")
            print(f"  PDB file: {generator.pdb_file}")
            assert generator.input_type == 'pdb', "Should detect PDB input type"
            print("  PASS")

            # Test 2: Read chain info
            print("\n[Test 2] Read Chain C (peptide) information")
            chain_info = generator._read_chain_info(str(test_pdb), 'C')
            print(f"  Chain ID: {chain_info['chain_id']}")
            print(f"  Residues: {chain_info['num_residues']}")
            print(f"  Atoms: {chain_info['num_atoms']}")
            print(f"  Residue range: {chain_info['residues'][0]} - {chain_info['residues'][-1]}")

            assert chain_info['chain_id'] == 'C', "Should read Chain C"
            assert chain_info['num_residues'] > 0, "Should have residues"
            assert chain_info['num_atoms'] > 0, "Should have atoms"
            print("  PASS")

            # Test 3: Generate index file
            print("\n[Test 3] Generate peptide index file")
            index_file = generator.generate_peptide_index(
                output_dir=temp_dir,
                peptide_chain_id='C',
                index_name="test_peptide.ndx"
            )
            print(f"  Index file: {index_file}")

            # Verify index file exists
            assert Path(index_file).exists(), "Index file should exist"
            print("  PASS")

            # Test 4: Verify index file content
            print("\n[Test 4] Verify index file content")
            with open(index_file) as f:
                content = f.read()
                print(f"  File size: {len(content)} bytes")

                # Check for peptide group
                assert "[ Peptide_Chain_C ]" in content, "Should contain peptide group"
                print("  Contains peptide group: YES")

                # Count atom numbers in Peptide_Chain_C group only
                lines = content.split('\n')
                in_peptide_group = False
                peptide_atoms = []

                for line in lines:
                    if '[ Peptide_Chain_C ]' in line:
                        in_peptide_group = True
                        continue
                    elif line.startswith('['):
                        in_peptide_group = False
                        continue

                    if in_peptide_group and line.strip():
                        peptide_atoms.extend(line.split())

                peptide_atom_count = len(peptide_atoms)
                print(f"  Atoms in Peptide_Chain_C group: {peptide_atom_count}")
                assert peptide_atom_count == chain_info['num_atoms'], "Peptide atom count should match"
                print("  PASS")

            # Test 5: Convenience function
            print("\n[Test 5] Test convenience function")
            index_file2 = generate_peptide_index_from_pdb(
                pdb_file=str(test_pdb),
                output_dir=temp_dir,
                peptide_chain_id='C'
            )
            assert Path(index_file2).exists(), "Convenience function should work"
            print(f"  Generated: {Path(index_file2).name}")
            print("  PASS")

            print("\n" + "=" * 80)
            print("All tests PASSED")
            print("=" * 80)
            return True

        except Exception as e:
            print(f"\n[ERROR] Test failed: {e}")
            import traceback
            traceback.print_exc()
            return False


def test_multiple_chains():
    """Test reading different chains from standardized PDB"""
    print("\n" + "=" * 80)
    print("Test: Multiple Chain Reading")
    print("=" * 80)

    test_pdb = project_root / "input/standardized_pdbs/1ao7_standardized.pdb"

    if not test_pdb.exists():
        print(f"Test PDB not found: {test_pdb}")
        return False

    generator = ChainBasedIndexGenerator(str(test_pdb))

    # Expected chains in standardized PDB
    expected_chains = {
        'A': 'HLA-alpha (longest)',
        'B': 'beta2-microglobulin',
        'C': 'Peptide (shortest)',
        'D': 'TCR-alpha',
        'E': 'TCR-beta'
    }

    print("\nReading all chains:")
    for chain_id, description in expected_chains.items():
        try:
            chain_info = generator._read_chain_info(str(test_pdb), chain_id)
            print(f"  Chain {chain_id} ({description}):")
            print(f"    Residues: {chain_info['num_residues']}")
            print(f"    Atoms: {chain_info['num_atoms']}")
        except ValueError as e:
            print(f"  Chain {chain_id}: Not found (may be missing in this PDB)")

    print("\n" + "=" * 80)
    print("Multiple chain reading test completed")
    print("=" * 80)
    return True


def test_error_handling():
    """Test error handling for invalid inputs"""
    print("\n" + "=" * 80)
    print("Test: Error Handling")
    print("=" * 80)

    # Test 1: Non-existent file
    print("\n[Test 1] Non-existent file")
    try:
        generator = ChainBasedIndexGenerator("/non/existent/file.pdb")
        print("  FAIL: Should raise FileNotFoundError")
        return False
    except FileNotFoundError:
        print("  PASS: Correctly raised FileNotFoundError")

    # Test 2: Invalid file type
    print("\n[Test 2] Invalid file type")
    with tempfile.NamedTemporaryFile(suffix='.txt', delete=False) as tmp:
        tmp_path = tmp.name

    try:
        generator = ChainBasedIndexGenerator(tmp_path)
        print("  FAIL: Should raise ValueError for invalid file type")
        os.unlink(tmp_path)
        return False
    except ValueError as e:
        print(f"  PASS: Correctly raised ValueError ({e})")
        os.unlink(tmp_path)

    # Test 3: Non-existent chain
    print("\n[Test 3] Non-existent chain ID")
    test_pdb = project_root / "input/standardized_pdbs/1ao7_standardized.pdb"
    if test_pdb.exists():
        generator = ChainBasedIndexGenerator(str(test_pdb))
        try:
            with tempfile.TemporaryDirectory() as temp_dir:
                generator.generate_peptide_index(
                    output_dir=temp_dir,
                    peptide_chain_id='Z',  # Non-existent chain
                    index_name="test.ndx"
                )
            print("  FAIL: Should raise ValueError for non-existent chain")
            return False
        except ValueError as e:
            print(f"  PASS: Correctly raised ValueError")

    print("\n" + "=" * 80)
    print("Error handling tests PASSED")
    print("=" * 80)
    return True


if __name__ == "__main__":
    print("\n")
    print("=" * 80)
    print("Immunex Chain Standardization - Test Suite")
    print("=" * 80)

    all_passed = True

    # Run tests
    all_passed &= test_chain_based_index_generator()
    all_passed &= test_multiple_chains()
    all_passed &= test_error_handling()

    # Summary
    print("\n" + "=" * 80)
    if all_passed:
        print("TEST SUITE: ALL TESTS PASSED")
    else:
        print("TEST SUITE: SOME TESTS FAILED")
    print("=" * 80)
    print()

    sys.exit(0 if all_passed else 1)
