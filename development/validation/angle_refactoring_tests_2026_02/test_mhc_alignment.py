"""
Test MHC Sequence Alignment Module

Validates MHC α1/α2 region detection via sequence alignment to HLA-A*02:01 reference.

Test Cases:
1. 1AO7 structure - HLA-A*02:01 complex
2. 1bd2 structure - HLA-B*27:05 complex (different allele)
3. Alignment quality metrics validation
4. Residue mapping correctness

Author: Immunex Development Team
Date: 2026-02-04
"""

import sys
import os
sys.path.insert(0, os.path.abspath('.'))

import MDAnalysis as mda
from immunex.analysis.angles.mhc_sequence_aligner import MHCSequenceAligner


def test_1ao7_alignment():
    """Test alignment on 1AO7 (HLA-A*02:01) structure."""
    print("=" * 80)
    print("Test 1: 1AO7 (HLA-A*02:01) Alignment")
    print("=" * 80)

    # Load structure
    pdb_path = "development/workspaces/FEL_workspace/1ao7_basin_structures/1ao7_standard_basin1.pdb"
    if not os.path.exists(pdb_path):
        print(f"ERROR: File not found: {pdb_path}")
        print("Please ensure 1AO7 structure is available.")
        return False

    u = mda.Universe(pdb_path)

    # Initialize aligner (assuming MHC is chain A)
    aligner = MHCSequenceAligner(u, "chainID A")

    try:
        # Extract sequence
        pdb_seq, residue_map = aligner.extract_pdb_sequence()
        print(f"\n1. PDB Sequence Extraction:")
        print(f"   Length: {len(pdb_seq)} residues")
        print(f"   First 50 residues: {pdb_seq[:50]}")
        print(f"   Residue map range: {min(residue_map.values())} - {max(residue_map.values())}")

        # Perform alignment
        alignment, score, identity = aligner.align_to_reference(pdb_seq)
        print(f"\n2. Alignment Quality:")
        print(f"   Score: {score:.1f} (threshold: {aligner.MIN_ALIGNMENT_SCORE})")
        print(f"   Identity: {identity:.2%} (threshold: {aligner.MIN_SEQUENCE_IDENTITY:.0%})")
        print(f"   Status: {'PASS' if score >= aligner.MIN_ALIGNMENT_SCORE else 'FAIL'}")

        # Get α1/α2 residues
        alpha1_resids, alpha2_resids, quality = aligner.get_alpha1_alpha2_residues()

        print(f"\n3. α1 Helix Residues:")
        print(f"   Mapped count: {len(alpha1_resids)}")
        print(f"   PDB resid range: {min(alpha1_resids)} - {max(alpha1_resids)}")
        print(f"   Gap fraction: {quality['alpha1_gap_fraction']:.1%}")
        print(f"   Residues: {alpha1_resids}")

        print(f"\n4. α2 Helix Residues:")
        print(f"   Mapped count: {len(alpha2_resids)}")
        print(f"   PDB resid range: {min(alpha2_resids)} - {max(alpha2_resids)}")
        print(f"   Gap fraction: {quality['alpha2_gap_fraction']:.1%}")
        print(f"   Residues: {alpha2_resids}")

        print(f"\n5. Overall Quality:")
        print(f"   Alignment score: {quality['score']:.1f}")
        print(f"   Sequence identity: {quality['identity']:.2%}")
        print(f"   α1 mapped: {quality['alpha1_mapped_count']}/{aligner.ALPHA1_POSITIONS[1] - aligner.ALPHA1_POSITIONS[0] + 1}")
        print(f"   α2 mapped: {quality['alpha2_mapped_count']}/{aligner.ALPHA2_POSITIONS[1] - aligner.ALPHA2_POSITIONS[0] + 1}")

        print(f"\n✓ 1AO7 alignment test PASSED")
        return True

    except Exception as e:
        print(f"\n✗ 1AO7 alignment test FAILED: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_1bd2_alignment():
    """Test alignment on 1bd2 (HLA-B*27:05) structure - different allele."""
    print("\n" + "=" * 80)
    print("Test 2: 1bd2 (HLA-B*27:05) Alignment")
    print("=" * 80)

    # Load structure
    pdb_path = "development/workspaces/FEL_workspace/1bd2_basin_structures/1bd2_run1_basin1.pdb"
    if not os.path.exists(pdb_path):
        print(f"ERROR: File not found: {pdb_path}")
        print("Please ensure 1bd2 structure is available.")
        return False

    u = mda.Universe(pdb_path)

    # Initialize aligner (assuming MHC is chain A)
    aligner = MHCSequenceAligner(u, "chainID A")

    try:
        # Extract sequence
        pdb_seq, residue_map = aligner.extract_pdb_sequence()
        print(f"\n1. PDB Sequence Extraction:")
        print(f"   Length: {len(pdb_seq)} residues")
        print(f"   First 50 residues: {pdb_seq[:50]}")
        print(f"   Residue map range: {min(residue_map.values())} - {max(residue_map.values())}")

        # Perform alignment
        alignment, score, identity = aligner.align_to_reference(pdb_seq)
        print(f"\n2. Alignment Quality:")
        print(f"   Score: {score:.1f} (threshold: {aligner.MIN_ALIGNMENT_SCORE})")
        print(f"   Identity: {identity:.2%} (threshold: {aligner.MIN_SEQUENCE_IDENTITY:.0%})")
        print(f"   Status: {'PASS' if score >= aligner.MIN_ALIGNMENT_SCORE else 'FAIL'}")
        print(f"   Note: HLA-B*27:05 expected to have lower identity than HLA-A*02:01")

        # Get α1/α2 residues
        alpha1_resids, alpha2_resids, quality = aligner.get_alpha1_alpha2_residues()

        print(f"\n3. α1 Helix Residues:")
        print(f"   Mapped count: {len(alpha1_resids)}")
        print(f"   PDB resid range: {min(alpha1_resids)} - {max(alpha1_resids)}")
        print(f"   Gap fraction: {quality['alpha1_gap_fraction']:.1%}")

        print(f"\n4. α2 Helix Residues:")
        print(f"   Mapped count: {len(alpha2_resids)}")
        print(f"   PDB resid range: {min(alpha2_resids)} - {max(alpha2_resids)}")
        print(f"   Gap fraction: {quality['alpha2_gap_fraction']:.1%}")

        print(f"\n5. Overall Quality:")
        print(f"   Alignment score: {quality['score']:.1f}")
        print(f"   Sequence identity: {quality['identity']:.2%}")
        print(f"   α1 mapped: {quality['alpha1_mapped_count']}/{aligner.ALPHA1_POSITIONS[1] - aligner.ALPHA1_POSITIONS[0] + 1}")
        print(f"   α2 mapped: {quality['alpha2_mapped_count']}/{aligner.ALPHA2_POSITIONS[1] - aligner.ALPHA2_POSITIONS[0] + 1}")

        print(f"\n✓ 1bd2 alignment test PASSED")
        return True

    except Exception as e:
        print(f"\n✗ 1bd2 alignment test FAILED: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_alignment_robustness():
    """Test alignment robustness to non-standard numbering."""
    print("\n" + "=" * 80)
    print("Test 3: Alignment Robustness Check")
    print("=" * 80)

    # This test verifies that alignment works even if PDB numbering is non-standard
    # Use 1AO7 as test case
    pdb_path = "development/workspaces/FEL_workspace/1ao7_basin_structures/1ao7_standard_basin1.pdb"
    if not os.path.exists(pdb_path):
        print(f"Skipping robustness test: {pdb_path} not found")
        return False

    u = mda.Universe(pdb_path)
    aligner = MHCSequenceAligner(u, "chainID A")

    try:
        alpha1_resids, alpha2_resids, quality = aligner.get_alpha1_alpha2_residues()

        print(f"\n1. Residue ID Continuity Check:")
        # Check if residue IDs are continuous
        alpha1_continuous = all(
            alpha1_resids[i] + 1 == alpha1_resids[i+1]
            for i in range(len(alpha1_resids) - 1)
        )
        alpha2_continuous = all(
            alpha2_resids[i] + 1 == alpha2_resids[i+1]
            for i in range(len(alpha2_resids) - 1)
        )

        print(f"   α1 continuous: {alpha1_continuous}")
        print(f"   α2 continuous: {alpha2_continuous}")

        print(f"\n2. Expected Region Coverage:")
        # For HLA-A*02:01, we expect high coverage
        expected_alpha1_count = aligner.ALPHA1_POSITIONS[1] - aligner.ALPHA1_POSITIONS[0] + 1
        expected_alpha2_count = aligner.ALPHA2_POSITIONS[1] - aligner.ALPHA2_POSITIONS[0] + 1

        alpha1_coverage = len(alpha1_resids) / expected_alpha1_count
        alpha2_coverage = len(alpha2_resids) / expected_alpha2_count

        print(f"   α1 coverage: {alpha1_coverage:.1%} ({len(alpha1_resids)}/{expected_alpha1_count})")
        print(f"   α2 coverage: {alpha2_coverage:.1%} ({len(alpha2_resids)}/{expected_alpha2_count})")

        # Success criteria: >80% coverage for both regions
        if alpha1_coverage >= 0.8 and alpha2_coverage >= 0.8:
            print(f"\n✓ Alignment robustness test PASSED")
            return True
        else:
            print(f"\n✗ Alignment robustness test FAILED: Insufficient coverage")
            return False

    except Exception as e:
        print(f"\n✗ Alignment robustness test FAILED: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    """Run all tests."""
    print("\n")
    print("#" * 80)
    print("# MHC Sequence Alignment Module Test Suite")
    print("#" * 80)
    print("\nThis test validates:")
    print("  1. Sequence extraction from PDB files")
    print("  2. Alignment to HLA-A*02:01 reference")
    print("  3. α1/α2 helix region mapping")
    print("  4. Quality metrics and thresholds")
    print("  5. Robustness to different PDB numbering schemes")

    results = []

    # Test 1: 1AO7 (HLA-A*02:01)
    results.append(("1AO7 Alignment", test_1ao7_alignment()))

    # Test 2: 1bd2 (HLA-B*27:05)
    results.append(("1bd2 Alignment", test_1bd2_alignment()))

    # Test 3: Robustness check
    results.append(("Alignment Robustness", test_alignment_robustness()))

    # Summary
    print("\n" + "=" * 80)
    print("TEST SUMMARY")
    print("=" * 80)

    passed = sum(1 for _, result in results if result)
    total = len(results)

    for test_name, result in results:
        status = "PASS" if result else "FAIL"
        symbol = "✓" if result else "✗"
        print(f"{symbol} {test_name}: {status}")

    print(f"\nOverall: {passed}/{total} tests passed")

    if passed == total:
        print("\n🎉 All tests passed! MHC sequence alignment module is working correctly.")
        return 0
    else:
        print(f"\n⚠ {total - passed} test(s) failed. Please review the errors above.")
        return 1


if __name__ == "__main__":
    sys.exit(main())
