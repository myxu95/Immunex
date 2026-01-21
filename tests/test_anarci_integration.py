#!/usr/bin/env python3
"""
Test ANARCI and Intelligent Chain Identification System
"""

import sys
from pathlib import Path

# Add aftermd to path
sys.path.insert(0, str(Path(__file__).parent.parent))

def test_anarci_basic():
    """Test basic ANARCI functionality."""
    print("="*80)
    print("Test 1: Basic ANARCI Functionality")
    print("="*80)

    try:
        from anarci import anarci

        # Real TCR beta sequence from PDB 1AO7
        tcr_beta_seq = (
            "NAGVTQTPKFQVLKTGQSMTLQCAQDMNHEYMSWYRQDPGMGLRLIHYSVGAGITDQGEVP"
            "NGYNVSRSTTEDFPLRLLSAAPSQTSVYFCASSYVGNTGELFFGEGSRLTVL"
        )

        print(f"\nTesting with TCR beta sequence ({len(tcr_beta_seq)} aa)...")
        print(f"Sequence: {tcr_beta_seq[:50]}...")

        results = anarci(
            [('tcr_beta', tcr_beta_seq)],
            scheme='imgt',
            output=False
        )

        if results and results[0] and results[0][0]:
            numbering, alignment, chain_type = results[0][0]
            print(f"\n✓ SUCCESS!")
            print(f"  Chain type: {chain_type}")
            print(f"  Alignment found: Yes")
            return True
        else:
            print("\n⚠ WARNING: No alignment found")
            print("  This may be normal if the sequence doesn't match TCR/Ig patterns")
            return False

    except Exception as e:
        print(f"\n✗ FAILED: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_anarci_wrapper():
    """Test AfterMD's ANARCIWrapper."""
    print("\n" + "="*80)
    print("Test 2: AfterMD ANARCIWrapper")
    print("="*80)

    try:
        from aftermd.utils.cdr_manager import ANARCIWrapper

        wrapper = ANARCIWrapper(allow_fallback=True)
        print(f"\n✓ ANARCIWrapper initialized")
        print(f"  ANARCI available: {wrapper.anarci_available}")

        if wrapper.anarci_available:
            # Test with TCR sequence
            tcr_seq = (
                "NAGVTQTPKFQVLKTGQSMTLQCAQDMNHEYMSWYRQDPGMGLRLIHYSVGAGITDQGEVP"
                "NGYNVSRSTTEDFPLRLLSAAPSQTSVYFCASSYVGNTGELFFGEGSRLTVL"
            )

            result = wrapper.run_anarci(tcr_seq, chain_type='TCR')
            print(f"\n  Result keys: {list(result.keys())}")

            if 'chain_type' in result:
                print(f"  Chain type identified: {result['chain_type']}")

            return True
        else:
            print("\n⚠ ANARCI not available, using fallback")
            return False

    except Exception as e:
        print(f"\n✗ FAILED: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_sequence_extractor():
    """Test PDBSequenceExtractor."""
    print("\n" + "="*80)
    print("Test 3: PDBSequenceExtractor")
    print("="*80)

    try:
        from aftermd.utils import PDBSequenceExtractor

        extractor = PDBSequenceExtractor()
        print(f"\n✓ PDBSequenceExtractor initialized")
        print(f"  Supported residues: {len(extractor.three_to_one)} types")

        # Test three-to-one conversion
        test_residues = ['ALA', 'GLY', 'VAL', 'MSE', 'HSE']
        print(f"\n  Testing residue conversion:")
        for res in test_residues:
            one_letter = extractor.three_to_one.get(res, 'X')
            print(f"    {res} → {one_letter}")

        return True

    except Exception as e:
        print(f"\n✗ FAILED: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_intelligent_identifier():
    """Test IntelligentChainIdentifier."""
    print("\n" + "="*80)
    print("Test 4: IntelligentChainIdentifier")
    print("="*80)

    try:
        from aftermd.utils import IntelligentChainIdentifier

        identifier = IntelligentChainIdentifier(use_anarci=True)
        print(f"\n✓ IntelligentChainIdentifier initialized")
        print(f"  ANARCI mode: {identifier.use_anarci}")

        if identifier.anarci_wrapper:
            print(f"  ANARCIWrapper available: {identifier.anarci_wrapper.anarci_available}")

        # Test length thresholds
        print(f"\n  Length thresholds:")
        print(f"    Peptide max: {identifier.PEPTIDE_MAX_LENGTH} aa")
        print(f"    Beta2m range: {identifier.BETA2M_MIN_LENGTH}-{identifier.BETA2M_MAX_LENGTH} aa")
        print(f"    HLA-alpha min: {identifier.HLA_ALPHA_MIN_LENGTH} aa")

        return True

    except Exception as e:
        print(f"\n✗ FAILED: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_pdb_chain_standardizer():
    """Test updated PDBChainStandardizer."""
    print("\n" + "="*80)
    print("Test 5: Updated PDBChainStandardizer")
    print("="*80)

    try:
        from aftermd.utils import PDBChainStandardizer

        # Test legacy mode
        std_legacy = PDBChainStandardizer(use_intelligent_identification=False)
        print(f"\n✓ Legacy mode initialized")
        print(f"  Intelligent mode: {std_legacy.use_intelligent_identification}")
        print(f"  Standard order: {std_legacy.standard_order}")

        # Test intelligent mode
        std_intelligent = PDBChainStandardizer(use_intelligent_identification=True)
        print(f"\n✓ Intelligent mode initialized")
        print(f"  Intelligent mode: {std_intelligent.use_intelligent_identification}")

        if std_intelligent.intelligent_identifier:
            print(f"  IntelligentChainIdentifier: Available")
        else:
            print(f"  IntelligentChainIdentifier: Not available (fallback to legacy)")

        return True

    except Exception as e:
        print(f"\n✗ FAILED: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    """Run all tests."""
    print("\n" + "="*80)
    print("AfterMD Intelligent Chain Identification - Test Suite")
    print("="*80)

    results = {
        'Basic ANARCI': test_anarci_basic(),
        'ANARCIWrapper': test_anarci_wrapper(),
        'SequenceExtractor': test_sequence_extractor(),
        'IntelligentIdentifier': test_intelligent_identifier(),
        'PDBChainStandardizer': test_pdb_chain_standardizer()
    }

    # Summary
    print("\n" + "="*80)
    print("Test Summary")
    print("="*80)

    for test_name, passed in results.items():
        status = "✓ PASS" if passed else "✗ FAIL"
        print(f"  {test_name:30s} {status}")

    total = len(results)
    passed = sum(results.values())
    print(f"\nTotal: {passed}/{total} tests passed")

    if passed == total:
        print("\n🎉 All tests passed!")
        return 0
    else:
        print(f"\n⚠ {total - passed} test(s) failed")
        return 1


if __name__ == "__main__":
    sys.exit(main())
