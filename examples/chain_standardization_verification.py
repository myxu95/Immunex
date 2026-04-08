#!/usr/bin/env python3
"""
Chain Standardization Verification Example

This example demonstrates how to verify that chains are properly standardized
before performing analysis with IndexManager.

Author: Immunex Development Team
Date: 2026-03-17
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from immunex.core import PipelineContext
from immunex.analysis import IndexManager


def example_chain_verification():
    """Example: Verify chain standardization status."""
    print("=" * 70)
    print("Chain Standardization Verification")
    print("=" * 70)

    # Create context
    context = PipelineContext(
        system_id="1ao7",
        topology="data/1ao7/md.tpr",
        trajectory_raw="data/1ao7/md_pbc.xtc",
        structure_pdb="data/1ao7/structure.pdb",
        output_dir="./output/1ao7"
    )

    # Access IndexManager (triggers chain validation)
    index_mgr = context.index_manager

    print("\n### Step 1: Check Chain Standardization Status ###\n")

    # Check if chains are standardized
    is_standardized = index_mgr.is_chain_standardized()
    print(f"Chains standardized: {is_standardized}")

    # Get detailed validation report
    report = index_mgr.get_chain_validation_report()

    print("\n### Step 2: Review Validation Report ###\n")
    print(f"Standardization performed: {report.get('standardized')}")
    print(f"Validation passed: {report.get('valid')}")

    if 'chains' in report:
        print("\nChain residue counts:")
        chains = report['chains']
        print(f"  Chain A (HLA-α):    {chains.get('A', 'N/A')} residues")
        print(f"  Chain B (β2m):      {chains.get('B', 'N/A')} residues")
        print(f"  Chain C (peptide):  {chains.get('C', 'N/A')} residues")
        print(f"  Chain D (TCR-α):    {chains.get('D', 'N/A')} residues")
        print(f"  Chain E (TCR-β):    {chains.get('E', 'N/A')} residues")

    if 'issues' in report:
        print("\n⚠️  Validation Issues:")
        for issue in report['issues']:
            print(f"  - {issue}")

    if 'reason' in report:
        print(f"\nReason: {report['reason']}")

    print("\n### Step 3: Decide How to Proceed ###\n")

    if is_standardized:
        print("✅ Chains are standardized and validated")
        print("   Safe to proceed with analysis")

        # Generate base index
        base_index = index_mgr.ensure_base_index()
        print(f"\nBase index generated: {base_index}")

    elif report.get('standardized') is False:
        print("❌ Chain standardization FAILED")
        print("   Analysis results will be INCORRECT!")
        print("\nRecommendations:")
        print("  1. Verify input PDB file has correct chain ordering")
        print("  2. Manually standardize chains using PDBChainStandardizer")
        print("  3. Provide pre-standardized PDB via context.structure_pdb")

    else:
        print("⚠️  Chain validation has warnings")
        print("   Proceed with caution and verify results manually")


def example_handling_unstandardized_chains():
    """Example: How to handle unstandardized chains."""
    print("\n" + "=" * 70)
    print("Handling Unstandardized Chains")
    print("=" * 70)

    context = PipelineContext(
        system_id="unknown_system",
        topology="data/unknown/md.tpr",
        trajectory_raw="data/unknown/md_pbc.xtc",
        output_dir="./output/unknown"
    )

    index_mgr = context.index_manager

    # Check standardization
    if not index_mgr.is_chain_standardized():
        print("\n⚠️  Chains are not standardized!")

        report = index_mgr.get_chain_validation_report()

        # Option 1: Abort analysis
        if report.get('standardized') is False:
            print("\nOption 1: Abort analysis (RECOMMENDED)")
            print("  Reason: Chain mapping is incorrect")
            print("  Action: Fix input data and retry")
            return

        # Option 2: Proceed with manual verification
        print("\nOption 2: Proceed with manual verification")
        print("  Action: Manually verify chain assignments")
        print("  Warning: Results may be incorrect!")

        # Check errors in context
        if context.has_errors():
            print("\nContext errors:")
            for error in context.errors:
                print(f"  - {error}")


def example_chain_standardization_in_pipeline():
    """Example: Chain verification integrated in pipeline."""
    print("\n" + "=" * 70)
    print("Chain Verification in Analysis Pipeline")
    print("=" * 70)

    context = PipelineContext(
        system_id="1ao7",
        topology="data/1ao7/md.tpr",
        trajectory_raw="data/1ao7/md_pbc.xtc",
        structure_pdb="data/1ao7/structure.pdb"
    )

    print("\n### Pipeline Step 1: Pre-flight Check ###\n")

    # Pre-flight check: Verify chains before starting analysis
    index_mgr = context.index_manager

    if not index_mgr.is_chain_standardized():
        print("❌ Pre-flight check FAILED")
        print("   Aborting analysis")

        # Log errors
        report = index_mgr.get_chain_validation_report()
        print(f"\nFailure reason: {report.get('reason', 'Unknown')}")

        # Add to context errors (already done by ensure_base_index)
        if context.has_errors():
            print("\nErrors logged to context:")
            for error in context.errors:
                print(f"  - {error}")

        return

    print("✅ Pre-flight check PASSED")

    print("\n### Pipeline Step 2: Generate Indices ###\n")

    # Generate base index (chain validation already done)
    base_index = index_mgr.ensure_base_index()
    print(f"Base index: {base_index}")

    # Get group IDs
    phla_id = index_mgr.get_group_id('pHLA')
    tcr_id = index_mgr.get_group_id('TCR')

    print(f"\nGroup IDs:")
    print(f"  pHLA: {phla_id}")
    print(f"  TCR:  {tcr_id}")

    print("\n### Pipeline Step 3: Perform Analysis ###\n")
    print("(RMSD calculation, RMSF analysis, etc.)")
    print("✅ Analysis completed successfully")


def main():
    """Run all examples."""
    print("\n" + "#" * 70)
    print("# Chain Standardization Verification Examples")
    print("#" * 70)

    try:
        example_chain_verification()
        # example_handling_unstandardized_chains()  # Requires data
        # example_chain_standardization_in_pipeline()  # Requires data

        print("\n" + "#" * 70)
        print("# Examples completed")
        print("#" * 70)

        print("\n📖 Summary:")
        print("  1. Always check index_mgr.is_chain_standardized() before analysis")
        print("  2. Review index_mgr.get_chain_validation_report() for details")
        print("  3. Abort analysis if standardization failed")
        print("  4. Chain validation is automatically performed on first access")

    except FileNotFoundError as e:
        print(f"\nNote: Examples require MD data files")
        print(f"Error: {e}")


if __name__ == "__main__":
    main()
