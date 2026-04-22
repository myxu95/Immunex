#!/usr/bin/env python3
"""
Index Manager Usage Example

This example demonstrates how to use the unified IndexManager for
GROMACS index file management and RMSD calculations.

Author: Immunex Development Team
Date: 2026-03-17
"""

import sys
from pathlib import Path

# Add parent directory to path for local imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from immunex.core import PipelineContext
from immunex.analysis import IndexManager, StandardComponent


def example_1_basic_usage():
    """Example 1: Basic IndexManager usage."""
    print("=" * 60)
    print("Example 1: Basic IndexManager Usage")
    print("=" * 60)

    # Create a context
    context = PipelineContext(
        system_id="1ao7",
        topology="data/1ao7/md.tpr",
        trajectory_raw="data/1ao7/md.xtc",
        structure_pdb="data/1ao7/structure.pdb",
        output_dir="./output/1ao7"
    )

    # Access IndexManager (lazily initialized)
    index_mgr = context.index_manager

    print(f"\nIndexManager initialized: {index_mgr}")

    # Ensure base index exists
    base_index = index_mgr.ensure_base_index()
    print(f"\nBase index file: {base_index}")

    # List available components
    components = index_mgr.list_available_components()
    print(f"\nAvailable components: {components}")

    # Get group IDs for components
    print("\nComponent group IDs:")
    for comp_name in ['pHLA', 'TCR', 'peptide', 'HLA_alpha', 'TCR_beta']:
        group_id = index_mgr.get_group_id(comp_name)
        group_info = index_mgr.get_group_info(comp_name)
        if group_info:
            print(f"  {comp_name:12s} -> Group {group_id:2s} ({group_info.atom_count} atoms)")


def example_2_standard_components():
    """Example 2: Using StandardComponent Enum."""
    print("\n" + "=" * 60)
    print("Example 2: StandardComponent Enum")
    print("=" * 60)

    print("\nStandardComponent definitions:")
    for comp in StandardComponent:
        print(f"  {comp.name:12s} -> chains: {comp.chains}, name: {comp.component_name}")


def example_3_rmsd_with_index_manager():
    """Example 3: RMSD calculation using IndexManager."""
    print("\n" + "=" * 60)
    print("Example 3: RMSD with IndexManager")
    print("=" * 60)

    # Create context
    context = PipelineContext(
        system_id="1ao7",
        topology="data/1ao7/md.tpr",
        trajectory_raw="data/1ao7/md_pbc.xtc",
        structure_pdb="data/1ao7/structure.pdb",
        output_dir="./output/1ao7"
    )

    # Get group IDs from IndexManager
    index_mgr = context.index_manager
    phla_group_id = index_mgr.get_group_id('pHLA')
    tcr_group_id = index_mgr.get_group_id('TCR')

    print(f"\nGroup IDs for RMSD calculation:")
    print(f"  pHLA (fit):  Group {phla_group_id}")
    print(f"  TCR (calc):  Group {tcr_group_id}")

    # Construct GROMACS stdin input
    stdin_input = f"{phla_group_id}\n{tcr_group_id}\n"

    print(f"\nGROMACS stdin input:")
    print(f"  {repr(stdin_input)}")

    print("\nThis can be used directly with gmx rms:")
    print(f"  gmx rms -f md_pbc.xtc -s md.tpr -n {index_mgr.base_index_file.name}")
    print(f"  stdin: {stdin_input.strip()}")


def example_4_cdr3_dynamic_index():
    """Example 4: CDR3 dynamic index generation."""
    print("\n" + "=" * 60)
    print("Example 4: CDR3 Dynamic Index")
    print("=" * 60)

    # Create context
    context = PipelineContext(
        system_id="1ao7",
        topology="data/1ao7/md.tpr",
        trajectory_raw="data/1ao7/md_pbc.xtc",
        structure_pdb="data/1ao7/structure.pdb",
        output_dir="./output/1ao7"
    )

    index_mgr = context.index_manager

    # Generate CDR3 index (example sequence)
    cdr3_beta_seq = "CASSLGQAYEQYF"

    print(f"\nGenerating CDR3 beta index for sequence: {cdr3_beta_seq}")

    cdr3_group_id = index_mgr.ensure_cdr3_index(
        chain='beta',
        sequence=cdr3_beta_seq,
        ca_only=True
    )

    if cdr3_group_id:
        print(f"CDR3 beta group ID: {cdr3_group_id}")

        # Get combined index with CDR3
        combined_index = index_mgr.get_combined_index('pHLA', 'TCR', 'CDR3_TCR_beta')
        print(f"\nCombined index file: {combined_index}")
    else:
        print("CDR3 sequence not found in structure")


def example_5_multiple_components():
    """Example 5: Working with multiple components."""
    print("\n" + "=" * 60)
    print("Example 5: Multiple Components")
    print("=" * 60)

    context = PipelineContext(
        system_id="1ao7",
        topology="data/1ao7/md.tpr",
        trajectory_raw="data/1ao7/md_pbc.xtc",
        structure_pdb="data/1ao7/structure.pdb",
        output_dir="./output/1ao7"
    )

    index_mgr = context.index_manager

    # Get group IDs for multiple analysis scenarios
    scenarios = {
        'pHLA stability': ['pHLA'],
        'TCR stability': ['TCR'],
        'TCR-pHLA interface': ['pHLA', 'TCR'],
        'Peptide-TCR contact': ['peptide', 'TCR'],
        'HLA alpha-TCR beta': ['HLA_alpha', 'TCR_beta']
    }

    print("\nAnalysis scenarios and their group IDs:")
    for scenario_name, components in scenarios.items():
        group_ids = [index_mgr.get_group_id(comp) for comp in components]
        print(f"\n  {scenario_name}:")
        for comp, gid in zip(components, group_ids):
            print(f"    {comp:12s} -> Group {gid}")


def main():
    """Run all examples."""
    print("\n" + "#" * 60)
    print("# IndexManager Usage Examples")
    print("#" * 60)

    try:
        example_1_basic_usage()
        example_2_standard_components()
        example_3_rmsd_with_index_manager()
        # example_4_cdr3_dynamic_index()  # Requires actual data
        example_5_multiple_components()

        print("\n" + "#" * 60)
        print("# All examples completed")
        print("#" * 60)
        print("\nNote: Examples 1, 3, 4, 5 require actual MD data files.")
        print("To run with real data:")
        print("  1. Prepare MD files: md.tpr, md_pbc.xtc, structure.pdb")
        print("  2. Update paths in the examples")
        print("  3. Run this script")

    except FileNotFoundError as e:
        print(f"\nNote: Skipping examples that require actual MD files")
        print(f"Error: {e}")


if __name__ == "__main__":
    main()
