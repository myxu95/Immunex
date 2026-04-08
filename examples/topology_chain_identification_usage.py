#!/usr/bin/env python3
"""
Topology-Based Chain Identification Usage Examples

This script demonstrates how to use topology-based chain identification
to avoid PDB-TPR chain ID mismatch issues.

Problem:
--------
When PDB structures are standardized (chains A/B/C/D/E in specific order),
but TPR/trajectory files maintain original chain IDs, index files generated
from standardized PDB will be incorrect when applied to TPR files.

Solution:
---------
Work directly with topology files (TPR/GRO) to identify chains based on
atom/residue counts, independent of chain ID labels.

Author: Immunex Development Team
Date: 2026-03-18
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from immunex.core import PipelineContext
from immunex.analysis.topology import TopologyChainIdentifier


def example1_direct_topology_identification():
    """
    Example 1: Direct topology-based chain identification

    Use TopologyChainIdentifier directly to identify chains from TPR file.
    """
    print("=" * 70)
    print("Example 1: Direct Topology Identification")
    print("=" * 70)

    # Initialize identifier
    identifier = TopologyChainIdentifier()

    # Identify chains from topology file
    result = identifier.identify_chains_from_topology(
        topology_file="data/1ao7/md.tpr",
        output_index_file="data/1ao7/topology_components.ndx"
    )

    if result.success:
        print("\nChain Identification Successful!")
        print("\nComponent Mapping:")
        for comp_name, group_id in sorted(result.component_map.items()):
            chain = result.chains[group_id]
            print(f"  {comp_name:12s} -> Group {group_id:2d} "
                  f"({chain.residue_count} residues, {chain.atom_count} atoms)")

        print("\nGenerated index file: topology_components.ndx")
        print("\nYou can now use group IDs in GROMACS commands:")
        print(f"  gmx trjconv -f md.xtc -s md.tpr -n topology_components.ndx")
    else:
        print(f"\nChain identification failed: {result.error_message}")


def example2_context_based_usage():
    """
    Example 2: Using topology identification through PipelineContext

    Recommended approach for integration with Immunex workflows.
    """
    print("\n" + "=" * 70)
    print("Example 2: Context-Based Usage (Recommended)")
    print("=" * 70)

    # Create pipeline context
    context = PipelineContext(
        system_id="1ao7",
        topology="data/1ao7/md.tpr",
        trajectory_raw="data/1ao7/md.xtc",
        output_dir="results/1ao7"
    )

    # Get index manager
    index_mgr = context.index_manager

    # Use topology-based method (instead of PDB-based)
    print("\nUsing topology-based chain identification...")
    result = index_mgr.identify_chains_from_topology()

    if result.success:
        print("\nChain Identification Successful!")
        print("\nComponent mapping stored in context metadata:")
        mapping = context.metadata.get('topology_chain_mapping', {})
        for comp_name, group_id in sorted(mapping.items()):
            print(f"  {comp_name:12s} -> Group {group_id}")

        # Generate base index using topology method
        print("\nGenerating base index from topology...")
        base_index = index_mgr.ensure_base_index_from_topology()
        print(f"Base index created: {base_index}")

        # Get group IDs for analysis
        phla_id = index_mgr.get_group_id('pHLA')
        tcr_id = index_mgr.get_group_id('TCR')
        peptide_id = index_mgr.get_group_id('peptide')

        if phla_id is not None:
            print(f"\nReady for analysis:")
            print(f"  pHLA group ID: {phla_id}")
            print(f"  TCR group ID: {tcr_id}")
            print(f"  Peptide group ID: {peptide_id}")
    else:
        print(f"\nFailed: {result.error_message}")


def main():
    """Run all examples."""
    print("\n")
    print("=" * 70)
    print("  Topology-Based Chain Identification Examples")
    print("=" * 70)

    print("\nNote: Examples use dummy file paths")
    print("Replace paths with your actual topology files to run.\n")

    print("=" * 70)
    print("Quick Reference")
    print("=" * 70)
    print("""
Basic Usage:
------------
from immunex.analysis.topology import TopologyChainIdentifier

identifier = TopologyChainIdentifier()
result = identifier.identify_chains_from_topology(
    topology_file="md.tpr",
    output_index_file="components.ndx"
)

Integration with PipelineContext:
---------------------------------
from immunex.core import PipelineContext

context = PipelineContext(system_id="1ao7", topology="md.tpr", ...)
index_mgr = context.index_manager
index_mgr.ensure_base_index_from_topology()

When to Use:
------------
Use topology-based when:
  - PDB and TPR have different chain orderings
  - Chain standardization warnings appear
  - You want deterministic, reproducible chain identification
    """)


if __name__ == "__main__":
    main()
