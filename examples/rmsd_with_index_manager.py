#!/usr/bin/env python3
"""
RMSD Calculation with IndexManager - Complete Example

This example demonstrates the complete workflow of using IndexManager
with RMSDCalculator for selective fit RMSD calculations.

Scenario: Calculate TCR RMSD while aligning on pHLA complex

Author: Immunex Development Team
Date: 2026-03-17
"""

import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from immunex.core import PipelineContext
from immunex.analysis import IndexManager, StandardComponent
from immunex.analysis.trajectory import RMSDCalculator
import numpy as np
import pandas as pd


def example_basic_rmsd_with_index_manager():
    """
    Example 1: Basic RMSD calculation with IndexManager.

    Align on pHLA, calculate TCR RMSD.
    """
    print("=" * 70)
    print("Example 1: Basic RMSD with IndexManager")
    print("=" * 70)

    # Step 1: Create PipelineContext
    context = PipelineContext(
        system_id="1ao7",
        topology="data/1ao7/md.tpr",
        trajectory_raw="data/1ao7/md_pbc.xtc",
        structure_pdb="data/1ao7/structure.pdb",
        output_dir="./output/1ao7"
    )

    print(f"\nTask: {context.system_id}")
    print(f"Topology: {context.topology}")
    print(f"Trajectory: {context.trajectory_raw}")

    # Step 2: Initialize RMSD calculator
    rmsd_calc = RMSDCalculator(
        topology=context.topology,
        trajectory=context.trajectory_raw
    )

    # Step 3: Use IndexManager for RMSD calculation
    print("\nCalculating RMSD with IndexManager...")
    print("  Alignment:    pHLA complex (chains A, B, C)")
    print("  Calculation:  TCR (chains D, E)")

    try:
        output_file = rmsd_calc.calculate_with_index_manager(
            context=context,
            fit_component='pHLA',
            calc_component='TCR',
            output_file=context.get_analysis_path('rmsd', 'tcr_rmsd.xvg')
        )

        print(f"\nRMSD calculation completed!")
        print(f"Output file: {output_file}")

        # Read and display statistics
        if Path(output_file).exists():
            # Parse XVG file (skip comments)
            data = np.loadtxt(output_file, comments=['#', '@'])
            times = data[:, 0]
            rmsd_values = data[:, 1]

            print(f"\nStatistics:")
            print(f"  Frames:      {len(rmsd_values)}")
            print(f"  Time range:  {times[0]:.1f} - {times[-1]:.1f} ps")
            print(f"  Mean RMSD:   {np.mean(rmsd_values):.3f} nm")
            print(f"  Std RMSD:    {np.std(rmsd_values):.3f} nm")
            print(f"  Min RMSD:    {np.min(rmsd_values):.3f} nm")
            print(f"  Max RMSD:    {np.max(rmsd_values):.3f} nm")

    except FileNotFoundError as e:
        print(f"\nNote: Example requires actual MD data files")
        print(f"Error: {e}")


def example_multiple_rmsd_calculations():
    """
    Example 2: Multiple RMSD calculations for different components.

    Demonstrates batch RMSD calculations for various scenarios:
    - pHLA stability (align & calc on pHLA)
    - TCR stability (align & calc on TCR)
    - Peptide stability (align on pHLA, calc on peptide)
    - TCR-peptide relative motion (align on pHLA, calc on TCR)
    """
    print("\n" + "=" * 70)
    print("Example 2: Multiple RMSD Calculations")
    print("=" * 70)

    context = PipelineContext(
        system_id="1ao7",
        topology="data/1ao7/md.tpr",
        trajectory_raw="data/1ao7/md_pbc.xtc",
        structure_pdb="data/1ao7/structure.pdb",
        output_dir="./output/1ao7"
    )

    rmsd_calc = RMSDCalculator(
        topology=context.topology,
        trajectory=context.trajectory_raw
    )

    # Define RMSD analysis scenarios
    scenarios = {
        'pHLA_stability': {
            'fit': 'pHLA',
            'calc': 'pHLA',
            'description': 'pHLA complex internal stability'
        },
        'TCR_stability': {
            'fit': 'TCR',
            'calc': 'TCR',
            'description': 'TCR complex internal stability'
        },
        'peptide_motion': {
            'fit': 'pHLA',
            'calc': 'peptide',
            'description': 'Peptide motion relative to pHLA'
        },
        'tcr_relative_motion': {
            'fit': 'pHLA',
            'calc': 'TCR',
            'description': 'TCR motion relative to pHLA'
        }
    }

    results = {}

    print("\nRunning multiple RMSD calculations...")
    print("-" * 70)

    for scenario_name, config in scenarios.items():
        print(f"\nScenario: {scenario_name}")
        print(f"  Description: {config['description']}")
        print(f"  Fit:  {config['fit']}")
        print(f"  Calc: {config['calc']}")

        try:
            output_file = rmsd_calc.calculate_with_index_manager(
                context=context,
                fit_component=config['fit'],
                calc_component=config['calc'],
                output_file=context.get_analysis_path('rmsd', f"{scenario_name}.xvg")
            )

            results[scenario_name] = output_file
            print(f"  Output: {output_file}")

        except FileNotFoundError as e:
            print(f"  Skipped (missing data files)")

    print("\n" + "-" * 70)
    print(f"Completed {len(results)} RMSD calculations")


def example_understanding_group_ids():
    """
    Example 3: Understanding group IDs and index file structure.

    Shows how IndexManager maps component names to group IDs.
    """
    print("\n" + "=" * 70)
    print("Example 3: Understanding Group IDs")
    print("=" * 70)

    context = PipelineContext(
        system_id="1ao7",
        topology="data/1ao7/md.tpr",
        trajectory_raw="data/1ao7/md_pbc.xtc",
        structure_pdb="data/1ao7/structure.pdb",
        output_dir="./output/1ao7"
    )

    # Access IndexManager
    index_mgr = context.index_manager

    print("\nStandardComponent Definitions:")
    print("-" * 70)
    for comp in StandardComponent:
        print(f"{comp.component_name:12s} | Chains: {', '.join(comp.chains)}")

    try:
        # Generate base index
        base_index = index_mgr.ensure_base_index()
        print(f"\nBase index file: {base_index}")

        # Show group registry
        print("\nGroup Registry (Component -> Group ID):")
        print("-" * 70)

        for comp_name in index_mgr.list_available_components():
            group_info = index_mgr.get_group_info(comp_name)
            if group_info:
                print(f"{comp_name:12s} | Group {group_info.group_id:2d} | "
                      f"{group_info.atom_count:5d} atoms")

        # Show how this translates to GROMACS stdin
        print("\nGROMACS stdin examples:")
        print("-" * 70)

        examples = [
            ('pHLA', 'TCR', 'Align pHLA, calc TCR RMSD'),
            ('TCR', 'peptide', 'Align TCR, calc peptide RMSD'),
            ('HLA_alpha', 'TCR_beta', 'Align HLA-alpha, calc TCR-beta RMSD')
        ]

        for fit_comp, calc_comp, description in examples:
            fit_id = index_mgr.get_group_id(fit_comp)
            calc_id = index_mgr.get_group_id(calc_comp)

            if fit_id and calc_id:
                stdin = f"{fit_id}\\n{calc_id}\\n"
                print(f"{description}")
                print(f"  Fit:  {fit_comp:12s} -> Group {fit_id}")
                print(f"  Calc: {calc_comp:12s} -> Group {calc_id}")
                print(f"  stdin: {repr(stdin)}")
                print()

    except FileNotFoundError as e:
        print(f"\nNote: Example requires actual MD data files")
        print(f"Error: {e}")


def example_comparison_old_vs_new():
    """
    Example 4: Code comparison - Old approach vs IndexManager.
    """
    print("\n" + "=" * 70)
    print("Example 4: Old vs New Approach Comparison")
    print("=" * 70)

    print("\n### OLD APPROACH (without IndexManager) ###\n")
    print("""
# 1. Manually generate index file
generator = IndexGenerator(topology="md.tpr")
generator.generate_multi_component_index(['pHLA', 'TCR'], 'phla_tcr.ndx')

# 2. Manually track group IDs (magic numbers!)
PHLA_GROUP = "20"  # Hard-coded, error-prone
TCR_GROUP = "21"   # Must match .ndx file order

# 3. Construct stdin manually
stdin_input = f"{PHLA_GROUP}\\n{TCR_GROUP}\\n"

# 4. Call GROMACS
subprocess.run(
    ["gmx", "rms", "-f", "md.xtc", "-s", "md.tpr",
     "-n", "phla_tcr.ndx", "-o", "rmsd.xvg"],
    input=stdin_input, text=True
)
""")

    print("\n### NEW APPROACH (with IndexManager) ###\n")
    print("""
# 1. Create context (automatic index management)
context = PipelineContext(
    system_id="1ao7",
    topology="md.tpr",
    trajectory_raw="md.xtc"
)

# 2. Calculate RMSD with one line (type-safe component names)
rmsd_calc = RMSDCalculator(topology="md.tpr", trajectory="md.xtc")
output = rmsd_calc.calculate_with_index_manager(
    context=context,
    fit_component='pHLA',   # No magic numbers!
    calc_component='TCR'    # Type-safe enums
)
""")

    print("\n### ADVANTAGES ###\n")
    print("✅ No manual index file generation")
    print("✅ No hard-coded group IDs")
    print("✅ Type-safe component names (StandardComponent enum)")
    print("✅ Automatic group registry building")
    print("✅ Reusable across all analysis modules")
    print("✅ Single source of truth for component definitions")


def main():
    """Run all examples."""
    print("\n" + "#" * 70)
    print("# RMSD with IndexManager - Complete Examples")
    print("#" * 70)

    try:
        example_basic_rmsd_with_index_manager()
        # example_multiple_rmsd_calculations()  # Requires MD data
        # example_understanding_group_ids()     # Requires MD data
        example_comparison_old_vs_new()

        print("\n" + "#" * 70)
        print("# Examples completed")
        print("#" * 70)

        print("\nNote: Some examples require actual MD data files:")
        print("  - md.tpr (topology)")
        print("  - md_pbc.xtc (PBC-corrected trajectory)")
        print("  - structure.pdb (standardized PDB)")
        print("\nTo run full examples, prepare these files and uncomment the")
        print("commented example function calls in main().")

    except Exception as e:
        print(f"\nExample error: {e}")
        print("\nThis is expected if MD data files are not available.")


if __name__ == "__main__":
    main()
