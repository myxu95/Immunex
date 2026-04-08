#!/usr/bin/env python3
"""
Group Selector usage examples for Immunex toolkit.

This script demonstrates how to use the GroupSelector for intelligent
GROMACS group selection in MD analysis workflows.
"""

import os
from immunex import RMSDCalculator, PBCProcessor
from immunex.utils import GroupSelector


def example_basic_group_selection():
    """Basic group selection examples."""
    print("=== Basic Group Selection ===")
    
    topology = "protein.tpr"
    
    if not os.path.exists(topology):
        print("Topology file not found - using placeholder examples")
        return
    
    # Initialize group selector
    selector = GroupSelector(topology)
    
    # List available groups
    print("Available groups:")
    groups = selector.list_available_groups()
    for group_id, info in groups.items():
        print(f"  {group_id}: {info['name']} ({info['atoms']} atoms)")
    
    # Auto-select groups for different purposes
    purposes = [
        "rmsd_backbone", "rmsd_protein", "rmsd_calpha",
        "center_protein", "output_system", "pbc_protein"
    ]
    
    print("\nAuto-selected groups:")
    for purpose in purposes:
        group = selector.select_group(purpose)
        print(f"  {purpose}: {group}")


def example_rmsd_with_group_selection():
    """RMSD calculation with intelligent group selection."""
    print("\n=== RMSD with Group Selection ===")
    
    topology = "protein.tpr"
    trajectory = "trajectory.xtc"
    
    if not (os.path.exists(topology) and os.path.exists(trajectory)):
        print("Files not found - using placeholder examples")
        return
    
    # Initialize RMSD calculator (automatically includes GroupSelector)
    rmsd_calc = RMSDCalculator(topology, trajectory)
    
    # Method 1: Use predefined RMSD types (auto-selects groups)
    print("1. Auto-selection based on RMSD type:")
    
    # Backbone RMSD - automatically selects appropriate backbone groups
    rmsd_file = rmsd_calc.calculate_gromacs(
        rmsd_type="backbone",
        output_file="rmsd_backbone.xvg"
    )
    print(f"  Backbone RMSD: {rmsd_file}")
    
    # Protein RMSD - automatically selects protein groups
    rmsd_file = rmsd_calc.calculate_gromacs(
        rmsd_type="protein", 
        output_file="rmsd_protein.xvg"
    )
    print(f"  Protein RMSD: {rmsd_file}")
    
    # C-alpha RMSD - automatically selects C-alpha groups
    rmsd_file = rmsd_calc.calculate_gromacs(
        rmsd_type="calpha",
        output_file="rmsd_calpha.xvg"
    )
    print(f"  C-alpha RMSD: {rmsd_file}")
    
    # Method 2: Custom group selection
    print("\n2. Custom group selection:")
    rmsd_file = rmsd_calc.calculate_gromacs_custom_groups(
        reference_group="Backbone",
        analysis_group="Protein",
        output_file="rmsd_custom.xvg"
    )
    print(f"  Custom RMSD: {rmsd_file}")


def example_pbc_processing_with_group_selection():
    """PBC processing with intelligent group selection."""
    print("\n=== PBC Processing with Group Selection ===")
    
    trajectory = "raw_trajectory.xtc"
    topology = "system.tpr"
    
    if not (os.path.exists(topology) and os.path.exists(trajectory)):
        print("Files not found - using placeholder examples")
        return
    
    # Initialize PBC processor
    pbc_processor = PBCProcessor()
    
    # Method 1: Auto-select groups (recommended)
    print("1. Auto-selection of groups:")
    processed_traj = pbc_processor.remove_pbc(
        trajectory=trajectory,
        topology=topology,
        output="processed_auto.xtc"
        # center_group and output_group will be auto-selected
    )
    print(f"  Processed trajectory: {processed_traj}")
    
    # Method 2: Manual group specification
    print("\n2. Manual group specification:")
    processed_traj = pbc_processor.remove_pbc(
        trajectory=trajectory,
        topology=topology,
        output="processed_manual.xtc",
        center_group="Protein",
        output_group="System"
    )
    print(f"  Processed trajectory: {processed_traj}")


def example_custom_group_mappings():
    """Creating and using custom group mappings."""
    print("\n=== Custom Group Mappings ===")
    
    topology = "protein.tpr"
    
    if not os.path.exists(topology):
        print("Topology file not found - using placeholder examples")
        return
    
    selector = GroupSelector(topology)
    
    # Add custom mappings for specific project needs
    print("Adding custom mappings:")
    
    # Custom mapping for membrane protein analysis
    selector.add_custom_mapping("membrane_protein", ["Protein_Membrane", "Protein", "1"])
    selector.add_custom_mapping("lipid_bilayer", ["POPC", "Lipids", "15"])
    selector.add_custom_mapping("water_ions", ["Water_Ion", "non-Protein", "11"])
    
    # Custom mapping for drug-protein interaction
    selector.add_custom_mapping("drug_molecule", ["LIG", "MOL", "Drug"])
    selector.add_custom_mapping("binding_site", ["r_100-120", "Binding_Site"])
    
    print("  Added custom mappings for membrane and drug analysis")
    
    # Use custom mappings
    print("\nUsing custom mappings:")
    purposes = ["membrane_protein", "lipid_bilayer", "drug_molecule"]
    for purpose in purposes:
        group = selector.select_group(purpose)
        print(f"  {purpose}: {group}")


def example_group_suggestions():
    """Getting group suggestions for specific purposes."""
    print("\n=== Group Suggestions ===")
    
    topology = "protein.tpr"
    
    if not os.path.exists(topology):
        print("Topology file not found - using placeholder examples")
        return
    
    selector = GroupSelector(topology)
    
    # Get suggestions for different analysis types
    analysis_types = [
        "protein_analysis",
        "backbone_fitting", 
        "water_analysis",
        "heavy_atoms"
    ]
    
    print("Group suggestions for different analysis types:")
    for analysis_type in analysis_types:
        suggestions = selector.suggest_groups_for_purpose(analysis_type)
        print(f"\n  {analysis_type}:")
        for group_name, reason in suggestions:
            print(f"    {group_name} ({reason})")


def example_selection_presets():
    """Creating and using selection presets for workflows."""
    print("\n=== Selection Presets ===")
    
    topology = "protein.tpr"
    
    if not os.path.exists(topology):
        print("Topology file not found - using placeholder examples")
        return
    
    selector = GroupSelector(topology)
    
    # Create preset for common protein analysis workflow
    protein_analysis_preset = {
        "rmsd_reference": "Backbone",
        "rmsd_analysis": "Protein", 
        "center_group": "Protein",
        "output_group": "System",
        "rdf_protein": "Protein",
        "rdf_water": "SOL"
    }
    
    selector.create_selection_preset("protein_analysis", protein_analysis_preset)
    print("Created 'protein_analysis' preset")
    
    # Create preset for membrane protein analysis
    membrane_analysis_preset = {
        "protein_group": "Protein",
        "membrane_group": "POPC",
        "water_group": "SOL",
        "center_group": "Protein_Membrane",
        "output_group": "System"
    }
    
    selector.create_selection_preset("membrane_analysis", membrane_analysis_preset)
    print("Created 'membrane_analysis' preset")
    
    # Use preset
    print("\nUsing protein_analysis preset:")
    preset_selections = selector.use_preset("protein_analysis")
    for purpose, group in preset_selections.items():
        print(f"  {purpose}: {group}")


def example_integration_workflow():
    """Complete workflow example with group selection integration."""
    print("\n=== Complete Workflow Example ===")
    
    topology = "protein.tpr"
    trajectory = "trajectory.xtc"
    
    if not (os.path.exists(topology) and os.path.exists(trajectory)):
        print("Files not found - simulating workflow")
        
        # Simulate the workflow steps
        print("Workflow steps (simulated):")
        print("1. Initialize GroupSelector with topology")
        print("2. PBC processing with auto-selected groups")
        print("3. RMSD analysis with different group types") 
        print("4. Multiple analysis types with consistent group selection")
        return
    
    print("Running complete analysis workflow:")
    
    # Step 1: PBC processing
    print("1. PBC Processing...")
    pbc_processor = PBCProcessor()
    processed_traj = pbc_processor.remove_pbc(
        trajectory=trajectory,
        topology=topology,
        output="workflow_processed.xtc"
    )
    
    # Step 2: Multiple RMSD analyses
    print("2. RMSD Analyses...")
    rmsd_calc = RMSDCalculator(topology, processed_traj)
    
    # Different RMSD types using auto-selection
    rmsd_types = ["backbone", "protein", "calpha"]
    rmsd_files = []
    
    for rmsd_type in rmsd_types:
        output_file = f"workflow_rmsd_{rmsd_type}.xvg"
        rmsd_file = rmsd_calc.calculate_gromacs(
            rmsd_type=rmsd_type,
            output_file=output_file
        )
        rmsd_files.append(rmsd_file)
        print(f"   - {rmsd_type.capitalize()} RMSD: {rmsd_file}")
    
    print(f"\nWorkflow completed! Generated files:")
    print(f"  - Processed trajectory: {processed_traj}")
    for rmsd_file in rmsd_files:
        print(f"  - RMSD analysis: {rmsd_file}")


def main():
    """Run all group selector examples."""
    print("Immunex Group Selector - Usage Examples")
    print("=" * 50)
    
    example_basic_group_selection()
    example_rmsd_with_group_selection()
    example_pbc_processing_with_group_selection()
    example_custom_group_mappings()
    example_group_suggestions()
    example_selection_presets()
    example_integration_workflow()
    
    print("\n" + "=" * 50)
    print("Group Selector examples completed!")
    print("\nKey Benefits:")
    print("- Automatic group selection based on analysis type")
    print("- Consistent group choices across analysis workflows")
    print("- Customizable mappings for specific project needs")
    print("- Intelligent fallbacks when preferred groups unavailable")
    print("- User configuration persistence across sessions")


if __name__ == "__main__":
    main()