#!/usr/bin/env python3
"""
Modular usage examples for Immunex toolkit.

This script demonstrates how to use the new modular structure 
with individual functional modules.
"""

import os
import numpy as np
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor

# Import individual modules for specific functionality
from immunex import (
    # Utils
    PlotManager, PathManager,
    # Trajectory analysis modules
    RMSDCalculator, RDFCalculator, RadiusGyrationCalculator,
    DistanceCalculator, HydrogenBondAnalyzer,
    # Structure analysis modules  
    BFactorAnalyzer, ContactMapCalculator,
    GeometryAnalyzer, AtomInfoExtractor
)
from immunex.analysis import PBCProcessor


def example_rmsd_analysis():
    """Example of using the dedicated RMSD module."""
    print("=== RMSD Analysis Module ===")
    
    topology = "protein.tpr"
    trajectory = "trajectory.xtc"
    
    if os.path.exists(topology) and os.path.exists(trajectory):
        # Initialize RMSD calculator
        rmsd_calc = RMSDCalculator(topology, trajectory)
        
        # Calculate RMSD using MDAnalysis
        times, rmsd_values = rmsd_calc.calculate_mdanalysis(
            selection="protein and name CA",
            output_file="rmsd_mdanalysis.csv"
        )
        
        # Calculate RMSD using GROMACS
        rmsd_file = rmsd_calc.calculate_gromacs(
            selection_group="Backbone",
            output_file="rmsd_gromacs.xvg"
        )
        
        print(f"RMSD analysis completed. Mean RMSD: {np.mean(rmsd_values):.3f} nm")
        
    else:
        print("Trajectory files not found - using placeholder")


def example_structure_analysis():
    """Example of using structure analysis modules."""
    print("\n=== Structure Analysis Modules ===")
    
    structure_file = "protein.pdb"
    
    if os.path.exists(structure_file):
        # B-factor analysis
        bfactor_analyzer = BFactorAnalyzer(structure_file)
        indices, bfactors = bfactor_analyzer.extract_bfactors(
            selection="protein",
            output_file="bfactors.csv"
        )
        
        residue_bfactors = bfactor_analyzer.analyze_by_residue(
            output_file="bfactor_by_residue.csv"
        )
        
        # Contact map calculation
        contact_calc = ContactMapCalculator(structure_file)
        contact_map = contact_calc.calculate(
            selection="name CA",
            cutoff=8.0,
            output_file="contact_map.csv"
        )
        
        # Geometry analysis
        geom_analyzer = GeometryAnalyzer(structure_file)
        geometry_props = geom_analyzer.analyze_basic_properties()
        print(f"Radius of gyration: {geometry_props['radius_of_gyration']:.3f} nm")
        
        # Atom information extraction
        atom_extractor = AtomInfoExtractor(structure_file)
        atom_info = atom_extractor.extract_comprehensive_info(
            output_file="atom_info.csv"
        )
        
        structure_summary = atom_extractor.get_structure_summary()
        print(f"Structure contains {structure_summary['n_atoms']} atoms")
        
    else:
        print("Structure file not found - using placeholder")


def example_trajectory_modules():
    """Example of using individual trajectory analysis modules."""
    print("\n=== Individual Trajectory Modules ===")
    
    topology = "protein.tpr"
    trajectory = "trajectory.xtc"
    
    if os.path.exists(topology) and os.path.exists(trajectory):
        # RDF calculation
        rdf_calc = RDFCalculator(topology, trajectory)
        bins, rdf_values = rdf_calc.calculate(
            selection1="name O and resname SOL",
            selection2="name O and resname SOL",
            output_file="water_rdf.csv"
        )
        
        # Radius of gyration
        rg_calc = RadiusGyrationCalculator(topology, trajectory)
        times, rg_values = rg_calc.calculate(
            selection="protein",
            output_file="radius_gyration.csv"
        )
        
        # Components analysis
        times, rgx, rgy, rgz = rg_calc.calculate_components(
            output_file="rg_components.csv"
        )
        
        # Distance calculations
        dist_calc = DistanceCalculator(topology, trajectory)
        times, distances = dist_calc.calculate_center_of_mass_distance(
            selection1="protein",
            selection2="resname SOL",
            output_file="protein_water_distance.csv"
        )
        
        # Hydrogen bond analysis
        hbond_analyzer = HydrogenBondAnalyzer(topology, trajectory)
        hbond_results = hbond_analyzer.calculate(
            donors_sel="protein",
            acceptors_sel="protein",
            output_file="hbonds.csv"
        )
        
        print(f"Mean H-bond count: {hbond_results['mean_count']:.2f}")
        
    else:
        print("Trajectory files not found - using placeholder")


def example_batch_processing_with_modules():
    """Example of batch processing with individual modules."""
    print("\n=== Batch Processing with Modules ===")

    def analyze_single_trajectory(traj_file, output_dir="batch_results"):
        """Analyze single trajectory using multiple modules."""
        try:
            Path(output_dir).mkdir(exist_ok=True)
            base_name = Path(traj_file).stem
            
            # Assume corresponding topology file exists
            topology = traj_file.replace('.xtc', '.tpr')
            
            if not os.path.exists(topology):
                return f"Topology not found for {traj_file}"
            
            results = {}
            
            # RMSD analysis
            rmsd_calc = RMSDCalculator(topology, traj_file)
            times, rmsd_vals = rmsd_calc.calculate_mdanalysis(
                output_file=f"{output_dir}/{base_name}_rmsd.csv"
            )
            results['rmsd_mean'] = float(np.mean(rmsd_vals))
            
            # Radius of gyration
            rg_calc = RadiusGyrationCalculator(topology, traj_file)
            times, rg_vals = rg_calc.calculate(
                output_file=f"{output_dir}/{base_name}_rg.csv"
            )
            results['rg_mean'] = float(np.mean(rg_vals))
            
            return results
            
        except Exception as e:
            return f"Error: {e}"

    # Find trajectory files
    trajectory_files = sorted(str(path) for path in Path(".").glob("*.xtc"))

    if trajectory_files:
        with ThreadPoolExecutor(max_workers=2) as executor:
            results = list(executor.map(analyze_single_trajectory, trajectory_files))
        print(f"Batch analysis completed: {results}")
    else:
        print("No trajectory files found for batch processing")


def example_plotting_integration():
    """Example of integrating plotting with analysis modules."""
    print("\n=== Plotting Integration ===")
    
    # Initialize plot manager
    plotter = PlotManager()
    
    # Example: plot RMSD data from modular calculation
    topology = "protein.tpr"
    trajectory = "trajectory.xtc"
    
    if os.path.exists(topology) and os.path.exists(trajectory):
        # Calculate RMSD
        rmsd_calc = RMSDCalculator(topology, trajectory)
        times, rmsd_values = rmsd_calc.calculate_mdanalysis()
        
        # Create DataFrame for plotting
        import pandas as pd
        df = pd.DataFrame({'Time': times, 'RMSD': rmsd_values})
        
        # Plot using PlotManager
        plotter.plot_time_series(
            df,
            title="RMSD Analysis from Modular Calculator",
            xlabel="Time (ps)",
            ylabel="RMSD (nm)",
            output_path="modular_rmsd_plot.png"
        )
        
        # Plot distribution
        plotter.plot_distribution(
            rmsd_values,
            title="RMSD Distribution",
            xlabel="RMSD (nm)",
            output_path="rmsd_distribution.png"
        )
        
    else:
        print("Files not found for plotting example")


def main():
    """Run all modular examples."""
    print("Immunex Toolkit - Modular Usage Examples")
    print("=" * 50)
    
    example_rmsd_analysis()
    example_structure_analysis()
    example_trajectory_modules()
    example_batch_processing_with_modules()
    example_plotting_integration()
    
    print("\n" + "=" * 50)
    print("Modular examples completed!")
    print("\nNew modular structure benefits:")
    print("- Focused, single-purpose modules")
    print("- Easier to test and maintain")
    print("- More flexible usage patterns")
    print("- Better code organization")


if __name__ == "__main__":
    main()
