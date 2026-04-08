#!/usr/bin/env python3
"""
Comprehensive Immunex Toolkit Usage Examples.

This script demonstrates the complete workflow and capabilities of Immunex,
including path management, trajectory analysis, structure analysis, plotting,
PBC preprocessing, and batch processing.
"""

import os
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from immunex import (
    PathManager,
    PlotManager,
    RMSDCalculator,
    BFactorAnalyzer
)
from immunex.analysis import PBCProcessor


def example_path_management():
    """Example of using PathManager."""
    print("=== Path Management Example ===")
    
    # Initialize path manager
    pm = PathManager("/path/to/project")
    
    # Create project structure
    pm.create_project_structure()
    
    # Get specific paths
    results_path = pm.get_path('results')
    plots_path = pm.get_path('plots')
    
    # Get timestamped output path
    timestamped_path = pm.get_timestamped_path('analysis', 'rmsd_analysis')
    
    print(f"Results path: {results_path}")
    print(f"Plots path: {plots_path}")
    print(f"Timestamped path: {timestamped_path}")
    
    # Print path structure
    print(pm)


def example_trajectory_analysis():
    """Example of trajectory analysis."""
    print("\n=== Trajectory Analysis Example ===")

    # Initialize RMSD calculator (replace with actual file paths)
    topology = "protein.tpr"
    trajectory = "trajectory.xtc"

    if os.path.exists(topology) and os.path.exists(trajectory):
        # Calculate RMSD using RMSDCalculator
        rmsd_calc = RMSDCalculator(topology, trajectory)

        times, rmsd = rmsd_calc.calculate_mdanalysis(
            selection="protein and name CA",
            output_file="rmsd_analysis.csv"
        )

        print(f"RMSD calculation completed")

    else:
        print("Trajectory files not found - using placeholder")


def example_structure_analysis():
    """Example of structure analysis."""
    print("\n=== Structure Analysis Example ===")

    # Initialize B-factor analyzer (replace with actual PDB file)
    structure_file = "protein.pdb"

    if os.path.exists(structure_file):
        # Extract B-factors using BFactorAnalyzer
        bfactor_analyzer = BFactorAnalyzer(structure_file)

        indices, bfactors = bfactor_analyzer.extract_bfactors(
            selection="protein",
            output_file="bfactors.csv"
        )

        # Analyze B-factors by residue
        residue_stats = bfactor_analyzer.analyze_bfactor_by_residue(
            output_file="bfactor_by_residue.csv"
        )

        print(f"B-factor analysis completed")

    else:
        print("Structure file not found - using placeholder")


def example_plotting():
    """Example of plotting capabilities."""
    print("\n=== Plotting Example ===")
    
    # Initialize plot manager
    plotter = PlotManager()
    
    # Example XVG file plotting (replace with actual file)
    xvg_file = "rmsd.xvg"
    
    if os.path.exists(xvg_file):
        # Plot time series
        plotter.plot_time_series(
            xvg_file,
            title="RMSD Over Time",
            xlabel="Time (ps)",
            ylabel="RMSD (nm)",
            output_path="rmsd_plot.png"
        )
        
        # Plot RMSD specifically
        plotter.plot_rmsd(xvg_file, output_path="rmsd_analysis.png")
        
    else:
        print("XVG file not found - skipping plotting example")


def example_preprocessing():
    """Example of PBC processing workflow."""
    print("\n=== PBC Processing Example ===")
    
    # Initialize PBC processor
    pbc_processor = PBCProcessor()
    
    # Example files (replace with actual paths)
    trajectory = "raw_trajectory.xtc"
    topology = "system.tpr"
    
    if os.path.exists(trajectory) and os.path.exists(topology):
        # Remove PBC with trajectory downsampling
        processed_traj = pbc_processor.remove_pbc(
            trajectory=trajectory,
            topology=topology,
            output="processed_trajectory.xtc",
            dt=10.0  # Sample every 10 ps
        )
        
        # Comprehensive PBC processing with downsampling
        results = pbc_processor.comprehensive_pbc_process(
            trajectory=trajectory,
            topology=topology,
            output_dir="processed_data",
            dt=20.0  # Sample every 20 ps
        )
        
        print(f"PBC processing results: {results}")
        
    else:
        print("Input files not found - using placeholder")


def example_batch_processing():
    """Example of batch processing."""
    print("\n=== Batch Processing Example ===")

    # Find trajectory files
    trajectory_files = sorted(str(path) for path in Path(".").glob("*.xtc"))
    print(f"Found {len(trajectory_files)} trajectory files")
    
    # Example processing function
    def analyze_trajectory(traj_file, output_dir="results"):
        """Example analysis function for batch processing."""
        try:
            # This would be replaced with actual analysis
            print(f"Processing: {traj_file}")
            return f"Processed {Path(traj_file).name}"
        except Exception as e:
            print(f"Error processing {traj_file}: {e}")
            return None
    
    if trajectory_files:
        with ThreadPoolExecutor(max_workers=4) as executor:
            results = list(executor.map(analyze_trajectory, trajectory_files))
        print(f"Batch processing summary: {results}")
    else:
        print("No trajectory files found for batch processing")


def main():
    """Run all examples."""
    print("Immunex Toolkit - Comprehensive Usage Examples")
    print("=" * 50)
    
    example_path_management()
    example_trajectory_analysis()
    example_structure_analysis()
    example_plotting()
    example_preprocessing()
    example_batch_processing()
    
    print("\n" + "=" * 50)
    print("Examples completed!")
    print("Note: Replace placeholder file paths with actual data files")


if __name__ == "__main__":
    main()
