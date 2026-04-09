#!/usr/bin/env python3
"""
Test the fixed angle calculation algorithm with reference structure strategy.

Should show:
1. Consistent cysteine detection (same residues in all frames)
2. Reduced geometric instability
3. More reasonable angle variations
"""

import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent))

from immunex.analysis.angles import DockingAngleAnalyzer, DockingAngleInput


def test_fixed_algorithm():
    """Test the fixed angle calculation algorithm."""
    print("=" * 70)
    print("Testing Fixed Angle Calculation Algorithm")
    print("=" * 70)

    # Test data
    base_dir = Path("development/workspaces/FEL_workspace")
    test_case = {
        'pdb': base_dir / '1bd2_basin_structures/1bd2_rest2_basin1.pdb',
        'xtc': base_dir / '1bd2_basin_trajectories/1bd2_rest2_basin1_frames0-151.xtc',
    }

    if not test_case['pdb'].exists() or not test_case['xtc'].exists():
        print("Test data not found. Exiting.")
        sys.exit(1)

    # Create analyzer
    analyzer = DockingAngleAnalyzer()

    # Configure input
    input_params = DockingAngleInput(
        topology=str(test_case['pdb']),
        trajectory=str(test_case['xtc']),
        stride=5,
        output_dir="test_output/angle_analysis_fixed",
        auto_identify_chains=True,
        use_anarci=True
    )

    print("\nRunning angle analysis with REFERENCE STRUCTURE strategy...")
    print("-" * 70)

    # Execute analysis
    result = analyzer.analyze(input_params)

    if not result.success:
        print(f"\nAnalysis failed: {result.error_message}")
        sys.exit(1)

    print("\n" + "=" * 70)
    print("Analysis Results")
    print("=" * 70)

    # Display statistics
    stats = result.statistics
    print(f"\nFrames analyzed: {stats['n_frames']}")

    print(f"\nCrossing Angle:")
    print(f"  Mean: {stats['crossing_mean']:.2f}°")
    print(f"  Std: {stats['crossing_std']:.2f}°")
    print(f"  Range: [{stats['crossing_min']:.2f}, {stats['crossing_max']:.2f}]°")
    print(f"  Variation: {stats['crossing_max'] - stats['crossing_min']:.2f}°")
    print(f"  CV: {stats['crossing_std'] / stats['crossing_mean']:.4f}")

    print(f"\nIncident Angle:")
    print(f"  Mean: {stats['incident_mean']:.2f}°")
    print(f"  Std: {stats['incident_std']:.2f}°")
    print(f"  Range: [{stats['incident_min']:.2f}, {stats['incident_max']:.2f}]°")
    print(f"  Variation: {stats['incident_max'] - stats['incident_min']:.2f}°")
    print(f"  CV: {stats['incident_std'] / stats['incident_mean']:.4f}")

    # Reference residue info
    if 'reference_residues' in result.metadata:
        ref_res = result.metadata['reference_residues']
        print(f"\nReference Residues (fixed for all frames):")
        print(f"  TCR Alpha Cys: {ref_res['tcr_alpha_cys']}")
        print(f"  TCR Beta Cys: {ref_res['tcr_beta_cys']}")
        print(f"  MHC α1: {ref_res['mhc_alpha1_count']} residues")
        print(f"  MHC α2: {ref_res['mhc_alpha2_count']} residues")

    # Output files
    print(f"\nOutput files:")
    for f in result.output_files:
        print(f"  {f}")

    # Compare with old results (if available)
    old_csv = Path("test_output/angle_analysis/docking_angles.csv")
    new_csv = Path("test_output/angle_analysis_fixed/docking_angles.csv")

    if old_csv.exists() and new_csv.exists():
        print("\n" + "=" * 70)
        print("Comparison with Old Algorithm")
        print("=" * 70)

        import pandas as pd
        old_df = pd.read_csv(old_csv)
        new_df = pd.read_csv(new_csv)

        print(f"\nCrossing Angle:")
        print(f"  Old: {old_df['Crossing(deg)'].mean():.2f} ± {old_df['Crossing(deg)'].std():.2f}° "
              f"(range: {old_df['Crossing(deg)'].max() - old_df['Crossing(deg)'].min():.2f}°)")
        print(f"  New: {new_df['Crossing(deg)'].mean():.2f} ± {new_df['Crossing(deg)'].std():.2f}° "
              f"(range: {new_df['Crossing(deg)'].max() - new_df['Crossing(deg)'].min():.2f}°)")

        print(f"\nIncident Angle:")
        print(f"  Old: {old_df['Incident(deg)'].mean():.2f} ± {old_df['Incident(deg)'].std():.2f}° "
              f"(range: {old_df['Incident(deg)'].max() - old_df['Incident(deg)'].min():.2f}°)")
        print(f"  New: {new_df['Incident(deg)'].mean():.2f} ± {new_df['Incident(deg)'].std():.2f}° "
              f"(range: {new_df['Incident(deg)'].max() - new_df['Incident(deg)'].min():.2f}°)")

        # Calculate improvement
        old_crossing_range = old_df['Crossing(deg)'].max() - old_df['Crossing(deg)'].min()
        new_crossing_range = new_df['Crossing(deg)'].max() - new_df['Crossing(deg)'].min()
        improvement = ((old_crossing_range - new_crossing_range) / old_crossing_range) * 100

        print(f"\nCrossing angle range reduction: {improvement:.1f}%")

    print("\n" + "=" * 70)
    print("Test Complete")
    print("=" * 70)

    return True


if __name__ == '__main__':
    try:
        success = test_fixed_algorithm()
        sys.exit(0 if success else 1)
    except Exception as e:
        print(f"\nTest failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
