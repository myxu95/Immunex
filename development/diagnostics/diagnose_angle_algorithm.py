#!/usr/bin/env python3
"""
Diagnose angle calculation algorithm stability.

Investigates why docking angles have large variations by examining:
1. TCR axis stability (cysteine positions, COM fluctuations)
2. MHC groove geometry stability (axis, normal vector consistency)
3. Frame-by-frame geometric quantity changes
"""

import sys
from pathlib import Path
import numpy as np
import MDAnalysis as mda

sys.path.insert(0, str(Path(__file__).parent.parent))

from immunex.analysis.angles import DockingAngleAnalyzer, DockingAngleInput


def diagnose_geometry_stability(topology, trajectory, output_dir="development/angle_diagnosis"):
    """
    Diagnose geometric stability issues in angle calculation.

    Parameters
    ----------
    topology : str
        Path to topology file
    trajectory : str
        Path to trajectory file
    output_dir : str
        Output directory for diagnostic plots
    """
    print("=" * 70)
    print("Angle Calculation Algorithm Diagnosis")
    print("=" * 70)

    # Load universe
    u = mda.Universe(topology, trajectory)

    # Create analyzer
    analyzer = DockingAngleAnalyzer()

    # Use auto-identification to get chains
    input_params = DockingAngleInput(
        topology=topology,
        trajectory=trajectory,
        stride=5,
        output_dir=output_dir,
        auto_identify_chains=True,
        use_anarci=True
    )

    # Get chain selections
    from immunex.utils import ChainIdentificationAdapter, SelectionStringBuilder

    adapter = ChainIdentificationAdapter(use_anarci=True)
    identifications = adapter.identify_chains(topology, None)
    builder = SelectionStringBuilder()
    selections = builder.build_all_selections(identifications)

    mhc_sel = selections['mhc_selection']
    tcr_alpha_sel = selections['tcr_alpha_selection']
    tcr_beta_sel = selections['tcr_beta_selection']

    print(f"\nChain selections:")
    print(f"  MHC: {mhc_sel}")
    print(f"  TCR Alpha: {tcr_alpha_sel}")
    print(f"  TCR Beta: {tcr_beta_sel}")

    # Initialize detectors (same as in analyzer)
    from immunex.analysis.angles.analyzer import _MHCGrooveDetector, _TCRAxisCalculator

    groove_detector = _MHCGrooveDetector(u, mhc_sel)
    tcr_calculator = _TCRAxisCalculator(u)

    # Arrays for tracking
    stride = 5
    n_frames = len(range(0, len(u.trajectory), stride))

    # TCR geometry tracking
    tcr_axis_vectors = np.zeros((n_frames, 3))
    alpha_com_positions = np.zeros((n_frames, 3))
    beta_com_positions = np.zeros((n_frames, 3))
    tcr_axis_lengths = np.zeros(n_frames)

    # MHC geometry tracking
    groove_axis_vectors = np.zeros((n_frames, 3))
    groove_normal_vectors = np.zeros((n_frames, 3))

    # Conserved cysteine tracking
    alpha_cys_resids = []
    beta_cys_resids = []

    print(f"\nAnalyzing {n_frames} frames (stride={stride})...")
    print("-" * 70)

    for i, ts in enumerate(u.trajectory[::stride]):
        # TCR axis calculation
        alpha_cys_info = tcr_calculator.detect_conserved_cysteines(tcr_alpha_sel, "TCR")
        beta_cys_info = tcr_calculator.detect_conserved_cysteines(tcr_beta_sel, "TCR")

        # Track which residues are detected
        alpha_cys_resids.append((
            alpha_cys_info['cys_23_resid'],
            alpha_cys_info['cys_104_resid']
        ))
        beta_cys_resids.append((
            beta_cys_info['cys_23_resid'],
            beta_cys_info['cys_104_resid']
        ))

        # Get CA atoms
        alpha_cys_ca = u.select_atoms(
            f"{alpha_cys_info['cys_23_ca_selection']} or {alpha_cys_info['cys_104_ca_selection']}"
        )
        beta_cys_ca = u.select_atoms(
            f"{beta_cys_info['cys_23_ca_selection']} or {beta_cys_info['cys_104_ca_selection']}"
        )

        com_alpha = alpha_cys_ca.center_of_mass()
        com_beta = beta_cys_ca.center_of_mass()

        tcr_axis_raw = com_beta - com_alpha
        tcr_axis_length = np.linalg.norm(tcr_axis_raw)
        tcr_axis = tcr_axis_raw / tcr_axis_length

        # MHC groove geometry
        groove_axis, _, _, _ = groove_detector.detect_groove_axis()
        _, groove_normal, _, _ = groove_detector.detect_groove_plane()

        # Store data
        tcr_axis_vectors[i] = tcr_axis
        alpha_com_positions[i] = com_alpha
        beta_com_positions[i] = com_beta
        tcr_axis_lengths[i] = tcr_axis_length

        groove_axis_vectors[i] = groove_axis
        groove_normal_vectors[i] = groove_normal

        if (i + 1) % 10 == 0:
            print(f"Frame {i+1}/{n_frames}: "
                  f"TCR_axis_length={tcr_axis_length:.2f}Å, "
                  f"Alpha_COM={com_alpha[:2]}, Beta_COM={com_beta[:2]}")

    print("-" * 70)
    print("\nDiagnostic Analysis:")
    print("=" * 70)

    # 1. Check conserved cysteine consistency
    print("\n1. Conserved Cysteine Detection Consistency:")
    print("-" * 70)

    # Check if same residues are detected in all frames
    alpha_cys_23_unique = set([r[0] for r in alpha_cys_resids])
    alpha_cys_104_unique = set([r[1] for r in alpha_cys_resids])
    beta_cys_23_unique = set([r[0] for r in beta_cys_resids])
    beta_cys_104_unique = set([r[1] for r in beta_cys_resids])

    print(f"  Alpha Cys23 unique resids: {alpha_cys_23_unique}")
    print(f"  Alpha Cys104 unique resids: {alpha_cys_104_unique}")
    print(f"  Beta Cys23 unique resids: {beta_cys_23_unique}")
    print(f"  Beta Cys104 unique resids: {beta_cys_104_unique}")

    cys_consistent = (
        len(alpha_cys_23_unique) == 1 and
        len(alpha_cys_104_unique) == 1 and
        len(beta_cys_23_unique) == 1 and
        len(beta_cys_104_unique) == 1
    )

    if cys_consistent:
        print("\n  ✓ CONSISTENT: Same cysteines detected in all frames")
    else:
        print("\n  ✗ PROBLEM: Different cysteines detected in different frames!")
        print("  This causes TCR axis instability.")

    # 2. TCR axis stability
    print("\n2. TCR Axis Stability:")
    print("-" * 70)

    # Length variation
    mean_length = np.mean(tcr_axis_lengths)
    std_length = np.std(tcr_axis_lengths)
    cv_length = std_length / mean_length

    print(f"  TCR axis length: {mean_length:.2f} ± {std_length:.2f} Å (CV={cv_length:.4f})")
    print(f"  Range: [{np.min(tcr_axis_lengths):.2f}, {np.max(tcr_axis_lengths):.2f}] Å")

    if cv_length > 0.05:
        print(f"  ⚠ WARNING: Large length variation (CV > 0.05)")
    else:
        print(f"  ✓ OK: Length is stable (CV < 0.05)")

    # Direction variation (angular deviation from mean)
    mean_tcr_axis = np.mean(tcr_axis_vectors, axis=0)
    mean_tcr_axis /= np.linalg.norm(mean_tcr_axis)

    angular_deviations = []
    for vec in tcr_axis_vectors:
        cos_angle = np.dot(vec, mean_tcr_axis)
        cos_angle = np.clip(cos_angle, -1.0, 1.0)
        angle_deg = np.degrees(np.arccos(cos_angle))
        angular_deviations.append(angle_deg)

    angular_deviations = np.array(angular_deviations)
    mean_dev = np.mean(angular_deviations)
    std_dev = np.std(angular_deviations)
    max_dev = np.max(angular_deviations)

    print(f"\n  TCR axis angular deviation from mean:")
    print(f"    Mean: {mean_dev:.2f}°")
    print(f"    Std: {std_dev:.2f}°")
    print(f"    Max: {max_dev:.2f}°")

    if max_dev > 20.0:
        print(f"  ✗ PROBLEM: Large angular deviation (max > 20°)")
        print(f"    The TCR axis direction changes significantly between frames.")
    elif max_dev > 10.0:
        print(f"  ⚠ WARNING: Moderate angular deviation (10° < max < 20°)")
    else:
        print(f"  ✓ OK: Stable TCR axis direction (max < 10°)")

    # 3. MHC groove geometry stability
    print("\n3. MHC Groove Geometry Stability:")
    print("-" * 70)

    # Groove axis angular deviation
    mean_groove_axis = np.mean(groove_axis_vectors, axis=0)
    mean_groove_axis /= np.linalg.norm(mean_groove_axis)

    groove_axis_deviations = []
    for vec in groove_axis_vectors:
        cos_angle = np.dot(vec, mean_groove_axis)
        cos_angle = np.clip(cos_angle, -1.0, 1.0)
        angle_deg = np.degrees(np.arccos(cos_angle))
        groove_axis_deviations.append(angle_deg)

    groove_axis_deviations = np.array(groove_axis_deviations)

    print(f"  Groove axis angular deviation from mean:")
    print(f"    Mean: {np.mean(groove_axis_deviations):.2f}°")
    print(f"    Std: {np.std(groove_axis_deviations):.2f}°")
    print(f"    Max: {np.max(groove_axis_deviations):.2f}°")

    if np.max(groove_axis_deviations) > 15.0:
        print(f"  ✗ PROBLEM: Large groove axis deviation (max > 15°)")
    elif np.max(groove_axis_deviations) > 10.0:
        print(f"  ⚠ WARNING: Moderate groove axis deviation")
    else:
        print(f"  ✓ OK: Stable groove axis")

    # Groove normal angular deviation
    mean_groove_normal = np.mean(groove_normal_vectors, axis=0)
    mean_groove_normal /= np.linalg.norm(mean_groove_normal)

    groove_normal_deviations = []
    for vec in groove_normal_vectors:
        cos_angle = np.dot(vec, mean_groove_normal)
        cos_angle = np.clip(cos_angle, -1.0, 1.0)
        angle_deg = np.degrees(np.arccos(cos_angle))
        groove_normal_deviations.append(angle_deg)

    groove_normal_deviations = np.array(groove_normal_deviations)

    print(f"\n  Groove normal angular deviation from mean:")
    print(f"    Mean: {np.mean(groove_normal_deviations):.2f}°")
    print(f"    Std: {np.std(groove_normal_deviations):.2f}°")
    print(f"    Max: {np.max(groove_normal_deviations):.2f}°")

    if np.max(groove_normal_deviations) > 15.0:
        print(f"  ✗ PROBLEM: Large groove normal deviation (max > 15°)")
    elif np.max(groove_normal_deviations) > 10.0:
        print(f"  ⚠ WARNING: Moderate groove normal deviation")
    else:
        print(f"  ✓ OK: Stable groove normal")

    # 4. Summary
    print("\n" + "=" * 70)
    print("Summary:")
    print("=" * 70)

    issues = []
    if not cys_consistent:
        issues.append("- Inconsistent cysteine detection → Use fixed reference structure")
    if cv_length > 0.05:
        issues.append("- Large TCR axis length variation")
    if max_dev > 20.0:
        issues.append(f"- Large TCR axis angular deviation ({max_dev:.1f}°)")
    if np.max(groove_axis_deviations) > 15.0:
        issues.append(f"- Large MHC groove axis deviation ({np.max(groove_axis_deviations):.1f}°)")
    if np.max(groove_normal_deviations) > 15.0:
        issues.append(f"- Large MHC groove normal deviation ({np.max(groove_normal_deviations):.1f}°)")

    if issues:
        print("\nProblems identified:")
        for issue in issues:
            print(f"  {issue}")

        print("\nRecommendations:")
        print("  1. Use a REFERENCE structure (e.g., first frame or average) for:")
        print("     - Conserved cysteine detection (do once, not per frame)")
        print("     - MHC sequence alignment (do once, not per frame)")
        print("  2. Calculate TCR axis and MHC groove geometry relative to:")
        print("     - Fixed cysteine residue IDs")
        print("     - Fixed α1/α2 region residue IDs")
        print("  3. Only update atomic POSITIONS per frame, not selections")
    else:
        print("\n✓ No major stability issues detected.")
        print("  The large angle variations may be REAL conformational changes.")

    print("\n" + "=" * 70)
    print("Diagnostic complete.")
    print("=" * 70)

    return {
        'cys_consistent': cys_consistent,
        'tcr_axis_deviation_max': max_dev,
        'groove_axis_deviation_max': np.max(groove_axis_deviations),
        'groove_normal_deviation_max': np.max(groove_normal_deviations),
        'tcr_axis_length_cv': cv_length
    }


if __name__ == '__main__':
    # Test data
    base_dir = Path("development/workspaces/FEL_workspace")

    test_case = {
        'pdb': base_dir / '1bd2_basin_structures/1bd2_rest2_basin1.pdb',
        'xtc': base_dir / '1bd2_basin_trajectories/1bd2_rest2_basin1_frames0-151.xtc',
    }

    if not test_case['pdb'].exists() or not test_case['xtc'].exists():
        print("Test data not found. Exiting.")
        sys.exit(1)

    results = diagnose_geometry_stability(
        str(test_case['pdb']),
        str(test_case['xtc'])
    )

    print(f"\nDiagnostic results: {results}")
