"""
Auto Chain Identification Usage Examples

Demonstrates the new automatic chain identification feature for docking angle analysis.

Author: Immunex Development Team
Date: 2026-03-10
"""

import logging
from immunex.analysis.angles import DockingAnglePrimaryAnalyzer

logging.basicConfig(level=logging.INFO)


def example_1_fully_automatic():
    """
    Example 1: Fully automatic - zero configuration required.

    The analyzer automatically:
    1. Detects file format (PDB vs TPR)
    2. Extracts protein chain sequences
    3. Identifies chains using ANARCI + length heuristics
    4. Generates MDAnalysis selection strings
    """
    print("\n" + "=" * 80)
    print("Example 1: Fully Automatic Chain Identification")
    print("=" * 80)

    # PDB file
    analyzer_pdb = DockingAnglePrimaryAnalyzer('1ao7.pdb')
    crossing, incident = analyzer_pdb.calculate_docking_angles()
    print(f"\nPDB Results:")
    print(f"  Crossing: {crossing:.2f}°")
    print(f"  Incident: {incident:.2f}°")

    # TPR file
    analyzer_tpr = DockingAnglePrimaryAnalyzer('md.tpr', 'md.xtc')
    times, crossing_traj, incident_traj = analyzer_tpr.calculate_docking_angles_trajectory(
        stride=10,
        output_file='docking_angles.csv'
    )
    print(f"\nTPR Trajectory Results:")
    print(f"  Frames: {len(times)}")
    print(f"  Crossing: {crossing_traj.mean():.2f} ± {crossing_traj.std():.2f}°")
    print(f"  Incident: {incident_traj.mean():.2f} ± {incident_traj.std():.2f}°")


def example_2_manual_override():
    """
    Example 2: Manual override for specific chains.

    You can still manually specify selections if needed.
    """
    print("\n" + "=" * 80)
    print("Example 2: Manual Selection Override (Backward Compatible)")
    print("=" * 80)

    # Disable auto identification
    analyzer = DockingAnglePrimaryAnalyzer(
        'md.tpr',
        'md.xtc',
        auto_identify_chains=False
    )

    # Manually specify selections (old API still works)
    crossing, incident = analyzer.calculate_docking_angles(
        mhc_selection='segname PROA',
        tcr_alpha_selection='segname PROD',
        tcr_beta_selection='segname PROE'
    )

    print(f"\nManual Selection Results:")
    print(f"  Crossing: {crossing:.2f}°")
    print(f"  Incident: {incident:.2f}°")


def example_3_mixed_mode():
    """
    Example 3: Mixed mode - auto identify + manual override.

    Auto identify chains, then optionally override specific selections.
    """
    print("\n" + "=" * 80)
    print("Example 3: Mixed Mode (Auto + Manual Override)")
    print("=" * 80)

    # Auto identify all chains
    analyzer = DockingAnglePrimaryAnalyzer('md.tpr', 'md.xtc')

    # Access auto-identified selections
    print("\nAuto-identified selections:")
    for key, value in analyzer.auto_selections.items():
        print(f"  {key}: {value}")

    # Use auto for most, but manually override MHC selection
    crossing, incident = analyzer.calculate_docking_angles(
        mhc_selection='segname PROA and resid 1:180',  # Custom MHC range
        # tcr_alpha_selection and tcr_beta_selection will use auto-identified values
    )

    print(f"\nMixed Mode Results:")
    print(f"  Crossing: {crossing:.2f}°")
    print(f"  Incident: {incident:.2f}°")


def example_4_graceful_degradation():
    """
    Example 4: Graceful degradation with missing components.

    The analyzer can still work even if peptide or beta2m are missing.
    Only requires: HLA-alpha + at least 1 TCR chain.
    """
    print("\n" + "=" * 80)
    print("Example 4: Graceful Degradation (Missing Components)")
    print("=" * 80)

    # This will work even if peptide or beta2m are missing
    # because validation is in graceful mode (strict=False)
    analyzer = DockingAnglePrimaryAnalyzer('incomplete_complex.pdb')

    # Check what was identified
    print("\nIdentified chains:")
    for chain_id, info in analyzer.identifications.items():
        print(f"  {chain_id}: {info.chain_type} ({info.length} AA, conf={info.confidence:.2f})")

    # Calculate angles (will succeed if HLA-α + 1 TCR chain present)
    try:
        crossing, incident = analyzer.calculate_docking_angles()
        print(f"\nResults (with missing components):")
        print(f"  Crossing: {crossing:.2f}°")
        print(f"  Incident: {incident:.2f}°")
    except ValueError as e:
        print(f"\nError: {e}")


def example_5_batch_processing():
    """
    Example 5: Batch processing with automatic identification.

    Process multiple structures with zero manual configuration.
    """
    print("\n" + "=" * 80)
    print("Example 5: Batch Processing with Auto Identification")
    print("=" * 80)

    structures = [
        '1ao7.pdb',
        '1bd2.pdb',
        '1oga.pdb',
    ]

    results = []

    for structure in structures:
        try:
            analyzer = DockingAnglePrimaryAnalyzer(structure)
            crossing, incident = analyzer.calculate_docking_angles()
            results.append({
                'structure': structure,
                'crossing': crossing,
                'incident': incident,
                'status': 'success'
            })
            print(f"✓ {structure}: Crossing={crossing:.2f}°, Incident={incident:.2f}°")
        except Exception as e:
            results.append({
                'structure': structure,
                'status': 'failed',
                'error': str(e)
            })
            print(f"✗ {structure}: Failed - {e}")

    # Summary
    success_count = sum(1 for r in results if r['status'] == 'success')
    print(f"\n{success_count}/{len(structures)} structures processed successfully")


def example_6_inspection_mode():
    """
    Example 6: Inspect auto-identification results before calculation.

    Use this to verify chain identification before running analysis.
    """
    print("\n" + "=" * 80)
    print("Example 6: Inspection Mode - Verify Identification")
    print("=" * 80)

    # Initialize analyzer
    analyzer = DockingAnglePrimaryAnalyzer('md.tpr', 'md.xtc')

    # Inspect identification results
    print(analyzer.chain_adapter.get_chain_summary(analyzer.identifications))

    # Validate identification
    is_valid, messages = analyzer.chain_adapter.validate_identification(
        analyzer.identifications,
        strict=False
    )

    print(f"\nValidation: {'✓ PASS' if is_valid else '✗ FAIL'}")
    for msg in messages:
        print(f"  {msg}")

    # Access individual chain info
    for chain_id, info in analyzer.identifications.items():
        if info.chain_type == 'HLA_alpha':
            print(f"\nHLA-alpha chain details:")
            print(f"  Chain ID: {info.chain_id}")
            print(f"  Length: {info.length} AA")
            print(f"  Confidence: {info.confidence:.2f}")
            print(f"  Residue range: {info.residue_range}")
            print(f"  Naming: {info.naming_convention}")


if __name__ == "__main__":
    # Run all examples
    import sys

    examples = [
        example_1_fully_automatic,
        example_2_manual_override,
        example_3_mixed_mode,
        example_4_graceful_degradation,
        example_5_batch_processing,
        example_6_inspection_mode,
    ]

    print("\n" + "=" * 80)
    print("AUTO CHAIN IDENTIFICATION USAGE EXAMPLES")
    print("=" * 80)

    for i, example_func in enumerate(examples, 1):
        try:
            example_func()
        except Exception as e:
            print(f"\nExample {i} encountered error: {e}")

    print("\n" + "=" * 80)
    print("All examples completed")
    print("=" * 80)
