#!/usr/bin/env python3
"""
Test PBC-RMSD Pipeline with Real MD Data

This script tests the complete preprocessing pipeline with real trajectory data:
1. PBC correction (2-step method)
2. Post-PBC validation
3. RMSD calculation
4. Convergence analysis
5. Quality grading and reporting

Author: Immunex Development Team
Date: 2026-03-18
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

from immunex.pipeline.pbc_rmsd_pipeline import PBCRMSDPipeline
import logging

logging.basicConfig(
    level=logging.INFO,
    format='%(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def test_pipeline_fast():
    """Test pipeline with downsampled trajectory (faster)"""

    # Test data location
    test_data_dir = Path("development/workspaces/FEL_workspace/input/1ao7/standard")

    if not test_data_dir.exists():
        logger.error(f"Test data directory not found: {test_data_dir}")
        logger.info("Expected location: development/workspaces/FEL_workspace/input/1ao7/standard")
        return 1

    trajectory = str(test_data_dir / "md.xtc")
    topology = str(test_data_dir / "md.tpr")

    # Check files exist
    if not Path(trajectory).exists():
        logger.error(f"Trajectory file not found: {trajectory}")
        return 1
    if not Path(topology).exists():
        logger.error(f"Topology file not found: {topology}")
        return 1

    # Get file size
    traj_size_mb = Path(trajectory).stat().st_size / (1024**2)
    logger.info(f"Trajectory size: {traj_size_mb:.1f} MB")

    # Output directory
    output_dir = Path("development/test_output/pipeline_test_real_data")
    output_dir.mkdir(parents=True, exist_ok=True)

    print()
    print("=" * 70)
    print("Testing PBC-RMSD Pipeline with Real MD Data")
    print("=" * 70)
    print()
    print(f"Trajectory: {trajectory}")
    print(f"Topology: {topology}")
    print(f"Output: {output_dir}")
    print(f"Mode: Fast (dt=100 ps downsampling)")
    print()
    print("This will take ~1-2 minutes...")
    print()

    # Initialize pipeline
    pipeline = PBCRMSDPipeline(
        gmx_executable="gmx",
        max_com_drift=1.0,
        max_rg_std_ratio=0.15,
        max_frame_jump_rmsd=0.5
    )

    # Process trajectory
    try:
        results = pipeline.process_single_trajectory(
            trajectory=trajectory,
            topology=topology,
            output_dir=str(output_dir),
            pbc_method="2step",
            dt=100.0,  # Sample every 100 ps for speed
            rmsd_selection="backbone",
            run_quality_check=True,
            validation_stride=1,
            generate_report=True
        )

        # Print results
        print()
        print("=" * 70)
        print("Pipeline Results")
        print("=" * 70)
        print()

        print(f"Overall Grade: {results['overall_grade']}")
        print(f"Qualified: {'YES' if results['is_qualified'] else 'NO'}")
        print()

        print(f"PBC Output: {results['pbc_output']}")
        print()

        # RMSD metrics
        if results.get('rmsd_metrics'):
            metrics = results['rmsd_metrics']
            print("RMSD Metrics:")
            print(f"  Mean RMSD: {metrics['mean_rmsd']:.3f} nm")
            print(f"  Std RMSD: {metrics['std_rmsd']:.3f} nm")
            print(f"  Min RMSD: {metrics['min_rmsd']:.3f} nm")
            print(f"  Max RMSD: {metrics['max_rmsd']:.3f} nm")
            print(f"  Converged: {'Yes' if metrics['is_converged'] else 'No'}")
            print(f"  Convergence Time: {metrics['convergence_time']:.1f}")
            print()

        # Post-PBC validation
        if results.get('post_pbc_validation'):
            pbc_val = results['post_pbc_validation']
            print("Post-PBC Validation:")
            print(f"  Status: {pbc_val['overall_status'].upper()}")
            print(f"  Grade: {pbc_val['overall_grade']}")

            if pbc_val.get('pbc_quality'):
                pbc_q = pbc_val['pbc_quality']
                print(f"  COM Drift: {pbc_q.get('com_drift_nm', 0):.3f} nm")

            if pbc_val.get('stability'):
                stab = pbc_val['stability']
                print(f"  Rg Mean: {stab.get('rg_mean_nm', 0):.3f} nm")
                print(f"  Rg Std/Mean: {stab.get('rg_std_ratio', 0):.3f}")

            if pbc_val.get('all_issues') and len(pbc_val['all_issues']) > 0:
                print(f"  Issues: {len(pbc_val['all_issues'])} found")
                for issue in pbc_val['all_issues'][:3]:  # Show first 3
                    print(f"    - {issue}")

            print()

        # Report files
        if results.get('report_files'):
            print("Generated Report Files:")
            for report_file in results['report_files']:
                print(f"  - {report_file}")
            print()

        print("=" * 70)
        print("Pipeline Test Completed Successfully!")
        print("=" * 70)
        print()

        # Summary interpretation
        print("Interpretation:")
        if results['overall_grade'] == 'A':
            print("  Excellent quality - ready for analysis")
        elif results['overall_grade'] == 'B':
            print("  Good quality - suitable for most analyses")
        elif results['overall_grade'] == 'C':
            print("  Acceptable quality - use with caution")
        else:
            print("  Poor quality - may need re-simulation")
        print()

        return 0

    except Exception as e:
        print()
        print("=" * 70)
        print("Pipeline Test Failed!")
        print("=" * 70)
        print(f"Error: {str(e)}")
        print()
        import traceback
        traceback.print_exc()
        return 1


def main():
    return test_pipeline_fast()


if __name__ == "__main__":
    sys.exit(main())
