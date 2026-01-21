#!/usr/bin/env python3
"""
Batch Contact Interface RMSD Calculation

Calculate RMSD for residues within contact distance of peptide.
This provides a focused view of interface stability.
"""

import sys
import logging
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import pandas as pd
import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent))

from aftermd.analysis.trajectory.interface_rmsd import ContactInterfaceRMSDCalculator

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def process_single_task(task_info):
    """
    Process single task for interface RMSD calculation.

    Args:
        task_info: Dict with keys:
            - name: task name
            - topology: path to topology file
            - trajectory: path to trajectory file
            - output_dir: output directory
            - cutoff: distance cutoff (default 10.0 A)

    Returns:
        Dict with task results
    """
    task_name = task_info['name']
    topology = task_info['topology']
    trajectory = task_info['trajectory']
    output_dir = Path(task_info['output_dir'])
    cutoff = task_info.get('cutoff', 10.0)

    logger.info(f"="*60)
    logger.info(f"Processing task: {task_name}")
    logger.info(f"="*60)

    try:
        # Initialize calculator
        calc = ContactInterfaceRMSDCalculator(topology, trajectory)

        # Calculate interface RMSD
        output_file = output_dir / f"{task_name}_interface_rmsd.csv"

        times, rmsd_values = calc.calculate_interface_rmsd(
            cutoff=cutoff,
            output_file=str(output_file),
            verbose=False
        )

        # Calculate statistics
        result = {
            'task': task_name,
            'status': 'success',
            'mean_rmsd_nm': float(np.mean(rmsd_values)),
            'std_rmsd_nm': float(np.std(rmsd_values)),
            'min_rmsd_nm': float(np.min(rmsd_values)),
            'max_rmsd_nm': float(np.max(rmsd_values)),
            'n_frames': len(rmsd_values),
            'cutoff_angstrom': cutoff,
            'output_file': str(output_file)
        }

        logger.info(f"  -> Mean interface RMSD: {result['mean_rmsd_nm']:.4f} nm")

        return result

    except Exception as e:
        logger.error(f"Task {task_name} failed: {e}")
        return {
            'task': task_name,
            'status': 'failed',
            'error': str(e),
            'cutoff_angstrom': cutoff
        }


def main():
    """Main batch processing function."""

    # Discover all processed trajectory tasks
    input_dirs = [
        Path("input/pbc_1000frames_2step"),
        Path("input/pbc_1000frames_2step_patch")
    ]

    tasks = []
    for input_dir in input_dirs:
        if not input_dir.exists():
            logger.warning(f"Directory not found: {input_dir}")
            continue

        for task_dir in sorted(input_dir.iterdir()):
            if not task_dir.is_dir():
                continue

            task_name = task_dir.name
            topology = task_dir / "md.tpr"
            trajectory_files = list(task_dir.glob("*_processed.xtc"))

            if not topology.exists():
                logger.warning(f"Missing topology: {task_name}")
                continue

            if not trajectory_files:
                logger.warning(f"Missing processed trajectory: {task_name}")
                continue

            trajectory = trajectory_files[0]

            tasks.append({
                'name': task_name,
                'topology': str(topology),
                'trajectory': str(trajectory),
                'output_dir': 'output/rmsd_contactsurf',
                'cutoff': 10.0
            })

    logger.info("="*80)
    logger.info("Batch Contact Interface RMSD Calculation")
    logger.info("="*80)
    logger.info(f"Total tasks: {len(tasks)}")
    logger.info(f"Distance cutoff: 10.0 A (peptide contact interface)")
    logger.info(f"Alignment: pHLA (chains A+B+C)")
    logger.info(f"Output directory: output/rmsd_contactsurf")
    logger.info("="*80)

    if len(tasks) == 0:
        logger.error("No valid tasks found!")
        return

    # Create output directory
    output_base = Path("output/rmsd_contactsurf")
    output_base.mkdir(parents=True, exist_ok=True)

    # Process tasks in parallel
    max_workers = 4
    logger.info(f"\nStarting parallel processing with {max_workers} workers...")

    all_results = []

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(process_single_task, task): task for task in tasks}

        for future in as_completed(futures):
            task = futures[future]
            try:
                result = future.result()
                all_results.append(result)

                status = result.get('status', 'unknown')
                logger.info(f"\nCompleted: {result['task']} - {status}")

            except Exception as e:
                logger.error(f"\nTask {task['name']} raised exception: {e}")
                all_results.append({
                    'task': task['name'],
                    'status': 'exception',
                    'error': str(e)
                })

    # Generate summary report
    logger.info("\n" + "="*80)
    logger.info("Generating summary report...")
    logger.info("="*80)

    # Save summary CSV
    summary_df = pd.DataFrame(all_results)
    summary_file = output_base / "batch_interface_rmsd_summary.csv"
    summary_df.to_csv(summary_file, index=False)
    logger.info(f"\nSummary saved: {summary_file}")

    # Calculate statistics
    success_df = summary_df[summary_df['status'] == 'success']
    failed_df = summary_df[summary_df['status'] != 'success']

    logger.info("\n" + "="*80)
    logger.info("Batch Processing Complete")
    logger.info("="*80)
    logger.info(f"Total tasks: {len(all_results)}")
    logger.info(f"Successful: {len(success_df)}")
    logger.info(f"Failed: {len(failed_df)}")

    if len(success_df) > 0:
        logger.info("\n" + "="*80)
        logger.info("Interface RMSD Statistics (All Tasks)")
        logger.info("="*80)
        logger.info(f"N = {len(success_df)}")
        logger.info(f"Mean: {success_df['mean_rmsd_nm'].mean():.4f} ± {success_df['mean_rmsd_nm'].std():.4f} nm")
        logger.info(f"Median: {success_df['mean_rmsd_nm'].median():.4f} nm")
        logger.info(f"Range: {success_df['mean_rmsd_nm'].min():.4f} - {success_df['mean_rmsd_nm'].max():.4f} nm")

        # Generate statistics summary
        stats_summary = {
            'Metric': 'Interface_RMSD',
            'N': len(success_df),
            'Mean_nm': f"{success_df['mean_rmsd_nm'].mean():.4f}",
            'Std_nm': f"{success_df['mean_rmsd_nm'].std():.4f}",
            'Median_nm': f"{success_df['mean_rmsd_nm'].median():.4f}",
            'Min_nm': f"{success_df['mean_rmsd_nm'].min():.4f}",
            'Max_nm': f"{success_df['mean_rmsd_nm'].max():.4f}",
            'Q25_nm': f"{success_df['mean_rmsd_nm'].quantile(0.25):.4f}",
            'Q75_nm': f"{success_df['mean_rmsd_nm'].quantile(0.75):.4f}"
        }

        stats_df = pd.DataFrame([stats_summary])
        stats_file = output_base / "INTERFACE_RMSD_STATISTICS.csv"
        stats_df.to_csv(stats_file, index=False)
        logger.info(f"\nStatistics summary saved: {stats_file}")

    logger.info("\n" + "="*80)


if __name__ == "__main__":
    main()
