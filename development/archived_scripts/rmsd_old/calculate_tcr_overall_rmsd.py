#!/usr/bin/env python3
"""
Calculate overall TCR RMSD aligned to pHLA.
Compare with CDR3β RMSD to assess relative flexibility.
"""

import subprocess
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import logging
import pandas as pd
import numpy as np

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def calculate_tcr_rmsd(task_info):
    """Calculate TCR overall RMSD for a single task."""
    task_name, task_dir = task_info

    result = {
        'task_name': task_name,
        'status': 'failed',
        'error': None
    }

    try:
        topology = task_dir / "md.tpr"
        trajectory = list(task_dir.glob("*_processed.xtc"))[0]
        index_file = task_dir / "combined_cdr3_index.ndx"
        rmsd_file = task_dir / "tcr_overall_rmsd.xvg"

        # Check files exist
        if not topology.exists():
            result['error'] = 'md.tpr not found'
            return result

        if not index_file.exists():
            result['error'] = 'combined_cdr3_index.ndx not found'
            return result

        # Run gmx rms: align to pHLA, calculate TCR RMSD
        cmd = [
            'gmx', 'rms',
            '-s', str(topology),
            '-f', str(trajectory),
            '-n', str(index_file),
            '-o', str(rmsd_file),
            '-tu', 'ns'
        ]

        # Provide selections: pHLA for alignment, TCR for RMSD
        process = subprocess.Popen(
            cmd,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )

        stdout, stderr = process.communicate(input="pHLA\nTCR\n", timeout=120)

        if process.returncode != 0:
            result['error'] = f'gmx rms failed: {stderr}'
            return result

        # Parse RMSD data
        rmsd_data = []
        with open(rmsd_file, 'r') as f:
            for line in f:
                if not line.startswith(('#', '@')):
                    parts = line.split()
                    if len(parts) >= 2:
                        time = float(parts[0])
                        rmsd = float(parts[1])
                        rmsd_data.append([time, rmsd])

        if not rmsd_data:
            result['error'] = 'No RMSD data parsed'
            return result

        # Calculate statistics
        rmsd_values = [r[1] for r in rmsd_data]
        result['status'] = 'success'
        result['rmsd_mean'] = sum(rmsd_values) / len(rmsd_values)
        result['rmsd_std'] = (sum((x - result['rmsd_mean'])**2 for x in rmsd_values) / len(rmsd_values))**0.5
        result['rmsd_min'] = min(rmsd_values)
        result['rmsd_max'] = max(rmsd_values)
        result['num_frames'] = len(rmsd_values)

        logger.info(f"✓ {task_name}: Mean TCR RMSD = {result['rmsd_mean']:.4f} nm")

        return result

    except Exception as e:
        result['error'] = str(e)
        logger.error(f"✗ {task_name}: {e}")
        return result


def main():
    """Main function to calculate TCR overall RMSD."""

    logger.info("=" * 80)
    logger.info("Calculate Overall TCR RMSD (pHLA aligned)")
    logger.info("=" * 80)

    # Paths
    trajectory_base = Path("/home/xumy/work/development/Immunex/input/pbc_1000frames_2step")
    output_dir = Path("/home/xumy/work/development/Immunex/output/cdr3_analysis")

    # Load successful tasks from CDR3β summary
    cdr3_summary = output_dir / "cdr3_beta_rmsd_complete_summary.csv"
    df_cdr3 = pd.read_csv(cdr3_summary)
    successful_tasks = df_cdr3['task_name'].tolist()

    logger.info(f"Found {len(successful_tasks)} tasks with successful CDR3β RMSD")

    # Prepare task info
    tasks_to_process = []
    for task_name in successful_tasks:
        task_dir = trajectory_base / task_name
        if task_dir.exists():
            tasks_to_process.append((task_name, task_dir))

    logger.info(f"Processing {len(tasks_to_process)} tasks\n")

    # Process in parallel
    with ProcessPoolExecutor(max_workers=8) as executor:
        results = list(executor.map(calculate_tcr_rmsd, tasks_to_process))

    # Summary
    successful = [r for r in results if r['status'] == 'success']
    failed = [r for r in results if r['status'] == 'failed']

    logger.info("\n" + "=" * 80)
    logger.info("TCR Overall RMSD Calculation Summary")
    logger.info("=" * 80)
    logger.info(f"Total tasks: {len(results)}")
    logger.info(f"Successful: {len(successful)}")
    logger.info(f"Failed: {len(failed)}")
    logger.info(f"Success rate: {len(successful)/len(results)*100:.1f}%")

    if successful:
        # Save summary CSV
        df_tcr = pd.DataFrame(successful)
        df_tcr = df_tcr[['task_name', 'rmsd_mean', 'rmsd_std', 'rmsd_min', 'rmsd_max', 'num_frames']]
        df_tcr.columns = ['task_name', 'tcr_rmsd_mean_nm', 'tcr_rmsd_std_nm', 'tcr_rmsd_min_nm', 'tcr_rmsd_max_nm', 'num_frames']

        tcr_summary_csv = output_dir / "tcr_overall_rmsd_summary.csv"
        df_tcr.to_csv(tcr_summary_csv, index=False)

        logger.info(f"\n✓ TCR RMSD summary saved to: {tcr_summary_csv}")

        # Statistics
        mean_rmsds = [r['rmsd_mean'] for r in successful]
        overall_mean = sum(mean_rmsds) / len(mean_rmsds)
        overall_std = (sum((x - overall_mean)**2 for x in mean_rmsds) / len(mean_rmsds))**0.5

        logger.info(f"\nTCR Overall RMSD Statistics:")
        logger.info(f"  Mean: {overall_mean:.4f} ± {overall_std:.4f} nm")
        logger.info(f"  Range: {min(mean_rmsds):.4f} - {max(mean_rmsds):.4f} nm")

        # Compare with CDR3β
        logger.info("\nComparison with CDR3β RMSD:")
        cdr3_mean = df_cdr3['rmsd_mean_nm'].mean()
        cdr3_std = df_cdr3['rmsd_mean_nm'].std()
        logger.info(f"  CDR3β Mean: {cdr3_mean:.4f} ± {cdr3_std:.4f} nm")
        logger.info(f"  TCR/CDR3β Ratio: {overall_mean/cdr3_mean:.2f}x")

        # Merge datasets for detailed comparison
        df_merged = df_cdr3.merge(df_tcr, on='task_name', how='inner')
        df_merged['rmsd_ratio'] = df_merged['tcr_rmsd_mean_nm'] / df_merged['rmsd_mean_nm']

        comparison_csv = output_dir / "tcr_cdr3_comparison.csv"
        df_merged.to_csv(comparison_csv, index=False)
        logger.info(f"\n✓ Detailed comparison saved to: {comparison_csv}")

    if failed:
        logger.info(f"\n✗ Failed tasks: {len(failed)}")
        for r in failed[:10]:
            logger.info(f"  {r['task_name']}: {r['error']}")

    logger.info("=" * 80)


if __name__ == "__main__":
    main()
