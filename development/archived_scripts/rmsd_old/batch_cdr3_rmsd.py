#!/usr/bin/env python3
"""
Batch CDR3 Beta RMSD Calculation Script

Calculates CDR3β RMSD for all tasks with TCR alignment.
"""

import os
import sys
import subprocess
import json
import numpy as np
import pandas as pd
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
from typing import Dict, Any, Tuple
import logging

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def process_single_task(task_info: Tuple[str, Path, Path]) -> Dict[str, Any]:
    """
    Process single task for CDR3β RMSD calculation.

    Args:
        task_info: (task_name, task_dir, cdr3_analysis_dir)

    Returns:
        Dictionary with RMSD statistics
    """
    task_name, task_dir, cdr3_analysis_dir = task_info

    result = {
        'task_name': task_name,
        'status': 'failed',
        'error': None
    }

    try:
        # Check required files
        trajectory = task_dir / f"{task_name}_2step_processed.xtc"
        topology = task_dir / "md.tpr"
        phla_index = task_dir / "phla_tcr.ndx"
        cdr3_index = cdr3_analysis_dir / task_name / "cdr3_only.ndx"

        if not trajectory.exists():
            result['error'] = 'Trajectory not found'
            return result

        if not topology.exists():
            result['error'] = 'Topology not found'
            return result

        if not phla_index.exists():
            result['error'] = 'pHLA-TCR index not found'
            return result

        if not cdr3_index.exists():
            result['error'] = 'CDR3 index not found'
            return result

        # Create combined index file
        combined_index = task_dir / "combined_cdr3_index.ndx"
        with open(combined_index, 'w') as outf:
            with open(phla_index, 'r') as inf:
                outf.write(inf.read())
            with open(cdr3_index, 'r') as inf:
                outf.write(inf.read())

        # Output RMSD file
        rmsd_file = task_dir / "cdr3_beta_rmsd.xvg"

        # Run gmx rms: align to TCR, calculate CDR3_beta_CA
        cmd = [
            'gmx', 'rms',
            '-s', str(topology),
            '-f', str(trajectory),
            '-n', str(combined_index),
            '-o', str(rmsd_file),
            '-tu', 'ns'
        ]

        # Provide selections: TCR (group 1) for alignment, CDR3_TCR_beta_CA (group 5) for RMSD
        process = subprocess.Popen(
            cmd,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )

        stdout, stderr = process.communicate(input="TCR\nCDR3_TCR_beta_CA\n", timeout=120)

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

        rmsd_data = np.array(rmsd_data)
        time_data = rmsd_data[:, 0]
        rmsd_values = rmsd_data[:, 1]

        # Calculate statistics
        result.update({
            'status': 'success',
            'n_frames': len(rmsd_values),
            'time_start_ns': float(time_data[0]),
            'time_end_ns': float(time_data[-1]),
            'rmsd_mean_nm': float(np.mean(rmsd_values)),
            'rmsd_std_nm': float(np.std(rmsd_values)),
            'rmsd_min_nm': float(np.min(rmsd_values)),
            'rmsd_max_nm': float(np.max(rmsd_values)),
            'rmsd_median_nm': float(np.median(rmsd_values)),
            'rmsd_range_nm': float(np.max(rmsd_values) - np.min(rmsd_values)),
            'rmsd_file': str(rmsd_file)
        })

        # Save per-task statistics
        stats_file = cdr3_analysis_dir / task_name / "cdr3_beta_rmsd_stats.json"
        with open(stats_file, 'w') as f:
            json.dump(result, f, indent=2)

        logger.info(f"✓ {task_name}: RMSD mean={result['rmsd_mean_nm']:.4f} nm")

        return result

    except Exception as e:
        result['error'] = str(e)
        logger.error(f"✗ {task_name}: {e}")
        return result


def main():
    """Main function for batch CDR3β RMSD calculation."""

    logger.info("="*60)
    logger.info("Batch CDR3β RMSD Calculation (TCR Alignment)")
    logger.info("="*60)

    # Paths
    base_dir = Path("/home/xumy/work/development/Immunex/input/pbc_1000frames_2step")
    cdr3_analysis_dir = Path("/home/xumy/work/development/Immunex/output/cdr3_analysis")

    # Find all tasks with required files
    tasks = []
    for task_dir in sorted(base_dir.iterdir()):
        if not task_dir.is_dir():
            continue

        task_name = task_dir.name

        # Check if CDR3 analysis exists
        if not (cdr3_analysis_dir / task_name).exists():
            logger.debug(f"Skipping {task_name}: no CDR3 analysis")
            continue

        # Check if trajectory exists
        trajectory = task_dir / f"{task_name}_2step_processed.xtc"
        if not trajectory.exists():
            logger.debug(f"Skipping {task_name}: no trajectory")
            continue

        tasks.append((task_name, task_dir, cdr3_analysis_dir))

    logger.info(f"Found {len(tasks)} tasks with complete data")

    if not tasks:
        logger.error("No tasks found for RMSD calculation")
        return

    # Process tasks in parallel
    logger.info(f"\nProcessing {len(tasks)} tasks with 8 workers...")

    with ProcessPoolExecutor(max_workers=8) as executor:
        results = list(executor.map(process_single_task, tasks))

    # Summarize results
    successful = [r for r in results if r['status'] == 'success']
    failed = [r for r in results if r['status'] == 'failed']

    logger.info("\n" + "="*60)
    logger.info("Batch RMSD Calculation Summary")
    logger.info("="*60)
    logger.info(f"Total tasks: {len(tasks)}")
    logger.info(f"Successful: {len(successful)}")
    logger.info(f"Failed: {len(failed)}")
    logger.info(f"Success rate: {len(successful)/len(tasks)*100:.1f}%")

    # Save results to CSV
    if successful:
        df = pd.DataFrame(successful)
        csv_file = cdr3_analysis_dir / "cdr3_beta_rmsd_summary.csv"
        df.to_csv(csv_file, index=False)
        logger.info(f"\nResults saved: {csv_file}")

        # Print statistics
        logger.info("\nRMSD Statistics (nm):")
        logger.info(f"  Mean of means: {df['rmsd_mean_nm'].mean():.4f} ± {df['rmsd_mean_nm'].std():.4f}")
        logger.info(f"  Min: {df['rmsd_min_nm'].min():.4f}")
        logger.info(f"  Max: {df['rmsd_max_nm'].max():.4f}")

    # Save failed tasks
    if failed:
        failed_file = cdr3_analysis_dir / "cdr3_beta_rmsd_failed.log"
        with open(failed_file, 'w') as f:
            f.write(f"Failed CDR3β RMSD Tasks: {len(failed)}\n")
            f.write("="*60 + "\n\n")
            for task in failed:
                f.write(f"Task: {task['task_name']}\n")
                f.write(f"Error: {task['error']}\n")
                f.write("-"*60 + "\n\n")
        logger.info(f"\nFailed tasks log: {failed_file}")

    logger.info("="*60)


if __name__ == "__main__":
    main()
