#!/usr/bin/env python3
"""
Regenerate pHLA-TCR index files for failed tasks using standardized PDB files.
"""

import os
import sys
import subprocess
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import logging

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def generate_phla_tcr_index(task_info):
    """Generate pHLA-TCR index file for a single task."""
    task_name, pdb_file, task_dir = task_info

    result = {
        'task_name': task_name,
        'status': 'failed',
        'error': None
    }

    try:
        if not pdb_file.exists():
            result['error'] = f'PDB file not found: {pdb_file}'
            return result

        # Create index file using gmx make_ndx
        index_file = task_dir / "phla_tcr.ndx"

        # Input commands for gmx make_ndx
        # chain A B C = pHLA
        # chain D E = TCR
        commands = "chain A B C\nname 20 pHLA\nchain D E\nname 21 TCR\nq\n"

        cmd = ['gmx', 'make_ndx', '-f', str(pdb_file), '-o', str(index_file)]

        process = subprocess.Popen(
            cmd,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )

        stdout, stderr = process.communicate(input=commands, timeout=30)

        if process.returncode != 0:
            result['error'] = f'gmx make_ndx failed: {stderr}'
            return result

        # Verify index file was created
        if not index_file.exists():
            result['error'] = 'Index file not created'
            return result

        # Check if pHLA and TCR groups are present
        with open(index_file, 'r') as f:
            content = f.read()
            if 'pHLA' not in content or 'TCR' not in content:
                result['error'] = 'pHLA or TCR group not found in index'
                return result

        result['status'] = 'success'
        result['index_file'] = str(index_file)
        logger.info(f"✓ {task_name}: Index file created")

        return result

    except Exception as e:
        result['error'] = str(e)
        logger.error(f"✗ {task_name}: {e}")
        return result


def main():
    """Main function to regenerate index files."""

    logger.info("=" * 60)
    logger.info("Regenerating pHLA-TCR Index Files")
    logger.info("=" * 60)

    # Paths
    protein_pdb_dir = Path("/home/xumy/work/development/Immunex/input/protein_only_pdbs")
    trajectory_base = Path("/home/xumy/work/development/Immunex/input/pbc_1000frames_2step")

    # Get failed tasks
    failed_log = Path("/home/xumy/work/development/Immunex/output/cdr3_analysis/cdr3_beta_rmsd_failed.log")

    failed_tasks = []
    with open(failed_log) as f:
        for line in f:
            if line.startswith("Task: "):
                task_name = line.strip().split(": ")[1]
                failed_tasks.append(task_name)

    logger.info(f"Found {len(failed_tasks)} failed tasks")

    # Prepare task info
    tasks_to_process = []
    for task_name in failed_tasks:
        pdb_id = task_name.replace('_1', '').replace('_run1', '').lower()
        pdb_file = protein_pdb_dir / f"{pdb_id}_protein_only.pdb"
        task_dir = trajectory_base / task_name

        if pdb_file.exists() and task_dir.exists():
            tasks_to_process.append((task_name, pdb_file, task_dir))
        else:
            logger.warning(f"Skipping {task_name}: Protein PDB or task dir not found")

    logger.info(f"Processing {len(tasks_to_process)} tasks with available PDB files")

    # Process tasks in parallel
    logger.info("\nGenerating index files...")

    with ProcessPoolExecutor(max_workers=8) as executor:
        results = list(executor.map(generate_phla_tcr_index, tasks_to_process))

    # Summarize results
    successful = [r for r in results if r['status'] == 'success']
    failed = [r for r in results if r['status'] == 'failed']

    logger.info("\n" + "=" * 60)
    logger.info("Index Generation Summary")
    logger.info("=" * 60)
    logger.info(f"Total tasks: {len(tasks_to_process)}")
    logger.info(f"Successful: {len(successful)}")
    logger.info(f"Failed: {len(failed)}")
    logger.info(f"Success rate: {len(successful)/len(tasks_to_process)*100:.1f}%")

    if failed:
        logger.info("\nFailed tasks:")
        for r in failed[:10]:
            logger.info(f"  {r['task_name']}: {r['error']}")

    logger.info("=" * 60)


if __name__ == "__main__":
    main()
