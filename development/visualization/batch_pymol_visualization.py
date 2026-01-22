#!/usr/bin/env python3
"""Batch PyMOL visualization generation for clustering results"""

import sys
from pathlib import Path
import argparse
import logging
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
import traceback

sys.path.insert(0, str(Path(__file__).parent.parent))

# Import the visualization generator functions
from scripts.generate_pymol_cluster_visualization import (
    extract_protein_pdb,
    generate_pymol_script
)

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('batch_pymol_viz.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)


def process_single_visualization(task_info):
    """Generate PyMOL visualization for a single task

    Args:
        task_info: Dict with task information

    Returns:
        Dict with visualization status
    """
    task_name = task_info['name']
    start_time = time.time()

    try:
        logger.info(f"[{task_name}] Generating PyMOL visualization...")

        # Output directory
        output_dir = Path(task_info['output_dir'])
        output_dir.mkdir(parents=True, exist_ok=True)

        # Extract PDB (frame 0)
        pdb_file = output_dir / "protein_structure.pdb"
        extract_protein_pdb(
            topology=task_info['topology'],
            trajectory=task_info['trajectory'],
            output_pdb=pdb_file,
            frame_index=0
        )

        # Generate PyMOL script
        pymol_script = output_dir / "cluster_visualization.pml"
        generate_pymol_script(
            cluster_assignments_csv=task_info['cluster_csv'],
            pdb_file=pdb_file,
            output_script=pymol_script,
            topology=task_info['topology'],
            top_n_clusters=15,
            min_cluster_size=5
        )

        elapsed = time.time() - start_time
        logger.info(f"[{task_name}] Completed in {elapsed:.2f}s")

        return {
            'task': task_name,
            'status': 'success',
            'elapsed': elapsed,
            'output_dir': str(output_dir)
        }

    except Exception as e:
        elapsed = time.time() - start_time
        error_msg = f"{str(e)}\n{traceback.format_exc()}"
        logger.error(f"[{task_name}] Failed after {elapsed:.2f}s: {e}")

        return {
            'task': task_name,
            'status': 'failed',
            'elapsed': elapsed,
            'error': error_msg
        }


def discover_visualization_tasks(input_dir, trajectory_dir):
    """Discover tasks that need PyMOL visualization

    Args:
        input_dir: Allostery analysis output directory
        trajectory_dir: Directory containing trajectories

    Returns:
        List of task info dictionaries
    """
    input_path = Path(input_dir)
    traj_path = Path(trajectory_dir)
    tasks = []

    for task_dir in sorted(input_path.iterdir()):
        if not task_dir.is_dir():
            continue

        task_name = task_dir.name
        cluster_csv = task_dir / "clustering" / "cluster_assignments.csv"

        # Check if clustering results exist
        if not cluster_csv.exists():
            logger.warning(f"[{task_name}] Missing clustering results, skipping")
            continue

        # Find trajectory files
        task_traj_dir = traj_path / task_name
        if not task_traj_dir.exists():
            logger.warning(f"[{task_name}] Trajectory directory not found, skipping")
            continue

        topology = task_traj_dir / "md.tpr"

        # Find trajectory file
        trajectory = None
        for pattern in ["*_processed.xtc", "*_2step_processed.xtc", "md.xtc"]:
            traj_files = list(task_traj_dir.glob(pattern))
            if traj_files:
                trajectory = traj_files[0]
                break

        if not topology.exists() or not trajectory:
            logger.warning(f"[{task_name}] Missing topology or trajectory, skipping")
            continue

        tasks.append({
            'name': task_name,
            'topology': str(topology),
            'trajectory': str(trajectory),
            'cluster_csv': str(cluster_csv),
            'output_dir': str(task_dir / "pymol_viz")
        })

    return tasks


def main():
    parser = argparse.ArgumentParser(
        description='Batch PyMOL visualization generation'
    )

    parser.add_argument(
        '--input-dir',
        default='./output/allostery_analysis',
        help='Input directory with clustering results'
    )
    parser.add_argument(
        '--trajectory-dir',
        default='./input/pbc_1000frames_2step',
        help='Directory containing trajectory files'
    )
    parser.add_argument(
        '--max-workers',
        type=int,
        default=8,
        help='Number of parallel workers (default: 8)'
    )
    parser.add_argument(
        '--limit',
        type=int,
        help='Limit number of tasks (for testing)'
    )

    args = parser.parse_args()

    logger.info("="*80)
    logger.info("Batch PyMOL Visualization Generation")
    logger.info("="*80)
    logger.info(f"Input directory: {args.input_dir}")
    logger.info(f"Trajectory directory: {args.trajectory_dir}")
    logger.info(f"Max workers: {args.max_workers}")

    # Discover tasks
    logger.info("\nDiscovering tasks...")
    tasks = discover_visualization_tasks(args.input_dir, args.trajectory_dir)

    if args.limit:
        tasks = tasks[:args.limit]
        logger.info(f"Limited to first {args.limit} tasks")

    logger.info(f"Found {len(tasks)} tasks to visualize")

    if len(tasks) == 0:
        logger.error("No valid tasks found!")
        return 1

    # Process in parallel
    logger.info(f"\nStarting parallel processing with {args.max_workers} workers...")
    logger.info(f"Estimated time: {len(tasks) * 5 / args.max_workers / 60:.1f} minutes")
    logger.info("="*80)

    start_time = time.time()
    results = []
    completed = 0
    failed = 0

    with ProcessPoolExecutor(max_workers=args.max_workers) as executor:
        future_to_task = {
            executor.submit(process_single_visualization, task): task
            for task in tasks
        }

        for future in as_completed(future_to_task):
            task = future_to_task[future]
            try:
                result = future.result()
                results.append(result)

                if result['status'] == 'success':
                    completed += 1
                else:
                    failed += 1

                total_done = completed + failed
                progress = total_done / len(tasks) * 100
                logger.info(
                    f"\nProgress: {total_done}/{len(tasks)} ({progress:.1f}%) | "
                    f"Success: {completed} | Failed: {failed}"
                )

            except Exception as e:
                logger.error(f"Task {task['name']} exception: {e}")
                failed += 1

    # Summary
    total_elapsed = time.time() - start_time
    logger.info("\n" + "="*80)
    logger.info("BATCH VISUALIZATION COMPLETED")
    logger.info("="*80)
    logger.info(f"Total time: {total_elapsed/60:.2f} minutes")
    logger.info(f"Total tasks: {len(tasks)}")
    logger.info(f"Successful: {completed}")
    logger.info(f"Failed: {failed}")
    logger.info(f"Success rate: {completed/len(tasks)*100:.1f}%")

    # Write summary
    summary_file = Path(args.input_dir) / "pymol_viz_summary.txt"
    with open(summary_file, 'w') as f:
        f.write("Batch PyMOL Visualization Summary\n")
        f.write("="*80 + "\n\n")
        f.write(f"Total time: {total_elapsed/60:.2f} minutes\n")
        f.write(f"Total tasks: {len(tasks)}\n")
        f.write(f"Successful: {completed}\n")
        f.write(f"Failed: {failed}\n\n")

        f.write("Successful Tasks:\n")
        f.write("-"*80 + "\n")
        for r in results:
            if r['status'] == 'success':
                f.write(f"{r['task']}: {r['output_dir']}\n")

        if failed > 0:
            f.write("\nFailed Tasks:\n")
            f.write("-"*80 + "\n")
            for r in results:
                if r['status'] == 'failed':
                    f.write(f"{r['task']}: {r['error'][:200]}\n")

    logger.info(f"\nSummary saved: {summary_file}")
    logger.info("="*80)

    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
