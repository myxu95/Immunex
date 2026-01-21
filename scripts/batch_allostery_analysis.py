"""
Batch Allostery Analysis for TCR-pMHC Complexes

This script performs contact correlation analysis on multiple MD trajectories
to identify allosteric communication pathways in TCR-pMHC complexes.
"""

import sys
from pathlib import Path
import logging
import json
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Dict, List

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from aftermd.analysis.allostery import ContactCorrelationAnalyzer


def analyze_single_task(task_info: Dict) -> Dict:
    """
    Analyze a single MD task for allosteric correlations

    Args:
        task_info: Dictionary with task information
            - name: Task name (e.g., "1ao7_run1")
            - topology: Path to md.tpr
            - trajectory: Path to processed trajectory
            - output_dir: Output directory for results

    Returns:
        Dictionary with analysis results and status
    """
    task_name = task_info['name']
    logger.info(f"Starting analysis for {task_name}...")

    try:
        # Initialize analyzer
        analyzer = ContactCorrelationAnalyzer(
            topology=task_info['topology'],
            trajectory=task_info['trajectory']
        )

        # Run analysis with optimized parameters for TCR-pMHC
        results = analyzer.run_full_analysis(
            selection="protein",           # Analyze entire complex
            output_dir=task_info['output_dir'],
            cutoff=4.5,                   # Standard contact distance
            seq_dist_cutoff=3,            # Exclude sequence neighbors
            min_frequency=0.15,           # More stringent for 1000 frames
            plot_heatmap=True,
            heatmap_figsize=(16, 14)
        )

        # Add task name to results
        results['task_name'] = task_name
        results['status'] = 'success'

        logger.info(f"Completed {task_name}: {results['n_contacts']} contacts, "
                   f"mean_corr={results['mean_correlation']:.3f}")

        return results

    except Exception as e:
        logger.error(f"Failed to analyze {task_name}: {e}")
        return {
            'task_name': task_name,
            'status': 'failed',
            'error': str(e)
        }


def collect_tasks(input_dir: str, max_tasks: int = None) -> List[Dict]:
    """
    Collect all tasks from input directory

    Args:
        input_dir: Base directory containing task subdirectories
        max_tasks: Maximum number of tasks to process (None = all)

    Returns:
        List of task information dictionaries
    """
    input_path = Path(input_dir)

    if not input_path.exists():
        raise ValueError(f"Input directory not found: {input_dir}")

    tasks = []

    # Iterate through task directories
    for task_dir in sorted(input_path.iterdir()):
        if not task_dir.is_dir():
            continue

        task_name = task_dir.name

        # Find required files
        topology_file = task_dir / "md.tpr"
        trajectory_file = task_dir / f"{task_name}_2step_processed.xtc"

        if not topology_file.exists():
            logger.warning(f"Skipping {task_name}: md.tpr not found")
            continue

        if not trajectory_file.exists():
            logger.warning(f"Skipping {task_name}: processed trajectory not found")
            continue

        # Create output directory
        output_dir = Path("output/allostery_analysis") / task_name

        tasks.append({
            'name': task_name,
            'topology': str(topology_file),
            'trajectory': str(trajectory_file),
            'output_dir': str(output_dir)
        })

        if max_tasks and len(tasks) >= max_tasks:
            break

    logger.info(f"Found {len(tasks)} tasks to process")
    return tasks


def batch_analysis(input_dir: str,
                   output_summary: str,
                   max_tasks: int = None,
                   max_workers: int = 4) -> Dict:
    """
    Run batch allostery analysis on multiple tasks

    Args:
        input_dir: Directory containing task subdirectories
        output_summary: Path to save summary JSON
        max_tasks: Maximum number of tasks (None = all)
        max_workers: Number of parallel workers

    Returns:
        Summary dictionary with all results
    """
    logger.info("=" * 80)
    logger.info("BATCH ALLOSTERY ANALYSIS")
    logger.info("=" * 80)

    # Collect tasks
    tasks = collect_tasks(input_dir, max_tasks)

    if len(tasks) == 0:
        logger.error("No valid tasks found!")
        return {'status': 'failed', 'reason': 'no_tasks'}

    # Run analysis in parallel
    results = []
    failed_tasks = []

    logger.info(f"Processing {len(tasks)} tasks with {max_workers} workers...")

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        # Submit all tasks
        future_to_task = {
            executor.submit(analyze_single_task, task): task['name']
            for task in tasks
        }

        # Collect results as they complete
        for idx, future in enumerate(as_completed(future_to_task), 1):
            task_name = future_to_task[future]

            try:
                result = future.result()
                results.append(result)

                if result['status'] == 'failed':
                    failed_tasks.append(task_name)

                logger.info(f"Progress: {idx}/{len(tasks)} tasks completed")

            except Exception as e:
                logger.error(f"Exception for task {task_name}: {e}")
                results.append({
                    'task_name': task_name,
                    'status': 'failed',
                    'error': str(e)
                })
                failed_tasks.append(task_name)

    # Compile summary
    successful_results = [r for r in results if r['status'] == 'success']

    summary = {
        'total_tasks': len(tasks),
        'successful': len(successful_results),
        'failed': len(failed_tasks),
        'failed_tasks': failed_tasks,
        'results': results
    }

    # Calculate statistics
    if successful_results:
        summary['statistics'] = {
            'mean_contacts': sum(r['n_contacts'] for r in successful_results) / len(successful_results),
            'mean_correlation': sum(r['mean_correlation'] for r in successful_results) / len(successful_results),
            'total_high_pos_corr': sum(r['n_high_pos_correlation'] for r in successful_results),
            'total_high_neg_corr': sum(r['n_high_neg_correlation'] for r in successful_results)
        }

    # Save summary
    output_path = Path(output_summary)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, 'w') as f:
        json.dump(summary, f, indent=2)

    logger.info("=" * 80)
    logger.info("BATCH ANALYSIS COMPLETED")
    logger.info("=" * 80)
    logger.info(f"Total tasks: {summary['total_tasks']}")
    logger.info(f"Successful: {summary['successful']}")
    logger.info(f"Failed: {summary['failed']}")

    if successful_results:
        logger.info(f"Average contacts per task: {summary['statistics']['mean_contacts']:.0f}")
        logger.info(f"Average correlation: {summary['statistics']['mean_correlation']:.3f}")

    logger.info(f"Summary saved to: {output_summary}")

    return summary


def test_single_task():
    """Test analysis on a single task"""
    logger.info("=" * 80)
    logger.info("TESTING SINGLE TASK")
    logger.info("=" * 80)

    input_dir = "/home/xumy/work/development/AfterMD/input/pbc_1000frames_2step"

    # Collect first task
    tasks = collect_tasks(input_dir, max_tasks=1)

    if not tasks:
        logger.error("No tasks found!")
        return

    task = tasks[0]
    logger.info(f"Testing task: {task['name']}")

    # Run analysis
    result = analyze_single_task(task)

    if result['status'] == 'success':
        logger.info("Test successful!")
        logger.info(f"Results:")
        logger.info(f"  Persistent contacts: {result['n_contacts']}")
        logger.info(f"  Mean correlation: {result['mean_correlation']:.3f}")
        logger.info(f"  High positive correlations (>0.7): {result['n_high_pos_correlation']}")
        logger.info(f"  High negative correlations (<-0.7): {result['n_high_neg_correlation']}")
        logger.info(f"  Output directory: {result['files']['correlation_heatmap']}")
    else:
        logger.error(f"Test failed: {result.get('error', 'Unknown error')}")


def main():
    """Main entry point"""
    import argparse

    parser = argparse.ArgumentParser(
        description="Batch allostery analysis for TCR-pMHC complexes"
    )
    parser.add_argument(
        '--input-dir',
        default='/home/xumy/work/development/AfterMD/input/pbc_1000frames_2step',
        help='Input directory containing task subdirectories'
    )
    parser.add_argument(
        '--output-summary',
        default='output/allostery_analysis/batch_summary.json',
        help='Path to save summary JSON'
    )
    parser.add_argument(
        '--max-tasks',
        type=int,
        default=None,
        help='Maximum number of tasks to process (default: all)'
    )
    parser.add_argument(
        '--max-workers',
        type=int,
        default=4,
        help='Number of parallel workers (default: 4)'
    )
    parser.add_argument(
        '--test',
        action='store_true',
        help='Test on a single task only'
    )

    args = parser.parse_args()

    if args.test:
        test_single_task()
    else:
        batch_analysis(
            input_dir=args.input_dir,
            output_summary=args.output_summary,
            max_tasks=args.max_tasks,
            max_workers=args.max_workers
        )


if __name__ == "__main__":
    main()
