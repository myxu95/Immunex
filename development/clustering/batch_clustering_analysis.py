#!/usr/bin/env python3
"""Batch clustering analysis for contact correlation matrices"""

import sys
from pathlib import Path
import argparse
import logging
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
import traceback
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent))

from aftermd.analysis.allostery import ContactCorrelationAnalyzer

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('batch_clustering.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)


def process_single_clustering(task_info):
    """Process clustering for a single task

    Args:
        task_info: Dict with 'name' and 'result_dir'

    Returns:
        Dict with clustering status and results
    """
    task_name = task_info['name']
    result_dir = Path(task_info['result_dir'])
    start_time = time.time()

    try:
        logger.info(f"[{task_name}] Starting clustering...")

        # Load correlation matrix and contact labels
        npz_file = result_dir / "correlation_matrix.npz"
        labels_file = result_dir / "contact_labels.csv"

        if not npz_file.exists():
            raise FileNotFoundError(f"Correlation matrix not found: {npz_file}")
        if not labels_file.exists():
            raise FileNotFoundError(f"Contact labels not found: {labels_file}")

        # Load data
        data = np.load(npz_file)
        correlation_matrix = data['correlation_matrix']

        labels_df = pd.read_csv(labels_file)
        contact_labels = list(labels_df[['ResA_ID', 'ResA_Name', 'ResB_ID', 'ResB_Name']].itertuples(index=False, name=None))

        logger.info(f"[{task_name}] Loaded {correlation_matrix.shape[0]} contacts")

        # Create temporary analyzer instance (just for clustering methods)
        analyzer = ContactCorrelationAnalyzer.__new__(ContactCorrelationAnalyzer)
        analyzer.correlation_matrix = correlation_matrix

        # Perform hierarchical clustering
        cluster_results = analyzer.cluster_hierarchical(
            correlation_matrix=correlation_matrix,
            method=task_info['linkage_method'],
            threshold=task_info['threshold'],
            min_cluster_size=task_info['min_cluster_size']
        )

        # Create output directory
        output_dir = result_dir / "clustering"
        output_dir.mkdir(exist_ok=True)

        # Visualize dendrogram
        analyzer.plot_cluster_dendrogram(
            linkage_matrix=cluster_results['linkage_matrix'],
            output_file=str(output_dir / "dendrogram.png"),
            threshold=task_info['threshold']
        )

        # Save clustering results
        saved_files = analyzer.save_cluster_results(
            cluster_labels=cluster_results['cluster_labels'],
            contact_labels=contact_labels,
            correlation_matrix=correlation_matrix,
            output_dir=str(output_dir)
        )

        elapsed = time.time() - start_time
        logger.info(f"[{task_name}] Completed in {elapsed:.2f}s - {cluster_results['n_clusters']} clusters")

        return {
            'task': task_name,
            'status': 'success',
            'elapsed': elapsed,
            'n_clusters': cluster_results['n_clusters'],
            'n_contacts': correlation_matrix.shape[0],
            'n_noise': int(np.sum(cluster_results['cluster_labels'] < 0)),
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


def discover_completed_tasks(output_dir):
    """Discover all completed allostery analysis tasks

    Args:
        output_dir: Path to allostery_analysis output directory

    Returns:
        List of task info dictionaries
    """
    output_path = Path(output_dir)
    tasks = []

    for task_dir in sorted(output_path.iterdir()):
        if not task_dir.is_dir():
            continue

        task_name = task_dir.name
        npz_file = task_dir / "correlation_matrix.npz"
        labels_file = task_dir / "contact_labels.csv"

        # Check if correlation analysis is complete
        if not npz_file.exists() or not labels_file.exists():
            logger.warning(f"[{task_name}] Missing correlation results, skipping")
            continue

        tasks.append({
            'name': task_name,
            'result_dir': str(task_dir)
        })

    return tasks


def main():
    parser = argparse.ArgumentParser(
        description='Batch clustering analysis for contact correlation matrices'
    )

    parser.add_argument(
        '--input-dir',
        default='./output/allostery_analysis',
        help='Input directory containing allostery analysis results'
    )
    parser.add_argument(
        '--method',
        choices=['hierarchical', 'network'],
        default='hierarchical',
        help='Clustering method (default: hierarchical)'
    )
    parser.add_argument(
        '--threshold',
        type=float,
        default=0.7,
        help='Distance/correlation threshold (default: 0.7)'
    )
    parser.add_argument(
        '--min-cluster-size',
        type=int,
        default=5,
        help='Minimum cluster size (default: 5)'
    )
    parser.add_argument(
        '--linkage-method',
        choices=['single', 'complete', 'average', 'ward'],
        default='average',
        help='Linkage method for hierarchical clustering (default: average)'
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
        help='Limit number of tasks to process (for testing)'
    )

    args = parser.parse_args()

    logger.info("="*80)
    logger.info("Batch Clustering Analysis - Parallel Processing")
    logger.info("="*80)
    logger.info(f"Input directory: {args.input_dir}")
    logger.info(f"Clustering method: {args.method}")
    logger.info(f"Threshold: {args.threshold}")
    logger.info(f"Min cluster size: {args.min_cluster_size}")
    logger.info(f"Linkage method: {args.linkage_method}")
    logger.info(f"Max workers: {args.max_workers}")

    # Discover tasks
    logger.info("\nDiscovering completed tasks...")
    tasks = discover_completed_tasks(args.input_dir)

    if args.limit:
        tasks = tasks[:args.limit]
        logger.info(f"Limited to first {args.limit} tasks")

    logger.info(f"Found {len(tasks)} tasks to cluster")

    if len(tasks) == 0:
        logger.error("No valid tasks found!")
        return 1

    # Add clustering parameters to each task
    for task in tasks:
        task['threshold'] = args.threshold
        task['min_cluster_size'] = args.min_cluster_size
        task['linkage_method'] = args.linkage_method

    # Process tasks in parallel
    logger.info(f"\nStarting parallel clustering with {args.max_workers} workers...")
    logger.info(f"Estimated total time: {len(tasks) * 2 / args.max_workers / 60:.1f} minutes")
    logger.info("="*80)

    start_time = time.time()
    results = []
    completed = 0
    failed = 0

    with ProcessPoolExecutor(max_workers=args.max_workers) as executor:
        # Submit all tasks
        future_to_task = {
            executor.submit(process_single_clustering, task): task
            for task in tasks
        }

        # Process completed tasks
        for future in as_completed(future_to_task):
            task = future_to_task[future]
            try:
                result = future.result()
                results.append(result)

                if result['status'] == 'success':
                    completed += 1
                else:
                    failed += 1

                # Print progress
                total_done = completed + failed
                progress = total_done / len(tasks) * 100
                logger.info(
                    f"\nProgress: {total_done}/{len(tasks)} ({progress:.1f}%) | "
                    f"Success: {completed} | Failed: {failed}"
                )

            except Exception as e:
                logger.error(f"Task {task['name']} generated an exception: {e}")
                failed += 1

    # Summary
    total_elapsed = time.time() - start_time
    logger.info("\n" + "="*80)
    logger.info("BATCH CLUSTERING COMPLETED")
    logger.info("="*80)
    logger.info(f"Total time: {total_elapsed/60:.2f} minutes ({total_elapsed/3600:.2f} hours)")
    logger.info(f"Total tasks: {len(tasks)}")
    logger.info(f"Successful: {completed}")
    logger.info(f"Failed: {failed}")
    logger.info(f"Success rate: {completed/len(tasks)*100:.1f}%")

    # Statistics
    if completed > 0:
        success_results = [r for r in results if r['status'] == 'success']
        avg_clusters = np.mean([r['n_clusters'] for r in success_results])
        avg_contacts = np.mean([r['n_contacts'] for r in success_results])
        avg_noise = np.mean([r['n_noise'] for r in success_results])

        logger.info("\nClustering Statistics:")
        logger.info(f"  Average clusters per task: {avg_clusters:.1f}")
        logger.info(f"  Average contacts per task: {avg_contacts:.1f}")
        logger.info(f"  Average noise contacts: {avg_noise:.1f} ({avg_noise/avg_contacts*100:.1f}%)")

    # Write summary report
    summary_file = Path(args.input_dir) / "clustering_summary.txt"
    with open(summary_file, 'w') as f:
        f.write("Batch Clustering Analysis Summary\n")
        f.write("="*80 + "\n\n")
        f.write(f"Total time: {total_elapsed/60:.2f} minutes\n")
        f.write(f"Total tasks: {len(tasks)}\n")
        f.write(f"Successful: {completed}\n")
        f.write(f"Failed: {failed}\n")
        f.write(f"Success rate: {completed/len(tasks)*100:.1f}%\n\n")

        f.write("Parameters:\n")
        f.write(f"  Method: {args.method}\n")
        f.write(f"  Threshold: {args.threshold}\n")
        f.write(f"  Min cluster size: {args.min_cluster_size}\n")
        f.write(f"  Linkage method: {args.linkage_method}\n\n")

        f.write("Successful Tasks:\n")
        f.write("-"*80 + "\n")
        for r in results:
            if r['status'] == 'success':
                f.write(f"{r['task']}: {r['n_clusters']} clusters, "
                       f"{r['n_contacts']} contacts, {r['n_noise']} noise, "
                       f"{r['elapsed']:.2f}s\n")

        f.write("\nFailed Tasks:\n")
        f.write("-"*80 + "\n")
        for r in results:
            if r['status'] == 'failed':
                f.write(f"{r['task']}: {r['error'][:200]}...\n")

    logger.info(f"\nSummary report saved: {summary_file}")
    logger.info("="*80)

    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
