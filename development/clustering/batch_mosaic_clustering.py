#!/usr/bin/env python3
"""Batch MoSAIC clustering for contact correlation matrices

This script performs MoSAIC clustering on all correlation matrices
in the allostery_analysis output directory.

Usage:
    python batch_mosaic_clustering.py \\
        --input output/allostery_analysis \\
        --mode CPM \\
        --resolution 0.75 \\
        --max-workers 4
"""

import sys
from pathlib import Path
import argparse
import logging
import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed
import time

sys.path.insert(0, str(Path(__file__).parent.parent))

from aftermd.analysis.allostery import ContactCorrelationAnalyzer

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def process_single_task(task_info):
    """Process a single task's MoSAIC clustering

    Args:
        task_info: Dictionary with task_name, input_dir, output_dir, and parameters

    Returns:
        Dictionary with task_name, status, and results
    """
    task_name = task_info['task_name']
    input_dir = task_info['input_dir']
    output_dir = task_info['output_dir']
    mode = task_info['mode']
    resolution_parameter = task_info['resolution_parameter']
    weighted = task_info['weighted']

    try:
        logger.info(f"Processing {task_name}...")

        # Load correlation matrix and contact labels
        npz_file = Path(input_dir) / "correlation_matrix.npz"
        labels_file = Path(input_dir) / "contact_labels.csv"

        if not npz_file.exists():
            return {
                'task_name': task_name,
                'status': 'failed',
                'error': 'correlation_matrix.npz not found'
            }

        data = np.load(npz_file)
        correlation_matrix = data['correlation_matrix']

        labels_df = pd.read_csv(labels_file)
        contact_labels = list(labels_df[['ResA_ID', 'ResA_Name', 'ResB_ID', 'ResB_Name']].itertuples(index=False, name=None))

        # Create temporary analyzer instance
        analyzer = ContactCorrelationAnalyzer.__new__(ContactCorrelationAnalyzer)

        # Perform MoSAIC clustering
        mosaic_result = analyzer.cluster_mosaic(
            correlation_matrix=correlation_matrix,
            mode=mode,
            resolution_parameter=resolution_parameter,
            weighted=weighted,
            seed=42
        )

        # Save results
        Path(output_dir).mkdir(parents=True, exist_ok=True)

        saved_files = analyzer.save_mosaic_results(
            mosaic_result=mosaic_result,
            contact_labels=contact_labels,
            correlation_matrix=correlation_matrix,
            output_dir=output_dir
        )

        # Plot block-diagonal matrix
        analyzer.plot_mosaic_matrix(
            matrix_reordered=mosaic_result['matrix_reordered'],
            ticks=mosaic_result['ticks'],
            output_file=str(Path(output_dir) / "mosaic_matrix.png"),
            figsize=(14, 12),
            title=f'{task_name} - MoSAIC Clustered Matrix'
        )

        # Count meaningful clusters (size > 1)
        cluster_sizes = [len(c) for c in mosaic_result['mosaic_clusters']]
        meaningful_clusters = sum(1 for s in cluster_sizes if s > 1)

        return {
            'task_name': task_name,
            'status': 'success',
            'n_clusters': mosaic_result['n_clusters'],
            'n_meaningful_clusters': meaningful_clusters,
            'n_contacts': len(contact_labels),
            'silhouette_score': mosaic_result['silhouette_score'],
            'resolution_parameter': mosaic_result['resolution_parameter']
        }

    except Exception as e:
        logger.error(f"Error processing {task_name}: {e}")
        return {
            'task_name': task_name,
            'status': 'failed',
            'error': str(e)
        }


def discover_tasks(base_dir, force_rerun=False):
    """Discover tasks with correlation matrices

    Args:
        base_dir: Base directory containing task subdirectories
        force_rerun: If True, rerun even if mosaic results exist

    Returns:
        List of task directories to process
    """
    base_path = Path(base_dir)

    tasks = []
    for task_dir in sorted(base_path.iterdir()):
        if not task_dir.is_dir():
            continue

        corr_matrix_file = task_dir / "correlation_matrix.npz"
        contact_labels_file = task_dir / "contact_labels.csv"
        mosaic_dir = task_dir / "mosaic"

        if not corr_matrix_file.exists() or not contact_labels_file.exists():
            continue

        # Check if already processed
        if not force_rerun and mosaic_dir.exists():
            mosaic_summary = mosaic_dir / "mosaic_cluster_summary.csv"
            if mosaic_summary.exists():
                logger.info(f"Skipping {task_dir.name} (already processed)")
                continue

        tasks.append(task_dir)

    return tasks


def main():
    parser = argparse.ArgumentParser(
        description='Batch MoSAIC clustering for contact correlation matrices',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        '--input', required=True,
        help='Input directory containing task subdirectories'
    )
    parser.add_argument(
        '--mode', choices=['CPM', 'modularity', 'linkage'],
        default='CPM',
        help='Clustering mode (default: CPM)'
    )
    parser.add_argument(
        '--resolution', type=float,
        help='Resolution parameter (default: 3rd quartile of matrix values)'
    )
    parser.add_argument(
        '--weighted', action='store_true', default=True,
        help='Use weighted edges (default: True)'
    )
    parser.add_argument(
        '--max-workers', type=int, default=4,
        help='Maximum parallel workers (default: 4)'
    )
    parser.add_argument(
        '--force-rerun', action='store_true',
        help='Force rerun even if mosaic results exist'
    )
    parser.add_argument(
        '--output-summary',
        default='./mosaic_batch_summary.csv',
        help='Output summary file (default: ./mosaic_batch_summary.csv)'
    )

    args = parser.parse_args()

    logger.info("=" * 80)
    logger.info("Batch MoSAIC Clustering Analysis")
    logger.info("=" * 80)
    logger.info(f"Input directory: {args.input}")
    logger.info(f"Mode: {args.mode}")
    logger.info(f"Resolution parameter: {args.resolution if args.resolution else 'auto (3rd quartile)'}")
    logger.info(f"Weighted: {args.weighted}")
    logger.info(f"Max workers: {args.max_workers}")
    logger.info(f"Force rerun: {args.force_rerun}")

    # Discover tasks
    tasks = discover_tasks(args.input, force_rerun=args.force_rerun)
    logger.info(f"\nFound {len(tasks)} tasks to process")

    if len(tasks) == 0:
        logger.info("No tasks to process. Exiting.")
        return

    # Prepare task info
    task_infos = []
    for task_dir in tasks:
        task_infos.append({
            'task_name': task_dir.name,
            'input_dir': str(task_dir),
            'output_dir': str(task_dir / 'mosaic'),
            'mode': args.mode,
            'resolution_parameter': args.resolution,
            'weighted': args.weighted
        })

    # Process tasks
    results = []
    start_time = time.time()

    if args.max_workers == 1:
        # Serial processing
        for task_info in task_infos:
            result = process_single_task(task_info)
            results.append(result)
    else:
        # Parallel processing
        with ProcessPoolExecutor(max_workers=args.max_workers) as executor:
            futures = {executor.submit(process_single_task, task_info): task_info
                      for task_info in task_infos}

            for i, future in enumerate(as_completed(futures), 1):
                task_info = futures[future]
                try:
                    result = future.result()
                    results.append(result)
                    logger.info(f"[{i}/{len(tasks)}] {result['task_name']}: {result['status']}")
                except Exception as e:
                    logger.error(f"Task {task_info['task_name']} failed with exception: {e}")
                    results.append({
                        'task_name': task_info['task_name'],
                        'status': 'failed',
                        'error': str(e)
                    })

    elapsed_time = time.time() - start_time

    # Generate summary
    results_df = pd.DataFrame(results)
    results_df.to_csv(args.output_summary, index=False)

    logger.info("\n" + "=" * 80)
    logger.info("BATCH PROCESSING SUMMARY")
    logger.info("=" * 80)
    logger.info(f"Total tasks: {len(tasks)}")
    logger.info(f"Successful: {sum(1 for r in results if r['status'] == 'success')}")
    logger.info(f"Failed: {sum(1 for r in results if r['status'] == 'failed')}")
    logger.info(f"Total time: {elapsed_time:.1f} seconds")
    logger.info(f"Average time per task: {elapsed_time / len(tasks):.1f} seconds")

    # Summary statistics for successful tasks
    successful = results_df[results_df['status'] == 'success']
    if len(successful) > 0:
        logger.info("\nClustering Statistics (successful tasks):")
        logger.info(f"  Mean clusters: {successful['n_clusters'].mean():.1f}")
        logger.info(f"  Mean meaningful clusters (size>1): {successful['n_meaningful_clusters'].mean():.1f}")
        logger.info(f"  Mean silhouette score: {successful['silhouette_score'].mean():.3f}")

    logger.info(f"\nSummary saved to: {args.output_summary}")
    logger.info("=" * 80)


if __name__ == "__main__":
    main()
