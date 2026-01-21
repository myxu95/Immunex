#!/usr/bin/env python3
"""
Batch Analysis Script for pHLA-TCR Complex Trajectories

Process multiple MD trajectories in parallel using the pHLATCRAnalyzer.

Author: AfterMD Development Team
Date: 2025-12-07
"""

import pandas as pd
import numpy as np
from pathlib import Path
from multiprocessing import Pool, cpu_count
import logging
import sys
import json
from datetime import datetime
from typing import List, Dict, Optional
import argparse

from aftermd.analysis.trajectory import pHLATCRAnalyzer

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class BatchpHLAAnalysis:
    """Batch processing manager for pHLA-TCR trajectory analysis."""

    def __init__(self,
                 processed_dir: str,
                 output_dir: str,
                 task_list_file: Optional[str] = None,
                 quality_filter: Optional[str] = None,
                 max_workers: int = 8):
        """
        Initialize batch analysis.

        Args:
            processed_dir: Directory containing processed trajectories
            output_dir: Output directory for analysis results
            task_list_file: Optional file with task names (one per line)
            quality_filter: Filter by quality ('excellent', 'good', 'acceptable')
            max_workers: Number of parallel workers
        """
        self.processed_dir = Path(processed_dir)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        self.max_workers = max_workers
        self.quality_filter = quality_filter

        # Load task list
        self.tasks = self._load_tasks(task_list_file)

        # Load quality information if filtering
        if quality_filter:
            self._apply_quality_filter()

        logger.info(f"Initialized batch analysis with {len(self.tasks)} tasks")
        logger.info(f"Max workers: {self.max_workers}")

    def _load_tasks(self, task_list_file: Optional[str]) -> List[Dict]:
        """Load list of tasks to process."""
        tasks = []

        if task_list_file:
            # Load from file
            with open(task_list_file, 'r') as f:
                task_names = [line.strip() for line in f if line.strip()]

            for task_name in task_names:
                task_info = self._find_task_files(task_name)
                if task_info:
                    tasks.append(task_info)
        else:
            # Auto-discover from processed_dir
            tasks = self._discover_tasks()

        return tasks

    def _discover_tasks(self) -> List[Dict]:
        """Auto-discover tasks from processed directory."""
        logger.info(f"Discovering tasks in {self.processed_dir}...")

        tasks = []

        # Iterate through PDB ID directories
        for pdb_dir in sorted(self.processed_dir.iterdir()):
            if not pdb_dir.is_dir():
                continue

            # Iterate through variant directories
            for variant_dir in sorted(pdb_dir.iterdir()):
                if not variant_dir.is_dir():
                    continue

                task_name = f"{pdb_dir.name}_{variant_dir.name}"

                # Check for required files
                xtc_file = variant_dir / f"{pdb_dir.name}_processed.xtc"
                pdb_file = variant_dir / f"{pdb_dir.name}_standardized.pdb"

                if xtc_file.exists() and pdb_file.exists():
                    tasks.append({
                        'task_name': task_name,
                        'pdb_id': pdb_dir.name,
                        'variant': variant_dir.name,
                        'trajectory': str(xtc_file),
                        'topology': str(pdb_file),
                        'task_dir': str(variant_dir)
                    })

        logger.info(f"Discovered {len(tasks)} tasks")
        return tasks

    def _find_task_files(self, task_name: str) -> Optional[Dict]:
        """Find trajectory and topology files for a specific task."""
        # Parse task name (format: pdbid_variant)
        parts = task_name.rsplit('_', 1)
        if len(parts) != 2:
            logger.warning(f"Invalid task name format: {task_name}")
            return None

        pdb_id, variant = parts

        # Construct paths
        task_dir = self.processed_dir / pdb_id / variant
        xtc_file = task_dir / f"{pdb_id}_processed.xtc"
        pdb_file = task_dir / f"{pdb_id}_standardized.pdb"

        if not xtc_file.exists():
            logger.warning(f"Trajectory not found: {xtc_file}")
            return None

        if not pdb_file.exists():
            logger.warning(f"Topology not found: {pdb_file}")
            return None

        return {
            'task_name': task_name,
            'pdb_id': pdb_id,
            'variant': variant,
            'trajectory': str(xtc_file),
            'topology': str(pdb_file),
            'task_dir': str(task_dir)
        }

    def _apply_quality_filter(self):
        """Filter tasks based on quality criteria."""
        summary_file = self.processed_dir.parent / 'logs' / 'pbc_summary_final.csv'

        if not summary_file.exists():
            logger.warning(f"Quality summary not found: {summary_file}")
            return

        logger.info(f"Applying quality filter: {self.quality_filter}")

        # Load quality data
        df_quality = pd.read_csv(summary_file)

        # Define quality thresholds
        quality_map = {
            'excellent': ['excellent'],
            'good': ['excellent', 'good'],
            'acceptable': ['excellent', 'good', 'moderate']
        }

        allowed_qualities = quality_map.get(self.quality_filter, [])

        # Filter tasks
        filtered_tasks = []
        for task in self.tasks:
            task_row = df_quality[df_quality['task_name'] == task['task_name']]
            if not task_row.empty:
                quality = task_row.iloc[0]['quality']
                if quality in allowed_qualities:
                    task['quality'] = quality
                    filtered_tasks.append(task)

        logger.info(f"Tasks after quality filter: {len(filtered_tasks)} / {len(self.tasks)}")
        self.tasks = filtered_tasks

    def process_single_task(self, task_info: Dict) -> Dict:
        """
        Process a single task.

        Args:
            task_info: Task information dictionary

        Returns:
            Result dictionary
        """
        task_name = task_info['task_name']
        logger.info(f"Processing {task_name}...")

        try:
            # Create analyzer
            analyzer = pHLATCRAnalyzer(
                task_name=task_name,
                trajectory=task_info['trajectory'],
                topology=task_info['topology'],
                output_dir=str(self.output_dir / task_name)
            )

            # Run analysis
            results = analyzer.run_full_analysis()

            return {
                'task_name': task_name,
                'status': 'success',
                'results': results
            }

        except Exception as e:
            logger.error(f"Failed to process {task_name}: {e}")
            return {
                'task_name': task_name,
                'status': 'failed',
                'error': str(e)
            }

    def run_batch_analysis(self) -> pd.DataFrame:
        """
        Run analysis on all tasks in parallel.

        Returns:
            DataFrame with summary results
        """
        logger.info(f"Starting batch analysis of {len(self.tasks)} tasks...")
        logger.info(f"Using {self.max_workers} parallel workers")

        start_time = datetime.now()

        # Process in parallel
        with Pool(processes=self.max_workers) as pool:
            results = pool.map(self.process_single_task, self.tasks)

        elapsed = (datetime.now() - start_time).total_seconds()

        # Compile results
        summary_data = []

        for result in results:
            task_name = result['task_name']
            status = result['status']

            if status == 'success':
                res = result['results']

                # Extract key metrics
                row = {
                    'task_name': task_name,
                    'status': status,
                    'global_rmsd_mean': res['global_rmsd']['mean'],
                    'global_rmsd_std': res['global_rmsd']['std'],
                    'peptide_rmsd_mean': res['chain_rmsd']['Peptide']['mean'],
                    'hla_alpha_rmsd_mean': res['chain_rmsd']['HLA_alpha']['mean'],
                    'hla_beta_rmsd_mean': res['chain_rmsd']['HLA_beta']['mean'],
                    'tcr_alpha_rmsd_mean': res['chain_rmsd']['TCR_alpha']['mean'],
                    'tcr_beta_rmsd_mean': res['chain_rmsd']['TCR_beta']['mean'],
                    'peptide_hla_dist_mean': res['interface_distances']['Peptide_HLA_COM']['mean'],
                    'tcr_peptide_dist_mean': res['interface_distances']['TCR_Peptide_min']['mean'],
                    'peptide_hla_contacts_mean': res['contacts']['Peptide_HLA']['mean_contacts'],
                    'tcr_peptide_contacts_mean': res['contacts']['TCR_Peptide']['mean_contacts'],
                }
            else:
                row = {
                    'task_name': task_name,
                    'status': status,
                    'error': result.get('error', 'Unknown error')
                }

            summary_data.append(row)

        # Create summary DataFrame
        df_summary = pd.DataFrame(summary_data)

        # Save summary
        summary_file = self.output_dir / f'batch_analysis_summary_{datetime.now().strftime("%Y%m%d_%H%M%S")}.csv'
        df_summary.to_csv(summary_file, index=False)

        # Print statistics
        n_success = sum(1 for r in results if r['status'] == 'success')
        n_failed = len(results) - n_success

        logger.info("\n" + "="*60)
        logger.info("BATCH ANALYSIS COMPLETE")
        logger.info("="*60)
        logger.info(f"Total tasks: {len(results)}")
        logger.info(f"Successful: {n_success}")
        logger.info(f"Failed: {n_failed}")
        logger.info(f"Total time: {elapsed:.2f} seconds")
        logger.info(f"Average time per task: {elapsed/len(results):.2f} seconds")
        logger.info(f"Summary saved to: {summary_file}")
        logger.info("="*60)

        return df_summary


def main():
    """Command-line interface."""
    parser = argparse.ArgumentParser(
        description='Batch analysis of pHLA-TCR MD trajectories',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process all tasks in directory
  python batch_phla_analysis.py /path/to/processed_once/output ./results

  # Process specific tasks from file
  python batch_phla_analysis.py /path/to/processed_once/output ./results --task-list tasks.txt

  # Process only excellent quality tasks
  python batch_phla_analysis.py /path/to/processed_once/output ./results --quality excellent

  # Use 16 parallel workers
  python batch_phla_analysis.py /path/to/processed_once/output ./results --workers 16
        """
    )

    parser.add_argument('processed_dir', help='Directory with processed trajectories')
    parser.add_argument('output_dir', help='Output directory for analysis results')
    parser.add_argument('--task-list', help='File with task names (one per line)')
    parser.add_argument('--quality', choices=['excellent', 'good', 'acceptable'],
                       help='Filter by quality level')
    parser.add_argument('--workers', type=int, default=8,
                       help='Number of parallel workers (default: 8)')

    args = parser.parse_args()

    # Create batch processor
    batch = BatchpHLAAnalysis(
        processed_dir=args.processed_dir,
        output_dir=args.output_dir,
        task_list_file=args.task_list,
        quality_filter=args.quality,
        max_workers=args.workers
    )

    # Run analysis
    df_summary = batch.run_batch_analysis()

    print("\n=== Summary Statistics ===")
    if 'global_rmsd_mean' in df_summary.columns:
        print(f"Mean global RMSD: {df_summary['global_rmsd_mean'].mean():.2f} ± {df_summary['global_rmsd_mean'].std():.2f} Å")
        print(f"Mean peptide RMSD: {df_summary['peptide_rmsd_mean'].mean():.2f} ± {df_summary['peptide_rmsd_mean'].std():.2f} Å")
        print(f"Mean peptide-HLA distance: {df_summary['peptide_hla_dist_mean'].mean():.2f} ± {df_summary['peptide_hla_dist_mean'].std():.2f} Å")


if __name__ == "__main__":
    main()
