#!/usr/bin/env python3
"""
Batch RMSD Calculation Script - Align to HLA_alpha

Process all MD trajectories with 5 RMSD calculations:
1. Align HLA_alpha -> Calculate TCR
2. Align HLA_alpha -> Calculate TCR_alpha
3. Align HLA_alpha -> Calculate TCR_beta
4. Align HLA_alpha -> Calculate peptide
5. Align HLA_alpha -> Calculate pHLA

Author: AfterMD Development Team
"""

import sys
import logging
from pathlib import Path
from typing import List, Dict, Optional
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed
import traceback

sys.path.insert(0, str(Path(__file__).parent.parent))

from aftermd.analysis.trajectory import RMSDInterface

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class BatchRMSDCalculator:
    """Batch RMSD calculator for multiple MD tasks."""

    CALCULATIONS = [
        {'align': 'HLA_alpha', 'calc': 'TCR', 'description': 'TCR overall motion relative to HLA-alpha'},
        {'align': 'HLA_alpha', 'calc': 'TCR_alpha', 'description': 'TCR-alpha chain motion relative to HLA-alpha'},
        {'align': 'HLA_alpha', 'calc': 'TCR_beta', 'description': 'TCR-beta chain motion relative to HLA-alpha'},
        {'align': 'HLA_alpha', 'calc': 'peptide', 'description': 'Peptide motion relative to HLA-alpha'},
        {'align': 'HLA_alpha', 'calc': 'pHLA', 'description': 'pMHC complex motion relative to HLA-alpha'},
    ]

    def __init__(self,
                 input_base_dir: str,
                 output_base_dir: str,
                 max_workers: int = 4):
        """
        Initialize batch RMSD calculator.

        Args:
            input_base_dir: Base directory containing all MD tasks
            output_base_dir: Base directory for output results
            max_workers: Number of parallel workers
        """
        self.input_base_dir = Path(input_base_dir)
        self.output_base_dir = Path(output_base_dir)
        self.max_workers = max_workers

        self.output_base_dir.mkdir(parents=True, exist_ok=True)

    def discover_tasks(self) -> List[Dict]:
        """
        Discover all valid MD tasks.

        Returns:
            List of task dictionaries with topology and trajectory paths
        """
        tasks = []

        for task_dir in sorted(self.input_base_dir.iterdir()):
            if not task_dir.is_dir():
                continue

            task_name = task_dir.name
            topology = task_dir / "md.tpr"

            # Find processed trajectory
            trajectory_files = list(task_dir.glob("*_processed.xtc"))

            if not topology.exists():
                logger.warning(f"[{task_name}] Missing md.tpr, skipping")
                continue

            if not trajectory_files:
                logger.warning(f"[{task_name}] No processed trajectory found, skipping")
                continue

            trajectory = trajectory_files[0]

            tasks.append({
                'name': task_name,
                'topology': str(topology),
                'trajectory': str(trajectory),
                'output_dir': str(self.output_base_dir / task_name)
            })

        logger.info(f"Discovered {len(tasks)} valid MD tasks")
        return tasks

    def process_single_task(self, task: Dict) -> Dict:
        """
        Process a single MD task with 5 RMSD calculations.

        Args:
            task: Task dictionary

        Returns:
            Results dictionary
        """
        task_name = task['name']

        try:
            logger.info(f"[{task_name}] Starting RMSD calculations...")

            # Initialize RMSD Interface
            rmsd = RMSDInterface(
                topology=task['topology'],
                trajectory=task['trajectory'],
                output_dir=task['output_dir'],
                auto_standardize=True
            )

            results = []

            for i, calc_config in enumerate(self.CALCULATIONS, 1):
                try:
                    logger.info(f"[{task_name}] Calculation {i}/5: {calc_config['align']} -> {calc_config['calc']}")

                    result = rmsd.calculate(
                        align=calc_config['align'],
                        calc=calc_config['calc']
                    )

                    results.append({
                        'task': task_name,
                        'align': calc_config['align'],
                        'calc': calc_config['calc'],
                        'description': calc_config['description'],
                        'mean_rmsd_nm': result['mean'],
                        'std_rmsd_nm': result['std'],
                        'min_rmsd_nm': result['min'],
                        'max_rmsd_nm': result['max'],
                        'n_frames': result['n_frames'],
                        'output_file': result['output_file'],
                        'status': 'success'
                    })

                except Exception as e:
                    logger.error(f"[{task_name}] Calculation {i} failed: {e}")
                    results.append({
                        'task': task_name,
                        'align': calc_config['align'],
                        'calc': calc_config['calc'],
                        'description': calc_config['description'],
                        'status': 'failed',
                        'error': str(e)
                    })

            logger.info(f"[{task_name}] Completed {len([r for r in results if r['status'] == 'success'])}/5 calculations")

            return {
                'task': task_name,
                'status': 'success',
                'results': results
            }

        except Exception as e:
            logger.error(f"[{task_name}] Task failed: {e}")
            logger.error(traceback.format_exc())
            return {
                'task': task_name,
                'status': 'failed',
                'error': str(e)
            }

    def run_batch(self, tasks: Optional[List[Dict]] = None) -> pd.DataFrame:
        """
        Run batch RMSD calculations.

        Args:
            tasks: Optional list of tasks (if None, auto-discover)

        Returns:
            DataFrame with all results
        """
        if tasks is None:
            tasks = self.discover_tasks()

        logger.info("="*80)
        logger.info(f"Starting batch RMSD calculation for {len(tasks)} tasks")
        logger.info(f"Using {self.max_workers} parallel workers")
        logger.info("="*80)

        all_results = []
        completed = 0
        failed = 0

        with ProcessPoolExecutor(max_workers=self.max_workers) as executor:
            future_to_task = {
                executor.submit(self.process_single_task, task): task
                for task in tasks
            }

            for future in as_completed(future_to_task):
                task = future_to_task[future]
                try:
                    task_result = future.result()

                    if task_result['status'] == 'success':
                        all_results.extend(task_result['results'])
                        completed += 1
                    else:
                        failed += 1

                except Exception as e:
                    logger.error(f"Task {task['name']} execution error: {e}")
                    failed += 1

        logger.info("="*80)
        logger.info(f"Batch processing complete!")
        logger.info(f"  Completed: {completed}/{len(tasks)} tasks")
        logger.info(f"  Failed: {failed}/{len(tasks)} tasks")
        logger.info("="*80)

        # Convert to DataFrame
        df = pd.DataFrame(all_results)

        # Save summary
        summary_file = self.output_base_dir / "batch_rmsd_summary.csv"
        df.to_csv(summary_file, index=False)
        logger.info(f"Summary saved to: {summary_file}")

        # Generate statistics report
        self._generate_statistics_report(df)

        return df

    def _generate_statistics_report(self, df: pd.DataFrame):
        """Generate statistics report."""

        report_file = self.output_base_dir / "batch_rmsd_statistics.txt"

        with open(report_file, 'w') as f:
            f.write("="*80 + "\n")
            f.write("Batch RMSD Calculation - Statistics Report\n")
            f.write("="*80 + "\n\n")

            # Overall statistics
            total_calculations = len(df)
            successful = len(df[df['status'] == 'success'])
            failed = len(df[df['status'] == 'failed'])

            f.write(f"Total calculations: {total_calculations}\n")
            f.write(f"Successful: {successful} ({successful/total_calculations*100:.1f}%)\n")
            f.write(f"Failed: {failed} ({failed/total_calculations*100:.1f}%)\n\n")

            # Statistics by calculation type
            if successful > 0:
                success_df = df[df['status'] == 'success']

                f.write("-"*80 + "\n")
                f.write("Statistics by Calculation Type\n")
                f.write("-"*80 + "\n\n")

                for calc_type in ['TCR', 'TCR_alpha', 'TCR_beta', 'peptide', 'pHLA']:
                    calc_df = success_df[success_df['calc'] == calc_type]

                    if len(calc_df) > 0:
                        f.write(f"Calculation: HLA_alpha -> {calc_type}\n")
                        f.write(f"  N = {len(calc_df)} tasks\n")
                        f.write(f"  Mean RMSD: {calc_df['mean_rmsd_nm'].mean():.4f} ± {calc_df['mean_rmsd_nm'].std():.4f} nm\n")
                        f.write(f"  Median RMSD: {calc_df['mean_rmsd_nm'].median():.4f} nm\n")
                        f.write(f"  Range: {calc_df['mean_rmsd_nm'].min():.4f} - {calc_df['mean_rmsd_nm'].max():.4f} nm\n\n")

        logger.info(f"Statistics report saved to: {report_file}")


def main():
    """Main function."""

    input_dir = "/home/xumy/work/development/AfterMD/input/pbc_1000frames_2step"
    output_dir = "/home/xumy/work/development/AfterMD/output/rmsd_batch_hla_alpha"

    calculator = BatchRMSDCalculator(
        input_base_dir=input_dir,
        output_base_dir=output_dir,
        max_workers=4
    )

    df_results = calculator.run_batch()

    logger.info(f"\nBatch processing complete!")
    logger.info(f"Results saved to: {output_dir}")


if __name__ == "__main__":
    main()
