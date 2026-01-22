#!/usr/bin/env python3
"""
Batch RMSD Calculation for Missing Tasks (patch directory)

Process the 54 tasks from pbc_1000frames_2step_patch directory.
"""

import sys
import logging
from pathlib import Path
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent))

from scripts.batch_rmsd_hla_alpha import BatchRMSDCalculator

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def main():
    # Read missing tasks list
    missing_csv = Path("output/rmsd_batch_hla_alpha/missing_tasks.csv")

    if not missing_csv.exists():
        logger.error(f"Missing tasks file not found: {missing_csv}")
        return

    missing_df = pd.read_csv(missing_csv)

    # Filter out 8qfy_1 (known problematic task)
    missing_df = missing_df[missing_df['name'] != '8qfy_1']

    logger.info("="*80)
    logger.info("RMSD Analysis for Missing Tasks")
    logger.info("="*80)
    logger.info(f"Total missing tasks: {len(missing_df)}")
    logger.info(f"Excluded: 8qfy_1 (problematic chain IDs)")
    logger.info(f"Tasks to process: {len(missing_df)}")

    # Show task sources
    source_counts = missing_df['source'].value_counts()
    for source, count in source_counts.items():
        logger.info(f"  From {source}: {count} tasks")

    # Prepare task list for BatchRMSDCalculator
    tasks = []
    for idx, row in missing_df.iterrows():
        task_name = row['name']
        task_path = Path(row['path'])

        topology = task_path / "md.tpr"
        trajectory_files = list(task_path.glob("*_processed.xtc"))

        if topology.exists() and trajectory_files:
            tasks.append({
                'name': task_name,
                'topology': str(topology),
                'trajectory': str(trajectory_files[0]),
                'output_dir': f"output/rmsd_batch_hla_alpha/{task_name}"
            })

    logger.info(f"\nValid tasks ready for analysis: {len(tasks)}")

    if len(tasks) == 0:
        logger.error("No valid tasks found!")
        return

    # Initialize calculator
    calculator = BatchRMSDCalculator(
        input_base_dir="input/pbc_1000frames_2step_patch",
        output_base_dir="output/rmsd_batch_hla_alpha",
        max_workers=4
    )

    # Run batch calculation
    logger.info("\n" + "="*80)
    logger.info("Starting batch RMSD calculation...")
    logger.info("="*80 + "\n")

    df_results = calculator.run_batch(tasks=tasks)

    # Merge with existing results
    existing_summary = Path("output/rmsd_batch_hla_alpha/batch_rmsd_summary.csv")

    if existing_summary.exists():
        logger.info("\n" + "="*80)
        logger.info("Merging with existing results...")
        logger.info("="*80)

        df_existing = pd.read_csv(existing_summary)
        df_combined = pd.concat([df_existing, df_results], ignore_index=True)

        # Save combined summary
        df_combined.to_csv(existing_summary, index=False)
        logger.info(f"Combined summary saved: {existing_summary}")
        logger.info(f"  Previous: {len(df_existing)} rows")
        logger.info(f"  New: {len(df_results)} rows")
        logger.info(f"  Total: {len(df_combined)} rows")

        # Update statistics
        df_success = df_combined[df_combined['status'] == 'success']

        stats_file = Path("output/rmsd_batch_hla_alpha/COMPLETE_STATISTICS.txt")
        with open(stats_file, 'w') as f:
            f.write("="*80 + "\n")
            f.write("Complete RMSD Analysis Statistics\n")
            f.write("="*80 + "\n\n")

            f.write(f"Total tasks analyzed: {len(df_combined['task'].unique())}\n")
            f.write(f"Total calculations: {len(df_combined)}\n")
            f.write(f"Successful: {len(df_success)} ({len(df_success)/len(df_combined)*100:.1f}%)\n")
            f.write(f"Failed: {len(df_combined) - len(df_success)}\n\n")

            f.write("="*80 + "\n")
            f.write("Statistics by Component\n")
            f.write("="*80 + "\n\n")

            for calc_type in ['TCR', 'TCR_alpha', 'TCR_beta', 'peptide', 'pHLA']:
                calc_df = df_success[df_success['calc'] == calc_type]
                if len(calc_df) > 0:
                    f.write(f"HLA_alpha -> {calc_type}:\n")
                    f.write(f"  N = {len(calc_df)}\n")
                    f.write(f"  Mean: {calc_df['mean_rmsd_nm'].mean():.4f} ± {calc_df['mean_rmsd_nm'].std():.4f} nm\n")
                    f.write(f"  Median: {calc_df['mean_rmsd_nm'].median():.4f} nm\n")
                    f.write(f"  Range: {calc_df['mean_rmsd_nm'].min():.4f} - {calc_df['mean_rmsd_nm'].max():.4f} nm\n\n")

        logger.info(f"Complete statistics saved: {stats_file}")

    logger.info("\n" + "="*80)
    logger.info("Batch processing complete!")
    logger.info("="*80)


if __name__ == "__main__":
    main()
