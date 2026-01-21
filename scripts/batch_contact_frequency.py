#!/usr/bin/env python3
"""
Batch Contact Frequency Calculation

Recalculate contact frequency using correct definition:
  frequency = (frames with contact) / (total frames)

A frame has contact if ANY atom pair is within cutoff distance.
"""

import sys
import logging
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent))

from aftermd.analysis.trajectory.contact_number import ContactNumberCalculator

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def process_single_task(task_info):
    """
    Process single task to calculate contact frequency.

    Args:
        task_info: Dict with keys:
            - name: task name
            - topology: path to topology file
            - trajectory: path to trajectory file
            - output_dir: output directory
            - selections: dict of selection pairs to analyze
            - cutoff: distance cutoff (default 4.5 A)

    Returns:
        Dict with task results
    """
    task_name = task_info['name']
    topology = task_info['topology']
    trajectory = task_info['trajectory']
    output_dir = Path(task_info['output_dir'])
    selections = task_info.get('selections', {})
    cutoff = task_info.get('cutoff', 4.5)

    logger.info(f"="*60)
    logger.info(f"Processing task: {task_name}")
    logger.info(f"="*60)

    try:
        # Initialize calculator
        calc = ContactNumberCalculator(topology, trajectory)

        results = {
            'task': task_name,
            'status': 'success',
            'contacts': {}
        }

        # Calculate frequency for each selection pair
        for contact_name, (sel1, sel2) in selections.items():
            logger.info(f"\n--- {contact_name} ---")

            output_file = output_dir / f"{task_name}_{contact_name}_freq.csv"

            try:
                times, contact_status, frequency = calc.calculate_contact_frequency(
                    selection1=sel1,
                    selection2=sel2,
                    cutoff=cutoff,
                    heavy_atoms_only=True,
                    output_file=str(output_file),
                    time_unit='ps'
                )

                results['contacts'][contact_name] = {
                    'frequency': frequency,
                    'frames_with_contact': int(contact_status.sum()),
                    'total_frames': len(contact_status),
                    'output_file': str(output_file)
                }

                logger.info(f"  -> Contact frequency: {frequency:.4f} ({frequency*100:.2f}%)")

            except Exception as e:
                logger.error(f"  -> Failed: {e}")
                results['contacts'][contact_name] = {
                    'error': str(e)
                }

        return results

    except Exception as e:
        logger.error(f"Task {task_name} failed: {e}")
        return {
            'task': task_name,
            'status': 'failed',
            'error': str(e)
        }


def main():
    """Main batch processing function."""

    # Define standard pHLA-TCR contact pairs
    # Chain mapping:
    #   A: HLA-alpha
    #   B: HLA-beta
    #   C: peptide
    #   D: TCR-alpha
    #   E: TCR-beta
    standard_selections = {
        'TCR_peptide': (
            'chainID D E',  # TCR (alpha + beta)
            'chainID C'     # peptide
        ),
        'TCR_pHLA': (
            'chainID D E',    # TCR
            'chainID A B C'   # pHLA (HLA-alpha + HLA-beta + peptide)
        ),
        'CDR3_peptide': (
            '(chainID D and resid 103-115) or (chainID E and resid 103-115)',  # CDR3 regions
            'chainID C'
        ),
        'TCRalpha_peptide': (
            'chainID D',  # TCR-alpha only
            'chainID C'   # peptide
        ),
        'TCRbeta_peptide': (
            'chainID E',  # TCR-beta only
            'chainID C'   # peptide
        )
    }

    # Discover all processed trajectory tasks
    input_dirs = [
        Path("input/pbc_1000frames_2step"),
        Path("input/pbc_1000frames_2step_patch")
    ]

    tasks = []
    for input_dir in input_dirs:
        if not input_dir.exists():
            logger.warning(f"Directory not found: {input_dir}")
            continue

        for task_dir in sorted(input_dir.iterdir()):
            if not task_dir.is_dir():
                continue

            task_name = task_dir.name
            topology = task_dir / "md.tpr"
            trajectory_files = list(task_dir.glob("*_processed.xtc"))

            if not topology.exists():
                logger.warning(f"Missing topology: {task_name}")
                continue

            if not trajectory_files:
                logger.warning(f"Missing processed trajectory: {task_name}")
                continue

            trajectory = trajectory_files[0]

            tasks.append({
                'name': task_name,
                'topology': str(topology),
                'trajectory': str(trajectory),
                'output_dir': 'output/contact_frequency_batch',
                'selections': standard_selections,
                'cutoff': 4.5
            })

    logger.info("="*80)
    logger.info("Batch Contact Frequency Calculation")
    logger.info("="*80)
    logger.info(f"Total tasks: {len(tasks)}")
    logger.info(f"Contact pairs per task: {len(standard_selections)}")
    logger.info(f"Distance cutoff: 4.5 A")
    logger.info(f"Output directory: output/contact_frequency_batch")
    logger.info("="*80)

    if len(tasks) == 0:
        logger.error("No valid tasks found!")
        return

    # Create output directory
    output_base = Path("output/contact_frequency_batch")
    output_base.mkdir(parents=True, exist_ok=True)

    # Process tasks in parallel
    max_workers = 4
    logger.info(f"\nStarting parallel processing with {max_workers} workers...")

    all_results = []

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(process_single_task, task): task for task in tasks}

        for future in as_completed(futures):
            task = futures[future]
            try:
                result = future.result()
                all_results.append(result)

                status = result.get('status', 'unknown')
                logger.info(f"\nCompleted: {result['task']} - {status}")

            except Exception as e:
                logger.error(f"\nTask {task['name']} raised exception: {e}")
                all_results.append({
                    'task': task['name'],
                    'status': 'exception',
                    'error': str(e)
                })

    # Generate summary report
    logger.info("\n" + "="*80)
    logger.info("Generating summary report...")
    logger.info("="*80)

    summary_rows = []
    for result in all_results:
        task_name = result['task']
        status = result['status']

        if status == 'success':
            for contact_name, contact_data in result.get('contacts', {}).items():
                if 'frequency' in contact_data:
                    summary_rows.append({
                        'task': task_name,
                        'contact_pair': contact_name,
                        'frequency': contact_data['frequency'],
                        'frames_with_contact': contact_data['frames_with_contact'],
                        'total_frames': contact_data['total_frames'],
                        'status': 'success'
                    })
                else:
                    summary_rows.append({
                        'task': task_name,
                        'contact_pair': contact_name,
                        'frequency': None,
                        'frames_with_contact': None,
                        'total_frames': None,
                        'status': 'failed',
                        'error': contact_data.get('error', 'Unknown error')
                    })
        else:
            summary_rows.append({
                'task': task_name,
                'contact_pair': 'all',
                'frequency': None,
                'frames_with_contact': None,
                'total_frames': None,
                'status': status,
                'error': result.get('error', 'Unknown error')
            })

    # Save summary CSV
    summary_df = pd.DataFrame(summary_rows)
    summary_file = output_base / "batch_contact_frequency_summary.csv"
    summary_df.to_csv(summary_file, index=False)
    logger.info(f"\nSummary saved: {summary_file}")

    # Statistics
    success_count = len([r for r in all_results if r['status'] == 'success'])
    failed_count = len(all_results) - success_count

    logger.info("\n" + "="*80)
    logger.info("Batch Processing Complete")
    logger.info("="*80)
    logger.info(f"Total tasks: {len(all_results)}")
    logger.info(f"Successful: {success_count}")
    logger.info(f"Failed: {failed_count}")

    # Contact frequency statistics
    if success_count > 0:
        logger.info("\nContact Frequency Statistics:")
        for contact_name in standard_selections.keys():
            contact_freqs = [
                r['contacts'][contact_name]['frequency']
                for r in all_results
                if r['status'] == 'success' and contact_name in r.get('contacts', {})
                and 'frequency' in r['contacts'][contact_name]
            ]

            if contact_freqs:
                import numpy as np
                logger.info(f"\n  {contact_name}:")
                logger.info(f"    N = {len(contact_freqs)}")
                logger.info(f"    Mean: {np.mean(contact_freqs):.4f} ({np.mean(contact_freqs)*100:.2f}%)")
                logger.info(f"    Std: {np.std(contact_freqs):.4f}")
                logger.info(f"    Range: {np.min(contact_freqs):.4f} - {np.max(contact_freqs):.4f}")

    logger.info("\n" + "="*80)


if __name__ == "__main__":
    main()
