#!/usr/bin/env python3
"""
Batch process raw_data_new tasks using 2-step PBC method
Output to pbc_1000frames_2step_patch directory
"""

import sys
import os
import logging
from pathlib import Path
from datetime import datetime
from multiprocessing import Pool
import json

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from immunex.core import PBCProcessor

# Configuration
INPUT_BASE = Path("/home/xumy/work/development/Immunex/input/raw_data_new")
OUTPUT_BASE = Path("/home/xumy/work/development/Immunex/input/pbc_1000frames_2step_patch")
LOGS_DIR = OUTPUT_BASE / "logs"
N_PROCESSES = 4
DT_PS = 100.0

# Setup logging
LOGS_DIR.mkdir(parents=True, exist_ok=True)
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
log_file = LOGS_DIR / f"batch_pbc_2step_{timestamp}.log"

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(log_file),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)


def find_all_md_tasks(base_dir):
    """
    Find all MD tasks with md.xtc and md.tpr files.

    Returns:
        List of tuples: (task_name, trajectory_path, topology_path)
    """
    tasks = []

    # Find all md.xtc files
    for xtc_file in base_dir.rglob("md.xtc"):
        prod_dir = xtc_file.parent
        tpr_file = prod_dir / "md.tpr"

        if not tpr_file.exists():
            logger.warning(f"Missing md.tpr for {xtc_file}")
            continue

        # Extract task name from directory structure
        # Example: .../12-29_run/parallel1/1g6r_run1/prod/md.xtc
        #          -> task_name: 1g6r_run1
        task_dir = prod_dir.parent
        task_name = task_dir.name

        tasks.append({
            'task_name': task_name,
            'trajectory': str(xtc_file),
            'topology': str(tpr_file),
            'source_dir': str(task_dir)
        })

    return tasks


def process_single_task(task_info):
    """
    Process a single MD task using 2-step PBC method.

    Args:
        task_info: Dictionary with task information

    Returns:
        Dictionary with processing results
    """
    task_name = task_info['task_name']
    start_time = datetime.now()

    logger.info(f"Processing {task_name}...")

    result = {
        'task_name': task_name,
        'status': 'failed',
        'error_message': '',
        'processing_time_s': 0.0
    }

    try:
        # Create output directory
        output_dir = OUTPUT_BASE / task_name
        output_dir.mkdir(parents=True, exist_ok=True)

        # Initialize PBC processor
        processor = PBCProcessor(gmx_executable="gmx", keep_temp_files=False)

        # Run 2-step PBC processing
        pbc_result = processor.comprehensive_pbc_process(
            trajectory=task_info['trajectory'],
            topology=task_info['topology'],
            output_dir=str(output_dir),
            method="2step",
            dt=DT_PS,
            auto_dt=True
        )

        # Rename output file to match reference format
        processed_xtc = Path(pbc_result['processed'])
        target_name = output_dir / f"{task_name}_2step_processed.xtc"

        if processed_xtc.exists() and processed_xtc != target_name:
            processed_xtc.rename(target_name)
            logger.info(f"  [{task_name}] Renamed to {target_name.name}")

        # Success
        result['status'] = 'success'
        result['output_dir'] = str(output_dir)
        result['processed_trajectory'] = str(target_name)

        processing_time = (datetime.now() - start_time).total_seconds()
        result['processing_time_s'] = processing_time

        logger.info(f"  [{task_name}] Completed in {processing_time:.1f}s")

    except Exception as e:
        result['error_message'] = str(e)
        logger.error(f"  [{task_name}] Failed: {e}")

    return result


def main():
    """Main processing function"""
    logger.info("=" * 80)
    logger.info("Batch 2-Step PBC Processing")
    logger.info("=" * 80)
    logger.info(f"Input directory: {INPUT_BASE}")
    logger.info(f"Output directory: {OUTPUT_BASE}")
    logger.info(f"Parallel processes: {N_PROCESSES}")
    logger.info(f"Downsampling: dt={DT_PS} ps")
    logger.info("")

    # Discover all tasks
    logger.info("Discovering MD tasks...")
    tasks = find_all_md_tasks(INPUT_BASE)
    logger.info(f"Found {len(tasks)} MD tasks to process")
    logger.info("")

    if len(tasks) == 0:
        logger.error("No tasks found!")
        return

    # Log task summary by source
    from collections import Counter
    sources = Counter([Path(t['source_dir']).parent.parent.name for t in tasks])
    for source, count in sources.items():
        logger.info(f"  {source}: {count} tasks")
    logger.info("")

    # Parallel processing
    start_time = datetime.now()

    with Pool(processes=N_PROCESSES) as pool:
        results = pool.map(process_single_task, tasks)

    total_time = (datetime.now() - start_time).total_seconds()

    # Summarize results
    logger.info("")
    logger.info("=" * 80)
    logger.info("Processing Summary")
    logger.info("=" * 80)

    success_count = sum(1 for r in results if r['status'] == 'success')
    failed_count = len(results) - success_count

    logger.info(f"Total tasks: {len(results)}")
    logger.info(f"Success: {success_count} ({success_count/len(results)*100:.1f}%)")
    logger.info(f"Failed: {failed_count}")
    logger.info(f"Total time: {total_time/60:.1f} minutes")
    logger.info("")

    # Save results to CSV
    csv_file = LOGS_DIR / f"batch_pbc_2step_summary_{timestamp}.csv"

    with open(csv_file, 'w') as f:
        f.write("task_name,status,processing_time_s,output_dir,error_message\n")

        for r in results:
            output_dir = r.get('output_dir', '')
            f.write(f"{r['task_name']},{r['status']},{r['processing_time_s']:.2f},"
                   f"{output_dir},\"{r['error_message']}\"\n")

    logger.info(f"Results saved to: {csv_file}")
    logger.info("")

    # Failed tasks list
    if failed_count > 0:
        logger.info("Failed tasks:")
        for r in results:
            if r['status'] != 'success':
                logger.info(f"  - {r['task_name']}: {r['error_message']}")
        logger.info("")

    # Save detailed results to JSON
    json_file = LOGS_DIR / f"batch_pbc_2step_results_{timestamp}.json"
    with open(json_file, 'w') as f:
        json.dump(results, f, indent=2)

    logger.info(f"Detailed results saved to: {json_file}")
    logger.info("")

    logger.info("=" * 80)
    logger.info("Processing completed!")
    logger.info("=" * 80)


if __name__ == "__main__":
    main()
