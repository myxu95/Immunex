#!/usr/bin/env python3
"""
Generic batch worker for Immunex SLURM processing.

This script reads task configurations from JSON files and executes
batch processing jobs. Supports different task types (PBC, RMSD, etc.).
"""

import sys
import json
import argparse
import logging
from pathlib import Path
from typing import Dict, Any, List
from datetime import datetime

# Add Immunex to path
IMMUNEX_PATH = "/public/home/xmy/tools/Immunex"
if Path(IMMUNEX_PATH).exists():
    sys.path.insert(0, IMMUNEX_PATH)

from immunex.core.pbc_processor import PBCProcessor
from immunex.analysis.trajectory.rmsd import RMSDCalculator


class BatchWorker:
    """Generic worker for batch task processing."""

    def __init__(self, config_file: str, log_file: str = None):
        """
        Initialize batch worker.

        Args:
            config_file: Path to JSON configuration file
            log_file: Optional log file path
        """
        self.config_file = config_file
        self.config = self._load_config()
        self._setup_logging(log_file)

    def _load_config(self) -> Dict[str, Any]:
        """Load configuration from JSON file."""
        config_path = Path(self.config_file)
        if not config_path.exists():
            raise FileNotFoundError(f"Config file not found: {self.config_file}")

        with open(config_path) as f:
            config = json.load(f)

        return config

    def _setup_logging(self, log_file: str = None):
        """Setup logging configuration."""
        log_format = '%(asctime)s - %(levelname)s - %(message)s'

        if log_file:
            logging.basicConfig(
                level=logging.INFO,
                format=log_format,
                handlers=[
                    logging.FileHandler(log_file),
                    logging.StreamHandler(sys.stdout)
                ]
            )
        else:
            logging.basicConfig(
                level=logging.INFO,
                format=log_format
            )

        self.logger = logging.getLogger(__name__)

    def process_pbc_task(self, task_name: str, task_data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Process a single PBC correction task.

        Args:
            task_name: Task identifier
            task_data: Task configuration including trajectory, topology, output_dir, etc.

        Returns:
            Result dictionary
        """
        try:
            trajectory = task_data['trajectory']
            topology = task_data['topology']
            output_dir = task_data['output_dir']
            dt = task_data.get('dt')

            self.logger.info(f"Processing PBC task: {task_name}")
            self.logger.info(f"  Trajectory: {trajectory}")
            self.logger.info(f"  Topology: {topology}")
            self.logger.info(f"  Output: {output_dir}")

            # Create output directory
            Path(output_dir).mkdir(parents=True, exist_ok=True)

            # Initialize PBC processor
            pbc_processor = PBCProcessor(gmx_executable='gmx')

            # Process task
            result = pbc_processor.comprehensive_pbc_process(
                trajectory=trajectory,
                topology=topology,
                output_dir=output_dir,
                dt=dt
            )

            self.logger.info(f"Task {task_name} completed successfully")

            return {
                'task_name': task_name,
                'status': 'success',
                'result': result
            }

        except Exception as e:
            self.logger.error(f"Task {task_name} failed: {e}")
            return {
                'task_name': task_name,
                'status': 'failed',
                'error': str(e)
            }

    def process_rmsd_task(self, task_name: str, task_data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Process a single RMSD calculation task.

        Args:
            task_name: Task identifier
            task_data: Task configuration including trajectory, topology, output_file, etc.

        Returns:
            Result dictionary
        """
        try:
            trajectory = task_data['trajectory']
            topology = task_data['topology']
            output_file = task_data['output_file']
            rmsd_type = task_data.get('rmsd_type', 'calpha')
            fit_group = task_data.get('fit_group')
            calc_group = task_data.get('calc_group')

            self.logger.info(f"Processing RMSD task: {task_name}")
            self.logger.info(f"  Trajectory: {trajectory}")
            self.logger.info(f"  Topology: {topology}")
            self.logger.info(f"  Output: {output_file}")
            self.logger.info(f"  RMSD type: {rmsd_type}")

            # Create output directory
            Path(output_file).parent.mkdir(parents=True, exist_ok=True)

            # Initialize RMSD calculator
            rmsd_calc = RMSDCalculator(topology, trajectory)

            # Calculate RMSD
            if fit_group and calc_group:
                result = rmsd_calc.calculate_gromacs_custom_groups(
                    reference_group=fit_group,
                    analysis_group=calc_group,
                    output_file=output_file
                )
            else:
                result = rmsd_calc.calculate_gromacs(
                    rmsd_type=rmsd_type,
                    output_file=output_file
                )

            self.logger.info(f"Task {task_name} completed successfully")

            return {
                'task_name': task_name,
                'status': 'success',
                'output_file': result
            }

        except Exception as e:
            self.logger.error(f"Task {task_name} failed: {e}")
            return {
                'task_name': task_name,
                'status': 'failed',
                'error': str(e)
            }

    def run(self) -> Dict[str, Any]:
        """
        Run batch processing based on configuration.

        Returns:
            Summary results dictionary
        """
        task_type = self.config.get('task_type', 'pbc')
        tasks = self.config.get('tasks', {})

        self.logger.info("=" * 60)
        self.logger.info(f"Batch Worker Started")
        self.logger.info(f"Task type: {task_type}")
        self.logger.info(f"Total tasks: {len(tasks)}")
        self.logger.info(f"Config file: {self.config_file}")
        self.logger.info("=" * 60)

        results = []

        # Process each task
        for task_name, task_data in tasks.items():
            if task_type == 'pbc':
                result = self.process_pbc_task(task_name, task_data)
            elif task_type == 'rmsd':
                result = self.process_rmsd_task(task_name, task_data)
            else:
                self.logger.error(f"Unknown task type: {task_type}")
                result = {
                    'task_name': task_name,
                    'status': 'failed',
                    'error': f'Unknown task type: {task_type}'
                }

            results.append(result)

        # Summary
        successful = len([r for r in results if r['status'] == 'success'])
        failed = len([r for r in results if r['status'] == 'failed'])

        self.logger.info("=" * 60)
        self.logger.info(f"Batch Processing Summary")
        self.logger.info(f"Total tasks: {len(results)}")
        self.logger.info(f"Successful: {successful}")
        self.logger.info(f"Failed: {failed}")
        self.logger.info("=" * 60)

        if failed > 0:
            self.logger.warning("Failed tasks:")
            for result in results:
                if result['status'] == 'failed':
                    self.logger.warning(f"  - {result['task_name']}: {result.get('error', 'Unknown error')}")

        return {
            'total': len(results),
            'successful': successful,
            'failed': failed,
            'results': results
        }


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Immunex batch worker for SLURM processing",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process PBC batch
  python batch_worker.py --config batch_001_config.json

  # Process RMSD batch with custom log file
  python batch_worker.py --config rmsd_batch_001.json --log-file batch.log

Configuration file format (JSON):
  {
    "task_type": "pbc",
    "tasks": {
      "task1": {
        "trajectory": "/path/to/md.xtc",
        "topology": "/path/to/md.tpr",
        "output_dir": "/path/to/output",
        "dt": null
      }
    }
  }
        """
    )

    parser.add_argument(
        "--config",
        required=True,
        help="Path to JSON configuration file"
    )

    parser.add_argument(
        "--log-file",
        help="Optional log file path (default: stdout only)"
    )

    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Enable verbose logging"
    )

    args = parser.parse_args()

    try:
        # Initialize worker
        worker = BatchWorker(
            config_file=args.config,
            log_file=args.log_file
        )

        # Run batch processing
        summary = worker.run()

        # Exit with appropriate code
        if summary['failed'] > 0:
            sys.exit(1)
        else:
            sys.exit(0)

    except Exception as e:
        logging.error(f"Batch worker failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
