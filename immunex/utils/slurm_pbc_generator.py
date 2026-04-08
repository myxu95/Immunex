#!/usr/bin/env python3
"""
SLURM batch processor for PBC corrections.

This module generates SLURM scripts for high-throughput PBC processing
of MD trajectories.
"""

from pathlib import Path
from typing import Dict, Any, Optional, Tuple
import logging
from .slurm_batch_base import SlurmBatchProcessor

logger = logging.getLogger(__name__)


class PBCBatchProcessor(SlurmBatchProcessor):
    """SLURM batch processor for PBC corrections."""

    def __init__(self,
                 template_file: Optional[str] = None,
                 default_slurm_params: Optional[Dict[str, Any]] = None):
        """
        Initialize PBC batch processor.

        Args:
            template_file: Path to custom SLURM template file
            default_slurm_params: Default SLURM parameters
        """
        super().__init__(template_file, default_slurm_params)

    def discover_tasks(self, input_dir: str, **kwargs) -> Dict[str, Tuple[str, str]]:
        """
        Discover PBC processing tasks.

        Looks for MD trajectories with .xtc and .tpr files.

        Args:
            input_dir: Directory containing MD tasks
            **kwargs: Additional parameters

        Returns:
            Dictionary of task_name -> (trajectory, topology) pairs
        """
        input_path = Path(input_dir)
        tasks = {}

        # Search for MD task directories
        for item in input_path.iterdir():
            if not item.is_dir():
                continue

            # Look for prod directory first
            prod_dir = item / "prod"
            search_dirs = [prod_dir] if prod_dir.exists() else [item]

            for search_dir in search_dirs:
                trajectory_file = None
                topology_file = None

                # Find trajectory file
                xtc_files = list(search_dir.glob("md.xtc"))
                if xtc_files:
                    trajectory_file = str(xtc_files[0])

                # Find topology file
                tpr_files = list(search_dir.glob("md.tpr"))
                if tpr_files:
                    topology_file = str(tpr_files[0])

                # If found both files, add task
                if trajectory_file and topology_file:
                    task_name = item.name
                    tasks[task_name] = (trajectory_file, topology_file)
                    logger.debug(f"Found PBC task: {task_name}")
                    break

        logger.info(f"Discovered {len(tasks)} PBC processing tasks")

        return tasks

    def build_task_config(self,
                         task_name: str,
                         task_data: Tuple[str, str],
                         output_dir: str,
                         **kwargs) -> Dict[str, Any]:
        """
        Build configuration for a single PBC task.

        Args:
            task_name: Task identifier
            task_data: Tuple of (trajectory, topology)
            output_dir: Output directory
            **kwargs: Additional parameters

        Returns:
            Task configuration dictionary
        """
        trajectory, topology = task_data
        dt = kwargs.get('dt')

        task_output_dir = f"{output_dir}/{task_name}_processed"

        config = {
            "trajectory": trajectory,
            "topology": topology,
            "output_dir": task_output_dir,
            "dt": dt
        }

        return config

    def generate_task_command(self,
                             task_name: str,
                             task_data: Tuple[str, str],
                             output_dir: str,
                             **kwargs) -> str:
        """
        Generate PBC processing command for a single task.

        This method is kept for backward compatibility but not used
        in the JSON config approach.

        Args:
            task_name: Task identifier
            task_data: Tuple of (trajectory, topology)
            output_dir: Output directory
            **kwargs: Additional parameters

        Returns:
            Command string (deprecated)
        """
        return f"# Task {task_name} - using JSON config instead"
