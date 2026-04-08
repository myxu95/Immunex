#!/usr/bin/env python3
"""
Base class for SLURM batch task processing.

This module provides a generic framework for generating SLURM scripts
to process large batches of tasks on HPC clusters.
"""

import os
import math
from pathlib import Path
from typing import List, Dict, Any, Optional, Tuple
from abc import ABC, abstractmethod
import logging
from datetime import datetime

logger = logging.getLogger(__name__)


class SlurmBatchProcessor(ABC):
    """
    Abstract base class for SLURM batch task processing.

    This class handles the common logic of:
    - Task partitioning into batches
    - SLURM script generation
    - Submission script creation

    Subclasses need to implement:
    - discover_tasks(): Find tasks to process
    - generate_task_command(): Generate command for a single task
    """

    def __init__(self,
                 template_file: Optional[str] = None,
                 default_slurm_params: Optional[Dict[str, Any]] = None):
        """
        Initialize SLURM batch processor.

        Args:
            template_file: Path to custom SLURM template file
            default_slurm_params: Default SLURM parameters
        """
        self.template_file = template_file
        self.default_params = default_slurm_params or self._get_default_slurm_params()

        # Load template
        if template_file and Path(template_file).exists():
            self.template = self._load_template(template_file)
        else:
            self.template = self._get_default_template()

    def _get_default_slurm_params(self) -> Dict[str, Any]:
        """Get default SLURM parameters."""
        return {
            "partition": "quick",
            "time": "12:00:00",
            "nodes": 1,
            "ntasks": 1,
            "cpus_per_task": 11,
            "gres": "gpu:1",
            "memory": None,
            "conda_env": "immunex",
            "gromacs_module": "/public/software/profile.d/apps_gromacs_2023.2.sh",
            "ld_library_path": "/public/software/lib/"
        }

    def _load_template(self, template_file: str) -> str:
        """Load SLURM template from file."""
        try:
            with open(template_file, 'r') as f:
                template = f.read()
            logger.info(f"Loaded SLURM template from: {template_file}")
            return template
        except Exception as e:
            logger.warning(f"Failed to load template file: {e}")
            return self._get_default_template()

    def _get_default_template(self) -> str:
        """Get default SLURM template."""
        return """#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH -p {partition}
#SBATCH --time={time}
#SBATCH --nodes={nodes}
#SBATCH --ntasks={ntasks}
#SBATCH --cpus-per-task={cpus_per_task}
{gres_line}
{memory_line}

# Environment setup
export LD_LIBRARY_PATH={ld_library_path}:$LD_LIBRARY_PATH
source {gromacs_module}

# Activate conda environment
source /public/home/xmy/app/miniconda/bin/activate {conda_env}

# Job information
echo "Start time: $(date)"
echo "SLURM_JOB_NODELIST: $SLURM_JOB_NODELIST"
echo "hostname: $(hostname)"
echo "CUDA_VISIBLE_DEVICES: $CUDA_VISIBLE_DEVICES"
echo "Job directory: $(pwd)"
echo "Processing {task_count} tasks in this batch"
echo "Config file: {config_file}"
echo ""

# Run batch worker
python {batch_worker_script} --config {config_file} {log_option}

echo ""
echo "End time: $(date)"
echo "Batch processing completed"
"""

    @abstractmethod
    def discover_tasks(self, input_dir: str, **kwargs) -> Dict[str, Any]:
        """
        Discover tasks to process.

        Args:
            input_dir: Input directory containing tasks
            **kwargs: Additional discovery parameters

        Returns:
            Dictionary of task_name -> task_data mappings
        """
        pass

    @abstractmethod
    def generate_task_command(self,
                             task_name: str,
                             task_data: Any,
                             output_dir: str,
                             **kwargs) -> str:
        """
        Generate command to process a single task.

        Args:
            task_name: Task identifier
            task_data: Task-specific data
            output_dir: Output directory for results
            **kwargs: Additional task parameters

        Returns:
            Command string to execute for this task
        """
        pass

    def partition_tasks(self,
                       tasks: Dict[str, Any],
                       tasks_per_batch: int) -> List[Dict[str, Any]]:
        """
        Partition tasks into batches for SLURM processing.

        Args:
            tasks: Dictionary of all tasks
            tasks_per_batch: Number of tasks per SLURM job

        Returns:
            List of task batches
        """
        if not tasks:
            logger.warning("No tasks provided for partitioning")
            return []

        if tasks_per_batch <= 0:
            raise ValueError("tasks_per_batch must be greater than 0")

        task_items = list(tasks.items())
        total_tasks = len(task_items)
        num_batches = math.ceil(total_tasks / tasks_per_batch)

        batches = []
        for i in range(num_batches):
            start_idx = i * tasks_per_batch
            end_idx = min(start_idx + tasks_per_batch, total_tasks)
            batch = dict(task_items[start_idx:end_idx])
            batches.append(batch)

        logger.info(f"Partitioned {total_tasks} tasks into {num_batches} batches "
                   f"({tasks_per_batch} tasks per batch)")

        return batches

    def generate_batch_config(self,
                             batch_id: int,
                             batch_tasks: Dict[str, Any],
                             output_dir: str,
                             task_type: str,
                             config_dir: str,
                             **kwargs) -> str:
        """
        Generate JSON configuration file for a batch of tasks.

        Args:
            batch_id: Batch identifier
            batch_tasks: Dictionary of tasks in this batch
            output_dir: Output directory
            task_type: Task type (pbc, rmsd, etc.)
            config_dir: Directory to save config file
            **kwargs: Additional task-specific parameters

        Returns:
            Path to generated config file
        """
        import json

        config_path = Path(config_dir)
        config_path.mkdir(parents=True, exist_ok=True)

        config_file = config_path / f"batch_{batch_id:03d}_config.json"

        # Build task configurations
        task_configs = {}
        for task_name, task_data in batch_tasks.items():
            task_configs[task_name] = self.build_task_config(
                task_name, task_data, output_dir, **kwargs
            )

        # Build full configuration
        config = {
            "task_type": task_type,
            "batch_id": batch_id,
            "output_dir": output_dir,
            "tasks": task_configs
        }

        # Add task-specific global parameters
        for key, value in kwargs.items():
            if key not in config and not key.startswith('_'):
                config[key] = value

        # Write config file
        with open(config_file, 'w') as f:
            json.dump(config, f, indent=2)

        logger.info(f"Generated config file: {config_file}")

        return str(config_file)

    @abstractmethod
    def build_task_config(self,
                         task_name: str,
                         task_data: Any,
                         output_dir: str,
                         **kwargs) -> Dict[str, Any]:
        """
        Build configuration for a single task.

        Args:
            task_name: Task identifier
            task_data: Task-specific data
            output_dir: Output directory
            **kwargs: Additional parameters

        Returns:
            Task configuration dictionary
        """
        pass

    def _generate_job_name(self,
                          input_dir: str,
                          batch_id: int,
                          total_batches: int,
                          task_type: str = "task") -> str:
        """
        Generate SLURM job name.

        Args:
            input_dir: Input directory path
            batch_id: Batch identifier
            total_batches: Total number of batches
            task_type: Task type identifier (e.g., 'pbc', 'rmsd')

        Returns:
            Job name string
        """
        input_path = Path(input_dir)
        dataset_name = input_path.name

        clean_dataset = "".join(c for c in dataset_name if c.isalnum() or c in "_")
        clean_dataset = clean_dataset[:12]

        job_name = f"amd_{task_type}_{clean_dataset}_{batch_id}of{total_batches}"

        if len(job_name) > 24:
            max_dataset_len = 24 - len(f"amd_{task_type}__{batch_id}of{total_batches}")
            if max_dataset_len > 0:
                clean_dataset = clean_dataset[:max_dataset_len]
                job_name = f"amd_{task_type}_{clean_dataset}_{batch_id}of{total_batches}"
            else:
                job_name = f"amd_{task_type}_{batch_id}of{total_batches}"

        return job_name

    def generate_slurm_script(self,
                             batch_id: int,
                             batch_tasks: Dict[str, Any],
                             input_dir: str,
                             output_dir: str,
                             total_batches: int,
                             task_type: str = "task",
                             config_dir: str = "./configs",
                             slurm_params: Optional[Dict[str, Any]] = None,
                             **kwargs) -> Tuple[str, str]:
        """
        Generate a single SLURM script for a batch of tasks.

        Args:
            batch_id: Batch identifier
            batch_tasks: Dictionary of tasks in this batch
            input_dir: Input directory
            output_dir: Output directory
            total_batches: Total number of batches
            task_type: Task type identifier
            config_dir: Directory to save config files
            slurm_params: Custom SLURM parameters
            **kwargs: Additional parameters for command generation

        Returns:
            Tuple of (SLURM script content, config file path)
        """
        params = self.default_params.copy()
        if slurm_params:
            params.update(slurm_params)

        job_name = params.get('job_name') or self._generate_job_name(
            input_dir, batch_id, total_batches, task_type
        )

        # Generate JSON config file
        config_file = self.generate_batch_config(
            batch_id=batch_id,
            batch_tasks=batch_tasks,
            output_dir=output_dir,
            task_type=task_type,
            config_dir=config_dir,
            **kwargs
        )

        # Batch worker script path
        batch_worker_script = params.get(
            'batch_worker_script',
            '/public/home/xmy/tools/AfrerMD/scripts/batch_worker.py'
        )

        # Optional log file
        log_file = f"{config_dir}/batch_{batch_id:03d}.log"
        log_option = f"--log-file {log_file}"

        template_vars = {
            'job_name': job_name,
            'partition': params['partition'],
            'time': params['time'],
            'nodes': params['nodes'],
            'ntasks': params['ntasks'],
            'cpus_per_task': params['cpus_per_task'],
            'gres_line': f"#SBATCH --gres={params['gres']}" if params.get('gres') else "",
            'memory_line': f"#SBATCH --mem={params['memory']}" if params.get('memory') else "",
            'ld_library_path': params['ld_library_path'],
            'gromacs_module': params['gromacs_module'],
            'conda_env': params.get('conda_env', 'immunex'),
            'task_count': len(batch_tasks),
            'config_file': config_file,
            'batch_worker_script': batch_worker_script,
            'log_option': log_option
        }

        script_content = self.template.format(**template_vars)
        return script_content, config_file

    def generate_batch_scripts(self,
                              tasks: Dict[str, Any],
                              input_dir: str,
                              output_dir: str,
                              tasks_per_batch: int,
                              output_script_dir: str,
                              task_type: str = "task",
                              slurm_params: Optional[Dict[str, Any]] = None,
                              script_prefix: str = "immunex_batch",
                              **kwargs) -> List[str]:
        """
        Generate multiple SLURM scripts for batch processing.

        Args:
            tasks: All tasks to process
            input_dir: Input directory
            output_dir: Output directory
            tasks_per_batch: Number of tasks per SLURM job
            output_script_dir: Directory to save generated scripts
            task_type: Task type identifier
            slurm_params: Custom SLURM parameters
            script_prefix: Prefix for script filenames
            **kwargs: Additional parameters

        Returns:
            List of generated script file paths
        """
        if not tasks:
            logger.error("No tasks provided for script generation")
            return []

        script_dir = Path(output_script_dir)
        script_dir.mkdir(parents=True, exist_ok=True)

        # Config directory alongside scripts
        config_dir = script_dir / "configs"
        config_dir.mkdir(parents=True, exist_ok=True)

        task_batches = self.partition_tasks(tasks, tasks_per_batch)
        total_batches = len(task_batches)

        generated_scripts = []

        for batch_id, batch_tasks in enumerate(task_batches, 1):
            script_content, config_file = self.generate_slurm_script(
                batch_id=batch_id,
                batch_tasks=batch_tasks,
                input_dir=input_dir,
                output_dir=output_dir,
                total_batches=total_batches,
                task_type=task_type,
                config_dir=str(config_dir),
                slurm_params=slurm_params,
                **kwargs
            )

            script_filename = f"{script_prefix}_{batch_id:03d}.sh"
            script_path = script_dir / script_filename

            try:
                with open(script_path, 'w') as f:
                    f.write(script_content)

                os.chmod(script_path, 0o755)
                generated_scripts.append(str(script_path))

                task_names = list(batch_tasks.keys())
                logger.info(f"Generated {script_filename} for {len(task_names)} tasks")
                logger.info(f"  Config: {config_file}")

            except Exception as e:
                logger.error(f"Failed to write script {script_filename}: {e}")

        logger.info(f"Generated {len(generated_scripts)} SLURM scripts in: {output_script_dir}")

        return generated_scripts

    def generate_submission_script(self,
                                  script_paths: List[str],
                                  output_dir: str,
                                  submit_delay: int = 5) -> str:
        """
        Generate a master script to submit all batch jobs.

        Args:
            script_paths: List of generated SLURM script paths
            output_dir: Directory to save the submission script
            submit_delay: Delay between job submissions (seconds)

        Returns:
            Path to the submission script
        """
        submit_script_path = Path(output_dir) / "submit_all_batches.sh"

        submit_content = f"""#!/bin/bash
# Immunex Batch Job Submission Script
# Generated on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

SCRIPT_DIR="$(cd "$(dirname "${{BASH_SOURCE[0]}}")" && pwd)"
cd "$SCRIPT_DIR"

echo "Current directory: $(pwd)"
echo "Submitting {len(script_paths)} batch jobs..."
echo ""

"""

        for i, script_path in enumerate(script_paths, 1):
            script_name = Path(script_path).name
            script_path_obj = Path(script_path)
            output_dir_obj = Path(output_dir)

            try:
                relative_path = script_path_obj.relative_to(output_dir_obj)
                script_ref = str(relative_path)
            except ValueError:
                script_ref = str(script_path_obj.resolve())

            submit_content += f"""# Submit batch {i}/{len(script_paths)}
echo "Submitting {script_name}..."
if [ -f "{script_ref}" ]; then
    sbatch "{script_ref}"
    if [ $? -eq 0 ]; then
        echo "  Successfully submitted {script_name}"
    else
        echo "  Failed to submit {script_name}"
    fi
else
    echo "  Script not found: {script_ref}"
fi
"""

            if i < len(script_paths) and submit_delay > 0:
                submit_content += f"sleep {submit_delay}\n"

            submit_content += "\n"

        submit_content += f"""echo "All {len(script_paths)} batch jobs submitted!"
echo "Use 'squeue -u $USER' to check job status"
"""

        with open(submit_script_path, 'w') as f:
            f.write(submit_content)

        os.chmod(submit_script_path, 0o755)

        logger.info(f"Generated submission script: {submit_script_path}")

        return str(submit_script_path)

    def process(self,
                input_dir: str,
                output_dir: str,
                tasks_per_batch: int,
                output_script_dir: str,
                task_type: str = "task",
                slurm_params: Optional[Dict[str, Any]] = None,
                script_prefix: str = "immunex_batch",
                **kwargs) -> Dict[str, Any]:
        """
        Main processing method: discover tasks, generate scripts, create submission script.

        Args:
            input_dir: Input directory
            output_dir: Output directory
            tasks_per_batch: Number of tasks per batch
            output_script_dir: Script output directory
            task_type: Task type identifier
            slurm_params: SLURM parameters
            script_prefix: Script filename prefix
            **kwargs: Additional parameters

        Returns:
            Results dictionary
        """
        try:
            logger.info(f"Discovering {task_type} tasks in: {input_dir}")
            tasks = self.discover_tasks(input_dir, **kwargs)

            if not tasks:
                return {
                    "success": False,
                    "error": f"No {task_type} tasks found",
                    "total_tasks": 0,
                    "scripts": [],
                    "submission_script": None
                }

            logger.info(f"Generating SLURM scripts for {len(tasks)} tasks...")
            script_paths = self.generate_batch_scripts(
                tasks=tasks,
                input_dir=input_dir,
                output_dir=output_dir,
                tasks_per_batch=tasks_per_batch,
                output_script_dir=output_script_dir,
                task_type=task_type,
                slurm_params=slurm_params,
                script_prefix=script_prefix,
                **kwargs
            )

            submission_script = self.generate_submission_script(
                script_paths=script_paths,
                output_dir=output_script_dir
            )

            return {
                "success": True,
                "total_tasks": len(tasks),
                "tasks_per_batch": tasks_per_batch,
                "num_batches": len(script_paths),
                "scripts": script_paths,
                "submission_script": submission_script,
                "output_script_dir": output_script_dir,
                "output_dir": output_dir
            }

        except Exception as e:
            logger.error(f"Processing failed: {e}")
            return {
                "success": False,
                "error": str(e),
                "total_tasks": 0,
                "scripts": [],
                "submission_script": None
            }
