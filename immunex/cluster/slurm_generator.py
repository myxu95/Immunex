#!/usr/bin/env python3
"""
SLURM Script Generator for Immunex Batch Processing

This module generates SLURM job scripts for processing multiple MD simulation
tasks on HPC clusters. It automatically partitions large numbers of MD tasks
into manageable chunks and creates individual SLURM scripts for each chunk.

Features:
- Automatic task partitioning based on user-specified batch size
- Customizable SLURM parameters (nodes, CPUs, time, partition, etc.)
- Integration with Immunex batch processing workflow
- Support for conda environment activation
- Flexible template system for different cluster configurations
"""

import os
import math
from pathlib import Path
from typing import List, Dict, Any, Optional, Tuple
import logging
from datetime import datetime

logger = logging.getLogger(__name__)


class SlurmScriptGenerator:
    """Generate SLURM job scripts for batch MD processing."""
    
    def __init__(self, 
                 template_file: Optional[str] = None,
                 default_slurm_params: Optional[Dict[str, Any]] = None):
        """
        Initialize SLURM script generator.
        
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
            "memory": None,  # Use default
            "conda_env": None,  # No conda environment by default
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
source /public/home/xmy/app/miniconda/bin/activate immunex

# Job information
echo "Start time: $(date)"
echo "SLURM_JOB_NODELIST: $SLURM_JOB_NODELIST"
echo "hostname: $(hostname)"
echo "CUDA_VISIBLE_DEVICES: $CUDA_VISIBLE_DEVICES"
echo "Job directory: $(pwd)"
echo "Processing {task_count} MD tasks in this batch"
echo ""

# Immunex batch processing
{processing_commands}

echo ""
echo "End time: $(date)"
echo "Batch processing completed"
"""
    
    def partition_tasks(self, 
                       md_tasks: Dict[str, Tuple[str, str]], 
                       tasks_per_batch: int) -> List[Dict[str, Tuple[str, str]]]:
        """
        Partition MD tasks into batches for SLURM processing.
        
        Args:
            md_tasks: Dictionary of task_name -> (trajectory, topology) pairs
            tasks_per_batch: Number of tasks per SLURM job
            
        Returns:
            List of task batches (each batch is a dict of tasks)
        """
        if not md_tasks:
            logger.warning("No MD tasks provided for partitioning")
            return []
        
        if tasks_per_batch <= 0:
            raise ValueError("tasks_per_batch must be greater than 0")
        
        # Convert to list for easier partitioning
        task_items = list(md_tasks.items())
        total_tasks = len(task_items)
        
        # Calculate number of batches needed
        num_batches = math.ceil(total_tasks / tasks_per_batch)
        
        batches = []
        for i in range(num_batches):
            start_idx = i * tasks_per_batch
            end_idx = min(start_idx + tasks_per_batch, total_tasks)
            
            # Create batch dictionary
            batch = dict(task_items[start_idx:end_idx])
            batches.append(batch)
        
        logger.info(f"Partitioned {total_tasks} tasks into {num_batches} batches "
                   f"({tasks_per_batch} tasks per batch)")
        
        return batches
    
    def generate_processing_commands(self, 
                                   batch_tasks: Dict[str, Tuple[str, str]],
                                   base_input_dir: str,
                                   base_output_dir: str,
                                   dt: Optional[float] = None,
                                   use_python_module: bool = True) -> str:
        """
        Generate Immunex processing commands for a batch of tasks.
        
        Args:
            batch_tasks: Dictionary of tasks in this batch
            base_input_dir: Base input directory containing all MD tasks
            base_output_dir: Base output directory for processed results
            dt: Time interval for frame sampling (optional)
            use_python_module: Use python module interface vs direct python script
            
        Returns:
            String containing the processing commands
        """
        if use_python_module:
            # Use specific task processing approach
            task_names = list(batch_tasks.keys())

            # Create python script that processes only specific tasks
            command = f"""python -c "
from immunex.core import discover_tasks_from_list
from immunex.pipeline import BatchExecutor, PreprocessOnlyPipeline
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

batch_tasks = {repr(batch_tasks)}
base_output_dir = '{base_output_dir}'

task_list = [
    {{
        'task_id': task_name,
        'topology': topology,
        'trajectory_raw': trajectory,
    }}
    for task_name, (trajectory, topology) in batch_tasks.items()
]

report = discover_tasks_from_list(
    task_list,
    required_files=['topology', 'trajectory'],
)
results = BatchExecutor(max_workers=len(batch_tasks)).execute_pipeline(
    report,
    PreprocessOnlyPipeline(dt={dt if dt else 'None'}),
    show_progress=False,
    output_base_dir=base_output_dir,
)

successful = len([result for result in results if not result.has_errors()])
failed = len(results) - successful
logger.info(f'Batch completed: {{successful}} successful, {{failed}} failed')
" """
            
        else:
            # Use direct python script approach
            task_names = list(batch_tasks.keys())
            task_list_str = " ".join(f'"{name}"' for name in task_names)
            
            command = f"""python -c "
from immunex import process_md_tasks

# Process specific tasks in this batch
task_names = [{', '.join(repr(name) for name in task_names)}]
print(f'Processing batch with tasks: {{task_names}}')

results = process_md_tasks(
    '{base_input_dir}',
    output_dir='{base_output_dir}',
    {'dt=' + str(dt) + ',' if dt else ''}
    max_workers={len(batch_tasks)}
)

print(f'Batch completed: {{results[\"successful\"]}}/{{results[\"total_tasks\"]}} tasks successful')
" """
        
        return command
    
    def _generate_simple_job_name(self, 
                                 base_input_dir: str, 
                                 batch_id: int, 
                                 total_batches: int) -> str:
        """
        Generate simple, recognizable job name for SLURM queue viewing.
        
        Args:
            base_input_dir: Base input directory path
            batch_id: Batch identifier
            total_batches: Total number of batches
            
        Returns:
            Simple job name string (format: amd_dataset_1of3)
        """
        # Extract dataset name from input directory
        input_path = Path(base_input_dir)
        dataset_name = input_path.name
        
        # Clean and shorten dataset name
        clean_dataset = "".join(c for c in dataset_name if c.isalnum() or c in "_")
        # Take first 12 characters to keep total length reasonable
        clean_dataset = clean_dataset[:12]
        
        # Generate simple job name: amd_dataset_XofY
        job_name = f"amd_{clean_dataset}_{batch_id}of{total_batches}"
        
        # Ensure reasonable length (should be < 25 characters for good squeue display)
        if len(job_name) > 24:
            # Further truncate dataset name if needed
            max_dataset_len = 24 - len(f"amd__{batch_id}of{total_batches}")
            if max_dataset_len > 0:
                clean_dataset = clean_dataset[:max_dataset_len]
                job_name = f"amd_{clean_dataset}_{batch_id}of{total_batches}"
            else:
                # Ultimate fallback
                job_name = f"amd_{batch_id}of{total_batches}"
        
        return job_name
    
    def generate_slurm_script(self, 
                             batch_id: int,
                             batch_tasks: Dict[str, Tuple[str, str]], 
                             base_input_dir: str,
                             base_output_dir: str,
                             total_batches: int,
                             slurm_params: Optional[Dict[str, Any]] = None,
                             dt: Optional[float] = None) -> str:
        """
        Generate a single SLURM script for a batch of tasks.
        
        Args:
            batch_id: Batch identifier number
            batch_tasks: Dictionary of tasks in this batch  
            base_input_dir: Base input directory
            base_output_dir: Base output directory
            total_batches: Total number of batches (for job naming)
            slurm_params: Custom SLURM parameters
            dt: Time sampling interval
            
        Returns:
            SLURM script content as string
        """
        # Merge parameters
        params = self.default_params.copy()
        if slurm_params:
            params.update(slurm_params)
        
        # Generate simple job name
        if 'job_name' in params:
            # Use explicitly provided job name
            job_name = params['job_name']
        else:
            # Generate simple, squeue-friendly job name
            job_name = self._generate_simple_job_name(
                base_input_dir, batch_id, total_batches
            )
        
        # Generate processing commands
        processing_commands = self.generate_processing_commands(
            batch_tasks, base_input_dir, base_output_dir, dt
        )
        
        # Prepare template variables
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
            'task_count': len(batch_tasks),
            'processing_commands': processing_commands
        }
        
        # Handle conda environment activation
        if params.get('conda_env'):
            conda_activation = f"""
# Activate conda environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate {params['conda_env']}
echo "Activated conda environment: {params['conda_env']}"
"""
        else:
            conda_activation = ""
        
        template_vars['conda_activation'] = conda_activation
        
        # Fill in the template
        script_content = self.template.format(**template_vars)
        
        return script_content
    
    def generate_batch_scripts(self, 
                              md_tasks: Dict[str, Tuple[str, str]],
                              base_input_dir: str,
                              base_output_dir: str,
                              tasks_per_batch: int,
                              output_script_dir: str,
                              slurm_params: Optional[Dict[str, Any]] = None,
                              dt: Optional[float] = None,
                              script_prefix: str = "immunex_batch") -> List[str]:
        """
        Generate multiple SLURM scripts for batch processing.
        
        Args:
            md_tasks: All MD tasks to process
            base_input_dir: Base input directory
            base_output_dir: Base output directory
            tasks_per_batch: Number of tasks per SLURM job
            output_script_dir: Directory to save generated scripts
            slurm_params: Custom SLURM parameters
            dt: Time sampling interval
            script_prefix: Prefix for script filenames
            
        Returns:
            List of generated script file paths
        """
        if not md_tasks:
            logger.error("No MD tasks provided for script generation")
            return []
        
        # Create output directory
        script_dir = Path(output_script_dir)
        script_dir.mkdir(parents=True, exist_ok=True)
        
        # Partition tasks into batches
        task_batches = self.partition_tasks(md_tasks, tasks_per_batch)
        total_batches = len(task_batches)
        
        generated_scripts = []
        
        for batch_id, batch_tasks in enumerate(task_batches, 1):
            # Generate script content
            script_content = self.generate_slurm_script(
                batch_id=batch_id,
                batch_tasks=batch_tasks,
                base_input_dir=base_input_dir,
                base_output_dir=base_output_dir,
                total_batches=total_batches,
                slurm_params=slurm_params,
                dt=dt
            )
            
            # Generate script filename
            script_filename = f"{script_prefix}_{batch_id:03d}.sh"
            script_path = script_dir / script_filename
            
            # Write script file
            try:
                with open(script_path, 'w') as f:
                    f.write(script_content)
                
                # Make script executable
                os.chmod(script_path, 0o755)
                
                generated_scripts.append(str(script_path))
                
                # Log batch information
                task_names = list(batch_tasks.keys())
                logger.info(f"Generated script {script_filename} for {len(task_names)} tasks:")
                for task_name in task_names:
                    logger.info(f"  • {task_name}")
                
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

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${{BASH_SOURCE[0]}}")" && pwd)"
cd "$SCRIPT_DIR"

echo "Current directory: $(pwd)"
echo "Submitting {len(script_paths)} Immunex batch jobs..."
echo ""

"""
        
        for i, script_path in enumerate(script_paths, 1):
            script_name = Path(script_path).name
            # Use relative path if script is in same directory as submission script
            script_path_obj = Path(script_path)
            output_dir_obj = Path(output_dir)
            
            try:
                # Try to get relative path from submission script directory
                relative_path = script_path_obj.relative_to(output_dir_obj)
                script_ref = str(relative_path)
            except ValueError:
                # Fall back to absolute path if relative path fails
                script_ref = str(script_path_obj.resolve())
            
            submit_content += f"""# Submit batch {i}/{len(script_paths)}
echo "Submitting {script_name}..."
if [ -f "{script_ref}" ]; then
    sbatch "{script_ref}"
    if [ $? -eq 0 ]; then
        echo "  ✅ Successfully submitted {script_name}"
    else
        echo "  ❌ Failed to submit {script_name}"
    fi
else
    echo "  ❌ Script not found: {script_ref}"
fi
"""
            
            if i < len(script_paths) and submit_delay > 0:
                submit_content += f"sleep {submit_delay}\n"
            
            submit_content += "\n"
        
        submit_content += f"""echo "All {len(script_paths)} batch jobs submitted!"
echo "Use 'squeue -u $USER' to check job status"
"""
        
        # Write submission script
        with open(submit_script_path, 'w') as f:
            f.write(submit_content)
        
        # Make executable
        os.chmod(submit_script_path, 0o755)
        
        logger.info(f"Generated submission script: {submit_script_path}")
        
        return str(submit_script_path)


def _generate_adaptive_output_dir(simulations_path: str, 
                                 custom_suffix: Optional[str] = None) -> str:
    """
    Generate adaptive output directory based on input path.
    
    Args:
        simulations_path: Input simulations directory path
        custom_suffix: Custom suffix for output directory
        
    Returns:
        Adaptive output directory path
    """
    input_path = Path(simulations_path)
    
    if custom_suffix:
        # Use custom suffix
        output_dir = input_path.parent / f"{input_path.name}_{custom_suffix}"
    else:
        # Default naming: add _processed suffix
        output_dir = input_path.parent / f"{input_path.name}_processed"
    
    return str(output_dir)


def generate_slurm_scripts_for_md_tasks(simulations_path: str,
                                       tasks_per_batch: int = 10,
                                       output_script_dir: str = "./slurm_scripts",
                                       base_output_dir: Optional[str] = None,
                                       output_suffix: Optional[str] = None,
                                       slurm_params: Optional[Dict[str, Any]] = None,
                                       dt: Optional[float] = None,
                                       template_file: Optional[str] = None,
                                       qualified_mds: Optional[List] = None) -> Dict[str, Any]:
    """
    Convenience function to generate SLURM scripts for MD batch processing.
    
    Args:
        simulations_path: Path containing MD task folders
        tasks_per_batch: Number of tasks per SLURM job (default: 10)
        output_script_dir: Directory to save generated scripts
        base_output_dir: Output directory for processed results (if None, auto-generated)
        output_suffix: Custom suffix for auto-generated output directory (default: "processed")
        slurm_params: Custom SLURM parameters
        dt: Time sampling interval  
        template_file: Custom SLURM template file
        qualified_mds: List of qualified MD paths (if None, process all MDs)
        
    Returns:
        Dictionary with generation results and file paths
        
    Example:
        >>> # Basic usage with simple job names and output paths
        >>> results = generate_slurm_scripts_for_md_tasks(
        ...     simulations_path="/data/antibody_simulations",
        ...     tasks_per_batch=5,
        ...     output_suffix="pbc_corrected",
        ...     slurm_params={
        ...         "partition": "gpu", 
        ...         "time": "12:00:00"
        ...     }
        ... )
        >>> print(f"Generated {len(results['scripts'])} SLURM scripts")
        >>> # Job names will be: amd_antibody_sim_1of4, amd_antibody_sim_2of4, etc.
        >>> # Output dir: /data/antibody_simulations_pbc_corrected
    """
    from ..pipeline.batch_workflow import discover_md_tasks
    
    # Discover MD tasks
    logger.info(f"Discovering MD tasks in: {simulations_path}")
    all_md_tasks = discover_md_tasks(simulations_path)
    
    # Filter MD tasks based on qualified list if provided
    if qualified_mds is not None:
        logger.info(f"Filtering MD tasks based on qualified list ({len(qualified_mds)} qualified)")
        md_tasks = {}
        qualified_paths = {str(Path(md_path).resolve()) for md_path in qualified_mds}
        
        for task_name, task_files in all_md_tasks.items():
            # Check if this task's directory is in the qualified list
            task_dir = Path(simulations_path) / task_name
            task_path = str(task_dir.resolve())
            
            if task_path in qualified_paths:
                md_tasks[task_name] = task_files
                logger.debug(f"Included qualified MD task: {task_name}")
            else:
                logger.debug(f"Skipped unqualified MD task: {task_name}")
        
        logger.info(f"Filtered to {len(md_tasks)} qualified MD tasks")
    else:
        logger.info("No quality filtering applied - processing all discovered MD tasks")
        md_tasks = all_md_tasks
    
    if not md_tasks:
        error_msg = "No qualified MD tasks found" if qualified_mds else "No valid MD tasks found"
        logger.error(f"{error_msg} for SLURM script generation")
        return {
            "success": False,
            "error": error_msg,
            "scripts": [],
            "submission_script": None
        }
    
    # Set output directory
    if base_output_dir is None:
        # Auto-generate output directory with custom suffix support
        base_output_dir = _generate_adaptive_output_dir(simulations_path, output_suffix)
    
    # Initialize generator
    generator = SlurmScriptGenerator(template_file=template_file)
    
    # Generate SLURM scripts
    try:
        script_paths = generator.generate_batch_scripts(
            md_tasks=md_tasks,
            base_input_dir=simulations_path,
            base_output_dir=base_output_dir,
            tasks_per_batch=tasks_per_batch,
            output_script_dir=output_script_dir,
            slurm_params=slurm_params,
            dt=dt
        )
        
        # Generate submission script
        submission_script = generator.generate_submission_script(
            script_paths=script_paths,
            output_dir=output_script_dir
        )
        
        results = {
            "success": True,
            "total_tasks": len(md_tasks),
            "tasks_per_batch": tasks_per_batch,
            "num_batches": len(script_paths),
            "scripts": script_paths,
            "submission_script": submission_script,
            "output_script_dir": output_script_dir,
            "base_output_dir": base_output_dir
        }
        
        logger.info(f"Successfully generated {len(script_paths)} SLURM scripts")
        logger.info(f"Use '{submission_script}' to submit all jobs")
        
        return results
        
    except Exception as e:
        logger.error(f"Failed to generate SLURM scripts: {e}")
        return {
            "success": False,
            "error": str(e),
            "scripts": [],
            "submission_script": None
        }


# Command line interface
def main():
    """Command line interface for SLURM script generation."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Generate SLURM scripts for Immunex batch processing",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python -m immunex.utils.slurm_generator /data/simulations --tasks-per-batch 5
    python -m immunex.utils.slurm_generator /data/simulations --partition gpu --time 24:00:00
    python -m immunex.utils.slurm_generator /data/simulations --template my_template.sh
        """
    )
    
    parser.add_argument(
        "simulations_path",
        help="Path containing MD task folders"
    )
    
    parser.add_argument(
        "--tasks-per-batch", "-t",
        type=int,
        default=10,
        help="Number of tasks per SLURM job (default: 10)"
    )
    
    parser.add_argument(
        "--output-scripts", "-o",
        default="./slurm_scripts",
        help="Directory to save generated scripts (default: ./slurm_scripts)"
    )
    
    parser.add_argument(
        "--output-data",
        help="Output directory for processed data (default: auto-generated)"
    )
    
    parser.add_argument(
        "--output-suffix",
        default="processed",
        help="Suffix for auto-generated output directory (default: 'processed')"
    )
    
    parser.add_argument(
        "--template",
        help="Custom SLURM template file"
    )
    
    parser.add_argument(
        "--dt",
        type=float,
        help="Time interval for frame sampling in ps"
    )
    
    # SLURM parameters
    parser.add_argument("--partition", help="SLURM partition")
    parser.add_argument("--time", help="Job time limit")
    parser.add_argument("--nodes", type=int, help="Number of nodes")
    parser.add_argument("--cpus-per-task", type=int, help="CPUs per task")
    parser.add_argument("--gres", help="Generic resources (e.g., gpu:1)")
    parser.add_argument("--memory", help="Memory requirement")
    parser.add_argument("--conda-env", help="Conda environment to activate")
    
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Enable verbose logging"
    )
    
    args = parser.parse_args()
    
    # Configure logging
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    
    # Build SLURM parameters
    slurm_params = {}
    if args.partition:
        slurm_params['partition'] = args.partition
    if args.time:
        slurm_params['time'] = args.time
    if args.nodes:
        slurm_params['nodes'] = args.nodes
    if args.cpus_per_task:
        slurm_params['cpus_per_task'] = args.cpus_per_task
    if args.gres:
        slurm_params['gres'] = args.gres
    if args.memory:
        slurm_params['memory'] = args.memory
    if args.conda_env:
        slurm_params['conda_env'] = args.conda_env
    
    try:
        results = generate_slurm_scripts_for_md_tasks(
            simulations_path=args.simulations_path,
            tasks_per_batch=args.tasks_per_batch,
            output_script_dir=args.output_scripts,
            base_output_dir=args.output_data,
            output_suffix=args.output_suffix,
            slurm_params=slurm_params if slurm_params else None,
            dt=args.dt,
            template_file=args.template
        )
        
        if results['success']:
            print()
            print("SLURM Script Generation Results:")
            print("=" * 40)
            print(f"Total MD tasks: {results['total_tasks']}")
            print(f"Tasks per batch: {results['tasks_per_batch']}")
            print(f"Number of batches: {results['num_batches']}")
            print(f"Generated scripts: {len(results['scripts'])}")
            print(f"Script directory: {results['output_script_dir']}")
            print(f"Submission script: {results['submission_script']}")
            print()
            print("To submit all jobs:")
            print(f"  bash {results['submission_script']}")
        else:
            print(f"Error: {results['error']}")
            return 1
    
    except Exception as e:
        logger.error(f"Error: {e}")
        return 1
    
    return 0


if __name__ == "__main__":
    exit(main())
