#!/usr/bin/env python3
"""
SLURM batch processor for RMSD calculations.

This module generates SLURM scripts for high-throughput RMSD calculations
on processed MD trajectories.
"""

from pathlib import Path
from typing import Dict, Any, Optional, Tuple
import logging
from .slurm_batch_base import SlurmBatchProcessor

logger = logging.getLogger(__name__)


class RMSDBatchProcessor(SlurmBatchProcessor):
    """SLURM batch processor for RMSD calculations."""

    def __init__(self,
                 template_file: Optional[str] = None,
                 default_slurm_params: Optional[Dict[str, Any]] = None):
        """
        Initialize RMSD batch processor.

        Args:
            template_file: Path to custom SLURM template file
            default_slurm_params: Default SLURM parameters
        """
        super().__init__(template_file, default_slurm_params)

    def discover_tasks(self, input_dir: str, **kwargs) -> Dict[str, Tuple[str, str]]:
        """
        Discover RMSD calculation tasks.

        Looks for processed trajectories with both .xtc and topology files.

        Args:
            input_dir: Directory containing processed trajectories
            **kwargs: Additional parameters
                - require_tpr: Require .tpr files (default: True)
                - processed_suffix: Suffix for processed trajectories (default: '_processed')

        Returns:
            Dictionary of task_name -> (trajectory, topology) pairs
        """
        require_tpr = kwargs.get('require_tpr', True)
        processed_suffix = kwargs.get('processed_suffix', '_processed')

        input_path = Path(input_dir)
        tasks = {}

        # Search for processed trajectory directories
        for item in input_path.iterdir():
            if not item.is_dir():
                continue

            # Look for processed trajectories
            trajectory_file = None
            topology_file = None

            # Find trajectory file (prioritize *_processed.xtc)
            for xtc_pattern in ['*_processed.xtc', '*.xtc']:
                xtc_files = list(item.glob(xtc_pattern))
                if xtc_files:
                    trajectory_file = str(xtc_files[0])
                    break

            if not trajectory_file:
                continue

            # Find topology file (prioritize .tpr > .gro > .pdb)
            for topo_ext in ['.tpr', '.gro', '.pdb']:
                topo_files = list(item.glob(f'md{topo_ext}'))
                if topo_files:
                    topology_file = str(topo_files[0])
                    break

            # If md.* not found, look for any topology file
            if not topology_file:
                for topo_ext in ['.tpr', '.gro', '.pdb']:
                    topo_files = list(item.glob(f'*{topo_ext}'))
                    if topo_files:
                        topology_file = str(topo_files[0])
                        break

            # Validate task
            if trajectory_file and topology_file:
                # Check if .tpr is required
                if require_tpr and not topology_file.endswith('.tpr'):
                    logger.warning(f"Skipping {item.name}: no .tpr file found (has {Path(topology_file).suffix})")
                    continue

                task_name = item.name
                tasks[task_name] = (trajectory_file, topology_file)
                logger.debug(f"Found RMSD task: {task_name}")

        logger.info(f"Discovered {len(tasks)} RMSD calculation tasks")

        return tasks

    def build_task_config(self,
                         task_name: str,
                         task_data: Tuple[str, str],
                         output_dir: str,
                         **kwargs) -> Dict[str, Any]:
        """
        Build configuration for a single RMSD task.

        Args:
            task_name: Task identifier
            task_data: Tuple of (trajectory, topology)
            output_dir: Output directory
            **kwargs: Additional parameters

        Returns:
            Task configuration dictionary
        """
        trajectory, topology = task_data
        rmsd_type = kwargs.get('rmsd_type', 'calpha')
        fit_group = kwargs.get('fit_group')
        calc_group = kwargs.get('calc_group')

        output_file = f"{output_dir}/{task_name}_rmsd.xvg"

        config = {
            "trajectory": trajectory,
            "topology": topology,
            "output_file": output_file,
            "rmsd_type": rmsd_type
        }

        if fit_group:
            config["fit_group"] = fit_group
        if calc_group:
            config["calc_group"] = calc_group

        return config

    def generate_task_command(self,
                             task_name: str,
                             task_data: Tuple[str, str],
                             output_dir: str,
                             **kwargs) -> str:
        """
        Generate RMSD calculation command for a single task.

        Args:
            task_name: Task identifier
            task_data: Tuple of (trajectory, topology)
            output_dir: Output directory for RMSD results
            **kwargs: Additional parameters
                - rmsd_type: Type of RMSD (default: 'calpha')
                - fit_group: GROMACS group for fitting
                - calc_group: GROMACS group for calculation

        Returns:
            Command string to calculate RMSD
        """
        trajectory, topology = task_data
        rmsd_type = kwargs.get('rmsd_type', 'calpha')
        fit_group = kwargs.get('fit_group')
        calc_group = kwargs.get('calc_group')

        # Create output filename
        output_file = f"{output_dir}/{task_name}_rmsd.xvg"

        # Build RMSD calculation command
        if fit_group and calc_group:
            # Custom groups
            command = f"""
# Task: {task_name} (custom groups: fit={fit_group}, calc={calc_group})
echo "Calculating RMSD for {task_name}..."
echo "{fit_group}" | gmx rms \\
    -s "{topology}" \\
    -f "{trajectory}" \\
    -o "{output_file}" \\
    -tu ns

if [ $? -eq 0 ]; then
    echo "  RMSD calculation completed: {task_name}"
else
    echo "  RMSD calculation failed: {task_name}"
fi
"""
        else:
            # Use Immunex python interface for intelligent group selection
            command = f"""
# Task: {task_name} (type: {rmsd_type})
echo "Calculating RMSD for {task_name}..."
python -c "
import sys
sys.path.append('/public/home/xmy/tools/AfrerMD')
from immunex.analysis.trajectory.rmsd import RMSDCalculator
import logging

logging.basicConfig(level=logging.INFO)

try:
    calc = RMSDCalculator('{topology}', '{trajectory}')
    result = calc.calculate_gromacs(
        rmsd_type='{rmsd_type}',
        output_file='{output_file}'
    )
    print(f'RMSD completed: {{result}}')
except Exception as e:
    print(f'RMSD failed: {{e}}')
    sys.exit(1)
"

if [ $? -eq 0 ]; then
    echo "  RMSD calculation completed: {task_name}"
else
    echo "  RMSD calculation failed: {task_name}"
fi
"""

        return command.strip()


def generate_rmsd_slurm_scripts(input_dir: str,
                                tasks_per_batch: int = 10,
                                output_dir: str = "./rmsd_results",
                                output_script_dir: str = "./slurm_scripts",
                                rmsd_type: str = "calpha",
                                fit_group: Optional[str] = None,
                                calc_group: Optional[str] = None,
                                require_tpr: bool = True,
                                slurm_params: Optional[Dict[str, Any]] = None,
                                template_file: Optional[str] = None) -> Dict[str, Any]:
    """
    Generate SLURM scripts for batch RMSD calculations.

    Args:
        input_dir: Directory containing processed trajectories
        tasks_per_batch: Number of RMSD calculations per SLURM job
        output_dir: Output directory for RMSD results
        output_script_dir: Directory to save SLURM scripts
        rmsd_type: Type of RMSD calculation ('calpha', 'backbone', etc.)
        fit_group: Custom GROMACS group for fitting
        calc_group: Custom GROMACS group for calculation
        require_tpr: Require .tpr files for RMSD reference
        slurm_params: Custom SLURM parameters
        template_file: Custom SLURM template file

    Returns:
        Results dictionary with script paths and statistics

    Example:
        >>> results = generate_rmsd_slurm_scripts(
        ...     input_dir="/data/processed_trajectories",
        ...     tasks_per_batch=20,
        ...     output_dir="./rmsd_results",
        ...     rmsd_type="calpha",
        ...     slurm_params={"partition": "quick", "cpus_per_task": 8}
        ... )
        >>> print(f"Generated {results['num_batches']} SLURM scripts")
    """
    processor = RMSDBatchProcessor(
        template_file=template_file,
        default_slurm_params=slurm_params
    )

    results = processor.process(
        input_dir=input_dir,
        output_dir=output_dir,
        tasks_per_batch=tasks_per_batch,
        output_script_dir=output_script_dir,
        task_type="rmsd",
        slurm_params=slurm_params,
        script_prefix="rmsd_batch",
        # Task-specific parameters
        rmsd_type=rmsd_type,
        fit_group=fit_group,
        calc_group=calc_group,
        require_tpr=require_tpr
    )

    return results


def main():
    """Command line interface for RMSD SLURM script generation."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Generate SLURM scripts for batch RMSD calculations",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage with C-alpha RMSD
  python -m immunex.utils.slurm_rmsd_generator /data/processed_trajectories

  # Custom batch size and partition
  python -m immunex.utils.slurm_rmsd_generator /data/processed \\
      --tasks-per-batch 20 --partition quick --cpus-per-task 8

  # Use custom GROMACS groups
  python -m immunex.utils.slurm_rmsd_generator /data/processed \\
      --fit-group 4 --calc-group 3

  # Backbone RMSD
  python -m immunex.utils.slurm_rmsd_generator /data/processed \\
      --rmsd-type backbone
        """
    )

    parser.add_argument(
        "input_dir",
        help="Directory containing processed trajectories"
    )

    parser.add_argument(
        "--tasks-per-batch", "-b",
        type=int,
        default=10,
        help="Number of RMSD calculations per SLURM job (default: 10)"
    )

    parser.add_argument(
        "--output-dir", "-o",
        default="./rmsd_results",
        help="Output directory for RMSD results (default: ./rmsd_results)"
    )

    parser.add_argument(
        "--output-scripts",
        default="./slurm_scripts",
        help="Directory to save SLURM scripts (default: ./slurm_scripts)"
    )

    parser.add_argument(
        "--rmsd-type", "-t",
        choices=["calpha", "backbone", "protein", "heavy"],
        default="calpha",
        help="Type of RMSD calculation (default: calpha)"
    )

    parser.add_argument(
        "--fit-group",
        help="Custom GROMACS group for fitting (e.g., 4 for Backbone)"
    )

    parser.add_argument(
        "--calc-group",
        help="Custom GROMACS group for calculation (e.g., 3 for C-alpha)"
    )

    parser.add_argument(
        "--allow-gro",
        action="store_true",
        help="Allow .gro files as topology (default: require .tpr)"
    )

    # SLURM parameters
    parser.add_argument("--partition", "-p", help="SLURM partition")
    parser.add_argument("--time", help="Job time limit")
    parser.add_argument("--cpus-per-task", type=int, help="CPUs per task")
    parser.add_argument("--memory", help="Memory requirement")

    parser.add_argument(
        "--template",
        help="Custom SLURM template file"
    )

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

    # Validate custom groups
    if bool(args.fit_group) != bool(args.calc_group):
        logger.error("Both --fit-group and --calc-group must be specified together")
        return 1

    # Build SLURM parameters
    slurm_params = {}
    if args.partition:
        slurm_params['partition'] = args.partition
    if args.time:
        slurm_params['time'] = args.time
    if args.cpus_per_task:
        slurm_params['cpus_per_task'] = args.cpus_per_task
    if args.memory:
        slurm_params['memory'] = args.memory

    try:
        results = generate_rmsd_slurm_scripts(
            input_dir=args.input_dir,
            tasks_per_batch=args.tasks_per_batch,
            output_dir=args.output_dir,
            output_script_dir=args.output_scripts,
            rmsd_type=args.rmsd_type,
            fit_group=args.fit_group,
            calc_group=args.calc_group,
            require_tpr=not args.allow_gro,
            slurm_params=slurm_params if slurm_params else None,
            template_file=args.template
        )

        if results['success']:
            print()
            print("RMSD SLURM Script Generation Results:")
            print("=" * 50)
            print(f"Total RMSD tasks: {results['total_tasks']}")
            print(f"Tasks per batch: {results['tasks_per_batch']}")
            print(f"Number of batches: {results['num_batches']}")
            print(f"Generated scripts: {len(results['scripts'])}")
            print(f"Script directory: {results['output_script_dir']}")
            print(f"Submission script: {results['submission_script']}")
            print(f"Output directory: {results['output_dir']}")
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
