#!/usr/bin/env python3
"""
RMSD Calculator Script for Batch Trajectory Analysis

This script automatically finds all .xtc trajectory files in a specified directory
and performs RMSD calculations using the AfterMD toolkit.

Default behavior:
- Uses C-alpha atoms (GROMACS group 3) for both fitting and calculation
- Suitable for pHLA-TCR complex analysis

GROMACS Group Reference:
Group 0 (System): 384489 elements
Group 1 (Protein): 11067 elements
Group 2 (Protein-H): 5679 elements
Group 3 (C-alpha): 708 elements
Group 4 (Backbone): 2124 elements
Group 5 (MainChain): 2827 elements
"""

import os
import sys
import argparse
import logging
from pathlib import Path
from typing import List, Tuple, Optional
import glob

# Add parent directory to path for imports
sys.path.append(str(Path(__file__).parent.parent))

from aftermd.analysis.trajectory.rmsd import RMSDCalculator
from aftermd.utils.batch_processor import BatchProcessor

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def find_trajectory_pairs(input_dir: str) -> List[Tuple[str, str]]:
    """
    Find all trajectory-topology pairs in the input directory.

    Args:
        input_dir: Directory to search for trajectory files

    Returns:
        List of (trajectory_file, topology_file) pairs
    """
    input_path = Path(input_dir)
    trajectory_pairs = []

    # Find all .xtc files
    xtc_files = list(input_path.rglob("*.xtc"))
    logger.info(f"Found {len(xtc_files)} .xtc files")

    for xtc_file in xtc_files:
        # Look for corresponding topology files in the same directory
        xtc_dir = xtc_file.parent
        topology_file = None

        # Priority order: .tpr -> .gro -> .pdb
        for topo_ext in [".tpr", ".gro", ".pdb"]:
            potential_topo = xtc_dir / f"md{topo_ext}"
            if potential_topo.exists():
                topology_file = str(potential_topo)
                break

        # If md.* not found, look for any topology file in the directory
        if not topology_file:
            for topo_ext in [".tpr", ".gro", ".pdb"]:
                topo_files = list(xtc_dir.glob(f"*{topo_ext}"))
                if topo_files:
                    topology_file = str(topo_files[0])
                    break

        if topology_file:
            trajectory_pairs.append((str(xtc_file), topology_file))
            logger.info(f"Paired: {xtc_file.name} with {Path(topology_file).name}")
        else:
            logger.warning(f"No topology file found for {xtc_file}")

    return trajectory_pairs


def calculate_rmsd_single(trajectory: str, topology: str,
                         output_dir: str, rmsd_type: str = "calpha",
                         custom_fit_group: Optional[str] = None,
                         custom_calc_group: Optional[str] = None) -> Optional[str]:
    """
    Calculate RMSD for a single trajectory.

    Args:
        trajectory: Path to trajectory file
        topology: Path to topology file
        output_dir: Output directory for results
        rmsd_type: Type of RMSD calculation
        custom_fit_group: Custom group for fitting (overrides rmsd_type)
        custom_calc_group: Custom group for calculation (overrides rmsd_type)

    Returns:
        Output file path if successful, None otherwise
    """
    try:
        # Create output directory
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        # Generate unique output file name including parent directory
        traj_path = Path(trajectory)
        parent_name = traj_path.parent.name
        traj_name = traj_path.stem

        # Create unique filename: {parent_dir}_{trajectory_name}_rmsd.xvg
        if parent_name and parent_name != '.':
            output_filename = f"{parent_name}_{traj_name}_rmsd.xvg"
        else:
            output_filename = f"{traj_name}_rmsd.xvg"

        output_file = output_path / output_filename

        logger.info(f"Calculating RMSD for {traj_name}")

        # Initialize RMSD calculator
        rmsd_calc = RMSDCalculator(topology, trajectory)

        # Calculate RMSD using GROMACS method (more reliable for batch processing)
        if custom_fit_group and custom_calc_group:
            # Use custom groups (e.g., --fit-group 3 --calc-group 3 for C-alpha)
            logger.info(f"Using custom groups: fit={custom_fit_group}, calc={custom_calc_group}")
            result_file = rmsd_calc.calculate_gromacs_custom_groups(
                reference_group=custom_fit_group,
                analysis_group=custom_calc_group,
                output_file=str(output_file)
            )
        else:
            # Use predefined rmsd_type (default: calpha = group 3 for both fit and calc)
            logger.info(f"Using predefined RMSD type: {rmsd_type}")
            result_file = rmsd_calc.calculate_gromacs(
                rmsd_type=rmsd_type,
                output_file=str(output_file)
            )

        logger.info(f"RMSD calculation completed: {result_file}")
        return result_file

    except Exception as e:
        logger.error(f"RMSD calculation failed for {trajectory}: {e}")
        return None


def calculate_rmsd_batch(input_dir: str, output_dir: str,
                        rmsd_type: str = "calpha",
                        max_workers: int = 4,
                        custom_fit_group: Optional[str] = None,
                        custom_calc_group: Optional[str] = None) -> dict:
    """
    Perform batch RMSD calculations.

    Args:
        input_dir: Directory containing trajectory files
        output_dir: Directory for output files
        rmsd_type: Type of RMSD calculation
        max_workers: Number of parallel workers

    Returns:
        Dictionary with processing results
    """
    # Find trajectory pairs
    trajectory_pairs = find_trajectory_pairs(input_dir)

    if not trajectory_pairs:
        logger.error("No valid trajectory-topology pairs found")
        return {"success": 0, "failed": 0, "results": []}

    logger.info(f"Found {len(trajectory_pairs)} trajectory pairs for processing")

    # Process trajectories
    results = {"success": 0, "failed": 0, "results": []}

    for trajectory, topology in trajectory_pairs:
        result = calculate_rmsd_single(trajectory, topology, output_dir, rmsd_type,
                                     custom_fit_group, custom_calc_group)

        if result:
            results["success"] += 1
            results["results"].append({
                "trajectory": trajectory,
                "topology": topology,
                "output": result,
                "status": "success"
            })
        else:
            results["failed"] += 1
            results["results"].append({
                "trajectory": trajectory,
                "topology": topology,
                "output": None,
                "status": "failed"
            })

    return results


def main():
    """Main function."""
    parser = argparse.ArgumentParser(
        description="Batch RMSD calculation for MD trajectories",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Calculate C-alpha RMSD (default: group 3 vs group 3)
  python rmsd_calculator.py ./

  # Calculate RMSD with specific output directory
  python rmsd_calculator.py /path/to/trajectories -o /path/to/results

  # Calculate Backbone RMSD instead of C-alpha
  python rmsd_calculator.py ./data -t backbone

  # Use 8 parallel workers
  python rmsd_calculator.py ./data -j 8

  # Use custom GROMACS groups (group 3=C-alpha for both fit and calc)
  python rmsd_calculator.py ./data --fit-group 3 --calc-group 3

  # Use different groups (fit on Backbone, calc on C-alpha)
  python rmsd_calculator.py ./data --fit-group 4 --calc-group 3
        """
    )

    parser.add_argument("input_dir",
                       help="Directory containing trajectory files (.xtc)")

    parser.add_argument("-o", "--output",
                       default="rmsd_results",
                       help="Output directory for RMSD results (default: rmsd_results)")

    parser.add_argument("-t", "--type",
                       choices=["backbone", "protein", "calpha", "heavy"],
                       default="calpha",
                       help="Type of RMSD calculation (default: calpha)")

    parser.add_argument("-j", "--jobs",
                       type=int, default=4,
                       help="Number of parallel jobs (default: 4)")

    parser.add_argument("--fit-group",
                       type=str,
                       help="Custom GROMACS group for fitting (e.g., '4' for Backbone)")

    parser.add_argument("--calc-group",
                       type=str,
                       help="Custom GROMACS group for RMSD calculation (e.g., '3' for C-alpha)")

    parser.add_argument("-v", "--verbose",
                       action="store_true",
                       help="Enable verbose logging")

    args = parser.parse_args()

    # Set logging level
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    # Validate input directory
    if not os.path.exists(args.input_dir):
        logger.error(f"Input directory does not exist: {args.input_dir}")
        sys.exit(1)

    logger.info(f"Starting RMSD calculation")
    logger.info(f"Input directory: {args.input_dir}")
    logger.info(f"Output directory: {args.output}")
    logger.info(f"RMSD type: {args.type}")
    logger.info(f"Parallel jobs: {args.jobs}")

    # Validate custom groups
    if bool(args.fit_group) != bool(args.calc_group):
        logger.error("Both --fit-group and --calc-group must be specified together")
        sys.exit(1)

    # Perform batch RMSD calculations
    results = calculate_rmsd_batch(
        input_dir=args.input_dir,
        output_dir=args.output,
        rmsd_type=args.type,
        max_workers=args.jobs,
        custom_fit_group=args.fit_group,
        custom_calc_group=args.calc_group
    )

    # Print summary
    logger.info("=" * 50)
    logger.info("RMSD Calculation Summary")
    logger.info("=" * 50)
    logger.info(f"Total trajectories processed: {results['success'] + results['failed']}")
    logger.info(f"Successful calculations: {results['success']}")
    logger.info(f"Failed calculations: {results['failed']}")

    if results['failed'] > 0:
        logger.info("\nFailed trajectories:")
        for result in results['results']:
            if result['status'] == 'failed':
                logger.info(f"  - {Path(result['trajectory']).name}")

    if results['success'] > 0:
        logger.info(f"\nResults saved to: {args.output}")
        logger.info("Successful calculations:")
        for result in results['results']:
            if result['status'] == 'success':
                logger.info(f"  - {Path(result['output']).name}")

    logger.info("RMSD calculation completed!")


if __name__ == "__main__":
    main()