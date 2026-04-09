#!/usr/bin/env python3
"""
Run complete CDR and TCR RMSD analysis with 20ns equilibration filtering
and generate combined boxplot
"""

import subprocess
import sys
import logging
import time
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def run_command(cmd, description, log_file=None):
    """Run a command and log output"""
    logger.info(f"Starting: {description}")
    logger.info(f"Command: {' '.join(cmd)}")

    if log_file:
        log_file = Path(log_file)
        log_file.parent.mkdir(parents=True, exist_ok=True)

        with open(log_file, 'w') as f:
            process = subprocess.Popen(
                cmd,
                stdout=f,
                stderr=subprocess.STDOUT,
                text=True
            )
            logger.info(f"Process started with PID: {process.pid}")
            logger.info(f"Log file: {log_file}")

            return process
    else:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True
        )

        if result.returncode != 0:
            logger.error(f"Failed: {description}")
            logger.error(f"Error: {result.stderr}")
            return None

        logger.info(f"Completed: {description}")
        return result


def main():
    python_bin = "/home/xumy/work/miniconda3/envs/immunex/bin/python"
    project_dir = Path("/home/xumy/work/development/Immunex")

    # Step 1: Run CDR RMSD analysis
    logger.info("=" * 80)
    logger.info("Step 1: Running CDR RMSD analysis with 20ns equilibration filter...")
    logger.info("=" * 80)

    cdr_script = project_dir / "scripts/batch_cdr_rmsd_exact.py"
    cdr_log = project_dir / "output/cdr_rmsd_exact_analysis/batch_run_equilibrated.log"

    cdr_process = run_command(
        [python_bin, str(cdr_script), "--max-workers", "8"],
        "CDR RMSD analysis",
        log_file=cdr_log
    )

    if cdr_process is None:
        logger.error("Failed to start CDR RMSD analysis")
        return 1

    # Step 2: Run TCR RMSD analysis
    logger.info("")
    logger.info("=" * 80)
    logger.info("Step 2: Running TCR RMSD analysis with 20ns equilibration filter...")
    logger.info("=" * 80)

    tcr_script = project_dir / "scripts/batch_tcr_rmsd.py"
    tcr_log = project_dir / "output/tcr_rmsd_phla_align/batch_run_equilibrated.log"

    tcr_process = run_command(
        [python_bin, str(tcr_script), "--max-workers", "8"],
        "TCR RMSD analysis",
        log_file=tcr_log
    )

    if tcr_process is None:
        logger.error("Failed to start TCR RMSD analysis")
        return 1

    # Wait for both analyses to complete
    logger.info("")
    logger.info("=" * 80)
    logger.info("Waiting for analyses to complete...")
    logger.info("=" * 80)
    logger.info(f"CDR RMSD PID: {cdr_process.pid}")
    logger.info(f"TCR RMSD PID: {tcr_process.pid}")

    logger.info("You can monitor progress with:")
    logger.info(f"  tail -f {cdr_log}")
    logger.info(f"  tail -f {tcr_log}")

    # Wait for CDR analysis
    logger.info("")
    logger.info("Waiting for CDR RMSD analysis...")
    cdr_returncode = cdr_process.wait()
    if cdr_returncode != 0:
        logger.error(f"CDR RMSD analysis failed with return code: {cdr_returncode}")
        logger.info("Check log file for details: {}".format(cdr_log))
    else:
        logger.info("CDR RMSD analysis completed successfully")

    # Wait for TCR analysis
    logger.info("")
    logger.info("Waiting for TCR RMSD analysis...")
    tcr_returncode = tcr_process.wait()
    if tcr_returncode != 0:
        logger.error(f"TCR RMSD analysis failed with return code: {tcr_returncode}")
        logger.info("Check log file for details: {}".format(tcr_log))
    else:
        logger.info("TCR RMSD analysis completed successfully")

    # Check if both analyses succeeded
    if cdr_returncode != 0 or tcr_returncode != 0:
        logger.error("One or both analyses failed. Skipping boxplot generation.")
        return 1

    # Step 3: Generate combined boxplot
    logger.info("")
    logger.info("=" * 80)
    logger.info("Step 3: Generating combined boxplot...")
    logger.info("=" * 80)

    plot_script = project_dir / "scripts/plot_cdr_tcr_combined_boxplot.py"

    plot_result = run_command(
        [python_bin, str(plot_script)],
        "Combined boxplot generation"
    )

    if plot_result is None:
        logger.error("Failed to generate combined boxplot")
        return 1

    logger.info("")
    logger.info("=" * 80)
    logger.info("All tasks completed successfully!")
    logger.info("=" * 80)
    logger.info("Output locations:")
    logger.info(f"  CDR RMSD: {project_dir}/output/cdr_rmsd_exact_analysis/")
    logger.info(f"  TCR RMSD: {project_dir}/output/tcr_rmsd_phla_align/")
    logger.info(f"  Combined plots: {project_dir}/output/cdr_tcr_combined_plots/")

    return 0


if __name__ == '__main__':
    exit(main())
