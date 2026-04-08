"""
PBC-RMSD Quality Pipeline

This module provides a comprehensive pipeline combining PBC correction and RMSD quality assessment.
It serves as a coordinator that orchestrates functional modules without reimplementing their logic.

Classes:
    PBCRMSDPipeline: Main pipeline for PBC processing + RMSD quality assessment

Author: Immunex Development Team
Date: 2026-03-15
"""

import os
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from concurrent.futures import ProcessPoolExecutor, as_completed
import json

from ..analysis.trajectory.pbc import PBCProcessor
from ..analysis.quality.post_pbc_validator import PostPBCValidator
from ..analysis.trajectory.rmsd import RMSDCalculator
from ..analysis.trajectory.rmsd_convergence import RMSDConvergenceAnalyzer
from ..utils.plotting import PlotManager


class PBCRMSDPipeline:
    """
    PBC Processing + RMSD Quality Assessment Pipeline

    This pipeline coordinates multiple functional modules to provide end-to-end
    processing from raw MD trajectory to quality-assessed results.

    Workflow:
    1. PBC correction (using PBCProcessor)
    2. Post-PBC validation (using PostPBCValidator)
    3. RMSD calculation (using RMSDCalculator)
    4. Convergence analysis (using RMSDConvergenceAnalyzer)
    5. Quality grading and reporting

    The pipeline acts as a coordinator and does NOT reimplement module functionality.
    All business logic is delegated to functional modules.

    Parameters
    ----------
    gmx_executable : str, default="gmx"
        GROMACS executable command
    max_com_drift : float, default=1.0
        Maximum allowed center-of-mass drift in nm (for PostPBCValidator)
    max_rg_std_ratio : float, default=0.15
        Maximum allowed Rg std/mean ratio (for PostPBCValidator)
    max_frame_jump_rmsd : float, default=0.5
        Maximum allowed frame-to-frame RMSD in nm (for PostPBCValidator)

    Examples
    --------
    >>> pipeline = PBCRMSDPipeline(gmx_executable="gmx")
    >>> results = pipeline.process_single_trajectory(
    ...     trajectory="md.xtc",
    ...     topology="md.tpr",
    ...     output_dir="./processed",
    ...     pbc_method="2step"
    ... )
    >>> print(f"Quality Grade: {results['overall_grade']}")
    """

    def __init__(self,
                 gmx_executable: str = "gmx",
                 max_com_drift: float = 1.0,
                 max_rg_std_ratio: float = 0.15,
                 max_frame_jump_rmsd: float = 0.5):
        """Initialize pipeline with functional module instances"""
        self.gmx_executable = gmx_executable
        self.logger = logging.getLogger(__name__)

        # Initialize functional modules
        self.pbc_processor = PBCProcessor(gmx_executable=gmx_executable)
        self.post_pbc_validator = PostPBCValidator(
            max_com_drift=max_com_drift,
            max_rg_std_ratio=max_rg_std_ratio,
            max_frame_jump_rmsd=max_frame_jump_rmsd
        )
        self.convergence_analyzer = RMSDConvergenceAnalyzer()
        self.plot_manager = PlotManager()

    def process_single_trajectory(self,
                                  trajectory: str,
                                  topology: str,
                                  output_dir: str,
                                  pbc_method: str = "2step",
                                  dt: Optional[float] = None,
                                  rmsd_selection: str = "backbone",
                                  run_quality_check: bool = True,
                                  validation_stride: int = 1,
                                  generate_report: bool = True) -> Dict:
        """
        Process single trajectory through complete pipeline

        This method coordinates all functional modules in sequence:
        1. Create output directory
        2. Run PBC correction (PBCProcessor)
        3. Optionally validate Post-PBC quality (PostPBCValidator)
        4. Calculate RMSD (RMSDCalculator)
        5. Analyze convergence (RMSDConvergenceAnalyzer)
        6. Assign overall quality grade
        7. Generate reports

        Parameters
        ----------
        trajectory : str
            Path to input trajectory file
        topology : str
            Path to topology file
        output_dir : str
            Output directory for processed files
        pbc_method : str, default="2step"
            PBC correction method: "2step" or "3step"
        dt : float, optional
            Time step for output trajectory in ps
        rmsd_selection : str, default="backbone"
            MDAnalysis selection string for RMSD calculation
        run_quality_check : bool, default=True
            Whether to run Post-PBC quality validation
        validation_stride : int, default=1
            Stride for validation (use >1 for faster checks on large trajectories)
        generate_report : bool, default=True
            Whether to generate quality report files

        Returns
        -------
        dict
            Processing results with keys:
            - pbc_output: path to PBC-corrected trajectory
            - post_pbc_validation: validation results (if run_quality_check=True)
            - rmsd_metrics: RMSD convergence metrics
            - overall_grade: quality grade ('A', 'B', 'C', 'D')
            - is_qualified: whether trajectory passes quality checks
            - report_files: paths to generated report files
        """
        self.logger.info(f"Starting pipeline for trajectory: {trajectory}")

        # Create output directory
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        results = {
            'pbc_output': None,
            'post_pbc_validation': None,
            'rmsd_metrics': None,
            'overall_grade': None,
            'is_qualified': False,
            'report_files': []
        }

        try:
            # Step 1: PBC correction
            self.logger.info("Step 1/5: Running PBC correction")
            pbc_output = output_path / "md_pbc.xtc"

            if pbc_method == "2step":
                self.pbc_processor.remove_pbc_2step(
                    trajectory=trajectory,
                    topology=topology,
                    output=str(pbc_output),
                    dt=dt
                )
            elif pbc_method == "3step":
                self.pbc_processor.remove_pbc_3step(
                    trajectory=trajectory,
                    topology=topology,
                    output=str(pbc_output),
                    dt=dt
                )
            else:
                raise ValueError(f"Unknown PBC method: {pbc_method}")

            results['pbc_output'] = str(pbc_output)
            self.logger.info(f"PBC-corrected trajectory: {pbc_output}")

            # Step 2: Post-PBC validation (optional)
            if run_quality_check:
                self.logger.info("Step 2/5: Running Post-PBC validation")
                validation_results = self.post_pbc_validator.comprehensive_validation(
                    trajectory=str(pbc_output),
                    topology=topology,
                    selection="protein",
                    stride=validation_stride
                )
                results['post_pbc_validation'] = validation_results

                if validation_results['overall_status'] == 'fail':
                    self.logger.warning(
                        f"Post-PBC validation failed: {validation_results['all_issues']}"
                    )
            else:
                self.logger.info("Step 2/5: Skipping Post-PBC validation")

            # Step 3: RMSD calculation
            self.logger.info("Step 3/5: Calculating RMSD")
            rmsd_calc = RMSDCalculator(topology, str(pbc_output))

            times, rmsd_values = rmsd_calc.calculate_mdanalysis(
                selection=rmsd_selection,
                output_file=str(output_path / "rmsd.csv")
            )

            # Step 4: Convergence analysis
            self.logger.info("Step 4/5: Analyzing RMSD convergence")
            rmsd_metrics = self.convergence_analyzer.calculate_convergence_metrics(
                times=times,
                rmsd_values=rmsd_values
            )
            results['rmsd_metrics'] = rmsd_metrics

            # Step 5: Overall quality grading
            self.logger.info("Step 5/5: Assigning quality grade")
            overall_grade = self._determine_overall_grade(
                post_pbc_validation=results.get('post_pbc_validation'),
                rmsd_metrics=rmsd_metrics
            )

            results['overall_grade'] = overall_grade
            results['is_qualified'] = overall_grade in ['A', 'B', 'C']

            # Generate reports
            if generate_report:
                report_files = self._generate_reports(
                    output_dir=output_path,
                    results=results,
                    trajectory_name=Path(trajectory).stem
                )
                results['report_files'] = report_files

            self.logger.info(
                f"Pipeline complete: Grade={overall_grade}, "
                f"Qualified={results['is_qualified']}"
            )

        except Exception as e:
            self.logger.error(f"Pipeline failed: {str(e)}")
            results['error'] = str(e)
            raise

        return results

    def _determine_overall_grade(self,
                                post_pbc_validation: Optional[Dict],
                                rmsd_metrics: Dict) -> str:
        """
        Determine overall quality grade from validation and RMSD metrics

        Grading logic:
        1. If Post-PBC validation grade is D, overall is D
        2. Otherwise, take the worse of Post-PBC grade and RMSD grade
        """
        rmsd_grade = self.convergence_analyzer.assign_quality_grade(rmsd_metrics)

        if post_pbc_validation is None:
            # No Post-PBC validation, use RMSD grade only
            return rmsd_grade

        pbc_grade = post_pbc_validation.get('overall_grade', 'D')

        # Grade D has highest priority
        if pbc_grade == 'D' or rmsd_grade == 'D':
            return 'D'

        # Take the worse grade
        grade_order = {'A': 1, 'B': 2, 'C': 3}
        pbc_score = grade_order.get(pbc_grade, 4)
        rmsd_score = grade_order.get(rmsd_grade, 4)

        worst_score = max(pbc_score, rmsd_score)
        grade_map = {1: 'A', 2: 'B', 3: 'C'}

        return grade_map.get(worst_score, 'D')

    def _generate_reports(self,
                         output_dir: Path,
                         results: Dict,
                         trajectory_name: str) -> List[str]:
        """
        Generate simplified quality report files

        Creates:
        1. Single TXT report (combines quality + convergence analysis)
        2. RMSD plot (PNG)
        """
        report_files = []

        # 1. Unified TXT report (combines quality + convergence)
        txt_file = output_dir / f"{trajectory_name}_quality_report.txt"
        txt_content = self._generate_unified_text_report(results, trajectory_name)
        with open(txt_file, 'w') as f:
            f.write(txt_content)
        report_files.append(str(txt_file))
        self.logger.info(f"Generated quality report: {txt_file}")

        # 2. RMSD plot
        rmsd_csv = output_dir / "rmsd.csv"
        if rmsd_csv.exists():
            rmsd_plot = output_dir / "rmsd.png"
            try:
                import pandas as pd
                import matplotlib.pyplot as plt

                # Read RMSD data
                df = pd.read_csv(rmsd_csv)

                # Create plot
                fig, ax = plt.subplots(figsize=(10, 6))
                ax.plot(df['Time (ps)'], df['RMSD (nm)'], linewidth=1.5, color='#2E86AB')

                # Add statistics
                if results.get('rmsd_metrics'):
                    metrics = results['rmsd_metrics']
                    mean_rmsd = metrics['mean_rmsd']
                    std_rmsd = metrics['std_rmsd']
                    grade = self.convergence_analyzer.assign_quality_grade(metrics)

                    # Add mean line
                    ax.axhline(y=mean_rmsd, color='red', linestyle='--',
                              linewidth=1, label=f'Mean: {mean_rmsd:.3f} nm')

                    # Add statistics text
                    stats_text = (
                        f'Grade: {grade}\n'
                        f'Mean: {mean_rmsd:.3f} ± {std_rmsd:.3f} nm\n'
                        f'CV: {metrics["cv_rmsd"]:.3f}\n'
                        f'Converged: {"Yes" if metrics["is_converged"] else "No"}'
                    )
                    ax.text(0.02, 0.98, stats_text, transform=ax.transAxes,
                           verticalalignment='top', bbox=dict(boxstyle='round',
                           facecolor='wheat', alpha=0.5))

                ax.set_xlabel('Time (ps)', fontsize=12)
                ax.set_ylabel('RMSD (nm)', fontsize=12)
                ax.set_title(f'RMSD Analysis - {trajectory_name}', fontsize=14, fontweight='bold')
                ax.grid(True, alpha=0.3)
                ax.legend(loc='upper right')

                plt.tight_layout()
                plt.savefig(rmsd_plot, dpi=300, bbox_inches='tight')
                plt.close()

                report_files.append(str(rmsd_plot))
                self.logger.info(f"Generated RMSD plot: {rmsd_plot}")

            except Exception as e:
                self.logger.warning(f"Failed to generate RMSD plot: {e}")

        return report_files

    def _generate_unified_text_report(self, results: Dict, trajectory_name: str) -> str:
        """Generate unified text report combining quality and convergence analysis"""
        lines = [
            "=" * 60,
            "MD Trajectory Quality Assessment Report",
            "=" * 60,
            "",
            f"Trajectory: {trajectory_name}",
            f"Overall Grade: {results['overall_grade']}",
            f"Quality Status: {'PASS' if results['is_qualified'] else 'FAIL'}",
            "",
            "=" * 60,
            "RMSD Quality Assessment",
            "=" * 60,
            ""
        ]

        # RMSD metrics
        if results.get('rmsd_metrics'):
            metrics = results['rmsd_metrics']
            rmsd_grade = self.convergence_analyzer.assign_quality_grade(metrics)

            lines.extend([
                f"Grade: {rmsd_grade}",
                "",
                "Overall Statistics:",
                f"  Mean RMSD: {metrics['mean_rmsd']:.3f} nm",
                f"  Std RMSD: {metrics['std_rmsd']:.3f} nm",
                f"  CV (Std/Mean): {metrics['cv_rmsd']:.3f}",
                f"  Min RMSD: {metrics['min_rmsd']:.3f} nm",
                f"  Max RMSD: {metrics['max_rmsd']:.3f} nm",
                "",
                "Convergence:",
                f"  Converged: {'Yes' if metrics['is_converged'] else 'No'}",
                f"  Convergence Time: {metrics['convergence_time']:.1f}",
                f"  Convergence Fraction: {metrics['convergence_fraction']:.1%}",
                ""
            ])

            # Block analysis
            if metrics.get('block_analysis'):
                block = metrics['block_analysis']
                lines.extend([
                    f"Block Analysis ({block['n_blocks']} blocks):",
                    f"  Block Stability (CV): {block['block_stability']:.3f}",
                    ""
                ])

                for info in block['block_info']:
                    t_range = f"{info['time_range'][0]:.0f}-{info['time_range'][1]:.0f}"
                    lines.append(
                        f"  Block {info['block_id']}: Mean={info['mean_rmsd']:.3f} nm, "
                        f"Std={info['std_rmsd']:.3f} nm"
                    )

                lines.append("")

            # Moving average
            if metrics.get('moving_avg'):
                ma = metrics['moving_avg']
                lines.extend([
                    "Moving Average Analysis:",
                    f"  Window Size: {ma['window_size']}",
                    f"  Trend Slope: {ma['trend_slope']:.6f} nm/time_unit",
                    f"  Has Drift: {'Yes' if ma['has_drift'] else 'No'}",
                    ""
                ])

        # Post-PBC validation
        if results.get('post_pbc_validation'):
            pbc_val = results['post_pbc_validation']
            lines.extend([
                "=" * 60,
                "Post-PBC Validation",
                "=" * 60,
                "",
                f"Status: {pbc_val['overall_status'].upper()}",
                f"Grade: {pbc_val['overall_grade']}",
                ""
            ])

            if pbc_val.get('pbc_quality'):
                pbc_q = pbc_val['pbc_quality']
                lines.extend([
                    "PBC Correction Quality:",
                    f"  COM Drift: {pbc_q.get('com_drift_nm', 0):.3f} nm",
                    ""
                ])

            if pbc_val.get('stability'):
                stab = pbc_val['stability']
                lines.extend([
                    "Structural Stability:",
                    f"  Rg Mean: {stab.get('rg_mean_nm', 0):.3f} nm",
                    f"  Rg Std: {stab.get('rg_std_nm', 0):.3f} nm",
                    f"  Rg Std/Mean: {stab.get('rg_std_ratio', 0):.3f}",
                    ""
                ])

            if pbc_val.get('all_issues'):
                lines.extend([
                    "Issues:",
                ])
                for issue in pbc_val['all_issues']:
                    lines.append(f"  - {issue}")
                lines.append("")

        lines.append("=" * 60)
        return "\n".join(lines)

    def batch_process(self,
                     tasks: List[Dict],
                     max_workers: int = 4,
                     fail_fast: bool = False) -> Dict:
        """
        Process multiple trajectories in parallel

        Parameters
        ----------
        tasks : List[Dict]
            List of task dictionaries, each containing:
            - trajectory: path to trajectory file
            - topology: path to topology file
            - output_dir: output directory
            - pbc_method: (optional) PBC method
            - dt: (optional) time step
            - rmsd_selection: (optional) RMSD selection
        max_workers : int, default=4
            Maximum number of parallel workers
        fail_fast : bool, default=False
            If True, stop on first failure; if False, continue with other tasks

        Returns
        -------
        dict
            Batch processing results with keys:
            - total_tasks: total number of tasks
            - successful: number of successful tasks
            - failed: number of failed tasks
            - task_results: list of individual task results
            - summary: grade distribution
        """
        self.logger.info(f"Starting batch processing of {len(tasks)} tasks")

        batch_results = {
            'total_tasks': len(tasks),
            'successful': 0,
            'failed': 0,
            'task_results': [],
            'summary': {'A': 0, 'B': 0, 'C': 0, 'D': 0}
        }

        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            # Submit all tasks
            future_to_task = {
                executor.submit(
                    self._process_task_wrapper,
                    task
                ): task for task in tasks
            }

            # Collect results
            for future in as_completed(future_to_task):
                task = future_to_task[future]
                task_name = Path(task['trajectory']).stem

                try:
                    result = future.result()
                    batch_results['successful'] += 1
                    batch_results['task_results'].append({
                        'task_name': task_name,
                        'status': 'success',
                        'result': result
                    })

                    # Update grade summary
                    grade = result.get('overall_grade', 'D')
                    batch_results['summary'][grade] += 1

                    self.logger.info(f"Task {task_name} completed: Grade {grade}")

                except Exception as e:
                    batch_results['failed'] += 1
                    batch_results['task_results'].append({
                        'task_name': task_name,
                        'status': 'failed',
                        'error': str(e)
                    })

                    self.logger.error(f"Task {task_name} failed: {str(e)}")

                    if fail_fast:
                        self.logger.error("Fail-fast enabled, stopping batch processing")
                        executor.shutdown(wait=False, cancel_futures=True)
                        break

        self.logger.info(
            f"Batch processing complete: {batch_results['successful']} successful, "
            f"{batch_results['failed']} failed"
        )

        return batch_results

    def _process_task_wrapper(self, task: Dict) -> Dict:
        """Wrapper for processing single task (used in batch processing)"""
        return self.process_single_trajectory(
            trajectory=task['trajectory'],
            topology=task['topology'],
            output_dir=task['output_dir'],
            pbc_method=task.get('pbc_method', '2step'),
            dt=task.get('dt'),
            rmsd_selection=task.get('rmsd_selection', 'backbone'),
            run_quality_check=task.get('run_quality_check', True),
            validation_stride=task.get('validation_stride', 1),
            generate_report=task.get('generate_report', True)
        )
