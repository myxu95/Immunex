"""
Standard Pipelines - Predefined pipeline templates.

This module provides commonly used pipeline configurations.
"""

from .base_pipeline import Pipeline
from .nodes import PreprocessNode, RMSDNode, RMSDPlotNode, RMSDQualityNode


class StandardTrajectoryPipeline(Pipeline):
    """
    Standard trajectory analysis pipeline.

    This pipeline performs:
    1. PBC preprocessing (3-step method)
    2. RMSD calculation (protein CA atoms)

    Example:
        >>> from immunex.pipeline import StandardTrajectoryPipeline
        >>> from immunex.core import PipelineContext
        >>>
        >>> context = PipelineContext(
        ...     system_id="1ao7",
        ...     topology="data/1ao7/md.tpr",
        ...     trajectory_raw="data/1ao7/md.xtc"
        ... )
        >>>
        >>> pipeline = StandardTrajectoryPipeline()
        >>> result = pipeline.execute(context)
        >>>
        >>> print(f"Mean RMSD: {result.results['rmsd']['mean_rmsd']:.3f} nm")
    """

    def __init__(self, preprocess_method: str = "3step", rmsd_selection: str = "protein and name CA"):
        """
        Initialize standard trajectory pipeline.

        Args:
            preprocess_method: PBC processing method ('2step' or '3step')
            rmsd_selection: Atom selection for RMSD calculation
        """
        nodes = [
            PreprocessNode(method=preprocess_method),
            RMSDNode(selection=rmsd_selection),
        ]
        super().__init__(nodes=nodes)


class QuickRMSDPipeline(Pipeline):
    """
    Quick RMSD pipeline (no preprocessing).

    Use this when trajectory is already preprocessed.

    Example:
        >>> pipeline = QuickRMSDPipeline()
        >>> context.trajectory_processed = "data/1ao7/md_processed.xtc"
        >>> result = pipeline.execute(context)
    """

    def __init__(self, rmsd_selection: str = "protein and name CA"):
        """
        Initialize quick RMSD pipeline.

        Args:
            rmsd_selection: Atom selection for RMSD calculation
        """
        nodes = [
            RMSDNode(selection=rmsd_selection),
        ]
        super().__init__(nodes=nodes)


class PreprocessOnlyPipeline(Pipeline):
    """
    Preprocessing-only pipeline.

    Use this to only perform PBC correction without analysis.

    Example:
        >>> pipeline = PreprocessOnlyPipeline(method="3step")
        >>> result = pipeline.execute(context)
        >>> print(result.trajectory_processed)
    """

    def __init__(self,
                 method: str = "3step",
                 dt: float = None,
                 gmx_executable: str = "gmx",
                 fit_group: str = None,
                 output_group: str = None,
                 center_group: str = None,
                 keep_temp_files: bool = False):
        """
        Initialize preprocess-only pipeline.

        Args:
            method: PBC processing method ('2step' or '3step')
            dt: Time interval for frame sampling in ps
            gmx_executable: GROMACS executable command
            fit_group: Optional fit group override
            output_group: Optional output group override
            center_group: Optional center group override
            keep_temp_files: Whether to keep temporary files
        """
        nodes = [
            PreprocessNode(
                method=method,
                dt=dt,
                gmx_executable=gmx_executable,
                fit_group=fit_group,
                output_group=output_group,
                center_group=center_group,
                keep_temp_files=keep_temp_files,
            ),
        ]
        super().__init__(nodes=nodes)


class PreprocessQualityPipeline(Pipeline):
    """
    预处理 + RMSD 质量评估流水线。

    该流水线以成熟的 pbc_rmsd 设计为参考，但保持在新架构中：
    1. PBC 预处理
    2. RMSD 计算
    3. RMSD 收敛分析与质量报告
    """

    def __init__(self,
                 method: str = "2step",
                 dt: float = None,
                 gmx_executable: str = "gmx",
                 fit_group: str = None,
                 output_group: str = None,
                 center_group: str = None,
                 keep_temp_files: bool = False,
                 rmsd_selection: str = "backbone"):
        nodes = [
            PreprocessNode(
                method=method,
                dt=dt,
                gmx_executable=gmx_executable,
                fit_group=fit_group,
                output_group=output_group,
                center_group=center_group,
                keep_temp_files=keep_temp_files,
            ),
            RMSDNode(
                selection=rmsd_selection,
                output_filename="rmsd.xvg",
                output_subdir="analysis/rmsd",
            ),
            RMSDPlotNode(),
            RMSDQualityNode(),
        ]
        super().__init__(nodes=nodes)
