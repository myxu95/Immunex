"""双体系结果对比分析模块。"""

from .comparison_schema import ComparisonArtifacts, SingleCaseArtifacts
from .single_case_loader import SingleCaseLoader
from .comparison_builder import (
    ComparisonBuildResult,
    SystemComparisonBuilder,
)
from .comparison_plotter import (
    write_flexibility_comparison_plot,
    write_interaction_family_comparison_plot,
    write_quality_interface_comparison_plot,
    write_rrcs_comparison_plot,
)

__all__ = [
    "ComparisonArtifacts",
    "SingleCaseArtifacts",
    "SingleCaseLoader",
    "ComparisonBuildResult",
    "SystemComparisonBuilder",
    "write_flexibility_comparison_plot",
    "write_interaction_family_comparison_plot",
    "write_quality_interface_comparison_plot",
    "write_rrcs_comparison_plot",
]
