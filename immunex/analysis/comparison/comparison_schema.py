"""对比模块的数据结构定义。"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd


@dataclass(slots=True)
class SingleCaseArtifacts:
    """单体系对比输入面。"""

    label: str
    case_root: Path
    overview_root: Path
    identity: dict
    quality: dict
    bsa: dict
    rmsf_summary: dict
    rmsf_regions: pd.DataFrame = field(default_factory=pd.DataFrame)
    rrcs_summary: dict = field(default_factory=dict)
    rrcs_regions: pd.DataFrame = field(default_factory=pd.DataFrame)
    interaction_overview: pd.DataFrame = field(default_factory=pd.DataFrame)
    source_paths: dict = field(default_factory=dict)


@dataclass(slots=True)
class ComparisonArtifacts:
    """对比模块结构化产物。"""

    summary_json: Path
    comparison_table_csv: Path
    identity_comparison_csv: Path
    rmsf_region_comparison_csv: Path
    rrcs_region_comparison_csv: Path
    interaction_family_comparison_csv: Path
    quality_interface_plot: Path
    flexibility_plot: Path
    rrcs_plot: Path
    interaction_plot: Path
