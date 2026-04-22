"""区域化 RMSF 分析节点。"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Optional

from ....analysis.trajectory import (
    ResidueRMSFAnalyzer,
    write_phla_rmsf_profile,
    write_region_rmsf_summary,
    write_tcr_rmsf_profile,
)
from ....core.base_node import PipelineNode
from ....core.context import PipelineContext
from ....core.exceptions import PipelineError


class ResidueRMSFNode(PipelineNode):
    """计算并落盘带语义注释的 residue-level RMSF。"""

    def __init__(
        self,
        selection: Optional[str] = None,
        stride: int = 1,
        time_unit: str = "ps",
        name: Optional[str] = None,
    ):
        super().__init__(name=name or "ResidueRMSFNode")
        self.selection = selection
        self.stride = stride
        self.time_unit = time_unit

    def validate_inputs(self, context: PipelineContext) -> None:
        missing = []
        if not context.topology:
            missing.append("topology")
        if not (context.trajectory_processed or context.trajectory_raw):
            missing.append("trajectory")
        if not context.structure_pdb:
            missing.append("structure_pdb")
        if "chain_mapping" not in context.metadata:
            missing.append("chain_mapping")
        if "cdr_detection" not in context.metadata:
            missing.append("cdr_detection")
        if missing:
            raise PipelineError(
                node_name=self.name,
                reason=f"Missing required inputs: {missing}",
                context_state={"system_id": context.system_id},
            )

    def execute(self, context: PipelineContext) -> PipelineContext:
        self.validate_inputs(context)

        try:
            analyzer = ResidueRMSFAnalyzer(
                topology_file=context.topology,
                trajectory_file=context.trajectory_processed or context.trajectory_raw,
                structure_pdb=context.structure_pdb,
            )
            result = analyzer.calculate(
                chain_mapping=context.metadata["chain_mapping"],
                cdr_detection=context.metadata["cdr_detection"],
                stride=self.stride,
                selection=self.selection,
                time_unit=self.time_unit,
            )

            analysis_dir = Path(context.get_analysis_path("rmsf", "residue_rmsf.csv")).parent
            analysis_dir.mkdir(parents=True, exist_ok=True)
            residue_csv = analysis_dir / "residue_rmsf.csv"
            summary_csv = analysis_dir / "region_rmsf_summary.csv"
            summary_json = analysis_dir / "rmsf_summary.json"
            tcr_plot = analysis_dir / "tcr_rmsf_profile.png"
            phla_plot = analysis_dir / "phla_rmsf_profile.png"
            region_plot = analysis_dir / "region_rmsf_summary.png"

            result.residue_frame.to_csv(residue_csv, index=False)
            result.region_summary.to_csv(summary_csv, index=False)
            summary_json.write_text(
                json.dumps(result.summary, indent=2, ensure_ascii=False),
                encoding="utf-8",
            )
            write_tcr_rmsf_profile(result.residue_frame, tcr_plot)
            write_phla_rmsf_profile(result.residue_frame, phla_plot)
            write_region_rmsf_summary(result.region_summary, region_plot)

            context.results["residue_rmsf"] = {
                "residue_csv": str(residue_csv),
                "region_summary_csv": str(summary_csv),
                "summary_json": str(summary_json),
                "tcr_plot": str(tcr_plot),
                "phla_plot": str(phla_plot),
                "region_plot": str(region_plot),
                "stride": self.stride,
                "time_unit": self.time_unit,
                "selection": result.summary.get("selection"),
            }
            return context
        except Exception as exc:
            raise PipelineError(
                node_name=self.name,
                reason=f"Residue RMSF analysis failed: {exc}",
                context_state={"system_id": context.system_id},
            ) from exc
