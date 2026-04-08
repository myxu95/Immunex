"""Normal mode / PRS 分析节点。"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Optional

from ...analysis.allostery.normal_mode import (
    NormalModeAnalyzer,
    write_hinge_profile,
    write_mode_mobility_profile,
    write_prs_ranking,
)
from ...core.base_node import PipelineNode
from ...core.context import PipelineContext
from ...core.exceptions import PipelineError


class NormalModeNode(PipelineNode):
    """执行 Cα 弹性网络的 normal mode / PRS 分析。"""

    def __init__(
        self,
        cutoff_angstrom: float = 10.0,
        n_low_modes: int = 10,
        prs_force_directions: int = 8,
        name: Optional[str] = None,
    ):
        super().__init__(name=name or "NormalModeNode")
        self.cutoff_angstrom = cutoff_angstrom
        self.n_low_modes = n_low_modes
        self.prs_force_directions = prs_force_directions

    def validate_inputs(self, context: PipelineContext) -> None:
        missing = []
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
            analyzer = NormalModeAnalyzer(context.structure_pdb)
            result = analyzer.calculate(
                chain_mapping=context.metadata["chain_mapping"],
                cdr_detection=context.metadata["cdr_detection"],
                cutoff_angstrom=self.cutoff_angstrom,
                n_low_modes=self.n_low_modes,
                prs_force_directions=self.prs_force_directions,
            )

            analysis_dir = Path(context.get_analysis_path("allostery", "normal_mode_residue_scores.csv")).parent
            analysis_dir.mkdir(parents=True, exist_ok=True)
            residue_csv = analysis_dir / "normal_mode_residue_scores.csv"
            region_csv = analysis_dir / "normal_mode_region_summary.csv"
            summary_json = analysis_dir / "normal_mode_summary.json"
            mobility_plot = analysis_dir / "mode_mobility_profile.png"
            hinge_plot = analysis_dir / "hinge_profile.png"
            prs_plot = analysis_dir / "prs_ranking.png"

            result.residue_frame.to_csv(residue_csv, index=False)
            result.region_summary.to_csv(region_csv, index=False)
            summary_json.write_text(
                json.dumps(result.summary, ensure_ascii=False, indent=2),
                encoding="utf-8",
            )
            write_mode_mobility_profile(result.residue_frame, mobility_plot)
            write_hinge_profile(result.residue_frame, hinge_plot)
            write_prs_ranking(result.residue_frame, prs_plot)

            context.results["normal_mode"] = {
                "residue_csv": str(residue_csv),
                "region_summary_csv": str(region_csv),
                "summary_json": str(summary_json),
                "mobility_plot": str(mobility_plot),
                "hinge_plot": str(hinge_plot),
                "prs_plot": str(prs_plot),
                "cutoff_angstrom": self.cutoff_angstrom,
                "n_low_modes": self.n_low_modes,
                "prs_force_directions": self.prs_force_directions,
            }
            return context
        except Exception as exc:
            raise PipelineError(
                node_name=self.name,
                reason=f"Normal mode analysis failed: {exc}",
                context_state={"system_id": context.system_id},
            ) from exc
