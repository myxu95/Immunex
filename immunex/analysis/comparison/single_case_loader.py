"""单体系结果读取与标准化。"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from .comparison_schema import SingleCaseArtifacts


class SingleCaseLoader:
    """从单体系结果根目录读取 compare 所需的最小结果面。"""

    def load(self, case_root: Path, label: str) -> SingleCaseArtifacts:
        case_root = Path(case_root).resolve()
        overview_root = self._resolve_overview_root(case_root)

        identity_path = self._first_existing(
            overview_root / "identity" / "biological_identity.json",
            case_root / "analysis" / "identity" / "biological_identity.json",
        )
        quality_path = self._first_existing(
            overview_root / "quality" / "preprocess_quality_report.json",
            case_root / "quality" / "preprocess_quality_report.json",
        )
        bsa_path = self._first_existing(
            overview_root / "bsa" / "interface_summary.json",
            case_root / "analysis" / "interface" / "interface_summary.json",
        )
        rmsf_summary_path = self._first_existing(
            overview_root / "rmsf" / "rmsf_summary.json",
            case_root / "analysis" / "rmsf" / "rmsf_summary.json",
        )
        rmsf_region_path = self._first_existing(
            overview_root / "rmsf" / "region_rmsf_summary.csv",
            case_root / "analysis" / "rmsf" / "region_rmsf_summary.csv",
        )
        rrcs_summary_path = self._first_existing(
            overview_root / "rrcs" / "rrcs_summary.json",
            case_root / "analysis" / "interactions" / "rrcs" / "rrcs_summary.json",
        )
        rrcs_region_path = self._first_existing(
            overview_root / "rrcs" / "rrcs_region_summary.csv",
            case_root / "analysis" / "interactions" / "rrcs" / "rrcs_region_summary.csv",
        )
        interaction_overview_path = self._first_existing(
            overview_root / "interaction_overview.csv",
        )

        return SingleCaseArtifacts(
            label=label,
            case_root=case_root,
            overview_root=overview_root,
            identity=self._load_json(identity_path),
            quality=self._load_json(quality_path),
            bsa=self._load_json(bsa_path),
            rmsf_summary=self._load_json(rmsf_summary_path),
            rmsf_regions=self._load_csv(rmsf_region_path),
            rrcs_summary=self._load_json(rrcs_summary_path),
            rrcs_regions=self._load_csv(rrcs_region_path),
            interaction_overview=self._load_csv(interaction_overview_path),
            source_paths={
                "identity": str(identity_path) if identity_path else "",
                "quality": str(quality_path) if quality_path else "",
                "bsa": str(bsa_path) if bsa_path else "",
                "rmsf_summary": str(rmsf_summary_path) if rmsf_summary_path else "",
                "rmsf_regions": str(rmsf_region_path) if rmsf_region_path else "",
                "rrcs_summary": str(rrcs_summary_path) if rrcs_summary_path else "",
                "rrcs_regions": str(rrcs_region_path) if rrcs_region_path else "",
                "interaction_overview": str(interaction_overview_path) if interaction_overview_path else "",
            },
        )

    def _resolve_overview_root(self, case_root: Path) -> Path:
        overview_root = case_root / "overview"
        if overview_root.exists():
            return overview_root
        return case_root

    def _first_existing(self, *candidates: Path) -> Path | None:
        for candidate in candidates:
            if candidate.exists():
                return candidate
        return None

    def _load_json(self, path: Path | None) -> dict:
        if not path:
            return {}
        return json.loads(path.read_text(encoding="utf-8"))

    def _load_csv(self, path: Path | None) -> pd.DataFrame:
        if not path:
            return pd.DataFrame()
        return pd.read_csv(path)
