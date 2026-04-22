"""查询型 reporter 的问题分类与结果路由。"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import List


@dataclass
class QueryRoute:
    """描述一个查询问题的路由结果。"""

    query_type: str
    sources: List[Path]
    preferred_missing: List[str]


class QueryRouter:
    """根据问题文本将查询路由到最小必要结果源。"""

    def __init__(self, case_dir: str | Path):
        self.case_dir = Path(case_dir).resolve()

    def classify(self, question: str) -> str:
        text = question.strip().lower()
        if any(token in text for token in ["cluster", "dominant", "state", "状态"]):
            return "cluster_status"
        if any(token in text for token in ["diagnose stability", "stability", "stable", "quality", "rmsd", "轨迹", "稳定性"]):
            return "quality_status"
        if any(token in text for token in ["summarize interface", "interface summary", "interface", "界面总结", "界面概述"]):
            return "interface_summary"
        if any(token in text for token in ["suggest mutation", "mutation", "mutate", "candidate site", "design", "engineer", "突变", "设计"]):
            return "design_hint"
        if any(token in text for token in ["identify hotspots", "hotspots", "hotspot", "关键位点", "热点"]):
            return "hotspot_summary"
        if any(token in text for token in ["rrcs", "hotspot pair", "最强配对", "最强热点对", "最强热点"]):
            return "top_pair"
        if any(token in text for token in ["最稳定", "persistent", "occupancy", "持续", "稳定配对"]):
            return "persistent_pair"
        if any(token in text for token in ["区域", "region", "活跃区域", "dominant region"]):
            return "top_region"
        if any(token in text for token in ["残基", "residue", "位点", "关键"]):
            return "top_residue"
        return "top_residue"

    def route(self, question: str) -> QueryRoute:
        query_type = self.classify(question)
        if query_type == "cluster_status":
            return QueryRoute(
                query_type=query_type,
                sources=self._existing([
                    "overview/cluster/interface_clustering_summary.json",
                    "overview/cluster/summary_table.csv",
                    "overview/cluster/cluster_feature_digest.csv",
                ]),
                preferred_missing=["overview/cluster/interface_clustering_summary.json", "overview/cluster/summary_table.csv"],
            )
        if query_type == "quality_status":
            return QueryRoute(
                query_type=query_type,
                sources=self._existing(self._glob_candidates(["*quality*summary*.json", "*rmsd*summary*.json", "*tail*summary*.json"]))
                + self._existing([
                    "overview/rmsf/rmsf_summary.json",
                ]),
                preferred_missing=["quality summary or RMSD summary"],
            )
        if query_type == "interface_summary":
            return QueryRoute(
                query_type=query_type,
                sources=self._existing([
                    "overview/bsa/interface_summary.json",
                    "overview/interaction_overview.json",
                    "overview/rrcs/rrcs_region_summary.csv",
                    "overview/rrcs/annotated_rrcs_pair_summary.csv",
                ]),
                preferred_missing=["overview/bsa/interface_summary.json", "overview/interaction_overview.json"],
            )
        if query_type == "hotspot_summary":
            return QueryRoute(
                query_type=query_type,
                sources=self._existing([
                    "overview/rrcs/annotated_rrcs_pair_summary.csv",
                    "overview/rrcs/rrcs_summary.json",
                    "overview/rrcs/rrcs_region_summary.csv",
                    "overview/cluster/interface_clustering_summary.json",
                    "overview/cluster/summary_table.csv",
                ]),
                preferred_missing=["overview/rrcs/annotated_rrcs_pair_summary.csv", "overview/rrcs/rrcs_summary.json"],
            )
        if query_type == "design_hint":
            return QueryRoute(
                query_type=query_type,
                sources=self._existing([
                    "overview/rrcs/annotated_rrcs_pair_summary.csv",
                    "overview/rrcs/rrcs_summary.json",
                    "overview/rrcs/rrcs_region_summary.csv",
                    "overview/bsa/interface_summary.json",
                    "overview/interaction_overview.json",
                    "overview/cluster/interface_clustering_summary.json",
                    "overview/cluster/summary_table.csv",
                ]),
                preferred_missing=["overview/rrcs/annotated_rrcs_pair_summary.csv", "overview/rrcs/rrcs_summary.json"],
            )
        if query_type == "top_pair":
            return QueryRoute(
                query_type=query_type,
                sources=self._existing([
                    "overview/rrcs/annotated_rrcs_pair_summary.csv",
                    "overview/rrcs/rrcs_summary.json",
                    "overview/rrcs/rrcs_region_summary.csv",
                ]),
                preferred_missing=["overview/rrcs/annotated_rrcs_pair_summary.csv", "overview/rrcs/rrcs_summary.json"],
            )
        if query_type == "persistent_pair":
            return QueryRoute(
                query_type=query_type,
                sources=self._existing(self._glob_candidates(["*occupancy*summary*.json", "*persistence*summary*.json", "*occupancy*region*.csv", "*persistent*pair*.csv"]))
                + self._existing([
                    "overview/rrcs/annotated_rrcs_pair_summary.csv",
                    "overview/rrcs/rrcs_summary.json",
                ]),
                preferred_missing=["occupancy or persistence summary"],
            )
        if query_type == "top_region":
            return QueryRoute(
                query_type=query_type,
                sources=self._existing([
                    "overview/rrcs/rrcs_region_summary.csv",
                    "overview/interaction_overview.json",
                    "overview/cluster/cluster_feature_digest.csv",
                ]),
                preferred_missing=["overview/rrcs/rrcs_region_summary.csv", "overview/interaction_overview.json"],
            )
        return QueryRoute(
            query_type="top_residue",
            sources=self._existing([
                "overview/rrcs/annotated_rrcs_pair_summary.csv",
                "overview/rrcs/rrcs_summary.json",
                "overview/rrcs/rrcs_region_summary.csv",
                "overview/interaction_overview.json",
            ]),
            preferred_missing=["RRCS annotated pair summary or RRCS summary"],
        )

    def _existing(self, candidates: list[str | Path]) -> list[Path]:
        paths: list[Path] = []
        seen: set[Path] = set()
        for candidate in candidates:
            path = candidate if isinstance(candidate, Path) else self.case_dir / candidate
            if path.exists() and path not in seen:
                paths.append(path)
                seen.add(path)
        return paths

    def _glob_candidates(self, patterns: list[str]) -> list[Path]:
        matched: list[Path] = []
        for pattern in patterns:
            matched.extend(sorted(self.case_dir.rglob(pattern)))
        return matched
