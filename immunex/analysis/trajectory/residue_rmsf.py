"""
区域化 residue-level RMSF 分析模块。
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import MDAnalysis as mda
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from MDAnalysis.analysis.rms import RMSF

from ..topology.complex_residue_semantics import ComplexResidueSemanticAnnotator


TCR_REGION_ORDER = [
    "CDR1_alpha",
    "CDR2_alpha",
    "CDR3_alpha",
    "non_cdr_alpha",
    "CDR1_beta",
    "CDR2_beta",
    "CDR3_beta",
    "non_cdr_beta",
]

PHLA_REGION_ORDER = [
    "peptide",
    "alpha1_helix",
    "alpha2_helix",
    "non_groove",
    "beta2m",
]


@dataclass
class ResidueRMSFResult:
    """区域化 RMSF 结果容器。"""

    residue_frame: pd.DataFrame
    region_summary: pd.DataFrame
    summary: dict
    stride: int
    time_unit: str


class ResidueRMSFAnalyzer:
    """计算并注释 pHLA-TCR residue-level RMSF。"""

    def __init__(
        self,
        topology_file: str,
        trajectory_file: str,
        structure_pdb: str,
    ):
        self.topology_file = str(topology_file)
        self.trajectory_file = str(trajectory_file)
        self.structure_pdb = str(structure_pdb)
        self.trajectory_universe = mda.Universe(self.topology_file, self.trajectory_file)
        self.structure_universe = mda.Universe(self.structure_pdb)

    def calculate(
        self,
        chain_mapping: dict,
        cdr_detection: Optional[dict] = None,
        stride: int = 1,
        selection: Optional[str] = None,
        time_unit: str = "ps",
    ) -> ResidueRMSFResult:
        """计算 residue-level RMSF 并进行区域注释。"""
        if stride < 1:
            raise ValueError("stride 必须 >= 1")

        selected_chain_ids = [chain_id for chain_id in chain_mapping.values() if chain_id]
        if selection is None:
            selected_chain_clause = " or ".join(f"chainID {chain_id}" for chain_id in selected_chain_ids)
            selection = f"name CA and ({selected_chain_clause})"

        structure_atoms = self.structure_universe.select_atoms(selection)
        if len(structure_atoms) == 0:
            raise ValueError(f"未找到 RMSF 选择原子: {selection}")

        trajectory_atoms = self.trajectory_universe.atoms[structure_atoms.indices]
        if len(trajectory_atoms) != len(structure_atoms):
            raise ValueError("结构参考与轨迹拓扑原子索引不一致，无法映射 RMSF")

        rmsf_values = RMSF(trajectory_atoms).run(step=stride).results.rmsf
        if len(rmsf_values) != len(structure_atoms):
            raise RuntimeError("RMSF 结果长度与原子数量不一致")

        annotator = ComplexResidueSemanticAnnotator(
            chain_mapping=chain_mapping,
            cdr_detection=cdr_detection,
        )

        records = []
        for atom, rmsf in zip(structure_atoms, rmsf_values):
            semantics = annotator.annotate_residue(atom.chainID, int(atom.resid))
            region_group = self._resolve_region_group(semantics)
            records.append(
                {
                    "chain_id": atom.chainID,
                    "resid": int(atom.resid),
                    "resname": str(atom.resname),
                    "atom_name": str(atom.name),
                    "component": semantics.component,
                    "complex_side": semantics.complex_side,
                    "phla_region": semantics.phla_region,
                    "mhc_region": semantics.mhc_region,
                    "mhc_subregion": semantics.mhc_subregion,
                    "tcr_chain": semantics.tcr_chain,
                    "tcr_region": semantics.tcr_region,
                    "tcr_region_detailed": semantics.tcr_region_detailed,
                    "region_group": region_group,
                    "rmsf_angstrom": float(rmsf),
                }
            )

        residue_frame = pd.DataFrame.from_records(records)
        region_summary = self._summarize_regions(residue_frame)

        n_frames = len(self.trajectory_universe.trajectory[::stride])
        start_time = float(self.trajectory_universe.trajectory[0].time)
        end_time = float(self.trajectory_universe.trajectory[-1].time)
        time_divisor = 1000.0 if time_unit == "ns" else 1.0
        summary = {
            "selection": selection,
            "stride": int(stride),
            "n_residues": int(len(residue_frame)),
            "n_frames": int(n_frames),
            "time_unit": time_unit,
            "time_start": start_time / time_divisor,
            "time_end": end_time / time_divisor,
            "mean_rmsf_angstrom": float(residue_frame["rmsf_angstrom"].mean()),
            "max_rmsf_angstrom": float(residue_frame["rmsf_angstrom"].max()),
            "tcr_mean_rmsf_angstrom": float(
                residue_frame.loc[residue_frame["component"].isin(["TCR_alpha", "TCR_beta"]), "rmsf_angstrom"].mean()
            ),
            "phla_mean_rmsf_angstrom": float(
                residue_frame.loc[residue_frame["component"].isin(["HLA_alpha", "peptide"]), "rmsf_angstrom"].mean()
            ),
        }

        return ResidueRMSFResult(
            residue_frame=residue_frame,
            region_summary=region_summary,
            summary=summary,
            stride=stride,
            time_unit=time_unit,
        )

    @staticmethod
    def _resolve_region_group(semantics) -> str:
        if semantics.component in {"TCR_alpha", "TCR_beta"}:
            if semantics.tcr_region_detailed:
                return semantics.tcr_region_detailed
            if semantics.tcr_chain:
                return f"non_cdr_{semantics.tcr_chain}"
            return "tcr_other"
        if semantics.component == "peptide":
            return "peptide"
        if semantics.component == "HLA_alpha":
            return semantics.mhc_subregion or "non_groove"
        if semantics.component == "beta2m":
            return "beta2m"
        return "unknown"

    @staticmethod
    def _summarize_regions(residue_frame: pd.DataFrame) -> pd.DataFrame:
        summary = (
            residue_frame.groupby(["region_group"], as_index=False)
            .agg(
                n_residues=("resid", "count"),
                mean_rmsf_angstrom=("rmsf_angstrom", "mean"),
                median_rmsf_angstrom=("rmsf_angstrom", "median"),
                max_rmsf_angstrom=("rmsf_angstrom", "max"),
                min_rmsf_angstrom=("rmsf_angstrom", "min"),
            )
        )
        order = {name: idx for idx, name in enumerate(TCR_REGION_ORDER + PHLA_REGION_ORDER)}
        summary["sort_order"] = summary["region_group"].map(lambda value: order.get(value, 999))
        summary = summary.sort_values(["sort_order", "region_group"]).drop(columns=["sort_order"]).reset_index(drop=True)
        return summary


def write_tcr_rmsf_profile(residue_frame: pd.DataFrame, output_file: Path) -> None:
    """绘制 TCR alpha / beta residue RMSF 曲线。"""
    fig, ax = plt.subplots(figsize=(9.8, 4.8))
    plotted = False
    colors = {"alpha": "#2f5d50", "beta": "#c37a3b"}
    for chain in ["alpha", "beta"]:
        chain_df = residue_frame[residue_frame["tcr_chain"] == chain].copy()
        if chain_df.empty:
            continue
        chain_df = chain_df.sort_values("resid")
        ax.plot(
            chain_df["resid"],
            chain_df["rmsf_angstrom"],
            linewidth=1.8,
            color=colors[chain],
            label=f"TCR {chain}",
        )
        plotted = True
    if not plotted:
        ax.text(0.5, 0.5, "No TCR RMSF profile", ha="center", va="center", fontsize=13)
        ax.axis("off")
    else:
        ax.set_title("TCR RMSF profile", fontsize=14, weight="bold")
        ax.set_xlabel("Residue ID")
        ax.set_ylabel("RMSF (Å)")
        ax.grid(alpha=0.18, linestyle="--")
        ax.legend(frameon=False)
        ax.set_axisbelow(True)
    fig.tight_layout()
    fig.savefig(output_file, dpi=220, facecolor="white")
    plt.close(fig)


def write_phla_rmsf_profile(residue_frame: pd.DataFrame, output_file: Path) -> None:
    """绘制 peptide 与 HLA alpha RMSF 曲线。"""
    fig, axes = plt.subplots(2, 1, figsize=(9.8, 6.4), sharex=False)
    mhc_df = residue_frame[residue_frame["component"] == "HLA_alpha"].copy().sort_values("resid")
    peptide_df = residue_frame[residue_frame["component"] == "peptide"].copy().sort_values("resid")

    if mhc_df.empty:
        axes[0].text(0.5, 0.5, "No HLA alpha RMSF profile", ha="center", va="center", fontsize=12)
        axes[0].axis("off")
    else:
        axes[0].plot(mhc_df["resid"], mhc_df["rmsf_angstrom"], color="#658c7c", linewidth=1.8)
        axes[0].set_title("MHC alpha RMSF profile", fontsize=13, weight="bold")
        axes[0].set_ylabel("RMSF (Å)")
        axes[0].grid(alpha=0.18, linestyle="--")
        axes[0].set_axisbelow(True)

    if peptide_df.empty:
        axes[1].text(0.5, 0.5, "No peptide RMSF profile", ha="center", va="center", fontsize=12)
        axes[1].axis("off")
    else:
        axes[1].bar(peptide_df["resid"].astype(str), peptide_df["rmsf_angstrom"], color="#e0a34d")
        axes[1].set_title("Peptide RMSF profile", fontsize=13, weight="bold")
        axes[1].set_ylabel("RMSF (Å)")
        axes[1].set_xlabel("Residue ID")
        axes[1].grid(axis="y", alpha=0.18, linestyle="--")
        axes[1].set_axisbelow(True)
    fig.tight_layout()
    fig.savefig(output_file, dpi=220, facecolor="white")
    plt.close(fig)


def write_region_rmsf_summary(region_summary: pd.DataFrame, output_file: Path) -> None:
    """绘制区域级 RMSF 概览柱状图。"""
    fig, ax = plt.subplots(figsize=(11.2, 5.2))
    if region_summary.empty:
        ax.text(0.5, 0.5, "No region RMSF summary", ha="center", va="center", fontsize=13)
        ax.axis("off")
    else:
        top = region_summary.copy()
        ax.bar(
            top["region_group"],
            top["mean_rmsf_angstrom"],
            color=["#2f5d50" if value in TCR_REGION_ORDER else "#c37a3b" for value in top["region_group"]],
        )
        ax.set_title("Region-level RMSF summary", fontsize=14, weight="bold")
        ax.set_ylabel("Mean RMSF (Å)")
        ax.grid(axis="y", alpha=0.18, linestyle="--")
        ax.set_axisbelow(True)
        ax.tick_params(axis="x", rotation=30)
    fig.tight_layout()
    fig.savefig(output_file, dpi=220, facecolor="white")
    plt.close(fig)
