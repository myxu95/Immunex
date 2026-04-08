"""
基于 Cα 弹性网络模型的 normal mode / PRS 分析。

第一版目标：
- 基于单个代表结构构建 ENM
- 输出 residue-level mobility / hinge / PRS 指标
- 结合现有 TCR / pHLA 语义注释，形成关键位点排名
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional

import matplotlib.pyplot as plt
import MDAnalysis as mda
import numpy as np
import pandas as pd
from scipy.linalg import eigh
from scipy.spatial.distance import cdist

from ..topology import ComplexResidueSemanticAnnotator


@dataclass
class NormalModeResult:
    """Normal mode 分析结果。"""

    residue_frame: pd.DataFrame
    region_summary: pd.DataFrame
    summary: Dict[str, object]
    mode_count: int
    cutoff_angstrom: float


class NormalModeAnalyzer:
    """基于 Cα 弹性网络模型的 normal mode 分析器。"""

    def __init__(self, structure_pdb: str):
        self.structure_pdb = structure_pdb
        self.universe = mda.Universe(structure_pdb)

    def calculate(
        self,
        chain_mapping: Dict[str, str],
        cdr_detection: Optional[Dict] = None,
        cutoff_angstrom: float = 10.0,
        n_low_modes: int = 10,
        prs_force_directions: int = 8,
        gamma: float = 1.0,
    ) -> NormalModeResult:
        """执行 normal mode / PRS 分析。"""
        ca_atoms = self._select_ca_atoms(chain_mapping)
        if len(ca_atoms) < 10:
            raise ValueError("Cα 节点数量过少，无法建立稳定的弹性网络模型")

        positions = ca_atoms.positions.astype(float)
        pairwise_distance = cdist(positions, positions)
        contact_mask = (pairwise_distance > 0.0) & (pairwise_distance <= cutoff_angstrom)

        kirchhoff = self._build_kirchhoff(contact_mask=contact_mask, gamma=gamma)
        hessian = self._build_hessian(
            positions=positions,
            pairwise_distance=pairwise_distance,
            contact_mask=contact_mask,
            gamma=gamma,
        )

        gnm_covariance, gnm_mobility, gnm_modes = self._solve_gnm(
            kirchhoff=kirchhoff,
            n_low_modes=n_low_modes,
        )
        anm_covariance, anm_mobility, anm_modes = self._solve_anm(
            hessian=hessian,
            n_low_modes=n_low_modes,
        )
        hinge_score = self._compute_hinge_score(
            anm_modes=anm_modes,
            n_residues=len(ca_atoms.residues),
        )
        prs_matrix, prs_effectiveness, prs_sensitivity = self._compute_prs(
            covariance=anm_covariance,
            n_residues=len(ca_atoms.residues),
            n_force_directions=prs_force_directions,
        )
        network_degree = contact_mask.sum(axis=1).astype(float)

        residue_frame = self._build_residue_frame(
            ca_atoms=ca_atoms,
            chain_mapping=chain_mapping,
            cdr_detection=cdr_detection,
            gnm_mobility=gnm_mobility,
            anm_mobility=anm_mobility,
            hinge_score=hinge_score,
            prs_effectiveness=prs_effectiveness,
            prs_sensitivity=prs_sensitivity,
            network_degree=network_degree,
        )
        region_summary = self._build_region_summary(residue_frame)
        summary = self._build_summary(
            residue_frame=residue_frame,
            region_summary=region_summary,
            cutoff_angstrom=cutoff_angstrom,
            n_low_modes=n_low_modes,
            prs_force_directions=prs_force_directions,
            prs_matrix=prs_matrix,
            contact_mask=contact_mask,
        )

        return NormalModeResult(
            residue_frame=residue_frame,
            region_summary=region_summary,
            summary=summary,
            mode_count=n_low_modes,
            cutoff_angstrom=cutoff_angstrom,
        )

    def _select_ca_atoms(self, chain_mapping: Dict[str, str]):
        chain_ids = [chain_id for chain_id in chain_mapping.values() if chain_id]
        if not chain_ids:
            raise ValueError("chain_mapping 中未提供有效链信息")

        chain_selection = " or ".join(f"chainID {chain_id}" for chain_id in chain_ids)
        atoms = self.universe.select_atoms(f"({chain_selection}) and name CA")
        if len(atoms) == 0:
            raise ValueError("在输入结构中未找到对应链的 Cα 原子")
        return atoms

    @staticmethod
    def _build_kirchhoff(contact_mask: np.ndarray, gamma: float) -> np.ndarray:
        adjacency = contact_mask.astype(float)
        kirchhoff = -gamma * adjacency
        np.fill_diagonal(kirchhoff, -kirchhoff.sum(axis=1))
        return kirchhoff

    @staticmethod
    def _build_hessian(
        positions: np.ndarray,
        pairwise_distance: np.ndarray,
        contact_mask: np.ndarray,
        gamma: float,
    ) -> np.ndarray:
        n_residues = positions.shape[0]
        hessian = np.zeros((3 * n_residues, 3 * n_residues), dtype=float)

        contact_indices = np.argwhere(np.triu(contact_mask, k=1))
        for i, j in contact_indices:
            delta = positions[j] - positions[i]
            distance = pairwise_distance[i, j]
            if distance <= 0.0:
                continue
            unit = delta / distance
            outer = gamma * np.outer(unit, unit)
            i_slice = slice(3 * i, 3 * i + 3)
            j_slice = slice(3 * j, 3 * j + 3)
            hessian[i_slice, i_slice] += outer
            hessian[j_slice, j_slice] += outer
            hessian[i_slice, j_slice] -= outer
            hessian[j_slice, i_slice] -= outer

        return hessian

    @staticmethod
    def _positive_modes(eigenvalues: np.ndarray, threshold: float = 1e-8) -> np.ndarray:
        return np.where(eigenvalues > threshold)[0]

    def _solve_gnm(self, kirchhoff: np.ndarray, n_low_modes: int):
        eigenvalues, eigenvectors = eigh(kirchhoff)
        positive = self._positive_modes(eigenvalues)
        if positive.size == 0:
            raise ValueError("GNM 特征分解未得到有效非零模式")

        inv_lambda = 1.0 / eigenvalues[positive]
        covariance = (eigenvectors[:, positive] * inv_lambda) @ eigenvectors[:, positive].T
        mobility = np.diag(covariance)

        selected = positive[: min(n_low_modes, positive.size)]
        modes = {
            "eigenvalues": eigenvalues[selected],
            "eigenvectors": eigenvectors[:, selected],
        }
        return covariance, mobility, modes

    def _solve_anm(self, hessian: np.ndarray, n_low_modes: int):
        eigenvalues, eigenvectors = eigh(hessian)
        positive = self._positive_modes(eigenvalues)
        if positive.size == 0:
            raise ValueError("ANM 特征分解未得到有效非零模式")

        inv_lambda = 1.0 / eigenvalues[positive]
        covariance = (eigenvectors[:, positive] * inv_lambda) @ eigenvectors[:, positive].T
        n_residues = hessian.shape[0] // 3
        mobility = np.array(
            [
                float(np.trace(covariance[3 * i : 3 * i + 3, 3 * i : 3 * i + 3]))
                for i in range(n_residues)
            ]
        )

        selected = positive[: min(n_low_modes, positive.size)]
        modes = {
            "eigenvalues": eigenvalues[selected],
            "eigenvectors": eigenvectors[:, selected],
        }
        return covariance, mobility, modes

    @staticmethod
    def _compute_hinge_score(anm_modes: Dict[str, np.ndarray], n_residues: int) -> np.ndarray:
        eigenvalues = anm_modes["eigenvalues"]
        eigenvectors = anm_modes["eigenvectors"]
        if eigenvectors.size == 0:
            return np.zeros(n_residues, dtype=float)

        amplitudes = np.zeros(n_residues, dtype=float)
        for mode_idx in range(eigenvectors.shape[1]):
            mode_vector = eigenvectors[:, mode_idx].reshape(n_residues, 3)
            residue_amplitude = np.linalg.norm(mode_vector, axis=1) ** 2
            amplitudes += residue_amplitude / max(float(eigenvalues[mode_idx]), 1e-8)

        return _normalize_01(amplitudes.max() - amplitudes)

    @staticmethod
    def _compute_prs(covariance: np.ndarray, n_residues: int, n_force_directions: int):
        rng = np.random.default_rng(20260406)
        force_directions = rng.normal(size=(n_force_directions, 3))
        force_directions /= np.linalg.norm(force_directions, axis=1, keepdims=True)

        response_matrix = np.zeros((n_residues, n_residues), dtype=float)
        for perturbed_residue in range(n_residues):
            block = covariance[:, 3 * perturbed_residue : 3 * perturbed_residue + 3]
            accumulated = np.zeros(n_residues, dtype=float)
            for direction in force_directions:
                displacement = block @ direction
                residue_displacement = displacement.reshape(n_residues, 3)
                accumulated += np.linalg.norm(residue_displacement, axis=1)
            response_matrix[perturbed_residue, :] = accumulated / n_force_directions

        effectiveness = response_matrix.mean(axis=1)
        sensitivity = response_matrix.mean(axis=0)
        return response_matrix, effectiveness, sensitivity

    def _build_residue_frame(
        self,
        ca_atoms,
        chain_mapping: Dict[str, str],
        cdr_detection: Optional[Dict],
        gnm_mobility: np.ndarray,
        anm_mobility: np.ndarray,
        hinge_score: np.ndarray,
        prs_effectiveness: np.ndarray,
        prs_sensitivity: np.ndarray,
        network_degree: np.ndarray,
    ) -> pd.DataFrame:
        annotator = ComplexResidueSemanticAnnotator(
            chain_mapping=chain_mapping,
            cdr_detection=cdr_detection,
        )

        records = []
        gnm_norm = _normalize_01(gnm_mobility)
        anm_norm = _normalize_01(anm_mobility)
        prs_eff_norm = _normalize_01(prs_effectiveness)
        prs_sens_norm = _normalize_01(prs_sensitivity)
        network_degree_norm = _normalize_01(network_degree)

        for idx, atom in enumerate(ca_atoms):
            semantics = annotator.annotate_residue(chain_id=atom.chainID, resid=int(atom.resid))
            region_group = _derive_region_group(semantics)
            combined_score = (
                0.32 * hinge_score[idx]
                + 0.30 * prs_eff_norm[idx]
                + 0.13 * prs_sens_norm[idx]
                + 0.25 * network_degree_norm[idx]
            )
            records.append(
                {
                    "chain_id": atom.chainID,
                    "resid": int(atom.resid),
                    "resname": atom.resname,
                    "residue_label": f"{atom.resname}{int(atom.resid)}",
                    "component": semantics.component,
                    "complex_side": semantics.complex_side,
                    "phla_region": semantics.phla_region,
                    "mhc_region": semantics.mhc_region,
                    "mhc_subregion": semantics.mhc_subregion,
                    "tcr_chain": semantics.tcr_chain,
                    "tcr_region": semantics.tcr_region,
                    "tcr_region_detailed": semantics.tcr_region_detailed,
                    "region_group": region_group,
                    "gnm_mobility": float(gnm_mobility[idx]),
                    "gnm_mobility_norm": float(gnm_norm[idx]),
                    "anm_mobility": float(anm_mobility[idx]),
                    "anm_mobility_norm": float(anm_norm[idx]),
                    "hinge_score": float(hinge_score[idx]),
                    "prs_effectiveness": float(prs_effectiveness[idx]),
                    "prs_effectiveness_norm": float(prs_eff_norm[idx]),
                    "prs_sensitivity": float(prs_sensitivity[idx]),
                    "prs_sensitivity_norm": float(prs_sens_norm[idx]),
                    "network_degree": float(network_degree[idx]),
                    "network_degree_norm": float(network_degree_norm[idx]),
                    "combined_key_residue_score": float(combined_score),
                }
            )

        frame = pd.DataFrame(records).sort_values(
            ["combined_key_residue_score", "prs_effectiveness", "hinge_score"],
            ascending=[False, False, False],
        )
        return frame.reset_index(drop=True)

    @staticmethod
    def _build_region_summary(residue_frame: pd.DataFrame) -> pd.DataFrame:
        summary = (
            residue_frame.groupby("region_group", dropna=False)
            .agg(
                n_residues=("residue_label", "count"),
                mean_gnm_mobility=("gnm_mobility", "mean"),
                mean_anm_mobility=("anm_mobility", "mean"),
                mean_hinge_score=("hinge_score", "mean"),
                mean_prs_effectiveness=("prs_effectiveness", "mean"),
                max_combined_score=("combined_key_residue_score", "max"),
            )
            .reset_index()
            .sort_values("max_combined_score", ascending=False)
        )
        return summary

    @staticmethod
    def _build_summary(
        residue_frame: pd.DataFrame,
        region_summary: pd.DataFrame,
        cutoff_angstrom: float,
        n_low_modes: int,
        prs_force_directions: int,
        prs_matrix: np.ndarray,
        contact_mask: np.ndarray,
    ) -> Dict[str, object]:
        top_key_residues = (
            residue_frame.nlargest(10, "combined_key_residue_score")[
                [
                    "chain_id",
                    "resid",
                    "resname",
                    "region_group",
                    "combined_key_residue_score",
                    "hinge_score",
                    "prs_effectiveness",
                ]
            ]
            .to_dict(orient="records")
        )
        top_hinges = (
            residue_frame.nlargest(10, "hinge_score")[
                ["chain_id", "resid", "resname", "region_group", "hinge_score"]
            ].to_dict(orient="records")
        )
        summary = {
            "n_residues": int(len(residue_frame)),
            "network_edges": int(np.count_nonzero(np.triu(contact_mask, k=1))),
            "cutoff_angstrom": float(cutoff_angstrom),
            "n_low_modes": int(n_low_modes),
            "prs_force_directions": int(prs_force_directions),
            "mean_gnm_mobility": float(residue_frame["gnm_mobility"].mean()),
            "mean_anm_mobility": float(residue_frame["anm_mobility"].mean()),
            "mean_hinge_score": float(residue_frame["hinge_score"].mean()),
            "mean_prs_effectiveness": float(residue_frame["prs_effectiveness"].mean()),
            "mean_prs_sensitivity": float(residue_frame["prs_sensitivity"].mean()),
            "mean_network_degree": float(residue_frame["network_degree"].mean()),
            "prs_matrix_shape": [int(prs_matrix.shape[0]), int(prs_matrix.shape[1])],
            "top_key_residues": top_key_residues,
            "top_hinges": top_hinges,
            "top_regions": region_summary.head(10).to_dict(orient="records"),
        }
        return summary


def write_mode_mobility_profile(residue_frame: pd.DataFrame, output_file: Path) -> None:
    """绘制 GNM/ANM mobility residue profile。"""
    ordered = residue_frame.sort_values(["chain_id", "resid"]).reset_index(drop=True)
    x = np.arange(1, len(ordered) + 1)

    fig, ax = plt.subplots(figsize=(11.2, 4.8))
    ax.plot(x, ordered["gnm_mobility_norm"], color="#2f5d50", linewidth=1.8, label="GNM mobility")
    ax.plot(x, ordered["anm_mobility_norm"], color="#c37a3b", linewidth=1.8, label="ANM mobility")
    ax.set_xlabel("Residue index")
    ax.set_ylabel("Normalized mobility")
    ax.set_title("Low-frequency mode mobility profile", fontsize=13, weight="bold")
    ax.grid(alpha=0.18, linestyle="--")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(output_file, dpi=220, facecolor="white")
    plt.close(fig)


def write_hinge_profile(residue_frame: pd.DataFrame, output_file: Path) -> None:
    """绘制 hinge score profile。"""
    ordered = residue_frame.sort_values(["chain_id", "resid"]).reset_index(drop=True)
    x = np.arange(1, len(ordered) + 1)
    fig, ax = plt.subplots(figsize=(11.2, 4.8))
    ax.fill_between(x, 0.0, ordered["hinge_score"], color="#86a897", alpha=0.55)
    ax.plot(x, ordered["hinge_score"], color="#2f5d50", linewidth=1.8)
    ax.set_xlabel("Residue index")
    ax.set_ylabel("Hinge score")
    ax.set_title("Hinge-like residue profile", fontsize=13, weight="bold")
    ax.grid(alpha=0.18, linestyle="--")
    fig.tight_layout()
    fig.savefig(output_file, dpi=220, facecolor="white")
    plt.close(fig)


def write_prs_ranking(residue_frame: pd.DataFrame, output_file: Path, top_n: int = 15) -> None:
    """绘制关键位点排名图。"""
    top_frame = residue_frame.nlargest(top_n, "combined_key_residue_score").copy()
    if top_frame.empty:
        top_frame = residue_frame.head(min(top_n, len(residue_frame))).copy()
    labels = [
        f"{row.resname}{int(row.resid)} ({row.region_group})"
        for row in top_frame.itertuples(index=False)
    ]
    y = np.arange(len(top_frame))
    fig, ax = plt.subplots(figsize=(11.2, max(5.2, len(top_frame) * 0.45)))
    ax.barh(y, top_frame["combined_key_residue_score"], color="#2f5d50", alpha=0.9)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel("Combined key residue score")
    ax.set_title("Top key residues by normal mode / PRS score", fontsize=13, weight="bold")
    ax.grid(axis="x", alpha=0.18, linestyle="--")
    fig.tight_layout()
    fig.savefig(output_file, dpi=220, facecolor="white")
    plt.close(fig)


def _normalize_01(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    if values.size == 0:
        return values
    value_min = float(np.min(values))
    value_max = float(np.max(values))
    if abs(value_max - value_min) < 1e-12:
        return np.zeros_like(values, dtype=float)
    return (values - value_min) / (value_max - value_min)


def _derive_region_group(semantics) -> str:
    if semantics.tcr_region_detailed:
        return semantics.tcr_region_detailed
    if semantics.mhc_subregion:
        return semantics.mhc_subregion
    if semantics.component != "unknown":
        return semantics.component
    return "unknown"
