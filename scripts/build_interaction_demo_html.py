#!/usr/bin/env python3
"""生成单体系 interaction 静态 HTML 报告。"""

from __future__ import annotations

import argparse
import csv
import html
import json
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from immunex.analysis.topology import BiologicalIdentityAnnotator


VIEWER_PDB_NAME = "phla_tcr_complex_viewer.pdb"
VIEWER_META_NAME = "phla_tcr_viewer_meta.json"
FIGURES_DIR_NAME = "figures"
QUALITY_DIR_NAME = "quality"
BSA_DIR_NAME = "bsa"
RMSF_DIR_NAME = "rmsf"
IDENTITY_DIR_NAME = "identity"
CLUSTER_DIR_NAME = "cluster"
ALL_REPORT_SECTIONS = (
    "overview",
    "quality",
    "interface",
    "flexibility",
    "cluster",
    "occupancy",
    "contact",
    "interactions",
    "downloads",
)
FAMILY_COLORS = {
    "contact": "#2d6f5c",
    "hbond": "#3f8a77",
    "saltbridge": "#b56b3a",
    "hydrophobic": "#d7922f",
    "pipi": "#8661c1",
    "cationpi": "#c05a84",
}
REPORT_PATHS: dict[str, Path] = {}
OCCUPANCY_LAYOUTS = {
    "contact": ("analysis/contacts/occupancy", "contact"),
    "hbond": ("analysis/interactions/hydrogen_bonds/occupancy", "hbond"),
    "saltbridge": ("analysis/interactions/salt_bridges/occupancy", "saltbridge"),
    "hydrophobic": ("analysis/interactions/hydrophobic_contacts/occupancy", "hydrophobic"),
    "pipi": ("analysis/interactions/pi_interactions/occupancy", "pipi"),
    "cationpi": ("analysis/interactions/cation_pi_interactions/occupancy", "cationpi"),
}
HEATMAP_CMAP = LinearSegmentedColormap.from_list(
    "immunex_report_occupancy",
    ["#fffaf2", "#eef3ed", "#c8d8cf", "#7da48f", "#2f5d50"],
)


def build_report_paths(base_dir: Path, system_id: str) -> dict[str, Path]:
    return {
        "contact": base_dir / "contact" / system_id / "analysis/contacts/contact_report.csv",
        "hbond": base_dir / "hbond" / system_id / "analysis/interactions/hydrogen_bonds/hbond_report.csv",
        "saltbridge": base_dir / "saltbridge" / system_id / "analysis/interactions/salt_bridges/salt_bridge_report.csv",
        "hydrophobic": base_dir / "hydrophobic" / system_id / "analysis/interactions/hydrophobic_contacts/hydrophobic_report.csv",
        "pipi": base_dir / "pipi" / system_id / "analysis/interactions/pi_interactions/pi_pi_report.csv",
        "cationpi": base_dir / "cationpi" / system_id / "analysis/interactions/cation_pi_interactions/cation_pi_report.csv",
    }


def ensure_symlink(dst: Path, src: Path) -> None:
    if dst.exists() or dst.is_symlink():
        dst.unlink()
    dst.symlink_to(src.resolve())


def choose_source_pdb(base_dir: Path, system_id: str, cli_value: Path | None) -> Path:
    candidates = []
    if cli_value is not None:
        candidates.append(cli_value)
    candidates.extend(
        [
            Path(f"/home/xumy/work/development/Immunex/output/parallel2_preprocess/{system_id}/md_processed_converted.pdb"),
            base_dir / "contact" / system_id / "structure.pdb",
            base_dir / "contact" / system_id / "md_converted.pdb",
            base_dir / "contact" / system_id / "md_processed_converted.pdb",
            Path(f"/home/xumy/work/development/Immunex/output/.interaction_batch_staging/single_raw/{system_id}/md_converted.pdb"),
            Path(f"/home/xumy/work/development/Immunex/output/.interaction_case_staging/{system_id}_raw/{system_id}/md_processed_converted.pdb"),
        ]
    )
    for candidate in candidates:
        if candidate.exists():
            return candidate
    raise FileNotFoundError(f"未找到 viewer 结构文件，可候选路径: {[str(p) for p in candidates]}")


def choose_quality_root(cli_value: Path | None, source_pdb: Path) -> Path | None:
    candidates: list[Path] = []
    if cli_value is not None:
        candidates.append(cli_value)
    candidates.append(source_pdb.parent)
    for candidate in candidates:
        if not candidate:
            continue
        if (
            (candidate / "quality/preprocess_quality_report.json").exists()
            and (candidate / "analysis/rmsd/rmsd.xvg").exists()
        ):
            return candidate
    return None


def choose_bsa_root(cli_value: Path | None, source_pdb: Path, system_id: str) -> Path | None:
    candidates: list[Path] = []
    if cli_value is not None:
        candidates.append(cli_value)
    candidates.extend(
        [
            source_pdb.parent,
            Path(f"/home/xumy/work/development/Immunex/output/bsa_demo_{system_id}"),
            Path(f"/home/xumy/work/development/Immunex/output/bsa_demo_{system_id.split('_')[0]}"),
            Path(f"/home/xumy/work/development/Immunex/output/bsa_batch/{source_pdb.parent.parent.name}/{system_id}"),
            Path(f"/home/xumy/work/development/Immunex/output/bsa_batch_smoke/{system_id}"),
        ]
    )
    for candidate in candidates:
        if not candidate:
            continue
        if (
            (candidate / "analysis/interface/interface_summary.json").exists()
            and (candidate / "analysis/interface/bsa_timeseries.csv").exists()
        ):
            return candidate
    return None


def choose_rmsf_root(cli_value: Path | None, source_pdb: Path, system_id: str) -> Path | None:
    candidates: list[Path] = []
    if cli_value is not None:
        candidates.append(cli_value)
    candidates.extend(
        [
            source_pdb.parent,
            Path(f"/home/xumy/work/development/Immunex/output/rmsf_demo_{system_id}"),
            Path(f"/home/xumy/work/development/Immunex/output/rmsf_demo_{system_id.split('_')[0]}"),
            Path(f"/home/xumy/work/development/Immunex/output/rmsf_batch_smoke/{system_id}"),
        ]
    )
    for candidate in candidates:
        if not candidate:
            continue
        if (
            (candidate / "analysis/rmsf/residue_rmsf.csv").exists()
            and (candidate / "analysis/rmsf/region_rmsf_summary.csv").exists()
            and (candidate / "analysis/rmsf/rmsf_summary.json").exists()
        ):
            return candidate
    return None


def choose_cluster_root(cli_value: Path | None, system_id: str) -> Path | None:
    candidates: list[Path] = []
    if cli_value is not None:
        candidates.append(cli_value)
    candidates.extend(
        [
            Path(f"/home/xumy/work/development/Immunex/output/interface_cluster_demo_{system_id}"),
            Path(f"/home/xumy/work/development/Immunex/output/interface_cluster_demo_{system_id.split('_')[0]}"),
            Path(f"/home/xumy/work/development/Immunex/output/interface_cluster_demo_{system_id.split('_')[0]}_s5"),
            Path(f"/home/xumy/work/development/Immunex/output/interface_cluster_demo_{system_id.split('_')[0]}_s5_v2"),
            Path(f"/home/xumy/work/development/Immunex/output/interface_cluster_demo_{system_id.split('_')[0]}_s5_v3"),
        ]
    )
    for candidate in candidates:
        if not candidate:
            continue
        analysis_dir = candidate / "analysis/conformation/interface_clustering"
        if (
            (analysis_dir / "summary_table.csv").exists()
            and (analysis_dir / "cluster_feature_digest.csv").exists()
            and (analysis_dir / "cluster_id_vs_time.png").exists()
            and (analysis_dir / "state_population_over_time.png").exists()
        ):
            return candidate
    return None


def load_overview_rows(csv_path: Path) -> list[dict[str, str]]:
    with csv_path.open("r", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def build_cards(rows: list[dict[str, str]]) -> str:
    cards = []
    for row in rows:
        family = html.escape(row["family"])
        label = html.escape(row["label"])
        pair_rows = html.escape(row["pair_rows"])
        report_rows = html.escape(row["report_rows"])
        report_href = html.escape(row["overview_report"])
        heatmap_count = html.escape(row["heatmap_count"])
        status = "empty" if row["pair_rows"] == "0" else "active"
        cards.append(
            f"""
            <article class="family-card {status}">
              <div class="family-head">
                <div>
                  <div class="family-name">{family}</div>
                  <div class="family-label">{label}</div>
                </div>
                <div class="family-badge">{pair_rows} pairs</div>
              </div>
              <div class="family-metrics">
                <div><span>Raw pairs</span><strong>{pair_rows}</strong></div>
                <div><span>Report rows</span><strong>{report_rows}</strong></div>
                <div><span>Heatmaps</span><strong>{heatmap_count}</strong></div>
              </div>
              <div class="family-links">
                <a href="{report_href}">Open report CSV</a>
              </div>
            </article>
            """
        )
    return "\n".join(cards)


def build_hero_summary(
    rows: list[dict[str, str]],
    quality_payload: dict[str, object] | None = None,
    bsa_payload: dict[str, object] | None = None,
    rmsf_payload: dict[str, object] | None = None,
) -> str:
    populated_rows = [row for row in rows if int(row["pair_rows"]) > 0]
    total_pairs = sum(int(row["pair_rows"]) for row in rows)
    top_family = max(rows, key=lambda row: int(row["pair_rows"]))["family"] if rows else "n/a"
    chips = [
        f"""
        <div class="hero-stat">
          <span>Active families</span>
          <strong>{len(populated_rows)}/{len(rows)}</strong>
        </div>
        """,
        f"""
        <div class="hero-stat">
          <span>Total residue pairs</span>
          <strong>{total_pairs}</strong>
        </div>
        """,
        f"""
        <div class="hero-stat">
          <span>Largest family</span>
          <strong>{html.escape(top_family)}</strong>
        </div>
        """,
    ]

    if quality_payload:
        metrics = quality_payload.get("metrics", {}) or {}
        chips.append(
            f"""
            <div class="hero-stat">
              <span>Tail-90% RMSD</span>
              <strong>{float(metrics.get('tail90_variation_nm', 0.0)):.3f} nm</strong>
            </div>
            """
        )
    if bsa_payload:
        summary = bsa_payload.get("summary", {}) or {}
        bsa_stats = summary.get("buried_surface_area", {}) or {}
        chips.append(
            f"""
            <div class="hero-stat">
              <span>Mean BSA</span>
              <strong>{float(bsa_stats.get('mean', 0.0)):.1f} Å²</strong>
            </div>
            """
        )
    if rmsf_payload:
        summary = rmsf_payload.get("summary", {}) or {}
        chips.append(
            f"""
            <div class="hero-stat">
              <span>Mean RMSF</span>
              <strong>{float(summary.get('mean_rmsf_angstrom', 0.0)):.2f} Å</strong>
            </div>
            """
        )

    return "\n".join(chips)


def choose_identity_root(cli_value: Path | None, source_pdb: Path) -> Path | None:
    candidates: list[Path] = []
    if cli_value is not None:
        candidates.append(cli_value)
    candidates.append(source_pdb.parent)
    for candidate in candidates:
        if not candidate:
            continue
        if (candidate / "analysis/identity/biological_identity.json").exists():
            return candidate
    return None


def build_biological_identity_assets(
    overview_dir: Path,
    source_pdb: Path,
    cdr_metadata_path: Path | None,
    identity_root: Path | None = None,
) -> dict[str, Any] | None:
    identity_dir = overview_dir / IDENTITY_DIR_NAME
    identity_dir.mkdir(exist_ok=True)

    json_path = identity_dir / "biological_identity.json"
    existing_root = choose_identity_root(identity_root, source_pdb)
    if existing_root is not None:
        source_json = existing_root / "analysis/identity/biological_identity.json"
        payload = json.loads(source_json.read_text(encoding="utf-8"))
        ensure_symlink(json_path, source_json)
    else:
        annotator = BiologicalIdentityAnnotator()
        payload = annotator.annotate(source_pdb, cdr_metadata_path)
        json_path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")

    payload["json_href"] = json_path.relative_to(overview_dir).as_posix()
    return payload


def normalize_visible_sections(
    include_sections: list[str] | None,
    exclude_sections: list[str] | None,
    quality_payload: dict[str, object] | None,
    bsa_payload: dict[str, object] | None,
    rmsf_payload: dict[str, object] | None,
    cluster_payload: dict[str, object] | None,
    occupancy_payload: list[dict[str, object]],
) -> set[str]:
    visible = set(include_sections or ALL_REPORT_SECTIONS)
    visible.add("overview")
    if exclude_sections:
        visible.difference_update(exclude_sections)
        visible.add("overview")

    if not quality_payload:
        visible.discard("quality")
    if not bsa_payload:
        visible.discard("interface")
    if not rmsf_payload:
        visible.discard("flexibility")
    if not cluster_payload:
        visible.discard("cluster")
    if not occupancy_payload:
        visible.discard("occupancy")
    return visible


def optional_block(enabled: bool, content: str) -> str:
    return content if enabled else ""


def build_hero_identity(identity_payload: dict[str, Any] | None) -> str:
    if not identity_payload:
        return ""

    complex_identity = identity_payload.get("complex_identity", {}) or {}
    peptide_identity = identity_payload.get("peptide_identity", {}) or {}
    tcr_identity = identity_payload.get("tcr_identity", {}) or {}
    hla_identity = identity_payload.get("hla_identity", {}) or {}

    hla_title = str(hla_identity.get("best_locus", "") or "MHC class I")
    hla_subtitle = str(hla_identity.get("best_candidate_allele", "") or "Best candidate unavailable")
    peptide_seq = str(peptide_identity.get("sequence", "") or "n/a")
    peptide_len = int(peptide_identity.get("length", 0) or 0)
    alpha_v = str(tcr_identity.get("alpha_v_gene", "") or "")
    alpha_j = str(tcr_identity.get("alpha_j_gene", "") or "")
    beta_v = str(tcr_identity.get("beta_v_gene", "") or "")
    beta_j = str(tcr_identity.get("beta_j_gene", "") or "")
    cdr3_alpha = str(tcr_identity.get("cdr3_alpha_sequence", "") or "n/a")
    cdr3_beta = str(tcr_identity.get("cdr3_beta_sequence", "") or "n/a")
    peptide_chain = str(complex_identity.get("peptide_chain", "") or "n/a")
    alpha_genotype = " / ".join(part for part in (alpha_v, alpha_j) if part) or "Genotype unavailable"
    beta_genotype = " / ".join(part for part in (beta_v, beta_j) if part) or "Genotype unavailable"

    cards = [
        f"""
        <article class="hero-identity-card">
          <span>HLA identity</span>
          <strong>{html.escape(hla_title)}</strong>
          <p>{html.escape(hla_subtitle)}</p>
        </article>
        """,
        f"""
        <article class="hero-identity-card">
          <span>Peptide</span>
          <strong>{html.escape(peptide_seq)}</strong>
          <p>{peptide_len} aa peptide</p>
        </article>
        """,
        f"""
        <article class="hero-identity-card">
          <span>TCR α</span>
          <strong>{html.escape(alpha_genotype)}</strong>
          <p>CDR3α: {html.escape(cdr3_alpha)}</p>
        </article>
        """,
        f"""
        <article class="hero-identity-card">
          <span>TCR β</span>
          <strong>{html.escape(beta_genotype)}</strong>
          <p>CDR3β: {html.escape(cdr3_beta)}</p>
        </article>
        """,
    ]
    return "\n".join(cards)


def build_hero_metric_strip(
    rows: list[dict[str, str]],
    quality_payload: dict[str, object] | None = None,
    bsa_payload: dict[str, object] | None = None,
) -> str:
    populated_rows = [row for row in rows if int(row["pair_rows"]) > 0]
    total_pairs = sum(int(row["pair_rows"]) for row in rows)
    chips = [
        f"""
        <div class="hero-stat">
          <span>Active families</span>
          <strong>{len(populated_rows)}/{len(rows)}</strong>
        </div>
        """,
        f"""
        <div class="hero-stat">
          <span>Total residue pairs</span>
          <strong>{total_pairs}</strong>
        </div>
        """,
    ]
    if quality_payload:
        metrics = quality_payload.get("metrics", {}) or {}
        chips.append(
            f"""
            <div class="hero-stat">
              <span>Tail-90% RMSD</span>
              <strong>{float(metrics.get('tail90_variation_nm', 0.0)):.3f} nm</strong>
            </div>
            """
        )
    if bsa_payload:
        summary = bsa_payload.get("summary", {}) or {}
        bsa_stats = summary.get("buried_surface_area", {}) or {}
        chips.append(
            f"""
            <div class="hero-stat">
              <span>Mean BSA</span>
              <strong>{float(bsa_stats.get('mean', 0.0)):.1f} Å²</strong>
            </div>
            """
        )
    return "\n".join(chips)


def build_table(rows: list[dict[str, str]]) -> str:
    lines = [
        "<table>",
        "<thead><tr><th>family</th><th>type</th><th>raw pairs</th><th>report rows</th><th>report</th></tr></thead>",
        "<tbody>",
    ]
    for row in rows:
        lines.append(
            "<tr>"
            f"<td>{html.escape(row['family'])}</td>"
            f"<td>{html.escape(row['label'])}</td>"
            f"<td>{html.escape(row['pair_rows'])}</td>"
            f"<td>{html.escape(row['report_rows'])}</td>"
            f"<td><a href=\"{html.escape(row['overview_report'])}\">Open</a></td>"
            "</tr>"
        )
    lines.append("</tbody></table>")
    return "\n".join(lines)


def build_heatmap_sections(rows: list[dict[str, str]]) -> str:
    sections = []
    for row in rows:
        if not row["heatmaps"]:
            continue
        images = []
        for item in row["heatmaps"].split("; "):
            label = Path(item).name.replace(".png", "")
            images.append(
                f"""
                <figure class="heatmap-card">
                  <img src="{html.escape(item)}" alt="{html.escape(label)}">
                  <figcaption>{html.escape(label)}</figcaption>
                </figure>
                """
            )
        sections.append(
            f"""
            <section class="heatmap-section">
              <div class="section-title">
                <h3>{html.escape(row['family'])}</h3>
                <p>{html.escape(row['label'])}</p>
              </div>
              <div class="heatmap-grid">
                {''.join(images)}
              </div>
            </section>
            """
        )
    return "\n".join(sections)


def load_rmsd_xvg(xvg_path: Path) -> pd.DataFrame:
    records: list[tuple[float, float]] = []
    with xvg_path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#") or stripped.startswith("@"):
                continue
            parts = stripped.split()
            if len(parts) < 2:
                continue
            try:
                time_ps = float(parts[0])
                rmsd_nm = float(parts[1])
            except ValueError:
                continue
            records.append((time_ps, rmsd_nm))
    return pd.DataFrame(records, columns=["time_ps", "rmsd_nm"])


def build_quality_assets(
    overview_dir: Path,
    figures_dir: Path,
    quality_root: Path | None,
) -> dict[str, object] | None:
    if quality_root is None:
        return None

    quality_json = quality_root / "quality/preprocess_quality_report.json"
    quality_md = quality_root / "quality/preprocess_quality_report.md"
    rmsd_xvg = quality_root / "analysis/rmsd/rmsd.xvg"
    rmsd_png = quality_root / "plots/rmsd/rmsd.png"
    if not (quality_json.exists() and rmsd_xvg.exists()):
        return None

    quality_dir = overview_dir / QUALITY_DIR_NAME
    quality_dir.mkdir(exist_ok=True)

    json_link = quality_dir / "preprocess_quality_report.json"
    xvg_link = quality_dir / "rmsd.xvg"
    ensure_symlink(json_link, quality_json)
    ensure_symlink(xvg_link, rmsd_xvg)
    md_link = None
    png_link = None
    if quality_md.exists():
        md_link = quality_dir / "preprocess_quality_report.md"
        ensure_symlink(md_link, quality_md)
    if rmsd_png.exists():
        png_link = quality_dir / "rmsd.png"
        ensure_symlink(png_link, rmsd_png)

    metrics = json.loads(quality_json.read_text(encoding="utf-8"))
    rmsd_df = load_rmsd_xvg(rmsd_xvg)
    if rmsd_df.empty:
        return None

    trim_fraction = float(metrics.get("trim_fraction", 0.1) or 0.1)
    trim_index = min(max(int(len(rmsd_df) * trim_fraction), 0), max(len(rmsd_df) - 1, 0))
    tail_df = rmsd_df.iloc[trim_index:].copy()
    tail_start_time = float(metrics.get("tail90_start_time_ps", tail_df["time_ps"].iloc[0]))
    tail_end_time = float(metrics.get("tail90_end_time_ps", tail_df["time_ps"].iloc[-1]))

    # 主图：全程 RMSD，突出后 90%
    full_plot = figures_dir / "quality_rmsd_overview.png"
    fig, ax = plt.subplots(figsize=(10.2, 4.8))
    ax.axvspan(rmsd_df["time_ps"].iloc[0], tail_start_time, color="#efe2c7", alpha=0.7, linewidth=0)
    ax.axvspan(tail_start_time, tail_end_time, color="#dfece5", alpha=0.95, linewidth=0)
    ax.plot(rmsd_df["time_ps"], rmsd_df["rmsd_nm"], color="#658c7c", linewidth=1.4, alpha=0.9)
    ax.plot(tail_df["time_ps"], tail_df["rmsd_nm"], color="#2f5d50", linewidth=1.8)
    ax.set_title("RMSD trajectory", fontsize=14, weight="bold")
    ax.set_xlabel("Time (ps)")
    ax.set_ylabel("RMSD (nm)")
    ax.grid(alpha=0.18, linestyle="--")
    ax.set_axisbelow(True)
    fig.tight_layout()
    fig.savefig(full_plot, dpi=220, facecolor="white")
    plt.close(fig)

    # 辅助图：后 90% 局部放大
    tail_plot = figures_dir / "quality_rmsd_tail90.png"
    fig, ax = plt.subplots(figsize=(10.2, 4.6))
    ax.plot(tail_df["time_ps"], tail_df["rmsd_nm"], color="#2f5d50", linewidth=1.7)
    ax.fill_between(tail_df["time_ps"], tail_df["rmsd_nm"], color="#c8d8cf", alpha=0.45)
    ax.set_title("Tail-90% RMSD focus", fontsize=14, weight="bold")
    ax.set_xlabel("Time (ps)")
    ax.set_ylabel("RMSD (nm)")
    ax.grid(alpha=0.18, linestyle="--")
    ax.set_axisbelow(True)
    fig.tight_layout()
    fig.savefig(tail_plot, dpi=220, facecolor="white")
    plt.close(fig)

    return {
        "metrics": metrics,
        "json_href": json_link.relative_to(overview_dir).as_posix(),
        "xvg_href": xvg_link.relative_to(overview_dir).as_posix(),
        "md_href": md_link.relative_to(overview_dir).as_posix() if md_link else "",
        "png_href": png_link.relative_to(overview_dir).as_posix() if png_link else "",
        "overview_plot": full_plot.name,
        "tail_plot": tail_plot.name,
    }


def build_nav(visible_sections: set[str]) -> str:
    items = []
    for section_id, label in [
        ("overview", "Overview"),
        ("quality", "Quality"),
        ("interface", "Interface"),
        ("flexibility", "Flexibility"),
        ("cluster", "Cluster"),
        ("occupancy", "Persistence"),
        ("contact", "Landscape"),
        ("interactions", "Residue Pairs"),
        ("downloads", "Downloads"),
    ]:
        if section_id in visible_sections:
            items.append(f'<a href="#{section_id}">{label}</a>')
    return "".join(items)


def build_bsa_assets(
    overview_dir: Path,
    bsa_root: Path | None,
) -> dict[str, object] | None:
    if bsa_root is None:
        return None

    summary_json = bsa_root / "analysis/interface/interface_summary.json"
    timeseries_csv = bsa_root / "analysis/interface/bsa_timeseries.csv"
    plot_png = bsa_root / "analysis/interface/bsa_timeseries.png"
    if not (summary_json.exists() and timeseries_csv.exists()):
        return None

    bsa_dir = overview_dir / BSA_DIR_NAME
    bsa_dir.mkdir(exist_ok=True)

    summary_link = bsa_dir / "interface_summary.json"
    csv_link = bsa_dir / "bsa_timeseries.csv"
    ensure_symlink(summary_link, summary_json)
    ensure_symlink(csv_link, timeseries_csv)

    plot_href = ""
    if plot_png.exists():
        plot_link = bsa_dir / "bsa_timeseries.png"
        ensure_symlink(plot_link, plot_png)
        plot_href = plot_link.relative_to(overview_dir).as_posix()

    summary = json.loads(summary_json.read_text(encoding="utf-8"))
    return {
        "summary": summary,
        "json_href": summary_link.relative_to(overview_dir).as_posix(),
        "csv_href": csv_link.relative_to(overview_dir).as_posix(),
        "plot_href": plot_href,
    }


def build_bsa_section(bsa_payload: dict[str, object] | None) -> str:
    if not bsa_payload:
        return ""

    summary = bsa_payload["summary"]
    bsa_stats = summary.get("buried_surface_area", {}) or {}
    ratio_stats = summary.get("interface_ratio", {}) or {}
    selection_a = str(summary.get("selection_a", ""))
    selection_b = str(summary.get("selection_b", ""))
    n_frames = int(summary.get("n_frames", 0) or 0)
    time_end = float(summary.get("time_end", 0.0) or 0.0)

    return f"""
    <section class="section" id="interface">
      <h2>Interface Burial</h2>
      <p class="note">This section summarizes buried surface area between the pHLA-side assembly and the TCR-side assembly, highlighting whether the interface stays broadly compact across the sampled trajectory.</p>
      <div class="quality-grid">
        <div class="quality-card">
          <span>Frames analysed</span>
          <strong>{n_frames}</strong>
        </div>
        <div class="quality-card">
          <span>Time span</span>
          <strong>{time_end/1000:.1f} ns</strong>
        </div>
        <div class="quality-card">
          <span>Mean buried surface area</span>
          <strong>{float(bsa_stats.get('mean', 0.0)):.1f} Å²</strong>
        </div>
        <div class="quality-card">
          <span>Mean interface ratio</span>
          <strong>{float(ratio_stats.get('mean', 0.0)):.3f}</strong>
        </div>
      </div>
      <div class="quality-note-box">
        <p><strong>Selection A.</strong> <code>{html.escape(selection_a)}</code></p>
        <p><strong>Selection B.</strong> <code>{html.escape(selection_b)}</code></p>
        <p><strong>Reading guide.</strong> Higher buried surface area and a tighter interface-ratio band usually indicate a more consistently packed interface, while broader oscillations suggest larger interface breathing.</p>
        <div class="quality-downloads">
          <a href="{html.escape(str(bsa_payload['json_href']))}">Open BSA summary JSON</a>
          <a href="{html.escape(str(bsa_payload['csv_href']))}">Open BSA timeseries CSV</a>
        </div>
      </div>
      <div class="insight-grid quality-figures">
        <figure class="insight-card">
          <img src="{html.escape(str(bsa_payload['plot_href']))}" alt="BSA timeseries">
          <figcaption>
            <strong>Buried surface area and interface ratio</strong>
            <span>The upper panel tracks buried surface area in Å², while the lower panel shows the normalized interface ratio over the same time axis.</span>
          </figcaption>
        </figure>
      </div>
    </section>
    """


def build_rmsf_assets(
    overview_dir: Path,
    rmsf_root: Path | None,
) -> dict[str, object] | None:
    if rmsf_root is None:
        return None

    residue_csv = rmsf_root / "analysis/rmsf/residue_rmsf.csv"
    region_csv = rmsf_root / "analysis/rmsf/region_rmsf_summary.csv"
    summary_json = rmsf_root / "analysis/rmsf/rmsf_summary.json"
    tcr_plot = rmsf_root / "analysis/rmsf/tcr_rmsf_profile.png"
    phla_plot = rmsf_root / "analysis/rmsf/phla_rmsf_profile.png"
    region_plot = rmsf_root / "analysis/rmsf/region_rmsf_summary.png"
    if not (residue_csv.exists() and region_csv.exists() and summary_json.exists()):
        return None

    rmsf_dir = overview_dir / RMSF_DIR_NAME
    rmsf_dir.mkdir(exist_ok=True)
    residue_link = rmsf_dir / "residue_rmsf.csv"
    region_link = rmsf_dir / "region_rmsf_summary.csv"
    summary_link = rmsf_dir / "rmsf_summary.json"
    ensure_symlink(residue_link, residue_csv)
    ensure_symlink(region_link, region_csv)
    ensure_symlink(summary_link, summary_json)

    tcr_href = ""
    if tcr_plot.exists():
        tcr_link = rmsf_dir / "tcr_rmsf_profile.png"
        ensure_symlink(tcr_link, tcr_plot)
        tcr_href = tcr_link.relative_to(overview_dir).as_posix()
    phla_href = ""
    if phla_plot.exists():
        phla_link = rmsf_dir / "phla_rmsf_profile.png"
        ensure_symlink(phla_link, phla_plot)
        phla_href = phla_link.relative_to(overview_dir).as_posix()
    region_href = ""
    if region_plot.exists():
        region_plot_link = rmsf_dir / "region_rmsf_summary.png"
        ensure_symlink(region_plot_link, region_plot)
        region_href = region_plot_link.relative_to(overview_dir).as_posix()

    residue_df = pd.read_csv(residue_csv)
    region_df = pd.read_csv(region_csv)
    summary = json.loads(summary_json.read_text(encoding="utf-8"))
    return {
        "residue_href": residue_link.relative_to(overview_dir).as_posix(),
        "region_href": region_link.relative_to(overview_dir).as_posix(),
        "summary_href": summary_link.relative_to(overview_dir).as_posix(),
        "tcr_plot_href": tcr_href,
        "phla_plot_href": phla_href,
        "region_plot_href": region_href,
        "residue_df": residue_df,
        "region_df": region_df,
        "summary": summary,
    }


def build_rmsf_section(rmsf_payload: dict[str, object] | None) -> str:
    if not rmsf_payload:
        return ""
    summary = rmsf_payload["summary"]
    region_df = rmsf_payload["region_df"]
    focus_regions = region_df[region_df["region_group"].isin([
        "CDR1_alpha", "CDR2_alpha", "CDR3_alpha", "CDR1_beta", "CDR2_beta", "CDR3_beta",
        "peptide", "alpha1_helix", "alpha2_helix",
    ])].head(9).copy()
    region_rows = []
    for _, row in focus_regions.iterrows():
        region_rows.append(
            "<tr>"
            f"<td>{html.escape(str(row['region_group']))}</td>"
            f"<td>{int(row['n_residues'])}</td>"
            f"<td>{float(row['mean_rmsf_angstrom']):.2f}</td>"
            f"<td>{float(row['max_rmsf_angstrom']):.2f}</td>"
            "</tr>"
        )
    region_table = (
        '<div class="mini-table-wrap"><div class="mini-table-title">Region RMSF snapshot</div>'
        '<table class="mini-table"><thead><tr><th>region</th><th>n</th><th>mean Å</th><th>max Å</th></tr></thead>'
        f"<tbody>{''.join(region_rows)}</tbody></table></div>"
    )
    return f"""
    <section class="section" id="flexibility">
      <h2>Flexibility</h2>
      <p class="note">This section captures where the pHLA–TCR complex remains rigid or flexible, complementing global RMSD with residue- and region-level fluctuation.</p>
      <div class="quality-grid">
        <div class="quality-card">
          <span>Residues analysed</span>
          <strong>{int(summary.get('n_residues', 0))}</strong>
        </div>
        <div class="quality-card">
          <span>Frames analysed</span>
          <strong>{int(summary.get('n_frames', 0))}</strong>
        </div>
        <div class="quality-card">
          <span>Mean RMSF</span>
          <strong>{float(summary.get('mean_rmsf_angstrom', 0.0)):.2f} Å</strong>
        </div>
        <div class="quality-card">
          <span>Max RMSF</span>
          <strong>{float(summary.get('max_rmsf_angstrom', 0.0)):.2f} Å</strong>
        </div>
      </div>
      <div class="insight-grid quality-figures">
        <figure class="insight-card">
          <img src="{html.escape(str(rmsf_payload['tcr_plot_href']))}" alt="TCR RMSF profile">
          <figcaption>
            <strong>TCR RMSF profile</strong>
            <span>TCR alpha and beta flexibility profiles reveal whether motion is concentrated in a specific receptor arm.</span>
          </figcaption>
        </figure>
        <figure class="insight-card">
          <img src="{html.escape(str(rmsf_payload['phla_plot_href']))}" alt="pHLA RMSF profile">
          <figcaption>
            <strong>pHLA RMSF profile</strong>
            <span>The peptide and MHC alpha chain are shown separately so groove-side fluctuation is not hidden by peptide motion.</span>
          </figcaption>
        </figure>
        <figure class="insight-card">
          <img src="{html.escape(str(rmsf_payload['region_plot_href']))}" alt="Region RMSF summary">
          <figcaption>
            <strong>Region RMSF summary</strong>
            <span>CDR and pHLA subregions are summarized side by side to connect flexibility with interaction hotspots and occupancy.</span>
          </figcaption>
        </figure>
      </div>
      <div class="stability-table-grid">
        {region_table}
      </div>
      <div class="quality-downloads">
        <a href="{html.escape(str(rmsf_payload['summary_href']))}">Open RMSF summary JSON</a>
        <a href="{html.escape(str(rmsf_payload['residue_href']))}">Open residue RMSF CSV</a>
        <a href="{html.escape(str(rmsf_payload['region_href']))}">Open region RMSF CSV</a>
      </div>
    </section>
    """


def build_cluster_assets(
    overview_dir: Path,
    cluster_root: Path | None,
) -> dict[str, object] | None:
    if cluster_root is None:
        return None

    analysis_dir = cluster_root / "analysis/conformation/interface_clustering"
    summary_csv = analysis_dir / "summary_table.csv"
    digest_csv = analysis_dir / "cluster_feature_digest.csv"
    difference_csv = analysis_dir / "cluster_difference_matrix.csv"
    assignments_csv = analysis_dir / "frame_cluster_assignments.csv"
    summary_json = analysis_dir / "interface_clustering_summary.json"
    cluster_id_png = analysis_dir / "cluster_id_vs_time.png"
    population_png = analysis_dir / "state_population_over_time.png"
    structures_dir = analysis_dir / "structures"
    if not (
        summary_csv.exists()
        and digest_csv.exists()
        and summary_json.exists()
        and cluster_id_png.exists()
        and population_png.exists()
    ):
        return None

    cluster_dir = overview_dir / CLUSTER_DIR_NAME
    cluster_dir.mkdir(exist_ok=True)

    summary_link = cluster_dir / "summary_table.csv"
    digest_link = cluster_dir / "cluster_feature_digest.csv"
    difference_link = cluster_dir / "cluster_difference_matrix.csv"
    assignments_link = cluster_dir / "frame_cluster_assignments.csv"
    summary_json_link = cluster_dir / "interface_clustering_summary.json"
    cluster_id_link = cluster_dir / "cluster_id_vs_time.png"
    population_link = cluster_dir / "state_population_over_time.png"
    structures_link = cluster_dir / "structures"

    ensure_symlink(summary_link, summary_csv)
    ensure_symlink(digest_link, digest_csv)
    if difference_csv.exists():
        ensure_symlink(difference_link, difference_csv)
    if assignments_csv.exists():
        ensure_symlink(assignments_link, assignments_csv)
    ensure_symlink(summary_json_link, summary_json)
    ensure_symlink(cluster_id_link, cluster_id_png)
    ensure_symlink(population_link, population_png)
    if structures_link.exists() or structures_link.is_symlink():
        structures_link.unlink()
    structures_link.symlink_to(structures_dir.resolve(), target_is_directory=True)

    summary_df = pd.read_csv(summary_csv)
    digest_df = pd.read_csv(digest_csv)
    summary = json.loads(summary_json.read_text(encoding="utf-8"))
    return {
        "summary_df": summary_df,
        "digest_df": digest_df,
        "summary": summary,
        "summary_href": summary_link.relative_to(overview_dir).as_posix(),
        "digest_href": digest_link.relative_to(overview_dir).as_posix(),
        "difference_href": difference_link.relative_to(overview_dir).as_posix() if difference_csv.exists() else "",
        "assignments_href": assignments_link.relative_to(overview_dir).as_posix() if assignments_csv.exists() else "",
        "summary_json_href": summary_json_link.relative_to(overview_dir).as_posix(),
        "cluster_id_plot_href": cluster_id_link.relative_to(overview_dir).as_posix(),
        "population_plot_href": population_link.relative_to(overview_dir).as_posix(),
        "structures_href": structures_link.relative_to(overview_dir).as_posix(),
    }


def build_cluster_section(cluster_payload: dict[str, object] | None) -> str:
    if not cluster_payload:
        return ""

    summary_df = cluster_payload["summary_df"]
    digest_df = cluster_payload["digest_df"]
    summary = cluster_payload["summary"]
    largest = summary.get("largest_cluster", {}) or {}

    metric_cards = [
        ("Sampled frames", int(summary.get("n_frames", 0) or 0)),
        ("Clusters", int(summary.get("n_clusters", 0) or 0)),
        ("Largest cluster", f"{float(largest.get('population_percent', 0.0)):.1f}%"),
        ("Contact stride", int(summary.get("contact_stride", 0) or 0)),
    ]
    metric_html = "".join(
        f"""
        <div class="quality-card">
          <span>{html.escape(str(label))}</span>
          <strong>{html.escape(str(value))}</strong>
        </div>
        """
        for label, value in metric_cards
    )

    summary_head = summary_df.copy().head(8) if not summary_df.empty else pd.DataFrame()
    summary_rows = []
    for _, row in summary_head.iterrows():
        summary_rows.append(
            "<tr>"
            f"<td>{int(row['cluster_id'])}</td>"
            f"<td>{float(row['population_percent']):.1f}</td>"
            f"<td>{int(row['frames'])}</td>"
            f"<td>{int(row['representative_frame'])}</td>"
            f"<td>{int(row['medoid_frame'])}</td>"
            f"<td>{float(row['mean_dwell_time_ps']):.0f}</td>"
            f"<td>{html.escape(str(row['main_structural_descriptor']))}</td>"
            "</tr>"
        )
    summary_table = (
        '<div class="mini-table-wrap"><div class="mini-table-title">Cluster summary</div>'
        '<table class="mini-table"><thead><tr>'
        '<th>ID</th><th>Pop. %</th><th>Frames</th><th>Rep.</th><th>Medoid</th><th>Dwell (ps)</th><th>Main descriptor</th>'
        f"</tr></thead><tbody>{''.join(summary_rows)}</tbody></table></div>"
        if summary_rows
        else '<div class="mini-table-empty">No cluster summary records.</div>'
    )

    digest_blocks = []
    if not digest_df.empty:
        for cluster_id, group in digest_df.groupby("cluster_id"):
            key_rows = []
            sig_rows = []
            for _, row in group[group["digest_type"] == "key_feature"].sort_values("rank").head(5).iterrows():
                pair_display = str(row.get("pair_display", row.get("pair_label", "-")))
                pair_description = str(row.get("pair_description", ""))
                key_rows.append(
                    "<tr>"
                    f"<td><div>{html.escape(pair_display)}</div><small>{html.escape(pair_description)}</small></td>"
                    f"<td>{float(row['contact_occupancy']):.2f}</td>"
                    f"<td>{float(row['mean_pair_distance_angstrom']):.2f}</td>"
                    "</tr>"
                )
            for _, row in group[group["digest_type"] == "state_signature"].sort_values("rank").head(5).iterrows():
                pair_display = str(row.get("pair_display", row.get("pair_label", "-")))
                pair_description = str(row.get("pair_description", ""))
                sig_rows.append(
                    "<tr>"
                    f"<td><div>{html.escape(pair_display)}</div><small>{html.escape(pair_description)}</small></td>"
                    f"<td>{float(row['score']):.2f}</td>"
                    f"<td>{float(row['contact_occupancy']):.2f}</td>"
                    "</tr>"
                )
            digest_blocks.append(
                f"""
                <article class="occupancy-card">
                  <div class="family-detail-head">
                    <div>
                      <h3>Cluster {int(cluster_id)}</h3>
                      <p>Compact interface signature digest</p>
                    </div>
                  </div>
                  <div class="stability-table-grid">
                    <div class="mini-table-wrap">
                      <div class="mini-table-title">Key interface features</div>
                      <table class="mini-table">
                        <thead><tr><th>pair</th><th>occ.</th><th>dist. Å</th></tr></thead>
                        <tbody>{''.join(key_rows) or '<tr><td colspan="3">No key features.</td></tr>'}</tbody>
                      </table>
                    </div>
                    <div class="mini-table-wrap">
                      <div class="mini-table-title">Cluster-specific signatures</div>
                      <table class="mini-table">
                        <thead><tr><th>pair</th><th>score</th><th>occ.</th></tr></thead>
                        <tbody>{''.join(sig_rows) or '<tr><td colspan="3">No signatures.</td></tr>'}</tbody>
                      </table>
                    </div>
                  </div>
                </article>
                """
            )

    return f"""
    <section class="section" id="cluster">
      <h2>Interface Clustering</h2>
      <p class="note">This section summarizes the dominant pHLA–TCR interface states, how they evolve over time, and which residue pairs best distinguish one state from another.</p>
      <div class="quality-grid">
        {metric_html}
      </div>
      <div class="insight-grid quality-figures">
        <figure class="insight-card">
          <img src="{html.escape(str(cluster_payload['cluster_id_plot_href']))}" alt="Cluster ID vs time">
          <figcaption>
            <strong>Cluster ID vs time</strong>
            <span>Shows when the trajectory occupies each interface state and whether switches are rare or frequent.</span>
          </figcaption>
        </figure>
        <figure class="insight-card">
          <img src="{html.escape(str(cluster_payload['population_plot_href']))}" alt="State population over time">
          <figcaption>
            <strong>State population trends</strong>
            <span>Rolling and cumulative population curves make the dominant state and slower redistribution easier to read.</span>
          </figcaption>
        </figure>
      </div>
      <div class="stability-table-grid">
        {summary_table}
      </div>
      <div class="stability-grid">
        {''.join(digest_blocks)}
      </div>
      <div class="quality-downloads">
        <a href="{html.escape(str(cluster_payload['summary_href']))}">Open summary table CSV</a>
        <a href="{html.escape(str(cluster_payload['digest_href']))}">Open cluster digest CSV</a>
        <a href="{html.escape(str(cluster_payload['difference_href']))}">Open difference matrix CSV</a>
        <a href="{html.escape(str(cluster_payload['assignments_href']))}">Open frame assignments CSV</a>
        <a href="{html.escape(str(cluster_payload['summary_json_href']))}">Open clustering summary JSON</a>
        <a href="{html.escape(str(cluster_payload['structures_href']))}">Open cluster structures</a>
      </div>
    </section>
    """


def _download_entry(href: str, title: str, subtitle: str) -> str:
    return (
        f'<a href="{html.escape(href)}">{html.escape(title)}'
        f"<small>{html.escape(subtitle)}</small></a>"
    )


def build_downloads_section(
    identity_payload: dict[str, object] | None,
    bsa_payload: dict[str, object] | None,
    rmsf_payload: dict[str, object] | None,
    cluster_payload: dict[str, object] | None,
) -> str:
    groups: list[tuple[str, str, list[str]]] = [
        (
            "Overview",
            "Unified report-level index files and summary payloads.",
            [
                _download_entry(
                    "interaction_overview.csv",
                    "interaction_overview.csv",
                    "Unified table across all interaction families",
                ),
                _download_entry(
                    "interaction_overview.json",
                    "interaction_overview.json",
                    "Machine-readable overview payload",
                ),
            ],
        ),
        (
            "Interaction Families",
            "Primary coarse-contact and typed-interaction reports.",
            [
                _download_entry("reports/contact_report.csv", "contact_report.csv", "Primary coarse-contact report"),
                _download_entry("reports/hbond_report.csv", "hbond_report.csv", "Primary hydrogen-bond report"),
                _download_entry("reports/saltbridge_report.csv", "saltbridge_report.csv", "Primary salt-bridge report"),
                _download_entry("reports/hydrophobic_report.csv", "hydrophobic_report.csv", "Primary hydrophobic-contact report"),
                _download_entry("reports/pipi_report.csv", "pipi_report.csv", "Primary pi-pi report"),
                _download_entry("reports/cationpi_report.csv", "cationpi_report.csv", "Primary cation-pi report"),
            ],
        ),
    ]

    if identity_payload:
        groups.append(
            (
                "Biological Identity",
                "Peptide, HLA and TCR identity annotations.",
                [
                    _download_entry(
                        str(identity_payload["json_href"]),
                        "biological_identity.json",
                        "Biological identity annotation for peptide, HLA and TCR chains",
                    )
                ],
            )
        )
    if bsa_payload:
        groups.append(
            (
                "Interface Burial",
                "Buried surface area and interface packing outputs.",
                [
                    _download_entry(
                        str(bsa_payload["json_href"]),
                        "interface_summary.json",
                        "BSA and interface-ratio summary for the pHLA–TCR interface",
                    ),
                    _download_entry(
                        str(bsa_payload["csv_href"]),
                        "bsa_timeseries.csv",
                        "Frame-wise buried surface area and interface ratio",
                    ),
                ],
            )
        )
    if rmsf_payload:
        groups.append(
            (
                "Flexibility",
                "Residue- and region-level RMSF outputs.",
                [
                    _download_entry(
                        str(rmsf_payload["summary_href"]),
                        "rmsf_summary.json",
                        "Residue- and region-level RMSF summary",
                    ),
                    _download_entry(
                        str(rmsf_payload["residue_href"]),
                        "residue_rmsf.csv",
                        "Annotated residue-level RMSF table",
                    ),
                    _download_entry(
                        str(rmsf_payload["region_href"]),
                        "region_rmsf_summary.csv",
                        "Grouped RMSF summary across CDRs, peptide and groove helices",
                    ),
                ],
            )
        )
    if cluster_payload:
        groups.append(
            (
                "Interface Clustering",
                "Cluster assignments, state summaries and cluster-specific signatures.",
                [
                    _download_entry(
                        str(cluster_payload["summary_href"]),
                        "cluster_summary_table.csv",
                        "Interface-state summary across all detected clusters",
                    ),
                    _download_entry(
                        str(cluster_payload["digest_href"]),
                        "cluster_feature_digest.csv",
                        "Compact per-cluster digest of key interface features and signatures",
                    ),
                    _download_entry(
                        str(cluster_payload["summary_json_href"]),
                        "interface_clustering_summary.json",
                        "Interface clustering metadata, weights and dominant-state summary",
                    ),
                ],
            )
        )

    blocks = []
    for idx, (title, desc, entries) in enumerate(groups):
        blocks.append(
            f"""
            <details class="download-group" {'open' if idx < 2 else ''}>
              <summary>
                <span>{html.escape(title)}</span>
                <small>{html.escape(desc)}</small>
              </summary>
              <div class="download-group-grid">
                {''.join(entries)}
              </div>
            </details>
            """
        )

    return f"""
    <section class="section" id="downloads">
      <h2>Downloads</h2>
      <p class="note">Outputs are grouped by analysis domain so you can scan the directory like an index instead of reading a long flat list.</p>
      <div class="download-directory">
        {''.join(blocks)}
      </div>
    </section>
    """


def build_quality_section(quality_payload: dict[str, object] | None) -> str:
    if not quality_payload:
        return ""
    metrics = quality_payload["metrics"]
    frames = int(metrics.get("n_frames", 0) or 0)
    time_end = float(metrics.get("tail90_end_time_ps", 0.0) or 0.0)
    full_var = float(metrics.get("full_variation_nm", 0.0) or 0.0)
    tail_var = float(metrics.get("tail90_variation_nm", 0.0) or 0.0)
    tail_start = float(metrics.get("tail90_start_time_ps", 0.0) or 0.0)
    tail_mean = float(metrics.get("tail90_mean_rmsd_nm", 0.0) or 0.0)
    tail_std = float(metrics.get("tail90_std_rmsd_nm", 0.0) or 0.0)
    tail_min = float(metrics.get("tail90_min_rmsd_nm", 0.0) or 0.0)
    tail_max = float(metrics.get("tail90_max_rmsd_nm", 0.0) or 0.0)
    ratio = (tail_var / full_var) if full_var > 0 else 0.0
    if ratio <= 0.45:
        verdict = "The trajectory relaxes early and becomes substantially narrower in the final 90%."
        badge = "Stable tail"
    elif ratio <= 0.7:
        verdict = "The tail segment is clearly tighter than the full trace, but a visible fluctuation band remains."
        badge = "Moderate tail variation"
    else:
        verdict = "The tail segment still spans a broad RMSD range, so downstream interaction changes should be read cautiously."
        badge = "Broad tail variation"
    md_link = quality_payload.get("md_href", "")
    md_anchor = (
        f'<a href="{html.escape(str(md_link))}">Open quality note</a>' if md_link else ""
    )
    return f"""
    <section class="section" id="quality">
      <h2>Trajectory Quality</h2>
      <p class="note">This section sits before contact interpretation so the reader can judge whether the downstream interaction patterns are grounded in a reasonably behaved trajectory.</p>
      <div class="quality-banner">
        <div class="quality-banner-badge">{badge}</div>
        <div class="quality-banner-body">
          <strong>Quality reading</strong>
          <span>{verdict}</span>
        </div>
      </div>
      <div class="quality-grid">
        <div class="quality-card">
          <span>Frames</span>
          <strong>{frames}</strong>
        </div>
        <div class="quality-card">
          <span>Time span</span>
          <strong>{time_end/1000:.1f} ns</strong>
        </div>
        <div class="quality-card">
          <span>Full RMSD variation</span>
          <strong>{full_var:.3f} nm</strong>
        </div>
        <div class="quality-card">
          <span>Tail-90% variation</span>
          <strong>{tail_var:.3f} nm</strong>
        </div>
      </div>
      <div class="insight-grid quality-figures">
        <figure class="insight-card">
          <img src="{html.escape(FIGURES_DIR_NAME + '/' + str(quality_payload['overview_plot']))}" alt="RMSD trajectory overview">
          <figcaption>
            <strong>RMSD trajectory overview</strong>
            <span>The first 10% is shaded separately to make the later, more informative 90% easier to read.</span>
          </figcaption>
        </figure>
        <figure class="insight-card">
          <img src="{html.escape(FIGURES_DIR_NAME + '/' + str(quality_payload['tail_plot']))}" alt="Tail-90% RMSD focus">
          <figcaption>
            <strong>Tail-90% RMSD focus</strong>
            <span>Downstream trust should be driven primarily by the fluctuation range after the initial relaxation period.</span>
          </figcaption>
        </figure>
      </div>
      <div class="quality-note-box">
        <p><strong>Reading guide.</strong> This report does not rank trajectories by an arbitrary grade. Instead, it foregrounds the RMSD fluctuation range and the post-relaxation tail segment.</p>
        <p><strong>Current observation.</strong> The final 90% starts at <strong>{tail_start/1000:.1f} ns</strong>, with a mean RMSD of <strong>{tail_mean:.3f} nm</strong>, a standard deviation of <strong>{tail_std:.3f} nm</strong>, and a tail envelope from <strong>{tail_min:.3f}</strong> to <strong>{tail_max:.3f} nm</strong>.</p>
        <div class="quality-downloads">
          <a href="{html.escape(str(quality_payload['json_href']))}">Open quality JSON</a>
          <a href="{html.escape(str(quality_payload['xvg_href']))}">Open RMSD xvg</a>
          {md_anchor}
        </div>
      </div>
    </section>
    """


def discover_occupancy_payloads(base_dir: Path, system_id: str, overview_dir: Path) -> list[dict[str, object]]:
    occupancy_dir = overview_dir / "occupancy"
    occupancy_dir.mkdir(exist_ok=True)
    rows: list[dict[str, object]] = []
    for family, (subdir, prefix) in OCCUPANCY_LAYOUTS.items():
        candidates = [
            base_dir / family / system_id / subdir,
            base_dir.parent / f"{base_dir.name}_occ_{family}" / system_id / subdir,
        ]
        occ_root = next((path for path in candidates if path.exists()), None)
        if occ_root is None:
            continue
        pair_csv = occ_root / "pair_stability.csv"
        contrib_csv = occ_root / "interaction_class_contribution.csv"
        ranking_csv = occ_root / "persistent_interaction_ranking.csv"
        matrix_csv = occ_root / "residue_pair_occupancy_matrix.csv"
        peptide_mhc_csv = occ_root / "peptide_vs_mhc_summary.csv"
        region_csv = occ_root / "region_occupancy_summary.csv"
        summary_json = occ_root / "occupancy_summary.json"
        scatter_png = occ_root / f"{prefix}_occupancy_scatter.png"
        contrib_png = occ_root / f"{prefix}_interaction_class_contribution.png"
        timeline_png = occ_root / f"{prefix}_persistent_timeline.png"
        matrix_preview_png = occupancy_dir / f"{family}_occupancy_matrix_preview.png"
        if not (
            pair_csv.exists()
            and contrib_csv.exists()
            and ranking_csv.exists()
            and matrix_csv.exists()
            and peptide_mhc_csv.exists()
            and region_csv.exists()
            and summary_json.exists()
        ):
            continue

        pair_link = occupancy_dir / f"{family}_pair_stability.csv"
        contrib_link = occupancy_dir / f"{family}_interaction_class_contribution.csv"
        ranking_link = occupancy_dir / f"{family}_persistent_interaction_ranking.csv"
        matrix_link = occupancy_dir / f"{family}_residue_pair_occupancy_matrix.csv"
        peptide_mhc_link = occupancy_dir / f"{family}_peptide_vs_mhc_summary.csv"
        region_link = occupancy_dir / f"{family}_region_occupancy_summary.csv"
        summary_link = occupancy_dir / f"{family}_occupancy_summary.json"
        ensure_symlink(pair_link, pair_csv)
        ensure_symlink(contrib_link, contrib_csv)
        ensure_symlink(ranking_link, ranking_csv)
        ensure_symlink(matrix_link, matrix_csv)
        ensure_symlink(peptide_mhc_link, peptide_mhc_csv)
        ensure_symlink(region_link, region_csv)
        ensure_symlink(summary_link, summary_json)

        scatter_link = ""
        if scatter_png.exists():
            scatter_dst = occupancy_dir / scatter_png.name
            ensure_symlink(scatter_dst, scatter_png)
            scatter_link = scatter_dst.relative_to(overview_dir).as_posix()
        contrib_link_png = ""
        if contrib_png.exists():
            contrib_dst = occupancy_dir / contrib_png.name
            ensure_symlink(contrib_dst, contrib_png)
            contrib_link_png = contrib_dst.relative_to(overview_dir).as_posix()
        timeline_link_png = ""
        if timeline_png.exists():
            timeline_dst = occupancy_dir / timeline_png.name
            ensure_symlink(timeline_dst, timeline_png)
            timeline_link_png = timeline_dst.relative_to(overview_dir).as_posix()

        pair_df = pd.read_csv(pair_csv) if pair_csv.exists() else pd.DataFrame()
        contrib_df = pd.read_csv(contrib_csv) if contrib_csv.exists() else pd.DataFrame()
        ranking_df = pd.read_csv(ranking_csv) if ranking_csv.exists() else pd.DataFrame()
        matrix_df = pd.read_csv(matrix_csv, index_col=0) if matrix_csv.exists() else pd.DataFrame()
        peptide_mhc_df = pd.read_csv(peptide_mhc_csv) if peptide_mhc_csv.exists() else pd.DataFrame()
        region_df = pd.read_csv(region_csv) if region_csv.exists() else pd.DataFrame()
        summary = json.loads(summary_json.read_text(encoding="utf-8"))
        if not matrix_df.empty:
            matrix_long = (
                matrix_df.reset_index(names="tcr_residue")
                .melt(id_vars="tcr_residue", var_name="phla_residue", value_name="occupancy")
            )
            matrix_head_df = (
                matrix_long[matrix_long["occupancy"] > 0]
                .sort_values("occupancy", ascending=False)
                .head(8)
                .reset_index(drop=True)
            )
        else:
            matrix_head_df = pd.DataFrame(columns=["tcr_residue", "phla_residue", "occupancy"])
        write_occupancy_matrix_preview(matrix_df, matrix_preview_png, family)

        rows.append(
            {
                "family": family,
                "label": next((row["label"] for row in load_overview_rows(overview_dir / "interaction_overview.csv") if row["family"] == family), family),
                "pair_df": pair_df,
                "contrib_df": contrib_df,
                "ranking_df": ranking_df,
                "matrix_head_df": matrix_head_df,
                "peptide_mhc_df": peptide_mhc_df,
                "region_df": region_df,
                "summary": summary,
                "pair_href": pair_link.relative_to(overview_dir).as_posix(),
                "contrib_href": contrib_link.relative_to(overview_dir).as_posix(),
                "ranking_href": ranking_link.relative_to(overview_dir).as_posix(),
                "matrix_href": matrix_link.relative_to(overview_dir).as_posix(),
                "matrix_preview_href": matrix_preview_png.relative_to(overview_dir).as_posix(),
                "peptide_mhc_href": peptide_mhc_link.relative_to(overview_dir).as_posix(),
                "region_href": region_link.relative_to(overview_dir).as_posix(),
                "summary_href": summary_link.relative_to(overview_dir).as_posix(),
                "scatter_href": scatter_link,
                "contribution_href": contrib_link_png,
                "timeline_href": timeline_link_png,
            }
        )
    return rows


def _build_small_pair_table(frame: pd.DataFrame, mode: str) -> str:
    if frame.empty:
        return '<div class="mini-table-empty">No pair-level stability records.</div>'
    work = frame.copy()
    for col in ["contact_frequency", "max_consecutive_fraction", "n_segments"]:
        if col in work.columns:
            work[col] = pd.to_numeric(work[col], errors="coerce").fillna(0.0)
    if mode == "stable":
        top = work.sort_values(
            ["max_consecutive_fraction", "contact_frequency", "n_segments"],
            ascending=[False, False, True],
        ).head(5)
        title = "Most persistent pairs"
    else:
        transient = work[work.get("stability_label", "").astype(str) == "transient"]
        source = transient if not transient.empty else work
        top = source.sort_values(
            ["n_segments", "contact_frequency", "max_consecutive_fraction"],
            ascending=[False, True, True],
        ).head(5)
        title = "Most intermittent pairs"

    rows = []
    for _, row in top.iterrows():
        rows.append(
            "<tr>"
            f"<td>{html.escape(str(row.get('pair_label', '-')))}</td>"
            f"<td>{float(row.get('contact_frequency', 0.0)):.2f}</td>"
            f"<td>{float(row.get('max_consecutive_fraction', 0.0)):.2f}</td>"
            f"<td>{int(float(row.get('n_segments', 0) or 0))}</td>"
            "</tr>"
        )
    return (
        f'<div class="mini-table-wrap"><div class="mini-table-title">{html.escape(title)}</div>'
        '<table class="mini-table"><thead><tr><th>pair</th><th>occ.</th><th>longest</th><th>segments</th></tr></thead>'
        f"<tbody>{''.join(rows)}</tbody></table></div>"
    )


def _build_compact_summary_table(frame: pd.DataFrame, columns: list[tuple[str, str]]) -> str:
    if frame.empty:
        return '<div class="mini-table-empty">No summary records.</div>'
    top = frame.head(8).copy()
    header = "".join(f"<th>{html.escape(label)}</th>" for _, label in columns)
    rows = []
    for _, row in top.iterrows():
        cells = []
        for key, _ in columns:
            value = row.get(key, "")
            if isinstance(value, float):
                cells.append(f"<td>{value:.2f}</td>")
            else:
                cells.append(f"<td>{html.escape(str(value))}</td>")
        rows.append(f"<tr>{''.join(cells)}</tr>")
    return f'<div class="mini-table-wrap"><table class="mini-table"><thead><tr>{header}</tr></thead><tbody>{"".join(rows)}</tbody></table></div>'


def write_occupancy_matrix_preview(
    matrix_df: pd.DataFrame,
    out_path: Path,
    family: str,
) -> None:
    fig, ax = plt.subplots(figsize=(5.1, 3.6))
    if matrix_df.empty or matrix_df.shape[0] == 0 or matrix_df.shape[1] == 0:
        ax.text(0.5, 0.5, "No occupancy matrix", ha="center", va="center", fontsize=12)
        ax.axis("off")
    else:
        matrix = matrix_df.copy()
        row_order = np.argsort(-matrix.max(axis=1).to_numpy())
        col_order = np.argsort(-matrix.max(axis=0).to_numpy())
        matrix = matrix.iloc[row_order, col_order]
        trimmed = matrix.iloc[:14, :18]
        image = ax.imshow(trimmed.to_numpy(dtype=float), cmap=HEATMAP_CMAP, vmin=0.0, vmax=1.0, aspect="auto")
        ax.set_title(f"{family} occupancy matrix", fontsize=12, weight="bold")
        ax.set_xlabel("pHLA residues")
        ax.set_ylabel("TCR residues")
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        cbar = fig.colorbar(image, ax=ax, fraction=0.045, pad=0.02)
        cbar.set_label("Occupancy", fontsize=9)
        cbar.ax.tick_params(labelsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=220, facecolor="white")
    plt.close(fig)


def build_occupancy_section(occupancy_rows: list[dict[str, object]]) -> str:
    if not occupancy_rows:
        return ""
    cards = []
    for row in occupancy_rows:
        ranking_table = _build_compact_summary_table(
            row["ranking_df"],
            [
                ("pair_label", "pair"),
                ("contact_frequency", "occ."),
                ("max_consecutive_fraction", "longest"),
            ],
        )
        peptide_mhc_table = _build_compact_summary_table(
            row["peptide_mhc_df"],
            [
                ("group", "group"),
                ("n_pairs", "pairs"),
                ("mean_occupancy", "mean occ."),
                ("stable_pair_fraction", "stable frac."),
            ],
        )
        region_table = _build_compact_summary_table(
            row["region_df"],
            [
                ("tcr_region_group", "TCR side"),
                ("partner_region_group", "partner side"),
                ("mean_occupancy", "mean occ."),
            ],
        )
        matrix_table = _build_compact_summary_table(
            row["matrix_head_df"],
            [("tcr_residue", "TCR residue"), ("phla_residue", "pHLA residue"), ("occupancy", "occ.")],
        )
        cards.append(
            f"""
            <article class="occupancy-card">
              <div class="family-detail-head">
                <div>
                  <h3>{html.escape(str(row['family']))}</h3>
                  <p>{html.escape(str(row['label']))}</p>
                </div>
                <div class="occupancy-links">
                  <a href="{html.escape(str(row['ranking_href']))}">Ranking CSV</a>
                  <a href="{html.escape(str(row['matrix_href']))}">Matrix CSV</a>
                </div>
              </div>
              <div class="stability-metrics">
                <span><strong>{int(row['summary'].get('n_pairs', 0))}</strong> pairs</span>
                <span><strong>{int(row['summary'].get('n_stable_pairs', 0))}</strong> persistent</span>
                <span><strong>{int(row['summary'].get('n_transient_pairs', 0))}</strong> intermittent</span>
              </div>
              <div class="occupancy-figure-grid">
                <figure class="insight-card">
                  <img src="{html.escape(str(row['timeline_href']))}" alt="{html.escape(str(row['family']))} timeline">
                  <figcaption>
                    <strong>{html.escape(str(row['family']))} interaction timeline</strong>
                    <span>Top persistent pairs shown as compressed segment tracks across the sampled frames.</span>
                  </figcaption>
                </figure>
                <figure class="insight-card">
                  <img src="{html.escape(str(row['contribution_href']))}" alt="{html.escape(str(row['family']))} peptide vs mhc contribution">
                  <figcaption>
                    <strong>{html.escape(str(row['family']))} peptide vs MHC contribution</strong>
                    <span>Mean occupancy separated into TCR–peptide and TCR–MHC interaction classes.</span>
                  </figcaption>
                </figure>
                <figure class="insight-card">
                  <img src="{html.escape(str(row['matrix_preview_href']))}" alt="{html.escape(str(row['family']))} occupancy matrix">
                  <figcaption>
                    <strong>{html.escape(str(row['family']))} occupancy matrix preview</strong>
                    <span>A compact residue-pair occupancy matrix sorted by dominant pairs for quick scanning.</span>
                  </figcaption>
                </figure>
              </div>
              <details class="occupancy-details">
                <summary>Open detailed occupancy tables</summary>
                <div class="stability-table-grid">
                  {ranking_table}
                  {peptide_mhc_table}
                </div>
                <div class="occupancy-region-block">
                  <div class="mini-table-title">Matrix head: strongest residue pairs</div>
                  {matrix_table}
                </div>
                <div class="occupancy-region-block">
                  <div class="mini-table-title">CDR α/β vs peptide / α1-helix / α2-helix</div>
                  {region_table}
                  <div class="occupancy-download-row">
                    <a href="{html.escape(str(row['matrix_href']))}">Open matrix CSV</a>
                    <a href="{html.escape(str(row['peptide_mhc_href']))}">Open peptide-vs-MHC CSV</a>
                    <a href="{html.escape(str(row['region_href']))}">Open region occupancy CSV</a>
                    <a href="{html.escape(str(row['pair_href']))}">Open pair stability CSV</a>
                  </div>
                </div>
              </details>
            </article>
            """
        )
    return f"""
    <section class="section" id="occupancy">
      <h2>Persistence Patterns</h2>
      <p class="note">This section focuses on which residue pairs persist over time, how peptide-facing and MHC-facing contributions differ, and where long-lived interaction support concentrates across CDR and groove regions.</p>
      <div class="stability-grid">
        {''.join(cards)}
      </div>
    </section>
    """


def build_stability_section(occupancy_rows: list[dict[str, object]]) -> str:
    if not occupancy_rows:
        return ""
    cards = []
    for row in occupancy_rows:
        summary = row["summary"]
        contrib_df = row["contrib_df"]
        peptide_mean = 0.0
        hla_mean = 0.0
        if not contrib_df.empty and "interaction_class" in contrib_df.columns:
            peptide = contrib_df[contrib_df["interaction_class"] == "peptide_tcr"]
            hla = contrib_df[contrib_df["interaction_class"] == "hla_tcr"]
            if not peptide.empty:
                peptide_mean = float(peptide["mean_occupancy"].iloc[0])
            if not hla.empty:
                hla_mean = float(hla["mean_occupancy"].iloc[0])
        stable_table = _build_small_pair_table(row["pair_df"], "stable")
        intermittent_table = _build_small_pair_table(row["pair_df"], "intermittent")
        cards.append(
            f"""
            <article class="stability-card">
              <div class="family-detail-head">
                <div>
                  <h3>{html.escape(str(row['family']))}</h3>
                  <p>{html.escape(str(row['label']))}</p>
                </div>
                <a href="{html.escape(str(row['pair_href']))}">Open pair stability CSV</a>
              </div>
              <div class="stability-metrics">
                <span><strong>{int(summary.get('n_pairs', 0))}</strong> pairs</span>
                <span><strong>{int(summary.get('n_stable_pairs', 0))}</strong> stable</span>
                <span><strong>{int(summary.get('n_transient_pairs', 0))}</strong> transient</span>
                <span><strong>{peptide_mean:.2f}</strong> peptide mean occ.</span>
                <span><strong>{hla_mean:.2f}</strong> MHC mean occ.</span>
              </div>
              <div class="stability-figure-grid">
                <figure class="insight-card">
                  <img src="{html.escape(str(row['scatter_href']))}" alt="{html.escape(str(row['family']))} occupancy scatter">
                  <figcaption>
                    <strong>{html.escape(str(row['family']))} occupancy scatter</strong>
                    <span>Occupancy versus longest continuous segment fraction. Stable interactions accumulate toward the upper-right corner.</span>
                  </figcaption>
                </figure>
                <figure class="insight-card">
                  <img src="{html.escape(str(row['contribution_href']))}" alt="{html.escape(str(row['family']))} class contribution">
                  <figcaption>
                    <strong>{html.escape(str(row['family']))} class contribution</strong>
                    <span>Mean occupancy and stable-pair fraction separated into peptide-facing and MHC-facing interaction classes.</span>
                  </figcaption>
                </figure>
              </div>
              <div class="stability-table-grid">
                {stable_table}
                {intermittent_table}
              </div>
            </article>
            """
        )
    return f"""
    <section class="section" id="stability">
      <h2>Interaction Stability</h2>
      <p class="note">This section separates long-lived interactions from fragmented encounters, and compares peptide-facing and MHC-facing contributions across families.</p>
      <div class="stability-grid">
        {''.join(cards)}
      </div>
    </section>
    """


def load_report_frames() -> dict[str, pd.DataFrame]:
    frames: dict[str, pd.DataFrame] = {}
    for family, path in REPORT_PATHS.items():
        if not path.exists():
            frames[family] = pd.DataFrame()
            continue
        frame = pd.read_csv(path)
        if "contact_frequency" in frame.columns:
            frame["contact_frequency"] = pd.to_numeric(frame["contact_frequency"], errors="coerce").fillna(0.0)
        frames[family] = frame
    return frames


def write_family_pair_counts(rows: list[dict[str, str]], figures_dir: Path) -> str:
    families = [row["family"] for row in rows]
    values = [int(row["pair_rows"]) for row in rows]
    colors = [FAMILY_COLORS.get(f, "#7b8b84") for f in families]
    fig, ax = plt.subplots(figsize=(8.5, 4.8))
    bars = ax.bar(families, values, color=colors, edgecolor="#26453c", linewidth=0.8)
    ax.set_title("Interaction Family Size", fontsize=14, weight="bold")
    ax.set_ylabel("Residue pairs")
    ax.grid(axis="y", alpha=0.2, linestyle="--")
    ax.set_axisbelow(True)
    for bar, value in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, value + max(values) * 0.02 if max(values) else 0.1, str(value),
                ha="center", va="bottom", fontsize=10)
    fig.tight_layout()
    out = figures_dir / "family_pair_counts.png"
    fig.savefig(out, dpi=220, facecolor="white")
    plt.close(fig)
    return out.name


def write_interaction_class_comparison(frames: dict[str, pd.DataFrame], figures_dir: Path) -> str:
    class_order = ["peptide_tcr", "hla_tcr", "groove_tcr"]
    family_order = list(REPORT_PATHS.keys())
    data = []
    for family in family_order:
        frame = frames[family]
        counts = {label: 0 for label in class_order}
        if not frame.empty and "interaction_class" in frame.columns:
            value_counts = frame["interaction_class"].value_counts()
            for label in class_order:
                counts[label] = int(value_counts.get(label, 0))
        data.append(counts)
    fig, ax = plt.subplots(figsize=(9, 5))
    bottom = [0] * len(family_order)
    palette = {
        "peptide_tcr": "#f1b447",
        "hla_tcr": "#7ea68d",
        "groove_tcr": "#e06c54",
    }
    for label in class_order:
        values = [row[label] for row in data]
        ax.bar(family_order, values, bottom=bottom, label=label, color=palette[label], edgecolor="white")
        bottom = [b + v for b, v in zip(bottom, values)]
    ax.set_title("Interaction Class Composition", fontsize=14, weight="bold")
    ax.set_ylabel("Annotated pairs")
    ax.legend(frameon=False, ncol=3, loc="upper right")
    ax.grid(axis="y", alpha=0.2, linestyle="--")
    ax.set_axisbelow(True)
    fig.tight_layout()
    out = figures_dir / "interaction_class_composition.png"
    fig.savefig(out, dpi=220, facecolor="white")
    plt.close(fig)
    return out.name


def write_tcr_region_comparison(frames: dict[str, pd.DataFrame], figures_dir: Path) -> str:
    region_order = ["CDR1", "CDR2", "CDR3", "non_cdr"]
    family_order = [f for f in REPORT_PATHS.keys() if not frames[f].empty]
    fig, ax = plt.subplots(figsize=(9.5, 5))
    bottom = [0] * len(family_order)
    palette = {
        "CDR1": "#2a7f62",
        "CDR2": "#4ca786",
        "CDR3": "#7bc8a4",
        "non_cdr": "#c6d6cf",
    }
    for region in region_order:
        values = []
        for family in family_order:
            frame = frames[family]
            counts = frame["tcr_region"].value_counts() if "tcr_region" in frame.columns else pd.Series(dtype=int)
            values.append(int(counts.get(region, 0)))
        ax.bar(family_order, values, bottom=bottom, label=region, color=palette[region], edgecolor="white")
        bottom = [b + v for b, v in zip(bottom, values)]
    ax.set_title("TCR Region Composition", fontsize=14, weight="bold")
    ax.set_ylabel("Annotated pairs")
    ax.legend(frameon=False, ncol=4, loc="upper right")
    ax.grid(axis="y", alpha=0.2, linestyle="--")
    ax.set_axisbelow(True)
    fig.tight_layout()
    out = figures_dir / "tcr_region_composition.png"
    fig.savefig(out, dpi=220, facecolor="white")
    plt.close(fig)
    return out.name


def write_family_region_bubble(frames: dict[str, pd.DataFrame], figures_dir: Path) -> str:
    region_order = ["CDR1", "CDR2", "CDR3", "non_cdr"]
    family_order = list(REPORT_PATHS.keys())
    xs, ys, sizes, colors = [], [], [], []
    for x, family in enumerate(family_order):
        frame = frames[family]
        counts = frame["tcr_region"].value_counts() if (not frame.empty and "tcr_region" in frame.columns) else pd.Series(dtype=int)
        for y, region in enumerate(region_order):
            value = int(counts.get(region, 0))
            xs.append(x)
            ys.append(y)
            sizes.append(max(value, 1) * 22)
            colors.append(FAMILY_COLORS.get(family, "#7b8b84"))
    fig, ax = plt.subplots(figsize=(9.5, 4.8))
    ax.scatter(xs, ys, s=sizes, c=colors, alpha=0.75, edgecolors="white", linewidths=0.8)
    ax.set_xticks(range(len(family_order)))
    ax.set_xticklabels(family_order)
    ax.set_yticks(range(len(region_order)))
    ax.set_yticklabels(region_order)
    ax.set_title("Family × TCR Region Map", fontsize=14, weight="bold")
    ax.set_xlabel("Interaction family")
    ax.set_ylabel("TCR region")
    ax.grid(alpha=0.18, linestyle="--")
    ax.set_axisbelow(True)
    fig.tight_layout()
    out = figures_dir / "family_region_map.png"
    fig.savefig(out, dpi=220, facecolor="white")
    plt.close(fig)
    return out.name


def write_contact_top_pairs(contact_frame: pd.DataFrame, figures_dir: Path) -> str:
    fig, ax = plt.subplots(figsize=(9.5, 5.5))
    if contact_frame.empty:
        ax.text(0.5, 0.5, "No contact pairs detected", ha="center", va="center", fontsize=14)
        ax.axis("off")
    else:
        top = (
            contact_frame.sort_values(["contact_frequency", "contact_frames"], ascending=[False, False])
            .head(12)
            .copy()
        )
        top["pair_label"] = top["phla_residue"] + " ↔ " + top["tcr_residue"] + " (" + top["tcr_region_detailed"] + ")"
        top = top.iloc[::-1]
        ax.barh(top["pair_label"], top["contact_frequency"], color="#2d6f5c")
        ax.set_title("Top Contact Pairs", fontsize=14, weight="bold")
        ax.set_xlabel("Contact frequency")
        ax.set_xlim(0, 1.05)
        ax.grid(axis="x", alpha=0.2, linestyle="--")
        ax.set_axisbelow(True)
    fig.tight_layout()
    out = figures_dir / "contact_top_pairs.png"
    fig.savefig(out, dpi=220, facecolor="white")
    plt.close(fig)
    return out.name


def write_contact_hotspots(contact_frame: pd.DataFrame, figures_dir: Path) -> str:
    fig, axes = plt.subplots(1, 2, figsize=(11, 5.5))
    if contact_frame.empty:
        for ax in axes:
            ax.text(0.5, 0.5, "No hotspot data", ha="center", va="center", fontsize=13)
            ax.axis("off")
    else:
        tcr = (
            contact_frame.groupby(["tcr_chain", "tcr_residue"], as_index=False)["contact_frequency"]
            .sum()
            .sort_values("contact_frequency", ascending=False)
            .head(10)
        )
        tcr["label"] = tcr["tcr_chain"].str.upper() + ":" + tcr["tcr_residue"]
        tcr = tcr.iloc[::-1]
        axes[0].barh(tcr["label"], tcr["contact_frequency"], color="#4ca786")
        axes[0].set_title("Top TCR Hotspots", fontsize=13, weight="bold")
        axes[0].set_xlabel("Summed frequency")
        axes[0].grid(axis="x", alpha=0.2, linestyle="--")
        axes[0].set_axisbelow(True)

        phla = (
            contact_frame.groupby(["phla_region", "phla_residue"], as_index=False)["contact_frequency"]
            .sum()
            .sort_values("contact_frequency", ascending=False)
            .head(10)
        )
        phla["label"] = phla["phla_region"].str.replace("_", " ", regex=False) + ":" + phla["phla_residue"]
        phla = phla.iloc[::-1]
        axes[1].barh(phla["label"], phla["contact_frequency"], color="#e08c54")
        axes[1].set_title("Top pHLA Hotspots", fontsize=13, weight="bold")
        axes[1].set_xlabel("Summed frequency")
        axes[1].grid(axis="x", alpha=0.2, linestyle="--")
        axes[1].set_axisbelow(True)
    fig.tight_layout()
    out = figures_dir / "contact_hotspots.png"
    fig.savefig(out, dpi=220, facecolor="white")
    plt.close(fig)
    return out.name


def write_family_top_pairs(family: str, frame: pd.DataFrame, figures_dir: Path) -> str:
    out = figures_dir / f"{family}_top_pairs.png"
    fig, ax = plt.subplots(figsize=(9.6, 5.4))
    if frame.empty:
        ax.text(0.5, 0.5, f"No {family} residue pairs detected", ha="center", va="center", fontsize=14)
        ax.axis("off")
    else:
        work = frame.copy()
        work["pair_label"] = (
            work["phla_residue"].astype(str)
            + " ↔ "
            + work["tcr_residue"].astype(str)
            + " ["
            + work["interaction_class"].astype(str)
            + "]"
        )
        top = (
            work.sort_values(["contact_frequency", "contact_frames"], ascending=[False, False])
            .head(10)
            .iloc[::-1]
        )
        ax.barh(top["pair_label"], top["contact_frequency"], color=FAMILY_COLORS.get(family, "#7b8b84"))
        ax.set_xlim(0, 1.05)
        ax.set_xlabel("Interaction frequency")
        ax.set_title(f"{family} · Top residue pairs", fontsize=14, weight="bold")
        ax.grid(axis="x", alpha=0.2, linestyle="--")
        ax.set_axisbelow(True)
    fig.tight_layout()
    fig.savefig(out, dpi=220, facecolor="white")
    plt.close(fig)
    return out.name


def build_pair_table_html(frame: pd.DataFrame) -> str:
    if frame.empty:
        return '<div class="pair-table-empty">No residue-pair records for this family.</div>'

    work = frame.copy()
    work["contact_frequency"] = pd.to_numeric(work["contact_frequency"], errors="coerce").fillna(0.0)
    if "contact_frames" in work.columns:
        work["contact_frames"] = pd.to_numeric(work["contact_frames"], errors="coerce").fillna(0).astype(int)
    top = work.sort_values(["contact_frequency", "contact_frames"], ascending=[False, False]).head(6)

    table_rows = []
    for _, row in top.iterrows():
        interaction_class = str(row.get("interaction_class", "-")).replace("_", " ")
        tcr_region = str(row.get("tcr_region_detailed", row.get("tcr_region", "-"))).replace("_", " ")
        table_rows.append(
            "<tr>"
            f"<td>{html.escape(str(row.get('phla_residue', '-')))}</td>"
            f"<td>{html.escape(str(row.get('tcr_residue', '-')))}</td>"
            f"<td>{html.escape(interaction_class)}</td>"
            f"<td>{html.escape(tcr_region)}</td>"
            f"<td>{float(row.get('contact_frequency', 0.0)):.2f}</td>"
            "</tr>"
        )

    return (
        '<div class="pair-table-wrap">'
        '<table class="pair-table">'
        '<thead><tr><th>pHLA residue</th><th>TCR residue</th><th>interface</th><th>TCR region</th><th>freq.</th></tr></thead>'
        f"<tbody>{''.join(table_rows)}</tbody>"
        "</table>"
        "</div>"
    )


def build_analytical_figures(rows: list[dict[str, str]], figures_dir: Path) -> list[dict[str, str]]:
    figures_dir.mkdir(exist_ok=True)
    frames = load_report_frames()
    outputs = [
        {
            "title": "Top Contact Pairs",
            "caption": "The highest-frequency residue pairs from the coarse contact layer. This is the fastest way to see what actually dominates the interface.",
            "file": write_contact_top_pairs(frames["contact"], figures_dir),
        },
        {
            "title": "Contact Hotspots",
            "caption": "A compact hotspot view on both sides of the interface, much easier to read than a full matrix when presenting key residues.",
            "file": write_contact_hotspots(frames["contact"], figures_dir),
        },
    ]
    return outputs


def build_key_takeaways(
    rows: list[dict[str, str]],
    bsa_payload: dict[str, object] | None = None,
    rmsf_payload: dict[str, object] | None = None,
) -> str:
    pair_ranking = sorted(rows, key=lambda row: int(row["pair_rows"]), reverse=True)
    top_family = pair_ranking[0]["family"] if pair_ranking else "n/a"
    nonzero = [row for row in rows if int(row["pair_rows"]) > 0]
    empty = [row["family"] for row in rows if int(row["pair_rows"]) == 0]
    lines = [
        f"<li><strong>{html.escape(top_family)}</strong> is currently the largest interaction family by residue-pair count.</li>",
        f"<li><strong>{len(nonzero)}</strong> out of <strong>{len(rows)}</strong> families produced non-empty interaction tables for this system.</li>",
    ]
    if empty:
        lines.append(
            f"<li>The following family remains empty in this demo system: <strong>{html.escape(', '.join(empty))}</strong>.</li>"
        )
    if bsa_payload:
        summary = bsa_payload.get("summary", {})
        bsa_stats = summary.get("buried_surface_area", {}) or {}
        ratio_stats = summary.get("interface_ratio", {}) or {}
        mean_bsa = float(bsa_stats.get('mean', 0.0))
        mean_ratio = float(ratio_stats.get('mean', 0.0))
        lines.append(
            "<li>"
            f"Interface burial remains centered around <strong>{mean_bsa:.1f} Å²</strong> "
            f"with a mean interface ratio of <strong>{mean_ratio:.3f}</strong>, "
            "providing a compact readout of pHLA–TCR packing alongside the interaction families."
            "</li>"
        )
        rich_families = [row for row in rows if int(row["pair_rows"]) >= 3]
        if mean_bsa >= 700 and len(rich_families) >= 3:
            lines.append(
                "<li>"
                "The interface appears consistently buried while multiple interaction families remain populated, "
                "which supports reading the current contact landscape as a structured interface rather than a sparse collision pattern."
                "</li>"
            )
        elif mean_bsa < 400 or len(rich_families) <= 1:
            lines.append(
                "<li>"
                "Interface burial is comparatively limited or only a narrow subset of interaction families is populated, "
                "so individual contacts should be interpreted cautiously and checked against the structural context."
                "</li>"
            )
        else:
            lines.append(
                "<li>"
                "Interface burial and interaction-family richness are broadly aligned, but the system still benefits from reading pair-level and occupancy summaries together."
                "</li>"
            )
    if rmsf_payload:
        summary = rmsf_payload.get("summary", {})
        region_df = rmsf_payload.get("region_df", pd.DataFrame())
        mean_rmsf = float(summary.get("mean_rmsf_angstrom", 0.0))
        max_rmsf = float(summary.get("max_rmsf_angstrom", 0.0))
        lines.append(
            "<li>"
            f"Residue-level flexibility stays centered around <strong>{mean_rmsf:.2f} Å</strong> "
            f"with a maximum observed fluctuation of <strong>{max_rmsf:.2f} Å</strong>, "
            "adding a regional motion layer beyond global RMSD."
            "</li>"
        )
        cdr3_rows = region_df[region_df["region_group"].astype(str).str.contains("CDR3", na=False)]
        peptide_row = region_df[region_df["region_group"] == "peptide"]
        if not cdr3_rows.empty and not peptide_row.empty:
            cdr3_mean = float(cdr3_rows["mean_rmsf_angstrom"].mean())
            peptide_mean = float(peptide_row["mean_rmsf_angstrom"].iloc[0])
            relation = "more flexible than" if cdr3_mean > peptide_mean else "less flexible than"
            lines.append(
                "<li>"
                f"Across this system, CDR3 segments are <strong>{relation}</strong> the peptide on average "
                f"(<strong>{cdr3_mean:.2f} Å</strong> vs <strong>{peptide_mean:.2f} Å</strong>), "
                "which helps interpret whether persistent interactions come from a tight clamp or a still-mobile interface."
                "</li>"
            )
    lines.append(
        "<li>Start with the analytical figures for scale and hotspots, then move to heatmaps only for spatial detail.</li>"
    )
    return "\n".join(lines)


def build_figure_gallery(figures: list[dict[str, str]]) -> str:
    cards = []
    for figure in figures:
        cards.append(
            f"""
            <figure class="insight-card">
              <img src="{html.escape(FIGURES_DIR_NAME + '/' + figure['file'])}" alt="{html.escape(figure['title'])}">
              <figcaption>
                <strong>{html.escape(figure['title'])}</strong>
                <span>{html.escape(figure['caption'])}</span>
              </figcaption>
            </figure>
            """
        )
    return "\n".join(cards)


def build_family_pair_gallery(rows: list[dict[str, str]], figures_dir: Path) -> str:
    frames = load_report_frames()
    cards = []
    for row in rows:
        family = row["family"]
        frame = frames.get(family, pd.DataFrame())
        figure_name = write_family_top_pairs(family, frame, figures_dir)
        pair_table = build_pair_table_html(frame)
        heatmaps = [item for item in row["heatmaps"].split("; ") if item][:2]
        thumbs = "".join(
            f'<img src="{html.escape(item)}" alt="{html.escape(Path(item).name)}">' for item in heatmaps
        )
        cards.append(
            f"""
            <article class="family-detail-card">
              <div class="family-detail-head">
                <div>
                  <h3>{html.escape(family)}</h3>
                  <p>{html.escape(row['label'])}</p>
                </div>
                <a href="{html.escape(row['overview_report'])}">Open CSV</a>
              </div>
              <div class="family-detail-metrics">
                <span><strong>{html.escape(str(row['pair_rows']))}</strong> raw pairs</span>
                <span><strong>{html.escape(str(row['report_rows']))}</strong> report rows</span>
                <span><strong>{html.escape(str(row['heatmap_count']))}</strong> heatmaps</span>
              </div>
              <img src="{html.escape(FIGURES_DIR_NAME + '/' + figure_name)}" alt="{html.escape(family + ' top pairs')}">
              {pair_table}
              <div class="family-detail-thumbs">{thumbs}</div>
            </article>
            """
        )
    return "\n".join(cards)


def clean_stale_figure_files(figures_dir: Path) -> None:
    valid_names = {
        "contact_top_pairs.png",
        "contact_hotspots.png",
        "contact_top_pairs.png",
        "hbond_top_pairs.png",
        "saltbridge_top_pairs.png",
        "hydrophobic_top_pairs.png",
        "pipi_top_pairs.png",
        "cationpi_top_pairs.png",
    }
    if not figures_dir.exists():
        return
    for file_path in figures_dir.glob("*.png"):
        if file_path.name not in valid_names:
            file_path.unlink()


def build_viewer_structure(src_pdb: Path, dst_pdb: Path) -> None:
    """构建仅保留 pHLA-TCR 复合物本体的轻量 PDB。"""
    keep_chains = {"A", "B", "C", "D", "E"}
    kept = []
    with src_pdb.resolve().open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            record = line[:6].strip()
            if record in {"ATOM", "HETATM"} and len(line) > 21 and line[21] in keep_chains:
                kept.append(line.rstrip("\n"))
            elif record in {"TER", "END"}:
                kept.append(line.rstrip("\n"))
    dst_pdb.write_text("\n".join(kept) + "\n", encoding="utf-8")


def build_viewer_metadata(cdr_json: Path, contact_csv: Path, dst_json: Path) -> None:
    cdr_meta = json.loads(cdr_json.read_text(encoding="utf-8"))

    cdr_ranges: dict[str, dict[str, list[int]]] = {}
    chain_map = {
        "A": "HLA_alpha",
        "B": "beta2m",
        "C": "peptide",
        "D": "TCR_alpha",
        "E": "TCR_beta",
    }
    groove_residues: set[int] = set()
    hotspot_tcr: dict[tuple[str, int], float] = {}
    hotspot_phla: dict[tuple[str, int], float] = {}

    for chain_name, chain_info in cdr_meta.get("chains", {}).items():
        selection = chain_info.get("selection", "")
        chain_id = selection.split()[-1] if selection else ""
        cdr_ranges[chain_id] = {}
        for key in ("cdr1", "cdr2", "cdr3"):
            if key in chain_info:
                cdr_ranges[chain_id][key.upper()] = chain_info[key]["residue_range"]

    with contact_csv.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            try:
                freq = float(row.get("contact_frequency", "0") or 0)
            except ValueError:
                freq = 0.0
            if row.get("mhc_region") == "groove" and row.get("phla_chain_id") == "A":
                try:
                    groove_residues.add(int(row["phla_resid"]))
                except (KeyError, ValueError):
                    pass
            try:
                t_key = (row["tcr_chain_id"], int(row["tcr_resid"]))
                p_key = (row["phla_chain_id"], int(row["phla_resid"]))
            except (KeyError, ValueError):
                continue
            hotspot_tcr[t_key] = max(freq, hotspot_tcr.get(t_key, 0.0))
            hotspot_phla[p_key] = max(freq, hotspot_phla.get(p_key, 0.0))

    top_tcr = sorted(hotspot_tcr.items(), key=lambda x: (-x[1], x[0][0], x[0][1]))[:8]
    top_phla = sorted(hotspot_phla.items(), key=lambda x: (-x[1], x[0][0], x[0][1]))[:8]
    family_hotspots: dict[str, dict[str, list[dict[str, object]]]] = {}
    family_pairs: dict[str, list[dict[str, object]]] = {}
    for family, report_path in REPORT_PATHS.items():
        family_hotspots[family] = {"tcr": [], "phla": []}
        family_pairs[family] = []
        if not report_path.exists():
            continue
        try:
            frame = pd.read_csv(report_path)
        except Exception:
            continue
        if frame.empty or "contact_frequency" not in frame.columns:
            continue
        frame["contact_frequency"] = pd.to_numeric(frame["contact_frequency"], errors="coerce").fillna(0.0)
        if {"tcr_chain_id", "tcr_resid"}.issubset(frame.columns):
            tcr_group = (
                frame.groupby(["tcr_chain_id", "tcr_resid"], as_index=False)["contact_frequency"]
                .max()
                .sort_values("contact_frequency", ascending=False)
                .head(6)
            )
            family_hotspots[family]["tcr"] = [
                {"chain": row["tcr_chain_id"], "resi": int(row["tcr_resid"]), "score": float(row["contact_frequency"])}
                for _, row in tcr_group.iterrows()
            ]
        if {"phla_chain_id", "phla_resid"}.issubset(frame.columns):
            phla_group = (
                frame.groupby(["phla_chain_id", "phla_resid"], as_index=False)["contact_frequency"]
                .max()
                .sort_values("contact_frequency", ascending=False)
                .head(6)
            )
            family_hotspots[family]["phla"] = [
                {"chain": row["phla_chain_id"], "resi": int(row["phla_resid"]), "score": float(row["contact_frequency"])}
                for _, row in phla_group.iterrows()
            ]
        if {
            "phla_chain_id",
            "phla_resid",
            "phla_residue",
            "tcr_chain_id",
            "tcr_resid",
            "tcr_residue",
            "interaction_class",
            "contact_frequency",
        }.issubset(frame.columns):
            pair_top = (
                frame.sort_values(["contact_frequency", "contact_frames"], ascending=[False, False])
                .copy()
            )
            family_pairs[family] = [
                {
                    "pair_id": f"{family}_{idx}",
                    "phla_chain": row["phla_chain_id"],
                    "phla_resi": int(row["phla_resid"]),
                    "phla_residue": str(row["phla_residue"]),
                    "tcr_chain": row["tcr_chain_id"],
                    "tcr_resi": int(row["tcr_resid"]),
                    "tcr_residue": str(row["tcr_residue"]),
                    "interaction_class": str(row["interaction_class"]),
                    "tcr_region": str(row.get("tcr_region_detailed", row.get("tcr_region", ""))),
                    "score": float(row["contact_frequency"]),
                }
                for idx, (_, row) in enumerate(pair_top.iterrows(), start=1)
            ]

    payload = {
        "chain_map": chain_map,
        "cdr_ranges": cdr_ranges,
        "groove_residues": sorted(groove_residues),
        "top_tcr_hotspots": [
            {"chain": chain, "resi": resid, "score": score}
            for (chain, resid), score in top_tcr
        ],
        "top_phla_hotspots": [
            {"chain": chain, "resi": resid, "score": score}
            for (chain, resid), score in top_phla
        ],
        "family_hotspots": family_hotspots,
        "family_pairs": family_pairs,
    }
    dst_json.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="生成单体系 interaction HTML 报告。")
    parser.add_argument(
        "--base-dir",
        type=Path,
        default=Path("/home/xumy/work/development/Immunex/output/interaction_full_demo"),
        help="包含各 interaction family 子目录的结果根目录。",
    )
    parser.add_argument(
        "--system-id",
        type=str,
        default="1ao7_run1",
        help="需要生成报告的体系 ID。",
    )
    parser.add_argument(
        "--source-pdb",
        type=Path,
        default=None,
        help="可选：指定 viewer 使用的复合物 PDB。",
    )
    parser.add_argument(
        "--bsa-root",
        type=Path,
        default=None,
        help="可选：指定 BSA 结果根目录（包含 analysis/interface）。",
    )
    parser.add_argument(
        "--rmsf-root",
        type=Path,
        default=None,
        help="可选：指定 RMSF 结果根目录（包含 analysis/rmsf）。",
    )
    parser.add_argument(
        "--identity-root",
        type=Path,
        default=None,
        help="可选：指定 identity 结果根目录（包含 analysis/identity）。",
    )
    parser.add_argument(
        "--cluster-root",
        type=Path,
        default=None,
        help="可选：指定关键界面聚类结果根目录（包含 analysis/conformation/interface_clustering）。",
    )
    parser.add_argument(
        "--include-sections",
        type=str,
        default=None,
        help="可选：逗号分隔，仅保留指定 section（overview,quality,interface,flexibility,cluster,occupancy,contact,interactions,downloads）。",
    )
    parser.add_argument(
        "--exclude-sections",
        type=str,
        default=None,
        help="可选：逗号分隔，排除指定 section。",
    )
    return parser.parse_args()


def parse_section_list(raw_value: str | None) -> list[str] | None:
    if not raw_value:
        return None
    values = [item.strip().lower() for item in raw_value.split(",") if item.strip()]
    values = [item for item in values if item in ALL_REPORT_SECTIONS]
    return values or None


def build_interaction_html_report(
    base_dir: Path,
    system_id: str,
    source_pdb: Path | None = None,
    bsa_root: Path | None = None,
    rmsf_root: Path | None = None,
    identity_root: Path | None = None,
    cluster_root: Path | None = None,
    include_sections: list[str] | None = None,
    exclude_sections: list[str] | None = None,
) -> Path:
    overview_dir = base_dir / "overview"
    global REPORT_PATHS
    REPORT_PATHS = build_report_paths(base_dir, system_id)
    overview_csv = overview_dir / "interaction_overview.csv"
    overview_json = overview_dir / "interaction_overview.json"
    rows = load_overview_rows(overview_csv)
    meta = json.loads(overview_json.read_text(encoding="utf-8"))
    system_id = meta["system_id"]
    viewer_pdb = overview_dir / VIEWER_PDB_NAME
    viewer_meta = overview_dir / VIEWER_META_NAME
    figures_dir = overview_dir / FIGURES_DIR_NAME
    clean_stale_figure_files(figures_dir)
    source_pdb = choose_source_pdb(base_dir, system_id, source_pdb)
    quality_root = choose_quality_root(None, source_pdb)
    bsa_root = choose_bsa_root(bsa_root, source_pdb, system_id)
    rmsf_root = choose_rmsf_root(rmsf_root, source_pdb, system_id)
    cluster_root = choose_cluster_root(cluster_root, system_id)
    occupancy_payload = discover_occupancy_payloads(base_dir, system_id, overview_dir)
    contact_report = base_dir / "contact" / system_id / "analysis/contacts/contact_report.csv"
    cdr_metadata = base_dir / "contact" / system_id / "cdr_metadata.json"
    build_viewer_structure(source_pdb, viewer_pdb)
    build_viewer_metadata(cdr_metadata, contact_report, viewer_meta)
    analytical_figures = build_analytical_figures(rows, figures_dir)
    quality_payload = build_quality_assets(overview_dir, figures_dir, quality_root)
    bsa_payload = build_bsa_assets(overview_dir, bsa_root)
    rmsf_payload = build_rmsf_assets(overview_dir, rmsf_root)
    identity_payload = build_biological_identity_assets(overview_dir, source_pdb, cdr_metadata, identity_root=identity_root)
    cluster_payload = build_cluster_assets(overview_dir, cluster_root)
    visible_sections = normalize_visible_sections(
        include_sections=include_sections,
        exclude_sections=exclude_sections,
        quality_payload=quality_payload,
        bsa_payload=bsa_payload,
        rmsf_payload=rmsf_payload,
        cluster_payload=cluster_payload,
        occupancy_payload=occupancy_payload,
    )

    quality_section = optional_block("quality" in visible_sections, build_quality_section(quality_payload))
    interface_section = optional_block("interface" in visible_sections, build_bsa_section(bsa_payload))
    flexibility_section = optional_block("flexibility" in visible_sections, build_rmsf_section(rmsf_payload))
    cluster_section = optional_block("cluster" in visible_sections, build_cluster_section(cluster_payload))
    occupancy_section = optional_block("occupancy" in visible_sections, build_occupancy_section(occupancy_payload))
    contact_section = optional_block(
        "contact" in visible_sections,
        f"""
    <section class="section" id="contact">
      <h2>Interface Landscape</h2>
      <p class="note">This section is the global map of the interface: which families are populated, how large they are, and where the main hotspots sit before drilling down into specific residue pairs.</p>
      <div class="family-grid">
        {build_cards(rows)}
      </div>
    </section>

    <section class="section">
      <h2>Unified Summary Table</h2>
      <p class="note">This table is kept compact on purpose. It gives one comparison layer across families before moving into persistence or pair-level interpretation.</p>
      {build_table(rows)}
    </section>

    <section class="section">
      <h2>Landscape Takeaways</h2>
      <p class="note">These are the shortest global conclusions that can be read before opening any family-specific detail.</p>
      <ul class="takeaways">
        {build_key_takeaways(rows, bsa_payload, rmsf_payload)}
      </ul>
    </section>

    <section class="section">
      <h2>Analytical Figures</h2>
      <p class="note">These figures summarize scale and hotspot structure at the landscape level. They are meant to orient the reader, not to replace pair-level interpretation.</p>
      <div class="insight-grid">
        {build_figure_gallery(analytical_figures)}
      </div>
    </section>

    <section class="section">
      <h2>Representative Heatmaps</h2>
      <p class="note">Only the most informative spatial views are shown here, so the landscape section stays global rather than turning into a full result browser.</p>
      {build_heatmap_sections(rows)}
    </section>
        """,
    )
    interactions_section = optional_block(
        "interactions" in visible_sections,
        f"""
    <section class="section" id="interactions">
      <h2>Residue-Pair Definitions</h2>
      <p class="note">This section shifts from family-level landscape to concrete residue pairs. Each family is shown through the amino-acid pairs that actually define it.</p>
      <div class="family-detail-grid">
        {build_family_pair_gallery(rows, figures_dir)}
      </div>
    </section>
        """,
    )
    downloads_section = optional_block(
        "downloads" in visible_sections,
        build_downloads_section(identity_payload, bsa_payload, rmsf_payload, cluster_payload),
    )

    html_text = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>Interaction Demo Report - {html.escape(system_id)}</title>
  <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
  <style>
    :root {{
      --bg: #f5f1e8;
      --paper: #fffdfa;
      --ink: #1e2a24;
      --muted: #6f746d;
      --line: #ddd5c7;
      --accent: #255b4b;
      --accent-soft: #dbe9e1;
      --warn: #a36b00;
      --empty: #8a8f88;
      --shadow: 0 12px 30px rgba(29, 37, 32, 0.08);
    }}
    * {{ box-sizing: border-box; }}
    body {{
      margin: 0;
      font-family: "Noto Sans", "PingFang SC", "Microsoft YaHei", sans-serif;
      color: var(--ink);
      background:
        radial-gradient(circle at top left, #f7efe0 0, transparent 30%),
        linear-gradient(180deg, #efe7d7 0%, var(--bg) 24%, #f8f5ef 100%);
    }}
    .page {{
      max-width: 1440px;
      margin: 0 auto;
      padding: 40px 28px 80px;
    }}
    .page-nav {{
      display: flex;
      flex-wrap: wrap;
      gap: 10px;
      margin-bottom: 14px;
      position: sticky;
      top: 14px;
      z-index: 8;
      padding: 10px 12px;
      border-radius: 18px;
      background: rgba(255, 253, 250, 0.88);
      border: 1px solid rgba(196, 184, 164, 0.75);
      box-shadow: 0 10px 24px rgba(29, 37, 32, 0.08);
      backdrop-filter: blur(10px);
    }}
    .page-nav a {{
      text-decoration: none;
      color: #174538;
      font-size: 13px;
      font-weight: 700;
      padding: 10px 14px;
      border-radius: 999px;
      background: linear-gradient(180deg, #ffffff 0%, #f3f7f4 100%);
      border: 1px solid #c7d7cf;
      backdrop-filter: blur(6px);
      box-shadow: 0 2px 0 rgba(255, 255, 255, 0.75) inset;
      transition: transform 120ms ease, box-shadow 120ms ease, background 120ms ease, color 120ms ease;
    }}
    .page-nav a:hover {{
      color: #0f2e25;
      background: linear-gradient(180deg, #ffffff 0%, #e7f0eb 100%);
      box-shadow: 0 6px 14px rgba(37, 91, 75, 0.14);
      transform: translateY(-1px);
    }}
    .page-nav a.active {{
      color: #f7f5ef;
      background: linear-gradient(180deg, #2f5d50 0%, #21453a 100%);
      border-color: #21453a;
      box-shadow: 0 10px 20px rgba(33, 69, 58, 0.22);
    }}
    .page-nav a:focus-visible {{
      outline: 2px solid #2f5d50;
      outline-offset: 2px;
    }}
    .hero {{
      padding: 22px 24px 24px;
      background: linear-gradient(135deg, #1c473c 0%, #2d6f5c 100%);
      color: #f7f5ef;
      border-radius: 24px;
      box-shadow: var(--shadow);
    }}
    .hero h1 {{
      margin: 0 0 8px;
      font-size: 42px;
      line-height: 1.02;
      letter-spacing: -0.02em;
    }}
    .hero p {{
      margin: 4px 0;
      color: rgba(247, 245, 239, 0.88);
      font-size: 15px;
    }}
    .hero-layout {{
      display: grid;
      grid-template-columns: minmax(280px, 0.72fr) minmax(460px, 1.28fr);
      gap: 18px;
      align-items: stretch;
    }}
    .hero-copy {{
      display: flex;
      flex-direction: column;
      justify-content: space-between;
      min-height: 100%;
      padding: 8px 4px;
    }}
    .hero-kicker {{
      display: inline-flex;
      align-items: center;
      gap: 8px;
      margin-bottom: 18px;
      padding: 7px 12px;
      border-radius: 999px;
      width: fit-content;
      background: rgba(255, 255, 255, 0.12);
      border: 1px solid rgba(255, 255, 255, 0.18);
      color: rgba(247, 245, 239, 0.92);
      font-size: 12px;
      font-weight: 700;
      letter-spacing: 0.08em;
      text-transform: uppercase;
    }}
    .hero-lead {{
      max-width: 28ch;
      margin-top: 12px;
      font-size: 17px;
      line-height: 1.5;
      color: rgba(247, 245, 239, 0.92);
    }}
    .hero-group {{
      margin-top: 22px;
    }}
    .hero-group-title {{
      margin: 0 0 10px;
      color: rgba(247, 245, 239, 0.72);
      font-size: 12px;
      font-weight: 800;
      letter-spacing: 0.08em;
      text-transform: uppercase;
    }}
    .hero-identity {{
      margin-top: 0;
      display: grid;
      grid-template-columns: repeat(2, minmax(0, 1fr));
      gap: 10px;
    }}
    .hero-identity-card {{
      margin: 0;
      padding: 14px 15px 13px;
      border-radius: 18px;
      background: linear-gradient(180deg, rgba(255, 255, 255, 0.13), rgba(255, 255, 255, 0.08));
      border: 1px solid rgba(255, 255, 255, 0.14);
      box-shadow: inset 0 1px 0 rgba(255, 255, 255, 0.06);
      min-height: 92px;
      position: relative;
      overflow: hidden;
    }}
    .hero-identity-card::before {{
      content: "";
      position: absolute;
      left: 0;
      top: 0;
      width: 100%;
      height: 3px;
      background: linear-gradient(90deg, rgba(255, 209, 102, 0.72), rgba(126, 194, 173, 0.42));
    }}
    .hero-identity-card span {{
      display: block;
      color: rgba(247, 245, 239, 0.66);
      font-size: 11px;
      font-weight: 700;
      letter-spacing: 0.04em;
      text-transform: uppercase;
      margin-bottom: 8px;
    }}
    .hero-identity-card strong {{
      display: block;
      color: #fffaf2;
      font-size: 17px;
      line-height: 1.18;
      word-break: break-word;
    }}
    .hero-identity-card p {{
      margin: 8px 0 0;
      font-size: 13px;
      color: rgba(247, 245, 239, 0.8);
    }}
    .hero-meta {{
      margin-top: 0;
      display: grid;
      grid-template-columns: repeat(2, minmax(0, 1fr));
      gap: 10px;
    }}
    .hero-stat {{
      padding: 11px 14px 12px;
      border-radius: 15px;
      background: rgba(18, 42, 35, 0.34);
      border: 1px solid rgba(255, 255, 255, 0.08);
      box-shadow: inset 0 1px 0 rgba(255, 255, 255, 0.03);
    }}
    .hero-stat span {{
      display: block;
      color: rgba(247, 245, 239, 0.58);
      font-size: 11px;
      font-weight: 700;
      letter-spacing: 0.05em;
      text-transform: uppercase;
      margin-bottom: 8px;
    }}
    .hero-stat strong {{
      display: block;
      color: #fffdf8;
      font-size: 24px;
      font-weight: 800;
      line-height: 1.1;
    }}
    .hero-figure {{
      background: rgba(255, 255, 255, 0.08);
      border: 1px solid rgba(255, 255, 255, 0.12);
      border-radius: 22px;
      padding: 16px;
      display: flex;
      flex-direction: column;
    }}
    .hero-figure svg {{
      display: block;
      width: 100%;
      height: auto;
    }}
    .viewer-shell {{
      border-radius: 20px;
      overflow: hidden;
      background: rgba(15, 29, 24, 0.28);
      border: 1px solid rgba(255, 255, 255, 0.16);
      position: relative;
    }}
    .viewer-stage {{
      width: 100%;
      height: 390px;
      position: relative;
      overflow: hidden;
      isolation: isolate;
    }}
    .viewer-toolbar {{
      display: flex;
      flex-wrap: wrap;
      gap: 8px;
      padding: 12px;
      background: rgba(255, 255, 255, 0.08);
      border-top: 1px solid rgba(255, 255, 255, 0.08);
    }}
    .viewer-toolbar button {{
      appearance: none;
      border: 0;
      border-radius: 999px;
      padding: 8px 12px;
      background: rgba(255, 255, 255, 0.14);
      color: #f7f5ef;
      font-size: 13px;
      cursor: pointer;
    }}
    .viewer-toolbar button:hover {{
      background: rgba(255, 255, 255, 0.2);
    }}
    .viewer-note {{
      margin-top: 10px;
      font-size: 13px;
      color: rgba(247, 245, 239, 0.82);
    }}
    .viewer-pair-panel {{
      margin-top: 12px;
      border-radius: 16px;
      border: 1px solid rgba(255, 255, 255, 0.12);
      background: rgba(255, 255, 255, 0.08);
      overflow: hidden;
    }}
    .viewer-pair-panel-head {{
      display: flex;
      align-items: center;
      justify-content: space-between;
      gap: 12px;
      padding: 12px 14px 10px;
      border-bottom: 1px solid rgba(255, 255, 255, 0.08);
    }}
    .viewer-pair-panel-head strong {{
      font-size: 13px;
      color: #f7f5ef;
    }}
    .viewer-pair-panel-head span {{
      font-size: 12px;
      color: rgba(247, 245, 239, 0.78);
    }}
    .viewer-pair-actions {{
      display: flex;
      flex-wrap: wrap;
      gap: 8px;
      padding: 10px 14px 0;
    }}
    .viewer-pair-actions button {{
      appearance: none;
      border: 0;
      border-radius: 999px;
      padding: 7px 11px;
      background: rgba(255, 255, 255, 0.14);
      color: #f7f5ef;
      font-size: 12px;
      cursor: pointer;
    }}
    .viewer-pair-list {{
      max-height: 208px;
      overflow-y: auto;
      padding: 10px 14px 14px;
      display: grid;
      gap: 8px;
    }}
    .viewer-pair-row {{
      display: grid;
      grid-template-columns: auto 1fr auto;
      gap: 10px;
      align-items: start;
      padding: 10px 12px;
      border-radius: 12px;
      background: rgba(255, 255, 255, 0.08);
    }}
    .viewer-pair-row input {{
      margin-top: 2px;
    }}
    .viewer-pair-main {{
      display: grid;
      gap: 3px;
    }}
    .viewer-pair-main strong {{
      font-size: 13px;
      color: #fffaf2;
      font-weight: 700;
    }}
    .viewer-pair-main span {{
      font-size: 12px;
      color: rgba(247, 245, 239, 0.76);
    }}
    .viewer-pair-score {{
      font-size: 12px;
      color: #fffaf2;
      font-variant-numeric: tabular-nums;
      padding-top: 1px;
    }}
    .viewer-pair-empty {{
      padding: 14px;
      color: rgba(247, 245, 239, 0.78);
      font-size: 13px;
    }}
    .viewer-subtoolbar {{
      display: flex;
      flex-wrap: wrap;
      gap: 8px;
      margin-top: 10px;
    }}
    .viewer-subtoolbar button {{
      appearance: none;
      border: 0;
      border-radius: 999px;
      padding: 7px 11px;
      background: rgba(255, 255, 255, 0.16);
      color: #f7f5ef;
      font-size: 12px;
      cursor: pointer;
    }}
    .section {{
      margin-top: 28px;
      background: var(--paper);
      border: 1px solid var(--line);
      border-radius: 22px;
      padding: 24px;
      box-shadow: var(--shadow);
    }}
    .quality-grid {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
      gap: 14px;
      margin-bottom: 16px;
    }}
    .quality-banner {{
      display: grid;
      grid-template-columns: auto 1fr;
      gap: 14px;
      align-items: center;
      margin-bottom: 16px;
      padding: 14px 16px;
      border-radius: 18px;
      border: 1px solid #cfe0d6;
      background: linear-gradient(135deg, #f8fcfa 0%, #eef5f0 100%);
    }}
    .quality-banner-badge {{
      padding: 8px 12px;
      border-radius: 999px;
      background: #2f5d50;
      color: #fffaf2;
      font-size: 12px;
      font-weight: 700;
      letter-spacing: 0.03em;
      white-space: nowrap;
    }}
    .quality-banner-body {{
      display: grid;
      gap: 4px;
    }}
    .quality-banner-body strong {{
      font-size: 15px;
      color: var(--ink);
    }}
    .quality-banner-body span {{
      color: var(--muted);
      font-size: 14px;
      line-height: 1.5;
    }}
    .quality-card {{
      padding: 16px 18px;
      border: 1px solid var(--line);
      border-radius: 18px;
      background: linear-gradient(180deg, #fffefb 0%, #f5faf7 100%);
    }}
    .quality-card span {{
      display: block;
      color: var(--muted);
      font-size: 13px;
      margin-bottom: 8px;
    }}
    .quality-card strong {{
      display: block;
      font-size: 28px;
      color: var(--ink);
      line-height: 1;
      font-variant-numeric: tabular-nums;
    }}
    .quality-note-box {{
      margin-top: 16px;
      padding: 16px 18px;
      border: 1px solid var(--line);
      border-radius: 18px;
      background: #fffef9;
      color: var(--ink);
    }}
    .quality-figures .insight-card {{
      border-color: #cfe0d6;
      background: linear-gradient(180deg, #ffffff 0%, #fbfdfc 100%);
    }}
    .quality-note-box p {{
      margin: 0 0 10px;
      line-height: 1.65;
    }}
    .quality-downloads {{
      display: flex;
      flex-wrap: wrap;
      gap: 10px;
      margin-top: 12px;
    }}
    .quality-downloads a {{
      text-decoration: none;
      color: var(--accent);
      font-weight: 700;
      font-size: 13px;
    }}
    .section h2 {{
      margin: 0 0 16px;
      font-size: 24px;
    }}
    .takeaways {{
      margin: 0;
      padding-left: 20px;
      display: grid;
      gap: 10px;
      color: var(--ink);
    }}
    .takeaways li {{
      line-height: 1.55;
    }}
    .section p.note {{
      margin-top: -4px;
      color: var(--muted);
      font-size: 14px;
    }}
    .family-grid {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(240px, 1fr));
      gap: 16px;
    }}
    .family-card {{
      border: 1px solid var(--line);
      border-radius: 18px;
      padding: 18px;
      background: #fffefb;
    }}
    .family-card.active {{
      border-color: #bfd4ca;
      background: linear-gradient(180deg, #ffffff 0%, #f4faf7 100%);
    }}
    .family-card.empty {{
      opacity: 0.78;
    }}
    .family-head {{
      display: flex;
      justify-content: space-between;
      gap: 12px;
      align-items: flex-start;
    }}
    .family-name {{
      font-size: 20px;
      font-weight: 700;
      letter-spacing: 0.02em;
    }}
    .family-label {{
      color: var(--muted);
      margin-top: 4px;
      font-size: 14px;
    }}
    .family-badge {{
      padding: 6px 10px;
      border-radius: 999px;
      background: var(--accent-soft);
      color: var(--accent);
      font-size: 13px;
      font-weight: 700;
      white-space: nowrap;
    }}
    .family-metrics {{
      display: grid;
      grid-template-columns: repeat(3, 1fr);
      gap: 10px;
      margin-top: 18px;
    }}
    .family-metrics div {{
      padding: 10px 12px;
      border-radius: 14px;
      background: #f7f3eb;
    }}
    .family-metrics span {{
      display: block;
      color: var(--muted);
      font-size: 12px;
      margin-bottom: 4px;
    }}
    .family-metrics strong {{
      font-size: 18px;
    }}
    .family-links {{
      margin-top: 14px;
    }}
    .family-links a,
    table a {{
      color: var(--accent);
      font-weight: 700;
      text-decoration: none;
    }}
    table {{
      width: 100%;
      border-collapse: collapse;
      font-size: 14px;
    }}
    th, td {{
      text-align: left;
      padding: 12px 10px;
      border-bottom: 1px solid var(--line);
    }}
    th {{
      color: var(--muted);
      font-weight: 700;
      font-size: 13px;
      text-transform: uppercase;
      letter-spacing: 0.05em;
    }}
    .heatmap-section + .heatmap-section {{
      margin-top: 28px;
      padding-top: 22px;
      border-top: 1px dashed var(--line);
    }}
    .insight-grid {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(320px, 1fr));
      gap: 18px;
    }}
    .insight-card {{
      margin: 0;
      border: 1px solid var(--line);
      border-radius: 18px;
      overflow: hidden;
      background: #fff;
    }}
    .insight-card img {{
      display: block;
      width: 100%;
      height: auto;
      background: #fff;
    }}
    .insight-card figcaption {{
      display: grid;
      gap: 6px;
      padding: 12px 14px 14px;
      font-size: 13px;
      color: var(--muted);
      border-top: 1px solid var(--line);
    }}
    .insight-card strong {{
      color: var(--ink);
      font-size: 15px;
    }}
    .section-title {{
      margin-bottom: 14px;
    }}
    .family-detail-grid {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(360px, 1fr));
      gap: 18px;
    }}
    .family-detail-card {{
      border: 1px solid var(--line);
      border-radius: 18px;
      overflow: hidden;
      background: #fff;
    }}
    .family-detail-head {{
      display: flex;
      justify-content: space-between;
      gap: 16px;
      align-items: flex-start;
      padding: 14px 16px 10px;
      border-bottom: 1px solid var(--line);
      background: #fffef9;
    }}
    .family-detail-head h3 {{
      margin: 0;
      font-size: 18px;
      text-transform: none;
    }}
    .family-detail-head p {{
      margin: 4px 0 0;
      color: var(--muted);
      font-size: 13px;
    }}
    .family-detail-head a {{
      color: var(--accent);
      font-weight: 700;
      text-decoration: none;
      white-space: nowrap;
      padding-top: 2px;
    }}
    .family-detail-card img {{
      display: block;
      width: 100%;
      height: auto;
      background: #fff;
    }}
    .pair-table-wrap {{
      margin: 14px 14px 0;
      border: 1px solid var(--line);
      border-radius: 14px;
      overflow: hidden;
      background: #fffefb;
    }}
    .pair-table {{
      width: 100%;
      border-collapse: collapse;
      font-size: 13px;
    }}
    .pair-table th,
    .pair-table td {{
      padding: 10px 12px;
      text-align: left;
      border-bottom: 1px solid #eee6d8;
    }}
    .pair-table th {{
      font-size: 12px;
      text-transform: uppercase;
      letter-spacing: 0.04em;
      color: var(--muted);
      background: #f9f4ea;
    }}
    .pair-table tbody tr:last-child td {{
      border-bottom: 0;
    }}
    .pair-table-empty {{
      margin: 14px 14px 0;
      padding: 14px 16px;
      border: 1px dashed var(--line);
      border-radius: 14px;
      color: var(--muted);
      font-size: 14px;
      background: #fffefb;
    }}
    .family-detail-metrics {{
      display: flex;
      flex-wrap: wrap;
      gap: 10px;
      padding: 10px 16px 0;
      color: var(--muted);
      font-size: 13px;
    }}
    .family-detail-metrics strong {{
      color: var(--ink);
    }}
    .family-detail-thumbs {{
      display: grid;
      grid-template-columns: repeat(2, minmax(0, 1fr));
      gap: 10px;
      padding: 10px 14px 14px;
      background: #fffef9;
      border-top: 1px solid var(--line);
    }}
    .stability-grid {{
      display: grid;
      gap: 18px;
    }}
    .stability-card {{
      border: 1px solid var(--line);
      border-radius: 18px;
      overflow: hidden;
      background: #fff;
    }}
    .stability-metrics {{
      display: flex;
      flex-wrap: wrap;
      gap: 10px;
      padding: 10px 16px 0;
      color: var(--muted);
      font-size: 13px;
    }}
    .stability-metrics strong {{
      color: var(--ink);
    }}
    .stability-figure-grid {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(280px, 1fr));
      gap: 14px;
      padding: 14px;
    }}
    .occupancy-figure-grid {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(280px, 1fr));
      gap: 14px;
      padding: 14px;
    }}
    .stability-table-grid {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(300px, 1fr));
      gap: 14px;
      padding: 0 14px 14px;
    }}
    .mini-table-wrap {{
      border: 1px solid var(--line);
      border-radius: 14px;
      overflow: hidden;
      background: #fffefb;
    }}
    .mini-table-title {{
      padding: 10px 12px;
      border-bottom: 1px solid var(--line);
      background: #f9f4ea;
      font-size: 13px;
      font-weight: 700;
      color: var(--ink);
    }}
    .mini-table {{
      width: 100%;
      border-collapse: collapse;
      font-size: 12px;
    }}
    .mini-table th,
    .mini-table td {{
      padding: 9px 10px;
      border-bottom: 1px solid #eee6d8;
      text-align: left;
      vertical-align: top;
    }}
    .mini-table th {{
      color: var(--muted);
      text-transform: uppercase;
      letter-spacing: 0.04em;
      font-size: 11px;
    }}
    .mini-table tbody tr:last-child td {{
      border-bottom: 0;
    }}
    .mini-table-empty {{
      padding: 14px;
      color: var(--muted);
      font-size: 13px;
      border: 1px dashed var(--line);
      border-radius: 12px;
      background: #fffefb;
    }}
    .occupancy-card {{
      border: 1px solid var(--line);
      border-radius: 18px;
      overflow: hidden;
      background: #fff;
    }}
    .occupancy-card .stability-metrics {{
      padding-bottom: 4px;
    }}
    .occupancy-links {{
      display: flex;
      flex-wrap: wrap;
      gap: 8px;
    }}
    .occupancy-links a,
    .occupancy-download-row a {{
      color: var(--accent);
      text-decoration: none;
      font-weight: 700;
      font-size: 13px;
    }}
    .occupancy-region-block {{
      padding: 0 14px 14px;
    }}
    .occupancy-details {{
      margin: 0 14px 14px;
      border: 1px solid var(--line);
      border-radius: 14px;
      background: #fffefb;
      overflow: hidden;
    }}
    .occupancy-details summary {{
      cursor: pointer;
      list-style: none;
      padding: 12px 14px;
      font-size: 13px;
      font-weight: 700;
      color: var(--accent);
      background: #f9f4ea;
      border-bottom: 1px solid transparent;
    }}
    .occupancy-details[open] summary {{
      border-bottom-color: var(--line);
    }}
    .occupancy-details summary::-webkit-details-marker {{
      display: none;
    }}
    .occupancy-details summary::after {{
      content: 'Expand';
      float: right;
      color: var(--muted);
      font-weight: 600;
    }}
    .occupancy-details[open] summary::after {{
      content: 'Collapse';
    }}
    .occupancy-download-row {{
      display: flex;
      flex-wrap: wrap;
      gap: 10px;
      margin-top: 10px;
    }}
    .family-detail-thumbs img {{
      display: block;
      width: 100%;
      height: auto;
      border-radius: 10px;
      border: 1px solid var(--line);
      background: #fff;
    }}
    .section-title h3 {{
      margin: 0;
      font-size: 20px;
    }}
    .section-title p {{
      margin: 4px 0 0;
      color: var(--muted);
      font-size: 14px;
    }}
    .heatmap-grid {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(260px, 1fr));
      gap: 16px;
    }}
    .heatmap-card {{
      margin: 0;
      border: 1px solid var(--line);
      border-radius: 18px;
      overflow: hidden;
      background: #fff;
    }}
    .heatmap-card img {{
      display: block;
      width: 100%;
      height: auto;
      background: #fff;
    }}
    .heatmap-card figcaption {{
      padding: 10px 12px;
      font-size: 13px;
      color: var(--muted);
      border-top: 1px solid var(--line);
    }}
    .download-directory {{
      display: grid;
      gap: 14px;
      margin-top: 18px;
    }}
    .download-group {{
      border: 1px solid var(--line);
      border-radius: 18px;
      background: linear-gradient(180deg, #fffefb 0%, #faf7f1 100%);
      overflow: hidden;
    }}
    .download-group summary {{
      list-style: none;
      cursor: pointer;
      display: grid;
      gap: 4px;
      padding: 16px 18px;
      border-bottom: 1px solid rgba(215, 205, 188, 0.75);
    }}
    .download-group summary::-webkit-details-marker {{
      display: none;
    }}
    .download-group summary span {{
      font-size: 15px;
      font-weight: 800;
      color: var(--ink);
    }}
    .download-group summary small {{
      color: var(--muted);
      font-size: 13px;
      line-height: 1.5;
    }}
    .download-group-grid {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(280px, 1fr));
      gap: 12px;
      padding: 16px 18px 18px;
    }}
    .download-group-grid a {{
      display: block;
      padding: 14px 16px;
      border: 1px solid var(--line);
      border-radius: 14px;
      color: var(--ink);
      text-decoration: none;
      background: #fffef9;
      transition: transform 120ms ease, box-shadow 120ms ease, border-color 120ms ease;
    }}
    .download-group-grid a:hover {{
      transform: translateY(-1px);
      border-color: #bfd4ca;
      box-shadow: 0 8px 20px rgba(37, 91, 75, 0.08);
    }}
    .download-group-grid small {{
      display: block;
      margin-top: 6px;
      color: var(--muted);
      line-height: 1.5;
    }}
    @media (max-width: 900px) {{
      .page {{ padding: 20px 14px 40px; }}
      .hero h1 {{ font-size: 28px; }}
      .hero-layout {{ grid-template-columns: 1fr; }}
      .hero-copy {{ padding: 0; }}
      .hero-lead, .hero-subcopy {{ max-width: none; }}
      .hero-identity {{ grid-template-columns: 1fr; }}
      .hero-meta {{ grid-template-columns: 1fr 1fr; }}
      .family-metrics {{ grid-template-columns: 1fr; }}
    }}
  </style>
</head>
<body>
  <main class="page">
    <nav class="page-nav">
      {build_nav(visible_sections)}
    </nav>
    <section class="hero" id="overview">
      <div class="hero-layout">
        <div class="hero-copy">
          <div>
            <div class="hero-kicker">pHLA-TCR trajectory briefing</div>
            <h1>Interaction Demo Report</h1>
            <p>System: <strong>{html.escape(system_id)}</strong></p>
            <p class="hero-lead">Trajectory quality, interface stability, flexibility, and interaction readouts in one report.</p>
            <div class="hero-group">
              <div class="hero-group-title">Biological Identity</div>
              <div class="hero-identity">
                {build_hero_identity(identity_payload)}
              </div>
            </div>
          </div>
          <div class="hero-group">
            <div class="hero-group-title">Quick Metrics</div>
            <div class="hero-meta">
              {build_hero_metric_strip(rows, quality_payload, bsa_payload)}
            </div>
          </div>
        </div>
        <div class="hero-figure">
          <div class="viewer-shell">
            <div id="mol-viewer" class="viewer-stage"></div>
            <div class="viewer-toolbar">
              <button type="button" data-style="chain">Color by chain</button>
              <button type="button" data-style="cdr">Color by CDR</button>
              <button type="button" data-style="groove">Highlight groove</button>
              <button type="button" data-style="hotspot">Highlight hotspots</button>
              <button type="button" data-style="sticks">Sticks</button>
              <button type="button" data-style="surface">Surface</button>
              <button type="button" data-action="reset">Reset view</button>
            </div>
            <div class="viewer-subtoolbar">
              <button type="button" data-family="contact">Show contact</button>
              <button type="button" data-family="hbond">Show hbond</button>
              <button type="button" data-family="saltbridge">Show saltbridge</button>
              <button type="button" data-family="hydrophobic">Show hydrophobic</button>
              <button type="button" data-family="pipi">Show pipi</button>
              <button type="button" data-family="cationpi">Show cationpi</button>
            </div>
            <div class="viewer-pair-panel">
              <div class="viewer-pair-panel-head">
                <div>
                  <strong id="pair-panel-title">Residue-pair browser</strong>
                  <span id="pair-panel-subtitle">Choose an interaction family, then select residue pairs to display.</span>
                </div>
              </div>
              <div class="viewer-pair-actions">
                <button type="button" id="pair-select-top">Select top 5</button>
                <button type="button" id="pair-clear-all">Clear</button>
                <button type="button" id="pair-apply">Show selected pairs</button>
              </div>
              <div id="viewer-pair-list" class="viewer-pair-list">
                <div class="viewer-pair-empty">No interaction family selected yet.</div>
              </div>
            </div>
          </div>
          <div class="viewer-note">
            Interactive pHLA-TCR complex view. Drag to rotate, scroll to zoom, then switch between interaction-family highlights.
          </div>
        </div>
      </div>
    </section>

    {quality_section}

    {interface_section}

    {flexibility_section}

    {cluster_section}

    {occupancy_section}

    {contact_section}

    {interactions_section}

    {downloads_section}
  </main>
</body>
<script>
  const viewerElement = document.getElementById('mol-viewer');
  const pdbUrl = './{VIEWER_PDB_NAME}';
  const metaUrl = './{VIEWER_META_NAME}';
  let viewer = null;
  let currentSurface = null;
  let viewerMeta = null;
  let framePending = false;
  let activeLabels = [];
  let activeShapes = [];
  let activeFamily = null;

  function frameModel() {{
    if (!viewer) return;
    viewer.resize();
    viewer.center();
    viewer.zoomTo();
    viewer.render();
  }}

  function scheduleFrame() {{
    if (!viewer || framePending) return;
    framePending = true;
    requestAnimationFrame(() => {{
      framePending = false;
      frameModel();
    }});
  }}

  function clearStyles() {{
    viewer.setStyle({{}}, {{}});
  }}

  function applyChainStyle() {{
    viewer.setStyle({{}}, {{ cartoon: {{ color: '#c9d5cf' }} }});
    viewer.setStyle({{ chain: 'A' }}, {{ cartoon: {{ color: '#c8a77e' }} }});
    viewer.setStyle({{ chain: 'B' }}, {{ cartoon: {{ color: '#9b7f63' }} }});
    viewer.setStyle({{ chain: 'C' }}, {{ cartoon: {{ color: '#f1b447' }} }});
    viewer.setStyle({{ chain: 'D' }}, {{ cartoon: {{ color: '#2b6c59' }} }});
    viewer.setStyle({{ chain: 'E' }}, {{ cartoon: {{ color: '#7ec2ad' }} }});
  }}

  function clearSurface() {{
    if (currentSurface !== null) {{
      viewer.removeSurface(currentSurface);
      currentSurface = null;
    }}
  }}

  function clearAnnotations() {{
    activeLabels.forEach(label => viewer.removeLabel(label));
    activeShapes.forEach(shape => viewer.removeShape(shape));
    activeLabels = [];
    activeShapes = [];
  }}

  function colorRange(chainId, range, color) {{
    if (!range || range.length !== 2) return;
    viewer.setStyle({{ chain: chainId, resi: `${{range[0]}}-${{range[1]}}` }}, {{
      cartoon: {{ color }},
      stick: {{ colorscheme: color, radius: 0.16 }}
    }});
  }}

  function applyCDRStyle() {{
    applyChainStyle();
    if (!viewerMeta) return;
    const alpha = viewerMeta.cdr_ranges.D || {{}};
    const beta = viewerMeta.cdr_ranges.E || {{}};
    colorRange('D', alpha.CDR1, '#1f6f5f');
    colorRange('D', alpha.CDR2, '#2e8a78');
    colorRange('D', alpha.CDR3, '#45a18e');
    colorRange('E', beta.CDR1, '#a9652d');
    colorRange('E', beta.CDR2, '#c7772f');
    colorRange('E', beta.CDR3, '#df8d2d');
  }}

  function applyGrooveStyle() {{
    applyChainStyle();
    if (!viewerMeta) return;
    const grooveList = viewerMeta.groove_residues || [];
    if (grooveList.length) {{
      viewer.setStyle({{ chain: 'A', resi: grooveList }}, {{
        cartoon: {{ color: '#f06b54' }},
        stick: {{ colorscheme: '#f06b54', radius: 0.18 }}
      }});
    }}
    viewer.setStyle({{ chain: 'C' }}, {{
      cartoon: {{ color: '#ffd166' }},
      stick: {{ colorscheme: '#ffd166', radius: 0.24 }}
    }});
  }}

  function applyHotspotStyle() {{
    applyChainStyle();
    if (!viewerMeta) return;
    (viewerMeta.top_tcr_hotspots || []).forEach(item => {{
      viewer.addStyle({{ chain: item.chain, resi: item.resi }}, {{
        stick: {{ colorscheme: '#f7f5ef', radius: 0.28 }},
        sphere: {{ color: '#f7f5ef', radius: 0.65 }}
      }});
    }});
    (viewerMeta.top_phla_hotspots || []).forEach(item => {{
      viewer.addStyle({{ chain: item.chain, resi: item.resi }}, {{
        stick: {{ colorscheme: '#ff7a5a', radius: 0.28 }},
        sphere: {{ color: '#ff7a5a', radius: 0.65 }}
      }});
    }});
  }}

  function applyFamilyHotspots(family) {{
    applyChainStyle();
    if (!viewerMeta || !viewerMeta.family_hotspots || !viewerMeta.family_hotspots[family]) return;
    const payload = viewerMeta.family_hotspots[family];
    (payload.tcr || []).forEach(item => {{
      viewer.addStyle({{ chain: item.chain, resi: item.resi }}, {{
        stick: {{ colorscheme: '#f7f5ef', radius: 0.26 }},
        sphere: {{ color: FAMILY_COLORS[family] || '#f7f5ef', radius: 0.6 }}
      }});
    }});
    (payload.phla || []).forEach(item => {{
      viewer.addStyle({{ chain: item.chain, resi: item.resi }}, {{
        stick: {{ colorscheme: FAMILY_COLORS[family] || '#ff7a5a', radius: 0.26 }},
        sphere: {{ color: FAMILY_COLORS[family] || '#ff7a5a', radius: 0.62 }}
      }});
    }});
  }}

  function getResidueAnchor(chain, resi) {{
    const caAtoms = viewer.selectedAtoms({{ chain, resi, atom: 'CA' }});
    if (caAtoms && caAtoms.length) return caAtoms[0];
    const atoms = viewer.selectedAtoms({{ chain, resi }});
    if (atoms && atoms.length) return atoms[0];
    return null;
  }}

  function addPairVisual(pair, color) {{
    viewer.addStyle({{ chain: pair.phla_chain, resi: pair.phla_resi }}, {{
      stick: {{ colorscheme: color, radius: 0.28 }},
      sphere: {{ color, radius: 0.62 }}
    }});
    viewer.addStyle({{ chain: pair.tcr_chain, resi: pair.tcr_resi }}, {{
      stick: {{ colorscheme: color, radius: 0.28 }},
      sphere: {{ color, radius: 0.62 }}
    }});

    const phlaAtom = getResidueAnchor(pair.phla_chain, pair.phla_resi);
    const tcrAtom = getResidueAnchor(pair.tcr_chain, pair.tcr_resi);
    if (!phlaAtom || !tcrAtom) return;

    activeShapes.push(
      viewer.addLine({{
        start: {{ x: phlaAtom.x, y: phlaAtom.y, z: phlaAtom.z }},
        end: {{ x: tcrAtom.x, y: tcrAtom.y, z: tcrAtom.z }},
        dashed: true,
        color,
        linewidth: 2.0,
      }})
    );
    activeLabels.push(
      viewer.addLabel(pair.phla_residue, {{
        position: {{ x: phlaAtom.x, y: phlaAtom.y, z: phlaAtom.z }},
        fontColor: '#10231d',
        backgroundColor: 'rgba(255,250,242,0.92)',
        borderColor: color,
        borderThickness: 1,
        fontSize: 11,
        inFront: true,
      }})
    );
    activeLabels.push(
      viewer.addLabel(pair.tcr_residue, {{
        position: {{ x: tcrAtom.x, y: tcrAtom.y, z: tcrAtom.z }},
        fontColor: '#10231d',
        backgroundColor: 'rgba(255,250,242,0.92)',
        borderColor: color,
        borderThickness: 1,
        fontSize: 11,
        inFront: true,
      }})
    );
  }}

  function applyFamilyPairs(family, pairs) {{
    applyChainStyle();
    if (!viewerMeta || !viewerMeta.family_pairs || !viewerMeta.family_pairs[family]) return;
    const color = FAMILY_COLORS[family] || '#f7f5ef';
    (pairs || []).forEach(pair => addPairVisual(pair, color));
  }}

  function updatePairPanelTitle(family, totalCount) {{
    const title = document.getElementById('pair-panel-title');
    const subtitle = document.getElementById('pair-panel-subtitle');
    if (!family) {{
      title.textContent = 'Residue-pair browser';
      subtitle.textContent = 'Choose an interaction family, then select residue pairs to display.';
      return;
    }}
    title.textContent = `${{family}} residue pairs`;
    subtitle.textContent = `${{totalCount}} residue pairs available. Tick the pairs you want to render in the 3D view.`;
  }}

  function renderPairList(family) {{
    activeFamily = family;
    const list = document.getElementById('viewer-pair-list');
    const pairs = (viewerMeta && viewerMeta.family_pairs && viewerMeta.family_pairs[family]) || [];
    updatePairPanelTitle(family, pairs.length);
    if (!pairs.length) {{
      list.innerHTML = '<div class="viewer-pair-empty">No residue pairs available for this interaction family.</div>';
      return;
    }}
    const rows = pairs.map((pair, idx) => {{
      const checked = idx < 5 ? 'checked' : '';
      const pairText = `${{pair.phla_residue}} ↔ ${{pair.tcr_residue}}`;
      const metaText = `${{pair.interaction_class.replaceAll('_', ' ')}} · ${{pair.tcr_region.replaceAll('_', ' ')}}`;
      return `
        <label class="viewer-pair-row">
          <input type="checkbox" data-pair-id="${{pair.pair_id}}" ${{checked}}>
          <div class="viewer-pair-main">
            <strong>${{pairText}}</strong>
            <span>${{metaText}}</span>
          </div>
          <div class="viewer-pair-score">${{pair.score.toFixed(2)}}</div>
        </label>
      `;
    }});
    list.innerHTML = rows.join('');
  }}

  function getSelectedPairs() {{
    if (!activeFamily || !viewerMeta || !viewerMeta.family_pairs) return [];
    const selectedIds = Array.from(document.querySelectorAll('#viewer-pair-list input[type=\"checkbox\"]:checked'))
      .map(input => input.dataset.pairId);
    const pairMap = new Map((viewerMeta.family_pairs[activeFamily] || []).map(pair => [pair.pair_id, pair]));
    return selectedIds.map(id => pairMap.get(id)).filter(Boolean);
  }}

  function selectTopPairs(limit) {{
    const inputs = Array.from(document.querySelectorAll('#viewer-pair-list input[type=\"checkbox\"]'));
    inputs.forEach((input, idx) => {{
      input.checked = idx < limit;
    }});
  }}

  const FAMILY_COLORS = {{
    contact: '#2d6f5c',
    hbond: '#3f8a77',
    saltbridge: '#b56b3a',
    hydrophobic: '#d7922f',
    pipi: '#8661c1',
    cationpi: '#c05a84'
  }};

  function setStyle(mode) {{
    clearSurface();
    clearAnnotations();
    clearStyles();
    if (mode === 'sticks') {{
      viewer.setStyle({{ chain: 'C' }}, {{ sticks: {{ colorscheme: '#f1b447', radius: 0.28 }} }});
      viewer.setStyle({{ chain: 'D' }}, {{ cartoon: {{ color: '#2b6c59' }}, sticks: {{ colorscheme: '#2b6c59', radius: 0.18 }} }});
      viewer.setStyle({{ chain: 'E' }}, {{ cartoon: {{ color: '#7ec2ad' }}, sticks: {{ colorscheme: '#7ec2ad', radius: 0.18 }} }});
      viewer.setStyle({{ chain: 'A' }}, {{ cartoon: {{ color: '#c8a77e' }} }});
      viewer.setStyle({{ chain: 'B' }}, {{ cartoon: {{ color: '#9b7f63' }} }});
    }} else if (mode === 'surface') {{
      applyChainStyle();
      currentSurface = viewer.addSurface($3Dmol.SurfaceType.VDW, {{ opacity: 0.2, color: '#d5c4ad' }}, {{ chain: 'A' }});
    }} else if (mode === 'cdr') {{
      applyCDRStyle();
    }} else if (mode === 'groove') {{
      applyGrooveStyle();
    }} else if (mode === 'hotspot') {{
      applyHotspotStyle();
    }} else {{
      applyChainStyle();
    }}
    scheduleFrame();
  }}

  Promise.all([
    fetch(pdbUrl).then(resp => resp.text()),
    fetch(metaUrl).then(resp => resp.json())
  ]).then(([pdbText, meta]) => {{
    viewerMeta = meta;
    viewer = $3Dmol.createViewer(viewerElement, {{
      backgroundColor: 'rgba(17, 24, 21, 0.08)'
    }});
    viewer.addModel(pdbText, 'pdb');
    setStyle('chain');
    scheduleFrame();
    setTimeout(() => {{
      scheduleFrame();
    }}, 80);
  }});

  document.querySelectorAll('.viewer-toolbar button').forEach(button => {{
    button.addEventListener('click', () => {{
      if (!viewer) return;
      const action = button.dataset.action;
      const style = button.dataset.style;
      if (action === 'reset') {{
        clearSurface();
        clearAnnotations();
        clearStyles();
        applyChainStyle();
        scheduleFrame();
        return;
      }}
      setStyle(style);
    }});
  }});

  document.querySelectorAll('.viewer-subtoolbar button').forEach(button => {{
    button.addEventListener('click', () => {{
      if (!viewer) return;
      const family = button.dataset.family;
      clearSurface();
      clearAnnotations();
      clearStyles();
      renderPairList(family);
      applyFamilyPairs(family, getSelectedPairs());
      scheduleFrame();
    }});
  }});

  document.getElementById('pair-select-top').addEventListener('click', () => {{
    selectTopPairs(5);
  }});

  document.getElementById('pair-clear-all').addEventListener('click', () => {{
    document.querySelectorAll('#viewer-pair-list input[type="checkbox"]').forEach(input => {{
      input.checked = false;
    }});
  }});

  document.getElementById('pair-apply').addEventListener('click', () => {{
    if (!viewer || !activeFamily) return;
    clearSurface();
    clearAnnotations();
    clearStyles();
    applyFamilyPairs(activeFamily, getSelectedPairs());
    scheduleFrame();
  }});

  window.addEventListener('resize', () => {{
    scheduleFrame();
  }});

  const navLinks = Array.from(document.querySelectorAll('.page-nav a'));
  const sectionPairs = Array.from(
    navLinks
      .map(link => {{
        const id = link.getAttribute('href');
        if (!id || !id.startsWith('#')) return null;
        const section = document.querySelector(id);
        return section ? {{ section, link }} : null;
      }})
      .filter(Boolean)
  );

  if (sectionPairs.length) {{
    const setActiveLink = link => {{
      navLinks.forEach(item => item.classList.remove('active'));
      if (link) link.classList.add('active');
    }};

    const updateActiveNav = () => {{
      const scrollY = window.scrollY;
      const navOffset = 96;
      const probeLine = scrollY + navOffset;

      let activePair = sectionPairs[0];
      for (const pair of sectionPairs) {{
        const top = pair.section.offsetTop;
        const bottom = top + pair.section.offsetHeight;
        if (probeLine >= top && probeLine < bottom) {{
          activePair = pair;
          break;
        }}
        if (probeLine >= top) {{
          activePair = pair;
        }}
      }}

      const nearBottom =
        window.innerHeight + window.scrollY >= document.documentElement.scrollHeight - 8;
      if (nearBottom) {{
        activePair = sectionPairs[sectionPairs.length - 1];
      }}

      setActiveLink(activePair.link);
    }};

    window.addEventListener('scroll', updateActiveNav, {{ passive: true }});
    window.addEventListener('resize', updateActiveNav);
    navLinks.forEach(link => {{
      link.addEventListener('click', () => setActiveLink(link));
    }});
    updateActiveNav();
  }}
</script>
</html>
"""

    output_path = overview_dir / "interaction_report_demo.html"
    output_path.write_text(html_text, encoding="utf-8")
    return output_path


def main() -> None:
    args = parse_args()
    build_interaction_html_report(
        base_dir=args.base_dir,
        system_id=args.system_id,
        source_pdb=args.source_pdb,
        bsa_root=args.bsa_root,
        rmsf_root=args.rmsf_root,
        identity_root=args.identity_root,
        include_sections=parse_section_list(args.include_sections),
        exclude_sections=parse_section_list(args.exclude_sections),
    )


if __name__ == "__main__":
    main()
