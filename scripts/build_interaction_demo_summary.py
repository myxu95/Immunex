#!/usr/bin/env python3
"""构建单体系 interaction 报告的紧凑展示入口。"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


FAMILIES = {
    "contact": {
        "label": "Coarse contact",
        "report": "analysis/contacts/contact_report.csv",
        "pairs": "analysis/contacts/residue_contact_frequencies.csv",
        "heatmaps": [
            "analysis/contacts/pep_TCRa_heatmap.png",
            "analysis/contacts/pep_TCRb_heatmap.png",
            "analysis/contacts/groove_TCRa_heatmap.png",
            "analysis/contacts/groove_TCRb_heatmap.png",
        ],
    },
    "hbond": {
        "label": "Hydrogen bond",
        "report": "analysis/interactions/hydrogen_bonds/hbond_report.csv",
        "pairs": "analysis/interactions/hydrogen_bonds/residue_pair_hbonds.csv",
        "heatmaps": [
            "analysis/interactions/hydrogen_bonds/pep_TCRa_heatmap.png",
            "analysis/interactions/hydrogen_bonds/pep_TCRb_heatmap.png",
            "analysis/interactions/hydrogen_bonds/groove_TCRa_heatmap.png",
            "analysis/interactions/hydrogen_bonds/groove_TCRb_heatmap.png",
        ],
    },
    "saltbridge": {
        "label": "Salt bridge",
        "report": "analysis/interactions/salt_bridges/salt_bridge_report.csv",
        "pairs": "analysis/interactions/salt_bridges/residue_pair_salt_bridges.csv",
        "heatmaps": [
            "analysis/interactions/salt_bridges/groove_TCRa_heatmap.png",
            "analysis/interactions/salt_bridges/groove_TCRb_heatmap.png",
            "analysis/interactions/salt_bridges/hla_TCRa_heatmap.png",
            "analysis/interactions/salt_bridges/hla_TCRb_heatmap.png",
        ],
    },
    "hydrophobic": {
        "label": "Hydrophobic contact",
        "report": "analysis/interactions/hydrophobic_contacts/hydrophobic_report.csv",
        "pairs": "analysis/interactions/hydrophobic_contacts/residue_pair_hydrophobic_contacts.csv",
        "heatmaps": [
            "analysis/interactions/hydrophobic_contacts/groove_TCRa_heatmap.png",
            "analysis/interactions/hydrophobic_contacts/groove_TCRb_heatmap.png",
            "analysis/interactions/hydrophobic_contacts/pep_TCRb_heatmap.png",
        ],
    },
    "pipi": {
        "label": "Pi-pi",
        "report": "analysis/interactions/pi_interactions/pi_pi_report.csv",
        "pairs": "analysis/interactions/pi_interactions/residue_pair_pi_pi.csv",
        "heatmaps": [],
    },
    "cationpi": {
        "label": "Cation-pi",
        "report": "analysis/interactions/cation_pi_interactions/cation_pi_report.csv",
        "pairs": "analysis/interactions/cation_pi_interactions/residue_pair_cation_pi.csv",
        "heatmaps": [
            "analysis/interactions/cation_pi_interactions/groove_TCRa_heatmap.png",
            "analysis/interactions/cation_pi_interactions/groove_TCRb_heatmap.png",
            "analysis/interactions/cation_pi_interactions/hla_TCRa_heatmap.png",
            "analysis/interactions/cation_pi_interactions/hla_TCRb_heatmap.png",
        ],
    },
}


def count_rows(path: Path) -> int:
    if not path.exists():
        return 0
    with path.open("r", encoding="utf-8") as handle:
        return max(sum(1 for _ in handle) - 1, 0)


def ensure_symlink(dst: Path, src: Path) -> None:
    if dst.exists() or dst.is_symlink():
        dst.unlink()
    dst.symlink_to(src.resolve())


def build_summary(base_dir: Path, system_id: str) -> None:
    overview_dir = base_dir / "overview"
    reports_dir = overview_dir / "reports"
    heatmaps_dir = overview_dir / "heatmaps"
    overview_dir.mkdir(exist_ok=True)
    reports_dir.mkdir(exist_ok=True)
    heatmaps_dir.mkdir(exist_ok=True)

    summary_rows = []

    for family, spec in FAMILIES.items():
        task_dir = base_dir / family / system_id
        report_path = task_dir / spec["report"]
        pair_path = task_dir / spec["pairs"]

        report_rows = count_rows(report_path)
        pair_rows = count_rows(pair_path)

        report_link = reports_dir / f"{family}_report.csv"
        if report_path.exists():
            ensure_symlink(report_link, report_path)

        linked_heatmaps = []
        for heatmap_rel in spec["heatmaps"]:
            heatmap_src = task_dir / heatmap_rel
            if not heatmap_src.exists():
                continue
            dst_name = f"{family}_{heatmap_src.name}"
            heatmap_link = heatmaps_dir / dst_name
            ensure_symlink(heatmap_link, heatmap_src)
            linked_heatmaps.append(heatmap_link.relative_to(overview_dir).as_posix())

        summary_rows.append(
            {
                "family": family,
                "label": spec["label"],
                "pair_rows": pair_rows,
                "report_rows": report_rows,
                "report_path": report_path.relative_to(base_dir).as_posix(),
                "overview_report": report_link.relative_to(overview_dir).as_posix()
                if report_path.exists()
                else "",
                "heatmap_count": len(linked_heatmaps),
                "heatmaps": "; ".join(linked_heatmaps),
            }
        )

    csv_path = overview_dir / "interaction_overview.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "family",
                "label",
                "pair_rows",
                "report_rows",
                "report_path",
                "overview_report",
                "heatmap_count",
                "heatmaps",
            ],
        )
        writer.writeheader()
        writer.writerows(summary_rows)

    json_path = overview_dir / "interaction_overview.json"
    json_path.write_text(
        json.dumps(
            {
                "system_id": system_id,
                "families": summary_rows,
            },
            ensure_ascii=False,
            indent=2,
        ),
        encoding="utf-8",
    )

    md_lines = [
        "# Interaction Demo Overview",
        "",
        f"- System: `{system_id}`",
        f"- Root directory: `{base_dir}`",
        "",
        "## Family Overview",
        "",
        "| family | type | raw pairs | report rows | overview report |",
        "| --- | --- | ---: | ---: | --- |",
    ]
    for row in summary_rows:
        md_lines.append(
            f"| `{row['family']}` | {row['label']} | {row['pair_rows']} | {row['report_rows']} | `{row['overview_report']}` |"
        )

    md_lines.extend(
        [
            "",
            "## Suggested Reading Order",
            "",
            "- Start with `interaction_overview.csv` to compare interaction-family scale at a glance.",
            "- Then open the primary report in `reports/` for each family.",
            "- Only dive into the deep region directories when you need `cdr1/cdr2/cdr3/non_cdr` detail.",
            "",
            "## Representative Heatmaps",
            "",
        ]
    )

    for row in summary_rows:
        if not row["heatmaps"]:
            continue
        md_lines.append(f"### {row['family']} / {row['label']}")
        for item in row["heatmaps"].split("; "):
            md_lines.append(f"- `{item}`")
        md_lines.append("")

    (overview_dir / "README.md").write_text("\n".join(md_lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="构建单体系 interaction overview 摘要。")
    parser.add_argument(
        "--base-dir",
        type=Path,
        default=Path("/home/xumy/work/development/Immunex/output/interaction_full_demo"),
        help="包含 contact/hbond/saltbridge 等 family 子目录的结果根目录。",
    )
    parser.add_argument(
        "--system-id",
        type=str,
        default="1ao7_run1",
        help="需要汇总的体系 ID。",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    build_summary(args.base_dir, args.system_id)
