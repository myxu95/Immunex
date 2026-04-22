#!/usr/bin/env python3
"""Inspect available summary and digest files for reporter-query."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

PREFERRED_NAMES = [
    "reporter_context.json",
    "interaction_overview.json",
    "rrcs_summary.json",
    "normal_mode_summary.json",
    "interface_clustering_summary.json",
    "summary_table.csv",
    "cluster_feature_digest.csv",
]
PREFERRED_SUFFIXES = (
    "_summary.json",
    "_digest.csv",
    "_digest.json",
    "_region_summary.csv",
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Inspect query-relevant sources inside an Immunex case directory")
    parser.add_argument("case_dir", help="Case root directory to inspect")
    parser.add_argument("-o", "--output", help="Optional JSON output path")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    case_dir = Path(args.case_dir).resolve()
    if not case_dir.exists():
        raise SystemExit(f"Case directory not found: {case_dir}")

    discovered: list[str] = []
    for path in case_dir.rglob("*"):
        if not path.is_file():
            continue
        name = path.name
        rel = str(path.relative_to(case_dir))
        if name in PREFERRED_NAMES or name.endswith(PREFERRED_SUFFIXES):
            discovered.append(rel)

    discovered = sorted(set(discovered))
    payload = {
        "case_dir": str(case_dir),
        "n_sources": len(discovered),
        "sources": discovered,
    }

    print(json.dumps(payload, indent=2))
    if args.output:
        Path(args.output).write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
