#!/usr/bin/env python3
"""构建 HLA 派生参考库。"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from immunex.analysis.topology import derive_hla_reference_library


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="从原始 HLA FASTA 构建派生参考库。")
    parser.add_argument(
        "--raw-dir",
        type=Path,
        default=Path("/home/xumy/work/development/Immunex/immunex/data/reference/hla/raw/class_i"),
        help="原始 HLA FASTA 目录。",
    )
    parser.add_argument(
        "--derived-dir",
        type=Path,
        default=Path("/home/xumy/work/development/Immunex/immunex/data/reference/hla/derived"),
        help="派生库输出目录。",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    fasta_path, metadata_path = derive_hla_reference_library(args.raw_dir, args.derived_dir)
    print(f"Derived FASTA: {fasta_path}")
    print(f"Derived metadata: {metadata_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
