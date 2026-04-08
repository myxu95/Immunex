#!/usr/bin/env python3
"""IMN NMA - 基于 ENM 的 normal mode / PRS 分析。"""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

from immunex.core import PipelineContext
from immunex.pipeline import NormalModePipeline


def create_parser():
    parser = argparse.ArgumentParser(
        prog="imn nma",
        description="Normal mode / ENM analysis for pHLA-TCR complexes",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  imn nma --structure md_processed_converted.pdb -o ./output/nma_demo

  imn nma --structure md_processed_converted.pdb -o ./output/nma_demo --cutoff 12.0 --low-modes 12
        """,
    )
    parser.add_argument("--structure", required=True, metavar="FILE", help="Representative structure PDB")
    parser.add_argument("-o", "--output", required=True, metavar="DIR", help="Output directory")
    parser.add_argument("--cutoff", type=float, default=10.0, metavar="A", help="ENM cutoff in Angstrom (default: 10.0)")
    parser.add_argument("--low-modes", type=int, default=10, metavar="N", help="Number of low-frequency modes (default: 10)")
    parser.add_argument("--prs-forces", type=int, default=8, metavar="N", help="Number of force directions for PRS (default: 8)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose output")
    return parser


def setup_logging(verbose: bool = False):
    logging.basicConfig(level=logging.DEBUG if verbose else logging.INFO, format="%(levelname)s: %(message)s")


def main(argv=None):
    parser = create_parser()
    args = parser.parse_args(argv)
    setup_logging(args.verbose)
    logger = logging.getLogger(__name__)

    structure_path = Path(args.structure)
    if not structure_path.exists():
        logger.error(f"Structure file not found: {structure_path}")
        return 1

    output_dir = Path(args.output).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    system_id = structure_path.stem

    context = PipelineContext(
        system_id=system_id,
        topology=str(structure_path.resolve()),
        trajectory_raw=str(structure_path.resolve()),
        structure_pdb=str(structure_path.resolve()),
        output_dir=str(output_dir),
    )

    pipeline = NormalModePipeline(
        cutoff_angstrom=args.cutoff,
        n_low_modes=args.low_modes,
        prs_force_directions=args.prs_forces,
        auto_identify_chains=True,
        auto_detect_cdr=True,
    )
    result = pipeline.execute(context)
    if result.has_errors():
        logger.error("Normal mode analysis failed")
        for error in result.errors:
            logger.error(error)
        return 1

    nma_result = result.results.get("normal_mode", {})
    logger.info("Normal mode analysis completed")
    logger.info(f"Residue CSV: {nma_result.get('residue_csv')}")
    logger.info(f"Region summary: {nma_result.get('region_summary_csv')}")
    logger.info(f"Summary JSON: {nma_result.get('summary_json')}")

    summary_path = nma_result.get("summary_json")
    if summary_path:
        summary = json.loads(Path(summary_path).read_text(encoding="utf-8"))
        logger.info(
            json.dumps(
                {
                    "n_residues": summary.get("n_residues"),
                    "network_edges": summary.get("network_edges"),
                    "top_key_residues": summary.get("top_key_residues", [])[:5],
                },
                ensure_ascii=False,
                indent=2,
            )
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
