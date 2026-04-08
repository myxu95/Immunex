#!/usr/bin/env python3
"""IMN Report - 单体系报告生成。"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

from scripts.build_interaction_demo_summary import build_summary
from scripts.build_interaction_demo_html import build_interaction_html_report, parse_section_list


def create_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="imn report",
        description="Generate single-system HTML reports from Immunex analysis outputs.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  imn report interaction --base-dir output/interaction_case_1OGA --system-id 1OGA_sd_run2

  imn report interaction \
    --base-dir output/interaction_case_1OGA \
    --system-id 1OGA_sd_run2 \
    --source-pdb output/parallel2_preprocess/1OGA_sd_run2/md_processed_converted.pdb \
    --bsa-root output/bsa_demo_1OGA
        """,
    )
    subparsers = parser.add_subparsers(dest="action", help="Report action")

    interaction = subparsers.add_parser("interaction", help="Generate interaction HTML report")
    interaction.add_argument("--base-dir", type=Path, required=True, help="Result root containing family subdirectories")
    interaction.add_argument("--system-id", type=str, required=True, help="System ID to report")
    interaction.add_argument("--source-pdb", type=Path, default=None, help="Optional structure PDB for the 3D viewer")
    interaction.add_argument("--bsa-root", type=Path, default=None, help="Optional BSA result root containing analysis/interface")
    interaction.add_argument("--rmsf-root", type=Path, default=None, help="Optional RMSF result root containing analysis/rmsf")
    interaction.add_argument("--identity-root", type=Path, default=None, help="Optional identity result root containing analysis/identity")
    interaction.add_argument("--cluster-root", type=Path, default=None, help="Optional interface clustering root containing analysis/conformation/interface_clustering")
    interaction.add_argument(
        "--include-sections",
        type=str,
        default=None,
        help="Optional comma-separated section list to include: overview,quality,interface,flexibility,cluster,occupancy,contact,interactions,downloads",
    )
    interaction.add_argument(
        "--exclude-sections",
        type=str,
        default=None,
        help="Optional comma-separated section list to exclude",
    )

    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose logging")
    return parser


def setup_logging(verbose: bool = False) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format="%(levelname)s: %(message)s")


def handle_interaction(args, logger: logging.Logger) -> int:
    base_dir = args.base_dir
    if not base_dir.exists():
        logger.error(f"Base directory not found: {base_dir}")
        return 1

    logger.info("=" * 60)
    logger.info("IMN Interaction Report")
    logger.info("=" * 60)
    logger.info(f"Base directory: {base_dir}")
    logger.info(f"System ID: {args.system_id}")
    if args.bsa_root:
        logger.info(f"BSA root: {args.bsa_root}")
    if args.rmsf_root:
        logger.info(f"RMSF root: {args.rmsf_root}")
    if args.identity_root:
        logger.info(f"Identity root: {args.identity_root}")
    if args.cluster_root:
        logger.info(f"Cluster root: {args.cluster_root}")
    if args.include_sections:
        logger.info(f"Include sections: {args.include_sections}")
    if args.exclude_sections:
        logger.info(f"Exclude sections: {args.exclude_sections}")
    logger.info("")

    build_summary(base_dir, args.system_id)
    html_path = build_interaction_html_report(
        base_dir=base_dir,
        system_id=args.system_id,
        source_pdb=args.source_pdb,
        bsa_root=args.bsa_root,
        rmsf_root=args.rmsf_root,
        identity_root=args.identity_root,
        cluster_root=args.cluster_root,
        include_sections=parse_section_list(args.include_sections),
        exclude_sections=parse_section_list(args.exclude_sections),
    )

    logger.info("Interaction report generated")
    logger.info(f"HTML: {html_path}")
    return 0


def main(argv=None) -> int:
    parser = create_parser()
    args = parser.parse_args(argv)

    if not args.action:
        parser.print_help()
        return 1

    setup_logging(args.verbose)
    logger = logging.getLogger(__name__)

    if args.action == "interaction":
        return handle_interaction(args, logger)

    logger.error(f"Unknown action: {args.action}")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
