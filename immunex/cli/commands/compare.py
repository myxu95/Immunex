#!/usr/bin/env python3
"""IMN Compare - 双体系结果对比。"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

from immunex.core import PipelineContext
from immunex.pipeline.analysis_pipelines import SystemComparisonPipeline
from scripts.build_comparison_report_html import build_comparison_html_report


def create_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="imn compare",
        description="Compare two analyzed Immunex conditions and generate a comparison report.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  imn compare systems \
    --case-a output/md_standard \
    --case-b output/md_enhanced \
    --label-a "Standard sampling" \
    --label-b "Enhanced sampling" \
    --comparison-mode sampling \
    --comparison-context "Pre/post enhanced sampling" \
    -o output/compare_sampling
        """,
    )
    subparsers = parser.add_subparsers(dest="action", help="Comparison action")

    systems = subparsers.add_parser("systems", help="Compare two analyzed condition roots")
    systems.add_argument("--case-a", type=Path, required=True, help="Condition A result root or overview root")
    systems.add_argument("--case-b", type=Path, required=True, help="Condition B result root or overview root")
    systems.add_argument("--label-a", type=str, default="Condition A", help="Display label for condition A")
    systems.add_argument("--label-b", type=str, default="Condition B", help="Display label for condition B")
    systems.add_argument(
        "--comparison-mode",
        type=str,
        choices=["generic", "mutation", "sampling", "replicate"],
        default="generic",
        help="Comparison framing mode used for report wording",
    )
    systems.add_argument(
        "--comparison-context",
        type=str,
        default="",
        help="Optional free-text comparison context shown on the report cover",
    )
    systems.add_argument("-o", "--output", type=Path, required=True, help="Comparison output directory")

    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose logging")
    return parser


def setup_logging(verbose: bool = False) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format="%(levelname)s: %(message)s")


def handle_systems(args, logger: logging.Logger) -> int:
    output_dir = args.output.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    context = PipelineContext(
        system_id=f"{args.label_a}_vs_{args.label_b}",
        topology="",
        trajectory_raw="",
        output_dir=str(output_dir),
    )
    pipeline = SystemComparisonPipeline(
        case_a_root=str(args.case_a),
        case_b_root=str(args.case_b),
        label_a=args.label_a,
        label_b=args.label_b,
        comparison_mode=args.comparison_mode,
        comparison_context=args.comparison_context,
    )

    logger.info("=" * 60)
    logger.info("IMN System Comparison")
    logger.info("=" * 60)
    logger.info("Condition A: %s", args.case_a)
    logger.info("Condition B: %s", args.case_b)
    logger.info("Comparison mode: %s", args.comparison_mode)
    if args.comparison_context:
        logger.info("Comparison context: %s", args.comparison_context)
    logger.info("Output: %s", output_dir)

    result_context = pipeline.execute(context)
    if result_context.has_errors():
        for error in result_context.errors:
            logger.error(error)
        return 1

    comparison_result = result_context.results["comparison"]
    html_path = build_comparison_html_report(output_dir=output_dir, comparison_result=comparison_result)

    logger.info("Comparison report generated")
    logger.info("HTML: %s", html_path)
    return 0


def main(argv=None) -> int:
    parser = create_parser()
    args = parser.parse_args(argv)
    if not args.action:
        parser.print_help()
        return 1

    setup_logging(args.verbose)
    logger = logging.getLogger(__name__)

    if args.action == "systems":
        return handle_systems(args, logger)

    logger.error("Unknown action: %s", args.action)
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
