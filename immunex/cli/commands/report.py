#!/usr/bin/env python3
"""IMN Report - 单体系报告生成与轻量问答服务。"""

from __future__ import annotations

import argparse
import json
import logging
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Type

from immunex.analysis.reporter import QueryAnswerBuilder
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

  imn report serve --overview-dir output/1OGA_demo_bundle/report/interaction_case_1OGA/overview
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
    interaction.add_argument("--rrcs-root", type=Path, default=None, help="Optional RRCS result root containing analysis/interactions/rrcs")
    interaction.add_argument("--cluster-root", type=Path, default=None, help="Optional interface clustering root containing analysis/conformation/interface_clustering")
    interaction.add_argument(
        "--include-sections",
        type=str,
        default=None,
        help="Optional comma-separated section list to include: overview,quality,interface,flexibility,rrcs,cluster,occupancy,contact,interactions,downloads",
    )
    interaction.add_argument(
        "--exclude-sections",
        type=str,
        default=None,
        help="Optional comma-separated section list to exclude",
    )

    serve = subparsers.add_parser("serve", help="Serve a generated report directory with a lightweight assistant API")
    serve.add_argument("--overview-dir", type=Path, required=True, help="Overview directory containing interaction_report_demo.html")
    serve.add_argument("--case-dir", type=Path, default=None, help="Optional case directory. Defaults to overview_dir.parent")
    serve.add_argument("--host", type=str, default="127.0.0.1", help="Bind address")
    serve.add_argument("--port", type=int, default=8770, help="Serve port")
    serve.add_argument("--config", type=Path, default=None, help="Optional config file path")
    serve.add_argument("--llm", action="store_true", help="Enable LLM-enhanced answers")
    serve.add_argument("--llm-provider", type=str, choices=["openai", "anthropic"], help="LLM provider")
    serve.add_argument("--llm-model", type=str, help="LLM model name")
    serve.add_argument("--llm-api-key", type=str, help="LLM API key (or use env var)")

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
    if args.rrcs_root:
        logger.info(f"RRCS root: {args.rrcs_root}")
    if args.cluster_root:
        logger.info(f"Cluster root: {args.cluster_root}")
    if args.include_sections:
        logger.info(f"Include sections: {args.include_sections}")
    if args.exclude_sections:
        logger.info(f"Exclude sections: {args.exclude_sections}")
    logger.info("")

    html_path = build_interaction_html_report(
        base_dir=base_dir,
        system_id=args.system_id,
        source_pdb=args.source_pdb,
        bsa_root=args.bsa_root,
        rmsf_root=args.rmsf_root,
        identity_root=args.identity_root,
        rrcs_root=args.rrcs_root,
        cluster_root=args.cluster_root,
        include_sections=parse_section_list(args.include_sections),
        exclude_sections=parse_section_list(args.exclude_sections),
    )

    logger.info("Interaction report generated")
    logger.info(f"HTML: {html_path}")
    return 0


def _create_report_handler(overview_dir: Path, case_dir: Path, logger: logging.Logger, llm_config: dict = None) -> Type[SimpleHTTPRequestHandler]:
    class ReportHandler(SimpleHTTPRequestHandler):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, directory=str(overview_dir), **kwargs)

        def log_message(self, format: str, *args) -> None:
            logger.info("report-server: " + format, *args)

        def do_POST(self) -> None:
            if self.path.rstrip("/") != "/api/assistant/query":
                self.send_error(404, "Unknown assistant endpoint")
                return
            content_length = int(self.headers.get("Content-Length", "0") or 0)
            try:
                raw = self.rfile.read(content_length) if content_length > 0 else b"{}"
                payload = json.loads(raw.decode("utf-8"))
            except json.JSONDecodeError:
                self._write_json({"error": "Invalid JSON payload"}, status=400)
                return
            question = str(payload.get("question", "")).strip()
            if not question:
                self._write_json({"error": "Question is required"}, status=400)
                return
            try:
                # Initialize QueryAnswerBuilder with LLM config if provided
                if llm_config:
                    builder = QueryAnswerBuilder(
                        case_dir,
                        use_llm=llm_config.get("enabled", False),
                        llm_provider=llm_config.get("provider", "openai"),
                        llm_model=llm_config.get("model"),
                        llm_api_key=llm_config.get("api_key"),
                    )
                else:
                    builder = QueryAnswerBuilder(case_dir)

                result = builder.answer(question).to_dict()
                result.setdefault("skill", "reporter-query")
                result.setdefault("sources", result.get("sources_used", []))
                self._write_json(result, status=200)
            except Exception as exc:  # pragma: no cover - runtime guard
                logger.exception("Assistant query failed")
                self._write_json({"error": f"Assistant query failed: {exc}"}, status=500)

        def _write_json(self, payload: dict, status: int = 200) -> None:
            body = json.dumps(payload, ensure_ascii=False).encode("utf-8")
            self.send_response(status)
            self.send_header("Content-Type", "application/json; charset=utf-8")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)

    return ReportHandler


def handle_serve(args, logger: logging.Logger) -> int:
    overview_dir = args.overview_dir.resolve()
    if not overview_dir.exists():
        logger.error(f"Overview directory not found: {overview_dir}")
        return 1
    html_path = overview_dir / "interaction_report_demo.html"
    if not html_path.exists():
        logger.error(f"Report HTML not found: {html_path}")
        return 1
    case_dir = args.case_dir.resolve() if args.case_dir else overview_dir.parent.resolve()
    if not case_dir.exists():
        logger.error(f"Case directory not found: {case_dir}")
        return 1

    # Load configuration
    from immunex.analysis.reporter.config import load_config

    config = load_config(args.config if hasattr(args, 'config') else None)

    # Override with CLI arguments
    if hasattr(args, 'llm') and args.llm:
        config.llm.enabled = True
    if hasattr(args, 'llm_provider') and args.llm_provider:
        config.llm.provider = args.llm_provider
    if hasattr(args, 'llm_model') and args.llm_model:
        config.llm.model = args.llm_model
    if hasattr(args, 'llm_api_key') and args.llm_api_key:
        config.llm.api_key = args.llm_api_key

    # Use config values for host/port if not overridden
    host = args.host if args.host != "127.0.0.1" else config.server.host
    port = args.port if args.port != 8770 else config.server.port

    # Prepare LLM config dict
    llm_config = None
    if config.llm.enabled:
        llm_config = {
            "enabled": True,
            "provider": config.llm.provider,
            "model": config.llm.model,
            "api_key": config.llm.api_key,
        }

    handler_cls = _create_report_handler(overview_dir, case_dir, logger, llm_config)
    server = ThreadingHTTPServer((host, port), handler_cls)
    logger.info("=" * 60)
    logger.info("IMN Report Server")
    logger.info("=" * 60)
    logger.info(f"Overview directory: {overview_dir}")
    logger.info(f"Case directory: {case_dir}")
    logger.info(f"Listening on: http://{host}:{port}/interaction_report_demo.html")
    logger.info("Assistant API: POST /api/assistant/query")
    if config.llm.enabled:
        logger.info(f"LLM: Enabled ({config.llm.provider} - {config.llm.model or 'default model'})")
    else:
        logger.info("LLM: Disabled (using rule-based answers)")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        logger.info("Stopping report server")
    finally:
        server.server_close()
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
    if args.action == "serve":
        return handle_serve(args, logger)

    logger.error(f"Unknown action: {args.action}")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
