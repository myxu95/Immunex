#!/usr/bin/env python3
"""
IMN Batch - Batch Processing Workflows

Usage:
    imn batch preprocess <base_dir> [options]
    imn batch quality <base_dir> [options]
    imn batch contact <base_dir> [options]
    imn batch hbond <base_dir> [options]
    imn batch saltbridge <base_dir> [options]
    imn batch hydrophobic <base_dir> [options]
    imn batch pipi <base_dir> [options]
    imn batch cationpi <base_dir> [options]
    imn batch angle <base_dir> [options]
    imn batch bsa <base_dir> [options]

Author: Immunex Development Team
Date: 2026-03-15
"""

import argparse
import json
import logging
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from immunex.core.task_discovery import discover_tasks
from immunex.pipeline import (
    BatchExecutor,
    ContactFrequencyPipeline,
    DockingAnglePipeline,
    HydrogenBondInteractionPipeline,
    SaltBridgeInteractionPipeline,
    HydrophobicInteractionPipeline,
    PiStackingInteractionPipeline,
    CationPiInteractionPipeline,
    BSAPipeline,
    PreprocessQualityPipeline,
    QualityAssessmentPipeline,
)


def create_parser():
    """Create argument parser for batch subcommand"""
    parser = argparse.ArgumentParser(
        prog='imn batch',
        description='Batch processing workflows for multiple MD tasks',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  Batch PBC processing:
    imn batch preprocess /data/md_tasks/ --workers 4

  Batch quality assessment:
    imn batch quality /data/processed_md/ --workers 4

  Batch contact analysis:
    imn batch contact /data/processed_md/ --workers 4 --cutoff 4.5

  Batch hydrogen-bond interaction analysis:
    imn batch hbond /data/processed_md/ --workers 4 --stride 1

  Batch salt-bridge interaction analysis:
    imn batch saltbridge /data/processed_md/ --workers 4 --stride 1

  Batch hydrophobic contact analysis:
    imn batch hydrophobic /data/processed_md/ --workers 4 --stride 1

  Batch pi-pi interaction analysis:
    imn batch pipi /data/processed_md/ --workers 4 --stride 1

  Batch cation-pi interaction analysis:
    imn batch cationpi /data/processed_md/ --workers 4 --stride 1

  Batch docking angle analysis:
    imn batch angle /data/md_tasks/ --workers 4 --stride 10 --limit 10

  Batch buried surface area analysis:
    imn batch bsa /data/processed_md/ --workers 4 --stride 10

  Combined workflow:
    imn batch workflow /data/md_tasks/ --workers 4
        """
    )

    subparsers = parser.add_subparsers(dest='action', help='Batch action')

    preproc_parser = subparsers.add_parser('preprocess', help='Batch PBC processing')
    preproc_parser.add_argument('base_dir', help='Base directory containing MD tasks')
    preproc_parser.add_argument('--workers', type=int, default=4, help='Number of parallel workers (default: 4)')
    preproc_parser.add_argument('--method', choices=['2step', '3step'], default='2step', help='PBC method (default: 2step)')
    preproc_parser.add_argument('--rmsd-selection', default='backbone', help='RMSD atom selection for preprocessing quality assessment (default: backbone)')
    preproc_parser.add_argument('-o', '--output-dir', help='Override output directory')

    quality_parser = subparsers.add_parser('quality', help='Batch quality assessment')
    quality_parser.add_argument('base_dir', help='Base directory containing processed trajectories')
    quality_parser.add_argument('--workers', type=int, default=4, help='Number of parallel workers (default: 4)')
    quality_parser.add_argument('-o', '--output-dir', help='Override output directory')
    quality_parser.add_argument('--summary', help='Save summary JSON file')
    quality_parser.add_argument('--stride', type=int, default=1, help='Validation stride (default: 1)')

    contact_parser = subparsers.add_parser('contact', help='Batch contact analysis')
    contact_parser.add_argument('base_dir', help='Base directory containing processed trajectories')
    contact_parser.add_argument('--workers', type=int, default=4, help='Number of parallel workers (default: 4)')
    contact_parser.add_argument('--cutoff', type=float, default=4.5, help='Contact cutoff in Angstrom (default: 4.5)')
    contact_parser.add_argument('--stride', type=int, default=1, help='Trajectory stride for contact analysis (default: 1)')
    contact_parser.add_argument('--min-frequency', type=float, default=0.0, help='Minimum contact frequency to retain (default: 0.0)')
    contact_parser.add_argument('-o', '--output-dir', help='Override output directory')
    contact_parser.add_argument('--summary', help='Save summary JSON file')

    hbond_parser = subparsers.add_parser('hbond', help='Batch hydrogen-bond interaction analysis')
    hbond_parser.add_argument('base_dir', help='Base directory containing processed trajectories')
    hbond_parser.add_argument('--workers', type=int, default=4, help='Number of parallel workers (default: 4)')
    hbond_parser.add_argument('--stride', type=int, default=1, help='Trajectory stride for hydrogen-bond analysis (default: 1)')
    hbond_parser.add_argument('--distance-cutoff', type=float, default=3.5, help='Hydrogen-bond donor-acceptor cutoff in Angstrom (default: 3.5)')
    hbond_parser.add_argument('--angle-cutoff', type=float, default=150.0, help='Hydrogen-bond angle cutoff in degrees (default: 150.0)')
    hbond_parser.add_argument('-o', '--output-dir', help='Override output directory')
    hbond_parser.add_argument('--summary', help='Save summary JSON file')

    saltbridge_parser = subparsers.add_parser('saltbridge', help='Batch salt-bridge interaction analysis')
    saltbridge_parser.add_argument('base_dir', help='Base directory containing processed trajectories')
    saltbridge_parser.add_argument('--workers', type=int, default=4, help='Number of parallel workers (default: 4)')
    saltbridge_parser.add_argument('--stride', type=int, default=1, help='Trajectory stride for salt-bridge analysis (default: 1)')
    saltbridge_parser.add_argument('--distance-cutoff', type=float, default=4.0, help='Salt-bridge cutoff in Angstrom (default: 4.0)')
    saltbridge_parser.add_argument('-o', '--output-dir', help='Override output directory')
    saltbridge_parser.add_argument('--summary', help='Save summary JSON file')

    hydrophobic_parser = subparsers.add_parser('hydrophobic', help='Batch hydrophobic contact analysis')
    hydrophobic_parser.add_argument('base_dir', help='Base directory containing processed trajectories')
    hydrophobic_parser.add_argument('--workers', type=int, default=4, help='Number of parallel workers (default: 4)')
    hydrophobic_parser.add_argument('--stride', type=int, default=1, help='Trajectory stride for hydrophobic analysis (default: 1)')
    hydrophobic_parser.add_argument('--distance-cutoff', type=float, default=4.5, help='Hydrophobic contact cutoff in Angstrom (default: 4.5)')
    hydrophobic_parser.add_argument('-o', '--output-dir', help='Override output directory')
    hydrophobic_parser.add_argument('--summary', help='Save summary JSON file')

    pipi_parser = subparsers.add_parser('pipi', help='Batch pi-pi interaction analysis')
    pipi_parser.add_argument('base_dir', help='Base directory containing processed trajectories')
    pipi_parser.add_argument('--workers', type=int, default=4, help='Number of parallel workers (default: 4)')
    pipi_parser.add_argument('--stride', type=int, default=1, help='Trajectory stride for pi-pi analysis (default: 1)')
    pipi_parser.add_argument('--distance-cutoff', type=float, default=6.5, help='pi-pi centroid cutoff in Angstrom (default: 6.5)')
    pipi_parser.add_argument('-o', '--output-dir', help='Override output directory')
    pipi_parser.add_argument('--summary', help='Save summary JSON file')

    cationpi_parser = subparsers.add_parser('cationpi', help='Batch cation-pi interaction analysis')
    cationpi_parser.add_argument('base_dir', help='Base directory containing processed trajectories')
    cationpi_parser.add_argument('--workers', type=int, default=4, help='Number of parallel workers (default: 4)')
    cationpi_parser.add_argument('--stride', type=int, default=1, help='Trajectory stride for cation-pi analysis (default: 1)')
    cationpi_parser.add_argument('--distance-cutoff', type=float, default=6.0, help='Cation-pi centroid cutoff in Angstrom (default: 6.0)')
    cationpi_parser.add_argument('--normal-angle-cutoff', type=float, default=60.0, help='Cation-pi ring-normal angle cutoff in degrees (default: 60.0)')
    cationpi_parser.add_argument('-o', '--output-dir', help='Override output directory')
    cationpi_parser.add_argument('--summary', help='Save summary JSON file')

    angle_parser = subparsers.add_parser('angle', help='Batch docking angle analysis')
    angle_parser.add_argument('base_dir', help='Base directory containing MD tasks')
    angle_parser.add_argument('--workers', type=int, default=4, help='Number of parallel workers (default: 4)')
    angle_parser.add_argument('--stride', type=int, default=1, help='Trajectory stride for angle analysis (default: 1)')
    angle_parser.add_argument('--limit', type=int, help='Only analyze the first N valid tasks')
    angle_parser.add_argument('--print-each-frame', action='store_true', help='Print docking angles for every analyzed frame')
    angle_parser.add_argument('-o', '--output-dir', help='Override output directory')
    angle_parser.add_argument('--summary', help='Save summary JSON file')

    bsa_parser = subparsers.add_parser('bsa', help='Batch buried surface area analysis')
    bsa_parser.add_argument('base_dir', help='Base directory containing processed trajectories')
    bsa_parser.add_argument('--workers', type=int, default=4, help='Number of parallel workers (default: 4)')
    bsa_parser.add_argument('--stride', type=int, default=1, help='Trajectory stride for BSA analysis (default: 1)')
    bsa_parser.add_argument('--probe-radius', type=float, default=1.4, help='SASA probe radius in Angstrom (default: 1.4)')
    bsa_parser.add_argument('--time-unit', choices=['ps', 'ns'], default='ps', help='Output time unit (default: ps)')
    bsa_parser.add_argument('-o', '--output-dir', help='Override output directory')
    bsa_parser.add_argument('--summary', help='Save summary JSON file')

    workflow_parser = subparsers.add_parser('workflow', help='Complete workflow (PBC + Quality)')
    workflow_parser.add_argument('base_dir', help='Base directory containing MD tasks')
    workflow_parser.add_argument('--workers', type=int, default=4, help='Number of parallel workers (default: 4)')
    workflow_parser.add_argument('--method', choices=['2step', '3step'], default='2step', help='PBC method (default: 2step)')
    workflow_parser.add_argument('-o', '--output-dir', help='Override output directory')
    workflow_parser.add_argument('--summary', help='Save workflow quality summary JSON file')
    workflow_parser.add_argument('--stride', type=int, default=1, help='Validation stride (default: 1)')

    parser.add_argument('-v', '--verbose', action='store_true', help='Verbose output')
    return parser


def setup_logging(verbose=False):
    """Setup logging configuration"""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format='%(levelname)s: %(message)s')


def _discover_execution_report(base_dir: Path, required_files: list[str]):
    """Discover tasks for a pipeline stage using the public discovery contract."""
    return discover_tasks(base_dir, required_files=required_files)


def _summarize_context_results(results):
    successful = [result for result in results if not result.has_errors()]
    failed = [result for result in results if result.has_errors()]
    return {
        'total_tasks': len(results),
        'successful': len(successful),
        'failed': len(failed),
        'results': [
            {
                'task_name': result.system_id,
                'status': 'success' if not result.has_errors() else 'failed',
                'processed': result.trajectory_processed,
                'output_dir': result.output_dir,
                'rmsd': result.get_result('rmsd'),
                'preprocess_quality': result.get_result('preprocess_quality'),
                'errors': list(result.errors),
            }
            for result in results
        ],
    }


def _write_preprocess_quality_summary(results, output_root: str):
    output_path = Path(output_root)
    output_path.mkdir(parents=True, exist_ok=True)

    rows = []
    for item in results['results']:
        quality = item.get('preprocess_quality') or {}
        rmsd = item.get('rmsd') or {}
        metrics = quality.get('metrics') or {}
        reports = quality.get('reports') or {}
        rows.append({
            'task_name': item['task_name'],
            'status': item['status'],
            'full_variation_nm': metrics.get('full_variation_nm'),
            'tail90_variation_nm': metrics.get('tail90_variation_nm'),
            'tail90_mean_rmsd_nm': metrics.get('tail90_mean_rmsd_nm'),
            'tail90_std_rmsd_nm': metrics.get('tail90_std_rmsd_nm'),
            'trimmed_start_frame': metrics.get('trimmed_start_frame'),
            'processed': item.get('processed'),
            'rmsd_file': rmsd.get('output_file'),
            'rmsd_plot': metrics.get('plot_file'),
            'quality_report_md': reports.get('markdown'),
            'quality_report_json': reports.get('json'),
            'errors': '; '.join(item.get('errors', [])),
        })

    df = pd.DataFrame(rows).sort_values(
        by=['status', 'tail90_variation_nm', 'task_name'],
        ascending=[True, False, True],
    )
    csv_path = output_path / 'preprocess_rmsd_variation_summary.csv'
    json_path = output_path / 'preprocess_rmsd_variation_summary.json'
    md_path = output_path / 'preprocess_rmsd_variation_summary.md'
    problematic_csv_path = output_path / 'problematic_tasks.csv'
    plot_path = output_path / 'rmsd_tail90_variation_rank.png'

    df.to_csv(csv_path, index=False)
    json_path.write_text(
        json.dumps(rows, indent=2, ensure_ascii=False, default=str),
        encoding='utf-8',
    )
    success_count = int((df['status'] == 'success').sum())
    failed_count = int((df['status'] == 'failed').sum())

    problem_df = df[
        (df['status'] == 'failed') |
        (
            df['status'].eq('success') &
            df['tail90_variation_nm'].notna()
        )
    ].copy()
    if not problem_df.empty:
        problem_df = problem_df.sort_values(by='tail90_variation_nm', ascending=False)
        problem_df.to_csv(problematic_csv_path, index=False)
    else:
        pd.DataFrame(columns=df.columns).to_csv(problematic_csv_path, index=False)

    top_rows = df[['task_name', 'status', 'tail90_variation_nm', 'full_variation_nm', 'tail90_mean_rmsd_nm']].head(20)
    lines = [
        "# Batch Preprocess RMSD Variation Summary",
        "",
        f"- Total tasks: {len(df)}",
        f"- Successful: {success_count}",
        f"- Failed: {failed_count}",
        f"- Tasks with RMSD variation data: {int(df['tail90_variation_nm'].notna().sum())}",
        "",
        "## Most Variable Tasks (Tail90)",
        "",
        "| Task | Status | Tail90 Variation (nm) | Full Variation (nm) | Tail90 Mean RMSD (nm) |",
        "|---|---|---:|---:|---:|",
    ]
    for _, row in top_rows.iterrows():
        tail_var = "" if pd.isna(row['tail90_variation_nm']) else f"{row['tail90_variation_nm']:.3f}"
        full_var = "" if pd.isna(row['full_variation_nm']) else f"{row['full_variation_nm']:.3f}"
        tail_mean = "" if pd.isna(row['tail90_mean_rmsd_nm']) else f"{row['tail90_mean_rmsd_nm']:.3f}"
        lines.append(f"| {row['task_name']} | {row['status']} | {tail_var} | {full_var} | {tail_mean} |")
    lines.append("")

    md_path.write_text('\n'.join(lines), encoding='utf-8')

    plot_df = df[(df['status'] == 'success') & df['tail90_variation_nm'].notna()].head(30)
    if not plot_df.empty:
        fig, ax = plt.subplots(figsize=(12, max(6, len(plot_df) * 0.3)))
        plot_df = plot_df.sort_values(by='tail90_variation_nm', ascending=True)
        ax.barh(plot_df['task_name'], plot_df['tail90_variation_nm'], color='#1f4e79', alpha=0.9)
        ax.set_xlabel('Tail90 RMSD Variation (nm)')
        ax.set_ylabel('Task')
        ax.set_title('RMSD Tail90 Variation Ranking')
        ax.grid(True, axis='x', alpha=0.2)
        fig.tight_layout()
        fig.savefig(plot_path, dpi=300, bbox_inches='tight')
        plt.close(fig)

    return {
        'csv': str(csv_path),
        'json': str(json_path),
        'markdown': str(md_path),
        'problematic_csv': str(problematic_csv_path),
        'plot': str(plot_path),
    }


def _load_angle_timeseries(csv_path: str | None):
    """Load docking angle time series from CSV and compute variation metrics."""
    if not csv_path:
        return {}

    csv_file = Path(csv_path)
    if not csv_file.exists():
        return {}

    df = pd.read_csv(csv_file)
    if df.empty:
        return {}

    total_frames = len(df)
    trim_start = min(total_frames - 1, int(total_frames * 0.1)) if total_frames > 1 else 0
    tail_df = df.iloc[trim_start:].copy()

    crossing_col = 'Crossing(deg)'
    incident_col = 'Incident(deg)'
    if crossing_col not in df.columns or incident_col not in df.columns:
        return {}

    def _variation(series_name: str, frame_df: pd.DataFrame):
        series = frame_df[series_name].dropna()
        if series.empty:
            return None
        return float(series.max() - series.min())

    metrics = {
        'n_frames': total_frames,
        'trimmed_start_frame': trim_start,
        'crossing_full_variation_deg': _variation(crossing_col, df),
        'crossing_tail90_variation_deg': _variation(crossing_col, tail_df),
        'incident_full_variation_deg': _variation(incident_col, df),
        'incident_tail90_variation_deg': _variation(incident_col, tail_df),
    }
    return metrics


def _write_angle_summary(results, output_root: str):
    """Write batch docking-angle summary files."""
    output_path = Path(output_root)
    output_path.mkdir(parents=True, exist_ok=True)

    rows = []
    for item in results['results']:
        angle_result = item.get('docking_angles') or {}
        stats = angle_result.get('statistics') or {}
        output_files = angle_result.get('output_files') or []
        csv_file = next(
            (path for path in output_files if str(path).endswith('docking_angles.csv')),
            None,
        )
        metrics = _load_angle_timeseries(csv_file)

        rows.append({
            'task_name': item['task_name'],
            'status': item['status'],
            'n_frames': metrics.get('n_frames'),
            'trimmed_start_frame': metrics.get('trimmed_start_frame'),
            'crossing_mean_deg': stats.get('crossing_mean'),
            'crossing_std_deg': stats.get('crossing_std'),
            'crossing_full_variation_deg': metrics.get('crossing_full_variation_deg'),
            'crossing_tail90_variation_deg': metrics.get('crossing_tail90_variation_deg'),
            'incident_mean_deg': stats.get('incident_mean'),
            'incident_std_deg': stats.get('incident_std'),
            'incident_full_variation_deg': metrics.get('incident_full_variation_deg'),
            'incident_tail90_variation_deg': metrics.get('incident_tail90_variation_deg'),
            'angles_csv': csv_file,
            'output_dir': item.get('output_dir'),
            'errors': '; '.join(item.get('errors', [])),
        })

    df = pd.DataFrame(rows).sort_values(
        by=['status', 'crossing_tail90_variation_deg', 'incident_tail90_variation_deg', 'task_name'],
        ascending=[True, False, False, True],
    )

    csv_path = output_path / 'docking_angle_variation_summary.csv'
    json_path = output_path / 'docking_angle_variation_summary.json'
    md_path = output_path / 'docking_angle_variation_summary.md'
    plot_crossing = output_path / 'crossing_tail90_variation_rank.png'
    plot_incident = output_path / 'incident_tail90_variation_rank.png'

    df.to_csv(csv_path, index=False)
    json_path.write_text(
        json.dumps(rows, indent=2, ensure_ascii=False, default=str),
        encoding='utf-8',
    )

    success_count = int((df['status'] == 'success').sum())
    failed_count = int((df['status'] == 'failed').sum())
    top_rows = df[
        [
            'task_name',
            'status',
            'crossing_tail90_variation_deg',
            'incident_tail90_variation_deg',
            'crossing_mean_deg',
            'incident_mean_deg',
        ]
    ].head(20)

    lines = [
        '# Batch Docking Angle Variation Summary',
        '',
        f'- Total tasks: {len(df)}',
        f'- Successful: {success_count}',
        f'- Failed: {failed_count}',
        f"- Tasks with crossing variation data: {int(df['crossing_tail90_variation_deg'].notna().sum())}",
        f"- Tasks with incident variation data: {int(df['incident_tail90_variation_deg'].notna().sum())}",
        '',
        '## Most Variable Tasks',
        '',
        '| Task | Status | Crossing Tail90 Variation (deg) | Incident Tail90 Variation (deg) | Crossing Mean (deg) | Incident Mean (deg) |',
        '|---|---|---:|---:|---:|---:|',
    ]

    for _, row in top_rows.iterrows():
        crossing_tail = '' if pd.isna(row['crossing_tail90_variation_deg']) else f"{row['crossing_tail90_variation_deg']:.2f}"
        incident_tail = '' if pd.isna(row['incident_tail90_variation_deg']) else f"{row['incident_tail90_variation_deg']:.2f}"
        crossing_mean = '' if pd.isna(row['crossing_mean_deg']) else f"{row['crossing_mean_deg']:.2f}"
        incident_mean = '' if pd.isna(row['incident_mean_deg']) else f"{row['incident_mean_deg']:.2f}"
        lines.append(
            f"| {row['task_name']} | {row['status']} | {crossing_tail} | {incident_tail} | {crossing_mean} | {incident_mean} |"
        )
    lines.append('')
    md_path.write_text('\n'.join(lines), encoding='utf-8')

    def _plot_rank(metric_name: str, title: str, output_file: Path, color: str):
        plot_df = df[(df['status'] == 'success') & df[metric_name].notna()].head(30)
        if plot_df.empty:
            return
        fig, ax = plt.subplots(figsize=(12, max(6, len(plot_df) * 0.3)))
        plot_df = plot_df.sort_values(by=metric_name, ascending=True)
        ax.barh(plot_df['task_name'], plot_df[metric_name], color=color, alpha=0.9)
        ax.set_xlabel('Angle Variation (deg)')
        ax.set_ylabel('Task')
        ax.set_title(title)
        ax.grid(True, axis='x', alpha=0.2)
        fig.tight_layout()
        fig.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close(fig)

    _plot_rank(
        'crossing_tail90_variation_deg',
        'Crossing Angle Tail90 Variation Ranking',
        plot_crossing,
        '#1f4e79',
    )
    _plot_rank(
        'incident_tail90_variation_deg',
        'Incident Angle Tail90 Variation Ranking',
        plot_incident,
        '#7a3b69',
    )

    return {
        'csv': str(csv_path),
        'json': str(json_path),
        'markdown': str(md_path),
        'crossing_plot': str(plot_crossing),
        'incident_plot': str(plot_incident),
    }


def handle_preprocess(args, logger):
    """Handle batch preprocessing action"""
    base_dir = Path(args.base_dir)
    if not base_dir.exists():
        logger.error(f"Base directory not found: {args.base_dir}")
        return 1

    logger.info("=" * 60)
    logger.info("IMN Batch Preprocess")
    logger.info("=" * 60)
    logger.info(f"Base directory: {args.base_dir}")
    logger.info(f"Workers: {args.workers}")
    logger.info(f"PBC method: {args.method}")
    logger.info("")

    output_root = args.output_dir or str(Path('./output') / 'preprocess_batch' / base_dir.name)

    try:
        report = _discover_execution_report(
            base_dir,
            required_files=['topology', 'trajectory'],
        )
    except Exception as exc:
        logger.error(f"Task discovery failed: {exc}")
        return 1

    if report.num_valid == 0:
        logger.error('No valid tasks found for batch preprocessing')
        return 1

    if report.num_invalid or report.num_ambiguous:
        logger.info(f"Skipped {report.num_invalid} invalid and {report.num_ambiguous} ambiguous tasks")

    pipeline = PreprocessQualityPipeline(
        method=args.method,
        rmsd_selection=args.rmsd_selection,
    )
    executor = BatchExecutor(max_workers=args.workers)

    try:
        contexts = executor.execute_pipeline(
            report,
            pipeline,
            show_progress=not args.verbose,
            output_base_dir=output_root,
        )
    except Exception as exc:
        logger.error(f"Batch preprocess failed: {exc}")
        return 1

    results = _summarize_context_results(contexts)
    results['output_directory'] = output_root
    results['summary_files'] = _write_preprocess_quality_summary(results, output_root)

    logger.info("")
    logger.info("Batch preprocess summary")
    logger.info(f"Total tasks: {results['total_tasks']}")
    logger.info(f"Successful: {results['successful']}")
    logger.info(f"Failed: {results['failed']}")
    logger.info(f"Output directory: {results['output_directory']}")
    logger.info(f"RMSD variation summary CSV: {results['summary_files']['csv']}")
    logger.info(f"Problematic tasks CSV: {results['summary_files']['problematic_csv']}")
    logger.info(f"Variation ranking plot: {results['summary_files']['plot']}")

    return 0 if results['failed'] == 0 else 1


def _run_quality_task(task, output_root: Path, stride: int):
    """Run quality assessment for a single discovered task."""
    pipeline = QualityAssessmentPipeline()
    task_output_dir = output_root / task.task_id
    task_output_dir.mkdir(parents=True, exist_ok=True)

    results = pipeline.run_comprehensive_assessment(
        trajectory=task.input_files.trajectory_path,
        topology=task.input_files.topology_path,
        output_dir=str(task_output_dir),
        validation_stride=stride,
    )

    markdown_report = task_output_dir / 'quality_report.md'
    pipeline.generate_quality_report(results=results, output_file=str(markdown_report), format='markdown')

    json_report = task_output_dir / 'quality_report.json'
    with open(json_report, 'w') as handle:
        json.dump(results, handle, indent=2, default=str)

    return {
        'task_id': task.task_id,
        'task_root': task.task_root,
        'trajectory': task.input_files.trajectory_path,
        'topology': task.input_files.topology_path,
        'output_dir': str(task_output_dir),
        'overall_grade': results['overall_grade'],
        'is_qualified': results['is_qualified'],
        'status': 'success',
        'reports': {
            'markdown': str(markdown_report),
            'json': str(json_report),
        },
    }


def handle_quality(args, logger):
    """Handle batch quality assessment action"""
    if not Path(args.base_dir).exists():
        logger.error(f"Base directory not found: {args.base_dir}")
        return 1

    logger.info("=" * 60)
    logger.info("IMN Batch Quality Assessment")
    logger.info("=" * 60)
    logger.info(f"Base directory: {args.base_dir}")
    logger.info(f"Workers: {args.workers}")
    logger.info("")

    output_root = Path(args.output_dir or Path('./output') / 'quality_batch' / Path(args.base_dir).name)
    output_root.mkdir(parents=True, exist_ok=True)

    logger.info('Discovering processed MD tasks...')
    report = discover_tasks(Path(args.base_dir), required_files=['topology', 'trajectory'])

    if report.num_valid == 0:
        logger.error('No valid tasks found for batch quality assessment')
        return 1

    logger.info(f'Discovered {report.num_valid} valid tasks')
    if report.num_invalid or report.num_ambiguous:
        logger.info(f'Skipped {report.num_invalid} invalid and {report.num_ambiguous} ambiguous tasks')

    task_results = []
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        future_map = {
            executor.submit(_run_quality_task, task, output_root, args.stride): task
            for task in report.valid_tasks
        }

        for future in as_completed(future_map):
            task = future_map[future]
            try:
                result = future.result()
                task_results.append(result)
                logger.info(f"[{result['overall_grade']}] {result['task_id']} {'PASS' if result['is_qualified'] else 'FAIL'}")
            except Exception as exc:
                logger.error(f'Task {task.task_id} failed: {exc}')
                task_results.append({
                    'task_id': task.task_id,
                    'task_root': task.task_root,
                    'trajectory': task.input_files.trajectory_path,
                    'topology': task.input_files.topology_path,
                    'status': 'failed',
                    'error': str(exc),
                })

    successful = [item for item in task_results if item['status'] == 'success']
    failed = [item for item in task_results if item['status'] == 'failed']
    qualified = [item for item in successful if item['is_qualified']]

    summary = {
        'base_dir': str(Path(args.base_dir).resolve()),
        'output_dir': str(output_root.resolve()),
        'total_discovered': report.total_tasks,
        'valid_tasks': report.num_valid,
        'invalid_tasks': report.num_invalid,
        'ambiguous_tasks': report.num_ambiguous,
        'successful': len(successful),
        'failed': len(failed),
        'qualified': len(qualified),
        'results': sorted(task_results, key=lambda item: item['task_id']),
    }

    if args.summary:
        summary_path = Path(args.summary)
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        with open(summary_path, 'w') as handle:
            json.dump(summary, handle, indent=2, default=str)
        logger.info(f'Summary JSON: {summary_path}')

    logger.info('')
    logger.info('Batch quality summary')
    logger.info(f'Valid tasks: {report.num_valid}')
    logger.info(f'Successful: {len(successful)}')
    logger.info(f'Failed: {len(failed)}')
    logger.info(f'Qualified: {len(qualified)}')
    logger.info(f'Output directory: {output_root}')

    return 0 if not failed else 1


def handle_contact(args, logger):
    """Handle batch contact analysis action."""
    base_dir = Path(args.base_dir)
    if not base_dir.exists():
        logger.error(f"Base directory not found: {args.base_dir}")
        return 1

    logger.info('=' * 60)
    logger.info('IMN Batch Contact Analysis')
    logger.info('=' * 60)
    logger.info(f'Base directory: {args.base_dir}')
    logger.info(f'Workers: {args.workers}')
    logger.info(f'Cutoff: {args.cutoff} A')
    logger.info(f'Stride: {args.stride}')
    logger.info(f'Min frequency: {args.min_frequency}')
    logger.info('')

    output_root = args.output_dir or str(Path('./output') / 'contact_batch' / base_dir.name)

    try:
        report = _discover_execution_report(
            base_dir,
            required_files=['topology', 'trajectory', 'structure'],
        )
    except Exception as exc:
        logger.error(f'Task discovery failed: {exc}')
        return 1

    if report.num_valid == 0:
        logger.error('No valid tasks found for batch contact analysis')
        return 1

    if report.num_invalid or report.num_ambiguous:
        logger.info(f'Skipped {report.num_invalid} invalid and {report.num_ambiguous} ambiguous tasks')

    pipeline = ContactFrequencyPipeline(
        cutoff=args.cutoff,
        stride=args.stride,
        min_frequency=args.min_frequency,
    )
    executor = BatchExecutor(max_workers=args.workers)

    try:
        contexts = executor.execute_pipeline(
            report,
            pipeline,
            show_progress=not args.verbose,
            output_base_dir=output_root,
        )
    except Exception as exc:
        logger.error(f'Batch contact analysis failed: {exc}')
        return 1

    results = _summarize_context_results(contexts)
    results['output_directory'] = output_root
    results['parameters'] = {
        'cutoff': args.cutoff,
        'stride': args.stride,
        'min_frequency': args.min_frequency,
    }
    for item, context in zip(results['results'], contexts):
        contact_result = context.results.get('contact_annotation', {})
        item['contact_report'] = contact_result.get('contact_report_file')

    if args.summary:
        summary_path = Path(args.summary)
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        with open(summary_path, 'w') as handle:
            json.dump(results, handle, indent=2, default=str)
        logger.info(f'Summary JSON: {summary_path}')

    logger.info('')
    logger.info('Batch contact summary')
    logger.info(f"Total tasks: {results['total_tasks']}")
    logger.info(f"Successful: {results['successful']}")
    logger.info(f"Failed: {results['failed']}")
    logger.info(f"Output directory: {results['output_directory']}")

    return 0 if results['failed'] == 0 else 1


def handle_hbond(args, logger):
    """Handle batch hydrogen-bond interaction analysis."""
    base_dir = Path(args.base_dir)
    if not base_dir.exists():
        logger.error(f"Base directory not found: {args.base_dir}")
        return 1

    logger.info('=' * 60)
    logger.info('IMN Batch Hydrogen-Bond Interaction Analysis')
    logger.info('=' * 60)
    logger.info(f'Base directory: {args.base_dir}')
    logger.info(f'Workers: {args.workers}')
    logger.info(f'Stride: {args.stride}')
    logger.info(f'Distance cutoff: {args.distance_cutoff} A')
    logger.info(f'Angle cutoff: {args.angle_cutoff} deg')
    logger.info('')

    output_root = args.output_dir or str(Path('./output') / 'hbond_batch' / base_dir.name)

    try:
        report = _discover_execution_report(
            base_dir,
            required_files=['topology', 'trajectory', 'structure'],
        )
    except Exception as exc:
        logger.error(f'Task discovery failed: {exc}')
        return 1

    if report.num_valid == 0:
        logger.error('No valid tasks found for batch hydrogen-bond analysis')
        return 1

    if report.num_invalid or report.num_ambiguous:
        logger.info(f'Skipped {report.num_invalid} invalid and {report.num_ambiguous} ambiguous tasks')

    pipeline = HydrogenBondInteractionPipeline(
        stride=args.stride,
        distance_cutoff=args.distance_cutoff,
        angle_cutoff=args.angle_cutoff,
    )
    executor = BatchExecutor(max_workers=args.workers)

    try:
        contexts = executor.execute_pipeline(
            report,
            pipeline,
            show_progress=not args.verbose,
            output_base_dir=output_root,
        )
    except Exception as exc:
        logger.error(f'Batch hydrogen-bond analysis failed: {exc}')
        return 1

    results = _summarize_context_results(contexts)
    results['output_directory'] = output_root
    results['parameters'] = {
        'stride': args.stride,
        'distance_cutoff': args.distance_cutoff,
        'angle_cutoff': args.angle_cutoff,
    }
    for item, context in zip(results['results'], contexts):
        hbond_result = context.results.get('hbond_annotation', {})
        item['hbond_report'] = hbond_result.get('hbond_report_file')

    if args.summary:
        summary_path = Path(args.summary)
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        with open(summary_path, 'w') as handle:
            json.dump(results, handle, indent=2, default=str)
        logger.info(f'Summary JSON: {summary_path}')

    logger.info('')
    logger.info('Batch hydrogen-bond summary')
    logger.info(f"Total tasks: {results['total_tasks']}")
    logger.info(f"Successful: {results['successful']}")
    logger.info(f"Failed: {results['failed']}")
    logger.info(f"Output directory: {results['output_directory']}")

    return 0 if results['failed'] == 0 else 1


def handle_saltbridge(args, logger):
    """Handle batch salt-bridge interaction analysis."""
    base_dir = Path(args.base_dir)
    if not base_dir.exists():
        logger.error(f"Base directory not found: {args.base_dir}")
        return 1

    logger.info('=' * 60)
    logger.info('IMN Batch Salt-Bridge Interaction Analysis')
    logger.info('=' * 60)
    logger.info(f'Base directory: {args.base_dir}')
    logger.info(f'Workers: {args.workers}')
    logger.info(f'Stride: {args.stride}')
    logger.info(f'Distance cutoff: {args.distance_cutoff} A')
    logger.info('')

    output_root = args.output_dir or str(Path('./output') / 'saltbridge_batch' / base_dir.name)

    try:
        report = _discover_execution_report(
            base_dir,
            required_files=['topology', 'trajectory', 'structure'],
        )
    except Exception as exc:
        logger.error(f'Task discovery failed: {exc}')
        return 1

    if report.num_valid == 0:
        logger.error('No valid tasks found for batch salt-bridge analysis')
        return 1

    if report.num_invalid or report.num_ambiguous:
        logger.info(f'Skipped {report.num_invalid} invalid and {report.num_ambiguous} ambiguous tasks')

    pipeline = SaltBridgeInteractionPipeline(
        stride=args.stride,
        distance_cutoff=args.distance_cutoff,
    )
    executor = BatchExecutor(max_workers=args.workers)

    try:
        contexts = executor.execute_pipeline(
            report,
            pipeline,
            show_progress=not args.verbose,
            output_base_dir=output_root,
        )
    except Exception as exc:
        logger.error(f'Batch salt-bridge analysis failed: {exc}')
        return 1

    results = _summarize_context_results(contexts)
    results['output_directory'] = output_root
    results['parameters'] = {'stride': args.stride, 'distance_cutoff': args.distance_cutoff}
    for item, context in zip(results['results'], contexts):
        pair_result = context.results.get('salt_bridge_annotation', {})
        item['salt_bridge_report'] = pair_result.get('salt_bridge_report_file')

    if args.summary:
        summary_path = Path(args.summary)
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        with open(summary_path, 'w') as handle:
            json.dump(results, handle, indent=2, default=str)
        logger.info(f'Summary JSON: {summary_path}')

    logger.info('')
    logger.info('Batch salt-bridge summary')
    logger.info(f"Total tasks: {results['total_tasks']}")
    logger.info(f"Successful: {results['successful']}")
    logger.info(f"Failed: {results['failed']}")
    logger.info(f"Output directory: {results['output_directory']}")
    return 0 if results['failed'] == 0 else 1


def handle_hydrophobic(args, logger):
    """Handle batch hydrophobic interaction analysis."""
    base_dir = Path(args.base_dir)
    if not base_dir.exists():
        logger.error(f"Base directory not found: {args.base_dir}")
        return 1

    logger.info('=' * 60)
    logger.info('IMN Batch Hydrophobic Interaction Analysis')
    logger.info('=' * 60)
    logger.info(f'Base directory: {args.base_dir}')
    logger.info(f'Workers: {args.workers}')
    logger.info(f'Stride: {args.stride}')
    logger.info(f'Distance cutoff: {args.distance_cutoff} A')
    logger.info('')

    output_root = args.output_dir or str(Path('./output') / 'hydrophobic_batch' / base_dir.name)

    try:
        report = _discover_execution_report(
            base_dir,
            required_files=['topology', 'trajectory', 'structure'],
        )
    except Exception as exc:
        logger.error(f'Task discovery failed: {exc}')
        return 1

    if report.num_valid == 0:
        logger.error('No valid tasks found for batch hydrophobic analysis')
        return 1

    if report.num_invalid or report.num_ambiguous:
        logger.info(f'Skipped {report.num_invalid} invalid and {report.num_ambiguous} ambiguous tasks')

    pipeline = HydrophobicInteractionPipeline(
        stride=args.stride,
        distance_cutoff=args.distance_cutoff,
    )
    executor = BatchExecutor(max_workers=args.workers)

    try:
        contexts = executor.execute_pipeline(
            report,
            pipeline,
            show_progress=not args.verbose,
            output_base_dir=output_root,
        )
    except Exception as exc:
        logger.error(f'Batch hydrophobic analysis failed: {exc}')
        return 1

    results = _summarize_context_results(contexts)
    results['output_directory'] = output_root
    results['parameters'] = {'stride': args.stride, 'distance_cutoff': args.distance_cutoff}
    for item, context in zip(results['results'], contexts):
        pair_result = context.results.get('hydrophobic_annotation', {})
        item['hydrophobic_report'] = pair_result.get('hydrophobic_report_file')

    if args.summary:
        summary_path = Path(args.summary)
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        with open(summary_path, 'w') as handle:
            json.dump(results, handle, indent=2, default=str)
        logger.info(f'Summary JSON: {summary_path}')

    logger.info('')
    logger.info('Batch hydrophobic summary')
    logger.info(f"Total tasks: {results['total_tasks']}")
    logger.info(f"Successful: {results['successful']}")
    logger.info(f"Failed: {results['failed']}")
    logger.info(f"Output directory: {results['output_directory']}")
    return 0 if results['failed'] == 0 else 1


def handle_pipi(args, logger):
    """Handle batch pi-pi interaction analysis."""
    base_dir = Path(args.base_dir)
    if not base_dir.exists():
        logger.error(f"Base directory not found: {args.base_dir}")
        return 1

    logger.info('=' * 60)
    logger.info('IMN Batch Pi-Pi Interaction Analysis')
    logger.info('=' * 60)
    logger.info(f'Base directory: {args.base_dir}')
    logger.info(f'Workers: {args.workers}')
    logger.info(f'Stride: {args.stride}')
    logger.info(f'Distance cutoff: {args.distance_cutoff} A')
    logger.info('')

    output_root = args.output_dir or str(Path('./output') / 'pipi_batch' / base_dir.name)

    try:
        report = _discover_execution_report(
            base_dir,
            required_files=['topology', 'trajectory', 'structure'],
        )
    except Exception as exc:
        logger.error(f'Task discovery failed: {exc}')
        return 1

    if report.num_valid == 0:
        logger.error('No valid tasks found for batch pi-pi analysis')
        return 1

    if report.num_invalid or report.num_ambiguous:
        logger.info(f'Skipped {report.num_invalid} invalid and {report.num_ambiguous} ambiguous tasks')

    pipeline = PiStackingInteractionPipeline(
        stride=args.stride,
        distance_cutoff=args.distance_cutoff,
    )
    executor = BatchExecutor(max_workers=args.workers)

    try:
        contexts = executor.execute_pipeline(
            report,
            pipeline,
            show_progress=not args.verbose,
            output_base_dir=output_root,
        )
    except Exception as exc:
        logger.error(f'Batch pi-pi analysis failed: {exc}')
        return 1

    results = _summarize_context_results(contexts)
    results['output_directory'] = output_root
    results['parameters'] = {'stride': args.stride, 'distance_cutoff': args.distance_cutoff}
    for item, context in zip(results['results'], contexts):
        pair_result = context.results.get('pi_pi_annotation', {})
        item['pi_pi_report'] = pair_result.get('pi_pi_report_file')

    if args.summary:
        summary_path = Path(args.summary)
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        with open(summary_path, 'w') as handle:
            json.dump(results, handle, indent=2, default=str)
        logger.info(f'Summary JSON: {summary_path}')

    logger.info('')
    logger.info('Batch pi-pi summary')
    logger.info(f"Total tasks: {results['total_tasks']}")
    logger.info(f"Successful: {results['successful']}")
    logger.info(f"Failed: {results['failed']}")
    logger.info(f"Output directory: {results['output_directory']}")
    return 0 if results['failed'] == 0 else 1


def handle_cationpi(args, logger):
    """Handle batch cation-pi interaction analysis."""
    base_dir = Path(args.base_dir)
    if not base_dir.exists():
        logger.error(f"Base directory not found: {args.base_dir}")
        return 1

    logger.info('=' * 60)
    logger.info('IMN Batch Cation-Pi Interaction Analysis')
    logger.info('=' * 60)
    logger.info(f'Base directory: {args.base_dir}')
    logger.info(f'Workers: {args.workers}')
    logger.info(f'Stride: {args.stride}')
    logger.info(f'Distance cutoff: {args.distance_cutoff} A')
    logger.info(f'Normal angle cutoff: {args.normal_angle_cutoff} deg')
    logger.info('')

    output_root = args.output_dir or str(Path('./output') / 'cationpi_batch' / base_dir.name)

    try:
        report = _discover_execution_report(
            base_dir,
            required_files=['topology', 'trajectory', 'structure'],
        )
    except Exception as exc:
        logger.error(f'Task discovery failed: {exc}')
        return 1

    if report.num_valid == 0:
        logger.error('No valid tasks found for batch cation-pi analysis')
        return 1

    if report.num_invalid or report.num_ambiguous:
        logger.info(f'Skipped {report.num_invalid} invalid and {report.num_ambiguous} ambiguous tasks')

    pipeline = CationPiInteractionPipeline(
        stride=args.stride,
        distance_cutoff=args.distance_cutoff,
        normal_angle_cutoff=args.normal_angle_cutoff,
    )
    executor = BatchExecutor(max_workers=args.workers)

    try:
        contexts = executor.execute_pipeline(
            report,
            pipeline,
            show_progress=not args.verbose,
            output_base_dir=output_root,
        )
    except Exception as exc:
        logger.error(f'Batch cation-pi analysis failed: {exc}')
        return 1

    results = _summarize_context_results(contexts)
    results['output_directory'] = output_root
    results['parameters'] = {
        'stride': args.stride,
        'distance_cutoff': args.distance_cutoff,
        'normal_angle_cutoff': args.normal_angle_cutoff,
    }
    for item, context in zip(results['results'], contexts):
        pair_result = context.results.get('cation_pi_annotation', {})
        item['cation_pi_report'] = pair_result.get('cation_pi_report_file')

    if args.summary:
        summary_path = Path(args.summary)
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        with open(summary_path, 'w') as handle:
            json.dump(results, handle, indent=2, default=str)
        logger.info(f'Summary JSON: {summary_path}')

    logger.info('')
    logger.info('Batch cation-pi summary')
    logger.info(f"Total tasks: {results['total_tasks']}")
    logger.info(f"Successful: {results['successful']}")
    logger.info(f"Failed: {results['failed']}")
    logger.info(f"Output directory: {results['output_directory']}")
    return 0 if results['failed'] == 0 else 1


def handle_angle(args, logger):
    """Handle batch docking angle analysis."""
    base_dir = Path(args.base_dir)
    if not base_dir.exists():
        logger.error(f'Base directory not found: {args.base_dir}')
        return 1

    logger.info('=' * 60)
    logger.info('IMN Batch Docking Angle Analysis')
    logger.info('=' * 60)
    logger.info(f'Base directory: {args.base_dir}')
    logger.info(f'Workers: {args.workers}')
    logger.info(f'Stride: {args.stride}')
    if args.limit:
        logger.info(f'Limit: {args.limit}')
    logger.info('')

    output_root = args.output_dir or str(Path('./output') / 'angle_batch' / base_dir.name)

    try:
        report = _discover_execution_report(
            base_dir,
            required_files=['structure', 'topology', 'trajectory'],
        )
    except Exception as exc:
        logger.error(f'Task discovery failed: {exc}')
        return 1

    if report.num_valid == 0:
        logger.error('No valid tasks found for batch docking angle analysis')
        return 1

    if args.limit:
        report.valid_tasks = report.valid_tasks[:args.limit]
        report.total_tasks = len(report.valid_tasks)

    if report.num_invalid or report.num_ambiguous:
        logger.info(f'Skipped {report.num_invalid} invalid and {report.num_ambiguous} ambiguous tasks')

    pipeline = DockingAnglePipeline(
        stride=args.stride,
        auto_identify_chains=True,
        print_each_frame=args.print_each_frame,
    )
    executor = BatchExecutor(max_workers=args.workers)

    try:
        contexts = executor.execute_pipeline(
            report,
            pipeline,
            show_progress=not args.verbose,
            output_base_dir=output_root,
        )
    except Exception as exc:
        logger.error(f'Batch docking angle analysis failed: {exc}')
        return 1

    results = _summarize_context_results(contexts)
    results['output_directory'] = output_root
    results['parameters'] = {
        'stride': args.stride,
        'limit': args.limit,
        'print_each_frame': args.print_each_frame,
    }
    for item, context in zip(results['results'], contexts):
        item['docking_angles'] = context.results.get('docking_angles')
        item['angle_analysis'] = context.metadata.get('angle_analysis')

    results['summary_files'] = _write_angle_summary(results, output_root)

    if args.summary:
        summary_path = Path(args.summary)
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        with open(summary_path, 'w') as handle:
            json.dump(results, handle, indent=2, default=str)
        logger.info(f'Summary JSON: {summary_path}')

    logger.info('')
    logger.info('Batch docking angle summary')
    logger.info(f"Total tasks: {results['total_tasks']}")
    logger.info(f"Successful: {results['successful']}")
    logger.info(f"Failed: {results['failed']}")
    logger.info(f"Output directory: {results['output_directory']}")
    logger.info(f"Summary CSV: {results['summary_files']['csv']}")
    logger.info(f"Crossing variation plot: {results['summary_files']['crossing_plot']}")
    logger.info(f"Incident variation plot: {results['summary_files']['incident_plot']}")

    return 0 if results['failed'] == 0 else 1


def handle_bsa(args, logger):
    """Handle batch buried surface area analysis."""
    base_dir = Path(args.base_dir)
    if not base_dir.exists():
        logger.error(f'Base directory not found: {args.base_dir}')
        return 1

    logger.info('=' * 60)
    logger.info('IMN Batch Buried Surface Area Analysis')
    logger.info('=' * 60)
    logger.info(f'Base directory: {args.base_dir}')
    logger.info(f'Workers: {args.workers}')
    logger.info(f'Stride: {args.stride}')
    logger.info(f'Probe radius: {args.probe_radius} A')
    logger.info(f'Time unit: {args.time_unit}')
    logger.info('')

    output_root = args.output_dir or str(Path('./output') / 'bsa_batch' / base_dir.name)

    try:
        report = _discover_execution_report(
            base_dir,
            required_files=['topology', 'trajectory', 'structure'],
        )
    except Exception as exc:
        logger.error(f'Task discovery failed: {exc}')
        return 1

    if report.num_valid == 0:
        logger.error('No valid tasks found for batch BSA analysis')
        return 1

    if report.num_invalid or report.num_ambiguous:
        logger.info(f'Skipped {report.num_invalid} invalid and {report.num_ambiguous} ambiguous tasks')

    pipeline = BSAPipeline(
        probe_radius=args.probe_radius,
        stride=args.stride,
        time_unit=args.time_unit,
        auto_identify_chains=True,
    )
    executor = BatchExecutor(max_workers=args.workers)

    try:
        contexts = executor.execute_pipeline(
            report,
            pipeline,
            show_progress=not args.verbose,
            output_base_dir=output_root,
        )
    except Exception as exc:
        logger.error(f'Batch BSA analysis failed: {exc}')
        return 1

    results = _summarize_context_results(contexts)
    results['output_directory'] = output_root
    results['parameters'] = {
        'stride': args.stride,
        'probe_radius': args.probe_radius,
        'time_unit': args.time_unit,
    }
    for item, context in zip(results['results'], contexts):
        item['bsa'] = context.results.get('bsa')

    if args.summary:
        summary_path = Path(args.summary)
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        with open(summary_path, 'w') as handle:
            json.dump(results, handle, indent=2, default=str)
        logger.info(f'Summary JSON: {summary_path}')

    logger.info('')
    logger.info('Batch BSA summary')
    logger.info(f"Total tasks: {results['total_tasks']}")
    logger.info(f"Successful: {results['successful']}")
    logger.info(f"Failed: {results['failed']}")
    logger.info(f"Output directory: {results['output_directory']}")
    return 0 if results['failed'] == 0 else 1


def handle_workflow(args, logger):
    """Handle complete workflow (PBC + Quality)"""
    logger.info('=' * 60)
    logger.info('IMN Batch Workflow')
    logger.info('=' * 60)
    logger.info(f'Base directory: {args.base_dir}')
    logger.info(f'Workers: {args.workers}')
    logger.info('')

    ret = handle_preprocess(args, logger)
    if ret != 0:
        logger.error('PBC processing failed')
        return ret

    quality_args = argparse.Namespace(
        base_dir=args.output_dir or str(Path('./output') / 'preprocess_batch' / Path(args.base_dir).name),
        workers=args.workers,
        output_dir=args.output_dir,
        summary=args.summary,
        stride=args.stride,
        verbose=args.verbose,
    )

    logger.info('')
    logger.info('Step 2: Batch quality assessment')
    logger.info('-' * 60)
    ret = handle_quality(quality_args, logger)
    if ret != 0:
        logger.error('Quality assessment failed')
        return ret

    logger.info('')
    logger.info('=' * 60)
    logger.info('Workflow completed successfully!')
    logger.info('=' * 60)
    return 0


def main(argv=None):
    """Main entry point for batch subcommand"""
    parser = create_parser()
    args = parser.parse_args(argv)

    if not args.action:
        parser.print_help()
        return 1

    setup_logging(args.verbose)
    logger = logging.getLogger(__name__)

    try:
        if args.action == 'preprocess':
            return handle_preprocess(args, logger)
        if args.action == 'quality':
            return handle_quality(args, logger)
        if args.action == 'contact':
            return handle_contact(args, logger)
        if args.action == 'hbond':
            return handle_hbond(args, logger)
        if args.action == 'saltbridge':
            return handle_saltbridge(args, logger)
        if args.action == 'hydrophobic':
            return handle_hydrophobic(args, logger)
        if args.action == 'pipi':
            return handle_pipi(args, logger)
        if args.action == 'cationpi':
            return handle_cationpi(args, logger)
        if args.action == 'angle':
            return handle_angle(args, logger)
        if args.action == 'bsa':
            return handle_bsa(args, logger)
        if args.action == 'workflow':
            return handle_workflow(args, logger)

        logger.error(f'Unknown action: {args.action}')
        return 1
    except KeyboardInterrupt:
        logger.info('Interrupted by user')
        return 130


if __name__ == '__main__':
    sys.exit(main())
