#!/usr/bin/env python3
"""Monitor batch MoSAIC clustering progress

This script checks the progress of the batch MoSAIC clustering
and provides real-time statistics.

Usage:
    python monitor_mosaic_batch.py
"""

import sys
from pathlib import Path
import pandas as pd
import time

def count_completed_tasks(base_dir):
    """Count tasks with completed MoSAIC results"""
    base_path = Path(base_dir)

    total = 0
    completed = 0

    for task_dir in sorted(base_path.iterdir()):
        if not task_dir.is_dir():
            continue

        corr_matrix = task_dir / "correlation_matrix.npz"
        if not corr_matrix.exists():
            continue

        total += 1

        mosaic_summary = task_dir / "mosaic" / "mosaic_cluster_summary.csv"
        if mosaic_summary.exists():
            completed += 1

    return total, completed


def print_progress(base_dir, summary_file, log_file):
    """Print batch processing progress"""
    total, completed = count_completed_tasks(base_dir)

    print("=" * 80)
    print("MoSAIC批量处理进度监控")
    print("=" * 80)
    print(f"总任务数: {total}")
    print(f"已完成: {completed}")
    print(f"剩余: {total - completed}")
    print(f"进度: {completed}/{total} ({completed/total*100:.1f}%)")

    # Check summary file
    if Path(summary_file).exists():
        print("\n汇总文件状态:")
        df = pd.read_csv(summary_file)
        print(f"  记录数: {len(df)}")

        if 'status' in df.columns:
            success = sum(df['status'] == 'success')
            failed = sum(df['status'] == 'failed')
            print(f"  成功: {success}")
            print(f"  失败: {failed}")

            if success > 0:
                print("\n聚类统计 (成功任务):")
                successful = df[df['status'] == 'success']
                print(f"  平均簇数: {successful['n_clusters'].mean():.1f}")
                print(f"  平均有意义簇数 (size>1): {successful['n_meaningful_clusters'].mean():.1f}")
                if 'silhouette_score' in successful.columns:
                    print(f"  平均轮廓系数: {successful['silhouette_score'].mean():.3f}")

    # Check log file
    if Path(log_file).exists():
        print("\n最新日志 (最后10行):")
        print("-" * 80)
        with open(log_file, 'r') as f:
            lines = f.readlines()
            for line in lines[-10:]:
                print(line.rstrip())

    print("=" * 80)


def watch_progress(base_dir, summary_file, log_file, interval=30):
    """Continuously monitor progress"""
    print("开始监控批量处理进度...")
    print(f"刷新间隔: {interval}秒")
    print("按 Ctrl+C 退出监控\n")

    try:
        while True:
            print_progress(base_dir, summary_file, log_file)
            time.sleep(interval)
            print("\n" * 2)
    except KeyboardInterrupt:
        print("\n\n监控已停止")


if __name__ == "__main__":
    base_dir = "/home/xumy/work/development/AfterMD/output/allostery_analysis"
    summary_file = f"{base_dir}/mosaic_batch_summary.csv"
    log_file = f"{base_dir}/mosaic_batch.log"

    if len(sys.argv) > 1 and sys.argv[1] == "--watch":
        # Continuous monitoring mode
        interval = int(sys.argv[2]) if len(sys.argv) > 2 else 30
        watch_progress(base_dir, summary_file, log_file, interval)
    else:
        # Single check mode
        print_progress(base_dir, summary_file, log_file)
