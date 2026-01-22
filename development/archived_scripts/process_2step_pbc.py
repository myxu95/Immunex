#!/usr/bin/env python3
"""
2-Step PBC Processing for 1000frames Trajectories

基于1000frames_extract目录中的轨迹进行简化2步PBC校正:
1. Step 1: gmx trjconv -pbc nojump (去除周期性跳变)
2. Step 2: gmx trjconv -fit rot+trans (结构对齐)
3. RMSD评估

输入: /public/home/xmy/dataset/workspace/1000frames_extract/
输出: /public/home/xmy/dataset/workspace/pbc_1000frames_2step/

作者: MD Analysis Workflow
日期: 2025-12-08
"""

import sys
import os
import logging
import subprocess
import json
from pathlib import Path
from typing import Dict, Optional, Tuple
from datetime import datetime
from multiprocessing import Pool
import numpy as np

# ==================== 配置路径 ====================
WORKSPACE = Path("/public/home/xmy/dataset/workspace")
INPUT_BASE = WORKSPACE / "1000frames_extract"
OUTPUT_BASE = WORKSPACE / "pbc_1000frames_2step"
LOGS_DIR = OUTPUT_BASE / "logs"

# 创建输出目录
OUTPUT_BASE.mkdir(parents=True, exist_ok=True)
LOGS_DIR.mkdir(parents=True, exist_ok=True)

# ==================== 处理参数 ====================
N_PROCESSES = 2  # 并行进程数（前台运行，避免资源占用过多）
EXPECTED_FRAMES = 1001  # 期望帧数
TIMEOUT = 600  # 超时时间(秒)

# ==================== 日志配置 ====================
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
log_file = LOGS_DIR / f"pbc_2step_processing_{timestamp}.log"

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(log_file),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)


def run_gmx_command(cmd: list, input_str: str = "", timeout: int = TIMEOUT) -> Tuple[bool, str, str]:
    """执行GROMACS命令"""
    try:
        result = subprocess.run(
            cmd,
            input=input_str,
            capture_output=True,
            text=True,
            timeout=timeout,
            shell=False
        )
        return result.returncode == 0, result.stdout, result.stderr
    except subprocess.TimeoutExpired:
        return False, "", f"Command timeout after {timeout}s"
    except Exception as e:
        return False, "", str(e)


def parse_rmsd_xvg(xvg_file: Path) -> Optional[Dict]:
    """解析RMSD XVG文件并计算统计量"""
    rmsd_values = []
    try:
        with open(xvg_file, 'r') as f:
            for line in f:
                if not line.startswith(('#', '@')):
                    parts = line.split()
                    if len(parts) >= 2:
                        rmsd_nm = float(parts[1])
                        rmsd_values.append(rmsd_nm * 10.0)  # 转换为埃

        if len(rmsd_values) == 0:
            return None

        rmsd_mean = float(np.mean(rmsd_values))
        rmsd_std = float(np.std(rmsd_values))

        # 计算帧间跳变
        jumps = [abs(rmsd_values[i+1] - rmsd_values[i]) for i in range(len(rmsd_values)-1)]
        max_jump = float(max(jumps)) if jumps else 0.0
        mean_jump = float(np.mean(jumps)) if jumps else 0.0

        # 质量评估
        if max_jump <= 2.0:
            quality = "excellent"
        elif max_jump <= 3.0:
            quality = "good"
        elif max_jump <= 5.0:
            quality = "moderate"
        else:
            quality = "poor"

        return {
            'mean': rmsd_mean,
            'std': rmsd_std,
            'max_jump': max_jump,
            'mean_jump': mean_jump,
            'quality': quality,
            'num_frames': len(rmsd_values)
        }
    except Exception as e:
        logger.error(f"Failed to parse RMSD file {xvg_file}: {e}")
        return None


def process_single_task(task_name: str) -> Dict:
    """处理单个任务"""
    start_time = datetime.now()
    logger.info(f"Processing {task_name}...")

    result = {
        'task_name': task_name,
        'status': 'Failed',
        'error_message': '',
        'step1_nojump': False,
        'step2_fit': False,
        'rmsd_calculated': False,
        'rmsd_mean': 0.0,
        'rmsd_std': 0.0,
        'max_jump': 0.0,
        'quality': 'unknown',
        'processing_time_s': 0.0
    }

    try:
        # ============ 1. 查找输入文件 ============
        input_dir = INPUT_BASE / task_name
        if not input_dir.exists():
            result['error_message'] = f"Input directory not found: {input_dir}"
            logger.error(result['error_message'])
            return result

        xtc_file = input_dir / f"{task_name}_1000frames.xtc"
        tpr_file = input_dir / "md.tpr"

        if not xtc_file.exists() or not tpr_file.exists():
            result['error_message'] = f"Missing input files: {xtc_file} or {tpr_file}"
            logger.error(result['error_message'])
            return result

        # ============ 2. 创建输出目录 ============
        output_dir = OUTPUT_BASE / task_name
        output_dir.mkdir(parents=True, exist_ok=True)

        # 输出文件路径
        nojump_xtc = output_dir / f"{task_name}_nojump.xtc"
        final_xtc = output_dir / f"{task_name}_2step_processed.xtc"
        rmsd_xvg = output_dir / "rmsd_2step.xvg"
        rmsd_json = output_dir / "rmsd_stats_2step.json"

        # ============ 3. Step 1: PBC nojump ============
        logger.info(f"  [{task_name}] Step 1: PBC nojump...")
        cmd_nojump = [
            "gmx", "trjconv",
            "-s", str(tpr_file),
            "-f", str(xtc_file),
            "-o", str(nojump_xtc),
            "-pbc", "nojump"
        ]

        success, stdout, stderr = run_gmx_command(cmd_nojump, input_str="0\n")
        if not success:
            result['error_message'] = f"Step 1 failed: {stderr}"
            logger.error(f"  [{task_name}] {result['error_message']}")
            return result

        result['step1_nojump'] = True
        logger.info(f"  [{task_name}] Step 1 completed")

        # ============ 4. Step 2: Fit rot+trans ============
        logger.info(f"  [{task_name}] Step 2: Fit rot+trans...")
        cmd_fit = [
            "gmx", "trjconv",
            "-s", str(tpr_file),
            "-f", str(nojump_xtc),
            "-o", str(final_xtc),
            "-fit", "rot+trans"
        ]

        success, stdout, stderr = run_gmx_command(cmd_fit, input_str="1\n0\n")
        if not success:
            result['error_message'] = f"Step 2 failed: {stderr}"
            logger.error(f"  [{task_name}] {result['error_message']}")
            return result

        result['step2_fit'] = True
        logger.info(f"  [{task_name}] Step 2 completed")

        # ============ 5. RMSD计算 ============
        logger.info(f"  [{task_name}] Calculating RMSD...")
        cmd_rmsd = [
            "gmx", "rms",
            "-s", str(tpr_file),
            "-f", str(final_xtc),
            "-o", str(rmsd_xvg)
        ]

        success, stdout, stderr = run_gmx_command(cmd_rmsd, input_str="4\n4\n")
        if not success:
            result['error_message'] = f"RMSD calculation failed: {stderr}"
            logger.warning(f"  [{task_name}] {result['error_message']}")
        else:
            result['rmsd_calculated'] = True

            # 解析RMSD统计
            rmsd_stats = parse_rmsd_xvg(rmsd_xvg)
            if rmsd_stats:
                result['rmsd_mean'] = rmsd_stats['mean']
                result['rmsd_std'] = rmsd_stats['std']
                result['max_jump'] = rmsd_stats['max_jump']
                result['quality'] = rmsd_stats['quality']

                # 保存RMSD统计到JSON
                with open(rmsd_json, 'w') as f:
                    json.dump(rmsd_stats, f, indent=2)

                logger.info(f"  [{task_name}] RMSD: {rmsd_stats['mean']:.2f}±{rmsd_stats['std']:.2f}Å, "
                           f"MaxJump: {rmsd_stats['max_jump']:.2f}Å, Quality: {rmsd_stats['quality']}")

        # ============ 6. 清理中间文件 ============
        if nojump_xtc.exists():
            nojump_xtc.unlink()
            logger.info(f"  [{task_name}] Cleaned up intermediate file: {nojump_xtc.name}")

        # ============ 7. 成功完成 ============
        result['status'] = 'Success'
        processing_time = (datetime.now() - start_time).total_seconds()
        result['processing_time_s'] = processing_time

        logger.info(f"  [{task_name}] ✓ Completed in {processing_time:.1f}s")

    except Exception as e:
        result['error_message'] = f"Unexpected error: {str(e)}"
        logger.error(f"  [{task_name}] {result['error_message']}")

    return result


def main():
    """主处理流程"""
    logger.info("=" * 80)
    logger.info("2-Step PBC Processing for 1000frames Trajectories")
    logger.info("=" * 80)
    logger.info(f"Input directory: {INPUT_BASE}")
    logger.info(f"Output directory: {OUTPUT_BASE}")
    logger.info(f"Parallel processes: {N_PROCESSES}")
    logger.info("")

    # ============ 1. 获取任务列表 ============
    task_dirs = [d for d in INPUT_BASE.iterdir() if d.is_dir()]
    task_names = [d.name for d in task_dirs]
    task_names.sort()

    logger.info(f"Found {len(task_names)} tasks to process")
    logger.info("")

    if len(task_names) == 0:
        logger.error("No tasks found!")
        return

    # ============ 2. 并行处理 ============
    start_time = datetime.now()

    with Pool(processes=N_PROCESSES) as pool:
        results = pool.map(process_single_task, task_names)

    total_time = (datetime.now() - start_time).total_seconds()

    # ============ 3. 汇总结果 ============
    logger.info("")
    logger.info("=" * 80)
    logger.info("Processing Summary")
    logger.info("=" * 80)

    success_count = sum(1 for r in results if r['status'] == 'Success')
    failed_count = len(results) - success_count

    # 质量统计
    quality_counts = {'excellent': 0, 'good': 0, 'moderate': 0, 'poor': 0, 'unknown': 0}
    for r in results:
        if r['status'] == 'Success':
            quality_counts[r['quality']] += 1

    logger.info(f"Total tasks: {len(results)}")
    logger.info(f"Success: {success_count} ({success_count/len(results)*100:.1f}%)")
    logger.info(f"Failed: {failed_count}")
    logger.info(f"Total time: {total_time/60:.1f} minutes")
    logger.info("")

    logger.info("Quality Distribution:")
    logger.info(f"  Excellent: {quality_counts['excellent']} ({quality_counts['excellent']/success_count*100:.1f}%)")
    logger.info(f"  Good:      {quality_counts['good']} ({quality_counts['good']/success_count*100:.1f}%)")
    logger.info(f"  Moderate:  {quality_counts['moderate']} ({quality_counts['moderate']/success_count*100:.1f}%)")
    logger.info(f"  Poor:      {quality_counts['poor']} ({quality_counts['poor']/success_count*100:.1f}%)")
    logger.info("")

    # ============ 4. 保存CSV结果 ============
    csv_file = LOGS_DIR / f"pbc_2step_summary_{timestamp}.csv"

    with open(csv_file, 'w') as f:
        # 写入表头
        f.write("task_name,status,rmsd_mean,rmsd_std,max_jump,quality,processing_time_s,error_message\n")

        # 写入数据
        for r in results:
            f.write(f"{r['task_name']},{r['status']},{r['rmsd_mean']:.6f},{r['rmsd_std']:.6f},"
                   f"{r['max_jump']:.6f},{r['quality']},{r['processing_time_s']:.2f},"
                   f"\"{r['error_message']}\"\n")

    logger.info(f"Results saved to: {csv_file}")
    logger.info("")

    # ============ 5. 失败任务列表 ============
    if failed_count > 0:
        logger.info("Failed tasks:")
        for r in results:
            if r['status'] != 'Success':
                logger.info(f"  - {r['task_name']}: {r['error_message']}")
        logger.info("")

    logger.info("=" * 80)
    logger.info("Processing completed!")
    logger.info("=" * 80)


if __name__ == "__main__":
    main()
