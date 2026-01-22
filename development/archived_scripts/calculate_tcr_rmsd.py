#!/usr/bin/env python3
"""
Calculate TCR (chain D, E) RMSD with pHLA (chain A, B, C) as fit group

用pHLA作为参考进行结构对齐，然后计算TCR的RMSD，
这样可以分析TCR相对于pHLA的运动和动力学特征。

输入: pbc_1000frames_2step/{task}/{task}_2step_processed.xtc
输出: pbc_1000frames_2step/{task}/tcr_rmsd.xvg
      pbc_1000frames_2step/{task}/tcr_rmsd_stats.json

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
INPUT_BASE = WORKSPACE / "pbc_1000frames_2step"
EXTRACT_BASE = WORKSPACE / "1000frames_extract"  # 获取tpr和pdb文件
LOGS_DIR = INPUT_BASE / "logs"

# ==================== 处理参数 ====================
N_PROCESSES = 4  # 并行进程数
TIMEOUT = 600  # 超时时间(秒)

# ==================== 日志配置 ====================
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
log_file = LOGS_DIR / f"tcr_rmsd_calculation_{timestamp}.log"

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


def create_custom_index(tpr_file: Path, pdb_file: Path, output_ndx: Path) -> bool:
    """使用MDAnalysis创建自定义索引文件，包含pHLA和TCR组"""
    try:
        import MDAnalysis as mda

        # 使用standardized PDB文件读取结构
        u = mda.Universe(str(tpr_file), str(pdb_file))

        # 选择pHLA: chain A + B + C (HLA-α, HLA-β, Peptide)
        phla = u.select_atoms("protein and (chainID A or chainID B or chainID C)")

        # 选择TCR: chain D + E (TCR-α, TCR-β)
        tcr = u.select_atoms("protein and (chainID D or chainID E)")

        if len(phla) == 0:
            logger.error(f"  pHLA selection is empty!")
            return False

        if len(tcr) == 0:
            logger.error(f"  TCR selection is empty!")
            return False

        # 写入索引文件
        with open(output_ndx, 'w') as f:
            # pHLA组
            f.write("[ pHLA ]\n")
            atoms = phla.indices + 1  # GROMACS uses 1-based indexing
            for i in range(0, len(atoms), 15):
                line = " ".join(f"{atom:6d}" for atom in atoms[i:i+15])
                f.write(line + "\n")

            # TCR组
            f.write("\n[ TCR ]\n")
            atoms = tcr.indices + 1
            for i in range(0, len(atoms), 15):
                line = " ".join(f"{atom:6d}" for atom in atoms[i:i+15])
                f.write(line + "\n")

        logger.debug(f"  Created index: pHLA={len(phla)} atoms, TCR={len(tcr)} atoms")
        return True

    except Exception as e:
        logger.error(f"  Error creating index: {e}")
        return False


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
        rmsd_min = float(np.min(rmsd_values))
        rmsd_max = float(np.max(rmsd_values))

        # 计算帧间跳变
        jumps = [abs(rmsd_values[i+1] - rmsd_values[i]) for i in range(len(rmsd_values)-1)]
        max_jump = float(max(jumps)) if jumps else 0.0
        mean_jump = float(np.mean(jumps)) if jumps else 0.0

        return {
            'mean': rmsd_mean,
            'std': rmsd_std,
            'min': rmsd_min,
            'max': rmsd_max,
            'max_jump': max_jump,
            'mean_jump': mean_jump,
            'num_frames': len(rmsd_values)
        }
    except Exception as e:
        logger.error(f"Failed to parse RMSD file {xvg_file}: {e}")
        return None


def process_single_task(task_name: str) -> Dict:
    """处理单个任务：计算TCR相对pHLA的RMSD"""
    start_time = datetime.now()
    logger.info(f"Processing {task_name}...")

    result = {
        'task_name': task_name,
        'status': 'Failed',
        'error_message': '',
        'index_created': False,
        'rmsd_calculated': False,
        'tcr_rmsd_mean': 0.0,
        'tcr_rmsd_std': 0.0,
        'tcr_rmsd_min': 0.0,
        'tcr_rmsd_max': 0.0,
        'tcr_max_jump': 0.0,
        'processing_time_s': 0.0
    }

    try:
        # ============ 1. 查找输入文件 ============
        task_dir = INPUT_BASE / task_name
        extract_dir = EXTRACT_BASE / task_name

        if not task_dir.exists():
            result['error_message'] = f"Task directory not found: {task_dir}"
            logger.error(result['error_message'])
            return result

        # 处理后的轨迹
        processed_xtc = task_dir / f"{task_name}_2step_processed.xtc"

        # TPR文件从extract目录获取
        tpr_file = extract_dir / "md.tpr"

        # Standardized PDB文件
        pdb_file = extract_dir / f"{task_name}_standardized.pdb"

        if not processed_xtc.exists():
            result['error_message'] = f"Processed XTC not found: {processed_xtc}"
            logger.error(result['error_message'])
            return result

        if not tpr_file.exists():
            result['error_message'] = f"TPR file not found: {tpr_file}"
            logger.error(result['error_message'])
            return result

        if not pdb_file.exists():
            result['error_message'] = f"Standardized PDB not found: {pdb_file}"
            logger.error(result['error_message'])
            return result

        # ============ 2. 创建自定义索引文件 ============
        logger.info(f"  [{task_name}] Creating custom index (pHLA, TCR)...")
        custom_ndx = task_dir / "phla_tcr.ndx"

        if not create_custom_index(tpr_file, pdb_file, custom_ndx):
            result['error_message'] = "Failed to create custom index"
            logger.error(f"  [{task_name}] {result['error_message']}")
            return result

        result['index_created'] = True

        # ============ 3. 计算TCR RMSD (fit on pHLA) ============
        logger.info(f"  [{task_name}] Calculating TCR RMSD (fit on pHLA)...")

        tcr_rmsd_xvg = task_dir / "tcr_rmsd.xvg"

        # gmx rms命令:
        # -s: TPR结构文件
        # -f: 轨迹文件
        # -n: 自定义索引文件
        # -o: 输出RMSD文件
        # 第一个选择: pHLA (用于fit)
        # 第二个选择: TCR (用于计算RMSD)

        cmd_rmsd = [
            "gmx", "rms",
            "-s", str(tpr_file),
            "-f", str(processed_xtc),
            "-n", str(custom_ndx),
            "-o", str(tcr_rmsd_xvg)
        ]

        # 选择: pHLA for fit, TCR for RMSD
        # 使用组名而不是组号
        success, stdout, stderr = run_gmx_command(cmd_rmsd, input_str="pHLA\nTCR\n")

        if not success:
            result['error_message'] = f"RMSD calculation failed: {stderr}"
            logger.error(f"  [{task_name}] {result['error_message']}")
            return result

        result['rmsd_calculated'] = True

        # ============ 4. 解析RMSD统计 ============
        rmsd_stats = parse_rmsd_xvg(tcr_rmsd_xvg)
        if rmsd_stats:
            result['tcr_rmsd_mean'] = rmsd_stats['mean']
            result['tcr_rmsd_std'] = rmsd_stats['std']
            result['tcr_rmsd_min'] = rmsd_stats['min']
            result['tcr_rmsd_max'] = rmsd_stats['max']
            result['tcr_max_jump'] = rmsd_stats['max_jump']

            # 保存统计到JSON
            json_file = task_dir / "tcr_rmsd_stats.json"
            with open(json_file, 'w') as f:
                json.dump(rmsd_stats, f, indent=2)

            logger.info(f"  [{task_name}] TCR RMSD: {rmsd_stats['mean']:.2f}±{rmsd_stats['std']:.2f}Å "
                       f"(range: {rmsd_stats['min']:.2f}-{rmsd_stats['max']:.2f}Å), "
                       f"MaxJump: {rmsd_stats['max_jump']:.2f}Å")

        # ============ 5. 成功完成 ============
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
    logger.info("TCR RMSD Calculation (fit on pHLA)")
    logger.info("=" * 80)
    logger.info(f"Input directory: {INPUT_BASE}")
    logger.info(f"Parallel processes: {N_PROCESSES}")
    logger.info("")

    # ============ 1. 获取任务列表 ============
    task_dirs = [d for d in INPUT_BASE.iterdir()
                 if d.is_dir() and not d.name.startswith('logs')]
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

    logger.info(f"Total tasks: {len(results)}")
    logger.info(f"Success: {success_count} ({success_count/len(results)*100:.1f}%)")
    logger.info(f"Failed: {failed_count}")
    logger.info(f"Total time: {total_time/60:.1f} minutes")
    logger.info("")

    # 统计TCR RMSD
    if success_count > 0:
        tcr_rmsd_values = [r['tcr_rmsd_mean'] for r in results if r['status'] == 'Success']
        logger.info("TCR RMSD Statistics:")
        logger.info(f"  Mean: {np.mean(tcr_rmsd_values):.2f} ± {np.std(tcr_rmsd_values):.2f} Å")
        logger.info(f"  Range: {np.min(tcr_rmsd_values):.2f} - {np.max(tcr_rmsd_values):.2f} Å")
        logger.info("")

    # ============ 4. 保存CSV结果 ============
    csv_file = LOGS_DIR / f"tcr_rmsd_summary_{timestamp}.csv"

    with open(csv_file, 'w') as f:
        # 写入表头
        f.write("task_name,status,tcr_rmsd_mean,tcr_rmsd_std,tcr_rmsd_min,tcr_rmsd_max,"
               "tcr_max_jump,processing_time_s,error_message\n")

        # 写入数据
        for r in results:
            f.write(f"{r['task_name']},{r['status']},{r['tcr_rmsd_mean']:.6f},"
                   f"{r['tcr_rmsd_std']:.6f},{r['tcr_rmsd_min']:.6f},{r['tcr_rmsd_max']:.6f},"
                   f"{r['tcr_max_jump']:.6f},{r['processing_time_s']:.2f},"
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
