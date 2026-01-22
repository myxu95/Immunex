#!/usr/bin/env python3
"""
Complete RMSF Analysis Pipeline for pHLA-TCR Trajectories

整合了三个层次的RMSF分析：
1. 全局RMSF (Global RMSF) - 整个蛋白质复合物
2. 按链RMSF (Per-Chain RMSF) - 5个链分别分析
3. CDR3 RMSF (CDR3-Specific RMSF) - TCR CDR3区域

Author: AfterMD Development Team
Date: 2025-12-09
"""

import pandas as pd
import numpy as np
from pathlib import Path
from multiprocessing import Pool
import logging
import sys
import argparse
from datetime import datetime
from typing import List, Dict, Optional

from phla_tcr_analyzer import pHLATCRAnalyzer
from visualization import pHLATCRVisualizer

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class RMSFAnalysisPipeline:
    """
    完整的RMSF分析Pipeline

    功能层次：
    Level 1: 全局RMSF - 整个复合物的柔性概览
    Level 2: 按链RMSF - 各链的柔性特征
    Level 3: CDR3 RMSF - 关键识别区域的柔性细节
    """

    def __init__(self,
                 trajectory_dir: str,
                 topology_dir: str,
                 output_dir: str,
                 cdr3_csv: Optional[str] = None,
                 max_workers: int = 8):
        """
        初始化RMSF分析pipeline

        Args:
            trajectory_dir: 轨迹目录
            topology_dir: 拓扑文件目录
            output_dir: 输出目录
            cdr3_csv: CDR3序列CSV文件路径（可选）
            max_workers: 并行worker数量
        """
        self.trajectory_dir = Path(trajectory_dir)
        self.topology_dir = Path(topology_dir)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.max_workers = max_workers

        # 加载CDR3序列（如果提供）
        self.cdr3_df = None
        if cdr3_csv:
            self.cdr3_csv = Path(cdr3_csv)
            if self.cdr3_csv.exists():
                self.cdr3_df = pd.read_csv(self.cdr3_csv)
                logger.info(f"加载CDR3序列: {len(self.cdr3_df)} 个任务")
            else:
                logger.warning(f"CDR3序列文件不存在: {cdr3_csv}")

        # 准备任务列表
        self.tasks = self._prepare_tasks()

        logger.info(f"初始化RMSF分析Pipeline: {len(self.tasks)} 个任务")
        logger.info(f"并行worker数: {self.max_workers}")
        logger.info(f"CDR3分析: {'启用' if self.cdr3_df is not None else '禁用'}")

    def _prepare_tasks(self) -> List[Dict]:
        """准备任务列表（智能匹配任务名）"""
        tasks = []

        # 扫描轨迹目录，建立PDB ID到任务名的映射
        pdb_to_tasks = {}

        for traj_dir in self.trajectory_dir.iterdir():
            if traj_dir.is_dir() and traj_dir.name != 'logs':
                task_name = traj_dir.name
                pdb_id = task_name.split('_')[0].lower()

                if pdb_id not in pdb_to_tasks:
                    pdb_to_tasks[pdb_id] = []
                pdb_to_tasks[pdb_id].append(task_name)

        logger.info(f"轨迹目录中发现 {sum(len(v) for v in pdb_to_tasks.values())} 个任务")

        # 如果提供了CDR3 CSV，使用它作为任务列表
        if self.cdr3_df is not None:
            for idx, row in self.cdr3_df.iterrows():
                cdr3_task_name = row['task_name']
                pdb_id = cdr3_task_name.split('_')[0].lower()

                # 智能匹配任务名
                matched_task_name = self._match_task_name(cdr3_task_name, pdb_id, pdb_to_tasks)

                if not matched_task_name:
                    continue

                # 构建任务信息
                task_info = self._build_task_info(matched_task_name, cdr3_task_name, row)

                if task_info:
                    tasks.append(task_info)
        else:
            # 没有CDR3 CSV，处理所有任务（仅Level 1和2）
            for pdb_id, task_names in pdb_to_tasks.items():
                for task_name in task_names:
                    task_info = self._build_task_info(task_name, None, None)

                    if task_info:
                        tasks.append(task_info)

        logger.info(f"准备了 {len(tasks)} 个任务")

        return tasks

    def _match_task_name(self, cdr3_task_name: str, pdb_id: str, pdb_to_tasks: Dict) -> Optional[str]:
        """智能匹配任务名"""
        if pdb_id not in pdb_to_tasks:
            return None

        candidates = pdb_to_tasks[pdb_id]

        # 策略1: 完全匹配
        if cdr3_task_name in candidates:
            return cdr3_task_name

        # 策略2: 只有一个候选
        if len(candidates) == 1:
            return candidates[0]

        # 策略3: 后缀匹配
        suffix = '_'.join(cdr3_task_name.split('_')[1:])
        for candidate in candidates:
            if candidate.endswith(suffix):
                return candidate

        # 策略4: 取第一个
        return candidates[0]

    def _build_task_info(self, task_name: str, cdr3_task_name: Optional[str], cdr3_row: Optional[pd.Series]) -> Optional[Dict]:
        """构建任务信息"""
        # 检查文件存在性
        trajectory_path = self.trajectory_dir / task_name / f"{task_name}_2step_processed.xtc"
        topology_path = self.topology_dir / task_name / f"{task_name}_standardized.pdb"

        if not trajectory_path.exists():
            logger.warning(f"轨迹文件不存在: {trajectory_path}")
            return None

        if not topology_path.exists():
            logger.warning(f"拓扑文件不存在: {topology_path}")
            return None

        task_info = {
            'task_name': task_name,
            'trajectory': str(trajectory_path),
            'topology': str(topology_path),
            'output_dir': str(self.output_dir / task_name),
            'cdr3_sequences': None
        }

        # 添加CDR3序列（如果有）
        if cdr3_row is not None:
            task_info['cdr3_sequences'] = {
                'cdr3_alpha': cdr3_row['tcra_cdr3'],
                'cdr3_beta': cdr3_row['tcrb_cdr3']
            }
            task_info['cdr3_task_name'] = cdr3_task_name

        return task_info

    def run_pipeline(self) -> pd.DataFrame:
        """运行完整的RMSF分析pipeline"""
        if not self.tasks:
            logger.error("没有任务可处理!")
            return pd.DataFrame()

        logger.info("="*60)
        logger.info("RMSF分析Pipeline")
        logger.info("="*60)
        logger.info(f"任务总数: {len(self.tasks)}")
        logger.info(f"分析层次:")
        logger.info(f"  Level 1: 全局RMSF - ✓")
        logger.info(f"  Level 2: 按链RMSF - ✓")
        logger.info(f"  Level 3: CDR3 RMSF - {'✓' if self.cdr3_df is not None else '✗'}")
        logger.info("="*60)

        start_time = datetime.now()

        # 并行处理
        with Pool(processes=self.max_workers) as pool:
            results = pool.map(_process_single_rmsf_task, self.tasks)

        # 汇总结果
        successful = [r for r in results if r is not None]
        failed = len(results) - len(successful)

        elapsed = (datetime.now() - start_time).total_seconds()

        logger.info("="*60)
        logger.info(f"Pipeline完成: {len(successful)} 成功, {failed} 失败")
        logger.info(f"总耗时: {elapsed:.1f} 秒 ({elapsed/60:.1f} 分钟)")
        logger.info("="*60)

        # 创建汇总DataFrame
        if successful:
            df_summary = pd.DataFrame(successful)

            # 保存汇总CSV
            timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
            summary_file = self.output_dir / f'rmsf_analysis_summary_{timestamp}.csv'
            df_summary.to_csv(summary_file, index=False)

            logger.info(f"汇总保存至: {summary_file}")

            # 打印统计信息
            self._print_statistics(df_summary)

            return df_summary
        else:
            logger.error("所有任务失败!")
            return pd.DataFrame()

    def _print_statistics(self, df: pd.DataFrame):
        """打印RMSF统计信息"""
        logger.info("\n" + "="*60)
        logger.info("RMSF分析统计")
        logger.info("="*60)

        # 基本统计
        logger.info(f"任务总数: {len(df)}")
        logger.info(f"平均处理时间: {df['processing_time_sec'].mean():.1f} 秒")

        # Level 1: 全局RMSF
        logger.info(f"\n【Level 1】全局RMSF统计:")
        if 'global_rmsf_mean' in df.columns:
            logger.info(f"  平均RMSF: {df['global_rmsf_mean'].mean():.2f} ± {df['global_rmsf_mean'].std():.2f} Å")
            logger.info(f"  范围: {df['global_rmsf_mean'].min():.2f} - {df['global_rmsf_mean'].max():.2f} Å")

        # Level 2: 按链RMSF
        logger.info(f"\n【Level 2】按链RMSF统计:")
        for chain in ['Peptide', 'HLA_alpha', 'HLA_beta', 'TCR_alpha', 'TCR_beta']:
            col = f'{chain}_rmsf_mean'
            if col in df.columns:
                logger.info(f"  {chain:12s}: {df[col].mean():.2f} ± {df[col].std():.2f} Å")

        # Level 3: CDR3 RMSF
        if 'CDR3_alpha_mean' in df.columns:
            logger.info(f"\n【Level 3】CDR3 RMSF统计:")

            # CDR3α
            valid_alpha = df[df['CDR3_alpha_mean'].notna()]
            if len(valid_alpha) > 0:
                logger.info(f"  CDR3α (n={len(valid_alpha)}):")
                logger.info(f"    平均: {valid_alpha['CDR3_alpha_mean'].mean():.2f} ± {valid_alpha['CDR3_alpha_mean'].std():.2f} Å")
                logger.info(f"    范围: {valid_alpha['CDR3_alpha_mean'].min():.2f} - {valid_alpha['CDR3_alpha_mean'].max():.2f} Å")

            # CDR3β
            valid_beta = df[df['CDR3_beta_mean'].notna()]
            if len(valid_beta) > 0:
                logger.info(f"  CDR3β (n={len(valid_beta)}):")
                logger.info(f"    平均: {valid_beta['CDR3_beta_mean'].mean():.2f} ± {valid_beta['CDR3_beta_mean'].std():.2f} Å")
                logger.info(f"    范围: {valid_beta['CDR3_beta_mean'].min():.2f} - {valid_beta['CDR3_beta_mean'].max():.2f} Å")

            # 序列比对成功率
            total = len(df)
            alpha_found = df['CDR3_alpha_mean'].notna().sum()
            beta_found = df['CDR3_beta_mean'].notna().sum()

            logger.info(f"\n  序列比对成功率:")
            logger.info(f"    CDR3α: {alpha_found}/{total} ({alpha_found/total*100:.1f}%)")
            logger.info(f"    CDR3β: {beta_found}/{total} ({beta_found/total*100:.1f}%)")

        logger.info("="*60 + "\n")


def _process_single_rmsf_task(task_info: Dict) -> Optional[Dict]:
    """
    处理单个RMSF分析任务

    Args:
        task_info: 任务信息字典

    Returns:
        RMSF分析结果字典
    """
    task_name = task_info['task_name']
    logger.info(f"处理任务: {task_name}")

    try:
        start_time = datetime.now()

        # 创建分析器
        analyzer = pHLATCRAnalyzer(
            task_name=task_name,
            trajectory=task_info['trajectory'],
            topology=task_info['topology'],
            output_dir=task_info['output_dir']
        )

        # Level 1 & 2: 全局和按链RMSF
        analyzer.analyze_global_rmsf()
        analyzer.analyze_chain_rmsf()

        # Level 3: CDR3 RMSF（如果提供了CDR3序列）
        if task_info['cdr3_sequences']:
            analyzer.analyze_cdr3_rmsf(task_info['cdr3_sequences'])

        # 生成可视化
        viz = pHLATCRVisualizer(task_info['output_dir'])
        viz.plot_rmsf_comparison()

        # 计算处理时间
        processing_time = (datetime.now() - start_time).total_seconds()

        # 提取结果
        summary = {
            'task_name': task_name,
            'processing_time_sec': processing_time,
        }

        # Level 1: 全局RMSF
        if 'global_rmsf' in analyzer.results:
            for key, val in analyzer.results['global_rmsf'].items():
                summary[f'global_rmsf_{key}'] = val

        # Level 2: 按链RMSF
        if 'chain_rmsf' in analyzer.results:
            for chain, stats in analyzer.results['chain_rmsf'].items():
                for key, val in stats.items():
                    summary[f'{chain}_rmsf_{key}'] = val

        # Level 3: CDR3 RMSF
        if 'cdr3_rmsf' in analyzer.results:
            for cdr, stats in analyzer.results['cdr3_rmsf'].items():
                for key, val in stats.items():
                    if key not in ['sequence', 'residue_range']:
                        summary[f'{cdr}_{key}'] = val
                    elif key == 'sequence':
                        summary[f'{cdr}_sequence'] = val

        logger.info(f"✓ 完成: {task_name} ({processing_time:.1f}秒)")

        return summary

    except Exception as e:
        logger.error(f"✗ 失败: {task_name} - {str(e)}")
        import traceback
        logger.debug(traceback.format_exc())
        return None


def main():
    parser = argparse.ArgumentParser(
        description='完整的RMSF分析Pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:

  # Level 1+2: 全局和按链RMSF（不包含CDR3）
  python rmsf_analysis_pipeline.py \\
      /public/home/xmy/dataset/workspace/pbc_1000frames_2step \\
      /public/home/xmy/dataset/workspace/1000frames_extract \\
      ./rmsf_results \\
      --workers 8

  # Level 1+2+3: 完整分析（包含CDR3）
  python rmsf_analysis_pipeline.py \\
      /public/home/xmy/dataset/workspace/pbc_1000frames_2step \\
      /public/home/xmy/dataset/workspace/1000frames_extract \\
      ./rmsf_results \\
      --cdr3-csv cdr3_sequences.csv \\
      --workers 8
        """
    )

    parser.add_argument('trajectory_dir', type=str,
                       help='轨迹目录 (pbc_1000frames_2step)')
    parser.add_argument('topology_dir', type=str,
                       help='拓扑文件目录 (1000frames_extract)')
    parser.add_argument('output_dir', type=str,
                       help='输出目录')
    parser.add_argument('--cdr3-csv', type=str, default=None,
                       help='CDR3序列CSV文件路径（可选，启用Level 3分析）')
    parser.add_argument('--workers', type=int, default=8,
                       help='并行worker数量 (默认: 8)')

    args = parser.parse_args()

    # 创建pipeline
    pipeline = RMSFAnalysisPipeline(
        trajectory_dir=args.trajectory_dir,
        topology_dir=args.topology_dir,
        output_dir=args.output_dir,
        cdr3_csv=args.cdr3_csv,
        max_workers=args.workers
    )

    # 运行分析
    df_summary = pipeline.run_pipeline()

    if not df_summary.empty:
        logger.info("RMSF分析Pipeline成功完成!")
        sys.exit(0)
    else:
        logger.error("RMSF分析Pipeline失败!")
        sys.exit(1)


if __name__ == '__main__':
    main()
