#!/usr/bin/env python3
"""
Immunex Analysis Pipeline

综合分析流水线，基于质量检测和PBC预处理后的规范文件进行后续分析。
包括RMSD、RDF、回转半径、距离分析、氢键分析等。
"""

import os
import logging
from pathlib import Path
from typing import Dict, List, Optional, Any, Union
import pandas as pd
import numpy as np

from ..analysis.trajectory.rmsd import RMSDCalculator
from ..analysis.trajectory.rdf import RDFCalculator
from ..analysis.trajectory.radius_gyration import RadiusGyrationCalculator
from ..analysis.trajectory.distance import DistanceCalculator
from ..analysis.trajectory.hydrogen_bonds import HydrogenBondAnalyzer
from ..utils.plotting import PlotManager
from ..utils.path_manager import PathManager

logger = logging.getLogger(__name__)


class AnalysisPipeline:
    """
    MD轨迹综合分析流水线
    
    基于标准化的输入文件格式进行各种轨迹分析:
    - 输入: PBC处理后的 .xtc 文件 + .gro/.tpr 参考结构
    - 输出: 分析结果 + 可视化图表 + 汇总报告
    """
    
    def __init__(self, 
                 processed_trajectory: str,
                 reference_structure: str,
                 output_dir: str,
                 gmx_executable: str = "gmx"):
        """
        初始化分析流水线
        
        Args:
            processed_trajectory: PBC处理后的轨迹文件路径 (.xtc)
            reference_structure: 参考结构文件路径 (.gro/.tpr/.pdb)
            output_dir: 分析结果输出目录
            gmx_executable: GROMACS可执行文件
        """
        self.processed_trajectory = Path(processed_trajectory)
        self.reference_structure = Path(reference_structure)
        self.output_dir = Path(output_dir)
        self.gmx = gmx_executable
        
        # 验证输入文件
        self._validate_inputs()
        
        # 创建输出目录结构
        self.path_manager = PathManager(base_dir=str(self.output_dir))
        self._setup_output_directories()
        
        # 初始化分析器
        self.analyzers = {}
        self._initialize_analyzers()
        
        # 初始化绘图管理器
        self.plot_manager = PlotManager(output_dir=str(self.output_dir / "plots"))
        
        # 存储分析结果
        self.results = {
            "input_files": {
                "trajectory": str(self.processed_trajectory),
                "reference": str(self.reference_structure),
                "trajectory_size_mb": self.processed_trajectory.stat().st_size / (1024*1024),
                "n_frames": None  # 将在首次分析时确定
            },
            "analysis_results": {},
            "plots_generated": {},
            "summary_statistics": {}
        }
    
    def _validate_inputs(self):
        """验证输入文件的有效性"""
        if not self.processed_trajectory.exists():
            raise FileNotFoundError(f"轨迹文件不存在: {self.processed_trajectory}")
        
        if not self.reference_structure.exists():
            raise FileNotFoundError(f"参考结构文件不存在: {self.reference_structure}")
        
        # 检查文件格式
        if self.processed_trajectory.suffix not in ['.xtc', '.trr']:
            raise ValueError(f"不支持的轨迹文件格式: {self.processed_trajectory.suffix}")
        
        if self.reference_structure.suffix not in ['.gro', '.tpr', '.pdb']:
            raise ValueError(f"不支持的参考结构格式: {self.reference_structure.suffix}")
        
        logger.info(f"输入文件验证通过:")
        logger.info(f"  轨迹: {self.processed_trajectory} ({self.processed_trajectory.stat().st_size/(1024*1024):.1f} MB)")
        logger.info(f"  结构: {self.reference_structure}")
    
    def _setup_output_directories(self):
        """创建分析输出目录结构"""
        directories = [
            "rmsd",
            "rdf", 
            "radius_gyration",
            "distances",
            "hydrogen_bonds",
            "plots",
            "data",
            "reports"
        ]
        
        for dir_name in directories:
            (self.output_dir / dir_name).mkdir(parents=True, exist_ok=True)
        
        logger.info(f"输出目录结构创建完成: {self.output_dir}")
    
    def _initialize_analyzers(self):
        """初始化各种分析器"""
        try:
            # 基础参数
            traj_str = str(self.processed_trajectory)
            ref_str = str(self.reference_structure)
            
            # RMSD分析器
            self.analyzers['rmsd'] = RMSDCalculator(
                topology=ref_str,
                trajectory=traj_str,
                gmx_executable=self.gmx
            )
            
            # RDF分析器  
            self.analyzers['rdf'] = RDFCalculator(
                topology=ref_str,
                trajectory=traj_str,
                gmx_executable=self.gmx
            )
            
            # 回转半径分析器
            self.analyzers['radius_gyration'] = RadiusGyrationCalculator(
                topology=ref_str,
                trajectory=traj_str,
                gmx_executable=self.gmx
            )
            
            # 距离分析器
            self.analyzers['distance'] = DistanceCalculator(
                topology=ref_str,
                trajectory=traj_str,
                gmx_executable=self.gmx
            )
            
            # 氢键分析器
            self.analyzers['hydrogen_bonds'] = HydrogenBondAnalyzer(
                topology=ref_str,
                trajectory=traj_str,
                gmx_executable=self.gmx
            )
            
            logger.info(f"成功初始化 {len(self.analyzers)} 个分析器")
            
        except Exception as e:
            logger.error(f"分析器初始化失败: {e}")
            raise
    
    def run_rmsd_analysis(self, 
                         selections: Optional[List[str]] = None,
                         reference_frame: int = 0,
                         **kwargs) -> Dict[str, Any]:
        """
        运行RMSD分析
        
        Args:
            selections: 原子选择列表，默认为蛋白质骨架原子
            reference_frame: 参考帧索引
            **kwargs: 其他RMSD分析参数
            
        Returns:
            RMSD分析结果字典
        """
        logger.info("开始RMSD分析...")
        
        if selections is None:
            selections = [
                "protein and name CA",        # α碳原子
                "protein and name CA and resname ALA VAL LEU ILE PHE TRP MET",  # 疏水残基
                "protein and name CA and resname SER THR ASN GLN",              # 极性残基
                "protein and backbone"        # 蛋白质骨架
            ]
        
        rmsd_results = {}
        
        for i, selection in enumerate(selections):
            try:
                selection_name = f"selection_{i+1}_{selection.replace(' ', '_').replace('and', '').strip('_')}"
                
                logger.info(f"计算RMSD: {selection}")
                
                # 使用MDAnalysis计算RMSD
                times, rmsd_values = self.analyzers['rmsd'].calculate_mdanalysis(
                    selection=selection,
                    reference_frame=reference_frame,
                    output_file=str(self.output_dir / "data" / f"rmsd_{selection_name}.xvg")
                )
                
                rmsd_results[selection_name] = {
                    "selection": selection,
                    "times": times,
                    "rmsd_values": rmsd_values,
                    "mean_rmsd": np.mean(rmsd_values),
                    "std_rmsd": np.std(rmsd_values),
                    "max_rmsd": np.max(rmsd_values),
                    "min_rmsd": np.min(rmsd_values)
                }
                
                # 生成RMSD图
                plot_file = self.output_dir / "plots" / f"rmsd_{selection_name}.png"
                self.plot_manager.plot_rmsd(
                    times=times,
                    rmsd_values=rmsd_values,
                    output_file=str(plot_file),
                    title=f"RMSD Analysis - {selection}",
                    xlabel="Time (ps)",
                    ylabel="RMSD (Å)"
                )
                
                logger.info(f"RMSD分析完成: {selection_name} (平均值: {np.mean(rmsd_values):.2f} Å)")
                
            except Exception as e:
                logger.error(f"RMSD分析失败 - {selection}: {e}")
                rmsd_results[selection_name] = {"error": str(e)}
        
        # 更新结果
        self.results["analysis_results"]["rmsd"] = rmsd_results
        if not self.results["input_files"]["n_frames"] and rmsd_results:
            # 从第一个成功的分析获取帧数
            for result in rmsd_results.values():
                if "times" in result:
                    self.results["input_files"]["n_frames"] = len(result["times"])
                    break
        
        logger.info(f"RMSD分析完成，共分析 {len([r for r in rmsd_results.values() if 'error' not in r])} 个选择")
        return rmsd_results
    
    def run_radius_gyration_analysis(self, 
                                   selections: Optional[List[str]] = None,
                                   **kwargs) -> Dict[str, Any]:
        """
        运行回转半径分析
        
        Args:
            selections: 原子选择列表
            **kwargs: 其他参数
            
        Returns:
            回转半径分析结果
        """
        logger.info("开始回转半径分析...")
        
        if selections is None:
            selections = [
                "protein",                    # 整个蛋白质
                "protein and name CA",        # α碳原子
            ]
        
        rg_results = {}
        
        for selection in selections:
            try:
                selection_name = selection.replace(' ', '_').replace('and', '')
                
                logger.info(f"计算回转半径: {selection}")
                
                times, rg_values = self.analyzers['radius_gyration'].calculate_mdanalysis(
                    selection=selection,
                    output_file=str(self.output_dir / "data" / f"radius_gyration_{selection_name}.xvg")
                )
                
                rg_results[selection_name] = {
                    "selection": selection,
                    "times": times,
                    "rg_values": rg_values,
                    "mean_rg": np.mean(rg_values),
                    "std_rg": np.std(rg_values),
                    "max_rg": np.max(rg_values),
                    "min_rg": np.min(rg_values)
                }
                
                # 生成回转半径图
                plot_file = self.output_dir / "plots" / f"radius_gyration_{selection_name}.png"
                self.plot_manager.plot_time_series(
                    times=times,
                    values=rg_values,
                    output_file=str(plot_file),
                    title=f"Radius of Gyration - {selection}",
                    xlabel="Time (ps)",
                    ylabel="Radius of Gyration (Å)"
                )
                
                logger.info(f"回转半径分析完成: {selection_name} (平均值: {np.mean(rg_values):.2f} Å)")
                
            except Exception as e:
                logger.error(f"回转半径分析失败 - {selection}: {e}")
                rg_results[selection_name] = {"error": str(e)}
        
        self.results["analysis_results"]["radius_gyration"] = rg_results
        return rg_results
    
    def run_distance_analysis(self,
                            distance_pairs: Optional[List[Dict]] = None,
                            **kwargs) -> Dict[str, Any]:
        """
        运行距离分析
        
        Args:
            distance_pairs: 距离对列表，每个元素为 {"name": "描述", "sel1": "选择1", "sel2": "选择2"}
            **kwargs: 其他参数
            
        Returns:
            距离分析结果
        """
        logger.info("开始距离分析...")
        
        if distance_pairs is None:
            # 默认分析一些常见的距离
            distance_pairs = [
                {
                    "name": "protein_center_of_mass", 
                    "sel1": "protein", 
                    "sel2": "protein",
                    "type": "center_of_mass"
                },
                {
                    "name": "protein_ends",
                    "sel1": "protein and resid 1", 
                    "sel2": "protein and resid 999",  # 自动调整为最后一个残基
                    "type": "minimum_distance"
                }
            ]
        
        distance_results = {}
        
        for pair in distance_pairs:
            try:
                pair_name = pair["name"]
                
                logger.info(f"计算距离: {pair_name}")
                
                times, distances = self.analyzers['distance'].calculate_distance(
                    selection1=pair["sel1"],
                    selection2=pair["sel2"],
                    distance_type=pair.get("type", "center_of_mass"),
                    output_file=str(self.output_dir / "data" / f"distance_{pair_name}.xvg")
                )
                
                distance_results[pair_name] = {
                    "selection1": pair["sel1"],
                    "selection2": pair["sel2"],
                    "distance_type": pair.get("type", "center_of_mass"),
                    "times": times,
                    "distances": distances,
                    "mean_distance": np.mean(distances),
                    "std_distance": np.std(distances),
                    "max_distance": np.max(distances),
                    "min_distance": np.min(distances)
                }
                
                # 生成距离图
                plot_file = self.output_dir / "plots" / f"distance_{pair_name}.png"
                self.plot_manager.plot_time_series(
                    times=times,
                    values=distances,
                    output_file=str(plot_file),
                    title=f"Distance Analysis - {pair_name}",
                    xlabel="Time (ps)",
                    ylabel="Distance (Å)"
                )
                
                logger.info(f"距离分析完成: {pair_name} (平均值: {np.mean(distances):.2f} Å)")
                
            except Exception as e:
                logger.error(f"距离分析失败 - {pair.get('name', 'unknown')}: {e}")
                distance_results[pair.get('name', 'unknown')] = {"error": str(e)}
        
        self.results["analysis_results"]["distances"] = distance_results
        return distance_results
    
    def run_comprehensive_analysis(self, 
                                 enable_rmsd: bool = True,
                                 enable_rg: bool = True,
                                 enable_distances: bool = True,
                                 enable_rdf: bool = False,  # RDF计算较慢，默认关闭
                                 enable_hbonds: bool = False,  # 氢键分析较慢，默认关闭
                                 **kwargs) -> Dict[str, Any]:
        """
        运行综合分析流水线
        
        Args:
            enable_rmsd: 是否启用RMSD分析
            enable_rg: 是否启用回转半径分析
            enable_distances: 是否启用距离分析
            enable_rdf: 是否启用RDF分析
            enable_hbonds: 是否启用氢键分析
            **kwargs: 传递给各分析模块的参数
            
        Returns:
            综合分析结果
        """
        logger.info("开始综合分析流水线...")
        
        # 记录开始时间
        import time
        start_time = time.time()
        
        analysis_summary = {
            "started_at": time.strftime("%Y-%m-%d %H:%M:%S"),
            "completed_analyses": [],
            "failed_analyses": [],
            "total_plots_generated": 0
        }
        
        # RMSD分析
        if enable_rmsd:
            try:
                self.run_rmsd_analysis(**kwargs.get('rmsd_params', {}))
                analysis_summary["completed_analyses"].append("RMSD")
            except Exception as e:
                logger.error(f"RMSD分析失败: {e}")
                analysis_summary["failed_analyses"].append(f"RMSD: {e}")
        
        # 回转半径分析
        if enable_rg:
            try:
                self.run_radius_gyration_analysis(**kwargs.get('rg_params', {}))
                analysis_summary["completed_analyses"].append("Radius of Gyration")
            except Exception as e:
                logger.error(f"回转半径分析失败: {e}")
                analysis_summary["failed_analyses"].append(f"Radius of Gyration: {e}")
        
        # 距离分析
        if enable_distances:
            try:
                self.run_distance_analysis(**kwargs.get('distance_params', {}))
                analysis_summary["completed_analyses"].append("Distance Analysis")
            except Exception as e:
                logger.error(f"距离分析失败: {e}")
                analysis_summary["failed_analyses"].append(f"Distance Analysis: {e}")
        
        # RDF分析 (可选)
        if enable_rdf:
            try:
                # 这里需要实现RDF分析调用
                logger.info("RDF分析已启用但尚未实现")
                analysis_summary["completed_analyses"].append("RDF (placeholder)")
            except Exception as e:
                logger.error(f"RDF分析失败: {e}")
                analysis_summary["failed_analyses"].append(f"RDF: {e}")
        
        # 氢键分析 (可选)
        if enable_hbonds:
            try:
                # 这里需要实现氢键分析调用
                logger.info("氢键分析已启用但尚未实现")
                analysis_summary["completed_analyses"].append("Hydrogen Bonds (placeholder)")
            except Exception as e:
                logger.error(f"氢键分析失败: {e}")
                analysis_summary["failed_analyses"].append(f"Hydrogen Bonds: {e}")
        
        # 计算总耗时
        end_time = time.time()
        analysis_summary["total_time_seconds"] = end_time - start_time
        analysis_summary["completed_at"] = time.strftime("%Y-%m-%d %H:%M:%S")
        
        # 统计生成的图片数量
        plots_dir = self.output_dir / "plots"
        if plots_dir.exists():
            analysis_summary["total_plots_generated"] = len(list(plots_dir.glob("*.png")))
        
        # 保存综合分析结果
        self.results["analysis_summary"] = analysis_summary
        
        # 生成分析报告
        self.generate_analysis_report()
        
        logger.info(f"综合分析完成! 耗时: {analysis_summary['total_time_seconds']:.1f}秒")
        logger.info(f"完成的分析: {', '.join(analysis_summary['completed_analyses'])}")
        if analysis_summary["failed_analyses"]:
            logger.warning(f"失败的分析: {', '.join(analysis_summary['failed_analyses'])}")
        
        return self.results
    
    def generate_analysis_report(self):
        """生成分析报告"""
        report_file = self.output_dir / "reports" / "analysis_report.md"
        
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write("# MD Trajectory Analysis Report\n\n")
            
            # 基本信息
            f.write("## Input Information\n")
            f.write(f"- **Trajectory**: {self.results['input_files']['trajectory']}\n")
            f.write(f"- **Reference Structure**: {self.results['input_files']['reference']}\n")
            f.write(f"- **Trajectory Size**: {self.results['input_files']['trajectory_size_mb']:.1f} MB\n")
            if self.results['input_files']['n_frames']:
                f.write(f"- **Number of Frames**: {self.results['input_files']['n_frames']}\n")
            f.write("\n")
            
            # 分析摘要
            summary = self.results.get("analysis_summary", {})
            f.write("## Analysis Summary\n")
            f.write(f"- **Started**: {summary.get('started_at', 'N/A')}\n")
            f.write(f"- **Completed**: {summary.get('completed_at', 'N/A')}\n")
            f.write(f"- **Total Time**: {summary.get('total_time_seconds', 0):.1f} seconds\n")
            f.write(f"- **Completed Analyses**: {', '.join(summary.get('completed_analyses', []))}\n")
            f.write(f"- **Generated Plots**: {summary.get('total_plots_generated', 0)}\n")
            f.write("\n")
            
            # 各分析结果摘要
            for analysis_type, results in self.results.get("analysis_results", {}).items():
                f.write(f"## {analysis_type.upper()} Results\n")
                
                for selection_name, data in results.items():
                    if "error" in data:
                        f.write(f"- **{selection_name}**: Failed - {data['error']}\n")
                    else:
                        if analysis_type == "rmsd":
                            f.write(f"- **{selection_name}**: Mean RMSD = {data.get('mean_rmsd', 0):.2f} ± {data.get('std_rmsd', 0):.2f} Å\n")
                        elif analysis_type == "radius_gyration":
                            f.write(f"- **{selection_name}**: Mean Rg = {data.get('mean_rg', 0):.2f} ± {data.get('std_rg', 0):.2f} Å\n")
                        elif analysis_type == "distances":
                            f.write(f"- **{selection_name}**: Mean Distance = {data.get('mean_distance', 0):.2f} ± {data.get('std_distance', 0):.2f} Å\n")
                
                f.write("\n")
            
            f.write("## Generated Files\n")
            f.write("- Data files: `data/`\n")
            f.write("- Plots: `plots/`\n")
            f.write("- This report: `reports/analysis_report.md`\n")
        
        logger.info(f"分析报告已生成: {report_file}")
    
    def save_results(self, output_file: Optional[str] = None):
        """保存分析结果到JSON文件"""
        import json
        
        if output_file is None:
            output_file = self.output_dir / "reports" / "analysis_results.json"
        
        # 转换numpy数组为列表以便JSON序列化
        def convert_numpy(obj):
            if isinstance(obj, np.ndarray):
                return obj.tolist()
            elif isinstance(obj, np.floating):
                return float(obj)
            elif isinstance(obj, np.integer):
                return int(obj)
            return obj
        
        def deep_convert(obj):
            if isinstance(obj, dict):
                return {k: deep_convert(v) for k, v in obj.items()}
            elif isinstance(obj, list):
                return [deep_convert(v) for v in obj]
            else:
                return convert_numpy(obj)
        
        results_serializable = deep_convert(self.results)
        
        with open(output_file, 'w', encoding='utf-8') as f:
            json.dump(results_serializable, f, indent=2, ensure_ascii=False)
        
        logger.info(f"分析结果已保存: {output_file}")


def create_analysis_pipeline(processed_trajectory: str,
                           reference_structure: str, 
                           output_dir: str,
                           **kwargs) -> AnalysisPipeline:
    """
    工厂函数：创建分析流水线实例
    
    Args:
        processed_trajectory: PBC处理后的轨迹文件
        reference_structure: 参考结构文件
        output_dir: 输出目录
        **kwargs: 其他参数
        
    Returns:
        AnalysisPipeline实例
    """
    return AnalysisPipeline(
        processed_trajectory=processed_trajectory,
        reference_structure=reference_structure,
        output_dir=output_dir,
        **kwargs
    )
