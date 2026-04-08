#!/usr/bin/env python3
"""
MD Product对比分析模块

收集所有MD Product的分析结果，进行统计对比和可视化。
包括RMSD对比、回转半径对比、质量统计等多维度比较。
"""

import json
import logging
from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from collections import defaultdict

logger = logging.getLogger(__name__)


class ComparativeAnalyzer:
    """
    MD Product对比分析器
    
    收集和比较多个MD Product的分析结果，生成对比图表和统计报告。
    """
    
    def __init__(self, output_dir: str):
        """
        初始化对比分析器
        
        Args:
            output_dir: 对比分析结果输出目录
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # 创建子目录
        self.plots_dir = self.output_dir / "plots"
        self.data_dir = self.output_dir / "data"
        self.reports_dir = self.output_dir / "reports"
        
        for dir_path in [self.plots_dir, self.data_dir, self.reports_dir]:
            dir_path.mkdir(parents=True, exist_ok=True)
        
        # 存储收集的数据
        self.md_products = {}
        self.comparative_results = {}
        
        # 配置matplotlib中文字体 (如果需要)
        plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'Liberation Sans']
        plt.rcParams['axes.unicode_minus'] = False
        
        logger.info(f"对比分析器初始化完成，输出目录: {self.output_dir}")
    
    def collect_analysis_results(self, processed_root: str) -> Dict[str, Dict]:
        """
        收集所有MD Product的分析结果
        
        Args:
            processed_root: 包含已处理MD Product的根目录
            
        Returns:
            收集到的分析结果字典
        """
        logger.info(f"开始收集分析结果: {processed_root}")
        
        processed_root_path = Path(processed_root)
        collected_data = {}
        
        # 查找所有包含分析结果的目录
        for analysis_dir in processed_root_path.rglob("analysis"):
            if not analysis_dir.is_dir():
                continue
            
            # 获取MD Product名称
            md_product_name = analysis_dir.parent.name
            
            # 查找分析结果JSON文件
            results_file = analysis_dir / "reports" / "analysis_results.json"
            if not results_file.exists():
                logger.warning(f"未找到分析结果文件: {results_file}")
                continue
            
            try:
                with open(results_file, 'r', encoding='utf-8') as f:
                    analysis_data = json.load(f)
                
                collected_data[md_product_name] = {
                    "analysis_data": analysis_data,
                    "analysis_dir": str(analysis_dir),
                    "input_files": analysis_data.get("input_files", {}),
                    "analysis_summary": analysis_data.get("analysis_summary", {}),
                    "analysis_results": analysis_data.get("analysis_results", {})
                }
                
                logger.debug(f"收集到 {md_product_name} 的分析结果")
                
            except Exception as e:
                logger.error(f"读取分析结果失败 {results_file}: {e}")
                continue
        
        self.md_products = collected_data
        logger.info(f"成功收集 {len(collected_data)} 个MD Product的分析结果")
        
        return collected_data
    
    def extract_rmsd_data(self) -> pd.DataFrame:
        """
        提取所有MD Product的RMSD数据
        
        Returns:
            包含RMSD统计数据的DataFrame
        """
        rmsd_data = []
        
        for md_name, md_data in self.md_products.items():
            rmsd_results = md_data["analysis_results"].get("rmsd", {})
            
            for selection_name, rmsd_info in rmsd_results.items():
                if "error" in rmsd_info:
                    continue
                
                rmsd_data.append({
                    "MD_Product": md_name,
                    "Selection": selection_name,
                    "Mean_RMSD": rmsd_info.get("mean_rmsd", 0),
                    "Std_RMSD": rmsd_info.get("std_rmsd", 0),
                    "Max_RMSD": rmsd_info.get("max_rmsd", 0),
                    "Min_RMSD": rmsd_info.get("min_rmsd", 0),
                    "Selection_Description": rmsd_info.get("selection", "")
                })
        
        df = pd.DataFrame(rmsd_data)
        logger.info(f"提取RMSD数据: {len(df)} 条记录")
        return df
    
    def extract_radius_gyration_data(self) -> pd.DataFrame:
        """
        提取所有MD Product的回转半径数据
        
        Returns:
            包含回转半径统计数据的DataFrame
        """
        rg_data = []
        
        for md_name, md_data in self.md_products.items():
            rg_results = md_data["analysis_results"].get("radius_gyration", {})
            
            for selection_name, rg_info in rg_results.items():
                if "error" in rg_info:
                    continue
                
                rg_data.append({
                    "MD_Product": md_name,
                    "Selection": selection_name,
                    "Mean_Rg": rg_info.get("mean_rg", 0),
                    "Std_Rg": rg_info.get("std_rg", 0),
                    "Max_Rg": rg_info.get("max_rg", 0),
                    "Min_Rg": rg_info.get("min_rg", 0),
                    "Selection_Description": rg_info.get("selection", "")
                })
        
        df = pd.DataFrame(rg_data)
        logger.info(f"提取回转半径数据: {len(df)} 条记录")
        return df
    
    def extract_distance_data(self) -> pd.DataFrame:
        """
        提取所有MD Product的距离数据
        
        Returns:
            包含距离统计数据的DataFrame
        """
        distance_data = []
        
        for md_name, md_data in self.md_products.items():
            distance_results = md_data["analysis_results"].get("distances", {})
            
            for distance_name, distance_info in distance_results.items():
                if "error" in distance_info:
                    continue
                
                distance_data.append({
                    "MD_Product": md_name,
                    "Distance_Type": distance_name,
                    "Mean_Distance": distance_info.get("mean_distance", 0),
                    "Std_Distance": distance_info.get("std_distance", 0),
                    "Max_Distance": distance_info.get("max_distance", 0),
                    "Min_Distance": distance_info.get("min_distance", 0),
                    "Selection1": distance_info.get("selection1", ""),
                    "Selection2": distance_info.get("selection2", "")
                })
        
        df = pd.DataFrame(distance_data)
        logger.info(f"提取距离数据: {len(df)} 条记录")
        return df
    
    def plot_rmsd_comparison(self, rmsd_df: pd.DataFrame) -> List[str]:
        """
        生成RMSD对比图表
        
        Args:
            rmsd_df: RMSD数据DataFrame
            
        Returns:
            生成的图表文件路径列表
        """
        plot_files = []
        
        # 1. 按选择类型分组的RMSD对比
        for selection in rmsd_df['Selection'].unique():
            selection_data = rmsd_df[rmsd_df['Selection'] == selection]
            
            if len(selection_data) < 2:
                continue
            
            plt.figure(figsize=(12, 8))
            
            # 创建条形图
            x_pos = np.arange(len(selection_data))
            bars = plt.bar(x_pos, selection_data['Mean_RMSD'], 
                          yerr=selection_data['Std_RMSD'],
                          capsize=5, alpha=0.7)
            
            # 设置颜色
            colors = plt.cm.Set3(np.linspace(0, 1, len(bars)))
            for bar, color in zip(bars, colors):
                bar.set_color(color)
            
            plt.xlabel('MD Products')
            plt.ylabel('RMSD (Å)')
            plt.title(f'RMSD Comparison - {selection}')
            plt.xticks(x_pos, selection_data['MD_Product'], rotation=45, ha='right')
            plt.grid(axis='y', alpha=0.3)
            plt.tight_layout()
            
            # 保存图表
            plot_file = self.plots_dir / f"rmsd_comparison_{selection}.png"
            plt.savefig(plot_file, dpi=300, bbox_inches='tight')
            plt.close()
            
            plot_files.append(str(plot_file))
            logger.info(f"生成RMSD对比图: {plot_file.name}")
        
        # 2. 所有选择类型的热力图
        if len(rmsd_df) > 0:
            pivot_data = rmsd_df.pivot(index='MD_Product', columns='Selection', values='Mean_RMSD')
            
            plt.figure(figsize=(14, 8))
            sns.heatmap(pivot_data, annot=True, fmt='.2f', cmap='YlOrRd', 
                       cbar_kws={'label': 'RMSD (Å)'})
            plt.title('RMSD Heatmap - All MD Products and Selections')
            plt.xlabel('Atom Selections')
            plt.ylabel('MD Products')
            plt.xticks(rotation=45, ha='right')
            plt.yticks(rotation=0)
            plt.tight_layout()
            
            heatmap_file = self.plots_dir / "rmsd_heatmap_all.png"
            plt.savefig(heatmap_file, dpi=300, bbox_inches='tight')
            plt.close()
            
            plot_files.append(str(heatmap_file))
            logger.info(f"生成RMSD热力图: {heatmap_file.name}")
        
        return plot_files
    
    def plot_radius_gyration_comparison(self, rg_df: pd.DataFrame) -> List[str]:
        """
        生成回转半径对比图表
        
        Args:
            rg_df: 回转半径数据DataFrame
            
        Returns:
            生成的图表文件路径列表
        """
        plot_files = []
        
        # 按选择类型分组的回转半径对比
        for selection in rg_df['Selection'].unique():
            selection_data = rg_df[rg_df['Selection'] == selection]
            
            if len(selection_data) < 2:
                continue
            
            plt.figure(figsize=(12, 8))
            
            # 创建条形图
            x_pos = np.arange(len(selection_data))
            bars = plt.bar(x_pos, selection_data['Mean_Rg'], 
                          yerr=selection_data['Std_Rg'],
                          capsize=5, alpha=0.7)
            
            # 设置颜色
            colors = plt.cm.Set2(np.linspace(0, 1, len(bars)))
            for bar, color in zip(bars, colors):
                bar.set_color(color)
            
            plt.xlabel('MD Products')
            plt.ylabel('Radius of Gyration (Å)')
            plt.title(f'Radius of Gyration Comparison - {selection}')
            plt.xticks(x_pos, selection_data['MD_Product'], rotation=45, ha='right')
            plt.grid(axis='y', alpha=0.3)
            plt.tight_layout()
            
            # 保存图表
            plot_file = self.plots_dir / f"rg_comparison_{selection}.png"
            plt.savefig(plot_file, dpi=300, bbox_inches='tight')
            plt.close()
            
            plot_files.append(str(plot_file))
            logger.info(f"生成回转半径对比图: {plot_file.name}")
        
        return plot_files
    
    def plot_distance_comparison(self, distance_df: pd.DataFrame) -> List[str]:
        """
        生成距离对比图表
        
        Args:
            distance_df: 距离数据DataFrame
            
        Returns:
            生成的图表文件路径列表
        """
        plot_files = []
        
        # 按距离类型分组的对比
        for distance_type in distance_df['Distance_Type'].unique():
            type_data = distance_df[distance_df['Distance_Type'] == distance_type]
            
            if len(type_data) < 2:
                continue
            
            plt.figure(figsize=(12, 8))
            
            # 创建条形图
            x_pos = np.arange(len(type_data))
            bars = plt.bar(x_pos, type_data['Mean_Distance'], 
                          yerr=type_data['Std_Distance'],
                          capsize=5, alpha=0.7)
            
            # 设置颜色
            colors = plt.cm.Pastel1(np.linspace(0, 1, len(bars)))
            for bar, color in zip(bars, colors):
                bar.set_color(color)
            
            plt.xlabel('MD Products')
            plt.ylabel('Distance (Å)')
            plt.title(f'Distance Comparison - {distance_type}')
            plt.xticks(x_pos, type_data['MD_Product'], rotation=45, ha='right')
            plt.grid(axis='y', alpha=0.3)
            plt.tight_layout()
            
            # 保存图表
            plot_file = self.plots_dir / f"distance_comparison_{distance_type}.png"
            plt.savefig(plot_file, dpi=300, bbox_inches='tight')
            plt.close()
            
            plot_files.append(str(plot_file))
            logger.info(f"生成距离对比图: {plot_file.name}")
        
        return plot_files
    
    def plot_quality_statistics(self) -> List[str]:
        """
        生成质量统计对比图表
        
        Returns:
            生成的图表文件路径列表
        """
        plot_files = []
        
        # 提取质量数据
        quality_data = []
        for md_name, md_data in self.md_products.items():
            input_files = md_data.get("input_files", {})
            analysis_summary = md_data.get("analysis_summary", {})
            
            quality_data.append({
                "MD_Product": md_name,
                "Trajectory_Size_MB": input_files.get("trajectory_size_mb", 0),
                "N_Frames": input_files.get("n_frames", 0),
                "Analysis_Time_Seconds": analysis_summary.get("total_time_seconds", 0),
                "Completed_Analyses": len(analysis_summary.get("completed_analyses", [])),
                "Failed_Analyses": len(analysis_summary.get("failed_analyses", []))
            })
        
        if not quality_data:
            return plot_files
        
        quality_df = pd.DataFrame(quality_data)
        
        # 1. 轨迹文件大小对比
        plt.figure(figsize=(12, 6))
        bars = plt.bar(quality_df['MD_Product'], quality_df['Trajectory_Size_MB'], alpha=0.7)
        plt.xlabel('MD Products')
        plt.ylabel('Trajectory Size (MB)')
        plt.title('Trajectory File Size Comparison')
        plt.xticks(rotation=45, ha='right')
        plt.grid(axis='y', alpha=0.3)
        plt.tight_layout()
        
        size_plot = self.plots_dir / "trajectory_size_comparison.png"
        plt.savefig(size_plot, dpi=300, bbox_inches='tight')
        plt.close()
        plot_files.append(str(size_plot))
        
        # 2. 分析时间对比
        plt.figure(figsize=(12, 6))
        bars = plt.bar(quality_df['MD_Product'], quality_df['Analysis_Time_Seconds'], alpha=0.7, color='orange')
        plt.xlabel('MD Products')
        plt.ylabel('Analysis Time (seconds)')
        plt.title('Analysis Time Comparison')
        plt.xticks(rotation=45, ha='right')
        plt.grid(axis='y', alpha=0.3)
        plt.tight_layout()
        
        time_plot = self.plots_dir / "analysis_time_comparison.png"
        plt.savefig(time_plot, dpi=300, bbox_inches='tight')
        plt.close()
        plot_files.append(str(time_plot))
        
        # 3. 分析成功率
        plt.figure(figsize=(12, 6))
        success_rate = quality_df['Completed_Analyses'] / (quality_df['Completed_Analyses'] + quality_df['Failed_Analyses'] + 0.001) * 100
        bars = plt.bar(quality_df['MD_Product'], success_rate, alpha=0.7, color='green')
        plt.xlabel('MD Products')
        plt.ylabel('Analysis Success Rate (%)')
        plt.title('Analysis Success Rate Comparison')
        plt.xticks(rotation=45, ha='right')
        plt.ylim(0, 105)
        plt.grid(axis='y', alpha=0.3)
        plt.tight_layout()
        
        success_plot = self.plots_dir / "analysis_success_rate.png"
        plt.savefig(success_plot, dpi=300, bbox_inches='tight')
        plt.close()
        plot_files.append(str(success_plot))
        
        logger.info(f"生成质量统计图表: {len(plot_files)} 个")
        return plot_files
    
    def plot_summary_dashboard(self) -> str:
        """
        生成综合仪表板图表
        
        Returns:
            仪表板图表文件路径
        """
        # 提取关键指标
        rmsd_df = self.extract_rmsd_data()
        rg_df = self.extract_radius_gyration_data()
        
        if rmsd_df.empty and rg_df.empty:
            logger.warning("无足够数据生成仪表板")
            return ""
        
        # 创建2x2子图
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
        
        # 1. RMSD汇总 (α碳原子)
        if not rmsd_df.empty:
            ca_rmsd = rmsd_df[rmsd_df['Selection'].str.contains('CA', na=False)]
            if not ca_rmsd.empty:
                ax1.bar(ca_rmsd['MD_Product'], ca_rmsd['Mean_RMSD'], alpha=0.7)
                ax1.set_title('RMSD Comparison (CA atoms)')
                ax1.set_ylabel('RMSD (Å)')
                ax1.tick_params(axis='x', rotation=45)
                ax1.grid(axis='y', alpha=0.3)
        
        # 2. 回转半径汇总
        if not rg_df.empty:
            protein_rg = rg_df[rg_df['Selection'].str.contains('protein', na=False)]
            if not protein_rg.empty:
                ax2.bar(protein_rg['MD_Product'], protein_rg['Mean_Rg'], alpha=0.7, color='orange')
                ax2.set_title('Radius of Gyration Comparison')
                ax2.set_ylabel('Rg (Å)')
                ax2.tick_params(axis='x', rotation=45)
                ax2.grid(axis='y', alpha=0.3)
        
        # 3. 质量指标散点图
        quality_data = []
        for md_name, md_data in self.md_products.items():
            input_files = md_data.get("input_files", {})
            analysis_summary = md_data.get("analysis_summary", {})
            
            quality_data.append({
                "MD_Product": md_name,
                "Size_MB": input_files.get("trajectory_size_mb", 0),
                "Time_s": analysis_summary.get("total_time_seconds", 0)
            })
        
        if quality_data:
            quality_df = pd.DataFrame(quality_data)
            scatter = ax3.scatter(quality_df['Size_MB'], quality_df['Time_s'], alpha=0.7, s=100)
            ax3.set_xlabel('Trajectory Size (MB)')
            ax3.set_ylabel('Analysis Time (s)')
            ax3.set_title('Size vs Analysis Time')
            ax3.grid(alpha=0.3)
            
            # 添加标签
            for i, row in quality_df.iterrows():
                ax3.annotate(row['MD_Product'], (row['Size_MB'], row['Time_s']), 
                           xytext=(5, 5), textcoords='offset points', fontsize=8)
        
        # 4. 分析类型统计
        analysis_counts = defaultdict(int)
        for md_data in self.md_products.values():
            completed = md_data.get("analysis_summary", {}).get("completed_analyses", [])
            for analysis_type in completed:
                analysis_counts[analysis_type] += 1
        
        if analysis_counts:
            labels = list(analysis_counts.keys())
            counts = list(analysis_counts.values())
            ax4.pie(counts, labels=labels, autopct='%1.1f%%', startangle=90)
            ax4.set_title('Analysis Type Distribution')
        
        plt.tight_layout()
        
        # 保存仪表板
        dashboard_file = self.plots_dir / "comparative_dashboard.png"
        plt.savefig(dashboard_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"生成综合仪表板: {dashboard_file}")
        return str(dashboard_file)
    
    def generate_comparative_report(self) -> str:
        """
        生成对比分析报告
        
        Returns:
            报告文件路径
        """
        report_file = self.reports_dir / "comparative_analysis_report.md"
        
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write("# MD Products Comparative Analysis Report\n\n")
            
            # 基本统计
            f.write("## Basic Statistics\n\n")
            f.write(f"- **Total MD Products Analyzed**: {len(self.md_products)}\n")
            
            # 质量统计
            total_size = sum(md_data.get("input_files", {}).get("trajectory_size_mb", 0) 
                           for md_data in self.md_products.values())
            total_time = sum(md_data.get("analysis_summary", {}).get("total_time_seconds", 0) 
                           for md_data in self.md_products.values())
            
            f.write(f"- **Total Trajectory Size**: {total_size:.1f} MB\n")
            f.write(f"- **Total Analysis Time**: {total_time:.1f} seconds\n")
            f.write(f"- **Average Analysis Time**: {total_time/len(self.md_products):.1f} seconds/MD\n\n")
            
            # RMSD统计
            rmsd_df = self.extract_rmsd_data()
            if not rmsd_df.empty:
                f.write("## RMSD Analysis Summary\n\n")
                
                # 按选择类型统计
                for selection in rmsd_df['Selection'].unique():
                    selection_data = rmsd_df[rmsd_df['Selection'] == selection]
                    f.write(f"### {selection}\n")
                    f.write(f"- **Number of MD Products**: {len(selection_data)}\n")
                    f.write(f"- **Mean RMSD Range**: {selection_data['Mean_RMSD'].min():.2f} - {selection_data['Mean_RMSD'].max():.2f} Å\n")
                    f.write(f"- **Average RMSD**: {selection_data['Mean_RMSD'].mean():.2f} ± {selection_data['Mean_RMSD'].std():.2f} Å\n\n")
            
            # 回转半径统计
            rg_df = self.extract_radius_gyration_data()
            if not rg_df.empty:
                f.write("## Radius of Gyration Analysis Summary\n\n")
                
                for selection in rg_df['Selection'].unique():
                    selection_data = rg_df[rg_df['Selection'] == selection]
                    f.write(f"### {selection}\n")
                    f.write(f"- **Number of MD Products**: {len(selection_data)}\n")
                    f.write(f"- **Mean Rg Range**: {selection_data['Mean_Rg'].min():.2f} - {selection_data['Mean_Rg'].max():.2f} Å\n")
                    f.write(f"- **Average Rg**: {selection_data['Mean_Rg'].mean():.2f} ± {selection_data['Mean_Rg'].std():.2f} Å\n\n")
            
            # 距离分析统计
            distance_df = self.extract_distance_data()
            if not distance_df.empty:
                f.write("## Distance Analysis Summary\n\n")
                
                for distance_type in distance_df['Distance_Type'].unique():
                    type_data = distance_df[distance_df['Distance_Type'] == distance_type]
                    f.write(f"### {distance_type}\n")
                    f.write(f"- **Number of MD Products**: {len(type_data)}\n")
                    f.write(f"- **Mean Distance Range**: {type_data['Mean_Distance'].min():.2f} - {type_data['Mean_Distance'].max():.2f} Å\n")
                    f.write(f"- **Average Distance**: {type_data['Mean_Distance'].mean():.2f} ± {type_data['Mean_Distance'].std():.2f} Å\n\n")
            
            # 生成的图表列表
            f.write("## Generated Visualizations\n\n")
            plot_files = list(self.plots_dir.glob("*.png"))
            for plot_file in sorted(plot_files):
                f.write(f"- {plot_file.name}\n")
            
            f.write(f"\n## Data Files\n\n")
            data_files = list(self.data_dir.glob("*.csv"))
            for data_file in sorted(data_files):
                f.write(f"- {data_file.name}\n")
        
        logger.info(f"生成对比分析报告: {report_file}")
        return str(report_file)
    
    def save_comparative_data(self):
        """保存对比分析数据到CSV文件"""
        # 保存RMSD数据
        rmsd_df = self.extract_rmsd_data()
        if not rmsd_df.empty:
            rmsd_file = self.data_dir / "rmsd_comparison_data.csv"
            rmsd_df.to_csv(rmsd_file, index=False)
            logger.info(f"保存RMSD对比数据: {rmsd_file}")
        
        # 保存回转半径数据
        rg_df = self.extract_radius_gyration_data()
        if not rg_df.empty:
            rg_file = self.data_dir / "radius_gyration_comparison_data.csv"
            rg_df.to_csv(rg_file, index=False)
            logger.info(f"保存回转半径对比数据: {rg_file}")
        
        # 保存距离数据
        distance_df = self.extract_distance_data()
        if not distance_df.empty:
            distance_file = self.data_dir / "distance_comparison_data.csv"
            distance_df.to_csv(distance_file, index=False)
            logger.info(f"保存距离对比数据: {distance_file}")
    
    def run_comprehensive_comparison(self, processed_root: str) -> Dict[str, Any]:
        """
        运行完整的对比分析流程
        
        Args:
            processed_root: 包含已处理MD Product的根目录
            
        Returns:
            对比分析结果摘要
        """
        logger.info("开始运行完整对比分析...")
        
        # 1. 收集数据
        collected_data = self.collect_analysis_results(processed_root)
        
        if not collected_data:
            logger.error("未收集到任何分析结果，无法进行对比分析")
            return {"status": "failed", "reason": "no_data"}
        
        # 2. 生成对比图表
        plot_files = []
        
        # RMSD对比图
        rmsd_df = self.extract_rmsd_data()
        if not rmsd_df.empty:
            rmsd_plots = self.plot_rmsd_comparison(rmsd_df)
            plot_files.extend(rmsd_plots)
        
        # 回转半径对比图
        rg_df = self.extract_radius_gyration_data()
        if not rg_df.empty:
            rg_plots = self.plot_radius_gyration_comparison(rg_df)
            plot_files.extend(rg_plots)
        
        # 距离对比图
        distance_df = self.extract_distance_data()
        if not distance_df.empty:
            distance_plots = self.plot_distance_comparison(distance_df)
            plot_files.extend(distance_plots)
        
        # 质量统计图
        quality_plots = self.plot_quality_statistics()
        plot_files.extend(quality_plots)
        
        # 综合仪表板
        dashboard_file = self.plot_summary_dashboard()
        if dashboard_file:
            plot_files.append(dashboard_file)
        
        # 3. 保存数据
        self.save_comparative_data()
        
        # 4. 生成报告
        report_file = self.generate_comparative_report()
        
        # 5. 返回结果摘要
        result = {
            "status": "success",
            "total_md_products": len(collected_data),
            "plots_generated": len(plot_files),
            "data_files_saved": len(list(self.data_dir.glob("*.csv"))),
            "report_file": report_file,
            "output_directory": str(self.output_dir),
            "plot_files": plot_files
        }
        
        logger.info(f"对比分析完成: {result['total_md_products']} MD Products, {result['plots_generated']} 图表")
        
        return result


def create_comparative_analyzer(output_dir: str) -> ComparativeAnalyzer:
    """
    工厂函数：创建对比分析器实例
    
    Args:
        output_dir: 输出目录
        
    Returns:
        ComparativeAnalyzer实例
    """
    return ComparativeAnalyzer(output_dir=output_dir)