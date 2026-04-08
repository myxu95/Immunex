#!/usr/bin/env python3
"""
前台MD质量检查脚本

在生成SLURM脚本之前，快速检查所有MD产品的质量，
只为通过质量检查的MD生成后续的处理脚本。

使用方法:
    python md_quality_check.py /path/to/md/simulations
    python md_quality_check.py /path/to/md/simulations --expected-chains 3
"""

import sys
import argparse
import json
import time
import numpy as np
from pathlib import Path
from datetime import datetime
from typing import List, Dict, Tuple

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent))

from immunex.analysis.quality import (
    MDCompletenessChecker,
    StructureValidator,
    BatchTracker,
    QualityReporter
)
from quality_analysis import find_pdb_files


class NumpyEncoder(json.JSONEncoder):
    """JSON encoder for numpy types."""
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        return super(NumpyEncoder, self).default(obj)


class ForegroundQualityChecker:
    """前台质量检查器"""
    
    def __init__(self,
                 expected_chains: int = 5,
                 min_traj_size_mb: float = 1.0,
                 min_sim_time_ps: float = 5000.0):
        """
        初始化质量检查器
        
        Args:
            expected_chains: 期望的蛋白质链数
            min_traj_size_mb: 最小轨迹文件大小(MB)
            min_sim_time_ps: 最小模拟时间(ps)
        """
        self.expected_chains = expected_chains
        self.min_traj_size_mb = min_traj_size_mb
        self.min_sim_time_ps = min_sim_time_ps
        
        # 初始化检查器
        self.md_checker = MDCompletenessChecker(
            min_trajectory_size_mb=min_traj_size_mb,
            min_simulation_time_ps=min_sim_time_ps
        )
        self.structure_validator = StructureValidator(
            expected_chain_count=expected_chains
        )
        self.batch_tracker = BatchTracker()
        self.quality_reporter = QualityReporter()
        
        # 结果存储
        self.results = {
            'start_time': datetime.now().isoformat(),
            'parameters': {
                'expected_chains': expected_chains,
                'min_traj_size_mb': min_traj_size_mb,
                'min_sim_time_ps': min_sim_time_ps
            },
            'total_found': 0,
            'qualified_count': 0,
            'failed_count': 0,
            'qualified_mds': [],
            'failed_mds': [],
            'quality_summary': {},
            'detailed_results': {}
        }
    
    def find_md_directories(self, root_path: Path) -> List[Path]:
        """
        查找所有MD Product目录
        
        MD Product结构应该是：
        MD_Product_Dir/
        ├── {pdb_id}.pdb          # PDB文件在根目录
        ├── prod/                 # 可选：MD文件在prod子目录
        │   ├── md.gro
        │   ├── md.xtc
        │   └── md.log
        └── md.gro                # 或者：MD文件直接在根目录
            md.xtc
            md.log
        """
        print("🔍 搜索MD Product目录...")
        
        md_dirs = []
        
        # 查找MD Product目录 (包含PDB文件的目录)
        for potential_dir in root_path.rglob("*"):
            if not potential_dir.is_dir():
                continue
                
            # 检查是否包含符合PDB ID命名规范的PDB文件
            pdb_files = self._find_pdb_files_in_product(potential_dir)
            if not pdb_files:
                continue
            
            # 检查MD文件位置
            md_files_location = self._find_md_files_location(potential_dir)
            
            if md_files_location:
                md_dirs.append(potential_dir)
        
        print(f"📁 发现 {len(md_dirs)} 个MD Product目录")
        return sorted(md_dirs)
    
    def _find_md_files_location(self, md_product_dir: Path) -> str:
        """
        在MD Product目录中查找MD文件的位置
        
        Returns:
            'root': MD文件在根目录
            'prod': MD文件在prod子目录  
            None: 未找到MD文件
        """
        # 1. 检查根目录是否有MD文件
        root_has_gro = any(md_product_dir.glob("md.gro")) or any(md_product_dir.glob("*.gro"))
        root_has_xtc = any(md_product_dir.glob("md.xtc")) or any(md_product_dir.glob("*.xtc")) or any(md_product_dir.glob("*.trr"))
        
        if root_has_gro and root_has_xtc:
            return 'root'
        
        # 2. 检查prod子目录是否有MD文件
        prod_dir = md_product_dir / 'prod'
        if prod_dir.exists() and prod_dir.is_dir():
            prod_has_gro = any(prod_dir.glob("md.gro")) or any(prod_dir.glob("*.gro"))
            prod_has_xtc = any(prod_dir.glob("md.xtc")) or any(prod_dir.glob("*.xtc")) or any(prod_dir.glob("*.trr"))
            
            if prod_has_gro and prod_has_xtc:
                return 'prod'
        
        return None
    
    def _find_pdb_files_in_product(self, md_product_dir: Path) -> List[Path]:
        """
        在MD Product根目录中查找PDB文件
        
        Args:
            md_product_dir: MD Product根目录
            
        Returns:
            PDB文件列表 (只包含符合PDB ID命名规范的文件)
        """
        from quality_analysis import is_pdb_id_filename
        
        pdb_files = []
        
        # 只在根目录查找PDB文件，不递归搜索
        for pdb_file in md_product_dir.glob("*.pdb"):
            if is_pdb_id_filename(pdb_file.name):
                pdb_files.append(pdb_file)
        
        return sorted(pdb_files)
    
    def check_single_md(self, md_dir: Path) -> Dict:
        """检查单个MD的质量"""
        result = {
            'md_path': str(md_dir),
            'md_name': md_dir.name,
            'is_qualified': False,
            'issues': [],
            'md_completeness': {},
            'structure_validation': {},
            'qualification_reasons': []
        }
        
        try:
            # 确定MD文件的实际位置
            md_files_location = self._find_md_files_location(md_dir)
            
            if not md_files_location:
                result['issues'].append('md_files_not_found')
                result['qualification_reasons'].append("MD文件缺失: 未找到md.gro和md.xtc文件")
                return result
            
            # 根据文件位置确定检查路径
            if md_files_location == 'prod':
                md_check_path = md_dir / 'prod'
            else:
                md_check_path = md_dir
            
            # 1. MD完整性检查
            print(f"   ├─ MD完整性检查... ({md_files_location})")
            md_result = self.md_checker.check_single_md(str(md_check_path))
            result['md_completeness'] = md_result
            result['md_completeness']['files_location'] = md_files_location
            
            # 检查MD是否完整
            if md_result.get('status') != 'complete':
                result['issues'].append('md_incomplete')
                result['qualification_reasons'].append(f"MD未完成: {md_result.get('message', 'Unknown')}")
                return result
            
            # 2. PDB结构验证 (PDB文件在MD Product根目录)
            print(f"   ├─ PDB结构验证...")
            pdb_files = self._find_pdb_files_in_product(md_dir)
            
            if pdb_files:
                # 只验证第一个PDB文件作为代表
                structure_result = self.structure_validator.validate_structure(pdb_files[0])
                result['structure_validation'] = structure_result
                
                # 检查链数
                chain_analysis = structure_result.get('chain_analysis', {})
                protein_chains = chain_analysis.get('protein_chains', 0)
                
                if protein_chains != self.expected_chains:
                    result['issues'].append('unexpected_chain_count')
                    result['qualification_reasons'].append(
                        f"链数异常: 发现{protein_chains}条链，期望{self.expected_chains}条链"
                    )
                    return result
            else:
                result['qualification_reasons'].append("未找到PDB文件进行结构验证")
            
            # 3. 如果所有检查都通过
            result['is_qualified'] = True
            result['qualification_reasons'].append("通过所有质量检查")
            
        except Exception as e:
            result['issues'].append('check_error')
            result['qualification_reasons'].append(f"检查过程出错: {str(e)}")
        
        return result
    
    def run_quality_check(self, md_root: Path) -> Dict:
        """运行完整的质量检查流程"""
        print("🚀 开始MD质量检查")
        print("=" * 60)
        print(f"📂 检查目录: {md_root}")
        print(f"🔧 期望链数: {self.expected_chains}")
        print(f"📏 最小轨迹: {self.min_traj_size_mb} MB")
        print(f"⏱️  最小时间: {self.min_sim_time_ps} ps")
        print("=" * 60)
        
        # 1. 查找MD目录
        md_dirs = self.find_md_directories(md_root)
        self.results['total_found'] = len(md_dirs)
        
        if not md_dirs:
            print("❌ 未找到任何MD目录！")
            return self.results
        
        # 2. 逐个检查MD质量
        print(f"\n📋 开始检查 {len(md_dirs)} 个MD目录...")
        
        for i, md_dir in enumerate(md_dirs, 1):
            print(f"\n[{i:3d}/{len(md_dirs)}] {md_dir.name}")
            
            # 检查单个MD
            check_result = self.check_single_md(md_dir)
            
            # 存储详细结果
            self.results['detailed_results'][str(md_dir)] = check_result
            
            # 分类结果
            if check_result['is_qualified']:
                self.results['qualified_mds'].append(str(md_dir))
                self.results['qualified_count'] += 1
                print(f"   ✅ 合格")
            else:
                self.results['failed_mds'].append({
                    'path': str(md_dir),
                    'reasons': check_result['qualification_reasons']
                })
                self.results['failed_count'] += 1
                print(f"   ❌ 不合格: {'; '.join(check_result['qualification_reasons'])}")
        
        # 3. 生成统计摘要
        self.results['end_time'] = datetime.now().isoformat()
        self.results['quality_summary'] = self._generate_summary()
        
        return self.results
    
    def _generate_summary(self) -> Dict:
        """生成质量检查摘要"""
        total = self.results['total_found']
        qualified = self.results['qualified_count']
        failed = self.results['failed_count']
        
        summary = {
            'total_checked': total,
            'qualified': qualified,
            'failed': failed,
            'qualification_rate': (qualified / total * 100) if total > 0 else 0,
            'main_failure_reasons': {}
        }
        
        # 统计失败原因
        reason_counts = {}
        for failed_md in self.results['failed_mds']:
            for reason in failed_md['reasons']:
                reason_counts[reason] = reason_counts.get(reason, 0) + 1
        
        summary['main_failure_reasons'] = dict(sorted(reason_counts.items(), 
                                                     key=lambda x: x[1], reverse=True))
        
        return summary
    
    def save_results(self, output_dir: Path):
        """保存检查结果"""
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # 1. 保存详细的JSON报告
        detailed_report = output_dir / "quality_check_report.json"
        with open(detailed_report, 'w', encoding='utf-8') as f:
            # Use custom encoder for numpy types
            json.dump(self.results, f, indent=2, ensure_ascii=False, cls=NumpyEncoder)
        
        # 2. 保存合格MD列表
        qualified_list = output_dir / "qualified_mds.txt"
        with open(qualified_list, 'w') as f:
            for md_path in self.results['qualified_mds']:
                f.write(f"{md_path}\n")
        
        # 3. 保存失败MD列表
        failed_list = output_dir / "failed_mds.txt"
        with open(failed_list, 'w', encoding='utf-8') as f:
            for failed_md in self.results['failed_mds']:
                f.write(f"{failed_md['path']}\t{'; '.join(failed_md['reasons'])}\n")
        
        # 4. 保存摘要报告
        summary_report = output_dir / "quality_summary.txt"
        with open(summary_report, 'w', encoding='utf-8') as f:
            self._write_summary_report(f)
        
        return {
            'detailed_report': detailed_report,
            'qualified_list': qualified_list,
            'failed_list': failed_list,
            'summary_report': summary_report
        }
    
    def _write_summary_report(self, file):
        """写入摘要报告"""
        summary = self.results['quality_summary']
        
        file.write("MD PRODUCTION QUALITY REPORT\n")
        file.write("=" * 80 + "\n")
        file.write(f"检查时间: {self.results['start_time']}\n")
        file.write(f"检查目录: {self.results.get('check_directory', 'N/A')}\n")
        file.write(f"期望链数: {self.results['parameters']['expected_chains']} | ")
        file.write(f"最小轨迹: {self.results['parameters']['min_traj_size_mb']}MB | ")
        file.write(f"最小时间: {self.results['parameters']['min_sim_time_ps']}ps\n")
        file.write("=" * 80 + "\n\n")
        
        # 1. 详细表格
        file.write("1. MD PRODUCTS DETAILED STATUS\n")
        file.write("-" * 80 + "\n")
        
        # 表头
        header = f"{'MD Product':<18} {'Chains':<8} {'Progress':<12} {'Traj(MB)':<10} {'Time(ps)':<10} {'Status':<12} {'Issues':<20}\n"
        file.write(header)
        file.write("-" * 80 + "\n")
        
        # 表格内容
        for md_path, details in self.results['detailed_results'].items():
            md_name = Path(md_path).name[:17]  # 截断过长名称
            
            # 获取链数
            chain_analysis = details.get('structure_validation', {}).get('chain_analysis', {})
            chains = chain_analysis.get('protein_chains', 'N/A')
            
            # 获取MD进度
            md_completeness = details.get('md_completeness', {})
            if md_completeness.get('status') == 'complete':
                progress = "Complete"
            elif md_completeness.get('status') == 'incomplete':
                progress = "Incomplete"
            else:
                progress = "Unknown"
            
            # 获取轨迹大小
            traj_size = md_completeness.get('trajectory_size_mb', 'N/A')
            if isinstance(traj_size, (int, float)):
                traj_size = f"{traj_size:.1f}"
            
            # 获取模拟时间
            sim_time = md_completeness.get('simulation_time_ps', 'N/A')
            if isinstance(sim_time, (int, float)):
                sim_time = f"{sim_time:.0f}"
            
            # 状态
            status = "PASS" if details['is_qualified'] else "FAIL"
            
            # 问题描述
            issues = "; ".join(details.get('qualification_reasons', ['']))[:19]
            
            line = f"{md_name:<18} {chains:<8} {progress:<12} {traj_size:<10} {sim_time:<10} {status:<12} {issues:<20}\n"
            file.write(line)
        
        file.write("-" * 80 + "\n\n")
        
        # 2. 统计摘要
        file.write("2. STATISTICS SUMMARY\n")
        file.write("-" * 40 + "\n")
        file.write(f"{'Total Checked:':<20} {summary['total_checked']:>8} MD Products\n")
        file.write(f"{'Qualified:':<20} {summary['qualified']:>8} ({summary['qualification_rate']:.1f}%)\n")
        file.write(f"{'Failed:':<20} {summary['failed']:>8} ({100-summary['qualification_rate']:.1f}%)\n\n")
        
        # 3. 异常统计
        if summary['main_failure_reasons']:
            file.write("3. FAILURE ANALYSIS\n")
            file.write("-" * 40 + "\n")
            for reason, count in summary['main_failure_reasons'].items():
                percentage = (count / summary['total_checked'] * 100)
                file.write(f"{'• ' + reason:<30} {count:>3} ({percentage:.1f}%)\n")
        
        file.write(f"\n4. RECOMMENDATIONS\n")
        file.write("-" * 40 + "\n")
        if summary['qualification_rate'] >= 80:
            file.write("✅ 大部分MD质量良好，可以继续批量处理\n")
        elif summary['qualification_rate'] >= 50:
            file.write("⚠️  部分MD存在质量问题，建议检查失败原因\n")
        else:
            file.write("❌ 大量MD存在质量问题，建议检查模拟参数和流程\n")
    
    def print_results(self):
        """打印检查结果"""
        summary = self.results['quality_summary']
        
        print(f"\n{'='*80}")
        print("🎯 MD PRODUCTION QUALITY REPORT")
        print(f"{'='*80}")
        
        # 简化的表格预览 (只显示前10个)
        print("TOP MD PRODUCTS STATUS:")
        print("-" * 80)
        header = f"{'MD Product':<18} {'Chains':<8} {'Progress':<12} {'Status':<10} {'Issues':<25}"
        print(header)
        print("-" * 80)
        
        # 显示前10个MD的状态
        count = 0
        for md_path, details in self.results['detailed_results'].items():
            if count >= 10:
                break
            
            md_name = Path(md_path).name[:17]
            
            # 获取链数
            chain_analysis = details.get('structure_validation', {}).get('chain_analysis', {})
            chains = chain_analysis.get('protein_chains', 'N/A')
            
            # 获取MD进度
            md_completeness = details.get('md_completeness', {})
            if md_completeness.get('status') == 'complete':
                progress = "Complete"
            elif md_completeness.get('status') == 'incomplete':
                progress = "Incomplete"
            else:
                progress = "Unknown"
            
            # 状态
            status = "PASS" if details['is_qualified'] else "FAIL"
            
            # 问题描述
            issues = "; ".join(details.get('qualification_reasons', ['']))[:24]
            
            print(f"{md_name:<18} {chains:<8} {progress:<12} {status:<10} {issues:<25}")
            count += 1
        
        if len(self.results['detailed_results']) > 10:
            remaining = len(self.results['detailed_results']) - 10
            print(f"... 以及其他 {remaining} 个MD (查看完整报告)")
        
        print("-" * 80)
        
        # 统计摘要
        print(f"\n📊 SUMMARY STATISTICS:")
        print(f"   Total Checked: {summary['total_checked']:>3} MD Products")
        print(f"   Qualified:     {summary['qualified']:>3} ({summary['qualification_rate']:.1f}%)")
        print(f"   Failed:        {summary['failed']:>3} ({100-summary['qualification_rate']:.1f}%)")
        
        # 失败原因统计
        if summary['main_failure_reasons']:
            print(f"\n🔍 FAILURE BREAKDOWN:")
            for reason, count in list(summary['main_failure_reasons'].items())[:3]:
                percentage = (count / summary['total_checked'] * 100)
                print(f"   • {reason}: {count} ({percentage:.1f}%)")
        
        print(f"\n💡 NEXT STEPS:")
        if summary['qualified'] > 0:
            print(f"   ✅ {summary['qualified']} 个合格MD可以进行批量处理")
            print(f"   📋 查看完整报告: quality_summary.txt")
            print(f"   🚀 生成SLURM脚本处理合格MD")
        else:
            print(f"   ⚠️  没有合格的MD，请检查质量问题")
            print(f"   📋 查看详细失败原因: quality_summary.txt")


def main():
    """主函数"""
    parser = argparse.ArgumentParser(
        description="MD质量前台检查工具",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例:
  基本检查:
    python md_quality_check.py /path/to/md/simulations
    
  自定义参数:
    python md_quality_check.py /path/to/md/simulations \\
        --expected-chains 3 \\
        --min-traj-size 5.0 \\
        --min-sim-time 10000
    
  指定输出目录:
    python md_quality_check.py /path/to/md/simulations \\
        --output ./quality_results
        """
    )
    
    # 必需参数
    parser.add_argument(
        "md_root",
        help="MD模拟根目录路径"
    )
    
    # 质量检查参数
    parser.add_argument(
        "--expected-chains",
        type=int,
        default=5,
        help="期望的蛋白质链数 (默认: 5)"
    )
    
    parser.add_argument(
        "--min-traj-size",
        type=float,
        default=1.0,
        help="最小轨迹文件大小(MB) (默认: 1.0)"
    )
    
    parser.add_argument(
        "--min-sim-time",
        type=float,
        default=5000.0,
        help="最小模拟时间(ps) (默认: 5000.0)"
    )
    
    # 输出参数
    parser.add_argument(
        "--output", "-o",
        type=Path,
        default="./quality_check_results",
        help="质量检查结果输出目录 (默认: ./quality_check_results)"
    )
    
    # 其他选项
    parser.add_argument(
        "--quiet", "-q",
        action="store_true",
        help="静默模式，减少输出信息"
    )
    
    args = parser.parse_args()
    
    # 验证输入目录
    md_root = Path(args.md_root)
    if not md_root.exists():
        print(f"❌ 错误: 目录不存在 - {md_root}")
        return 1
    
    if not md_root.is_dir():
        print(f"❌ 错误: 路径不是目录 - {md_root}")
        return 1
    
    try:
        # 创建质量检查器
        checker = ForegroundQualityChecker(
            expected_chains=args.expected_chains,
            min_traj_size_mb=args.min_traj_size,
            min_sim_time_ps=args.min_sim_time
        )
        
        # 存储检查目录信息
        checker.results['check_directory'] = str(md_root)
        
        # 运行质量检查
        start_time = time.time()
        results = checker.run_quality_check(md_root)
        end_time = time.time()
        
        # 保存结果
        saved_files = checker.save_results(args.output)
        
        # 显示结果
        if not args.quiet:
            checker.print_results()
        
        # 显示文件位置
        print(f"\n📁 结果已保存到: {args.output}")
        print(f"   • 详细报告: {saved_files['detailed_report'].name}")
        print(f"   • 合格列表: {saved_files['qualified_list'].name}")
        print(f"   • 失败列表: {saved_files['failed_list'].name}")
        print(f"   • 摘要报告: {saved_files['summary_report'].name}")
        
        print(f"\n⏱️  检查耗时: {end_time - start_time:.1f} 秒")
        
        # 返回状态码
        if results['qualified_count'] > 0:
            print(f"\n🎉 发现 {results['qualified_count']} 个合格MD，可以继续生成SLURM脚本！")
            return 0
        else:
            print(f"\n⚠️  没有发现合格的MD，请检查质量问题")
            return 1
            
    except Exception as e:
        print(f"\n❌ 检查过程出错: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    exit(main())