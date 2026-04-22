#!/usr/bin/env python3
"""
Immunex RMSF Quick Start Example

演示如何从Immunex包中使用RMSF分析功能。
"""

import sys
from pathlib import Path

def example_basic_rmsf():
    """示例1: 基础RMSF分析"""
    print("="*60)
    print("示例1: 基础RMSF分析（使用RMSFAnalyzer）")
    print("="*60)

    code = '''
import MDAnalysis as mda
from immunex.analysis import RMSFAnalyzer

# 加载轨迹
u = mda.Universe("topology.pdb", "trajectory.xtc")

# 创建RMSF分析器
analyzer = RMSFAnalyzer(u, output_dir="./rmsf_results")

# 计算全局RMSF
global_stats = analyzer.calculate_global_rmsf()
print(f"平均RMSF: {global_stats['mean']:.2f} Å")
print(f"最大RMSF: {global_stats['max']:.2f} Å")

# 计算按链RMSF
chain_selections = {
    'HLA_alpha': 'chainID A and protein',
    'TCR_alpha': 'chainID D and protein'
}
chain_stats = analyzer.calculate_chain_rmsf(chain_selections)
'''

    print(code)
    print()

def example_phla_tcr_analyzer():
    """示例2: 使用pHLATCRAnalyzer"""
    print("="*60)
    print("示例2: 完整pHLA-TCR分析（含RMSF）")
    print("="*60)

    code = '''
from immunex.analysis import pHLATCRAnalyzer
from immunex.utils import pHLATCRVisualizer

# 创建分析器
analyzer = pHLATCRAnalyzer(
    task_name="1ao7_run1",
    trajectory="path/to/1ao7_processed.xtc",
    topology="path/to/1ao7_standardized.pdb",
    output_dir="./results"
)

# 运行完整分析（包括RMSF）
analyzer.analyze_global_rmsf()
analyzer.analyze_chain_rmsf()

# 如果有CDR3序列
cdr3_sequences = {
    'cdr3_alpha': 'AVTTDSWGKLQ',
    'cdr3_beta': 'ASRPGLAGGRPEQY'
}
analyzer.analyze_cdr3_rmsf(cdr3_sequences)

# 生成RMSF可视化
viz = pHLATCRVisualizer(analyzer.output_dir)
viz.plot_rmsf_comparison()
print("RMSF图表已生成!")
'''

    print(code)
    print()

def example_cdr3_only():
    """示例3: 仅CDR3 RMSF分析"""
    print("="*60)
    print("示例3: 专注CDR3 RMSF分析")
    print("="*60)

    code = '''
import MDAnalysis as mda
from immunex.analysis import calculate_cdr3_rmsf

# 加载轨迹
u = mda.Universe("topology.pdb", "trajectory.xtc")

# 计算CDR3α RMSF
alpha_stats = calculate_cdr3_rmsf(
    universe=u,
    chain_selection="chainID D and protein",
    cdr3_sequence="AVTTDSWGKLQ",
    chain_name="CDR3_alpha",
    output_dir="./cdr3_results"
)

if alpha_stats:
    print(f"CDR3α:")
    print(f"  平均RMSF: {alpha_stats['mean']:.2f} Å")
    print(f"  标准差: {alpha_stats['std']:.2f} Å")
    print(f"  残基范围: {alpha_stats['residue_range']}")
    print(f"  序列: {alpha_stats['sequence']}")
'''

    print(code)
    print()

def example_batch_analysis():
    """示例4: 批量分析"""
    print("="*60)
    print("示例4: 批量RMSF分析（命令行）")
    print("="*60)

    code = '''
# 激活环境
source /public/home/xmy/app/miniconda/bin/activate immunex

# 进入Immunex目录
cd /public/home/xmy/tools/Immunex

# 运行批量分析（Level 1+2）
python scripts/rmsf_analysis_pipeline.py \\
    /path/to/trajectories \\
    /path/to/topologies \\
    ./rmsf_results \\
    --workers 8

# 运行批量分析（Level 1+2+3，含CDR3）
python scripts/rmsf_analysis_pipeline.py \\
    /path/to/trajectories \\
    /path/to/topologies \\
    ./rmsf_results \\
    --cdr3-csv cdr3_sequences.csv \\
    --workers 8

# 查看结果
ls -lh rmsf_results/
head rmsf_results/rmsf_analysis_summary_*.csv
'''

    print(code)
    print()

def example_advanced():
    """示例5: 高级用法"""
    print("="*60)
    print("示例5: 高级RMSF分析")
    print("="*60)

    code = '''
import MDAnalysis as mda
from immunex.analysis import analyze_phla_tcr_rmsf
from pathlib import Path
import pandas as pd

# 加载轨迹
u = mda.Universe("topology.pdb", "trajectory.xtc")

# 运行完整pHLA-TCR RMSF分析
cdr3_sequences = {
    'cdr3_alpha': 'AVTTDSWGKLQ',
    'cdr3_beta': 'ASRPGLAGGRPEQY'
}

results = analyze_phla_tcr_rmsf(
    universe=u,
    output_dir=Path("./rmsf_results"),
    cdr3_sequences=cdr3_sequences
)

# 分析结果
print("\\n=== 全局RMSF ===")
print(f"平均: {results['global_rmsf']['mean']:.2f} Å")

print("\\n=== 按链RMSF ===")
for chain, stats in results['chain_rmsf'].items():
    print(f"{chain:12s}: {stats['mean']:.2f} ± {stats['std']:.2f} Å")

if 'cdr3_rmsf' in results:
    print("\\n=== CDR3 RMSF ===")
    for cdr, stats in results['cdr3_rmsf'].items():
        print(f"{cdr}:")
        print(f"  RMSF: {stats['mean']:.2f} ± {stats['std']:.2f} Å")
        print(f"  序列: {stats['sequence']}")
        print(f"  残基: {stats['residue_range']}")

# 导出数据进行后续分析
import json
with open('rmsf_summary.json', 'w') as f:
    json.dump(results, f, indent=2)
'''

    print(code)
    print()

def main():
    """运行所有示例"""
    print("\n")
    print("*"*60)
    print("Immunex RMSF功能快速入门")
    print("*"*60)
    print()

    print("以下是5个常见使用场景的代码示例：")
    print()

    example_basic_rmsf()
    example_phla_tcr_analyzer()
    example_cdr3_only()
    example_batch_analysis()
    example_advanced()

    print("="*60)
    print("更多详细信息，请查看:")
    print("  - RMSF_INTEGRATION.md (集成文档)")
    print("  - RMSF_PIPELINE_REPORT.md (完整报告)")
    print("="*60)
    print()

    print("测试导入:")
    try:
        from immunex.analysis import RMSFAnalyzer, pHLATCRAnalyzer
        from immunex.utils import pHLATCRVisualizer
        print("  ✓ 所有模块导入成功!")
    except ImportError as e:
        print(f"  ✗ 导入失败: {e}")
        return 1

    return 0

if __name__ == "__main__":
    sys.exit(main())
