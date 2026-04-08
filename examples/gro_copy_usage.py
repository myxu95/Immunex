#!/usr/bin/env python3
"""
示例：md.gro 文件复制功能演示

这个示例展示了 Immunex 如何自动复制 md.gro 文件到处理结果目录中，
便于后续轨迹可视化和分析。
"""

from immunex.preprocessing import PBCProcessor
from pathlib import Path

def demo_gro_copy():
    """演示 .gro 文件复制功能"""
    
    print("🔬 Immunex .gro 文件复制功能演示")
    print("=" * 50)
    
    # 示例文件路径 (请根据实际情况修改)
    example_structure = """
    典型的MD模拟目录结构:
    
    task_directory/
    ├── md.xtc          # 轨迹文件
    ├── md.tpr          # 拓扑文件
    ├── md.gro          # 结构文件 (会被复制)
    └── prod/           # 或者在prod子目录
        ├── md.xtc
        ├── md.tpr
        └── md.gro
    """
    print(example_structure)
    
    print("\n📋 .gro 文件搜索优先级:")
    print("1. md.gro (最高优先级)")
    print("2. prod.gro")
    print("3. production.gro") 
    print("4. {trajectory_name}.gro")
    
    print("\n📁 搜索位置:")
    print("- 轨迹文件所在目录")
    print("- 拓扑文件所在目录")
    
    print("\n✨ 处理后的输出目录将包含:")
    print("- {task_name}_processed.xtc  # 处理后的轨迹")
    print("- md.gro                     # 复制的结构文件")
    print("- reference_*.gro            # 其他找到的结构文件")
    
    print("\n🎯 使用方法:")
    code_example = '''
from immunex.preprocessing import PBCProcessor

# 初始化处理器
pbc_processor = PBCProcessor()

# 运行完整的PBC处理 (自动复制.gro文件)
results = pbc_processor.comprehensive_pbc_process(
    trajectory="path/to/md.xtc",
    topology="path/to/md.tpr", 
    output_dir="output_directory"
)

# 检查复制的结构文件
if results["reference_structures"]:
    print(f"复制了 {len(results['reference_structures'])} 个结构文件:")
    for gro_file in results["reference_structures"]:
        print(f"  - {gro_file}")
else:
    print("未找到 .gro 结构文件")
'''
    print(code_example)
    
    print("\n💡 优势:")
    print("- 🔍 自动搜索常见的结构文件名")
    print("- 📂 在多个位置查找文件")
    print("- 📋 记录复制的文件信息")
    print("- 🖼️ 便于后续轨迹可视化")
    print("- ⚡ 无需手动操作")

if __name__ == "__main__":
    demo_gro_copy()