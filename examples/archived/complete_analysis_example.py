#!/usr/bin/env python3
"""
Immunex完整分析流程示例

演示从质量检测到轨迹分析的完整工作流程。
"""

import sys
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from immunex.pipeline import create_analysis_pipeline


def example_single_trajectory_analysis():
    """
    示例：分析单个已处理的轨迹
    """
    print("=== 单个轨迹分析示例 ===")
    
    # 假设的文件路径 (请替换为实际路径)
    processed_trajectory = "/path/to/md_product_1_processed/md_product_1_processed.xtc"
    reference_structure = "/path/to/md_product_1_processed/md.gro"
    output_dir = "/path/to/md_product_1_processed/analysis"
    
    try:
        # 创建分析流水线
        pipeline = create_analysis_pipeline(
            processed_trajectory=processed_trajectory,
            reference_structure=reference_structure,
            output_dir=output_dir
        )
        
        # 运行RMSD分析
        print("运行RMSD分析...")
        rmsd_results = pipeline.run_rmsd_analysis()
        
        # 运行回转半径分析
        print("运行回转半径分析...")
        rg_results = pipeline.run_radius_gyration_analysis()
        
        # 运行距离分析
        print("运行距离分析...")
        distance_results = pipeline.run_distance_analysis()
        
        # 生成报告
        pipeline.generate_analysis_report()
        pipeline.save_results()
        
        print("✅ 分析完成!")
        print(f"📁 结果位置: {output_dir}")
        
    except FileNotFoundError as e:
        print(f"❌ 文件未找到: {e}")
        print("💡 请确保文件路径正确且文件存在")
    except Exception as e:
        print(f"❌ 分析过程出错: {e}")


def example_custom_analysis():
    """
    示例：自定义分析参数
    """
    print("\n=== 自定义分析示例 ===")
    
    # 文件路径
    processed_trajectory = "/path/to/md_product_2_processed/md_product_2_processed.xtc"
    reference_structure = "/path/to/md_product_2_processed/md.gro"
    output_dir = "/path/to/md_product_2_processed/custom_analysis"
    
    try:
        pipeline = create_analysis_pipeline(
            processed_trajectory=processed_trajectory,
            reference_structure=reference_structure,
            output_dir=output_dir
        )
        
        # 自定义RMSD选择
        custom_rmsd_selections = [
            "protein and name CA",                          # 所有α碳
            "protein and name CA and resid 1-50",           # 前50个残基
            "protein and name CA and resname ALA VAL LEU",  # 疏水残基
        ]
        
        # 自定义距离对
        custom_distance_pairs = [
            {
                "name": "N_to_C_terminus",
                "sel1": "protein and resid 1 and name CA",
                "sel2": "protein and resid 200 and name CA",  # 假设有200个残基
                "type": "minimum_distance"
            },
            {
                "name": "domain_separation",
                "sel1": "protein and resid 1-100",
                "sel2": "protein and resid 101-200",
                "type": "center_of_mass"
            }
        ]
        
        # 运行自定义分析
        pipeline.run_rmsd_analysis(selections=custom_rmsd_selections)
        pipeline.run_distance_analysis(distance_pairs=custom_distance_pairs)
        
        print("✅ 自定义分析完成!")
        
    except Exception as e:
        print(f"❌ 自定义分析出错: {e}")


def example_comprehensive_analysis():
    """
    示例：运行所有分析类型
    """
    print("\n=== 综合分析示例 ===")
    
    processed_trajectory = "/path/to/md_product_3_processed/md_product_3_processed.xtc"
    reference_structure = "/path/to/md_product_3_processed/md.gro"
    output_dir = "/path/to/md_product_3_processed/comprehensive_analysis"
    
    try:
        pipeline = create_analysis_pipeline(
            processed_trajectory=processed_trajectory,
            reference_structure=reference_structure,
            output_dir=output_dir
        )
        
        # 运行所有类型的分析
        results = pipeline.run_comprehensive_analysis(
            enable_rmsd=True,
            enable_rg=True,
            enable_distances=True,
            enable_rdf=False,    # RDF分析较慢
            enable_hbonds=False  # 氢键分析较慢
        )
        
        # 打印分析摘要
        summary = results.get("analysis_summary", {})
        print(f"✅ 综合分析完成!")
        print(f"⏱️  总耗时: {summary.get('total_time_seconds', 0):.1f} 秒")
        print(f"📊 完成的分析: {', '.join(summary.get('completed_analyses', []))}")
        print(f"📈 生成图表: {summary.get('total_plots_generated', 0)} 个")
        
    except Exception as e:
        print(f"❌ 综合分析出错: {e}")


def example_batch_analysis_simulation():
    """
    示例：模拟批量分析
    """
    print("\n=== 批量分析模拟 ===")
    
    # 模拟多个已处理的轨迹目录
    processed_dirs = [
        "/path/to/md_products_processed/md_product_1_processed",
        "/path/to/md_products_processed/md_product_2_processed", 
        "/path/to/md_products_processed/md_product_3_processed",
    ]
    
    results_summary = []
    
    for processed_dir in processed_dirs:
        dir_path = Path(processed_dir)
        trajectory_file = dir_path / f"{dir_path.name}.xtc"
        reference_file = dir_path / "md.gro"
        output_dir = dir_path / "analysis"
        
        print(f"\n分析: {dir_path.name}")
        
        try:
            if trajectory_file.exists() and reference_file.exists():
                pipeline = create_analysis_pipeline(
                    processed_trajectory=str(trajectory_file),
                    reference_structure=str(reference_file),
                    output_dir=str(output_dir)
                )
                
                # 运行快速分析 (只包含基本类型)
                results = pipeline.run_comprehensive_analysis(
                    enable_rmsd=True,
                    enable_rg=True,
                    enable_distances=True,
                    enable_rdf=False,
                    enable_hbonds=False
                )
                
                summary = results.get("analysis_summary", {})
                results_summary.append({
                    "name": dir_path.name,
                    "status": "success",
                    "time": summary.get("total_time_seconds", 0),
                    "analyses": len(summary.get("completed_analyses", []))
                })
                
                print(f"   ✅ 成功 (耗时: {summary.get('total_time_seconds', 0):.1f}s)")
                
            else:
                print(f"   ⚠️  跳过: 缺少必要文件")
                results_summary.append({
                    "name": dir_path.name,
                    "status": "skipped",
                    "reason": "missing_files"
                })
                
        except Exception as e:
            print(f"   ❌ 失败: {e}")
            results_summary.append({
                "name": dir_path.name,
                "status": "failed",
                "error": str(e)
            })
    
    # 打印批量分析摘要
    print(f"\n📊 批量分析摘要:")
    successful = [r for r in results_summary if r["status"] == "success"]
    failed = [r for r in results_summary if r["status"] == "failed"]
    skipped = [r for r in results_summary if r["status"] == "skipped"]
    
    print(f"   ✅ 成功: {len(successful)} 个")
    print(f"   ❌ 失败: {len(failed)} 个")
    print(f"   ⚠️  跳过: {len(skipped)} 个")
    
    if successful:
        avg_time = sum(r["time"] for r in successful) / len(successful)
        print(f"   ⏱️  平均耗时: {avg_time:.1f} 秒/轨迹")


def main():
    """
    运行所有示例
    """
    print("🧬 Immunex轨迹分析流水线示例")
    print("=" * 50)
    
    print("\n💡 注意: 这些示例使用虚拟路径")
    print("请将路径替换为实际的文件位置后运行")
    
    # 运行示例 (使用虚拟路径，实际不会执行)
    try:
        example_single_trajectory_analysis()
        example_custom_analysis()
        example_comprehensive_analysis()
        example_batch_analysis_simulation()
        
    except Exception as e:
        print(f"\n💡 示例演示完成 (预期的错误: {e})")
    
    print("\n🎯 实际使用方法:")
    print("1. 替换示例中的路径为真实文件路径")
    print("2. 确保已完成PBC预处理步骤")
    print("3. 运行 python complete_analysis_example.py")
    
    print("\n📚 更多信息:")
    print("- 完整指南: TRAJECTORY_ANALYSIS_GUIDE.md")
    print("- 批量脚本: run_trajectory_analysis.py")
    print("- 工作流程: md_workflow.py")


if __name__ == "__main__":
    main()
