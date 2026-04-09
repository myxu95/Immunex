#!/usr/bin/env python3
"""
MD Products对比分析运行脚本

收集所有MD Product的分析结果，生成对比图表和统计报告。

使用方法:
    python run_comparative_analysis.py /path/to/processed/simulations
    python run_comparative_analysis.py /path/to/processed/simulations --output ./comparative_results
"""

import sys
import argparse
import logging
import time
from pathlib import Path
from datetime import datetime

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent))

from immunex.analysis.comparative_analysis import create_comparative_analyzer


def setup_logging(log_level: str = "INFO"):
    """设置日志配置"""
    logging.basicConfig(
        level=getattr(logging, log_level.upper()),
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler(),
            logging.FileHandler('comparative_analysis.log')
        ]
    )


def main():
    """主函数"""
    parser = argparse.ArgumentParser(
        description="MD Products对比分析工具 - 生成多MD对比图表和统计报告",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例:
  基本对比分析:
    python run_comparative_analysis.py /path/to/processed/simulations
    
  指定输出目录:
    python run_comparative_analysis.py /path/to/processed/simulations --output ./results
    
  详细日志:
    python run_comparative_analysis.py /path/to/processed/simulations --log-level DEBUG
    
  预览模式:
    python run_comparative_analysis.py /path/to/processed/simulations --dry-run

输出结果:
  - plots/: 对比图表 (PNG格式)
  - data/: 对比数据 (CSV格式)  
  - reports/: 分析报告 (Markdown格式)
        """
    )
    
    # 必需参数
    parser.add_argument(
        "processed_root",
        help="包含已处理MD Products的根目录"
    )
    
    # 输出选项
    parser.add_argument(
        "--output", "-o",
        type=Path,
        default="./comparative_analysis_results",
        help="对比分析结果输出目录 (默认: ./comparative_analysis_results)"
    )
    
    # 日志选项
    parser.add_argument(
        "--log-level",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        default="INFO",
        help="日志级别 (默认: INFO)"
    )
    
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="显示详细信息"
    )
    
    # 其他选项
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="干运行模式，只预览将要分析的MD Products"
    )
    
    args = parser.parse_args()
    
    # 设置日志
    log_level = "DEBUG" if args.verbose else args.log_level
    setup_logging(log_level)
    logger = logging.getLogger(__name__)
    
    # 验证输入目录
    processed_root = Path(args.processed_root)
    if not processed_root.exists():
        print(f"❌ 错误: 目录不存在 - {processed_root}")
        return 1
    
    if not processed_root.is_dir():
        print(f"❌ 错误: 路径不是目录 - {processed_root}")
        return 1
    
    print("🔍 MD Products对比分析")
    print("=" * 50)
    print(f"📂 输入目录: {processed_root}")
    print(f"📁 输出目录: {args.output}")
    print(f"🕐 开始时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    # 预览模式：扫描并显示可用的分析结果
    if args.dry_run:
        print("\n🔍 干运行模式 - 扫描可用的分析结果...")
        
        analysis_dirs = []
        for analysis_dir in processed_root.rglob("analysis"):
            if not analysis_dir.is_dir():
                continue
            
            results_file = analysis_dir / "reports" / "analysis_results.json"
            if results_file.exists():
                md_product_name = analysis_dir.parent.name
                analysis_dirs.append({
                    "name": md_product_name,
                    "path": str(analysis_dir),
                    "size_mb": results_file.stat().st_size / (1024*1024)
                })
        
        if analysis_dirs:
            print(f"\n📋 发现 {len(analysis_dirs)} 个MD Product的分析结果:")
            for i, item in enumerate(analysis_dirs, 1):
                print(f"  {i:2d}. {item['name']}")
                print(f"      路径: {item['path']}")
                print(f"      结果文件: {item['size_mb']:.3f} MB")
            
            print(f"\n📊 将生成的对比图表:")
            print("   • RMSD对比图 (按原子选择类型)")
            print("   • RMSD热力图 (所有MD Products)")
            print("   • 回转半径对比图")
            print("   • 距离分析对比图")
            print("   • 质量统计图表")
            print("   • 综合仪表板")
            
            print(f"\n💾 将保存的数据文件:")
            print("   • rmsd_comparison_data.csv")
            print("   • radius_gyration_comparison_data.csv")
            print("   • distance_comparison_data.csv")
            
            print(f"\n📄 将生成的报告:")
            print("   • comparative_analysis_report.md")
            
        else:
            print("❌ 未找到任何分析结果文件")
            print("💡 确保已运行轨迹分析: python run_trajectory_analysis.py")
        
        return 0
    
    try:
        # 创建对比分析器
        analyzer = create_comparative_analyzer(output_dir=str(args.output))
        
        # 运行对比分析
        print("\n🚀 开始对比分析...")
        start_time = time.time()
        
        results = analyzer.run_comprehensive_comparison(str(processed_root))
        
        end_time = time.time()
        total_time = end_time - start_time
        
        # 显示结果
        if results["status"] == "success":
            print(f"\n✅ 对比分析完成!")
            print("=" * 50)
            print(f"⏱️  总耗时: {total_time:.1f} 秒")
            print(f"📊 分析的MD Products: {results['total_md_products']} 个")
            print(f"🎨 生成的图表: {results['plots_generated']} 个")
            print(f"💾 保存的数据文件: {results['data_files_saved']} 个")
            
            print(f"\n📁 结果位置:")
            print(f"   • 输出目录: {results['output_directory']}")
            print(f"   • 图表目录: {results['output_directory']}/plots/")
            print(f"   • 数据目录: {results['output_directory']}/data/")
            print(f"   • 报告文件: {results['report_file']}")
            
            print(f"\n🎯 主要图表:")
            plot_files = results.get("plot_files", [])
            for plot_file in sorted(plot_files)[:8]:  # 显示前8个
                plot_name = Path(plot_file).name
                print(f"   • {plot_name}")
            if len(plot_files) > 8:
                print(f"   • ... 以及其他 {len(plot_files) - 8} 个图表")
            
            print(f"\n💡 下一步:")
            print(f"   1. 查看综合仪表板: {args.output}/plots/comparative_dashboard.png")
            print(f"   2. 阅读分析报告: {results['report_file']}")
            print(f"   3. 检查详细图表: {args.output}/plots/")
            print(f"   4. 使用数据文件进行进一步分析: {args.output}/data/")
            
            return 0
            
        else:
            print(f"\n❌ 对比分析失败: {results.get('reason', 'unknown')}")
            return 1
    
    except Exception as e:
        print(f"\n❌ 对比分析过程出错: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        return 1


if __name__ == "__main__":
    exit(main())