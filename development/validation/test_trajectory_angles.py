#!/usr/bin/env python3
"""
测试角度分析模块 - 实际轨迹分析

使用真实的TCR-pMHC轨迹数据测试重构后的角度分析模块。

Author: Immunex Development Team
Date: 2026-03-19
"""

import sys
from pathlib import Path
import numpy as np

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from immunex.analysis.angles import (
    DockingAngleAnalyzer,
    DockingAngleInput
)


def test_trajectory_angle_analysis():
    """测试真实轨迹的角度分析"""

    print("=" * 70)
    print("TCR-pMHC 对接角度轨迹分析测试")
    print("=" * 70)

    # 查找可用的测试数据
    base_dir = Path("development/workspaces/FEL_workspace")

    # 1bd2数据
    test_cases = [
        {
            'name': '1bd2_rest2_basin1',
            'pdb': base_dir / '1bd2_basin_structures/1bd2_rest2_basin1.pdb',
            'xtc': base_dir / '1bd2_basin_trajectories/1bd2_rest2_basin1_frames0-151.xtc',
        },
        {
            'name': '1ao7_rest2_basin1',
            'pdb': base_dir / '1ao7_basin_structures/1ao7_rest2_basin1.pdb',
            'xtc': base_dir / '1ao7_basin_trajectories/1ao7_rest2_basin1_frames942-1142.xtc',
        }
    ]

    # 查找第一个可用的数据
    selected_case = None
    for case in test_cases:
        if case['pdb'].exists() and case['xtc'].exists():
            selected_case = case
            print(f"\n✓ 找到可用数据: {case['name']}")
            print(f"  PDB: {case['pdb']}")
            print(f"  XTC: {case['xtc']}")
            break

    if not selected_case:
        print("\n✗ 未找到可用的测试数据")
        print("\n可选方案:")
        print("  1. 准备测试数据 (PDB + XTC)")
        print("  2. 使用示例数据进行测试")
        return False

    print("\n" + "=" * 70)
    print("开始角度分析")
    print("=" * 70)

    # 创建输出目录
    output_dir = Path("test_output/angle_analysis")
    output_dir.mkdir(parents=True, exist_ok=True)

    # 创建分析器
    analyzer = DockingAngleAnalyzer()

    # 设置进度回调
    def progress_callback(progress, message):
        print(f"[{progress*100:5.1f}%] {message}")

    analyzer.set_progress_callback(progress_callback)

    # 配置输入参数
    input_params = DockingAngleInput(
        topology=str(selected_case['pdb']),
        trajectory=str(selected_case['xtc']),
        stride=5,  # 每5帧采样一次
        output_dir=str(output_dir),
        auto_identify_chains=True,  # 自动识别链
        use_anarci=True
    )

    print("\n配置参数:")
    print(f"  Topology: {input_params.topology}")
    print(f"  Trajectory: {input_params.trajectory}")
    print(f"  Stride: {input_params.stride}")
    print(f"  Auto identify chains: {input_params.auto_identify_chains}")
    print(f"  Output directory: {input_params.output_dir}")

    # 执行分析
    print("\n" + "-" * 70)
    result = analyzer.analyze(input_params)
    print("-" * 70)

    # 检查结果
    if not result.success:
        print("\n✗ 分析失败")
        print(f"错误信息: {result.error_message}")
        return False

    print("\n" + "=" * 70)
    print("分析结果")
    print("=" * 70)

    # 显示统计信息
    stats = result.statistics
    print(f"\n✓ 分析成功!")
    print(f"\n帧数统计:")
    print(f"  总帧数: {stats['n_frames']}")

    print(f"\nCrossing角 (TCR轴 vs MHC沟槽轴):")
    print(f"  均值: {stats['crossing_mean']:.2f}°")
    print(f"  标准差: {stats['crossing_std']:.2f}°")
    print(f"  范围: [{stats['crossing_min']:.2f}, {stats['crossing_max']:.2f}]°")
    print(f"  波动: {stats['crossing_max'] - stats['crossing_min']:.2f}°")

    print(f"\nIncident角 (TCR轴 vs MHC沟槽法向量):")
    print(f"  均值: {stats['incident_mean']:.2f}°")
    print(f"  标准差: {stats['incident_std']:.2f}°")
    print(f"  范围: [{stats['incident_min']:.2f}, {stats['incident_max']:.2f}]°")
    print(f"  波动: {stats['incident_max'] - stats['incident_min']:.2f}°")

    # 分析角度稳定性
    print(f"\n角度稳定性评估:")
    crossing_cv = stats['crossing_std'] / stats['crossing_mean']
    incident_cv = stats['incident_std'] / stats['incident_mean']

    print(f"  Crossing角 CV: {crossing_cv:.3f} ", end="")
    if crossing_cv < 0.05:
        print("(非常稳定)")
    elif crossing_cv < 0.10:
        print("(稳定)")
    elif crossing_cv < 0.15:
        print("(中等波动)")
    else:
        print("(波动较大)")

    print(f"  Incident角 CV: {incident_cv:.3f} ", end="")
    if incident_cv < 0.05:
        print("(非常稳定)")
    elif incident_cv < 0.10:
        print("(稳定)")
    elif incident_cv < 0.15:
        print("(中等波动)")
    else:
        print("(波动较大)")

    # 输出文件
    print(f"\n输出文件:")
    for f in result.output_files:
        print(f"  ✓ {f}")

    # 显示前几帧的数据
    if result.times is not None and len(result.times) > 0:
        print(f"\n前5帧数据预览:")
        print(f"{'Time(ps)':>10s} {'Crossing(°)':>12s} {'Incident(°)':>12s}")
        print("-" * 36)
        for i in range(min(5, len(result.times))):
            print(f"{result.times[i]:10.1f} {result.crossing_angles[i]:12.2f} {result.incident_angles[i]:12.2f}")

    # 元数据
    if result.metadata:
        print(f"\n分析元数据:")
        print(f"  模块版本: {result.metadata.get('module_version', 'N/A')}")
        print(f"  分析时间: {result.metadata.get('timestamp', 'N/A')}")
        if 'chain_identifications' in result.metadata:
            print(f"  自动识别链: 成功")

    # 绘制简单的ASCII图
    print(f"\n" + "=" * 70)
    print("角度变化趋势 (ASCII图)")
    print("=" * 70)

    if result.crossing_angles is not None and len(result.crossing_angles) > 10:
        plot_ascii_trend(result.times, result.crossing_angles, "Crossing角")
        print()
        plot_ascii_trend(result.times, result.incident_angles, "Incident角")

    print("\n" + "=" * 70)
    print("测试完成!")
    print("=" * 70)

    return True


def plot_ascii_trend(times, angles, title, width=60, height=15):
    """绘制ASCII趋势图"""
    print(f"\n{title} 时间演化:")

    if len(angles) == 0:
        print("  (无数据)")
        return

    # 归一化数据到图表高度
    min_angle = np.min(angles)
    max_angle = np.max(angles)
    range_angle = max_angle - min_angle

    if range_angle == 0:
        print(f"  {min_angle:.2f}° (恒定)")
        return

    # 创建图表
    chart = [[' ' for _ in range(width)] for _ in range(height)]

    # 绘制数据点
    n_points = min(len(angles), width)
    for i in range(n_points):
        idx = int(i * len(angles) / n_points)
        normalized = (angles[idx] - min_angle) / range_angle
        y = int((1 - normalized) * (height - 1))
        x = int(i * width / n_points)
        if 0 <= y < height and 0 <= x < width:
            chart[y][x] = '●'

    # 打印图表
    print(f"  {max_angle:6.2f}° ┤", end="")
    for row in chart:
        print(''.join(row))
        if row != chart[-1]:
            print("           │", end="")
    print(f"  {min_angle:6.2f}° └" + "─" * width)
    print(f"           {times[0]:6.0f}ps" + " " * (width - 20) + f"{times[-1]:6.0f}ps")


if __name__ == '__main__':
    try:
        success = test_trajectory_angle_analysis()
        sys.exit(0 if success else 1)
    except Exception as e:
        print(f"\n✗ 测试过程中出错: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
