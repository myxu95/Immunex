#!/usr/bin/env python3
"""
角度数据可视化脚本

使用matplotlib绘制轨迹角度分析结果。
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

# 设置中文字体（如果需要）
plt.rcParams['font.sans-serif'] = ['DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

# 读取数据
data_file = Path("test_output/angle_analysis/docking_angles.csv")
if not data_file.exists():
    print(f"数据文件不存在: {data_file}")
    exit(1)

df = pd.read_csv(data_file)

print("=" * 70)
print("轨迹角度分析 - 数据可视化")
print("=" * 70)
print(f"\n数据文件: {data_file}")
print(f"数据点数: {len(df)}")
print(f"\n数据摘要:")
print(df.describe())

# 创建图表
fig, axes = plt.subplots(3, 2, figsize=(14, 12))
fig.suptitle('TCR-pMHC Docking Angle Trajectory Analysis (1bd2)',
             fontsize=16, fontweight='bold')

# 1. Crossing角时间演化
ax1 = axes[0, 0]
ax1.plot(df['Time(ps)'], df['Crossing(deg)'], 'b-', linewidth=1.5, alpha=0.7)
ax1.axhline(y=df['Crossing(deg)'].mean(), color='r', linestyle='--',
            label=f'Mean: {df["Crossing(deg)"].mean():.2f}°', alpha=0.5)
ax1.fill_between(df['Time(ps)'],
                  df['Crossing(deg)'].mean() - df['Crossing(deg)'].std(),
                  df['Crossing(deg)'].mean() + df['Crossing(deg)'].std(),
                  alpha=0.2, color='r', label=f'±1σ: {df["Crossing(deg)"].std():.2f}°')
ax1.set_xlabel('Time (ps)', fontsize=10)
ax1.set_ylabel('Crossing Angle (°)', fontsize=10)
ax1.set_title('Crossing Angle (TCR axis vs MHC groove axis)', fontsize=11)
ax1.legend(fontsize=9)
ax1.grid(True, alpha=0.3)

# 2. Incident角时间演化
ax2 = axes[0, 1]
ax2.plot(df['Time(ps)'], df['Incident(deg)'], 'g-', linewidth=1.5, alpha=0.7)
ax2.axhline(y=df['Incident(deg)'].mean(), color='r', linestyle='--',
            label=f'Mean: {df["Incident(deg)"].mean():.2f}°', alpha=0.5)
ax2.fill_between(df['Time(ps)'],
                  df['Incident(deg)'].mean() - df['Incident(deg)'].std(),
                  df['Incident(deg)'].mean() + df['Incident(deg)'].std(),
                  alpha=0.2, color='r', label=f'±1σ: {df["Incident(deg)"].std():.2f}°')
ax2.set_xlabel('Time (ps)', fontsize=10)
ax2.set_ylabel('Incident Angle (°)', fontsize=10)
ax2.set_title('Incident Angle (TCR axis vs MHC groove normal)', fontsize=11)
ax2.legend(fontsize=9)
ax2.grid(True, alpha=0.3)

# 3. Crossing角分布直方图
ax3 = axes[1, 0]
ax3.hist(df['Crossing(deg)'], bins=15, color='blue', alpha=0.7, edgecolor='black')
ax3.axvline(x=df['Crossing(deg)'].mean(), color='r', linestyle='--',
            linewidth=2, label='Mean')
ax3.set_xlabel('Crossing Angle (°)', fontsize=10)
ax3.set_ylabel('Frequency', fontsize=10)
ax3.set_title('Crossing Angle Distribution', fontsize=11)
ax3.legend(fontsize=9)
ax3.grid(True, alpha=0.3, axis='y')

# 4. Incident角分布直方图
ax4 = axes[1, 1]
ax4.hist(df['Incident(deg)'], bins=15, color='green', alpha=0.7, edgecolor='black')
ax4.axvline(x=df['Incident(deg)'].mean(), color='r', linestyle='--',
            linewidth=2, label='Mean')
ax4.set_xlabel('Incident Angle (°)', fontsize=10)
ax4.set_ylabel('Frequency', fontsize=10)
ax4.set_title('Incident Angle Distribution', fontsize=11)
ax4.legend(fontsize=9)
ax4.grid(True, alpha=0.3, axis='y')

# 5. 两个角度的相关性散点图
ax5 = axes[2, 0]
scatter = ax5.scatter(df['Crossing(deg)'], df['Incident(deg)'],
                      c=df['Time(ps)'], cmap='viridis', alpha=0.6, s=50)
ax5.set_xlabel('Crossing Angle (°)', fontsize=10)
ax5.set_ylabel('Incident Angle (°)', fontsize=10)
ax5.set_title('Crossing vs Incident Angle Correlation', fontsize=11)
ax5.grid(True, alpha=0.3)
cbar = plt.colorbar(scatter, ax=ax5)
cbar.set_label('Time (ps)', fontsize=9)

# 计算相关系数
corr = df['Crossing(deg)'].corr(df['Incident(deg)'])
ax5.text(0.05, 0.95, f'Correlation: {corr:.3f}',
         transform=ax5.transAxes, fontsize=10,
         verticalalignment='top',
         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

# 6. 滑动窗口标准差（波动性分析）
ax6 = axes[2, 1]
window_size = 5
crossing_rolling_std = df['Crossing(deg)'].rolling(window=window_size).std()
incident_rolling_std = df['Incident(deg)'].rolling(window=window_size).std()

ax6.plot(df['Time(ps)'], crossing_rolling_std, 'b-',
         label='Crossing (rolling std)', linewidth=1.5, alpha=0.7)
ax6.plot(df['Time(ps)'], incident_rolling_std, 'g-',
         label='Incident (rolling std)', linewidth=1.5, alpha=0.7)
ax6.set_xlabel('Time (ps)', fontsize=10)
ax6.set_ylabel('Rolling Std (°)', fontsize=10)
ax6.set_title(f'Angle Fluctuation (window={window_size} frames)', fontsize=11)
ax6.legend(fontsize=9)
ax6.grid(True, alpha=0.3)

plt.tight_layout()

# 保存图片
output_file = Path("test_output/angle_analysis/angle_trajectory_analysis.png")
plt.savefig(output_file, dpi=300, bbox_inches='tight')
print(f"\n✓ 图表已保存: {output_file}")

# 显示图片（如果在支持的环境中）
# plt.show()

print("\n" + "=" * 70)
print("统计分析")
print("=" * 70)

print("\nCrossing角:")
print(f"  均值: {df['Crossing(deg)'].mean():.2f}°")
print(f"  标准差: {df['Crossing(deg)'].std():.2f}°")
print(f"  中位数: {df['Crossing(deg)'].median():.2f}°")
print(f"  范围: [{df['Crossing(deg)'].min():.2f}, {df['Crossing(deg)'].max():.2f}]°")
print(f"  变异系数 (CV): {df['Crossing(deg)'].std() / df['Crossing(deg)'].mean():.3f}")

print("\nIncident角:")
print(f"  均值: {df['Incident(deg)'].mean():.2f}°")
print(f"  标准差: {df['Incident(deg)'].std():.2f}°")
print(f"  中位数: {df['Incident(deg)'].median():.2f}°")
print(f"  范围: [{df['Incident(deg)'].min():.2f}, {df['Incident(deg)'].max():.2f}]°")
print(f"  变异系数 (CV): {df['Incident(deg)'].std() / df['Incident(deg)'].mean():.3f}")

print("\n角度相关性:")
print(f"  Pearson相关系数: {corr:.3f}")

print("\n" + "=" * 70)
print("可视化完成!")
print("=" * 70)
