#!/usr/bin/env python3
"""
Visualize CDR3 beta flexibility analysis results.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import json

# Set style
plt.style.use('seaborn-v0_8-paper')
sns.set_palette("husl")

# Paths
summary_csv = Path("/home/xumy/work/development/Immunex/output/cdr3_analysis/cdr3_beta_rmsd_complete_summary.csv")
output_dir = Path("/home/xumy/work/development/Immunex/output/cdr3_analysis/figures")
output_dir.mkdir(exist_ok=True)

# Load data
df = pd.read_csv(summary_csv)

# Add stability category
def categorize_stability(rmsd):
    if rmsd < 0.10:
        return 'Very Stable'
    elif rmsd < 0.15:
        return 'Stable'
    elif rmsd < 0.20:
        return 'Moderate'
    elif rmsd < 0.25:
        return 'Flexible'
    else:
        return 'Highly Flexible'

df['stability'] = df['rmsd_mean_nm'].apply(categorize_stability)

# Define stability order
stability_order = ['Very Stable', 'Stable', 'Moderate', 'Flexible', 'Highly Flexible']
df['stability'] = pd.Categorical(df['stability'], categories=stability_order, ordered=True)

# Figure 1: Overall RMSD distribution histogram + KDE
fig, ax = plt.subplots(figsize=(10, 6))
ax.hist(df['rmsd_mean_nm'], bins=30, alpha=0.6, color='steelblue', edgecolor='black', density=True, label='Histogram')
df['rmsd_mean_nm'].plot(kind='kde', ax=ax, color='darkred', linewidth=2, label='KDE')
ax.axvline(df['rmsd_mean_nm'].mean(), color='red', linestyle='--', linewidth=2, label=f'Mean: {df["rmsd_mean_nm"].mean():.3f} nm')
ax.axvline(df['rmsd_mean_nm'].median(), color='green', linestyle='--', linewidth=2, label=f'Median: {df["rmsd_mean_nm"].median():.3f} nm')
ax.set_xlabel('RMSD (nm)', fontsize=12)
ax.set_ylabel('Density', fontsize=12)
ax.set_title('CDR3β RMSD Distribution (N=191)', fontsize=14, fontweight='bold')
ax.legend(fontsize=10)
ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(output_dir / 'cdr3_rmsd_distribution.png', dpi=300)
plt.close()

# Figure 2: Stability category pie chart
fig, ax = plt.subplots(figsize=(8, 8))
stability_counts = df['stability'].value_counts().reindex(stability_order)
colors = sns.color_palette("RdYlGn_r", len(stability_order))
wedges, texts, autotexts = ax.pie(
    stability_counts,
    labels=stability_counts.index,
    autopct='%1.1f%%',
    startangle=90,
    colors=colors,
    textprops={'fontsize': 11}
)
for autotext in autotexts:
    autotext.set_color('white')
    autotext.set_fontweight('bold')
ax.set_title('CDR3β Stability Classification', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig(output_dir / 'cdr3_stability_pie.png', dpi=300)
plt.close()

# Figure 3: Box plot by stability category
fig, ax = plt.subplots(figsize=(10, 6))
df_sorted = df.sort_values('stability')
bp = ax.boxplot(
    [df_sorted[df_sorted['stability'] == cat]['rmsd_mean_nm'].values for cat in stability_order],
    labels=stability_order,
    patch_artist=True,
    notch=True,
    showmeans=True,
    meanprops=dict(marker='D', markerfacecolor='red', markersize=6)
)
for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.7)
ax.set_xlabel('Stability Category', fontsize=12)
ax.set_ylabel('RMSD (nm)', fontsize=12)
ax.set_title('CDR3β RMSD Distribution by Stability Category', fontsize=14, fontweight='bold')
ax.grid(axis='y', alpha=0.3)
plt.xticks(rotation=15, ha='right')
plt.tight_layout()
plt.savefig(output_dir / 'cdr3_rmsd_boxplot.png', dpi=300)
plt.close()

# Figure 4: Violin plot
fig, ax = plt.subplots(figsize=(12, 6))
parts = ax.violinplot(
    [df_sorted[df_sorted['stability'] == cat]['rmsd_mean_nm'].values for cat in stability_order],
    positions=range(len(stability_order)),
    showmeans=True,
    showmedians=True
)
for pc, color in zip(parts['bodies'], colors):
    pc.set_facecolor(color)
    pc.set_alpha(0.7)
ax.set_xticks(range(len(stability_order)))
ax.set_xticklabels(stability_order, rotation=15, ha='right')
ax.set_xlabel('Stability Category', fontsize=12)
ax.set_ylabel('RMSD (nm)', fontsize=12)
ax.set_title('CDR3β RMSD Violin Plot by Stability Category', fontsize=14, fontweight='bold')
ax.grid(axis='y', alpha=0.3)
plt.tight_layout()
plt.savefig(output_dir / 'cdr3_rmsd_violin.png', dpi=300)
plt.close()

# Figure 5: Representative time series
representative_tasks = {
    'Most Stable': df.nsmallest(1, 'rmsd_mean_nm').iloc[0]['task_name'],
    'Moderate': df[df['stability'] == 'Moderate'].iloc[len(df[df['stability'] == 'Moderate'])//2]['task_name'],
    'Most Flexible': df.nlargest(1, 'rmsd_mean_nm').iloc[0]['task_name']
}

fig, axes = plt.subplots(3, 1, figsize=(12, 10))
trajectory_base = Path("/home/xumy/work/development/Immunex/input/pbc_1000frames_2step")

for ax, (label, task_name) in zip(axes, representative_tasks.items()):
    rmsd_file = trajectory_base / task_name / "cdr3_beta_rmsd.xvg"

    if rmsd_file.exists():
        # Parse XVG file
        time_data, rmsd_data = [], []
        with open(rmsd_file, 'r') as f:
            for line in f:
                if not line.startswith(('#', '@')):
                    parts = line.split()
                    if len(parts) >= 2:
                        time_data.append(float(parts[0]))
                        rmsd_data.append(float(parts[1]))

        ax.plot(time_data, rmsd_data, linewidth=1, alpha=0.8)
        ax.axhline(np.mean(rmsd_data), color='red', linestyle='--', linewidth=1.5, label=f'Mean: {np.mean(rmsd_data):.3f} nm')
        ax.fill_between(time_data,
                        np.mean(rmsd_data) - np.std(rmsd_data),
                        np.mean(rmsd_data) + np.std(rmsd_data),
                        alpha=0.2, color='gray', label=f'Std: {np.std(rmsd_data):.3f} nm')
        ax.set_ylabel('RMSD (nm)', fontsize=11)
        ax.set_title(f'{label}: {task_name}', fontsize=12, fontweight='bold')
        ax.legend(fontsize=9, loc='upper right')
        ax.grid(alpha=0.3)

axes[-1].set_xlabel('Time (ns)', fontsize=11)
plt.tight_layout()
plt.savefig(output_dir / 'cdr3_rmsd_timeseries.png', dpi=300)
plt.close()

# Figure 6: RMSD range plot (min-max-mean)
fig, ax = plt.subplots(figsize=(14, 8))
top_20 = df.nlargest(20, 'rmsd_mean_nm')
y_pos = np.arange(len(top_20))

ax.errorbar(
    top_20['rmsd_mean_nm'],
    y_pos,
    xerr=[top_20['rmsd_mean_nm'] - top_20['rmsd_min_nm'],
          top_20['rmsd_max_nm'] - top_20['rmsd_mean_nm']],
    fmt='o',
    markersize=8,
    capsize=5,
    capthick=2,
    elinewidth=2,
    alpha=0.7
)

ax.set_yticks(y_pos)
ax.set_yticklabels(top_20['task_name'])
ax.set_xlabel('RMSD (nm)', fontsize=12)
ax.set_ylabel('Task Name', fontsize=12)
ax.set_title('Top 20 Most Flexible CDR3β Complexes (RMSD Range)', fontsize=14, fontweight='bold')
ax.grid(axis='x', alpha=0.3)
plt.tight_layout()
plt.savefig(output_dir / 'cdr3_rmsd_range_top20.png', dpi=300)
plt.close()

# Generate summary statistics table
summary_stats = pd.DataFrame({
    'Statistic': ['Count', 'Mean', 'Std', 'Median', 'Min', 'Max', 'Q1', 'Q3', 'IQR'],
    'Value (nm)': [
        len(df),
        df['rmsd_mean_nm'].mean(),
        df['rmsd_mean_nm'].std(),
        df['rmsd_mean_nm'].median(),
        df['rmsd_mean_nm'].min(),
        df['rmsd_mean_nm'].max(),
        df['rmsd_mean_nm'].quantile(0.25),
        df['rmsd_mean_nm'].quantile(0.75),
        df['rmsd_mean_nm'].quantile(0.75) - df['rmsd_mean_nm'].quantile(0.25)
    ]
})
summary_stats.to_csv(output_dir / 'cdr3_rmsd_summary_stats.csv', index=False)

print("=" * 80)
print("CDR3β Flexibility Visualization Complete")
print("=" * 80)
print(f"\nGenerated figures:")
print(f"  1. RMSD distribution: {output_dir / 'cdr3_rmsd_distribution.png'}")
print(f"  2. Stability pie chart: {output_dir / 'cdr3_stability_pie.png'}")
print(f"  3. RMSD box plot: {output_dir / 'cdr3_rmsd_boxplot.png'}")
print(f"  4. RMSD violin plot: {output_dir / 'cdr3_rmsd_violin.png'}")
print(f"  5. Representative time series: {output_dir / 'cdr3_rmsd_timeseries.png'}")
print(f"  6. Top 20 flexible complexes: {output_dir / 'cdr3_rmsd_range_top20.png'}")
print(f"\nSummary statistics saved to: {output_dir / 'cdr3_rmsd_summary_stats.csv'}")
print("=" * 80)
