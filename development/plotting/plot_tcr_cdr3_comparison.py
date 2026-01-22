#!/usr/bin/env python3
"""
Compare TCR overall flexibility vs CDR3β flexibility.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy import stats

# Set style
plt.style.use('seaborn-v0_8-paper')
sns.set_palette("husl")

# Paths
comparison_csv = Path("/home/xumy/work/development/AfterMD/output/cdr3_analysis/tcr_cdr3_comparison.csv")
output_dir = Path("/home/xumy/work/development/AfterMD/output/cdr3_analysis/figures")
output_dir.mkdir(exist_ok=True)

# Load data
df = pd.read_csv(comparison_csv)

# Figure 1: Side-by-side distribution comparison
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# CDR3β distribution
axes[0].hist(df['rmsd_mean_nm'], bins=30, alpha=0.6, color='steelblue', edgecolor='black', density=True)
df['rmsd_mean_nm'].plot(kind='kde', ax=axes[0], color='darkblue', linewidth=2)
axes[0].axvline(df['rmsd_mean_nm'].mean(), color='red', linestyle='--', linewidth=2,
                label=f'Mean: {df["rmsd_mean_nm"].mean():.3f} nm')
axes[0].set_xlabel('RMSD (nm)', fontsize=12)
axes[0].set_ylabel('Density', fontsize=12)
axes[0].set_title('CDR3β RMSD Distribution', fontsize=13, fontweight='bold')
axes[0].legend(fontsize=10)
axes[0].grid(alpha=0.3)

# TCR overall distribution
axes[1].hist(df['tcr_rmsd_mean_nm'], bins=30, alpha=0.6, color='coral', edgecolor='black', density=True)
df['tcr_rmsd_mean_nm'].plot(kind='kde', ax=axes[1], color='darkred', linewidth=2)
axes[1].axvline(df['tcr_rmsd_mean_nm'].mean(), color='red', linestyle='--', linewidth=2,
                label=f'Mean: {df["tcr_rmsd_mean_nm"].mean():.3f} nm')
axes[1].set_xlabel('RMSD (nm)', fontsize=12)
axes[1].set_ylabel('Density', fontsize=12)
axes[1].set_title('TCR Overall RMSD Distribution', fontsize=13, fontweight='bold')
axes[1].legend(fontsize=10)
axes[1].grid(alpha=0.3)

plt.tight_layout()
plt.savefig(output_dir / 'tcr_cdr3_distributions.png', dpi=300)
plt.close()

# Figure 2: Overlaid distributions
fig, ax = plt.subplots(figsize=(12, 7))

ax.hist(df['rmsd_mean_nm'], bins=40, alpha=0.5, color='steelblue', label='CDR3β', density=True, edgecolor='black')
ax.hist(df['tcr_rmsd_mean_nm'], bins=40, alpha=0.5, color='coral', label='TCR Overall', density=True, edgecolor='black')

df['rmsd_mean_nm'].plot(kind='kde', ax=ax, color='darkblue', linewidth=2.5, linestyle='--')
df['tcr_rmsd_mean_nm'].plot(kind='kde', ax=ax, color='darkred', linewidth=2.5, linestyle='--')

ax.axvline(df['rmsd_mean_nm'].mean(), color='blue', linestyle='--', linewidth=2, alpha=0.7,
           label=f'CDR3β Mean: {df["rmsd_mean_nm"].mean():.3f} nm')
ax.axvline(df['tcr_rmsd_mean_nm'].mean(), color='red', linestyle='--', linewidth=2, alpha=0.7,
           label=f'TCR Mean: {df["tcr_rmsd_mean_nm"].mean():.3f} nm')

ax.set_xlabel('RMSD (nm)', fontsize=13)
ax.set_ylabel('Density', fontsize=13)
ax.set_title('TCR Overall vs CDR3β Flexibility Comparison', fontsize=14, fontweight='bold')
ax.legend(fontsize=11, loc='upper right')
ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig(output_dir / 'tcr_cdr3_overlay.png', dpi=300)
plt.close()

# Figure 3: Box plot comparison
fig, ax = plt.subplots(figsize=(10, 7))

data_to_plot = [df['rmsd_mean_nm'], df['tcr_rmsd_mean_nm']]
labels = ['CDR3β', 'TCR Overall']

bp = ax.boxplot(data_to_plot, labels=labels, patch_artist=True, notch=True,
                showmeans=True, meanprops=dict(marker='D', markerfacecolor='red', markersize=8))

colors = ['steelblue', 'coral']
for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.7)

ax.set_ylabel('RMSD (nm)', fontsize=13)
ax.set_title('TCR Overall vs CDR3β Flexibility', fontsize=14, fontweight='bold')
ax.grid(axis='y', alpha=0.3)

# Add statistics
ratio = df['tcr_rmsd_mean_nm'].mean() / df['rmsd_mean_nm'].mean()
ax.text(0.98, 0.98, f'Ratio: {ratio:.2f}x', transform=ax.transAxes,
        fontsize=12, verticalalignment='top', horizontalalignment='right',
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

plt.tight_layout()
plt.savefig(output_dir / 'tcr_cdr3_boxplot.png', dpi=300)
plt.close()

# Figure 4: Scatter plot - CDR3β vs TCR overall
fig, ax = plt.subplots(figsize=(10, 10))

ax.scatter(df['rmsd_mean_nm'], df['tcr_rmsd_mean_nm'], alpha=0.6, s=60, color='purple', edgecolor='black')

# Add diagonal reference line (1:1 ratio)
max_val = max(df['rmsd_mean_nm'].max(), df['tcr_rmsd_mean_nm'].max())
ax.plot([0, max_val], [0, max_val], 'k--', linewidth=2, alpha=0.5, label='1:1 Ratio')

# Linear regression
slope, intercept, r_value, p_value, std_err = stats.linregress(df['rmsd_mean_nm'], df['tcr_rmsd_mean_nm'])
x_fit = np.linspace(df['rmsd_mean_nm'].min(), df['rmsd_mean_nm'].max(), 100)
y_fit = slope * x_fit + intercept
ax.plot(x_fit, y_fit, 'r-', linewidth=2.5, label=f'Fit: y={slope:.2f}x+{intercept:.2f}\nR²={r_value**2:.3f}')

ax.set_xlabel('CDR3β RMSD (nm)', fontsize=13)
ax.set_ylabel('TCR Overall RMSD (nm)', fontsize=13)
ax.set_title('Correlation: CDR3β vs TCR Overall Flexibility', fontsize=14, fontweight='bold')
ax.legend(fontsize=11, loc='upper left')
ax.grid(alpha=0.3)
ax.set_aspect('equal', adjustable='box')

plt.tight_layout()
plt.savefig(output_dir / 'tcr_cdr3_scatter.png', dpi=300)
plt.close()

# Figure 5: Ratio distribution
fig, ax = plt.subplots(figsize=(10, 6))

ax.hist(df['rmsd_ratio'], bins=30, alpha=0.7, color='green', edgecolor='black')
ax.axvline(df['rmsd_ratio'].mean(), color='red', linestyle='--', linewidth=2,
           label=f'Mean Ratio: {df["rmsd_ratio"].mean():.2f}')
ax.axvline(df['rmsd_ratio'].median(), color='orange', linestyle='--', linewidth=2,
           label=f'Median Ratio: {df["rmsd_ratio"].median():.2f}')

ax.set_xlabel('TCR/CDR3β RMSD Ratio', fontsize=13)
ax.set_ylabel('Frequency', fontsize=13)
ax.set_title('Distribution of TCR/CDR3β Flexibility Ratio', fontsize=14, fontweight='bold')
ax.legend(fontsize=11)
ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig(output_dir / 'tcr_cdr3_ratio_distribution.png', dpi=300)
plt.close()

# Figure 6: Top 20 most flexible complexes - side by side comparison
fig, ax = plt.subplots(figsize=(14, 10))

# Select top 20 by TCR RMSD
top_20 = df.nlargest(20, 'tcr_rmsd_mean_nm')
y_pos = np.arange(len(top_20))

bar_width = 0.35
ax.barh(y_pos - bar_width/2, top_20['rmsd_mean_nm'], bar_width,
        label='CDR3β', color='steelblue', alpha=0.8)
ax.barh(y_pos + bar_width/2, top_20['tcr_rmsd_mean_nm'], bar_width,
        label='TCR Overall', color='coral', alpha=0.8)

ax.set_yticks(y_pos)
ax.set_yticklabels(top_20['task_name'], fontsize=9)
ax.set_xlabel('RMSD (nm)', fontsize=12)
ax.set_ylabel('Task Name', fontsize=12)
ax.set_title('Top 20 Most Flexible Complexes: CDR3β vs TCR Overall', fontsize=14, fontweight='bold')
ax.legend(fontsize=11, loc='lower right')
ax.grid(axis='x', alpha=0.3)

plt.tight_layout()
plt.savefig(output_dir / 'tcr_cdr3_top20_comparison.png', dpi=300)
plt.close()

# Statistical summary
print("=" * 80)
print("TCR Overall vs CDR3β Flexibility - Statistical Comparison")
print("=" * 80)

print("\nCDR3β Statistics:")
print(f"  Mean: {df['rmsd_mean_nm'].mean():.4f} nm")
print(f"  Std:  {df['rmsd_mean_nm'].std():.4f} nm")
print(f"  Range: {df['rmsd_mean_nm'].min():.4f} - {df['rmsd_mean_nm'].max():.4f} nm")

print("\nTCR Overall Statistics:")
print(f"  Mean: {df['tcr_rmsd_mean_nm'].mean():.4f} nm")
print(f"  Std:  {df['tcr_rmsd_mean_nm'].std():.4f} nm")
print(f"  Range: {df['tcr_rmsd_mean_nm'].min():.4f} - {df['tcr_rmsd_mean_nm'].max():.4f} nm")

print("\nComparison Metrics:")
print(f"  TCR/CDR3β Ratio (Mean): {df['tcr_rmsd_mean_nm'].mean() / df['rmsd_mean_nm'].mean():.2f}x")
print(f"  Correlation (Pearson): R = {r_value:.3f}, R² = {r_value**2:.3f}")
print(f"  P-value: {p_value:.2e}")

# Paired t-test
t_stat, t_pvalue = stats.ttest_rel(df['tcr_rmsd_mean_nm'], df['rmsd_mean_nm'])
print(f"\nPaired t-test:")
print(f"  t-statistic: {t_stat:.3f}")
print(f"  p-value: {t_pvalue:.2e}")
print(f"  Conclusion: TCR overall is {'significantly' if t_pvalue < 0.001 else 'not significantly'} more flexible than CDR3β")

print("\nFigures saved:")
print(f"  1. Side-by-side distributions: {output_dir / 'tcr_cdr3_distributions.png'}")
print(f"  2. Overlaid distributions: {output_dir / 'tcr_cdr3_overlay.png'}")
print(f"  3. Box plot comparison: {output_dir / 'tcr_cdr3_boxplot.png'}")
print(f"  4. Scatter plot (correlation): {output_dir / 'tcr_cdr3_scatter.png'}")
print(f"  5. Ratio distribution: {output_dir / 'tcr_cdr3_ratio_distribution.png'}")
print(f"  6. Top 20 comparison: {output_dir / 'tcr_cdr3_top20_comparison.png'}")

print("=" * 80)
