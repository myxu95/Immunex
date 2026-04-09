#!/usr/bin/env python3
"""
Visualize TCR overall RMSD distribution for trajectory quality assessment.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Set style
plt.style.use('seaborn-v0_8-paper')
sns.set_palette("husl")

# Paths
tcr_summary = Path("/home/xumy/work/development/Immunex/output/cdr3_analysis/tcr_overall_rmsd_summary.csv")
output_dir = Path("/home/xumy/work/development/Immunex/output/cdr3_analysis/figures")

# Load data
df = pd.read_csv(tcr_summary)

# Calculate statistics
mean_rmsd = df['tcr_rmsd_mean_nm'].mean()
std_rmsd = df['tcr_rmsd_mean_nm'].std()
median_rmsd = df['tcr_rmsd_mean_nm'].median()
q75 = df['tcr_rmsd_mean_nm'].quantile(0.75)
q95 = df['tcr_rmsd_mean_nm'].quantile(0.95)

# Create comprehensive figure
fig = plt.figure(figsize=(16, 10))
gs = fig.add_gridspec(3, 2, hspace=0.3, wspace=0.3)

# 1. Histogram + KDE
ax1 = fig.add_subplot(gs[0, :])
ax1.hist(df['tcr_rmsd_mean_nm'], bins=50, alpha=0.6, color='coral', edgecolor='black', density=True)
df['tcr_rmsd_mean_nm'].plot(kind='kde', ax=ax1, color='darkred', linewidth=2.5)
ax1.axvline(mean_rmsd, color='red', linestyle='--', linewidth=2, label=f'Mean: {mean_rmsd:.3f} nm')
ax1.axvline(median_rmsd, color='green', linestyle='--', linewidth=2, label=f'Median: {median_rmsd:.3f} nm')
ax1.axvline(q75, color='orange', linestyle='--', linewidth=2, alpha=0.7, label=f'Q75: {q75:.3f} nm')
ax1.axvline(q95, color='purple', linestyle='--', linewidth=2, alpha=0.7, label=f'Q95: {q95:.3f} nm')
ax1.axvline(mean_rmsd + 2*std_rmsd, color='blue', linestyle=':', linewidth=2, alpha=0.7,
            label=f'Mean+2σ: {mean_rmsd + 2*std_rmsd:.3f} nm')
ax1.set_xlabel('TCR Overall RMSD (nm)', fontsize=13)
ax1.set_ylabel('Density', fontsize=13)
ax1.set_title('TCR Overall RMSD Distribution (pHLA-aligned)', fontsize=14, fontweight='bold')
ax1.legend(fontsize=10, loc='upper right')
ax1.grid(alpha=0.3)

# 2. Box plot with outliers
ax2 = fig.add_subplot(gs[1, 0])
bp = ax2.boxplot([df['tcr_rmsd_mean_nm']], vert=False, patch_artist=True,
                  showmeans=True, meanprops=dict(marker='D', markerfacecolor='red', markersize=8),
                  widths=0.5)
bp['boxes'][0].set_facecolor('coral')
bp['boxes'][0].set_alpha(0.7)
ax2.set_xlabel('RMSD (nm)', fontsize=12)
ax2.set_yticklabels(['TCR Overall'])
ax2.set_title('Box Plot with Outliers', fontsize=13, fontweight='bold')
ax2.grid(axis='x', alpha=0.3)

# 3. Cumulative distribution
ax3 = fig.add_subplot(gs[1, 1])
sorted_rmsd = np.sort(df['tcr_rmsd_mean_nm'])
cumulative = np.arange(1, len(sorted_rmsd) + 1) / len(sorted_rmsd) * 100
ax3.plot(sorted_rmsd, cumulative, linewidth=2.5, color='darkred')
ax3.axhline(50, color='green', linestyle='--', linewidth=1.5, alpha=0.7, label='50%')
ax3.axhline(75, color='orange', linestyle='--', linewidth=1.5, alpha=0.7, label='75%')
ax3.axhline(95, color='purple', linestyle='--', linewidth=1.5, alpha=0.7, label='95%')
ax3.axvline(median_rmsd, color='green', linestyle='--', linewidth=1.5, alpha=0.7)
ax3.axvline(q75, color='orange', linestyle='--', linewidth=1.5, alpha=0.7)
ax3.axvline(q95, color='purple', linestyle='--', linewidth=1.5, alpha=0.7)
ax3.set_xlabel('TCR Overall RMSD (nm)', fontsize=12)
ax3.set_ylabel('Cumulative Percentage (%)', fontsize=12)
ax3.set_title('Cumulative Distribution', fontsize=13, fontweight='bold')
ax3.legend(fontsize=10)
ax3.grid(alpha=0.3)

# 4. Top 20 highest RMSD
ax4 = fig.add_subplot(gs[2, :])
top_20 = df.nlargest(20, 'tcr_rmsd_mean_nm')
colors_top20 = ['red' if x > q95 else 'orange' if x > q75 else 'green'
                for x in top_20['tcr_rmsd_mean_nm']]
bars = ax4.barh(range(len(top_20)), top_20['tcr_rmsd_mean_nm'], color=colors_top20, alpha=0.7, edgecolor='black')
ax4.set_yticks(range(len(top_20)))
ax4.set_yticklabels(top_20['task_name'], fontsize=9)
ax4.set_xlabel('TCR Overall RMSD (nm)', fontsize=12)
ax4.set_title('Top 20 Highest TCR RMSD (Potential Low Quality Trajectories)', fontsize=13, fontweight='bold')
ax4.axvline(q95, color='purple', linestyle='--', linewidth=2, alpha=0.7, label=f'95th percentile: {q95:.3f} nm')
ax4.axvline(mean_rmsd + 2*std_rmsd, color='blue', linestyle=':', linewidth=2, alpha=0.7,
            label=f'Mean+2σ: {mean_rmsd + 2*std_rmsd:.3f} nm')
ax4.legend(fontsize=10, loc='lower right')
ax4.grid(axis='x', alpha=0.3)
ax4.invert_yaxis()

plt.suptitle('TCR Overall RMSD Distribution Analysis\n(Aligned to pHLA, N=187 tasks)',
             fontsize=16, fontweight='bold', y=0.995)

plt.savefig(output_dir / 'tcr_overall_rmsd_distribution.png', dpi=300, bbox_inches='tight')
plt.close()

# Summary statistics
print("=" * 80)
print("TCR Overall RMSD Distribution Summary")
print("=" * 80)
print(f"\nAlignment: pHLA (reference)")
print(f"Measurement: TCR overall structure")
print(f"Total tasks: {len(df)}")

print(f"\nBasic Statistics:")
print(f"  Mean:   {mean_rmsd:.4f} nm")
print(f"  Median: {median_rmsd:.4f} nm")
print(f"  Std:    {std_rmsd:.4f} nm")
print(f"  Min:    {df['tcr_rmsd_mean_nm'].min():.4f} nm")
print(f"  Max:    {df['tcr_rmsd_mean_nm'].max():.4f} nm")

print(f"\nPercentiles:")
print(f"  25th (Q1):  {df['tcr_rmsd_mean_nm'].quantile(0.25):.4f} nm")
print(f"  50th (Q2):  {median_rmsd:.4f} nm")
print(f"  75th (Q3):  {q75:.4f} nm")
print(f"  95th:       {q95:.4f} nm")

print(f"\nQuality Thresholds (suggestions):")
print(f"  Conservative (Q75):      < {q75:.3f} nm  ({(df['tcr_rmsd_mean_nm'] < q75).sum()} tasks)")
print(f"  Moderate (Q95):          < {q95:.3f} nm  ({(df['tcr_rmsd_mean_nm'] < q95).sum()} tasks)")
print(f"  Lenient (Mean + 2σ):     < {mean_rmsd + 2*std_rmsd:.3f} nm  ({(df['tcr_rmsd_mean_nm'] < mean_rmsd + 2*std_rmsd).sum()} tasks)")

# Identify outliers
outliers_2sigma = df[df['tcr_rmsd_mean_nm'] > mean_rmsd + 2*std_rmsd]
outliers_q95 = df[df['tcr_rmsd_mean_nm'] > q95]

print(f"\nPotential Low Quality Trajectories:")
print(f"  Beyond Mean+2σ: {len(outliers_2sigma)} tasks")
if len(outliers_2sigma) > 0:
    print("    Tasks:", ", ".join(outliers_2sigma.nsmallest(10, 'tcr_rmsd_mean_nm')['task_name'].tolist()))

print(f"\n  Beyond 95th percentile: {len(outliers_q95)} tasks")
if len(outliers_q95) > 0:
    print("    Tasks:", ", ".join(outliers_q95.nsmallest(10, 'tcr_rmsd_mean_nm')['task_name'].tolist()))

print(f"\n✓ Figure saved to: {output_dir / 'tcr_overall_rmsd_distribution.png'}")
print("=" * 80)
