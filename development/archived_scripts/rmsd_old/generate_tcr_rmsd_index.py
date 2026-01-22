#!/usr/bin/env python3
"""
Generate sorted task index by TCR overall RMSD.
Identify and categorize tasks by quality levels.
"""

import pandas as pd
from pathlib import Path

# Load data
tcr_summary = Path("/home/xumy/work/development/AfterMD/output/cdr3_analysis/tcr_overall_rmsd_summary.csv")
output_dir = Path("/home/xumy/work/development/AfterMD/output/cdr3_analysis")

df = pd.read_csv(tcr_summary)

# Sort by RMSD descending
df_sorted = df.sort_values('tcr_rmsd_mean_nm', ascending=False).reset_index(drop=True)

# Calculate thresholds
mean_rmsd = df['tcr_rmsd_mean_nm'].mean()
std_rmsd = df['tcr_rmsd_mean_nm'].std()
q75 = df['tcr_rmsd_mean_nm'].quantile(0.75)
q95 = df['tcr_rmsd_mean_nm'].quantile(0.95)
threshold_2sigma = mean_rmsd + 2 * std_rmsd

# Categorize quality
def categorize_quality(rmsd):
    if rmsd > threshold_2sigma:
        return 'Extreme Outlier'
    elif rmsd > q95:
        return 'High RMSD'
    elif rmsd > q75:
        return 'Above Average'
    else:
        return 'Normal'

df_sorted['quality_category'] = df_sorted['tcr_rmsd_mean_nm'].apply(categorize_quality)

# Save full sorted list
full_output = output_dir / "tcr_rmsd_sorted_tasks.csv"
df_sorted[['task_name', 'tcr_rmsd_mean_nm', 'tcr_rmsd_std_nm', 'quality_category']].to_csv(
    full_output, index=True, index_label='rank'
)

# Generate summary report
print("=" * 80)
print("TCR Overall RMSD - Sorted Task Index")
print("=" * 80)

print(f"\nTotal tasks: {len(df_sorted)}")
print(f"\nQuality Thresholds:")
print(f"  Q75:           {q75:.3f} nm")
print(f"  Q95:           {q95:.3f} nm")
print(f"  Mean + 2σ:     {threshold_2sigma:.3f} nm")

# Category counts
category_counts = df_sorted['quality_category'].value_counts()
print(f"\nQuality Categories:")
for cat in ['Extreme Outlier', 'High RMSD', 'Above Average', 'Normal']:
    if cat in category_counts:
        count = category_counts[cat]
        pct = count / len(df_sorted) * 100
        print(f"  {cat:20s}: {count:3d} tasks ({pct:5.1f}%)")

# Extreme outliers (Mean + 2σ)
print("\n" + "=" * 80)
print("EXTREME OUTLIERS (RMSD > Mean+2σ = {:.3f} nm)".format(threshold_2sigma))
print("=" * 80)
extreme = df_sorted[df_sorted['quality_category'] == 'Extreme Outlier']
if len(extreme) > 0:
    print(f"\nTotal: {len(extreme)} tasks\n")
    for idx, row in extreme.iterrows():
        print(f"  Rank {idx+1:3d}. {row['task_name']:15s}  RMSD: {row['tcr_rmsd_mean_nm']:.4f} nm  (±{row['tcr_rmsd_std_nm']:.4f})")
else:
    print("  None")

# High RMSD (Q95 - Mean+2σ)
print("\n" + "=" * 80)
print("HIGH RMSD ({:.3f} < RMSD ≤ {:.3f} nm)".format(q95, threshold_2sigma))
print("=" * 80)
high = df_sorted[df_sorted['quality_category'] == 'High RMSD']
if len(high) > 0:
    print(f"\nTotal: {len(high)} tasks\n")
    for idx, row in high.iterrows():
        print(f"  Rank {idx+1:3d}. {row['task_name']:15s}  RMSD: {row['tcr_rmsd_mean_nm']:.4f} nm  (±{row['tcr_rmsd_std_nm']:.4f})")
else:
    print("  None")

# Above average (Q75 - Q95)
print("\n" + "=" * 80)
print("ABOVE AVERAGE ({:.3f} < RMSD ≤ {:.3f} nm)".format(q75, q95))
print("=" * 80)
above_avg = df_sorted[df_sorted['quality_category'] == 'Above Average']
print(f"\nTotal: {len(above_avg)} tasks")
print("\nShowing top 20:")
for idx, row in above_avg.head(20).iterrows():
    print(f"  Rank {idx+1:3d}. {row['task_name']:15s}  RMSD: {row['tcr_rmsd_mean_nm']:.4f} nm  (±{row['tcr_rmsd_std_nm']:.4f})")

# Normal range
print("\n" + "=" * 80)
print("NORMAL RANGE (RMSD ≤ {:.3f} nm)".format(q75))
print("=" * 80)
normal = df_sorted[df_sorted['quality_category'] == 'Normal']
print(f"\nTotal: {len(normal)} tasks")
print("\nShowing top 10 most stable:")
for idx, row in normal.tail(10).iterrows():
    print(f"  Rank {idx+1:3d}. {row['task_name']:15s}  RMSD: {row['tcr_rmsd_mean_nm']:.4f} nm  (±{row['tcr_rmsd_std_nm']:.4f})")

# Save filtered lists
extreme_output = output_dir / "tcr_rmsd_extreme_outliers.txt"
high_output = output_dir / "tcr_rmsd_high_rmsd.txt"

with open(extreme_output, 'w') as f:
    f.write("# Extreme Outlier Tasks (RMSD > {:.3f} nm)\n".format(threshold_2sigma))
    f.write("# These tasks should be excluded from analysis\n\n")
    for _, row in extreme.iterrows():
        f.write(f"{row['task_name']}\n")

with open(high_output, 'w') as f:
    f.write("# High RMSD Tasks ({:.3f} < RMSD ≤ {:.3f} nm)\n".format(q95, threshold_2sigma))
    f.write("# Consider excluding these tasks for stringent quality control\n\n")
    for _, row in high.iterrows():
        f.write(f"{row['task_name']}\n")

print("\n" + "=" * 80)
print("Files Generated")
print("=" * 80)
print(f"  Full sorted list:     {full_output}")
print(f"  Extreme outliers:     {extreme_output} ({len(extreme)} tasks)")
print(f"  High RMSD tasks:      {high_output} ({len(high)} tasks)")
print("=" * 80)
