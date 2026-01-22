#!/usr/bin/env python3
"""
Find the most representative sample from 240 MD tasks based on CDR RMSF patterns.
Representative = closest to the median RMSF values across all 6 CDR loops.
"""

import pandas as pd
import numpy as np
from pathlib import Path

# Read RMSF summary
csv_file = "/home/xumy/work/development/AfterMD/output/rmsf_cdr/batch_whole_protein_rmsf_summary.csv"
df = pd.read_csv(csv_file)

# Filter successful analyses
df = df[df['status'] == 'success'].copy()

print("=" * 80)
print("Finding Most Representative Sample for B-factor Projection")
print("=" * 80)
print(f"\nTotal samples: {df['task'].nunique()}")
print(f"Total CDR records: {len(df)}\n")

# Create CDR label
df['cdr_label'] = df['chain'].str.replace('TCR_', '') + '_' + df['cdr_region']

# Calculate median RMSF for each CDR across all samples
cdr_medians = {}
cdr_labels = ['alpha_CDR1', 'alpha_CDR2', 'alpha_CDR3',
              'beta_CDR1', 'beta_CDR2', 'beta_CDR3']

print("Median RMSF values across all samples:")
for cdr_label in cdr_labels:
    cdr_data = df[df['cdr_label'] == cdr_label]['mean_rmsf_nm']
    median_val = cdr_data.median()
    cdr_medians[cdr_label] = median_val
    print(f"  {cdr_label:15s}: {median_val:.4f} nm (n={len(cdr_data)})")

print("\n" + "=" * 80)
print("Calculating deviation from median for each sample...")
print("=" * 80)

# Calculate deviation for each task
task_deviations = {}

for task_name in df['task'].unique():
    task_data = df[df['task'] == task_name]

    deviations = []
    cdr_values = {}

    for cdr_label in cdr_labels:
        cdr_row = task_data[task_data['cdr_label'] == cdr_label]
        if len(cdr_row) > 0:
            value = cdr_row['mean_rmsf_nm'].values[0]
            cdr_values[cdr_label] = value
            deviation = abs(value - cdr_medians[cdr_label])
            deviations.append(deviation)

    if len(deviations) == 6:  # Has all 6 CDRs
        # Calculate total deviation (L1 norm)
        total_deviation = sum(deviations)
        # Calculate root mean square deviation
        rmsd = np.sqrt(np.mean([d**2 for d in deviations]))

        task_deviations[task_name] = {
            'total_deviation': total_deviation,
            'rmsd': rmsd,
            'mean_deviation': np.mean(deviations),
            'max_deviation': max(deviations),
            'cdr_values': cdr_values
        }

# Sort by RMSD (most representative = smallest RMSD)
sorted_tasks = sorted(task_deviations.items(),
                     key=lambda x: x[1]['rmsd'])

print("\nTop 10 Most Representative Samples (by RMSD from median):")
print("-" * 80)
print(f"{'Rank':<6} {'Task':<20} {'RMSD':>10} {'Mean Dev':>10} {'Max Dev':>10}")
print("-" * 80)

for rank, (task_name, stats) in enumerate(sorted_tasks[:10], 1):
    print(f"{rank:<6} {task_name:<20} {stats['rmsd']:>10.5f} {stats['mean_deviation']:>10.5f} {stats['max_deviation']:>10.5f}")

# Get the most representative sample
best_task = sorted_tasks[0][0]
best_stats = sorted_tasks[0][1]

print("\n" + "=" * 80)
print("MOST REPRESENTATIVE SAMPLE")
print("=" * 80)
print(f"Task name: {best_task}")
print(f"PDB ID: {best_task[:4].upper()}")
print(f"\nDeviation metrics:")
print(f"  RMSD from median: {best_stats['rmsd']:.5f} nm")
print(f"  Mean deviation: {best_stats['mean_deviation']:.5f} nm")
print(f"  Max deviation: {best_stats['max_deviation']:.5f} nm")
print(f"  Total L1 deviation: {best_stats['total_deviation']:.5f} nm")

print(f"\nCDR RMSF values:")
for cdr_label in cdr_labels:
    sample_value = best_stats['cdr_values'][cdr_label]
    median_value = cdr_medians[cdr_label]
    diff = sample_value - median_value
    print(f"  {cdr_label:15s}: {sample_value:.4f} nm (median: {median_value:.4f}, diff: {diff:+.4f})")

# Find the structure file
task_dir = Path(f"/home/xumy/work/development/AfterMD/output/rmsf_cdr/{best_task}")
if task_dir.exists():
    print(f"\nStructure files location:")
    print(f"  Task directory: {task_dir}")

    # Check for PDB and trajectory files
    pdb_candidates = list(task_dir.glob("*.pdb"))
    xtc_candidates = list(task_dir.glob("*.xtc"))

    if pdb_candidates:
        print(f"  PDB file: {pdb_candidates[0].name}")
    if xtc_candidates:
        print(f"  Trajectory: {xtc_candidates[0].name}")

# Also show top 5 for comparison
print("\n" + "=" * 80)
print("Top 5 Representative Samples:")
print("=" * 80)
for rank, (task_name, stats) in enumerate(sorted_tasks[:5], 1):
    pdb_id = task_name[:4].upper()
    print(f"\n{rank}. {task_name} (PDB: {pdb_id})")
    print(f"   RMSD: {stats['rmsd']:.5f} nm")
    print(f"   Mean dev: {stats['mean_deviation']:.5f} nm")

print("\n" + "=" * 80)
print(f"RECOMMENDATION: Use {best_task} for B-factor projection")
print("=" * 80)
