#!/usr/bin/env python3
"""
Analyze and quantify chain-chain contacts in MD trajectory.

This script calculates:
1. Number of contacts between chain pairs over time
2. Contact distance evolution
3. Contact residue identification
4. Generate contact time series plots
"""

import subprocess
from pathlib import Path
import sys
import matplotlib.pyplot as plt
import numpy as np

# Configuration
data_dir = Path("/home/xumy/work/development/Immunex/data/7n1e_1/prod")
processed_dir = data_dir / "pbc_comparison_nojump"  # Use nojump processed trajectory

tpr = processed_dir / "md.tpr"
xtc = processed_dir / "md_processed.xtc"
output_dir = data_dir / "chain_contact_analysis"
output_dir.mkdir(exist_ok=True)

# Contact distance threshold (nm)
CONTACT_CUTOFF = 0.35  # 3.5 Angstrom

print("=" * 80)
print("Chain Contact Analysis")
print("=" * 80)
print(f"Topology: {tpr}")
print(f"Trajectory: {xtc}")
print(f"Contact cutoff: {CONTACT_CUTOFF} nm")
print(f"Output: {output_dir}")
print("=" * 80)


def calculate_mindist(trajectory, topology, sel1, sel2, output_file):
    """
    Calculate minimum distance between two selections over time.

    Args:
        trajectory: Path to trajectory file
        topology: Path to topology file
        sel1: Selection string for group 1
        sel2: Selection string for group 2
        output_file: Output XVG file

    Returns:
        bool: Success status
    """
    cmd = [
        "gmx", "mindist",
        "-f", str(trajectory),
        "-s", str(topology),
        "-od", str(output_file),
        "-on", str(output_file.with_suffix('.num.xvg'))
    ]

    # Create temporary index file with custom selections
    index_file = output_dir / "temp_contacts.ndx"

    # Generate index for selections
    make_ndx_cmd = [
        "gmx", "make_ndx",
        "-f", str(topology),
        "-o", str(index_file)
    ]

    # Commands to create custom groups
    # This depends on your system - adjust chain selections as needed
    ndx_input = f"{sel1}\n{sel2}\nq\n"

    try:
        # Create index file
        result = subprocess.run(
            make_ndx_cmd,
            input=ndx_input,
            text=True,
            capture_output=True,
            timeout=60
        )

        if result.returncode != 0:
            print(f"Warning: make_ndx failed for {sel1} vs {sel2}")
            return False

        # Calculate mindist
        cmd.extend(["-n", str(index_file)])

        result = subprocess.run(
            cmd,
            input="0\n1\n",  # Use the two custom groups
            text=True,
            capture_output=True,
            timeout=120
        )

        if result.returncode != 0:
            print(f"Warning: mindist calculation failed for {sel1} vs {sel2}")
            return False

        # Clean up index file
        if index_file.exists():
            index_file.unlink()

        return True

    except Exception as e:
        print(f"Error calculating mindist: {e}")
        return False


def parse_mindist_xvg(xvg_file):
    """Parse minimum distance XVG file."""
    time = []
    distance = []

    with open(xvg_file, 'r') as f:
        for line in f:
            if line.startswith('#') or line.startswith('@'):
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    time.append(float(parts[0]) / 1000)  # ps to ns
                    distance.append(float(parts[1]))  # nm
                except ValueError:
                    pass

    return np.array(time), np.array(distance)


def count_contacts(distance_array, cutoff):
    """Count frames where distance is below cutoff (in contact)."""
    return np.sum(distance_array < cutoff)


# Define chain pairs to analyze
chain_pairs = [
    ("chain A", "chain B", "A-B"),
    ("chain A", "chain C", "A-C"),
    ("chain B", "chain C", "B-C"),
    ("chain A", "chain D", "A-D (Protein-Peptide)"),
    ("chain A", "chain E", "A-E (Protein-Peptide)"),
    ("chain B", "chain D", "B-D (Protein-Peptide)"),
    ("chain B", "chain E", "B-E (Protein-Peptide)"),
    ("chain C", "chain D", "C-D (Protein-Peptide)"),
    ("chain C", "chain E", "C-E (Protein-Peptide)"),
]

print("\nCalculating pairwise distances...")
print("-" * 80)

distance_data = {}

for chain1, chain2, label in chain_pairs:
    print(f"Analyzing {label}...")

    output_file = output_dir / f"mindist_{label.replace(' ', '_').replace('(', '').replace(')', '')}.xvg"

    # For GROMACS mindist, we use chain selections
    # Note: This is a simplified version - may need adjustment based on your system
    success = calculate_mindist(xtc, tpr, chain1, chain2, output_file)

    if success and output_file.exists():
        time, dist = parse_mindist_xvg(output_file)

        if len(time) > 0:
            distance_data[label] = {
                'time': time,
                'distance': dist,
                'min': np.min(dist),
                'mean': np.mean(dist),
                'max': np.max(dist),
                'contacts': count_contacts(dist, CONTACT_CUTOFF)
            }

            print(f"  Min distance: {distance_data[label]['min']:.3f} nm")
            print(f"  Mean distance: {distance_data[label]['mean']:.3f} nm")
            print(f"  Frames in contact (<{CONTACT_CUTOFF} nm): {distance_data[label]['contacts']}/{len(time)}")

print("\n" + "=" * 80)
print("Summary Statistics")
print("=" * 80)

for label, data in distance_data.items():
    total_frames = len(data['time'])
    contact_percent = (data['contacts'] / total_frames) * 100
    print(f"\n{label}:")
    print(f"  Contact percentage: {contact_percent:.1f}%")
    print(f"  Average distance: {data['mean']:.3f} nm")
    print(f"  Min distance: {data['min']:.3f} nm")

# Generate plots
print("\n" + "=" * 80)
print("Generating plots...")
print("=" * 80)

# Plot 1: All chain pair distances over time
fig, ax = plt.subplots(figsize=(14, 8))

colors = plt.cm.tab10(np.linspace(0, 1, len(distance_data)))

for idx, (label, data) in enumerate(distance_data.items()):
    ax.plot(data['time'], data['distance'], label=label,
            linewidth=1.5, alpha=0.7, color=colors[idx])

ax.axhline(y=CONTACT_CUTOFF, color='red', linestyle='--', linewidth=2,
           label=f'Contact cutoff ({CONTACT_CUTOFF} nm)')
ax.set_xlabel('Time (ns)', fontsize=12, fontweight='bold')
ax.set_ylabel('Minimum Distance (nm)', fontsize=12, fontweight='bold')
ax.set_title('Chain-Chain Minimum Distances Over Time', fontsize=14, fontweight='bold')
ax.legend(loc='best', fontsize=9, ncol=2)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plot1 = output_dir / "chain_distances_all.png"
plt.savefig(plot1, dpi=300, bbox_inches='tight')
plt.close()
print(f"Saved: {plot1}")

# Plot 2: Contact frequency bar chart
fig, ax = plt.subplots(figsize=(12, 6))

labels_list = list(distance_data.keys())
contact_percentages = [(data['contacts'] / len(data['time'])) * 100
                       for data in distance_data.values()]

bars = ax.bar(range(len(labels_list)), contact_percentages,
              color=colors[:len(labels_list)], alpha=0.7, edgecolor='black')

ax.set_xlabel('Chain Pair', fontsize=12, fontweight='bold')
ax.set_ylabel('Contact Percentage (%)', fontsize=12, fontweight='bold')
ax.set_title(f'Chain Contact Frequency (Distance < {CONTACT_CUTOFF} nm)',
             fontsize=14, fontweight='bold')
ax.set_xticks(range(len(labels_list)))
ax.set_xticklabels(labels_list, rotation=45, ha='right')
ax.grid(True, axis='y', alpha=0.3)

# Add percentage labels on bars
for bar, pct in zip(bars, contact_percentages):
    height = bar.get_height()
    ax.text(bar.get_x() + bar.get_width()/2., height,
            f'{pct:.1f}%', ha='center', va='bottom', fontsize=9)

plt.tight_layout()
plot2 = output_dir / "contact_frequency.png"
plt.savefig(plot2, dpi=300, bbox_inches='tight')
plt.close()
print(f"Saved: {plot2}")

# Plot 3: Individual chain pair plots (grid)
n_pairs = len(distance_data)
n_cols = 3
n_rows = (n_pairs + n_cols - 1) // n_cols

fig, axes = plt.subplots(n_rows, n_cols, figsize=(16, n_rows*3))
axes = axes.flatten() if n_pairs > 1 else [axes]

for idx, (label, data) in enumerate(distance_data.items()):
    ax = axes[idx]

    # Plot distance
    ax.plot(data['time'], data['distance'], linewidth=1.5, color=colors[idx])
    ax.axhline(y=CONTACT_CUTOFF, color='red', linestyle='--', linewidth=1.5, alpha=0.7)

    # Shade contact regions
    contact_mask = data['distance'] < CONTACT_CUTOFF
    ax.fill_between(data['time'], 0, CONTACT_CUTOFF, where=contact_mask,
                     alpha=0.3, color='green', label='In contact')

    ax.set_xlabel('Time (ns)', fontsize=10)
    ax.set_ylabel('Distance (nm)', fontsize=10)
    ax.set_title(label, fontsize=11, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend(loc='best', fontsize=8)

# Hide unused subplots
for idx in range(n_pairs, len(axes)):
    axes[idx].axis('off')

plt.tight_layout()
plot3 = output_dir / "chain_distances_individual.png"
plt.savefig(plot3, dpi=300, bbox_inches='tight')
plt.close()
print(f"Saved: {plot3}")

print("\n" + "=" * 80)
print("Analysis completed!")
print("=" * 80)
print(f"\nOutput directory: {output_dir}")
print(f"\nGenerated plots:")
print(f"  1. {plot1.name} - All chain distances")
print(f"  2. {plot2.name} - Contact frequency")
print(f"  3. {plot3.name} - Individual chain pairs")
print("=" * 80)
