#!/usr/bin/env python3
"""
Detailed terminus analysis with flexible terminus definition.
"""
import sys
sys.path.insert(0, '/home/xumy/work/development/Immunex')

from pathlib import Path
import MDAnalysis as mda
import numpy as np
from collections import defaultdict

from development.utilities.trim_pdb_by_distance import (
    find_peptide_chain,
    get_peptide_center_residue
)

input_dir = Path("/home/xumy/work/development/Immunex/input/standardizedpdbs/protein_only")
pdb_files = sorted(input_dir.glob("*.pdb"))

print("=" * 80)
print("Detailed Terminus Analysis")
print("=" * 80)

all_removals = []

for pdb_file in pdb_files:
    try:
        u = mda.Universe(str(pdb_file))
        peptide_segid = find_peptide_chain(u)
        center_residue = get_peptide_center_residue(u, peptide_segid)
        center_com = center_residue.atoms.center_of_mass()

        distance_cutoff_angstrom = 8.0 * 10.0

        chain_ranges = {}
        for seg in u.segments:
            protein_residues = u.select_atoms(f"segid {seg.segid} and protein").residues
            if len(protein_residues) > 0:
                resids = [r.resid for r in protein_residues]
                chain_ranges[seg.segid] = sorted(resids)

        all_residues = u.select_atoms("protein").residues

        for residue in all_residues:
            residue_com = residue.atoms.center_of_mass()
            distance = np.linalg.norm(residue_com - center_com)

            if distance > distance_cutoff_angstrom:
                chain = residue.segment.segid
                resid = residue.resid
                all_resids = chain_ranges[chain]

                pos_idx = all_resids.index(resid)
                total = len(all_resids)

                # Distance from ends
                dist_from_n = pos_idx
                dist_from_c = total - 1 - pos_idx

                all_removals.append({
                    'pdb': pdb_file.name,
                    'chain': chain,
                    'resid': resid,
                    'resname': residue.resname,
                    'distance_nm': distance / 10.0,
                    'position': pos_idx + 1,
                    'total': total,
                    'dist_from_n': dist_from_n,
                    'dist_from_c': dist_from_c
                })

    except Exception as e:
        print(f"ERROR: {pdb_file.name}: {e}")

# Analyze
print(f"\nTotal residues removed: {len(all_removals)}")
print("\n" + "=" * 80)
print("Distance from Chain Termini")
print("=" * 80)

# Group by distance from nearest terminus
terminus_dist = []
for r in all_removals:
    min_dist = min(r['dist_from_n'], r['dist_from_c'])
    terminus_dist.append(min_dist)
    terminus = "N-terminus" if r['dist_from_n'] < r['dist_from_c'] else "C-terminus"
    r['nearest_terminus'] = terminus
    r['terminus_distance'] = min_dist

# Statistics
terminus_dist.sort()
print(f"\nDistance from nearest terminus (residues):")
print(f"  Min:     {min(terminus_dist)}")
print(f"  Max:     {max(terminus_dist)}")
print(f"  Median:  {terminus_dist[len(terminus_dist)//2]}")
print(f"  Mean:    {sum(terminus_dist)/len(terminus_dist):.1f}")

# Histogram
print("\n" + "-" * 80)
print("Distribution:")
print("-" * 80)
bins = [0, 3, 6, 10, 15, 20, 999]
bin_labels = ["0-2 (immediate)", "3-5 (near)", "6-9 (close)", "10-14 (moderate)", "15-19 (far)", "20+ (very far)"]

for i in range(len(bins)-1):
    count = sum(1 for d in terminus_dist if bins[i] <= d < bins[i+1])
    pct = count / len(terminus_dist) * 100
    bar = "█" * int(pct / 2)
    print(f"{bin_labels[i]:20s}: {count:3d} ({pct:5.1f}%) {bar}")

# N vs C terminus
print("\n" + "=" * 80)
print("N-terminus vs C-terminus")
print("=" * 80)
n_term = sum(1 for r in all_removals if r['nearest_terminus'] == 'N-terminus')
c_term = sum(1 for r in all_removals if r['nearest_terminus'] == 'C-terminus')
print(f"N-terminus removals: {n_term} ({n_term/len(all_removals)*100:.1f}%)")
print(f"C-terminus removals: {c_term} ({c_term/len(all_removals)*100:.1f}%)")

# Far from terminus cases
print("\n" + "=" * 80)
print("Cases >10 residues from nearest terminus")
print("=" * 80)
far_cases = [r for r in all_removals if r['terminus_distance'] > 10]
if far_cases:
    print(f"\nFound {len(far_cases)} cases ({len(far_cases)/len(all_removals)*100:.1f}%):")

    # Group by PDB
    by_pdb = defaultdict(list)
    for r in far_cases:
        by_pdb[r['pdb']].append(r)

    for pdb in sorted(by_pdb.keys()):
        cases = by_pdb[pdb]
        print(f"\n{pdb}:")
        for r in cases:
            print(f"  Chain {r['chain']}: {r['resname']} {r['resid']}")
            print(f"    Position: {r['position']}/{r['total']}")
            print(f"    Distance from {r['nearest_terminus']}: {r['terminus_distance']} residues")
            print(f"    Distance from peptide center: {r['distance_nm']:.2f} nm")
else:
    print("\n✓ No cases found!")

print("\n" + "=" * 80)
print("CONCLUSION")
print("=" * 80)
within_10 = sum(1 for d in terminus_dist if d <= 10)
print(f"Removals within 10 residues of terminus: {within_10}/{len(all_removals)} ({within_10/len(all_removals)*100:.1f}%)")
print(f"Removals within 5 residues of terminus:  {sum(1 for d in terminus_dist if d <= 5)}/{len(all_removals)} ({sum(1 for d in terminus_dist if d <= 5)/len(all_removals)*100:.1f}%)")
print(f"Removals within 3 residues of terminus:  {sum(1 for d in terminus_dist if d <= 3)}/{len(all_removals)} ({sum(1 for d in terminus_dist if d <= 3)/len(all_removals)*100:.1f}%)")
print("=" * 80)
