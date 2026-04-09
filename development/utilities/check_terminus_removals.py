#!/usr/bin/env python3
"""
Check if removed residues are at chain termini.

Analyze whether removed residues are at N-terminus, C-terminus, or middle of chains.
"""
import sys
sys.path.insert(0, '/home/xumy/work/development/Immunex')

from pathlib import Path
import MDAnalysis as mda
import numpy as np
from collections import defaultdict

# Import trimming function components
from development.utilities.trim_pdb_by_distance import (
    find_peptide_chain,
    get_peptide_center_residue
)

input_dir = Path("/home/xumy/work/development/Immunex/input/standardizedpdbs/protein_only")
pdb_files = sorted(input_dir.glob("*.pdb"))

print("=" * 80)
print("Checking if removed residues are at chain termini")
print("=" * 80)

# Collect results
terminus_analysis = []
non_terminus_removals = []

for pdb_file in pdb_files:
    try:
        u = mda.Universe(str(pdb_file))
        peptide_segid = find_peptide_chain(u)
        center_residue = get_peptide_center_residue(u, peptide_segid)
        center_com = center_residue.atoms.center_of_mass()

        distance_cutoff_angstrom = 8.0 * 10.0

        # Get all protein chains and their residue ranges
        chain_ranges = {}
        for seg in u.segments:
            protein_residues = u.select_atoms(f"segid {seg.segid} and protein").residues
            if len(protein_residues) > 0:
                resids = [r.resid for r in protein_residues]
                chain_ranges[seg.segid] = {
                    'min': min(resids),
                    'max': max(resids),
                    'all': sorted(resids)
                }

        # Find removed residues
        all_residues = u.select_atoms("protein").residues
        residues_to_remove = []

        for residue in all_residues:
            residue_com = residue.atoms.center_of_mass()
            distance = np.linalg.norm(residue_com - center_com)

            if distance > distance_cutoff_angstrom:
                chain = residue.segment.segid
                resid = residue.resid

                # Determine position in chain
                chain_info = chain_ranges[chain]
                min_resid = chain_info['min']
                max_resid = chain_info['max']
                all_resids = chain_info['all']

                # Check if at terminus
                is_n_terminus = (resid <= min_resid + 2)  # Within 3 residues of N-terminus
                is_c_terminus = (resid >= max_resid - 2)  # Within 3 residues of C-terminus

                # Check if consecutive with chain end
                is_first = (resid == min_resid)
                is_last = (resid == max_resid)

                # Find position in ordered list
                try:
                    pos_idx = all_resids.index(resid)
                    is_n_end_segment = (pos_idx <= 2)
                    is_c_end_segment = (pos_idx >= len(all_resids) - 3)
                except ValueError:
                    is_n_end_segment = False
                    is_c_end_segment = False

                terminus_type = "MIDDLE"
                if is_first:
                    terminus_type = "N-terminus (first)"
                elif is_last:
                    terminus_type = "C-terminus (last)"
                elif is_n_end_segment:
                    terminus_type = "Near N-terminus"
                elif is_c_end_segment:
                    terminus_type = "Near C-terminus"

                removal_info = {
                    'pdb': pdb_file.name,
                    'chain': chain,
                    'resid': resid,
                    'resname': residue.resname,
                    'distance_nm': distance / 10.0,
                    'chain_min': min_resid,
                    'chain_max': max_resid,
                    'is_terminus': (is_n_end_segment or is_c_end_segment),
                    'terminus_type': terminus_type,
                    'position_from_start': pos_idx + 1,
                    'total_residues': len(all_resids)
                }

                residues_to_remove.append(removal_info)

                if not (is_n_end_segment or is_c_end_segment):
                    non_terminus_removals.append(removal_info)

        if residues_to_remove:
            terminus_analysis.extend(residues_to_remove)

    except Exception as e:
        print(f"ERROR processing {pdb_file.name}: {e}")

# Generate report
print("\n" + "=" * 80)
print("TERMINUS ANALYSIS SUMMARY")
print("=" * 80)

total_removed = len(terminus_analysis)
terminus_removed = sum(1 for r in terminus_analysis if r['is_terminus'])
middle_removed = len(non_terminus_removals)

print(f"Total residues removed:           {total_removed}")
print(f"At terminus (±3 residues):        {terminus_removed} ({terminus_removed/total_removed*100:.1f}%)")
print(f"In middle of chain:               {middle_removed} ({middle_removed/total_removed*100:.1f}%)")

# Categorize by terminus type
terminus_counts = defaultdict(int)
for r in terminus_analysis:
    terminus_counts[r['terminus_type']] += 1

print("\n" + "-" * 80)
print("Breakdown by position:")
print("-" * 80)
for ttype, count in sorted(terminus_counts.items(), key=lambda x: -x[1]):
    print(f"  {ttype:25s}: {count:3d} ({count/total_removed*100:.1f}%)")

if middle_removed > 0:
    print("\n" + "=" * 80)
    print(f"MIDDLE CHAIN REMOVALS (NOT at terminus) - {middle_removed} cases")
    print("=" * 80)

    for r in non_terminus_removals:
        print(f"\n{r['pdb']}:")
        print(f"  Chain {r['chain']}: {r['resname']} {r['resid']}")
        print(f"  Distance: {r['distance_nm']:.2f} nm")
        print(f"  Chain range: {r['chain_min']}-{r['chain_max']} ({r['total_residues']} residues)")
        print(f"  Position: {r['position_from_start']}/{r['total_residues']}")
        print(f"  Status: {r['terminus_type']}")
else:
    print("\n" + "=" * 80)
    print("✓ ALL removed residues are at chain termini!")
    print("=" * 80)

# Detailed breakdown by chain
print("\n" + "=" * 80)
print("REMOVALS BY CHAIN")
print("=" * 80)

chain_removals = defaultdict(list)
for r in terminus_analysis:
    chain_removals[r['chain']].append(r)

for chain in sorted(chain_removals.keys()):
    removals = chain_removals[chain]
    terminus_count = sum(1 for r in removals if r['is_terminus'])
    print(f"\nChain {chain}: {len(removals)} residues removed")
    print(f"  At terminus: {terminus_count} ({terminus_count/len(removals)*100:.1f}%)")

    # Check if all are C-terminus
    c_term = sum(1 for r in removals if 'C-terminus' in r['terminus_type'])
    n_term = sum(1 for r in removals if 'N-terminus' in r['terminus_type'])
    print(f"  C-terminus: {c_term}, N-terminus: {n_term}")

print("\n" + "=" * 80)
