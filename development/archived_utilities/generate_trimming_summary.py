#!/usr/bin/env python3
"""
Generate detailed summary report of trimmed residues.

Analyze all PDB files and create comprehensive report of removed residues.
"""
import sys
sys.path.insert(0, '/home/xumy/work/development/Immunex')

from pathlib import Path
from trim_pdb_by_distance import trim_pdb_by_distance
import json

# Process all files
input_dir = Path("/home/xumy/work/development/Immunex/input/standardizedpdbs/protein_only")
pdb_files = sorted(input_dir.glob("*.pdb"))

print("=" * 80)
print("PDB Trimming Summary Report (8nm cutoff)")
print("=" * 80)
print(f"Processing {len(pdb_files)} PDB files from protein_only/")
print("=" * 80)

# Collect all results
all_results = []
files_with_removals = []

for pdb_file in pdb_files:
    try:
        # Run trimming analysis (without actually saving, just get stats)
        from trim_pdb_by_distance import (
            find_peptide_chain,
            get_peptide_center_residue
        )
        import MDAnalysis as mda
        import numpy as np

        u = mda.Universe(str(pdb_file))
        peptide_segid = find_peptide_chain(u)
        center_residue = get_peptide_center_residue(u, peptide_segid)
        center_com = center_residue.atoms.center_of_mass()

        distance_cutoff_angstrom = 8.0 * 10.0  # 8nm = 80 Angstrom

        all_residues = u.select_atoms("protein").residues
        residues_to_remove = []

        for residue in all_residues:
            residue_com = residue.atoms.center_of_mass()
            distance = np.linalg.norm(residue_com - center_com)

            if distance > distance_cutoff_angstrom:
                residues_to_remove.append({
                    'segid': residue.segment.segid,
                    'resid': residue.resid,
                    'resname': residue.resname,
                    'distance_nm': distance / 10.0
                })

        result = {
            'pdb_file': pdb_file.name,
            'total_residues': len(all_residues),
            'removed_count': len(residues_to_remove),
            'removed_residues': residues_to_remove,
            'peptide_chain': peptide_segid,
            'peptide_center': f"{center_residue.resname}{center_residue.resid}"
        }

        all_results.append(result)

        if len(residues_to_remove) > 0:
            files_with_removals.append(result)

    except Exception as e:
        print(f"ERROR processing {pdb_file.name}: {e}")

# Generate report
print("\n" + "=" * 80)
print("SUMMARY STATISTICS")
print("=" * 80)
total_files = len(all_results)
files_with_trims = len(files_with_removals)
total_residues = sum(r['total_residues'] for r in all_results)
total_removed = sum(r['removed_count'] for r in all_results)

print(f"Total PDB files processed:    {total_files}")
print(f"Files with removed residues:  {files_with_trims} ({files_with_trims/total_files*100:.1f}%)")
print(f"Files without removals:       {total_files - files_with_trims} ({(total_files-files_with_trims)/total_files*100:.1f}%)")
print(f"\nTotal residues analyzed:      {total_residues}")
print(f"Total residues removed:       {total_removed} ({total_removed/total_residues*100:.2f}%)")
print(f"Average retention rate:       {(total_residues-total_removed)/total_residues*100:.2f}%")

if files_with_trims > 0:
    print("\n" + "=" * 80)
    print(f"DETAILED REMOVAL LIST ({files_with_trims} files with removals)")
    print("=" * 80)

    for result in files_with_removals:
        print(f"\n{result['pdb_file']}:")
        print(f"  Peptide: Chain {result['peptide_chain']} (center: {result['peptide_center']})")
        print(f"  Total residues: {result['total_residues']}")
        print(f"  Removed: {result['removed_count']} residues ({result['removed_count']/result['total_residues']*100:.1f}%)")

        # Group by chain
        by_chain = {}
        for res in result['removed_residues']:
            chain = res['segid']
            if chain not in by_chain:
                by_chain[chain] = []
            by_chain[chain].append(res)

        print(f"  Removed residues by chain:")
        for chain in sorted(by_chain.keys()):
            residues = by_chain[chain]
            print(f"    Chain {chain}: {len(residues)} residues")
            for res in residues:
                print(f"      {res['resname']:>3s} {res['resid']:>4d}  (distance: {res['distance_nm']:.2f} nm)")

# Save detailed JSON report
output_json = "/tmp/trimming_summary_protein_only_8nm.json"
with open(output_json, 'w') as f:
    json.dump({
        'summary': {
            'total_files': total_files,
            'files_with_removals': files_with_trims,
            'total_residues': int(total_residues),
            'total_removed': int(total_removed),
            'removal_rate': float(total_removed/total_residues*100)
        },
        'files_with_removals': files_with_removals,
        'all_results': all_results
    }, f, indent=2)

print("\n" + "=" * 80)
print(f"Detailed JSON report saved to: {output_json}")
print("=" * 80)
