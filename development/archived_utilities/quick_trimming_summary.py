"""
Generate clean summary report of 8nm distance trimming.
"""

import sys
sys.path.insert(0, '/home/xumy/work/development/Immunex')

from pathlib import Path
from trim_pdb_by_distance import find_peptide_chain, get_peptide_center_residue
import MDAnalysis as mda
import numpy as np

input_dir = Path("/home/xumy/work/development/Immunex/input/standardizedpdbs/protein_only")
pdb_files = sorted(input_dir.glob("*.pdb"))

all_results = []

for pdb_file in pdb_files:
    try:
        u = mda.Universe(str(pdb_file))
        peptide_segid = find_peptide_chain(u)
        center_residue = get_peptide_center_residue(u, peptide_segid)
        center_com = center_residue.atoms.center_of_mass()

        distance_cutoff_angstrom = 80.0  # 8nm

        all_residues = u.select_atoms("protein").residues
        removed_count = 0

        for residue in all_residues:
            residue_com = residue.atoms.center_of_mass()
            distance = np.linalg.norm(residue_com - center_com)
            if distance > distance_cutoff_angstrom:
                removed_count += 1

        if removed_count > 0:
            all_results.append({
                'pdb_id': pdb_file.stem.replace('_sd', ''),
                'removed': removed_count,
                'total': len(all_residues),
                'percentage': removed_count / len(all_residues) * 100
            })
    except:
        pass

# Sort by removed count
all_results.sort(key=lambda x: x['removed'], reverse=True)

# Print report
print("="*80)
print("8nm距离切割PDB统计报告")
print("="*80)
print(f"\n总共被切割的PDB数量: {len(all_results)}")
print(f"总共切除的残基数: {sum(r['removed'] for r in all_results)}")
print()
print("="*80)
print("被切割的PDB列表 (按切除残基数排序)")
print("="*80)
print(f"\n{'PDB ID':<10} {'切除残基数':>12} {'总残基数':>12} {'切除比例':>12}")
print("-"*80)

for r in all_results:
    print(f"{r['pdb_id']:<10} {r['removed']:>12} {r['total']:>12} {r['percentage']:>11.1f}%")

print()
print("="*80)
print("分类统计")
print("="*80)

bins = [(0, 2), (3, 5), (6, 10), (11, 20), (21, 100)]
for low, high in bins:
    count = len([r for r in all_results if low <= r['removed'] <= high])
    if count > 0:
        print(f"  切除{low:2d}-{high:2d}个残基: {count:2d} 个PDB")
