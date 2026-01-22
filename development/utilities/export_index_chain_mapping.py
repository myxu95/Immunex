#!/usr/bin/env python3
"""
Export chain and residue mapping for index groups.
"""

import subprocess
from pathlib import Path
from collections import defaultdict


def extract_pdb_from_tpr(task_dir, output_pdb):
    """Extract structure from trajectory to PDB - protein only."""
    tpr_file = task_dir / "md.tpr"

    # Find trajectory file
    xtc_files = list(task_dir.glob("*_processed.xtc"))
    if not xtc_files:
        return False

    traj_file = xtc_files[0]

    cmd = f"echo 'Protein' | gmx trjconv -s {tpr_file} -f {traj_file} -o {output_pdb} -dump 0 2>&1"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    return output_pdb.exists()


def read_pdb_atoms(pdb_file):
    """Read atom information from PDB file."""
    atoms = {}

    with open(pdb_file) as f:
        for line in f:
            if line.startswith("ATOM"):
                atom_num = int(line[6:11].strip())
                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()
                chain_id = line[21:22]
                res_num = int(line[22:26].strip())

                atoms[atom_num] = {
                    'atom_name': atom_name,
                    'res_name': res_name,
                    'chain_id': chain_id,
                    'res_num': res_num
                }

    return atoms


def read_index_group(index_file, group_name):
    """Read atom numbers from index group."""
    with open(index_file) as f:
        content = f.read()

    if f"[ {group_name} ]" not in content:
        return None

    # Extract atom numbers
    section = content.split(f"[ {group_name} ]")[1].split("[")[0].strip()
    atom_numbers = []
    for line in section.split("\n"):
        atom_numbers.extend([int(x) for x in line.split()])

    return atom_numbers


def analyze_group_mapping(atoms, atom_numbers):
    """Analyze which chains and residues are included in a group."""
    chain_info = defaultdict(lambda: {'residues': set(), 'atoms': []})

    for atom_num in atom_numbers:
        if atom_num in atoms:
            atom_data = atoms[atom_num]
            chain_id = atom_data['chain_id']
            res_num = atom_data['res_num']
            res_name = atom_data['res_name']

            chain_info[chain_id]['residues'].add((res_num, res_name))
            chain_info[chain_id]['atoms'].append(atom_num)

    # Convert to summary
    summary = {}
    for chain_id, info in chain_info.items():
        residues_sorted = sorted(info['residues'], key=lambda x: x[0])
        atoms_sorted = sorted(info['atoms'])

        summary[chain_id] = {
            'num_atoms': len(atoms_sorted),
            'num_residues': len(residues_sorted),
            'atom_range': (min(atoms_sorted), max(atoms_sorted)),
            'residue_range': (residues_sorted[0][0], residues_sorted[-1][0]),
            'first_residues': residues_sorted[:5],
            'last_residues': residues_sorted[-5:]
        }

    return summary


def main():
    task = "3KXF_1"
    task_dir = Path(f"/home/xumy/work/development/AfterMD/input/pbc_1000frames_2step/{task}")
    tpr_file = task_dir / "md.tpr"
    index_file = task_dir / "phla_tcr.ndx"
    temp_pdb = Path(f"/tmp/{task}_protein.pdb")

    print("=" * 80)
    print(f"Chain and Residue Mapping for Task: {task}")
    print("=" * 80)

    # Extract PDB
    print("\nExtracting structure from trajectory...")
    if not extract_pdb_from_tpr(task_dir, temp_pdb):
        print("Failed to extract PDB")
        return

    # Read PDB atoms
    print("Reading PDB atom information...")
    atoms = read_pdb_atoms(temp_pdb)
    print(f"Total atoms in PDB: {len(atoms)}")

    # Analyze each group
    for group_name in ['pHLA', 'TCR']:
        print(f"\n{'=' * 80}")
        print(f"Index Group: {group_name}")
        print(f"{'=' * 80}")

        # Read index group
        atom_numbers = read_index_group(index_file, group_name)
        if atom_numbers is None:
            print(f"Group {group_name} not found in index file")
            continue

        print(f"Total atoms in group: {len(atom_numbers)}")

        # Analyze mapping
        summary = analyze_group_mapping(atoms, atom_numbers)

        # Display results
        for chain_id in sorted(summary.keys()):
            info = summary[chain_id]
            print(f"\n  Chain {chain_id}:")
            print(f"    Atoms: {info['num_atoms']} atoms")
            print(f"    Atom range: {info['atom_range'][0]} - {info['atom_range'][1]}")
            print(f"    Residues: {info['num_residues']} residues")
            print(f"    Residue range: {info['residue_range'][0]} - {info['residue_range'][1]}")
            first_res = [f"{r[0]}:{r[1]}" for r in info['first_residues']]
            last_res = [f"{r[0]}:{r[1]}" for r in info['last_residues']]
            print(f"    First 5 residues: {first_res}")
            print(f"    Last 5 residues: {last_res}")

    # Cleanup
    temp_pdb.unlink()
    print("\n" + "=" * 80)


if __name__ == "__main__":
    main()
