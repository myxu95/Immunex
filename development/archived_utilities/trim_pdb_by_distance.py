#!/usr/bin/env python3
"""
Trim PDB file by removing residues beyond distance threshold from peptide center.

Remove all residues that are >8nm away from the central residue of the peptide chain.
"""
import sys
import argparse
from pathlib import Path
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array
import numpy as np

def find_peptide_chain(universe):
    """Find peptide chain (shortest chain, typically <=20 residues)."""
    chains = {}
    for seg in universe.segments:
        # Count protein residues only
        protein_residues = universe.select_atoms(
            f"segid {seg.segid} and protein"
        ).residues
        if len(protein_residues) > 0:
            chains[seg.segid] = len(protein_residues)

    # Find shortest chain
    if not chains:
        raise ValueError("No protein chains found!")

    peptide_segid = min(chains.items(), key=lambda x: x[1])[0]
    peptide_length = chains[peptide_segid]

    print(f"Peptide chain identified: {peptide_segid} ({peptide_length} residues)")
    return peptide_segid

def get_peptide_center_residue(universe, peptide_segid):
    """Get the central residue of peptide chain."""
    peptide_residues = universe.select_atoms(
        f"segid {peptide_segid} and protein"
    ).residues

    n_residues = len(peptide_residues)
    center_idx = n_residues // 2
    center_residue = peptide_residues[center_idx]

    print(f"Peptide center residue: {center_residue.resname}{center_residue.resid} (position {center_idx+1}/{n_residues})")
    return center_residue

def trim_pdb_by_distance(input_pdb, output_pdb, distance_cutoff_nm=8.0, verbose=True):
    """
    Remove residues beyond distance threshold from peptide center.

    Args:
        input_pdb: Input PDB file path
        output_pdb: Output PDB file path
        distance_cutoff_nm: Distance cutoff in nanometers (default: 8.0 nm)
        verbose: Print detailed information

    Returns:
        Dictionary with statistics
    """
    distance_cutoff_angstrom = distance_cutoff_nm * 10.0  # Convert nm to Angstrom

    print("=" * 80)
    print("PDB Trimming by Distance from Peptide Center")
    print("=" * 80)
    print(f"Input:  {input_pdb}")
    print(f"Output: {output_pdb}")
    print(f"Distance cutoff: {distance_cutoff_nm} nm ({distance_cutoff_angstrom} Å)")
    print("=" * 80)

    # Load structure
    u = mda.Universe(input_pdb)

    # Find peptide chain
    peptide_segid = find_peptide_chain(u)

    # Get peptide center residue
    center_residue = get_peptide_center_residue(u, peptide_segid)
    center_com = center_residue.atoms.center_of_mass()

    print(f"Center of mass: ({center_com[0]:.2f}, {center_com[1]:.2f}, {center_com[2]:.2f}) Å")
    print("=" * 80)

    # Calculate distances for all residues
    all_residues = u.select_atoms("protein").residues

    residues_to_keep = []
    residues_to_remove = []

    for residue in all_residues:
        residue_com = residue.atoms.center_of_mass()
        distance = np.linalg.norm(residue_com - center_com)

        residue_info = {
            'segid': residue.segment.segid,
            'resid': residue.resid,
            'resname': residue.resname,
            'distance': distance
        }

        if distance <= distance_cutoff_angstrom:
            residues_to_keep.append(residue)
        else:
            residues_to_remove.append(residue_info)

    # Print removed residues
    print(f"\nResidues to REMOVE (distance > {distance_cutoff_nm} nm):")
    print("-" * 80)
    if residues_to_remove:
        # Group by chain
        by_chain = {}
        for res in residues_to_remove:
            segid = res['segid']
            if segid not in by_chain:
                by_chain[segid] = []
            by_chain[segid].append(res)

        for segid in sorted(by_chain.keys()):
            print(f"\nChain {segid}:")
            for res in by_chain[segid]:
                print(f"  {res['resname']:>3s} {res['resid']:>4d}  (distance: {res['distance']/10:.2f} nm)")
    else:
        print("  None")

    print("\n" + "=" * 80)
    print("Statistics:")
    print("=" * 80)
    print(f"Total residues:         {len(all_residues)}")
    print(f"Residues to keep:       {len(residues_to_keep)}")
    print(f"Residues to remove:     {len(residues_to_remove)}")
    print(f"Retention rate:         {len(residues_to_keep)/len(all_residues)*100:.1f}%")

    # Create new selection with atoms from kept residues
    # Use atom indices for efficiency (avoids long selection strings)
    atom_indices = []
    for res in residues_to_keep:
        atom_indices.extend(res.atoms.indices)

    if not atom_indices:
        raise ValueError("No residues left after trimming!")

    # Select atoms by indices
    trimmed_atoms = u.atoms[atom_indices]

    print(f"Total atoms to keep:    {len(trimmed_atoms)}")
    print("=" * 80)

    # Write output PDB
    output_path = Path(output_pdb)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    trimmed_atoms.write(str(output_pdb))
    print(f"\nTrimmed PDB written to: {output_pdb}")

    # Return statistics
    return {
        'input_file': str(input_pdb),
        'output_file': str(output_pdb),
        'cutoff_nm': distance_cutoff_nm,
        'total_residues': len(all_residues),
        'kept_residues': len(residues_to_keep),
        'removed_residues': len(residues_to_remove),
        'removed_details': residues_to_remove
    }

def main():
    parser = argparse.ArgumentParser(
        description="Trim PDB by removing residues beyond distance from peptide center"
    )
    parser.add_argument(
        "input_pdb",
        help="Input PDB file"
    )
    parser.add_argument(
        "-o", "--output",
        help="Output PDB file (default: input_trimmed.pdb)",
        default=None
    )
    parser.add_argument(
        "-d", "--distance",
        type=float,
        default=8.0,
        help="Distance cutoff in nanometers (default: 8.0 nm)"
    )

    args = parser.parse_args()

    # Determine output filename
    if args.output is None:
        input_path = Path(args.input_pdb)
        output_pdb = input_path.parent / f"{input_path.stem}_trimmed.pdb"
    else:
        output_pdb = args.output

    # Run trimming
    try:
        stats = trim_pdb_by_distance(
            args.input_pdb,
            output_pdb,
            distance_cutoff_nm=args.distance
        )
        print("\nCompleted successfully!")
        return 0
    except Exception as e:
        print(f"\nERROR: {e}", file=sys.stderr)
        return 1

if __name__ == "__main__":
    sys.exit(main())
