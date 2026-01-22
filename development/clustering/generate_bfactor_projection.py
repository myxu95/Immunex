#!/usr/bin/env python3
"""
Map RMSF values to B-factor column of PDB structure for visualization.
B-factor = RMSF^2 * 100 (converting nm to Angstrom^2)
"""

import pandas as pd
import numpy as np
from pathlib import Path

print("=" * 80)
print("B-Factor Projection: Mapping RMSF to PDB Structure")
print("=" * 80)

# File paths for the representative sample (3gsn_1)
task_name = "3gsn_1"
rmsf_file = f"/home/xumy/work/development/AfterMD/output/rmsf_cdr/{task_name}/rmsf_whole_protein.xvg"
pdb_input = f"/home/xumy/work/development/AfterMD/input/pbc_1000frames_2step/{task_name}/md_converted.pdb"
pdb_output = f"/home/xumy/work/development/AfterMD/output/rmsf_cdr/{task_name}/{task_name}_bfactor.pdb"

print(f"\nTask: {task_name}")
print(f"PDB ID: {task_name[:4].upper()}")
print(f"\nInput files:")
print(f"  RMSF data: {rmsf_file}")
print(f"  PDB structure: {pdb_input}")
print(f"  Output PDB: {pdb_output}")

# Read RMSF data from XVG file
print("\n" + "-" * 80)
print("Step 1: Reading RMSF data from XVG file...")
print("-" * 80)

rmsf_data = []
with open(rmsf_file, 'r') as f:
    for line in f:
        if line.startswith('#') or line.startswith('@'):
            continue
        parts = line.split()
        if len(parts) >= 2:
            residue_id = int(parts[0])
            rmsf_nm = float(parts[1])
            rmsf_data.append((residue_id, rmsf_nm))

rmsf_df = pd.DataFrame(rmsf_data, columns=['residue_id', 'rmsf_nm'])

print(f"Loaded RMSF data for {len(rmsf_df)} residues")
print(f"RMSF range: {rmsf_df['rmsf_nm'].min():.4f} - {rmsf_df['rmsf_nm'].max():.4f} nm")
print(f"Mean RMSF: {rmsf_df['rmsf_nm'].mean():.4f} nm")

# Convert RMSF (nm) to B-factor (Angstrom^2)
# B-factor = RMSF^2 * 100  (since 1 nm = 10 Angstrom)
# This represents the mean square displacement
rmsf_df['bfactor'] = (rmsf_df['rmsf_nm'] * 10) ** 2  # Convert nm to Angstrom, then square

print(f"\nB-factor range: {rmsf_df['bfactor'].min():.2f} - {rmsf_df['bfactor'].max():.2f} Å²")
print(f"Mean B-factor: {rmsf_df['bfactor'].mean():.2f} Å²")

# Create residue_id to bfactor mapping
rmsf_dict = dict(zip(rmsf_df['residue_id'], rmsf_df['bfactor']))

# Read PDB and update B-factor column
print("\n" + "-" * 80)
print("Step 2: Mapping RMSF to PDB B-factor column...")
print("-" * 80)

output_lines = []
atoms_updated = 0
atoms_no_rmsf = 0
current_residue = None

with open(pdb_input, 'r') as f:
    for line in f:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            # Parse PDB ATOM/HETATM line
            # Columns: https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM
            atom_serial = line[6:11].strip()
            atom_name = line[12:16].strip()
            alt_loc = line[16:17]
            res_name = line[17:20].strip()
            chain_id = line[21:22]
            res_seq = line[22:26].strip()
            icode = line[26:27]
            x = line[30:38].strip()
            y = line[38:46].strip()
            z = line[46:54].strip()
            occupancy = line[54:60].strip()
            old_bfactor = line[60:66].strip()
            element = line[76:78].strip() if len(line) > 76 else ''
            charge = line[78:80].strip() if len(line) > 78 else ''

            # Get residue ID
            try:
                residue_id = int(res_seq)
            except:
                # Keep original line if residue ID parsing fails
                output_lines.append(line)
                continue

            # Map to RMSF B-factor
            if residue_id in rmsf_dict:
                new_bfactor = rmsf_dict[residue_id]
                atoms_updated += 1
            else:
                # No RMSF data for this residue, use 0.0
                new_bfactor = 0.0
                atoms_no_rmsf += 1

            # Reconstruct PDB line with new B-factor
            new_line = (
                f"{line[0:6]}"          # Record name (ATOM/HETATM)
                f"{atom_serial:>5} "    # Atom serial
                f"{atom_name:>4}"       # Atom name
                f"{alt_loc}"            # Alternate location
                f"{res_name:>3} "       # Residue name
                f"{chain_id}"           # Chain ID
                f"{res_seq:>4}"         # Residue sequence number
                f"{icode}   "           # Insertion code + 3 spaces
                f"{x:>8}"               # X coordinate
                f"{y:>8}"               # Y coordinate
                f"{z:>8}"               # Z coordinate
                f"{occupancy:>6}"       # Occupancy
                f"{new_bfactor:>6.2f}"  # NEW B-factor (RMSF-derived)
                f"          "           # 10 spaces
                f"{element:>2}"         # Element symbol
                f"{charge:>2}"          # Charge
                "\n"
            )

            output_lines.append(new_line)

        else:
            # Non-ATOM lines (headers, etc.), keep as is
            output_lines.append(line)

print(f"Atoms with RMSF data: {atoms_updated}")
print(f"Atoms without RMSF data: {atoms_no_rmsf} (set to 0.0)")

# Write output PDB
print("\n" + "-" * 80)
print("Step 3: Writing output PDB file...")
print("-" * 80)

Path(pdb_output).parent.mkdir(parents=True, exist_ok=True)

with open(pdb_output, 'w') as f:
    # Add header comment
    f.write(f"REMARK   B-factor column contains RMSF-derived values (Angstrom^2)\n")
    f.write(f"REMARK   Source: {rmsf_file}\n")
    f.write(f"REMARK   Conversion: B-factor = (RMSF_nm * 10)^2\n")
    f.write(f"REMARK   Representative sample from 240 MD tasks\n")
    f.write(f"REMARK   Task: {task_name}, PDB: {task_name[:4].upper()}\n")

    # Write all lines
    f.writelines(output_lines)

output_size = Path(pdb_output).stat().st_size / 1024 / 1024  # MB

print(f"Output PDB saved: {pdb_output}")
print(f"File size: {output_size:.2f} MB")

print("\n" + "=" * 80)
print("B-Factor Projection Complete!")
print("=" * 80)

print("\nVisualization Instructions:")
print("-" * 80)
print("1. PyMOL:")
print(f"   load {pdb_output}")
print("   spectrum b, rainbow, minimum=0, maximum=50")
print("   show cartoon")
print("   set cartoon_putty_scale_max, 5")
print("   set cartoon_putty_scale_min, 0.5")
print("   cartoon putty")
print("")
print("2. VMD:")
print(f"   mol new {pdb_output}")
print("   mol modcolor 0 0 Beta")
print("   mol modstyle 0 0 NewCartoon")
print("")
print("3. ChimeraX:")
print(f"   open {pdb_output}")
print("   color bfactor palette rainbow")
print("   cartoon")

print("\n" + "=" * 80)
print("Color Interpretation:")
print("=" * 80)
print("  Blue/Cyan:    Low flexibility (rigid regions)")
print("  Green/Yellow: Moderate flexibility")
print("  Red:          High flexibility (CDR loops expected here)")
print("=" * 80)
