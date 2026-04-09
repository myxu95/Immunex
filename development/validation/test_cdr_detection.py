#!/usr/bin/env python3
"""
Test CDR detection and add to index file

This script:
1. Detects CDR loops using ANARCI
2. Adds CDR regions to existing index file
3. Adds CDR backbone variants

Author: Immunex Development Team
Date: 2026-03-18
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import MDAnalysis as mda


def detect_cdr_with_anarci(sequence: str, chain_type: str):
    """
    Detect CDR regions using ANARCI directly.
    
    Args:
        sequence: Amino acid sequence
        chain_type: 'alpha' or 'beta'
    
    Returns:
        dict: CDR regions with positions
    """
    try:
        from anarci import anarci
        
        # Run ANARCI
        results = anarci([('seq', sequence)], scheme='imgt', output=False)

        if not results or not results[0]:
            print(f"  ANARCI failed for {chain_type} chain")
            return {}

        # Unpack ANARCI results (returns tuple of 3 elements)
        numbering_results, alignment_details, hit_table = results

        if not numbering_results or not numbering_results[0]:
            print(f"  No results for {chain_type} chain")
            return {}

        # Get first sequence result (numbering, start_idx, end_idx)
        numbering_info = numbering_results[0][0]
        numbering = numbering_info[0]
        
        if numbering is None:
            print(f"  No numbering for {chain_type} chain")
            return {}
        
        # IMGT CDR definitions (IMGT numbering)
        cdr_definitions = {
            'CDR1': (27, 38),    # IMGT 27-38
            'CDR2': (56, 65),    # IMGT 56-65
            'CDR3': (105, 117)   # IMGT 105-117
        }
        
        cdr_regions = {}
        
        # Build position mapping: imgt_number -> sequence_position
        imgt_to_seq = {}
        for seq_pos, (imgt_num, aa) in enumerate(numbering):
            # imgt_num is a tuple (number, insert_code)
            if imgt_num is not None and aa != '-':
                imgt_to_seq[imgt_num[0]] = seq_pos
        
        # Detect each CDR
        for cdr_name, (start_imgt, end_imgt) in cdr_definitions.items():
            # Find sequence positions
            seq_positions = []
            for imgt_num in range(start_imgt, end_imgt + 1):
                if imgt_num in imgt_to_seq:
                    seq_positions.append(imgt_to_seq[imgt_num])
            
            if seq_positions:
                start_pos = min(seq_positions)
                end_pos = max(seq_positions)
                cdr_regions[cdr_name] = {
                    'start': start_pos,
                    'end': end_pos + 1,  # Exclusive end
                    'length': end_pos - start_pos + 1
                }
                
                # Extract sequence
                cdr_seq = sequence[start_pos:end_pos+1]
                print(f"  {cdr_name}: pos {start_pos}-{end_pos}, seq={cdr_seq}")
        
        return cdr_regions
        
    except ImportError:
        print("  ANARCI not available")
        return {}
    except Exception as e:
        print(f"  ANARCI error: {e}")
        return {}


def add_cdr_to_index(pdb_file: str, index_file: str):
    """
    Detect CDR loops and add to index file.
    """
    print("Detecting CDR loops with ANARCI...")
    print("-" * 70)
    
    u = mda.Universe(pdb_file)
    
    # Get TCR alpha sequence (chain D)
    tcr_alpha = u.select_atoms("chainID D and protein and name CA")
    alpha_sequence = ''.join([
        mda.lib.util.convert_aa_code(r.resname)
        for r in tcr_alpha.residues
    ])
    
    print(f"TCR-alpha (chain D): {len(alpha_sequence)} residues")
    alpha_cdrs = detect_cdr_with_anarci(alpha_sequence, 'alpha')
    
    # Get TCR beta sequence (chain E)
    tcr_beta = u.select_atoms("chainID E and protein and name CA")
    beta_sequence = ''.join([
        mda.lib.util.convert_aa_code(r.resname)
        for r in tcr_beta.residues
    ])
    
    print(f"\nTCR-beta (chain E): {len(beta_sequence)} residues")
    beta_cdrs = detect_cdr_with_anarci(beta_sequence, 'beta')
    
    print()
    
    # Add to index file
    cdr_count = 0
    
    with open(index_file, 'a') as f:
        # Add alpha CDRs
        for cdr_name, cdr_info in alpha_cdrs.items():
            start_idx = cdr_info['start']
            end_idx = cdr_info['end'] - 1

            # Boundary check
            if end_idx >= len(tcr_alpha.residues):
                print(f"Warning: {cdr_name}_alpha end_idx {end_idx} exceeds sequence length {len(tcr_alpha.residues)}, truncating")
                end_idx = len(tcr_alpha.residues) - 1

            # Get residue IDs
            start_resid = tcr_alpha.residues[start_idx].resid
            end_resid = tcr_alpha.residues[end_idx].resid
            
            # Select all atoms in CDR
            cdr_atoms = u.select_atoms(
                f"chainID D and protein and resid {start_resid}:{end_resid}"
            )
            
            # Write group
            f.write(f"[ {cdr_name}_alpha ]\n")
            atom_indices = cdr_atoms.indices + 1
            for i in range(0, len(atom_indices), 15):
                line = ' '.join(f"{idx:4d}" for idx in atom_indices[i:i+15])
                f.write(f"{line}\n")
            f.write("\n")
            
            print(f"Added {cdr_name}_alpha: resid {start_resid}-{end_resid}, {len(cdr_atoms)} atoms")
            cdr_count += 1
            
            # Also add backbone variant
            cdr_backbone = u.select_atoms(
                f"chainID D and protein and resid {start_resid}:{end_resid} and (name CA C N O)"
            )
            
            f.write(f"[ {cdr_name}_alpha_backbone ]\n")
            atom_indices = cdr_backbone.indices + 1
            for i in range(0, len(atom_indices), 15):
                line = ' '.join(f"{idx:4d}" for idx in atom_indices[i:i+15])
                f.write(f"{line}\n")
            f.write("\n")
            
            print(f"Added {cdr_name}_alpha_backbone: {len(cdr_backbone)} atoms")
            cdr_count += 1
        
        # Add beta CDRs
        for cdr_name, cdr_info in beta_cdrs.items():
            start_idx = cdr_info['start']
            end_idx = cdr_info['end'] - 1

            # Boundary check
            if end_idx >= len(tcr_beta.residues):
                print(f"Warning: {cdr_name}_beta end_idx {end_idx} exceeds sequence length {len(tcr_beta.residues)}, truncating")
                end_idx = len(tcr_beta.residues) - 1

            start_resid = tcr_beta.residues[start_idx].resid
            end_resid = tcr_beta.residues[end_idx].resid
            
            cdr_atoms = u.select_atoms(
                f"chainID E and protein and resid {start_resid}:{end_resid}"
            )
            
            f.write(f"[ {cdr_name}_beta ]\n")
            atom_indices = cdr_atoms.indices + 1
            for i in range(0, len(atom_indices), 15):
                line = ' '.join(f"{idx:4d}" for idx in atom_indices[i:i+15])
                f.write(f"{line}\n")
            f.write("\n")
            
            print(f"Added {cdr_name}_beta: resid {start_resid}-{end_resid}, {len(cdr_atoms)} atoms")
            cdr_count += 1
            
            # Backbone variant
            cdr_backbone = u.select_atoms(
                f"chainID E and protein and resid {start_resid}:{end_resid} and (name CA C N O)"
            )
            
            f.write(f"[ {cdr_name}_beta_backbone ]\n")
            atom_indices = cdr_backbone.indices + 1
            for i in range(0, len(atom_indices), 15):
                line = ' '.join(f"{idx:4d}" for idx in atom_indices[i:i+15])
                f.write(f"{line}\n")
            f.write("\n")
            
            print(f"Added {cdr_name}_beta_backbone: {len(cdr_backbone)} atoms")
            cdr_count += 1
    
    print(f"\nTotal CDR groups added: {cdr_count}")
    
    return cdr_count


def main():
    pdb_file = "input/standardizedpdbs/raw_data/1AO7_SOL.pdb"
    index_file = "development/test_output/1ao7_comprehensive.ndx"
    
    if not Path(pdb_file).exists():
        print(f"PDB file not found: {pdb_file}")
        return 1
    
    if not Path(index_file).exists():
        print(f"Index file not found: {index_file}")
        print("Please run test_comprehensive_index.py first")
        return 1
    
    print("=" * 70)
    print("Adding CDR Loops to Index File")
    print("=" * 70)
    print()
    
    cdr_count = add_cdr_to_index(pdb_file, index_file)
    
    if cdr_count > 0:
        print("\n" + "=" * 70)
        print("Summary")
        print("=" * 70)
        
        # Count total groups
        with open(index_file, 'r') as f:
            content = f.read()
            total_groups = content.count('[')
        
        print(f"\nTotal groups in index: {total_groups}")
        print(f"Updated file: {index_file}")
        
        # Show CDR groups
        print("\nCDR groups added:")
        with open(index_file, 'r') as f:
            for line in f:
                if line.strip().startswith('[') and 'CDR' in line:
                    group_name = line.strip()[1:-1].strip()
                    print(f"  - {group_name}")
        
        print("\nUsage examples:")
        print("-" * 70)
        print("""
# Calculate RMSD of CDR3-beta backbone
gmx rms -f md.xtc -s md.tpr -n 1ao7_comprehensive.ndx \\
        -select CDR3_beta_backbone

# Calculate RMSF of CDR1-alpha
gmx rmsf -f md.xtc -s md.tpr -n 1ao7_comprehensive.ndx \\
         -select CDR1_alpha -res -o cdr1a_rmsf.xvg

# Extract CDR3-alpha structure
gmx trjconv -f md.xtc -s md.tpr -n 1ao7_comprehensive.ndx \\
            -select CDR3_alpha -o cdr3a.pdb
        """)
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
