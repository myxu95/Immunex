#!/usr/bin/env python3
"""
Comprehensive Index Generation Test

Generate complete index with:
1. Standard components (pHLA, TCR, peptide, etc.)
2. CDR loops (CDR1/2/3 for alpha and beta)
3. Backbone variants for each component

Author: Immunex Development Team
Date: 2026-03-18
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

from immunex.core.index_generation import (
    IndexGenerator,
    IndexGenerationInput,
    IndexGenerationMethod,
    ComponentDefinition
)
import MDAnalysis as mda


def detect_cdr_loops_anarci(pdb_file: str):
    """
    Detect CDR loops using ANARCI.
    
    Returns:
        dict: {
            'CDR1_alpha': (start, end),
            'CDR2_alpha': (start, end),
            'CDR3_alpha': (start, end),
            'CDR1_beta': (start, end),
            'CDR2_beta': (start, end),
            'CDR3_beta': (start, end)
        }
    """
    try:
        from immunex.utils.cdr_manager import CDRManager
        
        cdr_manager = CDRManager()
        
        # Load structure
        u = mda.Universe(pdb_file)
        
        # Get TCR alpha chain (D)
        tcr_alpha = u.select_atoms("chainID D and protein and name CA")
        alpha_sequence = ''.join([
            mda.lib.util.convert_aa_code(r.resname)
            for r in tcr_alpha.residues
        ])
        
        # Get TCR beta chain (E)
        tcr_beta = u.select_atoms("chainID E and protein and name CA")
        beta_sequence = ''.join([
            mda.lib.util.convert_aa_code(r.resname)
            for r in tcr_beta.residues
        ])
        
        print(f"TCR-alpha sequence length: {len(alpha_sequence)}")
        print(f"TCR-beta sequence length: {len(beta_sequence)}")
        
        # Detect CDR loops
        cdr_regions = {}
        
        # Alpha chain CDRs
        alpha_result = cdr_manager.detect_cdr_regions(alpha_sequence, 'alpha')
        if alpha_result['success']:
            for cdr_name, cdr_info in alpha_result['cdr_regions'].items():
                if cdr_info['detected']:
                    # Convert sequence position to residue IDs
                    start_resid = tcr_alpha.residues[cdr_info['start']].resid
                    end_resid = tcr_alpha.residues[cdr_info['end']-1].resid
                    cdr_regions[f"{cdr_name}_alpha"] = (start_resid, end_resid)
                    print(f"  {cdr_name}_alpha: resid {start_resid}-{end_resid}")
        
        # Beta chain CDRs
        beta_result = cdr_manager.detect_cdr_regions(beta_sequence, 'beta')
        if beta_result['success']:
            for cdr_name, cdr_info in beta_result['cdr_regions'].items():
                if cdr_info['detected']:
                    start_resid = tcr_beta.residues[cdr_info['start']].resid
                    end_resid = tcr_beta.residues[cdr_info['end']-1].resid
                    cdr_regions[f"{cdr_name}_beta"] = (start_resid, end_resid)
                    print(f"  {cdr_name}_beta: resid {start_resid}-{end_resid}")
        
        return cdr_regions
        
    except Exception as e:
        print(f"CDR detection failed: {e}")
        return {}


def generate_backbone_index(input_file: str, output_file: str, component_name: str, 
                           selection: str):
    """
    Generate backbone index for a specific component.
    
    Args:
        input_file: Input PDB file
        output_file: Output index file
        component_name: Name of component
        selection: MDAnalysis selection for the component
    """
    u = mda.Universe(input_file)
    
    # Select backbone atoms
    backbone_selection = f"({selection}) and (name CA C N O)"
    backbone_atoms = u.select_atoms(backbone_selection)
    
    if len(backbone_atoms) == 0:
        print(f"Warning: No backbone atoms found for {component_name}")
        return
    
    # Write to index file
    with open(output_file, 'a') as f:
        f.write(f"[ {component_name}_backbone ]\n")
        
        atom_indices = backbone_atoms.indices + 1  # 1-indexed
        for i in range(0, len(atom_indices), 15):
            line = ' '.join(f"{idx:4d}" for idx in atom_indices[i:i+15])
            f.write(f"{line}\n")
        
        f.write("\n")
    
    print(f"  Added {component_name}_backbone: {len(backbone_atoms)} atoms")


def main():
    print("=" * 70)
    print("Comprehensive Index Generation")
    print("=" * 70)
    print()
    
    # Input PDB file
    pdb_file = "input/standardizedpdbs/raw_data/1AO7_SOL.pdb"
    
    if not Path(pdb_file).exists():
        print(f"File not found: {pdb_file}")
        return 1
    
    output_dir = Path("development/test_output")
    output_dir.mkdir(exist_ok=True)
    output_file = str(output_dir / "1ao7_comprehensive.ndx")
    
    print(f"Input:  {pdb_file}")
    print(f"Output: {output_file}")
    print()
    
    # Step 1: Generate standard components
    print("Step 1: Generating standard components...")
    print("-" * 70)
    
    generator = IndexGenerator()
    
    input_params = IndexGenerationInput(
        topology=pdb_file,
        method=IndexGenerationMethod.PDB_BASED,
        components=[
            ComponentDefinition(name="pHLA", chains=['A', 'B', 'C']),
            ComponentDefinition(name="TCR", chains=['D', 'E']),
            ComponentDefinition(name="peptide", chains=['C']),
            ComponentDefinition(name="HLA_alpha", chains=['A']),
            ComponentDefinition(name="beta2m", chains=['B']),
            ComponentDefinition(name="TCR_alpha", chains=['D']),
            ComponentDefinition(name="TCR_beta", chains=['E']),
        ],
        output_file=output_file,
        auto_standardize=False
    )
    
    result = generator.generate(input_params)
    
    if not result.success:
        print(f"Failed: {result.error_message}")
        return 1
    
    print(f"Generated {len(result.components)} standard components")
    print()
    
    # Step 2: Detect and add CDR loops
    print("Step 2: Detecting CDR loops...")
    print("-" * 70)
    
    cdr_regions = detect_cdr_loops_anarci(pdb_file)
    
    if cdr_regions:
        u = mda.Universe(pdb_file)
        
        with open(output_file, 'a') as f:
            for cdr_name, (start_resid, end_resid) in cdr_regions.items():
                # Determine chain
                if 'alpha' in cdr_name:
                    chain_id = 'D'
                else:
                    chain_id = 'E'
                
                # Select CDR atoms
                cdr_atoms = u.select_atoms(
                    f"chainID {chain_id} and protein and resid {start_resid}:{end_resid}"
                )
                
                # Write to index
                f.write(f"[ {cdr_name} ]\n")
                atom_indices = cdr_atoms.indices + 1
                for i in range(0, len(atom_indices), 15):
                    line = ' '.join(f"{idx:4d}" for idx in atom_indices[i:i+15])
                    f.write(f"{line}\n")
                f.write("\n")
                
                print(f"  Added {cdr_name}: resid {start_resid}-{end_resid}, {len(cdr_atoms)} atoms")
        
        print(f"\nAdded {len(cdr_regions)} CDR regions")
    else:
        print("No CDR regions detected (ANARCI might not be available)")
    
    print()
    
    # Step 3: Add backbone variants
    print("Step 3: Adding backbone variants...")
    print("-" * 70)
    
    backbone_components = [
        ("pHLA", "chainID A B C and protein"),
        ("TCR", "chainID D E and protein"),
        ("peptide", "chainID C and protein"),
        ("HLA_alpha", "chainID A and protein"),
        ("beta2m", "chainID B and protein"),
        ("TCR_alpha", "chainID D and protein"),
        ("TCR_beta", "chainID E and protein"),
    ]
    
    for comp_name, selection in backbone_components:
        generate_backbone_index(pdb_file, output_file, comp_name, selection)
    
    # Also add CDR backbones
    if cdr_regions:
        for cdr_name, (start_resid, end_resid) in cdr_regions.items():
            if 'alpha' in cdr_name:
                chain_id = 'D'
            else:
                chain_id = 'E'
            
            selection = f"chainID {chain_id} and protein and resid {start_resid}:{end_resid}"
            generate_backbone_index(pdb_file, output_file, cdr_name, selection)
    
    print()
    
    # Summary
    print("=" * 70)
    print("Summary")
    print("=" * 70)
    
    # Count total groups
    with open(output_file, 'r') as f:
        content = f.read()
        groups = content.count('[')
    
    print(f"\nTotal groups in index: {groups}")
    print(f"Output file: {output_file}")
    
    # Show group list
    print("\nGenerated groups:")
    with open(output_file, 'r') as f:
        for line in f:
            if line.strip().startswith('['):
                group_name = line.strip()[1:-1].strip()
                print(f"  - {group_name}")
    
    print("\nUsage examples:")
    print("-" * 70)
    print("""
# RMSD of TCR backbone aligned to pHLA backbone
gmx rms -f md.xtc -s md.tpr -n 1ao7_comprehensive.ndx \\
        -a TCR_backbone -r pHLA_backbone

# RMSF of CDR3_alpha backbone
gmx rmsf -f md.xtc -s md.tpr -n 1ao7_comprehensive.ndx \\
         -res -o cdr3a_rmsf.xvg

# Extract peptide backbone
gmx trjconv -f md.xtc -s md.tpr -n 1ao7_comprehensive.ndx \\
            -select peptide_backbone -o peptide_bb.pdb
    """)
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
