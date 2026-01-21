#!/usr/bin/env python3
"""
Create a minimal test PDB file for testing intelligent chain identification
"""

def create_test_pdb():
    """Create a minimal 5-chain pHLA-TCR PDB structure."""

    pdb_content = """HEADER    IMMUNE SYSTEM                           01-JAN-26   TEST
TITLE     MINIMAL PHLA-TCR COMPLEX FOR TESTING
REMARK   1 CREATED FOR TESTING INTELLIGENT CHAIN IDENTIFICATION
"""

    # Chain V: HLA-alpha (longest, 50 residues)
    atom_num = 1
    for i in range(1, 51):
        x, y, z = 10.0 + i*0.1, 10.0, 10.0
        pdb_content += f"ATOM  {atom_num:5d}  CA  ALA V{i:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C\n"
        atom_num += 1

    # Chain Y: Beta2-microglobulin (30 residues)
    for i in range(1, 31):
        x, y, z = 20.0, 20.0 + i*0.1, 20.0
        pdb_content += f"ATOM  {atom_num:5d}  CA  GLY Y{i:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C\n"
        atom_num += 1

    # Chain X: Peptide (shortest, 10 residues)
    for i in range(1, 11):
        x, y, z = 30.0, 30.0, 30.0 + i*0.1
        pdb_content += f"ATOM  {atom_num:5d}  CA  VAL X{i:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C\n"
        atom_num += 1

    # Chain Z: TCR-alpha (40 residues)
    for i in range(1, 41):
        x, y, z = 40.0 + i*0.1, 40.0, 40.0
        pdb_content += f"ATOM  {atom_num:5d}  CA  LEU Z{i:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C\n"
        atom_num += 1

    # Chain W: TCR-beta (45 residues)
    for i in range(1, 46):
        x, y, z = 50.0, 50.0 + i*0.1, 50.0
        pdb_content += f"ATOM  {atom_num:5d}  CA  ILE W{i:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C\n"
        atom_num += 1

    pdb_content += "END\n"

    return pdb_content


if __name__ == "__main__":
    from pathlib import Path

    # Create test directory
    test_dir = Path("test_data")
    test_dir.mkdir(exist_ok=True)

    # Write test PDB
    test_pdb = test_dir / "test_complex.pdb"
    with open(test_pdb, 'w') as f:
        f.write(create_test_pdb())

    print(f"✓ Created test PDB: {test_pdb}")
    print(f"  Chains: V(50aa), Y(30aa), X(10aa), Z(40aa), W(45aa)")
    print(f"  Expected mapping:")
    print(f"    X → C (peptide, shortest)")
    print(f"    Y → B (beta2m, ~30aa)")
    print(f"    Z → D or E (TCR, 40aa)")
    print(f"    W → E or D (TCR, 45aa)")
    print(f"    V → A (HLA-alpha, longest)")
