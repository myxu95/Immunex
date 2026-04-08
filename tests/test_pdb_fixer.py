#!/usr/bin/env python3
"""
Test PDB Structure Fixer
"""
import sys
sys.path.insert(0, '/home/xumy/work/development/Immunex')

from immunex.analysis import PDBStructureFixer
from pathlib import Path

print("=" * 80)
print("Testing PDB Structure Fixer")
print("=" * 80)

# Test file
test_file = "input/standardizedpdbs/trimmed_sd_pdbs/1AO7_sd.pdb"

if not Path(test_file).exists():
    print(f"Test file not found: {test_file}")
    sys.exit(1)

print(f"\nTest file: {test_file}")
print("=" * 80)

# Initialize fixer
fixer = PDBStructureFixer(
    add_missing_residues=True,
    add_missing_atoms=True,
    add_hydrogens=True,
    remove_heterogens=False,
    pH=7.0
)

print("\nFixer settings:")
print(f"  Add missing residues: {fixer.add_missing_residues}")
print(f"  Add missing atoms: {fixer.add_missing_atoms}")
print(f"  Add hydrogens: {fixer.add_hydrogens}")
print(f"  Remove heterogens: {fixer.remove_heterogens}")
print(f"  pH: {fixer.pH}")

# Test fixing
output_file = "/tmp/test_fixed.pdb"

print("\n" + "=" * 80)
print("Running structure fix...")
print("=" * 80)

try:
    stats = fixer.fix_structure(test_file, output_file, verbose=True)

    print("\n" + "=" * 80)
    print("Fix Statistics:")
    print("=" * 80)
    print(f"Fixed successfully: {stats['fixed']}")
    print(f"Missing residues found: {stats['missing_residues']}")
    print(f"Missing heavy atoms found: {stats['missing_heavy_atoms']}")
    print(f"Hydrogens added: {stats.get('added_hydrogens', False)}")
    print(f"Output file: {stats['output_file']}")

    # Check output file size
    output_path = Path(output_file)
    if output_path.exists():
        size_kb = output_path.stat().st_size / 1024
        print(f"Output file size: {size_kb:.1f} KB")
        print("\n✓ Test PASSED")
    else:
        print("\n✗ Test FAILED: Output file not created")

except Exception as e:
    print(f"\n✗ Test FAILED: {e}")
    import traceback
    traceback.print_exc()

print("=" * 80)
