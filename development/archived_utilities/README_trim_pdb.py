#!/usr/bin/env python3
"""
Quick test of PDB trimming scripts

Usage examples for single and batch PDB trimming.
"""

print("=" * 80)
print("PDB Trimming Scripts - Usage Examples")
print("=" * 80)

print("\n1. Single PDB file trimming:")
print("-" * 80)
print("python development/utilities/trim_pdb_by_distance.py \\")
print("    input/standardizedpdbs/pdbs_raw/1ao7.pdb \\")
print("    -o /tmp/1ao7_trimmed.pdb \\")
print("    -d 8.0")

print("\n2. Batch PDB trimming (process entire directory):")
print("-" * 80)
print("python development/utilities/batch_trim_pdb_by_distance.py \\")
print("    input/standardizedpdbs/pdbs_raw/ \\")
print("    -o /tmp/trimmed_pdbs/ \\")
print("    -d 8.0 \\")
print("    -w 4 \\")
print("    --save-report /tmp/trimming_report.json")

print("\n3. Process only specific PDB files:")
print("-" * 80)
print("python development/utilities/batch_trim_pdb_by_distance.py \\")
print("    input/standardizedpdbs/pdbs_raw/ \\")
print("    -o /tmp/trimmed_pdbs/ \\")
print("    -d 8.0 \\")
print("    -p '1*.pdb'  # Only PDBs starting with '1'")

print("\n" + "=" * 80)
print("Parameters:")
print("=" * 80)
print("  -d, --distance    Distance cutoff in nanometers (default: 8.0 nm)")
print("  -o, --output      Output file/directory")
print("  -w, --workers     Number of parallel workers (batch only)")
print("  -p, --pattern     File pattern for batch processing")
print("  --save-report     Save detailed JSON report (batch only)")

print("\n" + "=" * 80)
print("Script Features:")
print("=" * 80)
print("✓ Automatically identifies peptide chain (shortest chain)")
print("✓ Finds central residue of peptide")
print("✓ Calculates distances from peptide center to all other residues")
print("✓ Removes residues beyond distance threshold")
print("✓ Prints detailed list of removed residues (by chain)")
print("✓ Shows statistics (total/kept/removed residues)")
print("✓ Batch processing with parallel workers")

print("\n" + "=" * 80)
