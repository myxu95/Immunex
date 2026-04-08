#!/usr/bin/env python3
"""
Multi-Model PDB Concatenation - Usage Examples

This example demonstrates how to use the MultiModelConcatenator to merge
multi-model PDB files (common in biological assemblies) into single PDB files
with continuous atom numbering.

Author: Immunex Development Team
Date: 2026-01
"""

from immunex.utils import MultiModelConcatenator

# ============================================================================
# Example 1: Basic Usage - Single File
# ============================================================================
print("=" * 70)
print("Example 1: Concatenate a Single Multi-Model PDB File")
print("=" * 70)

concatenator = MultiModelConcatenator()

# Check how many models are in the file
n_models = concatenator.count_models("input/standardizedpdbs/pdbs_raw/6VMA.pdb")
print(f"Found {n_models} models in 6VMA.pdb\n")

# Merge the models
result = concatenator.concatenate_models(
    input_pdb="input/standardizedpdbs/pdbs_raw/6VMA.pdb",
    output_pdb="output/6VMA_concatenated.pdb"
)

print(f"Results:")
print(f"  Input:  {result['input_file']}")
print(f"  Output: {result['output_file']}")
print(f"  Number of models merged: {result['n_models']}")
print(f"  Total atoms in output: {result['total_atoms']}")
print(f"  Atoms per model: {result['atoms_per_model']}")
print(f"  Status: {result['status']}")
print()


# ============================================================================
# Example 2: Batch Processing
# ============================================================================
print("=" * 70)
print("Example 2: Batch Process Multiple PDB Files")
print("=" * 70)

results = concatenator.batch_concatenate(
    input_dir="input/standardizedpdbs/pdbs_raw",
    output_dir="output/concatenated_pdbs",
    pattern="6VM*.pdb"
)

print(f"\nProcessed {len(results)} files:")
for result in results:
    if result.get('status') == 'success':
        print(f"  ✓ {result['input_file']}")
        print(f"    -> {result['n_models']} models, {result['total_atoms']} atoms")
    else:
        print(f"  ✗ {result['input_file']}: {result.get('error', 'Unknown error')}")

print()


# ============================================================================
# Example 3: Understanding the Output
# ============================================================================
print("=" * 70)
print("Example 3: What Changed in the Output?")
print("=" * 70)

print("""
Input PDB (multi-model):
  MODEL        1
  ATOM      1  N   GLY A   1      12.603   4.179  -5.862  ...
  ATOM      2  CA  GLY A   1      12.104   4.967  -6.974  ...
  ...
  ATOM   3177  ...
  ENDMDL
  MODEL        2
  ATOM      1  N   GLY A   1      12.100   4.100  -5.800  ...
  ATOM      2  CA  GLY A   1      11.900   4.900  -6.900  ...
  ...
  ATOM   3493  ...
  ENDMDL
  END

Output PDB (concatenated):
  ATOM      1  N   GLY A   1      12.603   4.179  -5.862  ...
  ATOM      2  CA  GLY A   1      12.104   4.967  -6.974  ...
  ...
  ATOM   3177  ...
  ATOM   3178  N   GLY A   1      12.100   4.100  -5.800  ...  <- MODEL 2 starts here
  ATOM   3179  CA  GLY A   1      11.900   4.900  -6.900  ...
  ...
  ATOM   6670  ...
  END

Key Changes:
  ✓ MODEL and ENDMDL records removed
  ✓ Atom numbers are continuous (1, 2, 3, ..., 6670)
  ✓ All coordinates preserved unchanged
  ✓ Residue numbering preserved (may have duplicates)
  ✓ HEADER/CRYST1 metadata preserved
""")


# ============================================================================
# Example 4: Using Command Line
# ============================================================================
print("=" * 70)
print("Example 4: Command-Line Usage")
print("=" * 70)

print("""
# Single file processing:
python scripts/concatenate_multimodel_pdbs.py \\
    input/standardizedpdbs/pdbs_raw/6VMA.pdb \\
    -o output/6VMA_concatenated.pdb

# Batch processing:
python scripts/concatenate_multimodel_pdbs.py \\
    input/standardizedpdbs/pdbs_raw/6VM*.pdb \\
    -o output/concatenated_pdbs/ \\
    --batch

# With JSON report:
python scripts/concatenate_multimodel_pdbs.py \\
    input/standardizedpdbs/pdbs_raw/6VMA.pdb \\
    -o output/6VMA_concatenated.pdb \\
    --report output/report.json

# Without metadata (HEADER/CRYST1):
python scripts/concatenate_multimodel_pdbs.py \\
    input/standardizedpdbs/pdbs_raw/6VMA.pdb \\
    -o output/6VMA_concatenated.pdb \\
    --no-remarks
""")

print("=" * 70)
print("All examples completed!")
print("=" * 70)
