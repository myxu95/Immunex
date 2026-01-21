# CDRSelector Usage Guide

## Overview

`CDRSelector` is a reusable utility for identifying TCR CDR loops in MD trajectories using exact sequence matching. This tool can be used across different analysis modules for consistent CDR selection.

## Installation

The `CDRSelector` class is part of the `aftermd.utils` module:

```python
from aftermd.utils import CDRSelector
```

## Basic Usage

### 1. Initialize the Selector

```python
from aftermd.utils import CDRSelector

# Initialize with CDR reference CSV file
selector = CDRSelector("input/standardizedpdbs/tcr_cdrs_output.csv")
```

### 2. Get CDR Regions for a Specific PDB

```python
# Get all CDR regions for PDB 1AO7
cdr_info = selector.get_cdr_regions(
    pdb_file="input/pbc_1000frames_2step/1ao7_run1/md_converted.pdb",
    pdb_id="1AO7"
)

# Access CDR information
alpha_cdr3 = cdr_info['alpha']['cdr3']
print(f"Chain: {alpha_cdr3['chain_id']}")
print(f"Sequence: {alpha_cdr3['sequence']}")
print(f"Residues: {alpha_cdr3['residue_range']}")  # [start, end]
print(f"CA atoms: {alpha_cdr3['ca_count']}")
```

### 3. Create Selection Strings for MDAnalysis

```python
import MDAnalysis as mda

u = mda.Universe("md.tpr", "md.xtc")
cdr_info = selector.get_cdr_regions("md_converted.pdb", "1AO7")

# Get selection string for alpha CDR3
sel_str = selector.create_selection_string(cdr_info, 'alpha', 3)
# Returns: "chainID D and name CA and resid 89:99"

# Select atoms
alpha_cdr3_atoms = u.select_atoms(sel_str)
```

### 4. Get All CDR AtomGroups at Once

```python
# Get all CDR AtomGroups in one call
cdr_atoms = selector.get_all_cdr_atoms(u, cdr_info)

# Access specific CDRs
alpha_cdr1 = cdr_atoms['alpha_cdr1']
beta_cdr3 = cdr_atoms['beta_cdr3']

# Calculate RMSD for each CDR
from MDAnalysis.analysis import rms

for cdr_name, atoms in cdr_atoms.items():
    R = rms.RMSD(atoms, atoms, select='name CA')
    R.run()
    print(f"{cdr_name}: mean RMSD = {R.results.rmsd[:, 2].mean():.3f} nm")
```

## Advanced Usage

### Example 1: CDR Flexibility Analysis

```python
from aftermd.utils import CDRSelector
import MDAnalysis as mda
import numpy as np

# Initialize
selector = CDRSelector("tcr_cdrs_output.csv")
u = mda.Universe("md.tpr", "md.xtc")

# Get CDR regions
cdr_info = selector.get_cdr_regions("md_converted.pdb", "1AO7")
cdr_atoms = selector.get_all_cdr_atoms(u, cdr_info)

# Calculate RMSF for each CDR
from MDAnalysis.analysis import rms

results = {}
for cdr_name, atoms in cdr_atoms.items():
    rmsf = rms.RMSF(atoms).run()
    results[cdr_name] = {
        'mean_rmsf': rmsf.results.rmsf.mean(),
        'max_rmsf': rmsf.results.rmsf.max()
    }

# Find most flexible CDR
most_flexible = max(results.items(), key=lambda x: x[1]['mean_rmsf'])
print(f"Most flexible CDR: {most_flexible[0]}")
```

### Example 2: CDR Contact Analysis

```python
from aftermd.utils import CDRSelector
import MDAnalysis as mda
from MDAnalysis.analysis import contacts

# Setup
selector = CDRSelector("tcr_cdrs_output.csv")
u = mda.Universe("md.tpr", "md.xtc")
cdr_info = selector.get_cdr_regions("md_converted.pdb", "1AO7")

# Select peptide (assuming chain C)
peptide = u.select_atoms("chainID C")

# Calculate contacts between each CDR and peptide
for chain_type in ['alpha', 'beta']:
    for cdr_num in [1, 2, 3]:
        sel_str = selector.create_selection_string(cdr_info, chain_type, cdr_num)
        if sel_str:
            cdr_atoms = u.select_atoms(sel_str)

            # Calculate distance-based contacts
            ca = contacts.Contacts(
                u,
                selection=(sel_str, "chainID C"),
                refgroup=(cdr_atoms, peptide),
                radius=4.5
            )
            ca.run()

            # Get contact frequency
            contact_freq = ca.results.timeseries[:, 1].mean()
            print(f"{chain_type.upper()} CDR{cdr_num}: {contact_freq:.2%} contact frequency")
```

### Example 3: Batch Processing Multiple Systems

```python
from aftermd.utils import CDRSelector
import MDAnalysis as mda
from pathlib import Path
import pandas as pd

# Initialize selector
selector = CDRSelector("tcr_cdrs_output.csv")

# Process multiple systems
input_dir = Path("input/pbc_1000frames_2step")
results = []

for task_dir in input_dir.glob("*_run1"):
    pdb_id = task_dir.name[:4].upper()
    pdb_file = task_dir / "md_converted.pdb"

    if not pdb_file.exists():
        continue

    try:
        # Get CDR regions
        cdr_info = selector.get_cdr_regions(str(pdb_file), pdb_id)

        # Count found CDRs
        n_cdrs = sum(1 for chain in ['alpha', 'beta']
                     for cdr in ['cdr1', 'cdr2', 'cdr3']
                     if cdr in cdr_info[chain])

        results.append({
            'pdb_id': pdb_id,
            'n_cdrs_found': n_cdrs,
            'status': 'success'
        })

    except Exception as e:
        results.append({
            'pdb_id': pdb_id,
            'n_cdrs_found': 0,
            'status': f'failed: {e}'
        })

# Save results
df = pd.DataFrame(results)
df.to_csv("cdr_detection_summary.csv", index=False)
print(f"Processed {len(results)} systems")
print(f"Average CDRs found: {df['n_cdrs_found'].mean():.1f}/6")
```

### Example 4: Generate Custom Index Files

```python
from aftermd.utils import CDRSelector
import MDAnalysis as mda

# Setup
selector = CDRSelector("tcr_cdrs_output.csv")
u = mda.Universe("md.tpr")
cdr_info = selector.get_cdr_regions("md_converted.pdb", "1AO7")

# Generate GROMACS index file for all CDRs
with open("cdr_loops.ndx", 'w') as f:
    for chain_type in ['alpha', 'beta']:
        for cdr_num in [1, 2, 3]:
            cdr_key = f'cdr{cdr_num}'
            if cdr_key in cdr_info[chain_type]:
                cdr_data = cdr_info[chain_type][cdr_key]

                # Write group header
                group_name = f"CDR{cdr_num}_{chain_type}_CA"
                f.write(f"[ {group_name} ]\n")

                # Write atom indices (10 per line)
                indices = cdr_data['ca_indices']
                for i in range(0, len(indices), 10):
                    f.write(" ".join(map(str, indices[i:i+10])) + "\n")
                f.write("\n")

print("Index file generated: cdr_loops.ndx")
```

## Return Data Structure

The `get_cdr_regions()` method returns a nested dictionary:

```python
{
    'alpha': {
        'cdr1': {
            'chain_id': 'D',                 # Chain where CDR was found
            'sequence': 'DRGSQS',            # CDR sequence
            'residue_range': [26, 31],       # 1-based residue numbers
            'residue_count': 6,              # Number of residues
            'ca_count': 6,                   # Number of CA atoms
            'ca_indices': [3201, 3209, ...]  # 1-based atom indices
        },
        'cdr2': {...},
        'cdr3': {...}
    },
    'beta': {
        'cdr1': {...},
        'cdr2': {...},
        'cdr3': {...}
    }
}
```

## Key Methods

### `get_cdr_regions(pdb_file, pdb_id)`
Main method to get all CDR information for a PDB.

### `create_selection_string(cdr_info, chain_type, cdr_number, atom_name='CA')`
Generate MDAnalysis selection string for a specific CDR.

### `get_all_cdr_atoms(universe, cdr_info, atom_name='CA')`
Get MDAnalysis AtomGroups for all CDRs at once.

### `get_chain_sequences(pdb_file)`
Extract sequences from all chains (useful for debugging).

### `find_sequence_in_chain(chain_sequence, target_sequence)`
Low-level method to locate a sequence within a chain.

## Error Handling

```python
from aftermd.utils import CDRSelector

selector = CDRSelector("tcr_cdrs_output.csv")

try:
    cdr_info = selector.get_cdr_regions("md.pdb", "XXXX")
except ValueError as e:
    print(f"PDB not found in reference: {e}")

# Check if specific CDR was found
cdr_info = selector.get_cdr_regions("md.pdb", "1AO7")

if 'cdr3' in cdr_info['alpha']:
    print("Alpha CDR3 found!")
else:
    print("Alpha CDR3 not detected")
```

## Integration with Existing AfterMD Modules

### With RMSD Analysis

```python
from aftermd.analysis.trajectory import RMSDCalculator
from aftermd.utils import CDRSelector

selector = CDRSelector("tcr_cdrs_output.csv")
cdr_info = selector.get_cdr_regions("md.pdb", "1AO7")

# Create selection for specific CDR
sel_str = selector.create_selection_string(cdr_info, 'beta', 3)

# Calculate RMSD
calc = RMSDCalculator("md.tpr", "md.xtc")
times, rmsd = calc.calculate_rmsd(selection=sel_str)
```

### With Contact Analysis

```python
from aftermd.analysis.interface import ContactResidueAnalyzer
from aftermd.utils import CDRSelector

selector = CDRSelector("tcr_cdrs_output.csv")
cdr_info = selector.get_cdr_regions("md.pdb", "1AO7")

# Analyze contacts between CDR3 and peptide
cdr3_sel = selector.create_selection_string(cdr_info, 'alpha', 3, atom_name='all')

analyzer = ContactResidueAnalyzer("md.tpr", "md.xtc")
contacts = analyzer.calculate_contact_frequency(
    selection_a=cdr3_sel,
    selection_b="chainID C",  # peptide
    cutoff=4.0
)
```

## Notes

1. **Chain Independence**: The selector searches all chains to find CDR sequences, so the actual chain ID in md_converted.pdb doesn't need to match the original PDB.

2. **Exact Matching**: Uses exact string matching for CDR sequences. Sequences must match perfectly (including any modifications in the reference CSV).

3. **Reference Data**: Requires a CSV file with columns: `PDB_ID`, `alpha_cdr1`, `alpha_cdr2`, `alpha_cdr3`, `beta_cdr1`, `beta_cdr2`, `beta_cdr3`.

4. **Performance**: Chain sequence extraction is cached per PDB file, making repeated queries efficient.

## See Also

- `scripts/batch_cdr_rmsd_exact.py` - Full batch analysis implementation
- `aftermd/utils/cdr_manager.py` - Original ANARCI-based CDR detection
- `docs/RMSD_CALCULATION_GUIDE.md` - RMSD analysis documentation
