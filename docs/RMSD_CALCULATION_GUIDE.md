# RMSD Calculation System

## Overview

AfterMD provides a unified RMSD calculation system with flexible component selection for alignment and analysis. The system automatically handles chain standardization and index file generation.

## Architecture

### Core Components

1. **IndexGenerator** (`aftermd/utils/index_generator.py`)
   - Automatically standardizes PDB chain ordering
   - Generates GROMACS index files for different components
   - Supports: HLA, pHLA, peptide, TCR, TCR_alpha, TCR_beta, CDR3_alpha, CDR3_beta

2. **PDBChainStandardizer** (`aftermd/utils/pdb_chain_standardizer.py`)
   - Standardizes chain IDs based on residue count
   - Standard order: C (peptide) < B (HLA-β) < D (TCR-α) < E (TCR-β) < A (HLA-α)
   - Handles TPR/GRO to PDB conversion

3. **RMSDInterface** (`aftermd/analysis/trajectory/rmsd_interface.py`)
   - High-level API for RMSD calculations
   - **Unified index strategy**: Generates one base index file with all fixed components
   - CDR3 components generated separately (require sequences)
   - Automatic index management and caching

### Index File Strategy

The system uses a **hybrid index approach** for optimal efficiency:

#### Fixed Components (Generated Once)
At initialization, a single unified index file is created containing:
- `HLA` (chains A, B)
- `pHLA` (chains A, B, C)
- `peptide` (chain C)
- `TCR` (chains D, E)
- `TCR_alpha` (chain D)
- `TCR_beta` (chain E)

**File**: `base_components.ndx`

#### Dynamic Components (Generated On-Demand)
CDR3 regions require amino acid sequences and are generated separately:
- `CDR3_alpha` (requires α chain sequence)
- `CDR3_beta` (requires β chain sequence)

**Files**: `cdr3_alpha_index.ndx`, `cdr3_beta_index.ndx`

#### Benefits of This Approach
✅ **Efficient**: Fixed components generated once, reused for all calculations
✅ **Flexible**: CDR3 components support different sequences
✅ **Clean**: Minimal file clutter (1-3 index files instead of 8+)
✅ **Fast**: No redundant index generation

### Standard Chain Ordering

After standardization, chains are ordered by residue count:

- **Chain C**: Peptide (shortest, 5-22 residues)
- **Chain B**: HLA-β/β2-microglobulin (~99-100 residues)
- **Chain D**: TCR-α (115-215 residues)
- **Chain E**: TCR-β (209-252 residues)
- **Chain A**: HLA-α (longest, ~270-280 residues)

## Usage

### 1. Command Line Interface

#### Basic Usage

```bash
# Calculate TCR RMSD aligned to pHLA (most common use case)
python scripts/calculate_rmsd.py md.tpr md_processed.xtc --align pHLA --calc TCR

# Calculate peptide RMSD aligned to pHLA
python scripts/calculate_rmsd.py md.tpr md_processed.xtc --align pHLA --calc peptide

# Calculate HLA RMSD aligned to TCR
python scripts/calculate_rmsd.py md.tpr md_processed.xtc --align TCR --calc HLA
```

#### CDR3 RMSD Calculation

```bash
# CDR3β RMSD aligned to TCR
python scripts/calculate_rmsd.py md.tpr md_processed.xtc \
    --align TCR \
    --calc CDR3_beta \
    --cdr3-beta CASSLGQAYEQYF

# CDR3α RMSD aligned to pHLA
python scripts/calculate_rmsd.py md.tpr md_processed.xtc \
    --align pHLA \
    --calc CDR3_alpha \
    --cdr3-alpha CAVRDDSNYQLIW
```

#### Multiple Calculations

```bash
# Calculate RMSD for multiple components
python scripts/calculate_rmsd.py md.tpr md_processed.xtc \
    --align pHLA \
    --calc TCR peptide TCR_alpha TCR_beta \
    --output-csv results.csv
```

#### Batch Mode with Config File

Create a JSON config file (`rmsd_config.json`):

```json
[
    {
        "align": "pHLA",
        "calc": "TCR",
        "output_file": "tcr_rmsd.xvg"
    },
    {
        "align": "pHLA",
        "calc": "peptide"
    },
    {
        "align": "TCR",
        "calc": "CDR3_beta",
        "cdr3_sequences": {"beta": "CASSLGQAYEQYF"}
    }
]
```

Run batch calculation:

```bash
python scripts/calculate_rmsd.py md.tpr md_processed.xtc \
    --batch-config rmsd_config.json \
    --output-csv batch_results.csv
```

### 2. Python API

#### Basic Usage

```python
from aftermd.analysis.trajectory import RMSDInterface

# Initialize interface
rmsd = RMSDInterface(
    topology="md.tpr",
    trajectory="md_processed.xtc",
    output_dir="rmsd_results"
)

# Calculate TCR RMSD aligned to pHLA
result = rmsd.calculate(align="pHLA", calc="TCR")

print(f"Mean RMSD: {result['mean']:.4f} nm")
print(f"Std RMSD: {result['std']:.4f} nm")
```

#### CDR3 Calculation

```python
# CDR3β RMSD aligned to TCR
result = rmsd.calculate(
    align="TCR",
    calc="CDR3_beta",
    cdr3_sequences={'beta': 'CASSLGQAYEQYF'}
)
```

#### Batch Calculations

```python
calculations = [
    {'align': 'pHLA', 'calc': 'TCR'},
    {'align': 'pHLA', 'calc': 'peptide'},
    {'align': 'TCR', 'calc': 'CDR3_beta',
     'cdr3_sequences': {'beta': 'CASSLGQAYEQYF'}}
]

df = rmsd.batch_calculate(
    calculations=calculations,
    output_csv="batch_results.csv"
)

print(df[['align_component', 'calc_component', 'rmsd_mean_nm']])
```

#### Available Components

```python
# List all available components
components = RMSDInterface.list_available_components()
print(components)
# Output: ['HLA', 'pHLA', 'peptide', 'TCR', 'TCR_alpha', 'TCR_beta', 'CDR3_alpha', 'CDR3_beta']
```

### 3. Direct Index Generator Usage

```python
from aftermd.utils import IndexGenerator

# Initialize with automatic chain standardization
generator = IndexGenerator(
    topology="md.tpr",
    auto_standardize=True
)

# Generate single component index
generator.generate_component_index(
    component="pHLA",
    output_file="phla.ndx"
)

# Generate multi-component index
generator.generate_multi_component_index(
    components=['pHLA', 'TCR', 'peptide'],
    output_file="combined.ndx"
)

# Generate CDR3 index
generator.generate_cdr3_index(
    chain='beta',
    cdr3_sequence='CASSLGQAYEQYF',
    output_file="cdr3_beta.ndx"
)
```

## Component Definitions

| Component | Chains | Description |
|-----------|--------|-------------|
| `HLA` | A, B | MHC alpha and beta chains |
| `pHLA` | A, B, C | pMHC complex (MHC + peptide) |
| `peptide` | C | Antigenic peptide |
| `TCR` | D, E | TCR alpha and beta chains |
| `TCR_alpha` | D | TCR alpha chain only |
| `TCR_beta` | E | TCR beta chain only |
| `CDR3_alpha` | D | CDR3 region of TCR alpha (sequence required) |
| `CDR3_beta` | E | CDR3 region of TCR beta (sequence required) |

## Common Use Cases

### 1. TCR Flexibility Analysis

Align to pHLA, calculate TCR RMSD to assess TCR flexibility relative to the pMHC:

```bash
python scripts/calculate_rmsd.py md.tpr md.xtc --align pHLA --calc TCR
```

### 2. CDR3 Flexibility Analysis

Align to TCR, calculate CDR3 RMSD to assess CDR3 loop flexibility:

```bash
python scripts/calculate_rmsd.py md.tpr md.xtc \
    --align TCR \
    --calc CDR3_alpha CDR3_beta \
    --cdr3-alpha CAVRDDSNYQLIW \
    --cdr3-beta CASSLGQAYEQYF
```

### 3. Peptide Stability Analysis

Align to pHLA, calculate peptide RMSD:

```bash
python scripts/calculate_rmsd.py md.tpr md.xtc --align pHLA --calc peptide
```

### 4. Interface Dynamics

Align to one component, calculate RMSD of the partner:

```bash
# Fix pHLA, see TCR movement
python scripts/calculate_rmsd.py md.tpr md.xtc --align pHLA --calc TCR

# Fix TCR, see pHLA movement
python scripts/calculate_rmsd.py md.tpr md.xtc --align TCR --calc pHLA
```

## Workflow Integration

### Typical Analysis Pipeline

```bash
# 1. PBC processing (prerequisite)
python scripts/pbc_process.py input/raw/md.xtc output/processed/

# 2. RMSD analysis with multiple alignments
python scripts/calculate_rmsd.py \
    input/md.tpr \
    output/processed/md_processed.xtc \
    --align pHLA \
    --calc TCR peptide \
    --output-dir output/rmsd_analysis \
    --output-csv rmsd_summary.csv

# 3. Visualization
python scripts/plot_rmsd.py output/rmsd_analysis/
```

## Output Files

### XVG Files

RMSD trajectory data in GROMACS XVG format:

```
# Time (ns)  RMSD (nm)
0.000       0.123
0.001       0.125
0.002       0.124
...
```

### CSV Summary

Batch calculation summary:

| align_component | calc_component | rmsd_mean_nm | rmsd_std_nm | rmsd_min_nm | rmsd_max_nm |
|----------------|----------------|--------------|-------------|-------------|-------------|
| pHLA | TCR | 0.234 | 0.045 | 0.156 | 0.389 |
| pHLA | peptide | 0.102 | 0.023 | 0.067 | 0.178 |

### Index Files

The system generates minimal index files:

**Basic Setup (Fixed Components Only)**:
```
rmsd_results/
├── base_components.ndx          # Unified index (HLA, pHLA, peptide, TCR, etc.)
├── md_converted.pdb             # Converted from TPR (if needed)
├── md_converted_standardized.pdb # Standardized chains (if needed)
├── rmsd_pHLA_to_TCR.xvg         # RMSD output
└── rmsd_pHLA_to_peptide.xvg     # RMSD output
```

**With CDR3 Analysis**:
```
rmsd_results/
├── base_components.ndx          # Fixed components
├── cdr3_alpha_index.ndx         # CDR3α (generated on-demand)
├── cdr3_beta_index.ndx          # CDR3β (generated on-demand)
├── combined_rmsd_index.ndx      # Merged (when using CDR3)
├── rmsd_pHLA_to_TCR.xvg
├── rmsd_TCR_to_CDR3_beta.xvg
└── rmsd_summary.csv
```

**Key Points**:
- `base_components.ndx`: Contains all 6 fixed components, generated once
- `cdr3_*_index.ndx`: Generated only when CDR3 analysis is requested
- `combined_rmsd_index.ndx`: Created only when mixing fixed + CDR3 components

## Troubleshooting

### Chain Order Issues

If you get incorrect chain assignments:

```python
# Check if PDB needs standardization
from aftermd.utils import PDBChainStandardizer

standardizer = PDBChainStandardizer()
chains = standardizer.analyze_chains("protein.pdb")

for chain in chains:
    print(f"Chain {chain.chain_id}: {chain.residue_count} residues")
```

### CDR3 Sequence Not Found

If CDR3 index generation fails:

- Verify the sequence matches the PDB structure
- Check for gaps or missing residues
- Use single-letter amino acid code

### Index Generation Errors

If index generation fails:

1. Check GROMACS installation: `gmx --version`
2. Verify topology file format
3. Enable verbose mode: `--verbose`

## Migration from Old Scripts

### Old Script Mapping

| Old Script | New Command |
|------------|-------------|
| `scripts/calculate_tcr_overall_rmsd.py` | `scripts/calculate_rmsd.py --align pHLA --calc TCR` |
| `scripts/batch_cdr3_rmsd.py` | `scripts/calculate_rmsd.py --align TCR --calc CDR3_beta --cdr3-beta SEQ` |
| `scripts/rmsd_calculator.py` | `scripts/calculate_rmsd.py` (with component flags) |

### Configuration File Migration

Old batch scripts using hard-coded paths can be replaced with JSON config files:

```json
{
    "tasks": [
        {
            "name": "1ao7_run1",
            "topology": "data/1ao7/md.tpr",
            "trajectory": "data/1ao7/md_processed.xtc",
            "calculations": [
                {"align": "pHLA", "calc": "TCR"},
                {"align": "TCR", "calc": "CDR3_beta",
                 "cdr3_sequences": {"beta": "CASSLGQAYEQYF"}}
            ]
        }
    ]
}
```

## Performance Considerations

### Index Caching

Index files are cached during a session. For large batch jobs, this significantly reduces overhead.

### Parallel Processing

For multiple independent calculations, use batch mode to parallelize:

```bash
# Process multiple tasks in parallel using GNU parallel
parallel -j 4 python scripts/calculate_rmsd.py {}/md.tpr {}/md.xtc --align pHLA --calc TCR ::: task_*
```

## References

- GROMACS `gmx rms` documentation
- MDAnalysis selection syntax
- AfterMD chain standardization protocol
