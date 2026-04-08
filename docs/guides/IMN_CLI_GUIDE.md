# IMN - Immunex CLI Usage Guide

**IMN** (Immunex) is the unified command-line interface for the Immunex MD Analysis Toolkit.

## Installation

```bash
cd /path/to/Immunex
pip install -e .
```

After installation, the `imn` command will be available globally.

## Quick Start

```bash
# Show help
imn --help

# Show version
imn --version

# Get help for specific subcommand
imn preprocess --help
imn quality --help
```

---

## Core Commands

### 1. Preprocess - PBC Correction

Remove periodic boundary condition artifacts from MD trajectories.

```bash
# Basic usage (interactive method selection)
imn preprocess -f md.xtc -s md.tpr -o processed.xtc

# Specify 2-step method directly
imn preprocess -f md.xtc -s md.tpr -o processed.xtc -m 2step

# Use 3-step method
imn preprocess -f md.xtc -s md.tpr -o processed.xtc -m 3step

# With custom time step
imn preprocess -f md.xtc -s md.tpr -o processed.xtc -m 2step --dt 100

# Verbose output
imn preprocess -f md.xtc -s md.tpr -o processed.xtc -v
```

**PBC Methods**:
- `2step`: `nojump` + `fit rot+trans`（推荐，速度更快）
- `3step`: `center` + `whole` + `fit`（更传统，也更彻底）

---

### 2. Quality - Quality Assessment

Assess trajectory quality with Post-PBC validation and RMSD convergence analysis.

```bash
# Basic quality assessment
imn quality -f processed.xtc -s md.tpr

# Save report to specific directory
imn quality -f processed.xtc -s md.tpr -o ./quality_reports

# Custom RMSD selection
imn quality -f processed.xtc -s md.tpr --selection "protein and name CA"

# Fast mode (skip validation)
imn quality -f processed.xtc -s md.tpr --skip-validation

# Choose report format
imn quality -f processed.xtc -s md.tpr --format markdown
```

**Quality Grades**:
- **Grade A** (Excellent): Mean RMSD < 0.3 nm, Std < 0.05 nm
- **Grade B** (Good): Mean RMSD < 0.5 nm, Std < 0.1 nm
- **Grade C** (Acceptable): Mean RMSD < 0.8 nm, Std < 0.15 nm
- **Grade D** (Poor): Mean RMSD > 0.8 nm or Std > 0.15 nm

---

### 3. RMSD - RMSD Calculation

Calculate Root Mean Square Deviation.

```bash
# Backbone RMSD
imn rmsd -f processed.xtc -s md.tpr --selection backbone

# C-alpha RMSD
imn rmsd -f processed.xtc -s md.tpr --selection "protein and name CA"

# Save to file
imn rmsd -f processed.xtc -s md.tpr --selection backbone -o rmsd.xvg

# Use GROMACS backend
imn rmsd -f processed.xtc -s md.tpr --selection backbone --backend gromacs
```

**Common Selections**:
- `backbone` - Protein backbone (N, CA, C)
- `protein and name CA` - C-alpha atoms only
- `resid 1-50` - Specific residue range
- `chainID A` - Specific chain

---

### 4. RMSF - RMSF Calculation

Calculate Root Mean Square Fluctuation.

```bash
# Basic RMSF
imn rmsf -f processed.xtc -s md.tpr --selection backbone

# C-alpha RMSF
imn rmsf -f processed.xtc -s md.tpr --selection "protein and name CA" -o rmsf.xvg
```

---

### 5. PDB - PDB File Processing

Download, fix, and standardize PDB files.

```bash
# Download PDB files
imn pdb download 1ao7 2ckb 6vma

# Download to specific directory
imn pdb download 1ao7 -o ./pdbs

# Fix PDB structure
imn pdb fix -f input.pdb -o fixed.pdb

# Standardize chain IDs
imn pdb standardize -f input.pdb -o standardized.pdb

# Batch processing
imn pdb batch --pdb-list pdbs.txt --workers 4
```

---

### 6. Batch - Batch Processing

Process multiple MD tasks in parallel.

```bash
# Batch PBC processing
imn batch preprocess /data/md_tasks/ --workers 4

# Batch quality assessment
imn batch quality /data/md_tasks/ --workers 4 --summary results.json

# Complete workflow (PBC + Quality)
imn batch workflow /data/md_tasks/ --workers 4
```

---

### 7. Contact - Contact Analysis

Analyze inter-molecular contacts (under development).

```bash
# Contact analysis
imn contact -f processed.xtc -s md.tpr --sel1 "chainID A" --sel2 "chainID B"
```

---

### 8. Angle - Docking Angle Analysis

Analyze TCR-pMHC docking angles (under development).

```bash
# Docking angle analysis
imn angle -f processed.xtc -s md.tpr --type docking
```

---

## Complete Workflow Example

```bash
# Step 1: PBC correction (interactive method selection)
imn preprocess -f md.xtc -s md.tpr -o processed.xtc

# Step 2: Quality assessment
imn quality -f processed.xtc -s md.tpr -o ./quality_reports

# Step 3: RMSD analysis
imn rmsd -f processed.xtc -s md.tpr --selection backbone -o rmsd.xvg

# Step 4: RMSF analysis
imn rmsf -f processed.xtc -s md.tpr --selection backbone -o rmsf.xvg
```

---

## Common Options

| Option | Description |
|--------|-------------|
| `-f`, `--trajectory` | Input trajectory file |
| `-s`, `--topology` | Structure/topology file (e.g., md.tpr) - GROMACS standard |
| `-o`, `--output` | Output file or directory |
| `-m`, `--method` | Processing method (for preprocess: 2step/3step) |
| `--selection` | Atom selection (MDAnalysis syntax) |
| `-v`, `--verbose` | Verbose output |
| `--help` | Show help message |

---

## Environment Variables

```bash
# Set GROMACS executable
export GMX=/path/to/gmx

# Set number of threads
export OMP_NUM_THREADS=4
```

---

## Troubleshooting

### Command not found: imn

Make sure Immunex is installed:
```bash
pip install -e /path/to/Immunex
```

Check installation:
```bash
which imn
imn --version
```

### Import errors

Ensure all dependencies are installed:
```bash
pip install -r requirements.txt
```

### GROMACS not found

Specify GROMACS executable:
```bash
imn preprocess -f md.xtc -s md.tpr -o processed.xtc --gmx /path/to/gmx
```

---

## Getting Help

```bash
# Main help
imn --help

# Subcommand help
imn preprocess --help
imn quality --help
imn rmsd --help
imn pdb --help
imn batch --help
```

---

## Version Information

```bash
imn --version
```

---

**Author**: Immunex Development Team
**Date**: 2026-03-15
**Version**: 0.1.0
