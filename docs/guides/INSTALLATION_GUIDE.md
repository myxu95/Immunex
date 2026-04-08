# Immunex Installation Guide

Complete installation guide for the Immunex MD analysis toolkit.

## Quick Start

### 1. Create Conda Environment

```bash
# Navigate to Immunex directory
cd /home/xumy/work/development/Immunex

# Run setup script
bash setup_immunex_env.sh

# Activate environment
conda activate immunex
```

### 2. Install Immunex Package

```bash
# Install in development mode
pip install -e .
```

### 3. Test Installation

```bash
# Run environment test
python development/validation/test_environment.py
```

---

## Detailed Installation Steps

### Prerequisites

- **Anaconda or Miniconda**: [Download here](https://docs.conda.io/en/latest/miniconda.html)
- **Git**: For cloning the repository
- **Linux/macOS**: Recommended (Windows via WSL2)

### Step 1: Clone Repository (if needed)

```bash
git clone https://github.com/your-org/Immunex.git
cd Immunex
```

### Step 2: Create Conda Environment

The environment includes all required dependencies:

- **Python 3.10**
- **MDAnalysis** - MD trajectory analysis
- **GROMACS** - MD simulation engine
- **ANARCI** - CDR loop detection (critical)
- **BioPython** - Structure manipulation
- **NumPy, Pandas, SciPy** - Scientific computing
- **Matplotlib, Seaborn, Plotly** - Visualization

#### Option A: Automated Setup (Recommended)

```bash
bash setup_immunex_env.sh
```

#### Option B: Manual Setup

```bash
# Create environment
conda env create -f environment.yml

# Activate environment
conda activate immunex

# Install Immunex package
pip install -e .
```

### Step 3: Verify Installation

```bash
conda activate immunex
python development/validation/test_environment.py
```

Expected output:
```
Testing Critical Dependencies
======================================================================

1. Core Scientific Libraries
----------------------------------------------------------------------
  NumPy                - OK              (v1.24.3)
  Pandas               - OK              (v2.0.2)
  SciPy                - OK              (v1.10.1)
  Matplotlib           - OK              (v3.7.1)

2. MD Analysis Tools
----------------------------------------------------------------------
  MDAnalysis           - OK              (v2.6.1)

3. Structure Analysis
----------------------------------------------------------------------
  BioPython            - OK              (v1.81)
  Biotite              - OK              (v0.37.0)

4. CDR Detection (CRITICAL)
----------------------------------------------------------------------
  ANARCI               - OK              (v1.3)
  ANARCI functional test: PASSED

...

✓ All critical dependencies installed successfully!
```

---

## Testing CDR Detection

Once the environment is set up, test CDR detection:

```bash
conda activate immunex

# Test with real PDB file
python development/validation/test_cdr_detection.py
```

Expected output:
```
Detecting CDR loops with ANARCI...
----------------------------------------------------------------------
TCR-alpha (chain D): 115 residues
  CDR1: pos 26-33, seq=SGDSAVN
  CDR2: pos 51-57, seq=DTQADS
  CDR3: pos 93-101, seq=AASIRSSYK

TCR-beta (chain E): 209 residues
  CDR1: pos 26-31, seq=SGHNS
  CDR2: pos 49-55, seq=FNNNVP
  CDR3: pos 92-102, seq=ASSLAPGTTN

Total CDR groups added: 12
```

---

## Generate Comprehensive Index

Test the full index generation with all components:

```bash
# Generate complete index with standard components + backbones
python development/validation/test_comprehensive_index.py

# Add CDR loops (requires ANARCI)
python development/validation/test_cdr_detection.py
```

Generated index will include:
- Standard components (pHLA, TCR, peptide, etc.)
- Backbone variants for all components
- CDR loops (CDR1/2/3 for alpha and beta chains)
- CDR backbone variants

---

## Troubleshooting

### ANARCI Not Found

If ANARCI installation fails:

```bash
# Install separately
conda activate immunex
conda install -c bioconda anarci

# Verify installation
python -c "from anarci import anarci; print('ANARCI OK')"
```

### GROMACS Not Found

If you have GROMACS installed system-wide:

```bash
# Remove GROMACS from environment.yml
# Then recreate environment

# Or install separately
conda install -c bioconda gromacs
```

### Import Errors

```bash
# Reinstall Immunex package
pip install -e . --force-reinstall

# Check installation
python -c "import immunex; print(immunex.__file__)"
```

### Environment Conflicts

```bash
# Remove and recreate environment
conda env remove -n immunex
conda env create -f environment.yml
conda activate immunex
pip install -e .
```

---

## Environment Management

### Update Environment

```bash
# After modifying environment.yml
conda env update -f environment.yml --prune
```

### Export Environment

```bash
# Export exact environment
conda env export > environment_exact.yml

# Export cross-platform environment
conda env export --from-history > environment_minimal.yml
```

### Remove Environment

```bash
conda env remove -n immunex
```

---

## Alternative Installation (pip only)

If you don't want to use conda:

```bash
# Create virtual environment
python -m venv immunex_venv
source immunex_venv/bin/activate

# Install dependencies
pip install -r requirements.txt

# Note: ANARCI requires HMMER which needs conda/system install
# Manual ANARCI installation: https://github.com/oxpig/ANARCI
```

---

## Verification Checklist

- [ ] Conda environment created successfully
- [ ] `conda activate immunex` works
- [ ] `python development/validation/test_environment.py` passes
- [ ] ANARCI import works: `python -c "from anarci import anarci"`
- [ ] GROMACS available: `gmx --version`
- [ ] Immunex imports: `python -c "import immunex"`
- [ ] CDR detection works: `python development/validation/test_cdr_detection.py`

---

## Next Steps

Once installation is complete:

1. **Test index generation**:
   ```bash
   python development/quick_index_test.py
   ```

2. **Try examples**:
   ```bash
   python examples/index_generation_usage.py
   ```

3. **Read documentation**:
   - `docs/INDEX_MANAGER_GUIDE.md`
   - `docs/architecture/NEW_ARCHITECTURE_QUICKSTART.md`

4. **Start analyzing**:
   ```bash
   # Your own analysis
   python your_analysis_script.py
   ```

---

## Support

If you encounter issues:

1. Check `development/validation/test_environment.py` output
2. Verify all dependencies installed correctly
3. Check ANARCI specifically: `python -c "from anarci import anarci"`
4. Create GitHub issue with error message and environment details

---

## System Requirements

- **OS**: Linux (Ubuntu 18.04+), macOS (10.14+), Windows (WSL2)
- **RAM**: 8GB minimum, 16GB recommended
- **Disk**: 10GB for conda environment + data
- **CPU**: Multi-core recommended for parallel processing

---

**Installation Date**: 2026-03-18  
**Immunex Version**: Development  
**Conda Environment**: immunex
