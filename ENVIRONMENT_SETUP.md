# AfterMD Environment Setup Guide

## Quick Start

### 1. Create Conda Environment

```bash
# Navigate to AfterMD directory
cd /home/xumy/work/development/AfterMD

# Create environment from yml file
conda env create -f environment.yml

# This will create an environment named 'aftermd'
```

### 2. Activate Environment

```bash
conda activate aftermd
```

### 3. Install AfterMD Package

```bash
# Install in development mode (recommended for development)
pip install -e .

# OR install in regular mode
pip install .
```

### 4. Verify Installation

```bash
# Test import
python -c "import aftermd; print('AfterMD version:', aftermd.__version__)"

# Test EnergyQualityChecker
python -c "from aftermd.analysis.quality import EnergyQualityChecker; print('Success!')"

# Run comprehensive test
python test_energy_standalone.py
```

---

## Environment Details

### Created Environment

**Name**: `aftermd`
**Python Version**: 3.9
**Location**: `~/work/miniconda3/envs/aftermd`

### Installed Packages

**Core Scientific Computing**:
- numpy >= 1.24.0
- pandas >= 2.0.0
- scipy >= 1.10.0
- matplotlib >= 3.7.0

**MD Analysis**:
- MDAnalysis >= 2.6.0 (for trajectory analysis)

**Visualization**:
- seaborn >= 0.12.0
- plotly >= 5.15.0

**Development Tools** (optional):
- pytest >= 7.0
- black >= 23.0
- flake8 >= 6.0

**System Tools**:
- git
- wget
- curl

---

## Usage

### Activate Environment

Every time you work with AfterMD, activate the environment first:

```bash
conda activate aftermd
```

### Run Scripts

```bash
# Activate environment
conda activate aftermd

# Run quality check
python scripts/md_quality_check.py /path/to/md/simulations

# Run PBC processing
python scripts/pbc_process.py /path/to/trajectory

# Use examples
python examples/energy_quality_check_usage.py
```

### Use in Python Scripts

```python
#!/usr/bin/env python3
# Make sure to run with: conda activate aftermd

from aftermd.analysis.quality import EnergyQualityChecker

# Your code here
checker = EnergyQualityChecker()
result = checker.comprehensive_energy_check("md.edr")
print(f"Energy Grade: {result['energy_grade']}")
```

---

## Running Without Activation

If you need to run scripts without manually activating the environment:

```bash
# Method 1: Use conda run
conda run -n aftermd python scripts/md_quality_check.py /path/to/data

# Method 2: Use full path to python
~/work/miniconda3/envs/aftermd/bin/python scripts/md_quality_check.py /path/to/data
```

---

## Environment Management

### List All Environments

```bash
conda env list
```

### Deactivate Environment

```bash
conda deactivate
```

### Remove Environment (if needed)

```bash
conda env remove -n aftermd
```

### Update Environment

If you modify `environment.yml`:

```bash
conda env update -f environment.yml --prune
```

### Export Environment

To share your exact environment:

```bash
# Export with exact versions
conda env export > environment_exact.yml

# Export with minimal versions
conda env export --from-history > environment_minimal.yml
```

---

## Requirements Files

### Two Installation Methods

**1. Using conda (recommended)**:
```bash
conda env create -f environment.yml
```

**2. Using pip only**:
```bash
pip install -r requirements.txt
pip install -e .
```

### Updating Requirements

If you add new dependencies:

1. Update `environment.yml` (conda packages)
2. Update `requirements.txt` (pip packages)
3. Update `setup.py` (package dependencies)

---

## Troubleshooting

### Issue: "ModuleNotFoundError: No module named 'aftermd'"

**Solution**:
```bash
# Make sure environment is activated
conda activate aftermd

# Reinstall package
pip install -e .
```

### Issue: "ModuleNotFoundError: No module named 'MDAnalysis'"

**Solution**:
```bash
# Activate environment
conda activate aftermd

# Reinstall MDAnalysis
conda install -c conda-forge mdanalysis
```

### Issue: GROMACS not found

AfterMD requires GROMACS for PBC processing. Install separately:

```bash
# On Ubuntu/Debian
sudo apt-get install gromacs

# Or use conda
conda install -c conda-forge gromacs

# Verify installation
gmx --version
```

### Issue: ImportError with Bio.Application

This is a deprecation warning from BioPython, not an error. Your code will work fine.
To suppress the warning:

```python
import warnings
warnings.filterwarnings('ignore', category=DeprecationWarning)
```

---

## Development Workflow

### Typical Development Session

```bash
# 1. Activate environment
conda activate aftermd

# 2. Make code changes
vim aftermd/analysis/quality/energy_quality.py

# 3. Test changes
python test_energy_standalone.py

# 4. Run your analysis
python scripts/md_quality_check.py /path/to/data

# 5. Deactivate when done
conda deactivate
```

### Running Tests

```bash
conda activate aftermd

# Run standalone test
python test_energy_standalone.py

# Run all examples (if implemented)
python examples/energy_quality_check_usage.py
```

---

## System Requirements

**Minimum**:
- Python 3.8+
- 4 GB RAM
- 1 GB disk space

**Recommended**:
- Python 3.9+
- 8 GB RAM
- 5 GB disk space (for conda cache)
- GROMACS 2020+ (for PBC processing)

---

## Quick Reference Commands

```bash
# Create environment
conda env create -f environment.yml

# Activate
conda activate aftermd

# Install AfterMD
pip install -e .

# Test installation
python -c "import aftermd; print(aftermd.__version__)"

# Run test
python test_energy_standalone.py

# Deactivate
conda deactivate

# List environments
conda env list

# Remove environment
conda env remove -n aftermd
```

---

## Environment Status: ✅ Ready

- ✅ Conda environment created: `aftermd`
- ✅ All dependencies installed
- ✅ AfterMD package installed in development mode
- ✅ Import tests passed
- ✅ EnergyQualityChecker module tested and working

**You can now start using AfterMD!**

```bash
conda activate aftermd
python scripts/md_quality_check.py /path/to/your/data
```
