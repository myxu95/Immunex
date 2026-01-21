# AfterMD Quick Start

## Environment Setup (One-Time)

```bash
# Create and activate environment
conda env create -f environment.yml
conda activate aftermd

# Install AfterMD
pip install -e .

# Verify installation
python check_environment.py
```

## Daily Usage

```bash
# Always activate environment first
conda activate aftermd

# Run your analysis
python scripts/md_quality_check.py /path/to/md/simulations

# When done
conda deactivate
```

## Common Commands

### Energy Quality Check (Pre-PBC)
```bash
conda activate aftermd
python -c "
from aftermd.analysis.quality import EnergyQualityChecker

checker = EnergyQualityChecker()
result = checker.comprehensive_energy_check('/path/to/md.edr')

print(f\"Grade: {result['energy_grade']}\")
print(f\"Score: {result['score']:.1f}/100\")
"
```

### PBC Processing
```bash
conda activate aftermd
python scripts/pbc_process.py /path/to/trajectory.xtc
```

### Batch Processing
```bash
conda activate aftermd
python scripts/md_quality_check.py /data/simulations \
    --expected-chains 5 \
    --min-traj-size 1.0 \
    --min-sim-time 5000
```

## Environment Status

✅ **Environment Name**: `aftermd`  
✅ **Python Version**: 3.9.23  
✅ **AfterMD Version**: 0.1.0  
✅ **GROMACS**: 2025.3  
✅ **All Dependencies**: Installed  

## Get Help

```bash
# Check environment status
python check_environment.py

# Test energy checker
python test_energy_standalone.py

# View documentation
cat ARCHITECTURE.md
cat ENVIRONMENT_SETUP.md
```
