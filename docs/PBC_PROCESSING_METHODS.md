# PBC Processing Methods Guide

## Overview

AfterMD now supports two PBC processing methods based on different use cases and empirical experience.

## Method Comparison

### 2-Step Method (Recommended, Default)

**Process Flow:**
```bash
Step 1: gmx trjconv -pbc nojump  # Prevent atom jumps across boundaries
Step 2: gmx trjconv -fit rot+trans  # Align structures
```

**Advantages:**
- Simple and effective for large complexes (e.g., TCR-pMHC)
- Avoids potential issues with centering and whole operations
- Based on empirical experience and validated workflows
- Faster processing time

**Recommended for:**
- Large protein complexes
- TCR-pMHC systems
- Systems where centering may cause artifacts
- Standard production MD analysis

---

### 3-Step Method (Standard GROMACS)

**Process Flow:**
```bash
Step 1: gmx trjconv -center -pbc nojump  # Center on specified group
Step 2: gmx trjconv -pbc whole  # Make molecules whole
Step 3: gmx trjconv -fit rot+trans  # Align structures
```

**Advantages:**
- Complete GROMACS standard workflow
- Better for systems requiring explicit centering
- Handles broken molecules across boundaries

**Recommended for:**
- Small to medium systems
- Systems with specific centering requirements
- When molecules are broken across boundaries

---

## Usage Examples

### Command Line (scripts/pbc_process.py)

#### Using 2-step method (default)
```bash
python scripts/pbc_process.py -f md.xtc -s md.tpr -o output/
```

#### Explicitly specify 2-step method
```bash
python scripts/pbc_process.py -f md.xtc -s md.tpr -o output/ --method 2step
```

#### Using 3-step method
```bash
python scripts/pbc_process.py -f md.xtc -s md.tpr -o output/ --method 3step
```

#### With custom parameters
```bash
# 2-step with custom dt
python scripts/pbc_process.py -f md.xtc -s md.tpr -o output/ \
    --method 2step --dt 50.0 --fit-group "Backbone"

# 3-step with custom groups
python scripts/pbc_process.py -f md.xtc -s md.tpr -o output/ \
    --method 3step --center-group "Protein" --fit-group "Backbone"
```

#### Batch processing
```bash
python scripts/pbc_process.py -d input_dir/ -o output_dir/ \
    --batch --method 2step --workers 4
```

---

### Python API

#### Using PBCProcessor directly

```python
from aftermd.core import PBCProcessor

# Initialize processor
processor = PBCProcessor(gmx_executable="gmx")

# 2-step method (recommended)
result = processor.comprehensive_pbc_process(
    trajectory="md.xtc",
    topology="md.tpr",
    output_dir="output/",
    method="2step",
    dt=100.0  # Optional: downsample to ~1000 frames
)

# 3-step method
result = processor.comprehensive_pbc_process(
    trajectory="md.xtc",
    topology="md.tpr",
    output_dir="output/",
    method="3step",
    center_group="Protein",
    fit_group="Backbone"
)
```

#### Using BatchProcessor for parallel processing

```python
from aftermd.utils import BatchProcessor

# Initialize batch processor
batch = BatchProcessor(max_workers=4)

# Batch PBC processing with 2-step method
results = batch.batch_pbc_processing(
    base_directory="input_tasks/",
    output_base_dir="output_tasks/",
    method="2step",
    dt=100.0
)

print(f"Processed {results['successful']}/{results['total_tasks']} tasks")
```

---

## Migration from Old Workflow

### Before (3-step only)
```python
processor.comprehensive_pbc_process(
    trajectory="md.xtc",
    topology="md.tpr",
    output_dir="output/"
)
```

### After (explicit 2-step, recommended)
```python
processor.comprehensive_pbc_process(
    trajectory="md.xtc",
    topology="md.tpr",
    output_dir="output/",
    method="2step"  # Now explicit
)
```

**Note:** The default method is now "2step", so existing code will automatically use the recommended method unless explicitly specified.

---

## Method Selection Guide

| Scenario | Recommended Method | Reason |
|----------|-------------------|--------|
| TCR-pMHC complexes | 2-step | Validated empirically, avoids centering artifacts |
| Antibody-antigen | 2-step | Large complex, similar to TCR-pMHC |
| Protein-protein (large) | 2-step | Simple and effective |
| Protein-ligand (small) | 3-step | Standard workflow sufficient |
| Membrane proteins | 3-step | May need explicit centering |
| Custom centering needed | 3-step | Provides center control |

---

## Default Parameters

| Parameter | 2-step Default | 3-step Default |
|-----------|---------------|---------------|
| `dt` | 100 ps (auto) | 100 ps (auto) |
| `fit_group` | Backbone | Backbone |
| `output_group` | System | System |
| `center_group` | N/A | Protein |
| `use_nojump` | Always on | User controlled |

---

## Output Files

Both methods produce the same output structure:

```
output/
├── md_processed.xtc           # Processed trajectory
├── md.tpr                      # Copied reference topology
├── md.gro                      # Copied reference structure
└── (temporary files cleaned)
```

**Result dictionary:**
```python
{
    "processed": "output/md_processed.xtc",
    "method": "2step",  # or "3step"
    "fit_group": "Backbone",
    "output_directory": "output/",
    "reference_structures": ["output/md.gro"],
    "reference_topology": "output/md.tpr"
}
```

---

## Performance Comparison

| Aspect | 2-step | 3-step |
|--------|--------|--------|
| Processing time | ~30-40s per 100ns | ~45-60s per 100ns |
| Temporary files | 1 intermediate | 2 intermediates |
| Memory usage | Lower | Moderate |
| RMSD quality | Excellent | Excellent |

---

## Troubleshooting

### Issue: High RMSD values with 3-step method
**Solution:** Try 2-step method, which avoids centering artifacts

### Issue: Molecules broken across boundaries
**Solution:** Use 3-step method with `use_nojump=True`

### Issue: Need custom centering
**Solution:** Use 3-step method with `center_group` parameter

---

## References

- Original 2-step implementation: `development/archived_scripts/process_2step_pbc.py`
- PBCProcessor source: `aftermd/core/pbc_processor.py`
- CLI script: `scripts/pbc_process.py`
