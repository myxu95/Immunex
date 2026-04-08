# Docking Angle Module Migration Guide

**Date**: 2026-02-04
**Version**: 3.0.0 (Major Breaking Changes)

## Overview

This guide helps users migrate from the old PCA-based docking angle module to the new geometric implementation following `docking_angles.md` specification.

---

## Summary of Changes

### ✅ What Changed

1. **Angle Definitions**:
   - **Old**: Twist/Tilt/Swing angles (PCA-based, relative to peptide)
   - **NEW**: Crossing/Incident angles (geometric, relative to MHC groove)

2. **MHC Region Detection**:
   - **Old**: Hardcoded residue ranges (resid 50:86 140:176)
   - **NEW**: Sequence alignment to HLA-A*02:01 reference (robust to non-standard numbering)

3. **TCR Axis Calculation**:
   - **Old**: Required manual specification of cysteine residues
   - **NEW**: Automatic via ANARCI (already implemented, no change)

4. **Module Names**:
   - **Old**: `DockingAngleAnalyzer` (deleted)
   - **NEW**: `DockingAnglePrimaryAnalyzer`

5. **Required Parameters**:
   - **Old**: `mhc_selection`, `tcr_alpha_selection`, `tcr_beta_selection`, `tcr_v_selection`, `peptide_selection`
   - **NEW**: `mhc_selection`, `tcr_alpha_selection`, `tcr_beta_selection` (simplified)

6. **Output Format**:
   - **Old**: Time, Twist(°), Tilt(°), Swing(°)
   - **NEW**: Time, Crossing(°), Incident(°)

---

## API Comparison

### Old API (❌ Deleted)

```python
from immunex.analysis.angles import DockingAngleAnalyzer

analyzer = DockingAngleAnalyzer('md.tpr', 'md.xtc')

times, twist, tilt, swing = analyzer.calculate_docking_angles_trajectory(
    mhc_selection='chainID A and resid 50:86 140:176 and name CA',  # Hardcoded
    tcr_alpha_selection='chainID D',
    tcr_beta_selection='chainID E',
    tcr_v_selection='chainID D E and resid 1:115 and name CA',     # Required
    peptide_selection='chainID C and name CA',                      # Required
    stride=10,
    output_file='docking_angles.csv'
)

# Output: Time, Twist(°), Tilt(°), Swing(°)
```

**Problems**:
- Hardcoded MHC residue ranges failed for non-standard PDB numbering
- Required specifying TCR V domain and peptide regions manually
- PCA-based angles lacked clear geometric interpretation
- Output included 3 angles (Twist/Tilt/Swing)

### New API (✅ Implemented)

```python
from immunex.analysis.angles import DockingAnglePrimaryAnalyzer

analyzer = DockingAnglePrimaryAnalyzer('md.tpr', 'md.xtc')

times, crossing, incident = analyzer.calculate_docking_angles_trajectory(
    mhc_selection='chainID A',           # Just chain ID - α1/α2 auto-detected
    tcr_alpha_selection='chainID D',     # ANARCI auto-detects cysteines
    tcr_beta_selection='chainID E',      # ANARCI auto-detects cysteines
    stride=10,
    output_file='docking_angles.csv'
)

# Output: Time, Crossing(°), Incident(°)
```

**Improvements**:
- MHC α1/α2 regions auto-detected via sequence alignment to HLA-A*02:01
- No need to specify TCR V domain or peptide regions
- Geometric angle definitions with clear physical meaning
- Simplified output (2 primary angles)

---

## Angle Name Changes

| Old Angle | New Angle | Definition |
|-----------|-----------|------------|
| **Twist** | **Crossing** | Lateral orientation of TCR across pHLA groove<br>`angle(TCR_axis, groove_axis)` |
| **Tilt** | **Incident** | Tilt/inclination of TCR relative to pHLA surface<br>`angle(TCR_axis, groove_normal)` |
| **Swing** | *(Removed)* | No longer calculated |

### Geometric Definitions

**Crossing Angle**:
```
θ_cross = arccos(v_TCR · v_groove / |v_TCR| |v_groove|)

where:
  v_TCR = COM(Vβ disulfide) - COM(Vα disulfide)
  v_groove = COM(α2, aligned) - COM(α1, aligned)
```

**Incident Angle**:
```
θ_inc = arccos(v_TCR · n_groove / |v_TCR| |n_groove|)

where:
  v_TCR = COM(Vβ disulfide) - COM(Vα disulfide)
  n_groove = normal from SVD plane fit of α1+α2 backbone
```

---

## Migration Examples

### Example 1: Single Frame Analysis

**Old Code**:
```python
from immunex.analysis.angles import DockingAngleAnalyzer

analyzer = DockingAngleAnalyzer('structure.pdb')
twist, tilt, swing = analyzer.calculate_docking_angles(
    mhc_selection='chainID A and resid 50:86 140:176 and name CA',
    tcr_alpha_selection='chainID D',
    tcr_beta_selection='chainID E',
    tcr_v_selection='chainID D E and resid 1:115 and name CA',
    peptide_selection='chainID C and name CA'
)
print(f"Twist: {twist:.2f}°, Tilt: {tilt:.2f}°, Swing: {swing:.2f}°")
```

**New Code**:
```python
from immunex.analysis.angles import DockingAnglePrimaryAnalyzer

analyzer = DockingAnglePrimaryAnalyzer('structure.pdb')
crossing, incident = analyzer.calculate_docking_angles(
    mhc_selection='chainID A',
    tcr_alpha_selection='chainID D',
    tcr_beta_selection='chainID E'
)
print(f"Crossing: {crossing:.2f}°, Incident: {incident:.2f}°")
```

### Example 2: MD Trajectory Analysis

**Old Code**:
```python
from immunex.analysis.angles import DockingAngleAnalyzer

analyzer = DockingAngleAnalyzer('md.tpr', 'md_pbc.xtc')
times, twist, tilt, swing = analyzer.calculate_docking_angles_trajectory(
    mhc_selection='chainID A and resid 50:86 140:176 and name CA',
    tcr_alpha_selection='chainID D',
    tcr_beta_selection='chainID E',
    tcr_v_selection='chainID D E and resid 1:115 and name CA',
    peptide_selection='chainID C and name CA',
    stride=10,
    output_file='docking_angles_old.csv'
)

import numpy as np
print(f"Twist: {np.mean(twist):.2f} ± {np.std(twist):.2f}°")
print(f"Tilt: {np.mean(tilt):.2f} ± {np.std(tilt):.2f}°")
print(f"Swing: {np.mean(swing):.2f} ± {np.std(swing):.2f}°")
```

**New Code**:
```python
from immunex.analysis.angles import DockingAnglePrimaryAnalyzer

analyzer = DockingAnglePrimaryAnalyzer('md.tpr', 'md_pbc.xtc')
times, crossing, incident = analyzer.calculate_docking_angles_trajectory(
    mhc_selection='chainID A',
    tcr_alpha_selection='chainID D',
    tcr_beta_selection='chainID E',
    stride=10,
    output_file='docking_angles_new.csv'
)

import numpy as np
print(f"Crossing: {np.mean(crossing):.2f} ± {np.std(crossing):.2f}°")
print(f"Incident: {np.mean(incident):.2f} ± {np.std(incident):.2f}°")
```

---

## CSV Output Format Changes

### Old Format (❌ Deleted)
```csv
Time(ps),Twist(deg),Tilt(deg),Swing(deg)
0.00,45.23,56.78,12.34
10.00,46.12,55.89,13.01
...
```

### New Format (✅ Implemented)
```csv
Time(ps),Crossing(deg),Incident(deg)
0.00,63.38,60.73
10.00,62.45,61.02
...
```

---

## Handling Non-Standard PDB Numbering

One of the key improvements is robust handling of non-standard PDB numbering.

### Old Approach (❌ Failed for Non-Standard Numbering)
```python
# Hardcoded residue ranges - breaks if PDB numbering starts from 1
mhc_selection = 'chainID A and resid 50:86 140:176 and name CA'
```

### New Approach (✅ Robust to Any Numbering)
```python
# Sequence alignment automatically maps α1/α2 regions
mhc_selection = 'chainID A'  # Just specify chain ID

# Behind the scenes:
# 1. Extract MHC sequence from PDB
# 2. Align to HLA-A*02:01 reference
# 3. Map reference positions (50-86, 140-176) to actual PDB resids
# 4. Use mapped residues for groove calculation
```

**Example**: If your PDB starts numbering from 1:
- Old code would look for resid 50-86 (wrong region)
- New code aligns sequence and finds correct α1/α2 regions (e.g., resid 26-62, 116-152)

---

## New Module Overview

### Core Modules (All NEW)

1. **`MHCSequenceAligner`** - Sequence alignment to HLA-A*02:01 reference
   - Automatic α1/α2 region detection
   - Robust to non-standard PDB numbering
   - Quality validation (alignment score, identity, gaps)

2. **`MHCGrooveDetector`** - Groove geometry calculation
   - Groove axis: COM(α2) - COM(α1)
   - Groove plane: SVD fit of α1/α2 backbone
   - Uses aligned residue positions

3. **`PlaneFitter`** - SVD-based plane fitting
   - Least-squares plane fit
   - RMS deviation validation
   - Point projection utilities

4. **`DockingAnglePrimaryAnalyzer`** - Main analyzer
   - Crossing angle: `angle(TCR_axis, groove_axis)`
   - Incident angle: `angle(TCR_axis, groove_normal)`
   - Single frame and trajectory support

5. **`ConservedCysteineDetector`** - ANARCI integration (no changes)
   - Automatic TCR disulfide detection
   - IMGT numbering

---

## Troubleshooting

### Error: "Alignment score too low"

**Problem**: PDB MHC sequence diverges significantly from HLA-A*02:01

**Solutions**:
1. Check if MHC is HLA Class I (HLA-A, HLA-B, HLA-C)
2. For HLA-B or HLA-C, alignment may have lower score but still work
3. Verify MHC chain ID is correct

### Error: "Too many gaps in α1/α2 region"

**Problem**: Structural regions missing in PDB

**Solutions**:
1. Check PDB completeness (missing residues)
2. Verify MHC chain selection includes full α1/α2 domains
3. Use structure fixing tools (e.g., `pdbfixer`)

### Warning: "Groove plane fit quality is poor"

**Note**: This is NORMAL for curved MHC grooves

**Typical RMS values**:
- 3-5 Å: Good (straight groove)
- 5-8 Å: Acceptable (typical curved groove)
- >10 Å: Poor (severe distortion)

---

## Expected Angle Ranges

Based on typical TCR-pMHC structures:

| Angle | Typical Range | Physical Meaning |
|-------|---------------|------------------|
| **Crossing** | 20-70° | Lateral orientation across groove |
| **Incident** | 40-80° | Tilt relative to binding surface |

**Note**: Values outside these ranges may indicate:
- Non-canonical docking modes
- Structural distortions during MD
- Unusual TCR-pMHC geometries

---

## Testing Your Migration

Run validation tests after migration:

```bash
# Test 1: MHC alignment
python development/test_mhc_alignment.py

# Test 2: Groove geometry
python development/test_groove_geometry.py

# Test 3: Angle calculation
python development/test_docking_angles_primary.py
```

All tests should pass (9/9 tests) before using in production.

---

## Additional Resources

- **Specification**: `docs/specs/docking_angles.md`
- **Implementation Plan**: `development/ANGLE_MODULE_REDESIGN.md` (if exists)
- **Example Scripts**: `examples/docking_angles_usage.py` (to be updated)
- **Batch Processing**: `scripts/batch_docking_angles.py` (to be updated)

---

## Support

If you encounter issues during migration:

1. Check that your PDB files contain complete MHC α1/α2 domains
2. Verify ANARCI is installed for TCR cysteine detection
3. Run validation tests to diagnose specific problems
4. Review alignment quality metrics in log output

---

## Conclusion

The new implementation provides:
- ✅ Robust MHC region detection (sequence-based)
- ✅ Simplified API (fewer required parameters)
- ✅ Clear geometric angle definitions
- ✅ Better handling of non-standard PDB files
- ✅ Improved scientific accuracy

While the API has changed significantly, the migration is straightforward:
just remove unnecessary parameters and rename angle variables.
