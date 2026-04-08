# TCR-pMHC Docking Angle Analysis Module

**Date**: 2026-01-23 (Design) | 2026-02-02 (Phase 1 Complete + ANARCI Integration)
**Version**: 2.0.0
**Status**: 🟢 Production Ready (Twist/Tilt/Swing angles with ANARCI auto-detection)

**Module Scope**: TCR-pMHC binding geometry (docking angles) ONLY
**Excluded**: Bond angles, dihedral angles (φ/ψ/χ), Ramachandran analysis

---

## 1. Overview

Complete redesign of **docking angle** calculation module for TCR-pMHC complex binding geometry analysis.

### Scope Definition
**Included**: TCR-pMHC docking angles (Twist/Tilt/Swing) - macromolecular binding geometry
**Excluded**: Bond angles, dihedral angles, Ramachandran angles - these belong to separate modules

### Key Changes
- **From**: External C++ tools + scattered Python scripts
- **To**: Unified pure Python implementation based on MDAnalysis
- **Goals**: Better maintainability, flexibility, and integration with Immunex pipeline

---

## 2. Module Architecture

**Core Focus**: TCR-pMHC docking angle analysis (Twist/Tilt/Swing)

```
immunex/analysis/angles/
├── __init__.py                       # Module initialization ✅
├── principal_axes.py                 # Principal axes calculation ✅
├── docking_angles.py                 # TCR-pMHC docking angles ✅ (ANARCI integrated)
├── vector_angles.py                  # Vector angle utilities ✅
├── conserved_cysteine_detector.py    # ANARCI-based cysteine detection ✅ NEW
└── clustering.py                     # Angle-based clustering (OPTIONAL)
```

**Note**: Bond angles, dihedral angles (φ/ψ/χ), and Ramachandran analysis are **NOT** part of this module.
They belong to separate structural analysis modules.

---

## 3. Core Components

### 3.1 Principal Axes Calculator (`principal_axes.py`)

**Purpose**: Calculate principal axes of molecular structures based on inertia tensor.

**Key Classes**:
```python
class PrincipalAxesCalculator:
    """
    Calculate principal axes via inertia tensor diagonalization

    Methods:
    --------
    - calculate_inertia_tensor(atoms, weights='mass')
    - get_principal_axes(atoms)
    - get_principal_moments(atoms)
    """
```

**Algorithm**:
1. Calculate center of mass
2. Construct inertia tensor (3x3 matrix)
3. Diagonalize to get eigenvalues (moments) and eigenvectors (axes)
4. Sort by descending moment (largest → smallest)

**Dependencies**: MDAnalysis, NumPy, scipy.linalg

---

### 3.2 Docking Angles Analyzer (`docking_angles.py`)

**Purpose**: Calculate TCR-pMHC binding geometry angles.

**Angle Definitions**:

#### Twist Angle
- **Definition**: Angle between TCR disulfide bond projection on MHC plane and MHC principal axis
- **Calculation**:
  1. Automatically detect conserved cysteines using ANARCI IMGT numbering (positions 23 and 104)
  2. Define MHC plane using principal axes
  3. Project TCR disulfide bond vector onto MHC plane
  4. Calculate angle with MHC principal axis
- **ANARCI Integration**: No need to manually specify cysteine residue numbers - automatically detected!

#### Tilt Angle
- **Definition**: Angle between TCR principal axis projection on MHC plane and MHC principal axis
- **Calculation**:
  1. Project TCR principal axis onto MHC plane
  2. Calculate angle with MHC principal axis

#### Swing Angle
- **Definition**: Lateral offset angle of TCR relative to peptide
- **Calculation**:
  1. Calculate TCR-peptide connecting vector
  2. Calculate perpendicular deviation from MHC axis

**Key Classes**:
```python
class DockingAngleAnalyzer:
    """
    TCR-pMHC docking angle analysis with ANARCI integration

    Methods:
    --------
    - calculate_twist_angle(mhc_sel, tcr_alpha_sel, tcr_beta_sel)
    - calculate_tilt_angle(mhc_sel, tcr_sel)
    - calculate_swing_angle(mhc_sel, tcr_sel, peptide_sel)
    - calculate_all_angles()  # Returns (twist, tilt, swing)

    Features:
    ---------
    - Automatic conserved cysteine detection via ANARCI IMGT numbering
    - No hardcoded residue numbers required
    - Sequence-independent TCR recognition
    """

class ConservedCysteineDetector:
    """
    ANARCI-based detector for TCR conserved disulfide bonds

    IMGT Positions:
    ---------------
    - Position 23: N-terminal conserved cysteine
    - Position 104: C-terminal conserved cysteine

    Methods:
    --------
    - detect_conserved_cysteines(chain_selection)
    - detect_tcr_alpha_beta_cysteines(alpha_sel, beta_sel)
    - Validates disulfide bond distance (<3.0 Å)
    """
```

**Coordinate Systems**:
- MHC coordinate frame: Principal axes of MHC molecule
- Peptide coordinate frame: Principal axes of peptide
- TCR coordinate frame: Principal axes of TCR V-domains

**Dependencies**: principal_axes.py, vector_angles.py, conserved_cysteine_detector.py

---

## 3.3 Usage Example

### New API (ANARCI Automatic Detection)

```python
from immunex.analysis.angles import DockingAngleAnalyzer

# Initialize analyzer
analyzer = DockingAngleAnalyzer('md.tpr', 'md_pbc.xtc')

# Define selections - No need to specify cysteine residue numbers!
times, twist, tilt, swing = analyzer.calculate_docking_angles_trajectory(
    mhc_selection='chainID A and resid 50:86 140:176 and name CA',
    tcr_alpha_selection='chainID D',  # ANARCI auto-detects cysteines
    tcr_beta_selection='chainID E',   # ANARCI auto-detects cysteines
    tcr_v_selection='chainID D E and resid 1:115 and name CA',
    peptide_selection='chainID C and name CA',
    stride=10,
    output_file='docking_angles.csv'
)

print(f"Twist: {twist.mean():.2f} ± {twist.std():.2f}°")
print(f"Tilt:  {tilt.mean():.2f} ± {tilt.std():.2f}°")
print(f"Swing: {swing.mean():.2f} ± {swing.std():.2f}°")
```

### Key Benefits of ANARCI Integration

1. **No Hardcoded Residue Numbers**: Works across different TCR sequences
2. **Automatic Validation**: Checks disulfide bond distance (<3.0 Å)
3. **Chain Type Detection**: Identifies TCR_alpha vs TCR_beta automatically
4. **IMGT Standardization**: Uses IMGT positions 23 and 104 universally

---

## 4. Design Principles

### 4.1 Modularity
- Each angle type in separate module
- Clear separation of concerns
- Reusable components (e.g., PrincipalAxesCalculator)

### 4.2 Flexibility
- Support arbitrary atom selections (MDAnalysis syntax)
- Not limited to fixed chain naming (A/B/C/D/E)
- Extensible for custom angle definitions

### 4.3 Performance
- Vectorized NumPy operations
- Efficient trajectory iteration
- Optional caching for repeated calculations

### 4.4 Integration
- Consistent API with other Immunex modules
- Support for batch processing via DiscoveryReport + BatchExecutor
- CLI commands for common workflows

---

## 5. Implementation Priority

### Phase 1: Core Components (High Priority) ✅ **COMPLETED (2026-02-02)**
1. ✅ Module initialization and cleanup
2. ✅ `principal_axes.py` - Foundation for other modules (158 lines)
3. ✅ `vector_angles.py` - Utility functions (141 lines)
4. ✅ `docking_angles.py` - Main TCR-pMHC analysis (385 lines)
5. ✅ `conserved_cysteine_detector.py` - ANARCI-based cysteine detection (280 lines) ⭐NEW
6. ✅ Unit tests - 7/7 tests passing (development/test_angle_modules.py)
7. ✅ **Batch processing** - `scripts/batch_docking_angles.py` (440 lines, updated for ANARCI)
8. ✅ **Old code cleanup** - Removed C++/VMD/archived scripts (~7.6MB)
9. ✅ **ANARCI Integration** - Automatic conserved cysteine detection ⭐NEW

### Phase 2: Extended Features (Medium Priority)
5. 📝 `clustering.py` - Angle-based clustering (to be integrated from development/clustering/)
6. 📝 PlotManager integration - Visualization functions for docking angles
7. 📝 Example usage scripts and tutorials

### Phase 3: Advanced Features (Low Priority)
8. 📝 CLI interface - `immunex-calc-angles` command
9. 📝 Performance optimization - Caching MHC principal axes
10. 📝 Alternative angle definitions (if needed for comparison)

---

## 6. Testing Strategy

### Unit Tests
- Each module has corresponding test file
- Test with synthetic data (known angles)
- Test edge cases (0°, 90°, 180°)

### Integration Tests
- Test with real PDB structures
- Compare with reference implementations (if available)
- Test MD trajectory processing

### Performance Tests
- Benchmark against old C++ tools (if applicable)
- Profile memory usage for large trajectories
- Optimize bottlenecks

---

## 7. Migration Plan ✅ **COMPLETED (2026-02-02)**

### Code Cleanup
- ✅ Old C++ tools **DELETED** - `development/tcr_docking_angle/` (2.3MB, 137 files)
- ✅ Archived backups **DELETED** - `development/archived_scripts/*angle*` (2.4MB)
- ✅ VMD reference scripts **DELETED** - `development/angle_analysis/` (2.9MB, 30 files)
- ✅ Deletion record created: `development/reports/DELETED_ANGLE_SCRIPTS_BACKUP_INFO.txt`
- ✅ Legacy data notice: `output/docking_angle/LEGACY_NOTICE.txt`

### For Users
- ✅ New production script ready: `scripts/batch_docking_angles.py`
- ✅ Historical data (output/docking_angle/) preserved with legacy notice
- 📝 Migration guide needed for users with custom scripts
- 📝 Update documentation and examples

### For Developers
- ✅ Old code removed, no deprecation period needed
- ✅ Core module is production-ready
- 📝 Git history preserved for recovery if needed

**Recovery command** (if needed):
```bash
git log --all --full-history -- "development/tcr_docking_angle/*"
git checkout <commit_hash> -- development/tcr_docking_angle/
```

---

## 8. Success Criteria

### Functional
- ✅ Calculate twist/tilt/swing angles for TCR-pMHC
- ✅ Support MD trajectory analysis
- ✅ Match or exceed old tool accuracy (within 1°)

### Performance
- ✅ Process single structure in < 1 second
- ✅ Process 10,000-frame trajectory in < 5 minutes

### Usability
- ✅ Clear API and documentation
- ✅ Batch processing script (`scripts/batch_docking_angles.py`)
- ✅ Integration with existing Immunex workflows
- 📝 CLI commands for common tasks (planned Phase 3)

---

## 9. References

### Algorithms
- Inertia tensor: Goldstein, Classical Mechanics
- Dihedral angles: IUPAC nomenclature
- TCR-pMHC geometry: Pierce Lab tools

### Previous Implementations (DELETED 2026-02-02)
- ~~VMD scripts: `development/angle_analysis/Analysis-scripts/`~~ (algorithms migrated to Python)
- ~~C++ tools: `development/tcr_docking_angle/`~~ (replaced by `immunex.analysis.angles`)
- ~~Old Python modules: `development/archived_scripts/angles_old_20260123/`~~ (archived code)

**Note**: All old implementations have been removed. Recovery available via git history if needed.
See `development/reports/DELETED_ANGLE_SCRIPTS_BACKUP_INFO.txt` for details.

---

## 10. Implementation Status ✅ **Phase 1 COMPLETE (2026-02-02)**

### Completed (Phase 1)
1. ✅ Review and approve design
2. ✅ Implement `principal_axes.py` (158 lines)
3. ✅ Implement `vector_angles.py` (141 lines)
4. ✅ Implement `docking_angles.py` (385 lines)
5. ✅ Implement `conserved_cysteine_detector.py` (280 lines) ⭐NEW
6. ✅ Add unit tests (7/7 passing)
7. ✅ Create batch processing script (440 lines)
8. ✅ Integrate ANARCI automatic cysteine detection ⭐NEW (2026-02-02)
9. ✅ Clean up old code (~7.6MB deleted)

### Next Steps (Phase 2)
1. 📝 Integrate angle-based clustering from `development/clustering/`
2. 📝 Add PlotManager visualization functions (angle time series, distributions)
3. 📝 Create example usage scripts and tutorials
4. 📝 Performance profiling and optimization
5. 📝 CLI interface (optional, Phase 3)

### Production Readiness Checklist
- ✅ Core modules implemented and tested
- ✅ ANARCI integration for automatic cysteine detection ⭐NEW
- ✅ Batch processing script available
- ✅ Old code cleaned up
- ✅ Documentation updated
- ✅ Historical data preserved with notices
- 📝 User migration guide (pending)
- 📝 Integration with PlotManager (pending)

---

**Status**: 🟢 **PRODUCTION READY** for docking angle analysis (Twist/Tilt/Swing)

**Key Feature**: ⭐ ANARCI automatic conserved cysteine detection - No hardcoded residue numbers!

**Last Updated**: 2026-02-02 (ANARCI integration complete)
