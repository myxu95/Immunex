# Archived Angle Calculation Modules

**Archive Date**: 2026-01-23

## Archived Contents

### 1. angles_old_20260123/
Original angle analysis modules from `immunex/analysis/angles/`:
- `angle_clustering.py` - Static angle clustering (K-means/hierarchical)
- `__init__.py` - Module initialization

**Reason for archiving**: Complete redesign of angle calculation architecture

### 2. tcr_docking_angle_backup_20260123/
External TCR docking angle calculation tools:
- C++ executable for docking/incident angle calculation
- Python wrappers (`tcr_angle_wrapper.py`, `batch_tcr_angles.py`)
- Requires GSL library and AHo numbering

**Reason for archiving**: Moving to pure Python implementation based on MDAnalysis

### 3. clustering_backup_20260123/
Static angle clustering analysis scripts:
- `step1_cluster_static_angles.py` - Clustering workflow script

**Reason for archiving**: Will be integrated into new angle analysis module

## New Implementation Plan

The new angle calculation module will include:
- Pure Python implementation (MDAnalysis + NumPy)
- Principal axes calculation (inertia tensor)
- Docking angles (Twist/Tilt/Swing for TCR-pMHC)
- Dihedral angles (backbone and sidechain)
- Vector angle utilities
- MD trajectory support
- Unified CLI interface

## Recovery

If you need to recover old functionality:
```bash
# View archived code
ls development/archived_scripts/angles_old_20260123/

# Copy back if needed
cp development/archived_scripts/angles_old_20260123/angle_clustering.py immunex/analysis/angles/
```

## Notes

- Old C++ tools remain in `development/tcr_docking_angle/` (not deleted)
- Archived copies are backups only
- New module design starts from clean slate
