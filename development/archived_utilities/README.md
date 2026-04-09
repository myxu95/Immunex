# Archived Utilities - 2026-03-12

This directory contains utilities that have been archived because their functionality is now available in the core Immunex modules.

## Archived Files

### Trim-related Scripts
**Reason**: Functionality now in `immunex/utils/pdb_distance_trimmer.py`

- `trim_pdb_by_distance.py` - Single PDB trimming
- `batch_trim_pdb_by_distance.py` - Batch PDB trimming
- `generate_trimming_summary.py` - Trimming summary generation
- `quick_trimming_summary.py` - Quick trimming summary
- `README_trim_pdb.py` - Documentation (should be .md)

### RMSD-related Scripts
**Reason**: Functionality now in `immunex/analysis/trajectory/rmsd*.py` and production scripts in `scripts/`

- `calculate_rmsd.py` - Single RMSD calculation
- `batch_rmsd_missing.py` - Batch RMSD for missing data

## Status

These scripts are preserved for reference but are no longer actively maintained. For current functionality, please use:

- **PDB trimming**: `immunex/utils/pdb_distance_trimmer.py`
- **RMSD calculation**: `immunex/analysis/trajectory/rmsd_analyzer.py`
- **Batch RMSD**: `scripts/batch_cdr_rmsd_exact.py`, `scripts/batch_tcr_rmsd.py`, etc.

## Archive Date
2026-03-12
