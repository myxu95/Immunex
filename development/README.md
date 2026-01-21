# Development Directory

This directory contains experimental code, temporary scripts, and case studies created during development.

## Directory Structure

### archived_scripts/
Deprecated scripts preserved for reference.
- `rmsd_old/` - Legacy RMSD implementations
- `rmsf_analysis_pipeline_phla.py` - Old RMSF analysis pipeline (merged into batch_whole_protein_rmsf.py)

### case_studies/
Scripts for specific case studies with hardcoded paths.
- `7n1e_chain_contacts_analysis.py` - Chain contact analysis for 7N1E complex

### workspaces/
Large workspace directories containing experimental data and intermediate results.

**Important Directories**:
- `FEL_workspace/` (23GB) - Free Energy Landscape analysis reference implementation
  - Purpose: Reference for developing aftermd.analysis.free_energy module
  - Contains complete FEL calculation workflow and validation data
  - See `FEL_README.md` for details

- `allosteric_workspace/` - Allosteric pathway analysis experiments
- `pattern_analysis/` - Pattern recognition experiments
- `pymol_visualization/` - PyMOL visualization scripts
- `cdr3_analysis_toolkit/` - CDR3 specific analysis tools

### logs/
Log file archives organized by date.
- `archive_2026_01/` - January 2026 log files

### Other Directories
- `clustering/`, `plotting/`, `utilities/`, `visualization/` - Various experimental scripts and tools
- `tcr_docking_angle/` - TCR docking angle analysis experiments

## Usage Notes

1. **Stability Not Guaranteed** - Code here may change or be removed without notice
2. **Hardcoded Paths** - case_studies scripts contain environment-specific paths
3. **Regular Cleanup** - Unused scripts should be moved to archived_scripts/
4. **Large Files** - workspaces/ directory is excluded from version control

## Migrating from Development to Production

To move a script from development to `scripts/` production:
- [ ] Remove hardcoded paths
- [ ] Add command-line argument support
- [ ] Follow naming conventions (batch_*.py)
- [ ] Add complete docstrings
- [ ] Pass code review

## Workspace Management

**FEL_workspace Importance**:
- Contains validated reference implementation for FEL analysis
- Used as blueprint for software module development
- Do not delete without extracting key algorithms

**General Workspaces**:
- Archive completed case studies to external storage
- Clean up temporary analysis results periodically
- Document important findings before archiving
