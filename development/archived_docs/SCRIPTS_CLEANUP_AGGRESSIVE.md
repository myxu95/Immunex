# Scripts Directory - Aggressive Cleanup Report

**Date**: 2026-03-16
**Action**: Remove obsolete scripts replaced by unified `imn` CLI

## Cleanup Script

Execute the following:

```bash
cd /home/xumy/work/development/Immunex
bash cleanup_scripts.sh
```

Or manually delete:

```bash
cd scripts/
rm -f imn_batch.py                    # Replaced by imn CLI
rm -f md_workflow.py                  # Wrapper (quality + slurm)
rm -f pbc_rmsd_workflow.py            # Single trajectory workflow
rm -f run_trajectory_analysis.py      # Old trajectory runner
rm -f run_batch_cdr_with_env.sh       # CDR batch wrapper
```

## Files to Delete (5 files)

| File | Reason | Replacement |
|------|--------|-------------|
| `imn_batch.py` | Prototype version | Use `imn` CLI |
| `md_workflow.py` | Wrapper script | Chain `md_quality_check.py` + `batch_pbc_slurm.py` |
| `pbc_rmsd_workflow.py` | Single trajectory workflow | Use `imn` subcommands |
| `run_trajectory_analysis.py` | Old trajectory runner | Use `imn tcr-rmsd`, `imn phla-rmsd`, etc. |
| `run_batch_cdr_with_env.sh` | CDR batch wrapper | Use `imn cdr-rmsf` |

## Expected Result

### Before Cleanup (12 files)
```
scripts/
├── imn                           # ✅ Keep
├── imn_batch.py                  # ❌ Delete
├── pbc_process.py                # ✅ Keep
├── md_quality_check.py           # ✅ Keep
├── batch_pbc_slurm.py            # ✅ Keep
├── batch_worker.py               # ✅ Keep
├── md_workflow.py                # ❌ Delete
├── pbc_rmsd_workflow.py          # ❌ Delete
├── run_trajectory_analysis.py    # ❌ Delete
├── md_task_organizer.py          # ✅ Keep
├── concatenate_multimodel_pdbs.py # ✅ Keep
├── run_batch_cdr_with_env.sh     # ❌ Delete
└── monitor_batch.sh              # ✅ Keep
```

### After Cleanup (8 files) ✨

```
scripts/
├── imn                           # Unified CLI (9 subcommands)
├── pbc_process.py                # Single trajectory PBC processing
├── md_quality_check.py           # MD quality check
├── batch_pbc_slurm.py            # SLURM script generator
├── batch_worker.py               # SLURM worker
├── md_task_organizer.py          # Task management utility
├── concatenate_multimodel_pdbs.py # PDB concatenation tool
└── monitor_batch.sh              # Batch monitoring script
```

## Benefits

- **File count**: 12 → 8 (33% reduction)
- **Clarity**: Each script has a clear, single purpose
- **No duplication**: All batch analysis goes through `imn`
- **Maintainability**: Fewer entry points, clearer responsibilities

## Core Workflow After Cleanup

### Batch Analysis (Use imn)
```bash
# All batch analysis through unified CLI
imn tcr-rmsd --input-dirs /data --output ./results
imn phla-rmsd --input-dirs /data --output ./results
imn docking-angles --input-dirs /data --output ./results
```

### Single Trajectory Processing
```bash
# PBC correction
python scripts/pbc_process.py -f md.xtc -s md.tpr -o processed/

# Quality check
python scripts/md_quality_check.py /data/simulations
```

### Cluster Deployment
```bash
# Generate SLURM scripts
python scripts/batch_pbc_slurm.py /data/simulations --batch-size 10

# Monitor jobs
bash scripts/monitor_batch.sh
```

### Utilities
```bash
# Organize incomplete tasks
python scripts/md_task_organizer.py /data/tasks --incomplete-dir /data/incomplete

# Concatenate multi-model PDBs
python scripts/concatenate_multimodel_pdbs.py input.pdb -o output.pdb
```

## Verification

After cleanup, verify the structure:

```bash
ls -lh scripts/
# Should show 7 .py files + 1 .sh file + 1 imn (no extension)
```

Test the unified CLI:

```bash
imn --help
# Should show 9 subcommands: tcr-rmsd, phla-rmsd, etc.
```

---

**Status**: Ready to execute
**Risk**: Low (all deleted scripts have replacements)
**Backup**: Git history preserves all files
