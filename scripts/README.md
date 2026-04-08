# Scripts Directory

This directory is an operational wrapper layer, not the core product layer.

The product surface of Immunex should be reached through the Python package and the `imn` CLI. Files in `scripts/` exist only for one of these reasons:

- HPC or shell-oriented wrappers around product workflows
- migration compatibility for older usage patterns
- one-off operational helpers that do not belong in the package API

## What Belongs Here

Current approved script categories:

### CLI Wrapper
- `imn`: thin launcher for the package CLI

### Operational Helpers
- `pbc_process.py`: single-trajectory operational wrapper around PBC processing
- `batch_pbc_slurm.py`: SLURM generation wrapper
- `batch_worker.py`: distributed worker wrapper
- `md_quality_check.py`: operational quality-check entrypoint
- `md_task_organizer.py`: task organization helper
- `monitor_batch.sh`: shell monitoring helper
- `concatenate_multimodel_pdbs.py`: standalone file utility
- `vmd_chain_contacts.tcl`: VMD integration helper

## What Does Not Belong Here

The following should not be introduced into `scripts/` going forward:

- core scientific logic
- new product workflows that bypass `immunex/cli`
- duplicate batch engines per analysis type
- modules that need to be imported as part of the stable Python API
- experimental scripts with hardcoded paths

## Ownership Rules

- If users should rely on it as part of the product, it belongs in `immunex/`.
- If it exists only to launch or adapt a product workflow, it may live in `scripts/`.
- If it is exploratory, put it in `development/`.

## Maintenance Rule

Any script added here must document:

- why it is not part of `immunex/`
- which package module it wraps or depends on
- whether it is temporary compatibility code or a long-term operational tool
