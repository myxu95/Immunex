# Scripts Directory

This directory is a thin tooling layer, not the core product layer.

The product surface of Immunex should be reached through the Python package and the `imn` CLI. Files in `scripts/` exist only for one of these reasons:

- report-build helpers that stay outside the package API
- standalone data-preparation helpers with narrow operational scope
- visualization utilities that are intentionally file-oriented

## What Belongs Here

Current approved script categories:

### Active Helpers
- `build_interaction_demo_html.py`: single-system interaction report builder
- `build_comparison_report_html.py`: two-condition comparison report builder
- `build_hla_reference_library.py`: HLA reference-library preparation utility
- `concatenate_multimodel_pdbs.py`: standalone file utility

## What Does Not Belong Here

The following should not be introduced into `scripts/` going forward:

- core scientific logic
- new product workflows that bypass `immunex/cli`
- duplicate CLI or batch entrypoints for workflows already exposed through `imn`
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
