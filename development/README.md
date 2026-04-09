# Development Directory

This directory is the experimental surface of the repository. Nothing here is part of the stable product contract unless it is explicitly promoted into `immunex/`, `examples/`, or `docs/`.

## Purpose

Use `development/` for:

- exploratory analysis scripts
- diagnostics and debugging tools
- case studies with environment-specific assumptions
- temporary validation code during refactors
- archived historical material that remains useful for reference

## Subdirectory Semantics

### `archived_docs/`, `archived_scripts/`, `archived_utilities/`
Historical material kept for reference only.

### `archived_assets/`
Bundled historical assets, external snapshots, or one-off source packages kept only for traceability.

### `reports/`
Implementation summaries, migration reports, and internal writeups that are useful to developers but are not part of the public product documentation set.

### `diagnostics/`
One-off debugging and inspection scripts for investigating algorithm behavior, geometry issues, or unexpected outputs.

### `validation/`
Developer-only validation and environment-check scripts kept outside the product-facing `tests/` suite.

### `case_studies/`
Hardcoded or sample-specific investigations.

### `clustering/`, `plotting/`, `visualization/`, `utilities/`
Experimental tooling and research support code.

### `validation/angle_refactoring_tests_2026_02/`
Focused refactor-validation work for the docking-angle module, kept as developer-only validation rather than product test coverage.

### `workspaces/`
Large research workspaces and embedded external experiments.

## Promotion Rules

A file should move out of `development/` only when all of the following are true:

- no hardcoded local paths remain
- the interface is stable enough for repeated use
- ownership is clear
- the output contract is defined
- at least a smoke-level verification path exists

Promotion targets:

- move product code to `immunex/`
- move supported examples to `examples/`
- move stable documentation to `docs/`
- move long-term operations wrappers to `scripts/`

## Boundary Rule

Do not make product code depend on files in `development/`.
