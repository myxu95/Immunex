# Immunex Package Layout

This package is the product core of the repository.

## Layer Map

- `core/`: contracts, models, context, task discovery, exceptions
- `analysis/`: scientific and MD-domain logic
- `pipeline/`: workflow orchestration and batch execution
- `cli/`: user-facing command entrypoints
- `cluster/`: HPC and remote-execution adapters
- `utils/`: compatibility helpers and generic support utilities only
- `legacy/`: retired interfaces kept only for explicit compatibility work

## Boundary Rule

If a module defines product behavior, users should be able to reach it through this package without going through `scripts/` or `development/`.

## Dependency Direction

Preferred dependency direction is:

`cli/` -> `pipeline/` -> `analysis/` and `core/`

`cluster/` -> `pipeline/` and `core/`

`analysis/` -> `core/` only when contracts or shared execution objects are needed

`utils/` should not grow into a second domain layer.

`legacy/` should never be treated as the default product surface.
