# Immunex Repository Restructuring Plan

**Status**: Proposed execution plan
**Goal**: Map the target architecture onto the current repository without destabilizing scientific functionality.

## Current Reality

The repository already contains the seeds of a better architecture:

- `immunex/core` for contracts and processing foundations
- `immunex/pipeline` for orchestration
- `immunex/analysis` for domain logic
- `immunex/cli` for entrypoints
- `immunex/cluster` for HPC execution

But the codebase still shows transitional symptoms:

- legacy and new discovery models coexist
- product workflows still leak into `scripts/`
- `development/` holds both notes and executable prototypes
- some `utils/` modules are actually domain modules
- public docs are ahead of implementation in some areas

## Restructuring Objective

Move the repository from a migration state to a stable product state without breaking useful scientific modules.

## Proposed Target Layout

```text
immunex/
  core/
  analysis/          # keep name if renaming is too costly
  pipeline/
  cli/
  cluster/
  integrations/      # optional future home for ANARCI, MoSAIC, external tools

docs/
examples/
tests/
scripts/             # non-product operational wrappers only
development/         # design notes, experiments, diagnostics, archived prototypes
```

## Package-Level Recommendations

### `immunex/core`
Keep and strengthen.

Should own:
- task models
- manifest and discovery contract
- pipeline context
- exceptions
- run metadata schemas

Should not own:
- domain-specific analysis logic
- CLI-like command behavior

### `immunex/analysis`
Keep as the domain layer.

Move here or keep here:
- RMSD and RMSF logic
- quality metrics
- angle logic
- contact logic
- chain and index logic that is scientifically meaningful

### `immunex/utils`
Shrink aggressively.

Current issue:
`utils/` mixes true utilities with domain-heavy modules such as chain identification and index generation.

Target rule:
- keep only generic helpers in `utils/`
- move scientific modules to `analysis/` or a future `domain/` package

### `immunex/pipeline`
Strengthen as the only orchestration layer.

Should own:
- standard workflows
- node composition
- batch execution behavior tied to pipeline semantics

### `immunex/cli`
Keep thin.

Should own:
- argument parsing
- help text
- user-facing summaries

Should not own:
- scientific logic
- manual script dispatching

### `scripts/`
Reduce to a compatibility and operations zone.

Allowed uses:
- temporary migration wrappers
- HPC environment bootstrap helpers
- maintenance scripts not part of the product API

Disallowed long-term use:
- core workflows that the CLI depends on
- duplicate batch engines

### `development/`
Split mentally into three buckets:
- design notes
- diagnostics and experiments
- archival material

If retained, add a small README that labels which files are historical, which are active diagnostics, and which are not product code.

## Execution Phases

### Phase 1: Contract Stabilization
Deliverables:
- stable product scope
- stable task and artifact spec
- stable top-level CLI workflow map

This phase is documentation-first and unblocks all later refactors.

### Phase 2: Product Surface Cleanup
Deliverables:
- CLI and README reflect actual implemented workflows
- no product workflow requires `scripts/`
- old entrypoints either removed or redirected cleanly

### Phase 3: Module Ownership Cleanup
Deliverables:
- audit `immunex/utils`
- move domain-heavy modules into clearer owners
- retire dead compatibility paths

### Phase 4: Output Standardization
Deliverables:
- canonical artifact names
- run metadata emitted by each workflow
- batch summary schema

### Phase 5: Verification Layer
Deliverables:
- smoke tests for CLI surfaces
- smoke tests for core task discovery
- minimal fixture-backed pipeline tests

## Immediate Refactor Candidates

These are the highest-value near-term moves.

1. Consolidate task discovery around one public contract.
2. Make `batch` a thin wrapper over real pipelines rather than separate script logic.
3. Decide which `utils/` modules are actually domain modules and move them deliberately.
4. Add `run_metadata.json` generation to core workflows.
5. Label or archive any development script that is no longer part of the product path.

## Decision Log To Maintain

During restructuring, keep a short decision log for each important move:

- why this module belongs where it does
- what compatibility path remains
- when the old path can be removed

Without that discipline, the repository will drift back into mixed ownership.

## Definition Of Done For The Restructure

The restructure is succeeding when:

- contributors can tell product code from experimental code immediately
- the CLI maps cleanly to pipelines
- task discovery and artifact layout are documented and enforced
- batch execution does not require duplicate logic per analysis
- README and docs describe the real product, not an aspirational hybrid state
