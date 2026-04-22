# Module Ownership Blueprint

**Status**: Target state
**Goal**: Define where code belongs so `core`, `analysis`, `pipeline`, and `utils` stop drifting.

## Decision Summary

Immunex should use this ownership model:

- `core`: contracts and execution primitives only
- `analysis`: scientific and MD-domain computation only
- `pipeline`: orchestration, workflow composition, and execution summaries
- `cli`: user-facing command translation only
- `cluster`: remote packaging of the same workflows
- `utils`: generic helper utilities only

## Layer Boundaries

### `core`

Owns:

- task discovery contract
- task and artifact models
- execution context objects
- exceptions
- minimal execution primitives shared by all workflows

Must not own:

- workflow wrappers
- scientific policy
- analysis-specific file naming conventions

Current files that fit:

- `immunex/core/context.py`
- `immunex/core/exceptions.py`
- `immunex/core/models.py`
- `immunex/core/task_discovery.py`
- `immunex/core/task_adapters.py`
- `immunex/core/base_node.py`

Current files that do not fit cleanly:

- none in the active `core/` surface

### `analysis`

Owns:

- RMSD, RMSF, RDF, Rg, hydrogen bonds
- angle and contact analysis
- quality metrics and validators
- scientifically meaningful topology, chain, and index interpretation
- preprocessing algorithms when they are domain logic rather than workflow glue

Must not own:

- CLI parsing
- batch task discovery
- job partitioning
- pipeline result summaries

Sub-groups that should exist conceptually:

- `analysis/trajectory`
- `analysis/structure`
- `analysis/quality`
- `analysis/angles`
- `analysis/support` or `analysis/topology`

### `pipeline`

Owns:

- workflow nodes
- workflow composition
- batch execution
- run-level summaries
- stage-specific input requirements

Must not own:

- low-level scientific implementation details
- ad hoc CLI behavior

Current files that fit:

- `immunex/pipeline/base_pipeline.py`
- `immunex/pipeline/batch_executor.py`
- `immunex/pipeline/standard_pipelines.py`
- `immunex/pipeline/analysis_pipelines.py`
- `immunex/pipeline/nodes/*`

Current files that should probably move here:

- future workflow-facing helpers extracted from topology/index generation

### `utils`

Owns only:

- plotting helpers
- generic path and cleanup helpers
- generic file/string helpers
- non-domain helper code with no workflow ownership

Must not own:

- chain identification policy
- topology interpretation
- index semantics tied to pHLA-TCR biology
- batch product entrypoints

## Final Ownership Decisions

### Keep In `core`

- `context.py`
- `exceptions.py`
- `models.py`
- `task_discovery.py`
- `task_adapters.py`
- `base_node.py`

### Move Out Of `core`

- `pbc_processor.py` -> `analysis/trajectory/pbc.py`
  Status: completed.
  Reason: this is domain processing, not a contract.

- `batch_process.py` -> `pipeline/batch_workflow.py`
  Status: completed.
  Reason: this is orchestration and compatibility wrapping, not core execution state.

- `index_generation.py` -> split:
  Status: partially completed.
  - workflow-facing pieces -> `pipeline/support/index_generation.py`
  - domain-facing pieces -> `analysis/topology/index_generation.py`

### Keep In `analysis`

- `analysis/trajectory/*`
- `analysis/structure/*`
- `analysis/quality/*`
- `analysis/angles/*`
- `analysis/allostery/*`
- `analysis/free_energy/*`

### Move Out Of `analysis`

- `analysis/analysis_pipeline.py` -> `pipeline/analysis_workflows.py`
  Status: completed.
  Reason: the name and responsibility are orchestration-shaped, not domain-shaped.

- `analysis/index_manager.py` -> `analysis/topology/index_manager.py`
  Status: completed.
  Reason: index ownership is domain-facing and belongs with topology semantics.

### Move Out Of `utils`

These modules are not generic helpers and should not remain in `utils/` long term.
The current high-priority domain-heavy modules have already been moved out.

Target home:

- chain and topology meaning -> `analysis/topology/`
- workflow helper wrappers -> `pipeline/support/`

### Keep In `utils`

- `plotting.py`
- `path_manager.py`
- `cleanup_manager.py`
- `group_selector.py`
- file conversion helpers that are truly generic

## Naming Corrections Still Needed

These names currently preserve history more than truth:

- `utils/group_selector.py`

If names remain misleading, contributors will continue putting code in the wrong layer.

## Non-Negotiable Rules For New Code

1. If a module owns a scientific calculation, it does not belong in `core` or `cli`.
2. If a module builds a workflow from multiple steps, it belongs in `pipeline`.
3. If a module exists only for backward compatibility, it should stay out of the active package tree and be archived explicitly.
4. If a module needs pHLA-TCR biological semantics, it is not a generic `utils` helper.
5. `utils` is for support code, not product ownership.

## Definition Of A Clean Structure

The structure is clean when:

- a contributor can choose the owning package without asking
- `core` contains no domain-heavy algorithms
- `analysis` contains no workflow orchestration
- `pipeline` contains no low-level science implementation
- `utils` contains no domain policy
- retired compatibility code stays outside the active package tree
