# Module Migration Inventory

**Status**: Execution inventory
**Goal**: Convert the ownership blueprint into a concrete module-by-module migration list.

## Priority Legend

- `P0`: do next, structural blocker
- `P1`: high value, should follow immediately after P0
- `P2`: worthwhile cleanup, not blocking

## P0 Moves

Completed:

- `immunex/core/pbc_processor.py` -> `immunex/analysis/trajectory/pbc.py`
- `immunex/core/batch_process.py` -> `immunex/pipeline/batch_workflow.py`
- `immunex/analysis/analysis_pipeline.py` -> `immunex/pipeline/analysis_workflows.py`
- `immunex/core/index_generation.py` -> `immunex/analysis/topology/index_generation.py` as the domain owner

### 1. `immunex/core/pbc_processor.py`

Status:

- completed

Current issue:

- domain-heavy preprocessing logic lives in `core`

Target:

- move to `immunex/analysis/trajectory/pbc.py`

Why:

- PBC correction is scientific trajectory processing, not a core contract

Follow-up:

- update `PreprocessNode`
- update any direct imports
- keep a temporary forwarding import only if absolutely necessary

### 2. `immunex/core/batch_process.py`

Status:

- completed

Current issue:

- workflow wrapper still lives in `core`

Target:

- move to `immunex/pipeline/batch_workflow.py`

Why:

- it is orchestration and compatibility wrapping over discovery + executor

Follow-up:

- re-export from `immunex.pipeline` only if you still want a Python convenience API
- keep `core` free of workflow surfaces

### 3. `immunex/analysis/analysis_pipeline.py`

Status:

- completed

Current issue:

- file name and role suggest orchestration inside `analysis`

Target:

- move to `immunex/pipeline/analysis_workflows.py` or delete if superseded

Why:

- “pipeline” is orchestration language, not domain language

### 4. `immunex/core/index_generation.py`

Status:

- partially completed

Current issue:

- mixed concerns: part domain semantics, part workflow helper

Current state:

- moved to `immunex/analysis/topology/index_generation.py`

Remaining target:

- extract any workflow-facing helpers into `immunex/pipeline/support/index_generation.py` only if they survive owner review

Why:

- domain contracts no longer belong in `core`
- any residual workflow logic should be isolated later under `pipeline/support`

## P1 Moves

### 5. `immunex/analysis/index_manager.py`

Status:

- completed as a topology-owned module

Current issue:

- mixed domain convenience and execution-facing behavior

Target:

- `immunex/analysis/topology/index_manager.py`

### 6. `immunex/utils/intelligent_chain_identifier.py`

Status:

- completed

Target:

- `immunex/analysis/topology/intelligent_chain_identifier.py`

Why:

- chain identity is domain interpretation

### 7. `immunex/utils/topology_chain_identifier.py`

Status:

- completed

Target:

- `immunex/analysis/topology/topology_chain_identifier.py`

### 8. `immunex/utils/tpr_chain_extractor.py`

Status:

- completed

Target:

- `immunex/analysis/topology/tpr_chain_extractor.py`

### 9. `immunex/utils/shortest_chain_detector.py`

Status:

- completed

Target:

- `immunex/analysis/topology/shortest_chain_detector.py`

### 10. `immunex/utils/strict_shortest_chain_detector.py`

Status:

- completed

Target:

- `immunex/analysis/topology/strict_shortest_chain_detector.py`

### 11. `immunex/utils/index_generator.py`

Status:

- completed

Target:

- `immunex/analysis/topology/component_index_generator.py`

### 12. `immunex/utils/chain_based_index_generator.py`

Status:

- completed

Target:

- `immunex/analysis/topology/chain_based_index_generator.py`

### 13. `immunex/utils/chain_identification_adapter.py`

Status:

- completed as a domain-facing module

Target:

- `immunex/analysis/topology/chain_identification_adapter.py`

Decision rule:

- if it transforms scientific chain identity models, keep in `analysis`
- if it adapts context or workflow execution, move to `pipeline`

### 14. `immunex/utils/cdr_manager.py`

Status:

- completed

Target:

- `immunex/analysis/topology/cdr_manager.py`

Why:

- CDR semantics are domain logic, not generic utility behavior

## P2 Moves

### 15. `immunex/utils/pdb_chain_standardizer.py`

Status:

- completed

Target:

- `immunex/analysis/structure/pdb_chain_standardizer.py`

### 16. `immunex/utils/pdb_sequence_extractor.py`

Status:

- completed

Target:

- `immunex/analysis/structure/pdb_sequence_extractor.py`

### 17. `immunex/utils/pdb_structure_fixer.py`

Status:

- completed

Target:

- `immunex/analysis/structure/pdb_structure_fixer.py`

### 18. `immunex/utils/pdb_distance_trimmer.py`

Status:

- completed

Target:

- `immunex/analysis/structure/pdb_distance_trimmer.py`

### 19. `immunex/utils/pdb_downloader.py`

Target:

- keep in `utils/` if treated as generic fetch helper
- move to `integrations/pdb/` if it grows network/provider logic

## Keep Where They Are

### Keep In `core`

- `context.py`
- `exceptions.py`
- `models.py`
- `task_discovery.py`
- `task_adapters.py`
- `base_node.py`

### Keep In `pipeline`

- `base_pipeline.py`
- `batch_executor.py`
- `standard_pipelines.py`
- `analysis_pipelines.py`
- `nodes/*`

### Keep In `cluster`

- `slurm_generator.py`

But rule:

- cluster code must only package pipeline-facing entrypoints

### Keep In `legacy`

- `legacy/task_discovery.py`
- `legacy/batch_processor.py`

## Migration Order

1. Move remaining topology utilities out of `utils/`
2. Review whether `analysis/topology/index_generation.py` still contains workflow-facing helpers worth extracting
3. Rebuild imports and examples for remaining topology ownership cleanup
4. Add package-level README notes after each move

## Verification Checklist Per Move

For each migrated module:

- update imports
- remove stale default exports
- run `py_compile`
- update active docs/examples
- add or update one smoke test

## Stop Conditions

Do not move modules blindly when:

- the file mixes two distinct responsibilities that should first be split
- the move would create a larger ambiguous module in the target location
- the new owner has not been defined clearly

Split first, move second.
