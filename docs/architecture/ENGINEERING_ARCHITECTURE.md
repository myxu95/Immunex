# Immunex Engineering Architecture

**Status**: Proposed target architecture
**Goal**: Define stable responsibilities for each layer so development can proceed without continuous structural drift.

## Architectural Principle

Immunex should be designed around this rule:

**Domain logic, workflow orchestration, and delivery interfaces must remain separate.**

When these concerns are mixed, the system degrades into script sprawl.

## Canonical Layer Model

```text
Layer 5: Documentation and Configuration
Layer 4: Delivery Interfaces
Layer 3: Workflow Orchestration
Layer 2: Domain Modules
Layer 1: Execution Contract
```

## Layer 1: Execution Contract

Current home:
- `immunex/core`

Responsibility:
- Task models
- pipeline context
- manifest formats
- exceptions
- path and artifact identity rules

This layer defines the objects the rest of the system speaks.

Target core objects:
- `TaskSpec`
- `TaskInputFiles`
- `TaskArtifacts`
- `PipelineContext`
- `DiscoveryReport`
- `ExecutionSummary`

Rule:
No scientific algorithm belongs here.

## Layer 2: Domain Modules

Current home:
- `immunex/analysis`
- part of `immunex/utils`

Responsibility:
- PBC processing primitives
- RMSD and RMSF calculations
- angle analysis
- contact analysis
- quality metrics
- domain-specific chain and index logic

Rule:
These modules do the actual work. They should not parse CLI arguments, discover batch tasks, or decide output folder policy.

## Layer 3: Workflow Orchestration

Current home:
- `immunex/pipeline`

Responsibility:
- compose nodes into workflows
- validate prerequisites between steps
- define workflow inputs and outputs
- create structured reports from lower-level results

Canonical workflow types:
- `PreprocessPipeline`
- `QualityPipeline`
- `AnalysisPipeline`
- `CompositeWorkflowPipeline`

Rule:
Pipelines should call domain modules and write artifacts through the execution contract, not via manual path concatenation.

## Layer 4: Delivery Interfaces

Current home:
- `immunex/cli`
- `immunex/cluster`

Responsibility:
- CLI parsing
- user-facing help and command structure
- HPC wrappers that package the same workflows for remote execution

Rule:
This layer must stay thin. It should translate user intent into pipeline execution, not contain core product logic.

## Layer 5: Documentation and Configuration

Current home:
- `docs/`
- `templates/`
- `scripts/configs/`

Responsibility:
- user documentation
- operator documentation
- workflow config templates
- environment descriptions

Rule:
Documentation must describe actual product entrypoints, not legacy execution paths that bypass the architecture.

## Current Repository Gaps

The current repository still has several cross-layer leaks:

- Some workflow behavior still depends on `scripts/`.
- Some `utils/` modules are domain modules but are named as generic utilities.
- Some topology and chain utilities still live outside `analysis/topology`.
- CLI, docs, and historical script surfaces still need periodic drift checks.

## Target Package Boundaries

```text
immunex/
  core/        # contracts, models, manifest, context, exceptions
  domain/      # scientific and MD algorithms
  pipeline/    # workflow composition and execution
  cli/         # local user interface
  cluster/     # remote execution adapters
```

If renaming `analysis/` to `domain/` is too disruptive, keep `analysis/` but treat it semantically as the domain layer.

## Engineering Rules

1. A CLI command may call a pipeline, but not a low-level script.
2. A pipeline may call a domain module, but not parse raw CLI arguments.
3. A domain module may return results, but should not decide final report layout.
4. `scripts/` may assist migration or operations, but should not be required for core product workflows.
5. Any new feature must declare:
   - task inputs
   - produced artifacts
   - owning module
   - pipeline surface

## Recommended Primary Workflows

### `preprocess`
Uses raw task inputs and emits processed trajectory artifacts.

### `quality`
Uses processed task inputs and emits machine-readable and human-readable reports.

### `analyze`
Uses qualified processed tasks and emits one or more analysis artifact sets.

### `batch`
Wraps the same workflows over multiple tasks and emits a structured summary.

## Observability Requirements

Each workflow execution should record:

- input task identity
- selected parameters
- software version
- warnings and errors
- produced artifacts
- execution timestamps

These should live in a run-level metadata artifact so results are auditable.
