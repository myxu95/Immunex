# Immunex Product Scope

**Status**: Proposed engineering baseline
**Audience**: Core developers, research users, future contributors
**Goal**: Define what Immunex is, what it is not, and where engineering effort should concentrate.

## Product Definition

Immunex is a task-oriented MD analysis platform for TCR-pMHC systems.

It is not just a loose collection of analysis scripts. The product should present a stable workflow that takes a simulation task as input, runs standardized preprocessing and quality control, executes analysis pipelines, and emits reproducible artifacts.

## Primary User Problems

Immunex should solve these concrete problems:

1. Convert heterogeneous MD task folders into a consistent execution model.
2. Enforce a correct analysis order: preprocess first, then quality assessment, then downstream analysis.
3. Package domain-specific analyses behind stable commands and Python APIs.
4. Produce outputs that are easy to review, compare, archive, and batch-execute on HPC systems.

## Product Pillars

### 1. Task-Centric Execution
Everything starts from a `Task` rather than ad hoc file paths.

### 2. Standardized Workflows
The product should expose a small set of canonical workflows instead of many partially overlapping scripts.

### 3. Reproducible Artifacts
Every pipeline run should leave behind machine-readable outputs, human-readable reports, and execution metadata.

### 4. Layered Extensibility
New analysis methods should be added as modules and pipelines, not as standalone scripts that duplicate discovery and orchestration logic.

## In Scope

These capabilities are in scope for the core product:

- Raw trajectory preprocessing and PBC correction
- Processed trajectory quality assessment
- Core structural and trajectory analyses
- Task discovery and manifest generation
- Batch execution on workstation and HPC environments
- Standard output structure and report generation
- Python API and CLI for the same workflows

## Out of Scope

These should not be part of the core product unless a dedicated subsystem is designed:

- General-purpose MD workflow management outside Immunex task formats
- Arbitrary workflow builders for unrelated biomolecular systems
- Long-lived database services
- Web frontend work before the task/artifact contract is stable
- Experimental scripts that bypass product APIs and write arbitrary outputs

## Core Product Surfaces

Immunex should have exactly four top-level workflow surfaces:

1. `preprocess`
   Input: raw MD task
   Output: processed trajectory artifacts
2. `quality`
   Input: processed MD task
   Output: quality metrics and qualification reports
3. `analyze`
   Input: processed and qualified MD task
   Output: scientific analysis artifacts such as RMSD, RMSF, angles, contacts
4. `batch`
   Input: a directory or manifest of tasks
   Output: aggregated execution results and per-task artifacts

## Design Constraints

The following constraints should drive design decisions:

- CLI and Python API must map to the same workflow semantics.
- `scripts/` is not a product layer; it is an operations and migration layer.
- A pipeline should never depend on README instructions or manual folder conventions that are not represented in code.
- Every output that matters downstream must have a stable file name or schema.
- Adding a new analysis should not require creating a new batch engine.

## Product Maturity Model

### Stage 1: Structural Stability
- Stable task model
- Stable pipeline boundaries
- Stable artifact layout
- Minimal smoke coverage

### Stage 2: Workflow Completeness
- Single-task and batch workflows aligned
- Documentation consistent with implementation
- HPC execution path integrated with the same task contract

### Stage 3: Scientific Expansion
- More analysis modules
- Better comparative reporting
- Higher-level dashboards or web UX if needed

## Success Criteria

Immunex is becoming a coherent product when the following are true:

- A new contributor can identify the product entrypoints in under five minutes.
- A user can run a canonical workflow without touching `scripts/`.
- Every batch run leaves behind structured per-task results and a summary artifact.
- The team can add a new analysis pipeline without copy-pasting discovery, logging, or batch logic.
