# Immunex Task And Artifact Specification

**Status**: Proposed stable contract
**Goal**: Standardize how Immunex discovers work and how downstream tooling consumes outputs.

## Core Entities

### Task
A `Task` is the smallest executable unit in Immunex.

A task contains:
- an identity
- input files
- optional metadata
- artifacts produced by workflows

### Artifact
An `Artifact` is any file or directory produced or consumed as part of a workflow contract.

Artifacts must be:
- typed
- named predictably
- attributable to a task and workflow step

## Canonical Task Directory

```text
task_root/
  task.yaml
  input/
    structure.pdb | md.gro
    md.tpr
    md.xtc
  artifacts/
    preprocess/
    quality/
    analysis/
  reports/
  logs/
```

For backward compatibility, the discovery layer may continue to support legacy layouts. But pipelines should write to the canonical structure.

## `task.yaml` Minimum Schema

```yaml
task_id: 1ao7_standard
metadata:
  system: 1ao7
  variant: standard
  source: experiment_a
input_files:
  structure: input/md.gro
  topology: input/md.tpr
  trajectory: input/md.xtc
tags:
  - tcr-pmhc
  - production-run
```

## Workflow States

A task should conceptually move through these states:

1. `raw`
2. `preprocessed`
3. `quality_checked`
4. `analyzed`
5. `batched`

These states do not need to be stored as a database field, but the filesystem and metadata should make them inferable.

## Required Artifacts By Workflow

### Preprocess Workflow
Input:
- raw trajectory
- topology
- structure if needed

Outputs:
- `artifacts/preprocess/trajectory_processed.xtc`
- `artifacts/preprocess/reference_topology.tpr`
- `artifacts/preprocess/reference_structure.pdb` or `.gro`
- `artifacts/preprocess/run_metadata.json`

### Quality Workflow
Input:
- processed trajectory
- topology

Outputs:
- `artifacts/quality/rmsd.csv`
- `artifacts/quality/quality_report.json`
- `reports/quality_report.md`
- `artifacts/quality/run_metadata.json`

### Analysis Workflow
Input:
- processed trajectory
- topology
- optional derived indices and chain mappings

Outputs:
- `artifacts/analysis/{analysis_name}/result.json`
- `artifacts/analysis/{analysis_name}/metrics.csv`
- `reports/{analysis_name}.md` if human-readable reporting is needed
- `artifacts/analysis/{analysis_name}/run_metadata.json`

### Batch Workflow
Input:
- a directory of tasks or a manifest

Outputs:
- `batch_summary.json`
- `batch_summary.csv` optional
- per-task workflow artifacts unchanged
- `logs/batch_execution.log`

## Artifact Metadata Requirements

Every `run_metadata.json` should include at least:

```json
{
  "task_id": "1ao7_standard",
  "workflow": "quality",
  "started_at": "2026-03-22T10:00:00",
  "finished_at": "2026-03-22T10:05:00",
  "parameters": {},
  "inputs": {},
  "outputs": {},
  "warnings": [],
  "errors": [],
  "immunex_version": "0.1.0"
}
```

## Discovery Rules

The discovery layer should support two modes.

### Mode 1: Explicit Task Mode
Use `task.yaml` when present.

### Mode 2: Legacy Discovery Mode
Infer files from historical layouts such as:
- files in task root
- files under `input/`
- files under `prod/` for raw trajectories

Legacy discovery is a compatibility feature. New pipelines should still emit canonical artifacts.

## Naming Rules

Use semantic names, not workflow-internal temporary names, for durable outputs.

Preferred:
- `trajectory_processed.xtc`
- `quality_report.json`
- `run_metadata.json`

Avoid as durable contract names:
- `md_pbc.xtc`
- `output_final_v2.csv`
- `tmp_quality.json`

Legacy filenames can still exist, but the canonical artifact name should be stable.

## Summary Contract

A downstream consumer should be able to answer these questions for any task without reading source code:

- What were the declared inputs?
- Which workflow ran?
- Which artifacts were produced?
- Did the task pass quality control?
- Which parameters were used?
- Where is the machine-readable summary?
