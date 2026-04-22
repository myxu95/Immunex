# Task Discovery and Manifest Generation Guide

## Overview

The Immunex task discovery module provides a structured way to:
1. Scan batch root directories for single-task directories
2. Validate task completeness (structure, topology, trajectory files)
3. Generate structured manifests (JSONL or CSV format)
4. Classify tasks as valid, invalid, or ambiguous

## Design Principles

- **Single is a special case of batch** - Batch operations are built on single-task foundations
- **Standard directory enables discoverability** - Conventional directory layout allows automatic discovery
- **Manifest enables executability** - Structured manifest files drive batch execution
- **Scanning discovers candidates** - Discovery identifies potential tasks, manifest drives execution

## Quick Start

### Basic Usage

```python
from pathlib import Path
from immunex.core import discover_tasks, discover_tasks_from_list, write_manifest

# Discover tasks
report = discover_tasks(Path("/data/batch"))

# Print summary
report.print_summary()

# Write manifests
write_manifest(report, Path("manifest.jsonl"), format="jsonl")
write_manifest(report, Path("manifest.csv"), format="csv")

# Access results
print(f"Found {report.num_valid} valid tasks")
for task in report.valid_tasks:
    print(f"  - {task.task_id}: {task.input_files.structure_path}")
```


### Adapting To Legacy Pipelines

```python
from pathlib import Path
from immunex.core import discover_tasks, discovery_report_to_contexts

report = discover_tasks(Path("/data/batch"))
contexts = discovery_report_to_contexts(report)

print(f"Converted {len(contexts)} valid tasks to PipelineContext objects")
```

### Running the Example

```bash
python examples/task_discovery_example.py
```

## Directory Structure Conventions

### Standard Input Directory

**Recommended structure** (files in `input/` subdirectory):

```
batch_root/
  system_001/
    input/
      structure.pdb
      topology.tpr
      trajectory.xtc
    meta/
      metadata.json
```

### Root Directory Files

**Alternative structure** (files in task root):

```
batch_root/
  system_002/
    structure.pdb
    topology.tpr
    trajectory.xtc
```

### Explicit Configuration

**With task.yaml** (explicit file paths):

```
batch_root/
  system_003/
    task.yaml           # Explicit configuration
    data/
      my_structure.pdb
      my_topology.tpr
      my_trajectory.xtc
```

**task.yaml example:**

```yaml
task_id: custom_task_id  # Optional: override directory name
input_files:
  structure: data/my_structure.pdb
  topology: data/my_topology.tpr
  trajectory: data/my_trajectory.xtc
metadata:
  system_type: protein-ligand
  temperature: 300
  pressure: 1.0
tags:
  - production
  - validated
```

## Supported File Extensions

### Structure Files
- `.pdb` - Protein Data Bank format
- `.gro` - GROMACS structure format

### Topology Files
- `.tpr` - GROMACS portable run input (preferred)
- `.top` - GROMACS topology

### Trajectory Files
- `.xtc` - GROMACS compressed trajectory (preferred)
- `.trr` - GROMACS trajectory

## Task Validation

### Valid Task
All three core files (structure, topology, trajectory) are present and unambiguous.

```
Status: valid
Validation messages: []
```

### Invalid Task
One or more core files are missing.

```
Status: invalid
Validation messages: ["Missing trajectory file"]
```

### Ambiguous Task
Multiple candidates exist for a file type, preventing unique identification.

```
Status: ambiguous
Validation messages: ["Ambiguous trajectory: found 2 candidates"]
```

## Discovery Report Structure

### DiscoveryReport

```python
@dataclass
class DiscoveryReport:
    source_root: str           # Batch root directory
    discovered_at: str         # ISO 8601 timestamp
    total_tasks: int           # Total tasks discovered
    valid_tasks: List[DiscoveredTask]
    invalid_tasks: List[DiscoveredTask]
    ambiguous_tasks: List[DiscoveredTask]
    all_tasks: List[DiscoveredTask]

    # Computed properties
    num_valid: int
    num_invalid: int
    num_ambiguous: int
    success_rate: float
```

### DiscoveredTask

```python
@dataclass
class DiscoveredTask:
    task_id: str               # Unique identifier
    task_root: str             # Absolute path to task directory
    source_root: str           # Absolute path to batch root
    input_files: TaskInputFiles
    task_file_path: Optional[str]  # Path to task.yaml if present
    validation_status: str     # "valid", "invalid", "ambiguous"
    validation_messages: List[str]
    tags: List[str]
    metadata: Dict[str, Any]
    discovered_at: str         # ISO 8601 timestamp
```

### TaskInputFiles

```python
@dataclass
class TaskInputFiles:
    structure_path: Optional[str]
    topology_path: Optional[str]
    trajectory_path: Optional[str]
```

## Manifest Formats

### JSONL Format

One task per line, JSON-encoded:

```jsonl
{"task_id": "system_001", "task_root": "/data/batch/system_001", ...}
{"task_id": "system_002", "task_root": "/data/batch/system_002", ...}
```

**Advantages:**
- Streamable (can process line-by-line)
- Preserves full data structure
- Easy to parse with standard JSON tools

**Usage:**
```python
write_manifest(report, Path("manifest.jsonl"), format="jsonl")
```

### CSV Format

Tabular format with JSON-encoded complex fields:

```csv
task_id,task_root,structure_path,validation_status,metadata
system_001,/data/batch/system_001,/path/to/struct.pdb,valid,"{""key"": ""value""}"
```

**Advantages:**
- Excel/spreadsheet compatible
- Easy to filter and sort
- Human-readable

**Usage:**
```python
write_manifest(report, Path("manifest.csv"), format="csv")
```

## Error Handling

### Exception Hierarchy

```
ImmunexError
├── TaskDiscoveryError        # Fatal discovery errors (e.g., root not found)
├── TaskValidationError       # Internal validation errors
├── ManifestWriteError        # Manifest write failures
└── TaskFileParseError        # task.yaml parsing errors
```

### Example Error Handling

```python
from immunex.core import (
    discover_tasks,
    TaskDiscoveryError,
    ManifestWriteError
)

try:
    report = discover_tasks(Path("/data/batch"))
except TaskDiscoveryError as e:
    print(f"Discovery failed: {e}")
    sys.exit(1)

try:
    write_manifest(report, Path("manifest.jsonl"))
except ManifestWriteError as e:
    print(f"Manifest write failed: {e}")
    sys.exit(1)
```

## Advanced Usage

### Filtering Tasks

```python
# Get only valid tasks with specific tags
valid_production_tasks = [
    task for task in report.valid_tasks
    if "production" in task.tags
]

# Get tasks with specific metadata
high_temp_tasks = [
    task for task in report.valid_tasks
    if task.metadata.get("temperature", 0) > 350
]
```

### Custom Processing

```python
# Process only valid tasks
for task in report.valid_tasks:
    print(f"Processing {task.task_id}...")
    # Run analysis on task.input_files.trajectory_path

# Report on invalid tasks
for task in report.invalid_tasks:
    print(f"Skipping {task.task_id}:")
    for msg in task.validation_messages:
        print(f"  - {msg}")
```

### Programmatic Manifest Reading

```python
import json

# Read JSONL manifest
tasks = []
with open("manifest.jsonl", 'r') as f:
    for line in f:
        task_data = json.loads(line)
        tasks.append(task_data)

print(f"Loaded {len(tasks)} tasks from manifest")
```

## Integration with Batch Execution

### Recommended Workflow

1. **Discovery Phase** - Scan batch root, generate manifest
   ```python
   report = discover_tasks(Path("/data/batch"))
   write_manifest(report, Path("manifest.jsonl"))
   ```

2. **Review Phase** - Inspect manifest, fix issues
   ```bash
   # Review CSV in spreadsheet
   libreoffice manifest.csv

   # Or inspect programmatically
   python -c "from immunex.core import DiscoveryReport; ..."
   ```

3. **Execution Phase** - Consume manifest, run pipeline
   ```python
   # Load manifest
   tasks = load_tasks_from_manifest("manifest.jsonl")

   # Filter to valid tasks only
   valid_tasks = [t for t in tasks if t['validation_status'] == 'valid']

   # Execute pipeline
   for task in valid_tasks:
       run_pipeline(task)
   ```

## Best Practices

1. **Always review discovery reports** - Check validation messages before executing
2. **Use task.yaml for complex setups** - Explicit configuration prevents ambiguity
3. **Version your manifests** - Keep manifest alongside results for reproducibility
4. **Filter ambiguous tasks** - Resolve ambiguities before batch execution
5. **Log discovery results** - Save reports for debugging and auditing

## Troubleshooting

### No tasks discovered

**Symptom:** `report.total_tasks == 0`

**Possible causes:**
- Incorrect batch root path
- Non-standard directory structure
- Required files have unexpected names

**Solution:** Check directory structure matches conventions, or use task.yaml

### All tasks marked ambiguous

**Symptom:** `report.num_ambiguous == report.total_tasks`

**Possible causes:**
- Multiple trajectory files in each directory
- Inconsistent naming conventions

**Solution:** Remove duplicate files or use task.yaml to specify exact paths

### Task.yaml not recognized

**Symptom:** Task ignored despite having task.yaml

**Possible causes:**
- YAML syntax errors
- Invalid file permissions

**Solution:** Validate YAML syntax, check file is readable

## API Reference

See module docstrings for detailed API documentation:

```python
from immunex.core import discover_tasks, discover_tasks_from_list, write_manifest
help(discover_tasks)
help(write_manifest)
```

## Examples

Complete examples are available in `examples/task_discovery_example.py`.


## Multi-Level Discovery

Use `task_depth` when tasks are nested more than one directory below the batch root.

```python
from pathlib import Path
from immunex.core import discover_tasks

report = discover_tasks(Path("/data/nested_batch"), task_depth=2)
```

## Stage-Specific Discovery Requirements

Some pipeline stages do not need the full structure/topology/trajectory triad. For example, preprocessing and quality assessment can validate tasks with only topology and trajectory present.

```python
from pathlib import Path
from immunex.core import discover_tasks

report = discover_tasks(
    Path("/data/raw_md"),
    required_files=["topology", "trajectory"],
)
```

## Building A DiscoveryReport From A Manual Task List

```python
from immunex.core import discover_tasks_from_list

report = discover_tasks_from_list([
    {
        "task_id": "1ao7_standard",
        "structure": "data/1ao7/standard/md.gro",
        "topology": "data/1ao7/standard/md.tpr",
        "trajectory_raw": "data/1ao7/standard/md.xtc"
    }
])
```
