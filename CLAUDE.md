# CLAUDE.md

This file provides guidance to Claude Code when working with the Immunex codebase.

## Project Overview

**Immunex** is a task-oriented analysis platform for TCR-pMHC molecular dynamics trajectories.

- **Language**: Python 3.9+
- **Main Dependencies**: MDAnalysis, GROMACS, NumPy, pandas
- **Entry Point**: `imn` CLI command (not scattered scripts)
- **Development Path**: `/home/xumy/work/development/Immunex`

## Communication Style

1. **Language**: Use Chinese for conversation, English for code/comments
2. **No emojis** in code or scripts unless explicitly requested
3. **Concise responses**: Focus on actionable information, avoid verbose summaries
4. **No temporary markdown reports**: Report results directly in conversation

## File Organization Rules

### Core Principles

1. **Test scripts go to `development/`**, not `scripts/` or project root
2. **No version suffixes**: Avoid `xxx_v2.py`, `xxx_new.py`, `xxx_exact.py`
3. **Consolidate duplicates**: Use parameters to control behavior, not multiple files
4. **Clean up after testing**: Remove temporary scripts promptly
5. **No temporary markdown files**: Don't create `UPDATE_2026-xx-xx.md`, `SUMMARY_xxx.md` unless explicitly requested

### Directory Structure

```
immunex/
├── core/           # Execution contracts, context, exceptions
├── analysis/       # Domain algorithms (15 submodules)
│   ├── allostery/
│   ├── angles/
│   ├── comparison/
│   ├── conformation/
│   ├── free_energy/
│   ├── interactions/
│   ├── interface/
│   ├── quality/
│   ├── reporter/
│   ├── structure/
│   ├── topology/
│   └── trajectory/
├── pipeline/       # Workflow nodes and orchestration
├── cli/            # User-facing CLI commands
├── utils/          # Utility modules
├── cluster/        # HPC integration
└── data/           # Reference data

scripts/            # Support scripts (minimal, stable)
development/        # Temporary development and testing
docs/               # Documentation
```

## Development Workflow

### 1. Primary Entry Point

Use `imn` CLI commands, not standalone scripts:

```bash
imn preprocess -f md.xtc -s md.tpr -o ./output
imn quality -f md_processed.xtc -s md.tpr -o ./output
imn rmsd / rmsf / contact / bsa / nma / inter_cluster ...
imn batch preprocess /path/to/tasks --workers 4
imn report interaction --base-dir ./output --system-id 1OGA
```

### 2. Adding New Features

New capabilities should land in:
- `immunex/analysis/` - Core algorithms
- `immunex/pipeline/nodes/` - Pipeline nodes
- `immunex/cli/commands/` - CLI commands

**Not** in `scripts/` or project root.

### 3. Testing and Development

- Temporary scripts → `development/`
- Unit tests → `test/` (if exists) or inline testing
- One-off analysis → `development/utilities/` or `development/case_studies/`

## Module Design Standards

All core modules must follow these principles:

### 1. Clear Inputs
- Use type hints for all parameters
- Use `dataclass` for structured inputs
- Implement `validate()` method
- Provide reasonable defaults

### 2. Clear Outputs
- Return structured data (dataclass), not strings/None
- Include `success: bool` field
- Include processing stats and metadata
- Include error messages if failed

### 3. Clear Side Effects
- Document all file operations
- Track commands executed
- Support dry-run mode (optional)

### 4. Clear Errors
- Define specific exception types
- Distinguish user errors from system errors
- Provide recovery suggestions

### 5. Testable
- Support unit testing
- Mock external dependencies (GROMACS, etc.)
- Test validation, success, and failure scenarios

### 6. Schedulable
- Support async/parallel execution
- Provide progress callbacks
- Support cancellation

## Code Quality Guidelines

### Naming Conventions

- **Functions/methods**: `snake_case`
- **Classes**: `PascalCase`
- **Constants**: `UPPER_SNAKE_CASE`
- **Private members**: `_leading_underscore`

### File Naming

- **Single task**: `pbc_process.py`
- **Batch processing**: `batch_pbc.py`
- **Cluster processing**: `batch_pbc_slurm.py`

Avoid: `process_v2.py`, `new_process.py`, `exact_process.py`

### Import Organization

```python
# Standard library
import os
from pathlib import Path

# Third-party
import numpy as np
import pandas as pd
import MDAnalysis as mda

# Local
from immunex.core import PipelineContext
from immunex.analysis.trajectory import RMSDAnalyzer
```

## Common Patterns

### 1. Analysis Module Template

```python
from dataclasses import dataclass
from typing import Optional
import MDAnalysis as mda

@dataclass
class AnalysisInput:
    topology: str
    trajectory: str
    selection: str
    output_file: str

    def validate(self):
        if not Path(self.topology).exists():
            raise FileNotFoundError(f"Topology not found: {self.topology}")

@dataclass
class AnalysisResult:
    success: bool
    output_file: str
    stats: dict
    error_message: Optional[str] = None

class MyAnalyzer:
    def analyze(self, input_params: AnalysisInput) -> AnalysisResult:
        try:
            input_params.validate()
            # Do analysis
            return AnalysisResult(
                success=True,
                output_file=output_path,
                stats={'mean': 2.34, 'std': 0.56}
            )
        except Exception as e:
            return AnalysisResult(
                success=False,
                output_file="",
                stats={},
                error_message=str(e)
            )
```

### 2. Pipeline Node Template

```python
from immunex.pipeline.base_pipeline import PipelineNode
from immunex.core import PipelineContext

class MyAnalysisNode(PipelineNode):
    def execute(self, context: PipelineContext) -> PipelineContext:
        # 1. Validate inputs
        if 'trajectory_processed' not in context:
            context.errors.append("Missing processed trajectory")
            return context

        # 2. Call core module
        analyzer = MyAnalyzer()
        result = analyzer.analyze(...)

        # 3. Record outputs
        context.results['my_analysis'] = {
            'output_file': result.output_file,
            'stats': result.stats
        }

        # 4. Handle errors
        if not result.success:
            context.errors.append(result.error_message)

        return context
```

### 3. CLI Command Template

```python
import click
from immunex.analysis.my_module import MyAnalyzer

@click.command()
@click.option('-f', '--trajectory', required=True)
@click.option('-s', '--topology', required=True)
@click.option('-o', '--output', required=True)
def my_command(trajectory, topology, output):
    """My analysis command description."""
    analyzer = MyAnalyzer()
    result = analyzer.analyze(...)

    if result.success:
        click.echo(f"Analysis completed: {result.output_file}")
    else:
        click.echo(f"Error: {result.error_message}", err=True)
```

## Key Documentation

Start with these documents:

- `README.md` - Project overview and quick start
- `docs/architecture/PRODUCT_SCOPE.md` - Product boundaries
- `docs/architecture/ENGINEERING_ARCHITECTURE.md` - Architecture design
- `docs/guides/IMN_CLI_GUIDE.md` - CLI usage guide
- `docs/guides/INSTALLATION_GUIDE.md` - Installation instructions
- `docs/guides/INTERACTION_MODULE_GUIDE.md` - Interaction analysis guide

## HTML Report System

### Overview

Immunex provides a comprehensive **single-system HTML report** with an embedded **AI assistant drawer** for interactive queries.

**Report Generation**:
```bash
imn report interaction \
  --base-dir output/interaction_case \
  --system-id 1OGA_sd_run2 \
  --source-pdb output/preprocess/1OGA/md_processed_converted.pdb \
  --bsa-root output/bsa_demo \
  --rmsf-root output/rmsf_demo \
  --identity-root output/identity_demo \
  --rrcs-root output/rrcs_demo \
  --cluster-root output/cluster_demo
```

**Report Server with AI Assistant**:
```bash
imn report serve \
  --overview-dir output/report/overview \
  --host 127.0.0.1 \
  --port 8770
```

### Report Sections (10 modules)

1. **overview** - System overview and biological identity
2. **quality** - Trajectory quality and RMSD convergence
3. **interface** - Buried surface area (BSA) analysis
4. **flexibility** - RMSF and flexibility analysis
5. **rrcs** - RRCS hotspot analysis
6. **cluster** - Interface clustering and dominant states
7. **occupancy** - Interaction occupancy statistics
8. **contact** - Coarse-grained contact analysis
9. **interactions** - Typed interactions (6 families)
10. **downloads** - Data download directory

### Interaction Families (6 types)

- **contact** - Residue-residue minimum heavy atom distance
- **hbond** - Hydrogen bonds (MDAnalysis HydrogenBondAnalysis)
- **saltbridge** - Salt bridges (charged atom distance)
- **hydrophobic** - Hydrophobic contacts (sidechain heavy atoms)
- **pipi** - Pi-Pi stacking (aromatic ring geometry)
- **cationpi** - Cation-Pi interactions (ring-charge geometry)

Each family provides:
- Residue-pair tables with occupancy
- Region-level summaries (CDR1/CDR2/CDR3/non-CDR)
- Interface views (peptide-TCR, HLA-TCR, groove-TCR)
- TCRα/TCRβ heatmaps

### AI Assistant Features

**Drawer UI**:
- Fixed floating action button (FAB) in bottom-right corner
- Slide-in drawer panel (420px width)
- Context-aware quick actions
- Evidence-based responses

**Query Types** (via `QueryRouter`):
- `cluster_status` - Dominant cluster identification
- `quality_status` - Trajectory stability diagnosis
- `interface_summary` - Interface overview
- `hotspot_summary` - RRCS hotspot identification
- `design_hint` - Mutation design suggestions
- `top_pair` - Strongest RRCS pairs
- `persistent_pair` - Most stable interface pairs
- `top_region` - Most active regions
- `top_residue` - Key residues

**Answer Structure** (via `QueryAnswerBuilder`):
- Direct answer sentence
- 3-5 evidence bullets
- Top 5 ranked items (if applicable)
- Confidence level
- Source citations

### Skills System

Located in `skills/` directory, providing specialized AI capabilities:

**Available Skills**:
1. **reporter-query** - Lightweight lookup questions
   - Identify key residues, regions, clusters
   - Find strongest hotspots and persistent pairs
   - Quick single-system queries

2. **reporter-diagnostic** - Full trajectory diagnosis
   - Stability assessment
   - Quality evaluation
   - Mechanistic interpretation

3. **reporter-design** - Mutation design suggestions
   - Hotspot-based design hints
   - Interface engineering recommendations

4. **reporter-compare** - Two-system comparison
   - Differential analysis
   - Comparative interpretation

Each skill includes:
- `SKILL.md` - Skill description and usage
- `references/` - Query routing, evidence rules, examples
- `agents/` - Agent configurations (if applicable)
- `scripts/` - Supporting scripts

## Common Tasks

### Running Analysis

```bash
# Preprocess trajectory (PBC correction)
imn preprocess -f md.xtc -s md.tpr -o ./output/preprocess

# Quality assessment
imn quality -f md_processed.xtc -s md.tpr -o ./output/quality

# RMSD analysis
imn rmsd -f md_processed.xtc -s md.tpr -o ./output/rmsd

# Batch interaction analysis (6 families)
imn batch contact /path/to/tasks --workers 4 -o ./output/contact
imn batch hbond /path/to/tasks --workers 4 -o ./output/hbond
imn batch saltbridge /path/to/tasks --workers 4 -o ./output/saltbridge
imn batch hydrophobic /path/to/tasks --workers 4 -o ./output/hydrophobic
imn batch pipi /path/to/tasks --workers 4 -o ./output/pipi
imn batch cationpi /path/to/tasks --workers 4 -o ./output/cationpi

# Generate HTML report
imn report interaction --base-dir ./output --system-id 1OGA

# Serve report with AI assistant
imn report serve --overview-dir ./output/report/overview --port 8770
```

### Development Commands

```bash
# Install in development mode
pip install -e .

# Run tests (if test suite exists)
pytest test/

# Format code
black immunex/
isort immunex/

# Check code quality
flake8 immunex/
```

## Important Notes

1. **PBC correction is mandatory**: All trajectory analysis should use PBC-corrected trajectories
2. **Quality assessment first**: Run quality checks before heavy analysis
3. **Use CLI commands**: Prefer `imn` commands over direct script execution
4. **Batch processing**: Use `imn batch` for multiple tasks
5. **Report generation**: Use `imn report interaction` for integrated reports

## What NOT to Do

1. Don't create scripts in project root or `scripts/` for testing
2. Don't create temporary markdown summary files
3. Don't add version suffixes to filenames
4. Don't duplicate functionality across multiple scripts
5. Don't skip PBC correction for trajectory analysis
6. Don't use Chinese or emojis in code/comments

## Current Development Priorities

Focus on:
- CLI consistency and documentation
- Report standardization
- Cross-module integration
- Batch summarization
- Cleanup of historical code

See `docs/architecture/PLATFORM_TODO_ROADMAP.md` for detailed roadmap.
