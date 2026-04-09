# Immunex Software Architecture

**Version:** 0.1.0
**Last Updated:** 2025-10-12

---

## Table of Contents

1. [Overview](#overview)
2. [Architecture Principles](#architecture-principles)
3. [Module Organization](#module-organization)
4. [Quality Control Workflow](#quality-control-workflow)
5. [Data Flow](#data-flow)
6. [Module Details](#module-details)
7. [Integration Points](#integration-points)
8. [Extensibility](#extensibility)

---

## Overview

Immunex is a **GROMACS MD Analysis Toolkit** designed for large-scale batch processing and HPC cluster deployment. The architecture follows a **modular, layered design** with clear separation between utilities, preprocessing, analysis, and quality control.

### Key Features

- **Staged Quality Control**: Pre-PBC and Post-PBC validation
- **PBC-Aware Architecture**: Clear dependency tracking
- **Batch Processing**: Parallel execution with configurable workers
- **HPC Integration**: SLURM script generation and cluster management
- **Modular Design**: Independent, reusable components

---

## Architecture Principles

### 1. Layered Architecture

```
┌─────────────────────────────────────────────────────────────┐
│                     Application Layer                        │
│  (Scripts, Examples, User Interfaces)                        │
└─────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────┐
│                    Quality Control Layer                     │
│  Pre-PBC Quality Check → Post-PBC Quality Check             │
└─────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────┐
│                     Processing Layer                         │
│  PBC Processing → Trajectory Analysis → Structure Analysis  │
└─────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────┐
│                       Utility Layer                          │
│  Batch Processing, Path Management, Plotting, etc.          │
└─────────────────────────────────────────────────────────────┘
```

### 2. Dependency Awareness

**PBC-Independent Modules** (fast screening):
- `EnergyQualityChecker`
- `MDCompletenessChecker`
- `StructureValidator`

**PBC-Dependent Modules** (require processed trajectory):
- `RMSDCalculator`
- `RadiusGyrationCalculator`
- `TrajectoryStabilityChecker` (future)

### 3. Design Patterns

- **Factory Pattern**: Module instantiation with configurable parameters
- **Strategy Pattern**: Pluggable analysis algorithms
- **Observer Pattern**: Quality reporting and issue tracking
- **Template Method**: Standardized batch processing workflow

---

## Module Organization

```
immunex/
├── __init__.py                           # Main package exports
│
├── utils/                                # UTILITY LAYER
│   ├── __init__.py
│   ├── batch_processor.py               # Parallel batch processing
│   ├── path_manager.py                  # File organization
│   ├── plotting.py                      # Visualization
│   ├── group_selector.py                # GROMACS group selection
│   ├── cleanup_manager.py               # Temporary file cleanup
│   ├── slurm_generator.py               # SLURM script generation
│   ├── slurm_pbc_generator.py           # PBC-specific SLURM
│   └── slurm_rmsd_generator.py          # RMSD-specific SLURM
│
├── preprocessing/                        # PREPROCESSING LAYER
│   ├── __init__.py
│   └── pbc_processor.py                 # PBC artifact removal
│
├── analysis/                             # ANALYSIS LAYER
│   ├── __init__.py
│   ├── analysis_pipeline.py             # Integrated analysis workflow
│   ├── comparative_analysis.py          # Multi-trajectory comparison
│   │
│   ├── quality/                          # Quality Control Modules
│   │   ├── __init__.py
│   │   ├── md_completeness.py           # [Pre-PBC] File integrity
│   │   ├── structure_validator.py       # [Pre-PBC] Structure validation
│   │   ├── energy_quality.py            # [Pre-PBC] Energy-based QC ⭐
│   │   ├── batch_tracker.py             # Batch processing tracking
│   │   └── quality_reporter.py          # Comprehensive reporting
│   │
│   ├── trajectory/                       # Trajectory Analysis (PBC-dependent)
│   │   ├── __init__.py
│   │   ├── rmsd.py                      # [Post-PBC] RMSD calculation
│   │   ├── rdf.py                       # [Post-PBC] Radial distribution
│   │   ├── radius_gyration.py           # [Post-PBC] Radius of gyration
│   │   ├── distance.py                  # [Post-PBC] Distance analysis
│   │   └── hydrogen_bonds.py            # [Post-PBC] H-bond analysis
│   │
│   └── structure/                        # Structure Analysis (PDB-based)
│       ├── __init__.py
│       ├── bfactor.py                   # B-factor extraction
│       ├── contact_map.py               # Contact map calculation
│       ├── geometry.py                  # Geometric properties
│       └── atom_info.py                 # Atom information extraction
│
└── batch_process.py                      # High-level batch interface

scripts/                                  # EXECUTABLE SCRIPTS
├── pbc_process.py                       # Single trajectory PBC
├── batch_pbc_slurm.py                   # Cluster batch PBC
├── md_quality_check.py                  # Quality check workflow
├── rmsd_calculator.py                   # RMSD calculation
├── run_trajectory_analysis.py           # Trajectory analysis pipeline
└── run_comparative_analysis.py          # Comparative analysis

examples/                                 # USAGE EXAMPLES
├── energy_quality_check_usage.py        # Energy QC examples ⭐
├── quality_analysis_usage.py            # Quality analysis demo
├── toolkit_comprehensive_example.py     # Full toolkit demo
└── modular_usage.py                     # Modular usage patterns

templates/                                # CONFIGURATION TEMPLATES
└── slurm_template.sh                    # SLURM job template
```

**Legend:**
- ⭐ = Newly implemented module
- [Pre-PBC] = Independent of PBC processing
- [Post-PBC] = Requires PBC-processed trajectory

---

## Quality Control Workflow

### Three-Stage Quality Assessment

```
┌──────────────────────────────────────────────────────────────────────┐
│                          STAGE 1: Pre-PBC Quality Check              │
│                          (Fast Screening - No PBC Required)          │
└──────────────────────────────────────────────────────────────────────┘
                                    ↓
        ┌───────────────────────────────────────────────────┐
        │  1.1 File Completeness Check                      │
        │      - md.gro exists                              │
        │      - md.xtc size > 1MB                          │
        │      - md.log contains completion marker          │
        │      Module: MDCompletenessChecker                │
        └───────────────────────────────────────────────────┘
                                    ↓
        ┌───────────────────────────────────────────────────┐
        │  1.2 Energy Quality Check ⭐                       │
        │      - Temperature stability (T ± 5K)             │
        │      - Pressure stability (P ± 50 bar, NPT)       │
        │      - Energy conservation (drift < 2%)           │
        │      - Potential/Kinetic balance                  │
        │      Module: EnergyQualityChecker                 │
        └───────────────────────────────────────────────────┘
                                    ↓
        ┌───────────────────────────────────────────────────┐
        │  1.3 Structure Integrity Check                    │
        │      - PDB chain count validation                 │
        │      - Coordinate range verification              │
        │      Module: StructureValidator                   │
        └───────────────────────────────────────────────────┘
                                    ↓
                        ┌───────────────────────┐
                        │  Decision Point 1     │
                        │  Pass Pre-QC?         │
                        └───────────────────────┘
                                    │
                    ┌───────────────┴──────────────┐
                    NO                            YES
                    ↓                              ↓
            ┌───────────────┐            ┌─────────────────┐
            │ REJECT        │            │ Proceed to      │
            │ Skip PBC      │            │ PBC Processing  │
            │ processing    │            └─────────────────┘
            └───────────────┘                      ↓

┌──────────────────────────────────────────────────────────────────────┐
│                       STAGE 2: PBC Processing                        │
└──────────────────────────────────────────────────────────────────────┘
                                    ↓
        ┌───────────────────────────────────────────────────┐
        │  2.1 Shortest Chain Detection                     │
        │      - gmx make_ndx + splitch                     │
        │      - Identify peptide (shortest chain)          │
        │      Module: GroupSelector                        │
        └───────────────────────────────────────────────────┘
                                    ↓
        ┌───────────────────────────────────────────────────┐
        │  2.2 Three-Step PBC Removal                       │
        │      Step 1: Center on peptide (pbc atom)         │
        │      Step 2: Make molecules whole (pbc whole)     │
        │      Step 3: Fit rotation+translation             │
        │      Module: PBCProcessor                         │
        └───────────────────────────────────────────────────┘
                                    ↓
        ┌───────────────────────────────────────────────────┐
        │  2.3 Reference File Management                    │
        │      - Copy md.tpr to output                      │
        │      - Copy md.gro to output                      │
        │      Module: PBCProcessor                         │
        └───────────────────────────────────────────────────┘
                                    ↓

┌──────────────────────────────────────────────────────────────────────┐
│                     STAGE 3: Post-PBC Quality Check                  │
│                     (Trajectory-Based Validation)                    │
└──────────────────────────────────────────────────────────────────────┘
                                    ↓
        ┌───────────────────────────────────────────────────┐
        │  3.1 RMSD Quality Assessment                      │
        │      - RMSD mean < 0.5 nm                         │
        │      - RMSD std < 0.3 nm                          │
        │      - No systematic drift                        │
        │      Module: RMSDCalculator                       │
        └───────────────────────────────────────────────────┘
                                    ↓
        ┌───────────────────────────────────────────────────┐
        │  3.2 Radius of Gyration Check                     │
        │      - Rg change < 15%                            │
        │      - Structure compactness maintained           │
        │      Module: RadiusGyrationCalculator             │
        └───────────────────────────────────────────────────┘
                                    ↓
        ┌───────────────────────────────────────────────────┐
        │  3.3 Convergence Analysis                         │
        │      - First half vs second half RMSD             │
        │      - Kolmogorov-Smirnov test (p > 0.05)         │
        │      Module: TrajectoryStabilityChecker (future)  │
        └───────────────────────────────────────────────────┘
                                    ↓
        ┌───────────────────────────────────────────────────┐
        │  3.4 Coordinate Continuity Check                  │
        │      - Max atom jump < 1.0 nm                     │
        │      - No PBC artifacts detected                  │
        │      Module: TrajectoryStabilityChecker (future)  │
        └───────────────────────────────────────────────────┘
                                    ↓
                        ┌───────────────────────┐
                        │  Decision Point 2     │
                        │  Pass Post-QC?        │
                        └───────────────────────┘
                                    │
                    ┌───────────────┴──────────────┐
                    NO                            YES
                    ↓                              ↓
            ┌───────────────┐            ┌─────────────────┐
            │ Grade: C/D    │            │ Grade: A/B      │
            │ Use with      │            │ Ready for       │
            │ caution       │            │ publication     │
            └───────────────┘            └─────────────────┘
```

### Quality Grading System

| Grade | Score    | Pre-PBC Criteria | Post-PBC Criteria | Recommendation |
|-------|----------|------------------|-------------------|----------------|
| **A** | 90-100   | Energy: A<br>Temp: ±2K<br>Drift: <0.5% | RMSD: <0.3nm<br>Rg: <5% change<br>KS test: p>0.1 | Publication-ready |
| **B** | 75-89    | Energy: B<br>Temp: ±5K<br>Drift: <1.0% | RMSD: <0.5nm<br>Rg: <10% change<br>KS test: p>0.05 | Good for analysis |
| **C** | 60-74    | Energy: C<br>Temp: ±10K<br>Drift: <2.0% | RMSD: <0.8nm<br>Rg: <15% change<br>KS test: p>0.01 | Usable with caution |
| **D** | <60      | Energy: D<br>Poor stability | RMSD: >0.8nm<br>Rg: >15% change<br>Not converged | Re-simulation recommended |

### Scoring Weights

**Pre-PBC Score (40%)**:
- File Completeness: 10%
- Energy Quality: 25%
- Structure Integrity: 5%

**Post-PBC Score (60%)**:
- RMSD Quality: 20%
- Rg Quality: 15%
- Convergence: 15%
- Coordinate Quality: 10%

---

## Data Flow

### Input → Processing → Output Pipeline

```
┌──────────────────────────────────────────────────────────────────────┐
│                            INPUT FILES                               │
└──────────────────────────────────────────────────────────────────────┘
                                    ↓
    ┌────────────────────────────────────────────────────────┐
    │  MD Simulation Directory Structure:                    │
    │                                                         │
    │  task_directory/                                        │
    │  ├── {pdbid}.pdb          (structure file)             │
    │  ├── md.gro               (final structure)            │
    │  ├── md.xtc               (trajectory)                 │
    │  ├── md.tpr               (topology)                   │
    │  ├── md.edr               (energy data) ⭐              │
    │  └── md.log               (simulation log)             │
    │                                                         │
    │  OR with prod/ subfolder:                              │
    │  task_directory/                                        │
    │  ├── {pdbid}.pdb                                       │
    │  └── prod/                                             │
    │      ├── md.gro                                        │
    │      ├── md.xtc                                        │
    │      ├── md.tpr                                        │
    │      ├── md.edr           ⭐                            │
    │      └── md.log                                        │
    └────────────────────────────────────────────────────────┘
                                    ↓
┌──────────────────────────────────────────────────────────────────────┐
│                       PROCESSING WORKFLOW                            │
└──────────────────────────────────────────────────────────────────────┘
                                    ↓
    ┌────────────────────────────────────────────────────────┐
    │  Stage 1: Pre-PBC Quality Check                        │
    │  Input:  md.edr, md.log, md.gro, {pdbid}.pdb          │
    │  Output: quality_report.json                           │
    │          qualified_mds.txt                             │
    │          failed_mds.txt                                │
    │  Time:   ~5 sec/task                                   │
    └────────────────────────────────────────────────────────┘
                                    ↓
    ┌────────────────────────────────────────────────────────┐
    │  Stage 2: PBC Processing                               │
    │  Input:  md.xtc, md.tpr (from qualified tasks)        │
    │  Output: {task}_processed.xtc                         │
    │          md.gro (copied)                              │
    │          md.tpr (copied)                              │
    │          chains.ndx (temporary, cleaned up)           │
    │  Time:   ~1-5 min/task (depends on trajectory size)   │
    └────────────────────────────────────────────────────────┘
                                    ↓
    ┌────────────────────────────────────────────────────────┐
    │  Stage 3: Post-PBC Analysis                            │
    │  Input:  {task}_processed.xtc, md.tpr                 │
    │  Output: rmsd.xvg                                      │
    │          rg.xvg                                        │
    │          distance.xvg                                  │
    │          analysis_plots.png                            │
    │  Time:   ~30 sec/task                                  │
    └────────────────────────────────────────────────────────┘
                                    ↓
┌──────────────────────────────────────────────────────────────────────┐
│                          OUTPUT STRUCTURE                            │
└──────────────────────────────────────────────────────────────────────┘
                                    ↓
    ┌────────────────────────────────────────────────────────┐
    │  output_directory/                                     │
    │  ├── quality_check_results/                           │
    │  │   ├── quality_check_report.json                    │
    │  │   ├── qualified_mds.txt                            │
    │  │   ├── failed_mds.txt                               │
    │  │   └── quality_summary.txt                          │
    │  │                                                     │
    │  ├── task1/                                           │
    │  │   ├── task1_processed.xtc       (PBC-corrected)   │
    │  │   ├── md.gro                    (reference)        │
    │  │   ├── md.tpr                    (topology)         │
    │  │   ├── rmsd.xvg                  (analysis)         │
    │  │   ├── rg.xvg                                       │
    │  │   └── processing_log.txt                           │
    │  │                                                     │
    │  ├── task2/                                           │
    │  │   └── ...                                          │
    │  │                                                     │
    │  └── slurm_scripts/                (if using cluster) │
    │      ├── immunex_batch_1.sh                           │
    │      ├── immunex_batch_2.sh                           │
    │      └── submit_all_batches.sh                        │
    └────────────────────────────────────────────────────────┘
```

### Parallel Batch Processing

```
                    ┌─────────────────────────┐
                    │   BatchProcessor        │
                    │   (max_workers=4)       │
                    └─────────────────────────┘
                              ↓
        ┌─────────────────────┼─────────────────────┐
        ↓                     ↓                     ↓
    Worker 1              Worker 2              Worker 3
    task1                 task2                 task3
    task5                 task6                 task7
    task9                 task10                task11
        ↓                     ↓                     ↓
    ┌──────────────────────────────────────────────────┐
    │         Results Collection & Summary             │
    └──────────────────────────────────────────────────┘
```

---

## Module Details

### 1. Quality Control Modules

#### 1.1 EnergyQualityChecker ⭐ (NEW)

**Purpose**: Pre-PBC energy-based quality assessment
**Dependencies**: None (uses md.edr)
**Performance**: Fast (~2-5 sec/task)

```python
from immunex.analysis.quality import EnergyQualityChecker

checker = EnergyQualityChecker(
    target_temperature=300.0,      # Target simulation temperature (K)
    temp_tolerance=5.0,            # Temperature deviation tolerance (K)
    target_pressure=1.0,           # Target pressure for NPT (bar)
    pressure_tolerance=50.0,       # Pressure deviation tolerance (bar)
    energy_drift_threshold=2.0     # Max energy drift percentage
)

# Comprehensive check
result = checker.comprehensive_energy_check("md.edr")

# Individual checks
temp_result = checker.check_temperature_stability("md.edr")
pressure_result = checker.check_pressure_stability("md.edr")
energy_result = checker.check_total_energy_conservation("md.edr")
balance_result = checker.check_potential_kinetic_balance("md.edr")
```

**Output Structure**:
```python
{
    'status': 'success',
    'energy_grade': 'A',           # A/B/C/D
    'score': 92.5,                 # 0-100
    'overall_status': 'excellent',  # excellent/good/acceptable/poor
    'issues': [],                   # List of detected issues
    'details': {
        'temperature': {
            'grade': 'A',
            'statistics': {'mean': 300.1, 'std': 1.8},
            'deviation': 0.1,
            'in_tolerance': True
        },
        'energy_conservation': {
            'grade': 'A',
            'drift': {'drift_percentage': 0.3}
        }
    }
}
```

#### 1.2 MDCompletenessChecker

**Purpose**: Verify MD simulation completeness
**Dependencies**: None (file system)
**Performance**: Very fast (~1 sec/task)

```python
from immunex.analysis.quality import MDCompletenessChecker

checker = MDCompletenessChecker(
    min_trajectory_size_mb=1.0,
    min_simulation_time_ps=5000.0
)

result = checker.check_directory("md_directory")
```

#### 1.3 StructureValidator

**Purpose**: Validate PDB structure integrity
**Dependencies**: MDAnalysis (for PDB parsing)
**Performance**: Fast (~2-3 sec/structure)

```python
from immunex.analysis.quality import StructureValidator

validator = StructureValidator(expected_chain_count=5)
result = validator.validate_structure("structure.pdb")
```

### 2. Preprocessing Modules

#### 2.1 PBCProcessor

**Purpose**: Remove PBC artifacts using three-step process
**Dependencies**: GROMACS, GroupSelector
**Performance**: Moderate (~1-5 min/trajectory)

**Three-Step PBC Removal**:
1. **Center**: `gmx trjconv -center -pbc atom` (using peptide)
2. **Whole**: `gmx trjconv -pbc whole`
3. **Fit**: `gmx trjconv -fit rot+trans`

```python
from immunex.preprocessing import PBCProcessor

processor = PBCProcessor(gmx_executable="gmx")

# Comprehensive processing
result = processor.comprehensive_pbc_process(
    trajectory="md.xtc",
    topology="md.tpr",
    output_dir="processed_output",
    dt=10.0  # Optional: downsample to 10ps intervals
)
```

### 3. Analysis Modules

#### 3.1 Trajectory Analysis (Post-PBC)

**RMSDCalculator**:
```python
from immunex.analysis import RMSDCalculator

calc = RMSDCalculator(topology="md.tpr", trajectory="processed.xtc")

# MDAnalysis method
times, rmsd = calc.calculate_mdanalysis(selection="protein and name CA")

# GROMACS method
xvg_file = calc.calculate_gromacs(rmsd_type="backbone")
```

**RadiusGyrationCalculator**:
```python
from immunex.analysis import RadiusGyrationCalculator

calc = RadiusGyrationCalculator(topology="md.tpr", trajectory="processed.xtc")
times, rg = calc.calculate_radius_gyration(selection="protein")
```

#### 3.2 Structure Analysis (PDB-based)

**BFactorAnalyzer**:
```python
from immunex.analysis import BFactorAnalyzer

analyzer = BFactorAnalyzer(structure="protein.pdb")
bfactors = analyzer.extract_bfactors_by_residue()
```

### 4. Utility Modules

#### 4.1 BatchProcessor

**Purpose**: Parallel batch processing framework
**Performance**: Configurable workers (default: CPU count - 1)

```python
from immunex.utils import BatchProcessor

processor = BatchProcessor(max_workers=4)

# Discover MD tasks
tasks = processor.discover_batch_tasks("/data/simulations")

# Batch PBC processing
results = processor.batch_pbc_processing(
    base_directory="/data/simulations",
    output_base_dir="/data/processed",
    dt=10.0
)
```

#### 4.2 PlotManager

**Purpose**: Publication-ready visualization
**Supported Formats**: XVG, CSV, NumPy arrays

```python
from immunex.utils import PlotManager

plotter = PlotManager()
plotter.plot_xvg("rmsd.xvg", output="rmsd_plot.png")
```

---

## Integration Points

### High-Level Interfaces

#### 1. Simple Batch Processing

```python
from immunex import process_md_tasks

# One-line batch processing
results = process_md_tasks(
    simulations_path="/data/md_simulations",
    output_dir="/data/processed",
    dt=10.0,
    max_workers=4
)
```

#### 2. Quality-Aware Workflow

```python
from immunex.analysis.quality import (
    EnergyQualityChecker,
    MDCompletenessChecker
)
from immunex.preprocessing import PBCProcessor

# Stage 1: Pre-PBC screening
energy_checker = EnergyQualityChecker()
md_checker = MDCompletenessChecker()

qualified_tasks = []
for md_dir in all_md_dirs:
    # File completeness
    if md_checker.check_directory(md_dir)['status'] != 'complete':
        continue

    # Energy quality
    edr_file = f"{md_dir}/md.edr"
    energy_result = energy_checker.comprehensive_energy_check(edr_file)

    if energy_result['energy_grade'] in ['A', 'B', 'C']:
        qualified_tasks.append(md_dir)

# Stage 2: PBC processing (only qualified tasks)
processor = PBCProcessor()
for task in qualified_tasks:
    processor.comprehensive_pbc_process(...)
```

#### 3. SLURM Cluster Deployment

```python
from immunex import generate_slurm_scripts_for_md_tasks

# Generate SLURM scripts
results = generate_slurm_scripts_for_md_tasks(
    simulations_path="/data/simulations",
    tasks_per_batch=10,
    slurm_params={
        'partition': 'gpu',
        'time': '24:00:00',
        'cpus_per_task': 8
    }
)

# Output:
# slurm_scripts/immunex_batch_1.sh
# slurm_scripts/immunex_batch_2.sh
# slurm_scripts/submit_all_batches.sh
```

---

## Extensibility

### Adding New Quality Checkers

**Example: TrajectoryStabilityChecker (Post-PBC)**

```python
# immunex/analysis/quality/trajectory_stability.py

class TrajectoryStabilityChecker:
    """Post-PBC trajectory quality checker"""

    def __init__(self,
                 rmsd_max_threshold: float = 1.5,
                 rmsd_std_threshold: float = 0.3):
        self.rmsd_max = rmsd_max_threshold
        self.rmsd_std = rmsd_std_threshold

    def check_rmsd_quality(self, processed_trajectory: str,
                          topology: str) -> Dict:
        """
        RMSD-based quality check (requires PBC processing)

        Returns:
            Dictionary with RMSD quality metrics
        """
        # Implementation
        pass

    def check_convergence(self, processed_trajectory: str,
                         topology: str) -> Dict:
        """
        Trajectory convergence analysis

        Uses KS test to compare first/second half distributions
        """
        pass
```

### Plugin Architecture

```python
# Custom analysis module
from immunex.analysis.trajectory import BaseTrajectoryAnalyzer

class MyCustomAnalyzer(BaseTrajectoryAnalyzer):
    def __init__(self, trajectory, topology):
        super().__init__(trajectory, topology)

    def calculate_my_property(self):
        # Custom analysis logic
        return results
```

---

## Performance Characteristics

### Typical Processing Times (per trajectory)

| Stage | Module | Time | Bottleneck |
|-------|--------|------|------------|
| Pre-QC | File Check | 1 sec | File I/O |
| Pre-QC | Energy Check ⭐ | 2-5 sec | `gmx energy` |
| Pre-QC | Structure Check | 2-3 sec | MDAnalysis parsing |
| PBC | Chain Detection | 5-10 sec | `gmx make_ndx` |
| PBC | Three-Step Process | 1-5 min | Trajectory size |
| Post-QC | RMSD Calculation | 10-30 sec | Trajectory frames |
| Post-QC | Rg Calculation | 10-20 sec | Trajectory frames |

### Scalability

**Local Processing** (4 workers):
- 100 tasks: ~2-6 hours (with PBC)
- 100 tasks: ~10 minutes (Pre-QC only)

**Cluster Processing** (10 nodes, 10 tasks/node):
- 1000 tasks: ~2-6 hours (parallel)

---

## Version History

### v0.1.0 (Current)
- ✅ Core PBC processing pipeline
- ✅ Batch processing framework
- ✅ SLURM cluster integration
- ✅ Pre-PBC quality control (File, Energy ⭐, Structure)
- ✅ Post-PBC trajectory analysis (RMSD, RDF, Rg, Distance, H-bonds)
- ✅ Structure analysis (B-factor, Contact map, Geometry)

### Future Enhancements
- 🔄 TrajectoryStabilityChecker (Post-PBC convergence analysis)
- 🔄 QualityCheckWorkflow (integrated Pre/Post pipeline)
- 🔄 Real-time monitoring dashboard
- 🔄 Machine learning-based quality prediction

---

## References

### Key Algorithms

1. **Shortest Chain Detection**: Based on `gmx make_ndx -splitch`
2. **PBC Removal**: Three-step protocol (center → whole → fit)
3. **Energy Analysis**: Linear regression for drift detection
4. **Quality Grading**: Weighted scoring system (Pre: 40%, Post: 60%)

### Dependencies

- **GROMACS**: 2020+ (PBC processing, energy extraction)
- **MDAnalysis**: 2.6.0+ (trajectory/structure analysis)
- **NumPy**: 1.24+ (numerical operations)
- **Pandas**: 2.0+ (data management)
- **Matplotlib**: 3.7+ (visualization)

---

**Document Version**: 1.0
**Generated**: 2025-10-12
**Architecture Status**: Stable (Production-ready for Pre-PBC QC)
