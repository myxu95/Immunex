# AfterMD Scripts Directory

This directory contains core production scripts for MD trajectory analysis. All scripts are stable, well-tested, and intended for regular use.

## Organization Principle

- **scripts/**: Production-ready tools (18 core scripts)
- **development/**: Experimental, testing, one-time analysis scripts
- **examples/**: Example usage scripts for documentation

## Script Categories

### Core Processing (3 scripts)

#### pbc_process.py
Single trajectory PBC (Periodic Boundary Conditions) correction.
```bash
python scripts/pbc_process.py -f trajectory.xtc -s topology.tpr -o output/
```

#### batch_pbc_slurm.py
Generate SLURM batch scripts for cluster processing.
```bash
python scripts/batch_pbc_slurm.py --input_dir data/ --output_dir results/
```

#### batch_worker.py
Generic SLURM worker for distributed task execution.
```bash
python scripts/batch_worker.py --config task_config.json
```

### Quality Control (3 scripts)

#### md_quality_check.py
MD simulation quality validation (completeness, energy, structure).
```bash
python scripts/md_quality_check.py --input_dir md_results/ --report quality_report.json
```

#### md_task_organizer.py
Organize and categorize MD tasks for batch processing.
```bash
python scripts/md_task_organizer.py --tasks_dir raw_data/ --output organized_tasks/
```

#### md_workflow.py
Complete MD analysis workflow orchestration.
```bash
python scripts/md_workflow.py --config workflow_config.yaml
```

### Batch Analysis (9 scripts)

#### batch_allostery_analysis.py
Dynamic contact correlation analysis for allostery studies.
```bash
python scripts/batch_allostery_analysis.py --tasks tasks.json --workers 4
```

#### batch_cdr_rmsd_exact.py
CDR (Complementarity Determining Region) RMSD calculation.
**Note**: Consolidated version, replaces older variants.
```bash
python scripts/batch_cdr_rmsd_exact.py --input_dir cdr_data/ --align_strategy exact
```

#### batch_cdr_rmsf.py
CDR RMSF (Root Mean Square Fluctuation) analysis.
```bash
python scripts/batch_cdr_rmsf.py --tasks tasks.json
```

#### batch_contact_frequency.py
Residue contact frequency calculation across trajectories.
```bash
python scripts/batch_contact_frequency.py --selection_a "segname PROA" --selection_b "segname PROC"
```

#### batch_interface_rmsd.py
Protein-protein interface RMSD analysis.
```bash
python scripts/batch_interface_rmsd.py --tasks interface_tasks.json
```

#### batch_phla_analysis.py
pMHC (peptide-MHC) complex specific analysis.
```bash
python scripts/batch_phla_analysis.py --input_dir phla_data/
```

#### batch_rmsd_hla_alpha.py
HLA alpha chain RMSD calculation.
```bash
python scripts/batch_rmsd_hla_alpha.py --tasks hla_tasks.json
```

#### batch_tcr_rmsd.py
T Cell Receptor (TCR) RMSD analysis.
```bash
python scripts/batch_tcr_rmsd.py --tasks tcr_tasks.json
```

#### batch_whole_protein_rmsf.py
Full protein RMSF analysis (all residues).
```bash
python scripts/batch_whole_protein_rmsf.py --input_dir trajectories/
```

### Analysis Pipelines (3 scripts)

#### analyze_chain_contacts.py
Detailed chain-chain contact analysis with VMD integration.
```bash
python scripts/analyze_chain_contacts.py --topology md.tpr --trajectory md_pbc.xtc
```

#### rmsf_analysis_pipeline.py
Complete RMSF analysis pipeline with visualization.
```bash
python scripts/rmsf_analysis_pipeline.py --tasks rmsf_tasks.json --output results/
```

#### run_trajectory_analysis.py
General trajectory analysis launcher.
```bash
python scripts/run_trajectory_analysis.py --config analysis_config.yaml
```

### VMD Integration

#### vmd_chain_contacts.tcl
VMD script for chain contact analysis.
```bash
vmd -dispdev text -e scripts/vmd_chain_contacts.tcl
```
See: `docs/VMD_CHAIN_CONTACTS_GUIDE.md`

## Shell Scripts

#### monitor_batch.sh
Monitor SLURM batch job progress.
```bash
bash scripts/monitor_batch.sh job_ids.txt
```

#### run_batch_cdr_with_env.sh
Run CDR batch analysis with environment setup.
```bash
bash scripts/run_batch_cdr_with_env.sh tasks.json
```

## Usage Guidelines

### 1. Production Scripts Only
- Only stable, tested scripts belong here
- No test_*.py or experimental scripts
- No version suffixes (_v2, _new, _exact should be consolidated)

### 2. Prefer Module Functions
When possible, import from aftermd package instead of duplicating code:
```python
from aftermd.core import PBCProcessor
from aftermd.analysis.trajectory import RMSDCalculator
from aftermd.utils import BatchProcessor
```

### 3. Consistent CLI Interface
All scripts should follow these conventions:
- Use argparse for command-line arguments
- Support --help flag with detailed usage
- Accept --config for JSON/YAML configuration
- Provide --verbose for debugging output
- Use --workers for parallel processing control

### 4. Logging Standards
```python
import logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)
```

## Development Workflow

### Adding New Scripts

1. **Test first**: Create script in `development/` directory
2. **Validate**: Run on multiple datasets, ensure robustness
3. **Document**: Add docstrings and --help text
4. **Review**: Get code review before moving to `scripts/`
5. **Update**: Add entry to this README

### Removing Scripts

1. **Deprecate**: Move to `development/archived_scripts/`
2. **Document**: Note reason in `docs/PROJECT_CLEANUP_*.md`
3. **Update**: Remove from this README
4. **Wait**: Keep archived for 1 month before deletion

## Related Documentation

- `ARCHITECTURE.md` - Overall project architecture
- `docs/CDR3_RMSD_PLOTTING_GUIDE.md` - CDR3 RMSD analysis
- `docs/PBC_PROCESSING_METHODS.md` - PBC correction methods
- `docs/VMD_CHAIN_CONTACTS_GUIDE.md` - VMD integration guide
- `docs/PROJECT_CLEANUP_2026-01.md` - Recent cleanup report

## Common Workflows

### Complete Analysis Pipeline
```bash
# 1. Quality check
python scripts/md_quality_check.py --input_dir raw_md_data/

# 2. PBC correction
python scripts/batch_pbc_slurm.py --input_dir qualified_tasks/

# 3. RMSD analysis
python scripts/batch_cdr_rmsd_exact.py --input_dir pbc_processed/

# 4. Contact analysis
python scripts/batch_contact_frequency.py --tasks contact_tasks.json

# 5. Allostery analysis
python scripts/batch_allostery_analysis.py --tasks allostery_tasks.json
```

### TCR-pMHC Specific Analysis
```bash
# CDR RMSD
python scripts/batch_cdr_rmsd_exact.py --tasks tcr_tasks.json

# TCR RMSD
python scripts/batch_tcr_rmsd.py --tasks tcr_tasks.json

# Interface RMSD
python scripts/batch_interface_rmsd.py --tasks interface_tasks.json

# Chain contacts
python scripts/analyze_chain_contacts.py --topology md.tpr --trajectory md_pbc.xtc
```

## Maintenance

**Last Cleanup**: 2026-01-20
**Script Count**: 18 Python scripts + 3 shell/TCL scripts
**Next Review**: 2026-04-20 (quarterly)

For historical changes, see `docs/PROJECT_CLEANUP_2026-01.md`
