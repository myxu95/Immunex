# Immunex

Immunex is a task-oriented analysis platform for **TCR-pMHC molecular dynamics trajectories**.  
Its goal is not to provide a loose collection of scripts, but to provide stable product-level entry points that connect:

- preprocessing
- quality assessment
- interface and interaction analysis
- biological identity annotation
- dynamical interpretation
- report generation

into one coherent workflow.

## Current Repository Status

The repository already exposes a complete product-facing command surface:

- `imn preprocess`: PBC handling and trajectory preprocessing
- `imn quality`: trajectory quality assessment
- `imn rmsd` / `imn rmsf`: basic dynamics analysis
- `imn contact`: coarse contact analysis
- `imn identity`: biological identity annotation
- `imn bsa`: buried surface area analysis
- `imn nma`: Normal Mode / ENM / PRS analysis
- `imn inter_cluster`: key interface-state clustering
- `imn batch ...`: batch workflows
- `imn report interaction`: integrated single-system HTML report

## Installation

Conda is the recommended setup path:

```bash
git clone <your-repo-url> Immunex
cd Immunex
conda env create -f environment.yml
conda activate immunex
pip install -e .
```

Basic verification:

```bash
python -c "import immunex; print(immunex.__version__)"
imn --help
gmx --version
```

For the full environment and validation guide, see:

- [docs/guides/INSTALLATION_GUIDE.md](docs/guides/INSTALLATION_GUIDE.md)

## Recommended Workflow

### 1. Preprocess

All trajectory analysis should be based on **PBC-corrected trajectories**.

```bash
imn preprocess -f md.xtc -s md.tpr -o ./output/preprocess_case
```

### 2. Quality Assessment

```bash
imn quality -f ./output/preprocess_case/md_processed.xtc -s md.tpr -o ./output/quality_case
```

### 3. Single-Module Analysis

```bash
imn rmsd -f ./output/preprocess_case/md_processed.xtc -s md.tpr -o ./output/rmsd_case
imn rmsf -f ./output/preprocess_case/md_processed.xtc -s md.tpr --structure ./output/preprocess_case/md_processed_converted.pdb -o ./output/rmsf_case
imn contact -f ./output/preprocess_case/md_processed.xtc -s md.tpr --structure ./output/preprocess_case/md_processed_converted.pdb -o ./output/contact_case
imn identity --structure ./output/preprocess_case/md_processed_converted.pdb -o ./output/identity_case
imn bsa -f ./output/preprocess_case/md_processed.xtc -s md.tpr --structure ./output/preprocess_case/md_processed_converted.pdb -o ./output/bsa_case
imn nma --structure ./output/preprocess_case/md_processed_converted.pdb -o ./output/nma_case
imn inter_cluster -f ./output/preprocess_case/md_processed.xtc -s md.tpr --structure ./output/preprocess_case/md_processed_converted.pdb -o ./output/inter_cluster_case
```

### 4. Batch Processing

```bash
imn batch preprocess /path/to/raw_tasks --workers 4 --output-dir ./output/preprocess_batch
imn batch quality ./output/preprocess_batch --workers 4 --output-dir ./output/quality_batch
imn batch bsa ./output/preprocess_batch --workers 4 --output-dir ./output/bsa_batch
```

### 5. Generate an Integrated Report

```bash
imn report interaction \
  --base-dir ./output/interaction_case \
  --system-id 1OGA_sd_run2 \
  --bsa-root ./output/bsa_case \
  --rmsf-root ./output/rmsf_case \
  --identity-root ./output/identity_case \
  --cluster-root ./output/inter_cluster_case
```

## Current Analysis Coverage

### Preprocessing and Quality

- PBC handling
- RMSD convergence and quality evaluation
- preprocessing quality reports

### Structure and Trajectory Analysis

- RMSD
- RMSF
- contact frequency
- docking angle

### Interface and Interaction Analysis

- coarse contact
- hydrogen bond
- salt bridge
- hydrophobic contact
- pi-pi
- cation-pi
- occupancy / persistence

### Biological Interpretation Layer

- chain mapping
- CDR detection
- biological identity
- basic HLA / peptide / TCR annotation
- BSA
- normal mode / PRS
- interface-aware clustering

### Reporting and Delivery

- integrated single-system HTML report
- 3D structure viewer
- section-based organization
- directory-style downloads

## Repository Structure

The repository should be understood in the following way:

```text
immunex/
  core/        execution contracts, context, exceptions, task discovery
  analysis/    domain algorithms and analysis modules
  pipeline/    nodes and workflow orchestration
  cli/         user-facing entry points
  cluster/     HPC / cluster integration
  data/        reference data

docs/
  architecture/  architecture boundaries and engineering rules
  guides/        user and developer guides
  design/        module design notes
  specs/         output and data contracts
  archive/       historical materials
```

## Recommended Documentation Entry Points

Start with:

- [docs/architecture/PRODUCT_SCOPE.md](docs/architecture/PRODUCT_SCOPE.md)
- [docs/architecture/ENGINEERING_ARCHITECTURE.md](docs/architecture/ENGINEERING_ARCHITECTURE.md)
- [docs/architecture/PLATFORM_TODO_ROADMAP.md](docs/architecture/PLATFORM_TODO_ROADMAP.md)
- [docs/guides/IMN_CLI_GUIDE.md](docs/guides/IMN_CLI_GUIDE.md)
- [docs/guides/INSTALLATION_GUIDE.md](docs/guides/INSTALLATION_GUIDE.md)
- [docs/guides/INTERACTION_MODULE_GUIDE.md](docs/guides/INTERACTION_MODULE_GUIDE.md)
- [docs/guides/NORMAL_MODE_MODULE_GUIDE.md](docs/guides/NORMAL_MODE_MODULE_GUIDE.md)
- [docs/design/INTERFACE_CLUSTERING_PIPELINE_DESIGN.md](docs/design/INTERFACE_CLUSTERING_PIPELINE_DESIGN.md)

## Current Development Priorities

The main priority is not to keep adding isolated scripts, but to keep tightening these areas:

- root entry documents and CLI consistency
- report-chain standardization
- stronger cross-module interpretation
- better batch summarization
- cleanup of historical leftovers

The detailed execution roadmap is:

- [docs/architecture/PLATFORM_TODO_ROADMAP.md](docs/architecture/PLATFORM_TODO_ROADMAP.md)

## Notes

- `scripts/` is a support layer, not the primary product entry point
- new capabilities should by default land in:
  - `analysis`
  - `pipeline/nodes`
  - `pipeline`
  - `cli`
- public-facing usage and demos should prefer the `imn` CLI

## License

MIT
