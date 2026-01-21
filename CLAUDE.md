# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Status
**Development**: Local development at `/home/xumy/work/development/AfterMD`
**Testing**: All tests run locally - no remote server deployment needed
**Deployment**: User will manually deploy to server for production testing

## 语言
1.脚本中避免出现中文和emoji
2.对话采取中文对话

## 文件操作规则
1. 及时清理测试脚本，保持项目整洁性
2. markdown只保留对用户必要的，对于某个步骤修改的解释可以放在一个临时文件夹里或者直接以对话的形式告诉我
3. 创建文件时的命名要兼顾在整体程序中的功能和职责，也要保持精简
4. 测试脚本和一次性分析脚本应放在development/目录，不要污染scripts/
5. 重复功能的脚本应该整合，通过参数控制不同策略，而不是创建多个版本

## 文件命名规范
### 处理脚本命名规则
对于PBC处理相关的脚本，按以下格式命名：
- **单任务处理**: `pbc_process.py` - 处理单个MD轨迹的PBC校正
- **批量处理**: `batch_pbc.py` - 本地批量处理多个MD任务
- **集群处理**: `batch_pbc_slurm.py` - 生成SLURM集群批处理脚本

### 命名原则
1. **功能优先**: 体现核心处理步骤(如pbc_process)
2. **规模区分**: single < batch < batch_slurm
3. **执行环境**: 本地处理 vs 集群调度(slurm)
4. **简洁明确**: 避免冗余词汇，保持功能性描述
5. **避免版本后缀**: 不要创建 xxx_v2.py, xxx_new.py, xxx_exact.py 等版本变体

## 项目目录结构

### scripts/ - 核心生产脚本 (18个)
仅保留稳定的、常用的生产环境脚本：

**核心处理:**
- `pbc_process.py` - 单轨迹PBC处理
- `batch_pbc_slurm.py` - SLURM集群批处理
- `batch_worker.py` - 通用SLURM worker

**质量控制:**
- `md_quality_check.py` - MD质量检查
- `md_task_organizer.py` - 任务组织
- `md_workflow.py` - 完整工作流

**批量分析:**
- `batch_allostery_analysis.py` - 变构分析
- `batch_cdr_rmsd_exact.py` - CDR RMSD (整合版本)
- `batch_cdr_rmsf.py` - CDR RMSF
- `batch_contact_frequency.py` - 接触频率
- `batch_interface_rmsd.py` - 界面RMSD
- `batch_phla_analysis.py` - pHLA分析
- `batch_rmsd_hla_alpha.py` - HLA alpha RMSD
- `batch_tcr_rmsd.py` - TCR RMSD
- `batch_whole_protein_rmsf.py` - 全蛋白RMSF

**分析流程:**
- `analyze_chain_contacts.py` - 链接触分析
- `rmsf_analysis_pipeline.py` - RMSF流程
- `run_trajectory_analysis.py` - 轨迹分析

**辅助脚本:**
- `vmd_chain_contacts.tcl` - VMD接触分析

### development/ - 开发和临时脚本
所有测试、一次性分析、实验性脚本都应放在这里：

**development/case_studies/** - 特定案例研究
**development/utilities/** - 临时工具脚本
**development/plotting/** - 可视化实验脚本
**development/clustering/** - 聚类分析实验
**development/visualization/** - 高级可视化
**development/archived_scripts/** - 已废弃脚本

## log板块
对于文件运行的反馈要抓住重点，只在关键的地方去做反馈


## Getting Started

Since this is an empty repository, you'll need to:

1. Initialize the project structure based on the intended technology stack
2. Set up package management files (package.json, requirements.txt, Cargo.toml, etc.)
3. Configure build and development tools
4. Establish testing framework

## Development Commands

### Installation and Setup
```bash
# Install in development mode
pip install -e .

# Install dependencies
pip install -r requirements.txt
```

### Basic Usage
```bash
# Run example script
python examples/basic_usage.py
  - 分层架构设计（Application → Quality Control → Processing → Utility）
  - PBC依赖性分析（Pre-PBC vs Post-PBC模块）
  - 设计模式应用

  2. 模块组织

  - 完整的目录树结构
# Run specific test module
pytest tests/test_trajectory.py
```

### Code Quality
```bash
# Format code
black aftermd/
isort aftermd/

# Lint code
flake8 aftermd/
pylint aftermd/
```

## Architecture

### Project Structure
```
aftermd/
├── __init__.py                 # Main package imports
├── utils/                      # Utility modules (tools)
│   ├── __init__.py
│   ├── batch_processor.py     # BatchProcessor class
│   ├── plotting.py            # PlotManager class  
│   └── path_manager.py        # PathManager class
├── analysis/                   # Analysis modules (core functionality)
│   ├── __init__.py
│   ├── quality/               # Quality control modules
│   │   ├── __init__.py
│   │   ├── md_completeness.py    # MDCompletenessChecker
│   │   ├── structure_validator.py # StructureValidator
│   │   ├── batch_tracker.py      # BatchTracker
│   │   └── quality_reporter.py   # QualityReporter
│   ├── trajectory/            # Trajectory analysis modules
│   │   ├── __init__.py
│   │   ├── rmsd.py           # RMSDCalculator
│   │   ├── rdf.py            # RDFCalculator
│   │   ├── radius_gyration.py # RadiusGyrationCalculator
│   │   ├── distance.py       # DistanceCalculator
│   │   └── hydrogen_bonds.py # HydrogenBondAnalyzer
│   ├── structure/             # Structure analysis modules
│   │   ├── __init__.py
│   │   ├── bfactor.py        # BFactorAnalyzer
│   │   ├── contact_map.py    # ContactMapCalculator
│   │   ├── geometry.py       # GeometryAnalyzer
│   │   └── atom_info.py      # AtomInfoExtractor
│   └── allostery/             # Allostery analysis modules ⭐NEW
│       ├── __init__.py
│       └── contact_correlation.py # ContactCorrelationAnalyzer
└── preprocessing/              # Preprocessing modules
    ├── __init__.py
    └── pbc_processor.py       # PBCProcessor class
```

### Module Organization

**Utils (Utility Modules)**:
- `BatchProcessor`: Parallel processing of multiple MD files
- `PathManager`: File organization and directory structure management
- `PlotManager`: Publication-ready plotting from XVG and analysis data
  - Standard RMSD plotting with statistics and convergence analysis
  - CDR3 RMSD specialized visualizations:
    - `plot_cdr3_rmsd_summary()`: 4-panel comprehensive summary
    - `plot_cdr3_stability_distribution()`: Stability category analysis
    - `plot_cdr3_rmsd_comparison()`: Multi-task trajectory comparison
    - `plot_cdr3_rmsd_heatmap()`: Statistical heatmap across tasks
  - Interactive Plotly support for time series data
  - See: `docs/CDR3_RMSD_PLOTTING_GUIDE.md` for detailed usage
- `GroupSelector`: GROMACS group selection and shortest chain detection
- `ShortestChainDetector`: Detect shortest protein chain from GRO files for PBC centering
- `PDBChainStandardizer`: Standardize PDB chain IDs based on residue count ordering
- `IndexGenerator`: Generate GROMACS index files for pHLA-TCR components
- `CleanupManager`: Temporary file cleanup management

**Analysis/Quality (Quality Control Modules)**:
- `MDCompletenessChecker`: [Pre-PBC] MD simulation completeness verification
- `StructureValidator`: [Pre-PBC] PDB structure quality validation and chain analysis
- `EnergyQualityChecker`: [Pre-PBC] Energy-based quality assessment (NEW)
- `BatchTracker`: Batch processing tracking and duplicate detection
- `QualityReporter`: Comprehensive quality analysis reporting

**Analysis/Trajectory (Trajectory Analysis Modules - Post-PBC)**:
- `RMSDCalculator`: RMSD calculation (MDAnalysis + GROMACS)
- `RDFCalculator`: Radial distribution function analysis
- `RadiusGyrationCalculator`: Radius of gyration and components
- `DistanceCalculator`: Distance measurements (COM, minimum, atom-atom)
- `HydrogenBondAnalyzer`: Hydrogen bond analysis

**Analysis/Structure (Structure Analysis Modules)**:
- `BFactorAnalyzer`: B-factor extraction and analysis by residue/atom
- `ContactMapCalculator`: Contact map calculation and analysis
- `GeometryAnalyzer`: Geometric properties and shape analysis
- `AtomInfoExtractor`: Comprehensive atom information extraction

**Analysis/Allostery (Allostery Analysis Modules)** ⭐NEW:
- `ContactCorrelationAnalyzer`: Dynamic cross-correlation of residue contacts for allosteric analysis

**Preprocessing (Preprocessing Modules)**:
- `PBCProcessor`: PBC artifact removal and trajectory preparation (three-step process)

### Data Flow (Three-Stage Quality Control)

**Stage 1: Pre-PBC Quality Check** (Fast screening, PBC-independent):
1. Raw MD files (md.edr, md.log, md.gro) → File completeness check
2. Energy file (md.edr) → EnergyQualityChecker (temperature, pressure, energy conservation)
3. Structure file (pdb) → StructureValidator (chain count, coordinate validation)
4. Decision: Pass Pre-QC? → Proceed to PBC or Reject

**Stage 2: PBC Processing** (Only for qualified tasks):
5. Qualified MD data → GroupSelector (shortest chain detection)
6. Raw trajectory (md.xtc) → PBCProcessor (center → whole → fit)
7. Output: processed_trajectory.xtc + reference files (md.tpr, md.gro)

**Stage 3: Post-PBC Analysis** (Trajectory-dependent):
8. Processed trajectories → Trajectory analysis (RMSD, Rg, RDF, Distance, H-bonds)
9. PDB structures → Structure analysis (B-factor, Contact map, Geometry)
10. Analysis results → PlotManager (visualization)
11. All quality data → QualityReporter (comprehensive reporting and grading)

**Parallel Processing**:
- BatchProcessor coordinates parallel execution across all stages
- Configurable workers for optimal performance

## Notes

-脚本中避免出现中文和emoji
-对于功能模块我们要综合应用的考虑来设计怎样的绘图可以呈现出直观的结果

## 主要语言
Python 

## 简介
本项目的目的为对gromacs MD production的结果做分析，主要利用MDAnalysis和gromacs内置的程序来进行处理分析

## 初步框架的构建
项目的功能会逐渐的增加，功能主要分为静态的结构的分析和动态的轨迹的分析
## 设计的模块：
### 批量处理模块
我们每个功能都有可能需要对一批文件进行处理，这里需要自动对批量文件处理的能力

### 路径管理模块
我们中间的路径错综复杂，需要集成一个路径管理的模块来专门负责这个部分的内容，比如在处理文件后，结果存放的位置等等

### 绘图模块
一般结果产生的都是xvg之类的文件，我们需要针对这些文件进行可视化的绘图，我们可以集成一些专业的包来处理绘图的问题

## 核心功能部分：
### MD Production质量分析模块

MD Production质量分析模块是AfterMD的核心质量控制功能，用于自动检测MD模拟的完整性、正确性和数据质量。

#### 主要功能

1. **MD完整性检测**
   - 检测MD模拟是否完整运行完成
   - 验证关键输出文件的存在性（md.gro, md.xtc/trr, md.log等）
   - 分析模拟时长和预期时间的匹配度

2. **结构质量验证** (仅限PDB文件)
   - PDB复合物结构分析和链数检测
   - 异常结构识别（链数异常、原子坐标异常等）
   - 体系完整性验证
   - 注意：GRO文件不包含链信息，无法进行链分析

3. **批次追踪和统计**
   - PDB处理次数统计和重复检测
   - 批次间数据对比分析
   - 遗漏数据识别和报告

4. **质量报告生成**
   - 生成详细的质量分析报告
   - 批次统计汇总和可视化
   - 异常数据标记和建议

#### 模块架构

```
analysis/
└── quality/
    ├── __init__.py
    ├── md_completeness.py      # MDCompletenessChecker [Pre-PBC]
    ├── structure_validator.py  # StructureValidator [Pre-PBC]
    ├── energy_quality.py       # EnergyQualityChecker [Pre-PBC] ⭐NEW
    ├── batch_tracker.py        # BatchTracker
    └── quality_reporter.py     # QualityReporter
```

#### 核心类设计

**MDCompletenessChecker** [Pre-PBC]: MD completeness checker
- Detect MD simulation completion status
- Verify output file integrity
- Analyze simulation duration and quality metrics

**StructureValidator** [Pre-PBC]: Structure quality validator (PDB files only)
- PDB file structure analysis and chain count detection
- Protein chain anomaly identification (GRO files lack chain info)
- Complex integrity verification

**EnergyQualityChecker** [Pre-PBC] ⭐NEW: Energy-based quality checker
- Temperature stability check (target ± tolerance)
- Pressure stability check (NPT ensemble)
- Total energy conservation analysis (drift detection)
- Potential/Kinetic energy balance verification
- Comprehensive energy grading (A/B/C/D system)

**BatchTracker**: Batch processing tracker
- PDB processing record management
- Duplicate and missing data detection
- Inter-batch statistics

**QualityReporter**: Quality report generator
- Comprehensive quality analysis reports
- Statistical data visualization
- Issue diagnosis and recommendations

#### 质量检查标准

1. **MD完整性标准**
   - md.gro文件存在且非空
   - 轨迹文件大小合理（> 1MB）
   - 日志文件显示正常结束
   - 模拟时长达到预期

2. **结构质量标准**
   - 标准体系：5条蛋白质链
   - 异常体系：链数 != 5
   - 原子坐标范围合理
   - 缺失残基比例 < 5%

3. **数据完整性标准**
   - 关键文件完整存在
   - 文件大小在合理范围
   - 时间戳逻辑正确

#### 报告输出格式

**质量检查报告**:
```json
{
  "summary": {
    "total_jobs": 150,
    "completed": 145,
    "failed": 3,
    "incomplete": 2,
    "success_rate": 96.7
  },
  "pdb_statistics": {
    "unique_pdbs": 75,
    "multiple_runs": 12,
    "missing_pdbs": 1
  },
  "issues": [
    {
      "pdb_id": "1ABC",
      "issue_type": "incomplete_md",
      "description": "Missing md.gro file",
      "severity": "high"
    }
  ]
}
```

#### 使用示例

```python
from aftermd.analysis.quality import (
    MDCompletenessChecker, StructureValidator,
    EnergyQualityChecker, BatchTracker, QualityReporter
)

# Stage 1: Pre-PBC Quality Screening

# 1.1 MD completeness check
completeness_checker = MDCompletenessChecker(
    min_simulation_time_ps=5000.0
)
completeness_results = completeness_checker.batch_check("/path/to/md/results")

# 1.2 Energy quality check (NEW - Fast screening)
energy_checker = EnergyQualityChecker(
    target_temperature=300.0,
    temp_tolerance=5.0,
    energy_drift_threshold=2.0
)
energy_result = energy_checker.comprehensive_energy_check("/path/to/md.edr")

# Decision: Only proceed with Grade A/B/C
if energy_result['energy_grade'] in ['A', 'B', 'C']:
    print(f"Qualified for PBC processing: Grade {energy_result['energy_grade']}")
else:
    print(f"Rejected: Poor energy quality (Grade {energy_result['energy_grade']})")

# 1.3 PDB structure validation (PDB files only)
structure_validator = StructureValidator(expected_chain_count=5)
pdb_files = ["/path/to/structure1.pdb", "/path/to/structure2.pdb"]
validation_results = structure_validator.batch_validate(pdb_files)

# Batch tracking
batch_tracker = BatchTracker()
batch_data = batch_tracker.analyze_batch("/path/to/batch")

# Comprehensive quality report
reporter = QualityReporter()
comprehensive_report = reporter.generate_comprehensive_report(
    completeness_results, validation_results, batch_data
)
```

#### 可能面临的问题及解决方案

1. **PDB结构包含多个复合物**
   - **问题**: 复合物链数异常，正常体系为5条链
   - **检测**: 通过MDAnalysis分析PDB文件链数
   - **处理**: 标记异常结构，生成警告报告

2. **MD模拟未完整运行**  
   - **问题**: MD因异常中断，输出不完整
   - **检测**: 检查md.gro存在性和日志文件状态
   - **处理**: 标记为失败任务，建议重新运行

3. **批次数据重复和遗漏**
   - **问题**: 同一PDB多次处理或某些PDB被遗漏
   - **检测**: 基于PDB ID（文件名前4字符）统计
   - **处理**: 生成重复/遗漏报告，建议数据清理

### Preprocessing Workflow

**PBC Artifact Removal** - Three-step process:
1. **Center**: Center trajectory on shortest chain (peptide) using `gmx trjconv -center -pbc atom`
2. **Whole**: Make molecules whole using `gmx trjconv -pbc whole`
3. **Fit**: Fit rotational and translational motion using `gmx trjconv -fit rot+trans`

**Key Features**:
- Automatic shortest chain detection using `gmx make_ndx -splitch`
- Strict validation (no fallback, ensures scientific correctness)
- Reference file copying (md.tpr, md.gro) for downstream analysis

### Trajectory Analysis (Post-PBC)

**Available Modules**:
1. RMSD calculation and plotting
2. Radius of gyration analysis
3. Radial distribution function (RDF)
4. Distance measurements (COM, minimum, atom-atom)
5. Hydrogen bond analysis

### Structure Analysis

**Available Modules**:
1. B-factor extraction and visualization
2. Contact map calculation
3. Geometric property analysis
4. Atom information extraction

### Quality Control Strategy

**Three-Stage Approach**:
- **Stage 1** [Pre-PBC]: Fast screening (file check, energy quality, structure validation)
- **Stage 2** [PBC Processing]: Only for qualified tasks
- **Stage 3** [Post-PBC]: Trajectory-based validation (RMSD, Rg, convergence)

**Grading System**: A (90-100) / B (75-89) / C (60-74) / D (<60)

See `ARCHITECTURE.md` for detailed workflow diagrams and module specifications.

## 新增模块 (2025-12)

### Analysis/Angles (角度分析模块)

专门用于分子动力学轨迹中的角度相关分析，特别针对TCR-pMHC复合物对接姿态研究。

**核心类**:
- `PrincipalAxesCalculator`: 主轴计算基础类，基于惯性张量对角化
- `DockingAngleAnalyzer`: TCR-pMHC对接角度分析 (Twist/Tilt/Swing三种角度)
- `DihedralCalculator`: 通用二面角计算器
- `VectorAngleCalculator`: 向量夹角计算工具类

**模块架构**:
```
analysis/angles/
├── __init__.py
├── principal_axes.py      # 主轴计算基础类 (~150行)
├── docking_angles.py      # TCR-pMHC对接角度 (~250行)
├── dihedral.py            # 通用二面角 (~100行)
└── vector_angles.py       # 向量夹角 (~100行)
```

**主要功能**:

1. **主轴计算** (principal_axes.py)
   - 基于惯性张量的主轴计算 (复用GeometryAnalyzer方法)
   - 支持单帧和整条轨迹的主轴演化
   - 返回特征值和特征向量 (按惯性矩降序排列)

2. **TCR-pMHC对接角度** (docking_angles.py)
   - **Twist角**: TCR二硫键连线在MHC平面上投影与MHC主轴的夹角
   - **Tilt角**: TCR主轴在MHC平面上投影与MHC主轴的夹角
   - **Swing角**: TCR连线相对peptide的侧向偏移角度
   - 源自VMD脚本，完全Python实现

3. **二面角计算** (dihedral.py)
   - 支持任意4个原子定义的二面角
   - 主链φ/ψ角和侧链旋转角
   - 返回-180到180度范围

4. **向量夹角** (vector_angles.py)
   - 主轴-主轴夹角
   - 质心连线-参考方向夹角
   - 灵活的向量角度计算接口

**使用示例**:
```python
from aftermd.analysis.angles import DockingAngleAnalyzer

analyzer = DockingAngleAnalyzer("md.tpr", "md_pbc.xtc")
times, twist, tilt, swing = analyzer.calculate_docking_angles_trajectory(
    mhc_selection="segname PROA and resid 50:86 140:176",
    tcr_alpha_cys="segname PROD and resid 89:94 20:25 and resname CYS and name CA",
    tcr_beta_cys="segname PROE and resid 89:94 20:25 and resname CYS and name CA",
    tcr_v_selection="segname PROD PROE and resid 3:115",
    peptide_selection="segname PROC",
    output_file="docking_angles.csv"
)
```

**技术特点**:
- 纯Python实现，基于MDAnalysis和NumPy
- 复用现有GeometryAnalyzer的惯性张量计算 (geometry.py:55-68)
- 向量投影和几何计算优化
- 输出XVG/CSV格式，兼容GROMACS工具链

---

### Analysis/Interface (界面分析模块)

用于分析蛋白质-蛋白质相互作用界面的物理化学性质。

**核心类**:
- `BuriedSurfaceCalculator`: 埋藏表面积计算
- `ContactResidueAnalyzer`: 接触残基对统计
- `ContactAtomAnalyzer`: 接触原子详细信息
- `ContactHeatmapGenerator`: 界面接触热力图

**模块架构**:
```
analysis/interface/
├── __init__.py
├── buried_surface.py      # 埋藏表面积 (~120行)
├── contact_residues.py    # 接触残基统计 (~150行)
├── contact_atoms.py       # 接触原子详情 (~130行)
└── contact_heatmap.py     # 界面热力图 (~120行)
```

**主要功能**:

1. **埋藏表面积 (BSA)** (buried_surface.py)
   - 公式: BSA = (SASA_A + SASA_B - SASA_AB) / 2
   - 使用MDAnalysis.analysis.sasa模块
   - 探针半径1.4 Å (水分子)
   - 支持轨迹时间演化

2. **接触残基对统计** (contact_residues.py)
   - 基于距离截断 (默认4.0 Å) 识别接触
   - 统计接触频率 (接触帧数/总帧数)
   - 输出DataFrame: ResA_ID, ResA_Name, ResB_ID, ResB_Name, Frequency

3. **接触原子详细信息** (contact_atoms.py)
   - 逐帧输出所有接触原子对
   - 包含原子ID、名称、残基信息、距离
   - 使用距离矩阵加速计算
   - 支持stride参数减少计算量

4. **界面热力图** (contact_heatmap.py)
   - 生成残基-残基接触频率矩阵
   - 使用seaborn绘制热力图
   - 输出300 dpi高分辨率图像
   - 集成PlotManager样式

**使用示例**:
```python
from aftermd.analysis.interface import (
    BuriedSurfaceCalculator,
    ContactResidueAnalyzer,
    ContactHeatmapGenerator
)

# BSA分析
bsa_calc = BuriedSurfaceCalculator("md.tpr", "md.xtc")
times, bsa_values = bsa_calc.calculate_bsa_trajectory(
    selection_a="segname PROA PROB PROC",  # pMHC
    selection_b="segname PROD PROE",       # TCR
    output_file="bsa_trajectory.csv"
)

# 接触频率统计
contact_analyzer = ContactResidueAnalyzer("md.tpr", "md.xtc")
contact_df = contact_analyzer.calculate_contact_frequency(
    selection_a="segname PROD",  # TCR alpha
    selection_b="segname PROC",  # Peptide
    cutoff=4.0
)

# 热力图生成
heatmap_gen = ContactHeatmapGenerator("md.tpr", "md.xtc")
contact_matrix = heatmap_gen.generate_and_plot(
    selection_a="segname PROD",
    selection_b="segname PROC",
    output_file="contact_heatmap.png"
)
```

**性能优化**:
- 距离矩阵向量化计算 (MDAnalysis.lib.distances.distance_array)
- 轨迹分块处理 (每次1000帧)
- 流式写入大文件
- 并行计算多任务

---

### Utils/BatchAnalyzer (统一批量分析工具)

集成角度分析和界面分析的统一批量处理接口。

**文件**: `aftermd/utils/batch_analyzer.py` (~400行)

**核心类**:
```python
class BatchAnalyzer:
    def batch_docking_angles()        # 批量对接角度分析
    def batch_bsa_analysis()          # 批量BSA分析
    def batch_contact_analysis()      # 批量接触分析
    def comprehensive_analysis()      # 综合批量分析
```

**主要功能**:

1. **批量角度分析**
   - 并行处理多个MD任务的对接角度
   - 自动创建任务输出目录
   - 错误隔离：单任务失败不影响其他

2. **批量界面分析**
   - 支持frequency/atoms/heatmap三种模式
   - 统一的任务配置格式
   - 生成综合分析报告

3. **综合分析**
   - 灵活的模块组合 (docking_angles, bsa, contact_frequency, contact_heatmap)
   - 基于ProcessPoolExecutor的多进程并行
   - 详细的处理日志和错误报告

**使用示例**:
```python
from aftermd.utils import BatchAnalyzer

tasks = [
    {
        'name': '1AO7',
        'topology': 'data/1AO7/md.tpr',
        'trajectory': 'data/1AO7/md_pbc.xtc',
        'mhc_selection': '...',
        'tcr_alpha_cys': '...',
        # ... 其他选择参数
    },
    # 更多任务...
]

batch_analyzer = BatchAnalyzer()
results = batch_analyzer.comprehensive_analysis(
    task_list=tasks,
    analysis_modules=['docking_angles', 'bsa', 'contact_frequency', 'contact_heatmap'],
    max_workers=4,
    output_dir="./results"
)
```

**设计特点**:
- 基于现有BatchProcessor模式 (batch_processor.py)
- 完善的错误处理和日志记录
- 支持YAML/JSON配置文件
- 生成汇总统计报告

**输出目录结构**:
```
results/
└── task_name/
    ├── angles/
    │   ├── docking_angles.csv
    │   ├── twist_angle.xvg
    │   ├── tilt_angle.xvg
    │   └── swing_angle.xvg
    ├── interface/
    │   ├── bsa_trajectory.csv
    │   ├── contact_frequency.csv
    │   ├── contact_atoms.csv
    │   └── contact_heatmap.png
    └── plots/
        ├── docking_angles.png
        └── bsa_evolution.png
```

---

### 模块依赖关系

```
principal_axes.py (基础)
    ├── docking_angles.py (继承)
    ├── vector_angles.py (继承)
    └── 独立使用

dihedral.py (独立)

buried_surface.py (独立)
contact_residues.py (独立)
contact_atoms.py (独立)
contact_heatmap.py (依赖contact_residues算法)

batch_analyzer.py (调用所有上述模块)
```

### 数据流扩展

```
Processed trajectories
    ├── Angle Analysis
    │   ├── Principal axes calculation
    │   ├── Docking angles (Twist/Tilt/Swing)
    │   ├── Dihedral angles
    │   └── Vector angles
    ├── Interface Analysis
    │   ├── Buried surface area (BSA)
    │   ├── Contact residue statistics
    │   ├── Contact atom details
    │   └── Contact heatmap visualization
    └── Allostery Analysis ⭐NEW
        ├── Contact time series matrix (N_contacts × T_frames)
        ├── Correlation matrix calculation (Pearson)
        └── Cooperative motion visualization (red-blue heatmap)
            ↓
    PlotManager (统一可视化)
            ↓
    Batch analysis reports
```

---

### Analysis/Allostery (变构分析模块) ⭐NEW 2026-01

用于研究蛋白质内部的协同运动和变构通讯路径。

**核心类**:
- `ContactCorrelationAnalyzer`: 残基接触动态相关性分析

**模块架构**:
```
analysis/allostery/
├── __init__.py
└── contact_correlation.py     # ContactCorrelationAnalyzer (~500行)
```

**主要功能**:

**动态接触相关性分析** (contact_correlation.py)
- 遍历轨迹计算残基对之间的最短重原子距离
- 构建接触时间序列二值矩阵 (N_contacts × T_frames)
- 过滤条件:
  - 序列相邻残基 (|i-j| < 3)
  - 稀有接触 (频率 < 10%)
- 计算接触对之间的皮尔逊相关系数 (N_contacts × N_contacts)
- 识别协同运动 (正相关) 和反协同运动 (负相关)

**使用示例**:
```python
from aftermd.analysis.allostery import ContactCorrelationAnalyzer

# 初始化分析器
analyzer = ContactCorrelationAnalyzer("md.tpr", "md_pbc.xtc")

# 运行完整分析
results = analyzer.run_full_analysis(
    selection="protein",           # 分析整个蛋白
    output_dir="./allostery_analysis",
    cutoff=4.5,                    # 接触距离阈值 (Angstrom)
    seq_dist_cutoff=3,             # 排除 |i-j| < 3
    min_frequency=0.1,             # 仅保留 >= 10% 帧数中的接触
    plot_heatmap=True              # 生成红蓝热力图
)

# 查看结果
print(f"持久接触对数量: {results['n_contacts']}")
print(f"高度正相关对 (>0.7): {results['n_high_pos_correlation']}")
print(f"高度负相关对 (<-0.7): {results['n_high_neg_correlation']}")
```

**输出文件**:
```
allostery_analysis/
├── contact_labels.csv          # 接触对定义和频率
├── correlation_matrix.csv      # 相关系数矩阵 (CSV)
├── correlation_matrix.npz      # 相关系数矩阵 (NumPy压缩)
├── contact_time_series.npz     # 接触时间序列 (N_contacts × T_frames)
└── correlation_heatmap.png     # 红蓝热力图 (红=正相关, 蓝=负相关)
```

**技术特点**:
- 两遍算法：第一遍识别持久接触，第二遍构建二值矩阵 (内存高效)
- 重原子过滤 (`not name H*`)
- 向量化距离计算 (MDAnalysis.lib.distances)
- 红蓝发散色图 (`RdBu_r`): 红色=协同运动，蓝色=反协同运动
- 基于NumPy的高效皮尔逊相关计算

**应用场景**:
- 变构通讯路径识别
- 协同运动区域检测
- 功能域动态耦合分析
- 药物设计中的变构位点预测

**性能优化**:
- 使用 `uint8` 存储二值矩阵 (节省内存)
- 字典索引加速第二遍扫描
- 支持大规模体系 (300残基, 10000+帧)

---

### 参考来源

这些模块的算法源自：
- VMD脚本: `development/angle_analysis/Analysis-scripts/`
  - mkvmd_twistangle.sh
  - mkvmd_tiltangle.sh
  - mkvmd_swingangle.sh
  - mkvmd_buriedarea.sh
  - mkvmd_PairContact.sh

完全用Python/MDAnalysis重新实现，无需VMD依赖。

