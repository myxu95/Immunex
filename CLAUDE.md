# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Status
**Development**: Local development at `/home/xumy/work/development/Immunex`
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
6. **禁止创建临时markdown报告文件**：如无必要，不要写markdown文件作为短期任务汇报（如UPDATE_2026-xx-xx.md、SUMMARY_xxx.md等），这会让项目变得混乱。应该直接在对话中汇报结果，或将重要信息更新到已有文档中

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
black immunex/
isort immunex/

# Lint code
flake8 immunex/
pylint immunex/
```

## Architecture

### 四层架构设计 (Layered Architecture)

Immunex 采用分层架构设计，确保模块职责清晰、可复用性高、易于维护。

#### 第一层：Core Modules（原子功能模块）

**定义**: 最小可复用分析单元，提供单一、独立的功能。

**核心模块示例**:
- `preprocess_trajectory` - 轨迹预处理（PBC校正）
- `compute_rmsd` - RMSD计算
- `compute_rmsf` - RMSF计算
- `compute_contacts` - 接触分析
- `compute_distance` - 距离测量
- `cluster_trajectory` - 轨迹聚类
- `extract_frames` - 帧提取

**设计原则**:
1. **单一职责**: 每个模块只做一件事，并做好
2. **输入输出标准化**: 使用 dataclass 定义输入/输出结构
3. **上下文无关**: 不关心上游是谁、下游是谁、是否批处理、是否来自REST API

**示例**: `compute_rmsd` 模块
```python
class RMSDCalculator:
    def calculate(self,
                  topology: str,
                  trajectory: str,
                  selection: str,
                  reference: Optional[str] = None) -> RMSDResult:
        """
        职责：计算RMSD

        它不该关心：
        - 是不是 batch 跑的
        - 是不是来自 REST2 模拟
        - 后面要不要画图
        - 是不是 agent 调用的

        只关心：
        - 接收处理好的轨迹/拓扑/selection
        - 计算 RMSD
        - 返回结果文件路径、统计值、元数据
        """
        pass
```

---

#### 第二层：Pipeline Steps / Workflow Nodes（工作流节点）

**定义**: 对原子模块的轻量封装，负责流程控制而非实际分析。

**节点示例**:
- `PreprocessNode` - 预处理节点
- `RMSDNode` - RMSD计算节点
- `ContactNode` - 接触分析节点
- `QualityCheckNode` - 质量检查节点

**节点职责**:
1. **校验输入上下文**: 检查必需的输入是否存在
2. **调用底层模块**: 委托给 Core Module 执行实际工作
3. **记录输出到 context**: 将结果写入 pipeline context
4. **处理异常**: 失败处理、跳过逻辑、缓存机制

**示例**: RMSDNode
```python
class RMSDNode:
    """RMSD计算节点 - 不亲自做分析，而是编排"""

    def execute(self, context: PipelineContext) -> PipelineContext:
        # 1. 校验输入
        if 'trajectory_processed' not in context:
            raise ValueError("Missing processed trajectory")

        # 2. 调用底层模块
        calculator = RMSDCalculator()
        result = calculator.calculate(
            topology=context['topology'],
            trajectory=context['trajectory_processed'],
            selection=context['selections']['protein']
        )

        # 3. 记录输出
        context['results']['rmsd'] = {
            'output_file': result.output_file,
            'mean': result.mean_rmsd,
            'std': result.std_rmsd
        }

        # 4. 处理失败（可选）
        if not result.success:
            context['errors'].append(f"RMSD failed: {result.error_message}")

        return context
```

---

#### 第三层：Workflow / Pipeline Orchestration（编排层）

**定义**: 决定节点执行顺序、数据流动、依赖关系和并行策略。

**编排职责**:
1. **定义执行顺序**: 哪些节点按什么顺序执行
2. **数据流管理**: 节点之间的数据如何流动
3. **失败策略**: 某个节点失败后怎么办
4. **并行调度**: 哪些节点可以并行执行
5. **依赖解析**: 哪些节点依赖前置结果

**线性流程示例**:
```python
class StandardTrajectoryPipeline:
    """标准轨迹分析流程"""

    def __init__(self):
        self.nodes = [
            LoadInputNode(),
            PreprocessNode(),
            AlignNode(),
            RMSDNode(),
            RMSFNode(),
            ContactNode(),
            ReportNode()
        ]

    def execute(self, context: PipelineContext) -> PipelineContext:
        for node in self.nodes:
            context = node.execute(context)
            if context.should_stop:
                break
        return context
```

**分叉流程示例**:
```python
class AdvancedAnalysisPipeline:
    """高级分析流程 - 支持分叉和并行"""

    def execute(self, context: PipelineContext) -> PipelineContext:
        # 预处理（必需）
        context = PreprocessNode().execute(context)

        # 并行执行多个分析
        parallel_nodes = [
            RMSDNode(),
            RMSFNode(),
            ClusteringNode(),
            ContactMapNode()
        ]

        with ProcessPoolExecutor(max_workers=4) as executor:
            futures = [
                executor.submit(node.execute, context.copy())
                for node in parallel_nodes
            ]
            results = [f.result() for f in futures]

        # 合并结果
        context = self.merge_results(results)
        return context
```

**动态流程示例**:
```python
class AdaptivePipeline:
    """自适应流程 - 根据质量检查结果动态决定"""

    def execute(self, context: PipelineContext) -> PipelineContext:
        # 质量检查
        context = QualityCheckNode().execute(context)

        # 根据质量等级选择不同流程
        if context['quality_grade'] == 'A':
            # 高质量：完整分析
            return self.full_analysis(context)
        elif context['quality_grade'] in ['B', 'C']:
            # 中等质量：基础分析
            return self.basic_analysis(context)
        else:
            # 低质量：仅报告问题
            return self.report_issues(context)
```

**关键设计原则**:
- 不要把流程逻辑写死在每个模块里
- 流程配置应该独立存在（代码或配置文件）
- 支持流程的动态组合和重用

---

#### 第四层：Configuration / User Interface（配置与入口层）

**定义**: 用户如何驱动 pipeline，提供多种入口方式。

**入口方式**:

1. **Python API** (编程接口)
```python
from immunex.pipeline import StandardTrajectoryPipeline
from immunex.core import PipelineContext

context = PipelineContext(
    system_id="1ao7",
    topology="data/1ao7/md.tpr",
    trajectory_raw="data/1ao7/md.xtc"
)

pipeline = StandardTrajectoryPipeline()
result = pipeline.execute(context)
```

2. **CLI** (命令行接口)
```bash
immunex run-pipeline \
    --config config.yaml \
    --system 1ao7 \
    --topology data/1ao7/md.tpr \
    --trajectory data/1ao7/md.xtc
```

3. **YAML Config** (配置文件驱动)
```yaml
# config.yaml
pipeline: standard_trajectory

input:
  system_id: 1ao7
  topology: data/1ao7/md.tpr
  trajectory_raw: data/1ao7/md.xtc

nodes:
  - name: preprocess
    enabled: true
  - name: rmsd
    params:
      selection: "protein"
  - name: rmsf
    params:
      selection: "backbone"

output:
  directory: ./results/1ao7
```

4. **Web UI** (Web界面，未来扩展)
```
用户通过网页上传文件、选择分析模块、查看结果
```

5. **Agent调用** (AI Agent，未来扩展)
```python
# Agent可以根据用户自然语言指令动态构建pipeline
agent.run("分析1ao7的RMSD和RMSF，如果RMSD不收敛就做聚类")
```

**关键原则**:
- 这一层不要和底层分析逻辑耦合
- 每种入口方式都是对 Pipeline 的薄封装
- 配置验证在这一层完成

---

### 上下文传递机制 (Context Passing)

**问题**: 避免"文件路径地狱"

**错误做法**:
```python
# 模块 A 输出 xxx_processed.xtc
# 模块 B 手工去读这个路径
# 模块 C 再拼字符串找 xxx_rmsd.xvg
# 模块 D 再假设文件名规则
# 结果：路径硬编码、命名冲突、维护困难
```

**正确做法**: 使用 **PipelineContext** 统一管理数据

#### PipelineContext 设计

```python
from dataclasses import dataclass, field
from typing import Dict, Any, Optional
from pathlib import Path

@dataclass
class PipelineContext:
    """Pipeline上下文 - 统一管理所有数据和状态"""

    # 系统标识
    system_id: str

    # 输入文件
    topology: str
    trajectory_raw: str
    structure_pdb: Optional[str] = None

    # 处理后的文件
    trajectory_processed: Optional[str] = None
    trajectory_aligned: Optional[str] = None

    # 原子选择
    selections: Dict[str, str] = field(default_factory=dict)
    # 示例: {'protein': 'protein', 'backbone': 'name CA C N O'}

    # 分析结果
    results: Dict[str, Any] = field(default_factory=dict)
    # 示例: {'rmsd': {...}, 'rmsf': {...}}

    # 元数据
    metadata: Dict[str, Any] = field(default_factory=dict)
    # 示例: {'timestamp': '2026-03-16', 'version': '1.0.0'}

    # 输出目录
    output_dir: Optional[str] = None

    # 错误和警告
    errors: list = field(default_factory=list)
    warnings: list = field(default_factory=list)

    # 流程控制
    should_stop: bool = False

    def get_output_path(self, filename: str) -> str:
        """获取输出文件路径"""
        if self.output_dir is None:
            self.output_dir = f"./results/{self.system_id}"
        Path(self.output_dir).mkdir(parents=True, exist_ok=True)
        return str(Path(self.output_dir) / filename)
```

#### 使用示例

```python
# 初始化 context
context = PipelineContext(
    system_id="1ao7",
    topology="data/1ao7/md.tpr",
    trajectory_raw="data/1ao7/md.xtc",
    selections={
        'protein': 'protein',
        'backbone': 'name CA C N O',
        'tcr': 'chainID D E'
    }
)

# 预处理节点
preprocessor = PreprocessNode()
context = preprocessor.execute(context)
# context.trajectory_processed = "./results/1ao7/md_processed.xtc"

# RMSD节点（无需手动指定路径）
rmsd_node = RMSDNode()
context = rmsd_node.execute(context)
# context.results['rmsd'] = {
#     'output_file': './results/1ao7/rmsd.xvg',
#     'mean': 2.34,
#     'std': 0.56
# }

# RMSF节点（自动从context获取所需数据）
rmsf_node = RMSFNode()
context = rmsf_node.execute(context)
# context.results['rmsf'] = {...}

# 访问结果
print(f"RMSD mean: {context.results['rmsd']['mean']} nm")
print(f"Output: {context.results['rmsd']['output_file']}")
```

#### 优势总结

1. **无路径硬编码**: 所有路径由 context 统一管理
2. **类型安全**: 使用 dataclass，IDE 可以自动补全
3. **易于测试**: 可以构造 mock context 进行单元测试
4. **状态透明**: 所有中间结果都在 context 中可见
5. **易于序列化**: 可以保存/恢复 context 状态
6. **错误追踪**: 统一收集错误和警告信息

---

### Project Structure
```
immunex/
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

---

## 模块设计标准

所有核心功能模块（如PBCProcessor、RMSDCalculator等）必须遵循以下6大设计原则，确保模块的健壮性和可维护性：

### 1️⃣ 输入明确 (Clear Inputs)

**原则**:
- 所有参数必须有类型注解 (type hints)
- 必需参数和可选参数明确区分
- 提供输入验证机制
- 参数有合理的默认值

**推荐实现**:
```python
from dataclasses import dataclass
from typing import Optional
from pathlib import Path

@dataclass
class ProcessingInput:
    """标准化输入参数"""
    # 必需参数
    trajectory: str
    topology: str
    output: str

    # 可选参数（有默认值）
    method: str = "2step"
    dt: Optional[float] = None
    verbose: bool = False

    def validate(self) -> None:
        """验证输入参数"""
        if not Path(self.trajectory).exists():
            raise FileNotFoundError(f"Trajectory not found: {self.trajectory}")
        if not Path(self.topology).exists():
            raise FileNotFoundError(f"Topology not found: {self.topology}")
        if self.method not in ["2step", "3step"]:
            raise ValueError(f"Invalid method: {self.method}")
```

### 2️⃣ 输出明确 (Clear Outputs)

**原则**:
- 返回结构化数据而非简单字符串/None
- 包含处理统计信息
- 包含所有生成的文件路径
- 包含处理元数据（时间、版本等）
- 包含成功/失败状态

**推荐实现**:
```python
@dataclass
class ProcessingResult:
    """标准化输出结果"""
    success: bool                          # 是否成功
    output_file: str                       # 主输出文件
    temporary_files: List[str] = None      # 临时文件列表

    # 统计信息
    processing_stats: Dict[str, Any] = None
    # 示例: {'n_frames': 1000, 'processing_time_sec': 45.2}

    # 元数据
    metadata: Dict[str, Any] = None
    # 示例: {'timestamp': '2026-03-16', 'version': '0.1.0'}

    # 错误信息（如果失败）
    error_message: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        """支持序列化"""
        return asdict(self)
```

### 3️⃣ 副作用明确 (Clear Side Effects)

**原则**:
- 明确记录所有文件操作（创建、修改、删除）
- 明确记录所有外部命令调用
- 提供副作用报告
- 支持dry-run模式（可选）

**推荐实现**:
```python
@dataclass
class SideEffectTracker:
    """副作用跟踪器"""
    files_created: List[str] = field(default_factory=list)
    files_deleted: List[str] = field(default_factory=list)
    commands_executed: List[Dict] = field(default_factory=list)

    def track_file_creation(self, filepath: str):
        self.files_created.append(filepath)

    def track_command(self, command: List[str]):
        self.commands_executed.append({
            'command': ' '.join(command),
            'timestamp': datetime.now().isoformat()
        })
```

### 4️⃣ 错误明确 (Clear Errors)

**原则**:
- 定义具体的异常类型
- 提供详细的错误信息和上下文
- 区分用户错误和系统错误
- 提供错误恢复建议

**推荐实现**:
```python
class ImmunexError(Exception):
    """Immunex基础异常"""
    pass

class InputValidationError(ImmunexError):
    """输入验证错误（用户错误）"""
    def __init__(self, param_name: str, value: Any, reason: str):
        super().__init__(f"Invalid {param_name}: {reason} (got: {value})")

class ProcessingError(ImmunexError):
    """处理错误（系统错误）"""
    def __init__(self, step: str, reason: str, suggestion: str = None):
        msg = f"Failed at {step}: {reason}"
        if suggestion:
            msg += f"\nSuggestion: {suggestion}"
        super().__init__(msg)
```

**错误处理模式**:
```python
def process(self, input_params: ProcessingInput) -> ProcessingResult:
    try:
        input_params.validate()  # 可能抛出InputValidationError
        result = self._do_processing(...)
        return ProcessingResult(success=True, output_file=result)

    except InputValidationError as e:
        # 用户错误 - 返回失败结果
        return ProcessingResult(
            success=False,
            error_message=str(e)
        )
    except Exception as e:
        # 系统错误 - 记录日志并返回失败结果
        logger.exception("Unexpected error")
        return ProcessingResult(
            success=False,
            error_message=f"System error: {str(e)}"
        )
```

### 5️⃣ 可独立测试 (Testable)

**原则**:
- 每个模块必须有单元测试
- 测试覆盖率 >= 80%
- 测试输入验证、成功场景、失败场景
- 支持mock外部依赖（GROMACS等）

**推荐测试结构**:
```python
class TestProcessor:
    @pytest.fixture
    def processor(self):
        return Processor()

    def test_input_validation_missing_file(self):
        """测试：缺失文件应抛出异常"""
        with pytest.raises(FileNotFoundError):
            input_params.validate()

    def test_processing_success(self, monkeypatch):
        """测试：成功处理场景"""
        # Mock外部依赖
        monkeypatch.setattr(...)
        result = processor.process(valid_input)
        assert result.success is True

    def test_processing_failure(self):
        """测试：失败场景处理"""
        result = processor.process(invalid_input)
        assert result.success is False
        assert result.error_message is not None
```

### 6️⃣ 可被调度 (Schedulable)

**原则**:
- 支持异步/并行调用
- 提供进度回调机制
- 支持取消操作
- 可被pipeline统一调度

**推荐实现**:
```python
class SchedulableProcessor:
    """可调度处理器基类"""

    def __init__(self):
        self._cancel_flag = threading.Event()
        self.progress_callback: Optional[Callable] = None

    def set_progress_callback(self, callback: Callable[[float, str], None]):
        """设置进度回调"""
        self.progress_callback = callback

    def _report_progress(self, progress: float, message: str):
        """报告进度 (0.0-1.0)"""
        if self.progress_callback:
            self.progress_callback(progress, message)

    def cancel(self):
        """取消操作"""
        self._cancel_flag.set()

    def is_cancelled(self) -> bool:
        return self._cancel_flag.is_set()

    def process(self, input_params: Any) -> Any:
        """处理入口 - 支持进度报告和取消"""
        if self.is_cancelled():
            return ProcessingResult(success=False, error_message="Cancelled")

        self._report_progress(0.0, "Starting...")
        # ... 处理逻辑 ...
        self._report_progress(0.5, "Half done...")
        # ... 继续处理 ...
        self._report_progress(1.0, "Completed")
```

**Pipeline调度示例**:
```python
from concurrent.futures import ProcessPoolExecutor

class PipelineScheduler:
    def submit_task(self, task_id: str, processor: SchedulableProcessor,
                   input_params: Any) -> Future:
        """提交任务到调度器"""
        with ProcessPoolExecutor() as executor:
            future = executor.submit(processor.process, input_params)
            return future
```

---

### 设计检查清单

在开发新模块或重构现有模块时，使用此清单确保符合标准：

- [ ] **输入明确**
  - [ ] 使用dataclass定义输入结构
  - [ ] 所有参数有类型注解
  - [ ] 实现validate()方法
  - [ ] 有合理的默认值

- [ ] **输出明确**
  - [ ] 使用dataclass定义输出结构
  - [ ] 包含success字段
  - [ ] 包含processing_stats和metadata
  - [ ] 支持to_dict()序列化

- [ ] **副作用明确**
  - [ ] 记录文件操作
  - [ ] 记录命令执行
  - [ ] 文档中列出所有副作用

- [ ] **错误明确**
  - [ ] 定义专用异常类
  - [ ] 错误信息包含上下文和建议
  - [ ] 返回结构化错误信息

- [ ] **可独立测试**
  - [ ] 单元测试覆盖率 >= 80%
  - [ ] 测试输入验证
  - [ ] 测试成功和失败场景
  - [ ] 支持mock外部依赖

- [ ] **可被调度**
  - [ ] 支持进度回调
  - [ ] 支持取消操作
  - [ ] 可异步/并行调用

---

### 适用范围

本标准适用于所有Immunex核心模块，包括：
- PBCProcessor
- RMSDCalculator
- QualityAssessmentPipeline
- ContactAnalyzer
- AngleAnalyzer
- 以及未来新增的所有处理模块

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

MD Production质量分析模块是Immunex的核心质量控制功能，用于自动检测MD模拟的完整性、正确性和数据质量。

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
from immunex.analysis.quality import (
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

### Analysis/Angles (角度分析模块) ✅ **已完成 (2026-01-23)**

**状态**: Phase 1 核心功能已实现并测试通过

**实施总结**: `development/reports/ANGLES_IMPLEMENTATION_SUMMARY.md`
**设计文档**: `docs/design/angles/ANGLE_MODULE_REDESIGN.md`
**归档代码**: `development/archived_scripts/README_ANGLES_ARCHIVE.md`

专门用于分子动力学轨迹中的角度相关分析，特别针对TCR-pMHC复合物对接姿态研究。

**已实现的核心类**:
- `PrincipalAxesCalculator`: 主轴计算基础类，基于惯性张量对角化和PCA
- `DockingAngleAnalyzer`: TCR-pMHC对接角度分析 (Twist/Tilt/Swing三种角度)
- **工具函数**: `angle_between_vectors`, `project_vector_onto_plane`, `signed_angle_between_vectors`, `dihedral_angle`

**未来扩展** (可选):
- `TrajectoryAngleAnalyzer`: 高层轨迹分析接口 (优先级低)
- 角度聚类和分布分析
- 批处理脚本和CLI命令

**模块架构**:
```
analysis/angles/
├── __init__.py              # ✅ 已完成 - 导入所有模块
├── principal_axes.py        # ✅ 已实现 - 主轴计算 (~150行)
├── docking_angles.py        # ✅ 已实现 - TCR-pMHC对接角度 (~380行)
├── vector_angles.py         # ✅ 已实现 - 向量夹角工具 (~140行)
└── trajectory_angles.py     # ⏳ 待实现 - 轨迹分析封装 (可选)
```

**设计目标** (✅ 已实现):

1. **纯Python实现**
   - ✅ 基于MDAnalysis + NumPy
   - ✅ 不依赖外部C++工具
   - ✅ 易于维护和扩展

2. **灵活的原子选择**
   - ✅ 支持MDAnalysis选择语法
   - ✅ 不限于固定链命名 (A/B/C/D/E)
   - ✅ 适用于各种蛋白复合物

3. **完整的轨迹支持**
   - ✅ 单结构角度计算
   - ✅ MD轨迹时间演化
   - ✅ 角度统计和稳定性分析

4. **集成Immunex流程**
   - ✅ 统一的输出格式 (CSV/XVG)
   - ⏳ CLI命令接口 (待实现)
   - ⏳ 批处理脚本 (待实现)

**已实现的主要功能**:

1. **主轴计算** (principal_axes.py) ✅
   - 质量加权惯性张量对角化
   - 主轴和主惯性矩提取
   - 单帧和轨迹演化支持

2. **TCR-pMHC对接角度** (docking_angles.py) ✅
   - **Twist角**: TCR二硫键连线投影与MHC主轴夹角
   - **Tilt角**: TCR V域主轴投影与MHC主轴夹角
   - **Swing角**: TCR相对peptide的侧向偏移角度
   - 基于PCA的MHC参考坐标系
   - 二硫键锚点定义TCR取向

3. **向量夹角工具** (vector_angles.py) ✅
   - 向量间夹角计算 (unsigned)
   - 向量投影到平面
   - 平面内有符号角度 (signed)
   - 二面角计算

**使用示例**:
```python
# MD轨迹对接角度分析
from immunex.analysis.angles import DockingAngleAnalyzer

# 初始化分析器
analyzer = DockingAngleAnalyzer('md.tpr', 'md_pbc.xtc')

# 定义原子选择
selections = {
    'mhc': 'chainID A and resid 50:86 140:176 and name CA',
    'tcr_alpha_cys': 'chainID D and resid 22 92 and name CA',
    'tcr_beta_cys': 'chainID E and resid 23 89 and name CA',
    'tcr_v': 'chainID D E and resid 1:115 and name CA',
    'peptide': 'chainID C and name CA'
}

# 计算整条轨迹
times, twist, tilt, swing = analyzer.calculate_docking_angles_trajectory(
    mhc_selection=selections['mhc'],
    tcr_alpha_cys_selection=selections['tcr_alpha_cys'],
    tcr_beta_cys_selection=selections['tcr_beta_cys'],
    tcr_v_selection=selections['tcr_v'],
    peptide_selection=selections['peptide'],
    stride=10,
    output_file='docking_angles.csv'
)

# 统计分析
import numpy as np
print(f"Twist: {np.mean(twist):.2f} ± {np.std(twist):.2f}°")
print(f"Tilt: {np.mean(tilt):.2f} ± {np.std(tilt):.2f}°")
print(f"Swing: {np.mean(swing):.2f} ± {np.std(swing):.2f}°")
```

**更多示例**: `examples/docking_angles_usage.py`

**实施状态**:
1. ✅ `principal_axes.py` - 基础模块 (已完成)
2. ✅ `vector_angles.py` - 工具函数 (已完成)
3. ✅ `docking_angles.py` - 主要功能 (已完成)
4. ⏳ `trajectory_angles.py` - 高层接口 (可选，待实现)

**测试状态**: 7/7 单元测试通过 (100%)

**文档**:
- 实施总结: `development/reports/ANGLES_IMPLEMENTATION_SUMMARY.md`
- 使用示例: `examples/docking_angles_usage.py`
- 单元测试: `development/test_angle_modules.py`

**迁移说明**:
- 旧的C++工具已归档: `development/archived_scripts/tcr_docking_angle_backup_20260123/`
- 旧的Python模块已归档: `development/archived_scripts/angles_old_20260123/`
- 归档代码仍可访问，但建议使用新模块

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
from immunex.analysis.interface import (
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

**文件**: `immunex/utils/batch_analyzer.py` (~400行)

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
from immunex.utils import BatchAnalyzer

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
from immunex.analysis.allostery import ContactCorrelationAnalyzer

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

---

## 🚧 重构计划 (2026-03)

### 当前架构问题

**现有batch实现的主要问题**:
1. **职责混乱**: 业务逻辑和批处理逻辑紧耦合
2. **代码重复**: 10+ batch脚本重复实现任务发现、并行执行
3. **路径管理混乱**: 字符串拼接路径，模块间靠文件名约定传递数据
4. **难以组合**: 无法灵活组合多个分析步骤
5. **缺少中间层**: 直接从Core跳到Batch，缺少Pipeline Nodes和Orchestration层

### 重构目标

**目标**: 将现有系统重构为完整的四层架构

```
Layer 4: Configuration / UI      ← CLI, YAML, Python API
          ↓
Layer 3: Pipeline Orchestration  ← 流程编排，依赖管理，并行调度
          ↓
Layer 2: Pipeline Nodes          ← 轻量封装，上下文管理
          ↓
Layer 1: Core Modules            ← 原子功能，单一职责
```

### 实施策略

**核心原则**:
- ✅ **渐进式重构**: 不推倒重来，逐步迁移
- ✅ **向后兼容**: 保留现有脚本，新旧并存
- ✅ **先核心后边缘**: 先重构高频使用的模块
- ✅ **测试驱动**: 每个阶段都有单元测试和集成测试

**实施阶段** (共12周):
1. **Phase 1** (第1周): 基础设施层 - PipelineContext, PipelineNode基类
2. **Phase 2** (第2-3周): Core Modules重构 - 标准化输入/输出
3. **Phase 3** (第4周): Pipeline Nodes层 - 创建节点封装
4. **Phase 4** (第5周): Pipeline Orchestration - 编排器和BatchExecutor
5. **Phase 5** (第6周): Configuration Layer - YAML配置，CLI命令
6. **Phase 6** (第7-10周): 迁移现有脚本 - 逐步迁移10+ batch脚本
7. **Phase 7** (第11周): 文档和培训
8. **Phase 8** (第12周): 性能优化和稳定性

### 详细方案

**完整文档**:
- 📄 `docs/archive/plans/REFACTORING_PLAN_2026_03.md` - 详细实施方案
- ✅ `docs/archive/plans/REFACTORING_CHECKLIST.md` - 进度检查清单

**核心改进**:

1. **PipelineContext统一数据传递**
```python
@dataclass
class PipelineContext:
    system_id: str
    topology: str
    trajectory_raw: str
    trajectory_processed: Optional[str] = None
    selections: Dict[str, str] = field(default_factory=dict)
    results: Dict[str, Any] = field(default_factory=dict)
    metadata: Dict[str, Any] = field(default_factory=dict)
    output_dir: Optional[str] = None
    errors: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    should_stop: bool = False
```

2. **PipelineNode标准化业务逻辑**
```python
class RMSDNode(PipelineNode):
    def execute(self, context: PipelineContext) -> PipelineContext:
        # 1. 校验输入
        # 2. 调用底层模块
        # 3. 记录输出到context
        # 4. 处理错误
        return context
```

3. **通用BatchExecutor**
```python
executor = BatchExecutor(max_workers=4)
results = executor.execute_pipeline(
    tasks=discover_tasks(base_dir),
    pipeline=StandardTrajectoryPipeline()
)
# ✅ 一套批处理逻辑支持所有分析
```

4. **灵活的Pipeline组合**
```python
class StandardTrajectoryPipeline(Pipeline):
    nodes = [
        QualityCheckNode(),
        PreprocessNode(),
        RMSDNode(),
        RMSFNode(),
        ContactAnalysisNode(),
        ReportNode()
    ]
```

5. **多种配置方式**
```bash
# CLI
immunex pipeline run --config standard_trajectory.yaml

# Python API
pipeline = StandardTrajectoryPipeline()
result = pipeline.execute(context)

# YAML配置
pipeline: standard_trajectory
nodes:
  - name: preprocess
    params: {method: "3step", dt: 10.0}
  - name: rmsd
    params: {selection: "protein"}
```

### 当前状态

**状态**: 📋 计划阶段 - 未开始实施
**优先级**: 🔴 高 - 影响项目长期可维护性
**预计时间**: 12周 (2026-03-16 至 2026-06-01)

**下一步行动**:
1. Review重构方案，确认技术路线
2. 开始Phase 1: 创建PipelineContext和基础设施
3. 为Phase 1编写单元测试

### 参与方式

**问题反馈**: 创建GitHub Issue标记 `refactoring` 标签
**技术讨论**: 参与项目讨论区
**贡献代码**: 按照 `docs/archive/plans/REFACTORING_CHECKLIST.md` 认领任务

