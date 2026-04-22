# Immunex 四层架构重构方案

**制定日期**: 2026-03-16
**目标**: 将现有batch处理系统重构为符合四层架构设计的模块化系统

---

## 📋 目录

1. [整改目标](#整改目标)
2. [整改策略](#整改策略)
3. [实施阶段](#实施阶段)
4. [技术方案](#技术方案)
5. [风险控制](#风险控制)
6. [验收标准](#验收标准)

---

## 🎯 整改目标

### 核心问题
- ❌ **职责混乱**: 业务逻辑和批处理逻辑紧耦合
- ❌ **代码重复**: 10+ batch脚本重复实现任务发现、并行执行
- ❌ **路径管理混乱**: 字符串拼接路径，模块间靠文件名约定传递数据
- ❌ **难以组合**: 无法灵活组合多个分析步骤
- ❌ **测试困难**: 缺少统一的接口，难以mock和单元测试

### 期望成果
- ✅ **清晰的四层架构**: Core → Nodes → Pipeline → Interface
- ✅ **统一的Context管理**: PipelineContext传递数据，消除路径字符串
- ✅ **可组合的Pipeline**: 灵活组合分析步骤
- ✅ **通用的BatchExecutor**: 一套批处理逻辑支持所有分析
- ✅ **向后兼容**: 旧脚本继续工作，逐步迁移

---

## 🔧 整改策略

### 核心原则
1. **渐进式重构**: 不推倒重来，逐步迁移
2. **向后兼容**: 保留现有脚本，新旧并存
3. **先核心后边缘**: 先重构高频使用的模块
4. **测试驱动**: 每个阶段都有单元测试和集成测试

### 迁移路径
```
现有系统 (保留)
    ↓
新架构核心组件 (Phase 1-2)
    ↓
示例Pipeline (Phase 3)
    ↓
逐步迁移现有脚本 (Phase 4-5)
    ↓
弃用旧脚本 (Phase 6)
```

---

## 📅 实施阶段

### **Phase 1: 基础设施层** (优先级: 🔴 最高)
**时间**: 第1周
**目标**: 建立核心数据结构和基础工具

#### 任务清单
- [ ] 1.1 创建 `PipelineContext` 数据类
- [ ] 1.2 创建 `ProcessingResult` 标准化输出
- [ ] 1.3 创建自定义异常类 `ImmunexError`
- [ ] 1.4 创建 `PipelineNode` 抽象基类
- [ ] 1.5 编写单元测试

#### 产出物
```
immunex/core/
├── __init__.py
├── context.py           # PipelineContext, ProcessingResult
├── exceptions.py        # ImmunexError, InputValidationError, etc.
└── base_node.py         # PipelineNode抽象基类

tests/core/
├── test_context.py
├── test_exceptions.py
└── test_base_node.py
```

---

### **Phase 2: Core Modules 重构** (优先级: 🔴 最高)
**时间**: 第2-3周
**目标**: 将现有分析模块改造为符合六大设计原则的原子模块

#### 任务清单
- [ ] 2.1 审计现有模块（PBCProcessor, RMSDCalculator等）
- [ ] 2.2 重构 `PBCProcessor` - 标准化输入/输出
- [ ] 2.3 重构 `RMSDCalculator` - 移除上下文依赖
- [ ] 2.4 重构 `RMSFCalculator`
- [ ] 2.5 重构 `ContactAnalyzer`
- [ ] 2.6 为每个模块编写单元测试

#### 产出物
```
immunex/analysis/trajectory/
├── rmsd.py              # 重构后的RMSDCalculator
├── rmsf.py              # 重构后的RMSFCalculator
└── ...

immunex/core/
└── pbc_processor.py     # 重构后的PBCProcessor

tests/analysis/
├── test_rmsd.py
├── test_rmsf.py
└── ...
```

#### 重构示例
**重构前**:
```python
class RMSDCalculator:
    def calculate_rmsd(self, topology, trajectory, output_dir):
        # 硬编码文件路径
        output_file = os.path.join(output_dir, "rmsd.xvg")
        # 返回字符串路径
        return output_file
```

**重构后**:
```python
@dataclass
class RMSDInput:
    topology: str
    trajectory: str
    selection: str
    reference_frame: int = 0

    def validate(self):
        if not Path(self.topology).exists():
            raise FileNotFoundError(f"Topology not found: {self.topology}")

@dataclass
class RMSDResult:
    success: bool
    output_file: str
    mean_rmsd: float
    std_rmsd: float
    metadata: Dict[str, Any]
    error_message: Optional[str] = None

class RMSDCalculator:
    def calculate(self, input_params: RMSDInput) -> RMSDResult:
        """职责单一：只计算RMSD"""
        input_params.validate()
        # 计算逻辑
        return RMSDResult(...)
```

---

### **Phase 3: Pipeline Nodes 层** (优先级: 🟡 高)
**时间**: 第4周
**目标**: 创建Pipeline Nodes层，封装业务逻辑

#### 任务清单
- [ ] 3.1 创建 `PreprocessNode` (PBC处理)
- [ ] 3.2 创建 `RMSDNode`
- [ ] 3.3 创建 `RMSFNode`
- [ ] 3.4 创建 `ContactAnalysisNode`
- [ ] 3.5 创建 `QualityCheckNode`
- [ ] 3.6 编写节点单元测试

#### 产出物
```
immunex/pipeline/
├── __init__.py
├── nodes/
│   ├── __init__.py
│   ├── preprocess_node.py    # PreprocessNode
│   ├── rmsd_node.py           # RMSDNode
│   ├── rmsf_node.py           # RMSFNode
│   ├── contact_node.py        # ContactAnalysisNode
│   └── quality_node.py        # QualityCheckNode

tests/pipeline/nodes/
├── test_preprocess_node.py
├── test_rmsd_node.py
└── ...
```

#### 节点实现示例
```python
class RMSDNode(PipelineNode):
    """RMSD计算节点 - 不亲自计算，而是编排"""

    def __init__(self, selection: str = "protein"):
        self.selection = selection

    def execute(self, context: PipelineContext) -> PipelineContext:
        """执行RMSD计算并更新context"""

        # 1. 校验输入
        if 'trajectory_processed' not in context or not context.trajectory_processed:
            raise ValueError("Missing processed trajectory - run PreprocessNode first")

        # 2. 调用底层模块
        calculator = RMSDCalculator()
        input_params = RMSDInput(
            topology=context.topology,
            trajectory=context.trajectory_processed,
            selection=self.selection
        )

        result = calculator.calculate(input_params)

        # 3. 记录输出到context
        if result.success:
            output_path = context.get_output_path("rmsd.xvg")
            context.results['rmsd'] = {
                'output_file': result.output_file,
                'mean': result.mean_rmsd,
                'std': result.std_rmsd,
                'metadata': result.metadata
            }
        else:
            context.errors.append(f"RMSD calculation failed: {result.error_message}")
            context.should_stop = True

        # 4. 返回更新后的context
        return context
```

---

### **Phase 4: Pipeline Orchestration** (优先级: 🟡 高)
**时间**: 第5周
**目标**: 实现Pipeline编排器和通用BatchExecutor

#### 任务清单
- [ ] 4.1 创建 `Pipeline` 基类
- [ ] 4.2 创建预定义Pipeline (StandardTrajectoryPipeline等)
- [ ] 4.3 创建 `PipelineOrchestrator` (支持DAG依赖)
- [ ] 4.4 创建通用 `BatchExecutor`
- [ ] 4.5 创建 `TaskDiscovery` 工具
- [ ] 4.6 编写集成测试

#### 产出物
```
immunex/pipeline/
├── base_pipeline.py          # Pipeline基类
├── orchestrator.py           # PipelineOrchestrator
├── standard_pipelines.py     # 预定义Pipeline
└── batch_executor.py         # BatchExecutor

immunex/core/
└── task_discovery.py         # TaskDiscovery

tests/pipeline/
├── test_orchestrator.py
├── test_batch_executor.py
└── test_integration.py
```

#### Pipeline实现示例
```python
class StandardTrajectoryPipeline(Pipeline):
    """标准轨迹分析流程"""

    def __init__(self):
        self.nodes = [
            QualityCheckNode(),
            PreprocessNode(),
            RMSDNode(selection="protein"),
            RMSFNode(selection="backbone"),
            ContactAnalysisNode(cutoff=4.0),
            ReportNode()
        ]

    def execute(self, context: PipelineContext) -> PipelineContext:
        """顺序执行所有节点"""
        for node in self.nodes:
            if context.should_stop:
                break

            try:
                context = node.execute(context)
            except Exception as e:
                logger.error(f"Node {node.__class__.__name__} failed: {e}")
                context.errors.append(str(e))
                context.should_stop = True
                break

        return context


class ParallelAnalysisPipeline(Pipeline):
    """并行分析流程"""

    def execute(self, context: PipelineContext) -> PipelineContext:
        # 预处理（必需）
        context = PreprocessNode().execute(context)

        # 并行执行多个分析
        parallel_nodes = [
            RMSDNode(),
            RMSFNode(),
            ContactAnalysisNode()
        ]

        with ThreadPoolExecutor(max_workers=3) as executor:
            futures = [
                executor.submit(node.execute, context.copy())
                for node in parallel_nodes
            ]
            results = [f.result() for f in futures]

        # 合并结果
        context = self.merge_contexts(results)
        return context
```

#### BatchExecutor实现
```python
class BatchExecutor:
    """通用批处理执行器 - 可执行任意Pipeline"""

    def __init__(self, max_workers: int = 4):
        self.max_workers = max_workers

    def execute_pipeline(self,
                        tasks: List[PipelineContext],
                        pipeline: Pipeline,
                        show_progress: bool = True) -> List[PipelineContext]:
        """
        批量执行Pipeline

        Args:
            tasks: 任务列表（每个是一个PipelineContext）
            pipeline: 要执行的Pipeline
            show_progress: 是否显示进度条

        Returns:
            执行结果列表
        """

        with ProcessPoolExecutor(max_workers=self.max_workers) as executor:
            futures = {
                executor.submit(pipeline.execute, task): task
                for task in tasks
            }

            results = []
            for future in tqdm(as_completed(futures), total=len(tasks),
                              disable=not show_progress):
                try:
                    result = future.result()
                    results.append(result)
                except Exception as e:
                    task = futures[future]
                    logger.error(f"Task {task.system_id} failed: {e}")
                    task.errors.append(str(e))
                    results.append(task)

        return results

    def summarize_results(self, results: List[PipelineContext]) -> Dict:
        """汇总批处理结果"""
        successful = sum(1 for r in results if not r.errors)
        failed = len(results) - successful

        return {
            'total_tasks': len(results),
            'successful': successful,
            'failed': failed,
            'success_rate': successful / len(results) * 100 if results else 0,
            'results': results
        }
```

---

### **Phase 5: Configuration Layer** (优先级: 🟢 中)
**时间**: 第6周
**目标**: 提供多种入口方式

#### 任务清单
- [ ] 5.1 创建 YAML 配置解析器
- [ ] 5.2 创建统一的 CLI 命令 `immunex pipeline`
- [ ] 5.3 创建 Python API 接口
- [ ] 5.4 编写配置文件示例
- [ ] 5.5 编写使用文档

#### 产出物
```
immunex/config/
├── __init__.py
├── parser.py            # YAML/JSON配置解析
└── validator.py         # 配置验证

immunex/cli/commands/
└── pipeline.py          # pipeline子命令

examples/configs/
├── standard_trajectory.yaml
├── cdr_analysis.yaml
└── batch_processing.yaml

docs/
└── PIPELINE_USER_GUIDE.md
```

#### 配置文件示例
```yaml
# standard_trajectory.yaml
pipeline:
  name: standard_trajectory
  description: "Standard trajectory analysis pipeline"

input:
  base_directory: "/data/simulations"
  file_pattern: "*/md.xtc"
  topology_pattern: "*/md.tpr"

nodes:
  - name: quality_check
    enabled: true
    params:
      min_frames: 1000

  - name: preprocess
    enabled: true
    params:
      method: "3step"
      dt: 10.0

  - name: rmsd
    enabled: true
    params:
      selection: "protein"
      reference_frame: 0

  - name: rmsf
    enabled: true
    params:
      selection: "name CA"

  - name: contact_analysis
    enabled: true
    params:
      cutoff: 4.0
      selection_a: "chainID D E"
      selection_b: "chainID A B C"

output:
  base_directory: "./results"
  save_context: true

execution:
  max_workers: 4
  show_progress: true
  stop_on_error: false
```

#### CLI 使用
```bash
# 使用YAML配置文件
immunex pipeline run --config standard_trajectory.yaml

# 命令行参数
immunex pipeline run \
  --input /data/simulations \
  --output ./results \
  --pipeline standard \
  --workers 4

# 列出可用的Pipeline
immunex pipeline list

# 验证配置文件
immunex pipeline validate --config myconfig.yaml

# 仅发现任务，不执行
immunex pipeline discover --input /data/simulations
```

#### Python API 使用
```python
# 方式1: 直接使用Pipeline
from immunex.pipeline import StandardTrajectoryPipeline
from immunex.core import PipelineContext

context = PipelineContext(
    system_id="1ao7",
    topology="data/1ao7/md.tpr",
    trajectory_raw="data/1ao7/md.xtc"
)

pipeline = StandardTrajectoryPipeline()
result = pipeline.execute(context)

# 方式2: 使用BatchExecutor
from immunex.pipeline import BatchExecutor
from immunex.core import TaskDiscovery

# 发现任务
discovery = TaskDiscovery()
tasks = discovery.discover("/data/simulations")

# 批量执行
executor = BatchExecutor(max_workers=4)
results = executor.execute_pipeline(tasks, StandardTrajectoryPipeline())

summary = executor.summarize_results(results)
print(f"成功: {summary['successful']}/{summary['total_tasks']}")

# 方式3: 从配置文件
from immunex.config import PipelineConfig

config = PipelineConfig.from_yaml("standard_trajectory.yaml")
pipeline = config.create_pipeline()
tasks = config.discover_tasks()

executor = BatchExecutor(max_workers=config.execution.max_workers)
results = executor.execute_pipeline(tasks, pipeline)
```

---

### **Phase 6: 迁移现有脚本** (优先级: 🟢 中)
**时间**: 第7-10周
**目标**: 逐步将现有batch_*.py脚本迁移到新架构

#### 迁移优先级
1. **高频使用** (第7周)
   - [ ] `batch_pbc_slurm.py` → `PBCPipeline`
   - [ ] `batch_cdr_rmsd_exact.py` → `CDRAnalysisPipeline`
   - [ ] `batch_phla_analysis.py` → `pHLAAnalysisPipeline`

2. **中频使用** (第8周)
   - [ ] `batch_contact_frequency.py` → `ContactPipeline`
   - [ ] `batch_interface_rmsd.py` → `InterfaceRMSDPipeline`
   - [ ] `batch_docking_angles.py` → `DockingAnglePipeline`

3. **低频使用** (第9周)
   - [ ] `batch_cdr_rmsf.py`
   - [ ] `batch_whole_protein_rmsf.py`
   - [ ] `batch_allostery_analysis.py`

4. **特殊脚本** (第10周)
   - [ ] `batch_chain_identification.py`
   - [ ] `batch_pdb_pipeline.py`

#### 迁移步骤
```
对于每个脚本：

1. 分析现有逻辑
   - 识别Core Module（分析逻辑）
   - 识别配置参数
   - 识别任务发现规则

2. 创建对应的Nodes
   - 将分析逻辑封装为Node
   - 使用PipelineContext传递数据

3. 创建Pipeline
   - 定义节点执行顺序
   - 处理依赖关系

4. 创建配置文件
   - YAML配置示例
   - 文档说明

5. 创建兼容脚本
   - 保留旧脚本名称
   - 内部调用新Pipeline
   - 添加deprecation警告

6. 测试
   - 对比新旧输出一致性
   - 性能测试
   - 集成测试
```

#### 迁移示例: batch_cdr_rmsd_exact.py
```python
# 新建 immunex/pipeline/cdr_analysis_pipeline.py
class CDRAnalysisPipeline(Pipeline):
    def __init__(self, cdr_csv_file: str):
        self.nodes = [
            LoadCDRReferenceNode(cdr_csv_file),
            ChainIdentificationNode(),
            CDRMappingNode(),
            CDRRMSDNode(),
            CDRReportNode()
        ]

# 兼容性包装 - scripts/batch_cdr_rmsd_exact.py
import warnings
from immunex.pipeline import CDRAnalysisPipeline, BatchExecutor
from immunex.core import TaskDiscovery

warnings.warn(
    "batch_cdr_rmsd_exact.py is deprecated. Use 'immunex pipeline run' instead.",
    DeprecationWarning
)

def main():
    # 保持原有命令行接口
    # 内部调用新Pipeline
    pipeline = CDRAnalysisPipeline(cdr_csv_file=args.cdr_csv)

    discovery = TaskDiscovery()
    tasks = discovery.discover(args.input_dirs)

    executor = BatchExecutor(max_workers=args.max_workers)
    results = executor.execute_pipeline(tasks, pipeline)

    # 保持原有输出格式
    # ...
```

---

### **Phase 7: 文档和培训** (优先级: 🟢 中)
**时间**: 第11周
**目标**: 完善文档，培训用户

#### 任务清单
- [ ] 7.1 编写架构设计文档
- [ ] 7.2 编写Pipeline开发指南
- [ ] 7.3 编写迁移指南
- [ ] 7.4 编写API文档
- [ ] 7.5 录制教学视频（可选）
- [ ] 7.6 更新README和QUICKSTART

#### 产出物
```
docs/
├── ARCHITECTURE_REDESIGN.md      # 架构设计文档
├── PIPELINE_DEVELOPMENT_GUIDE.md # Pipeline开发指南
├── MIGRATION_GUIDE.md            # 从旧脚本迁移指南
├── API_REFERENCE.md              # API文档
└── EXAMPLES.md                   # 使用示例

examples/
├── basic_pipeline_usage.py
├── custom_pipeline_example.py
├── batch_processing_example.py
└── advanced_orchestration.py
```

---

### **Phase 8: 性能优化和稳定性** (优先级: 🔵 低)
**时间**: 第12周
**目标**: 性能优化和稳定性提升

#### 任务清单
- [ ] 8.1 性能基准测试
- [ ] 8.2 内存优化（大规模批处理）
- [ ] 8.3 错误恢复机制（断点续传）
- [ ] 8.4 日志和监控增强
- [ ] 8.5 压力测试（1000+ 任务）
- [ ] 8.6 CI/CD集成

---

## 🛠️ 技术方案

### 核心数据结构

#### PipelineContext
```python
# immunex/core/context.py
from dataclasses import dataclass, field
from typing import Dict, Any, Optional, List
from pathlib import Path
import json

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

    # 分析结果
    results: Dict[str, Any] = field(default_factory=dict)

    # 元数据
    metadata: Dict[str, Any] = field(default_factory=dict)

    # 输出目录
    output_dir: Optional[str] = None

    # 错误和警告
    errors: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)

    # 流程控制
    should_stop: bool = False

    # 临时文件追踪（用于清理）
    temporary_files: List[str] = field(default_factory=list)

    def get_output_path(self, filename: str) -> str:
        """获取输出文件路径"""
        if self.output_dir is None:
            self.output_dir = f"./results/{self.system_id}"
        Path(self.output_dir).mkdir(parents=True, exist_ok=True)
        return str(Path(self.output_dir) / filename)

    def copy(self) -> 'PipelineContext':
        """深拷贝context（用于并行节点）"""
        import copy
        return copy.deepcopy(self)

    def to_dict(self) -> Dict[str, Any]:
        """序列化为字典"""
        from dataclasses import asdict
        return asdict(self)

    def save(self, filepath: str):
        """保存context到JSON文件"""
        with open(filepath, 'w') as f:
            json.dump(self.to_dict(), f, indent=2)

    @classmethod
    def load(cls, filepath: str) -> 'PipelineContext':
        """从JSON文件加载context"""
        with open(filepath, 'r') as f:
            data = json.load(f)
        return cls(**data)

    def add_result(self, key: str, value: Any):
        """添加分析结果"""
        self.results[key] = value

    def get_result(self, key: str, default: Any = None) -> Any:
        """获取分析结果"""
        return self.results.get(key, default)

    def has_errors(self) -> bool:
        """是否有错误"""
        return len(self.errors) > 0

    def add_error(self, error: str):
        """添加错误"""
        self.errors.append(error)

    def add_warning(self, warning: str):
        """添加警告"""
        self.warnings.append(warning)
```

#### PipelineNode基类
```python
# immunex/core/base_node.py
from abc import ABC, abstractmethod
from typing import Optional
import logging

logger = logging.getLogger(__name__)

class PipelineNode(ABC):
    """Pipeline节点抽象基类"""

    def __init__(self, name: Optional[str] = None):
        self.name = name or self.__class__.__name__
        self.logger = logging.getLogger(f"immunex.pipeline.{self.name}")

    @abstractmethod
    def execute(self, context: PipelineContext) -> PipelineContext:
        """
        执行节点逻辑

        Args:
            context: Pipeline上下文

        Returns:
            更新后的Pipeline上下文

        Raises:
            可以抛出任何异常，由Pipeline捕获处理
        """
        pass

    def validate_inputs(self, context: PipelineContext):
        """
        校验输入（子类可选重写）

        默认不做任何校验
        """
        pass

    def __repr__(self):
        return f"{self.__class__.__name__}(name='{self.name}')"
```

### 任务发现机制

```python
# immunex/core/task_discovery.py
from pathlib import Path
from typing import List, Dict, Optional, Callable
import logging

logger = logging.getLogger(__name__)

class TaskDiscovery:
    """任务发现工具 - 自动扫描目录发现MD任务"""

    def __init__(self):
        self.discovery_rules = []

    def add_rule(self, rule: Callable[[Path], Optional[Dict]]):
        """添加自定义发现规则"""
        self.discovery_rules.append(rule)

    def discover(self,
                base_directory: str,
                pattern: str = "*",
                require_files: Optional[List[str]] = None) -> List[PipelineContext]:
        """
        发现任务

        Args:
            base_directory: 基础目录
            pattern: 任务目录匹配模式
            require_files: 必需文件列表

        Returns:
            PipelineContext列表
        """
        base_path = Path(base_directory)

        if not base_path.exists():
            raise ValueError(f"Base directory not found: {base_directory}")

        tasks = []

        for task_dir in sorted(base_path.glob(pattern)):
            if not task_dir.is_dir():
                continue

            task_name = task_dir.name

            # 使用默认规则或自定义规则发现
            if self.discovery_rules:
                for rule in self.discovery_rules:
                    task_data = rule(task_dir)
                    if task_data:
                        tasks.append(self._create_context(task_name, task_data))
                        break
            else:
                # 默认发现规则
                task_data = self._default_discovery(task_dir, require_files)
                if task_data:
                    tasks.append(self._create_context(task_name, task_data))

        logger.info(f"Discovered {len(tasks)} tasks in {base_directory}")
        return tasks

    def _default_discovery(self, task_dir: Path, require_files: Optional[List[str]]) -> Optional[Dict]:
        """默认发现规则：查找md.xtc和md.tpr"""

        # 搜索位置
        search_locations = [
            task_dir,
            task_dir / "prod"
        ]

        for location in search_locations:
            if not location.exists():
                continue

            trajectory = location / "md.xtc"
            topology = location / "md.tpr"

            if trajectory.exists() and topology.exists():
                # 检查其他必需文件
                if require_files:
                    missing = [f for f in require_files if not (location / f).exists()]
                    if missing:
                        continue

                return {
                    'topology': str(topology),
                    'trajectory_raw': str(trajectory),
                    'task_dir': str(task_dir)
                }

        return None

    def _create_context(self, task_name: str, task_data: Dict) -> PipelineContext:
        """创建PipelineContext"""
        return PipelineContext(
            system_id=task_name,
            topology=task_data['topology'],
            trajectory_raw=task_data['trajectory_raw'],
            metadata={'task_dir': task_data.get('task_dir', '')}
        )
```

---

## ⚠️ 风险控制

### 风险识别

| 风险 | 影响 | 概率 | 缓解措施 |
|------|------|------|---------|
| **向后兼容性破坏** | 高 | 中 | 保留旧脚本，添加deprecation警告，逐步迁移 |
| **性能下降** | 中 | 低 | 性能基准测试，优化热点代码 |
| **学习曲线** | 中 | 高 | 详细文档，示例代码，培训支持 |
| **时间延期** | 中 | 中 | 分阶段实施，每阶段可独立交付 |
| **测试覆盖不足** | 高 | 中 | 强制单元测试，集成测试，覆盖率>80% |

### 回滚策略
- **Phase 1-2**: 如果基础设施有问题，可以完全回滚
- **Phase 3-5**: 新旧系统并存，可随时切换
- **Phase 6**: 每个脚本独立迁移，失败可单独回滚

### 测试策略
```
每个Phase都需要：
1. 单元测试 - 覆盖率 >= 80%
2. 集成测试 - 端到端测试
3. 性能测试 - 与旧版本对比
4. 兼容性测试 - 确保旧脚本继续工作
```

---

## ✅ 验收标准

### Phase 1-2 验收
- [ ] PipelineContext可序列化/反序列化
- [ ] 所有Core Modules有单元测试
- [ ] 测试覆盖率 >= 80%
- [ ] 文档完整

### Phase 3-4 验收
- [ ] 至少3个示例Pipeline可运行
- [ ] BatchExecutor支持并行执行
- [ ] 错误处理正确
- [ ] 性能不低于旧版本

### Phase 5 验收
- [ ] YAML配置可正常解析
- [ ] CLI命令正常工作
- [ ] Python API文档完整
- [ ] 至少5个配置文件示例

### Phase 6 验收
- [ ] 所有高频脚本已迁移
- [ ] 输出结果与旧版本一致
- [ ] 旧脚本有deprecation警告
- [ ] 迁移文档完整

### 最终验收
- [ ] 所有测试通过
- [ ] 文档完整
- [ ] 性能达标（不低于旧版本）
- [ ] 用户满意度调查 >= 4/5

---

## 📚 参考资料

### 设计原则
- CLAUDE.md - 四层架构设计
- CLAUDE.md - 六大模块设计标准
- CLAUDE.md - PipelineContext设计

### 相关文档
- ARCHITECTURE.md - 当前架构文档
- README.md - 项目介绍

### 外部参考
- Apache Airflow - DAG编排
- Prefect - 现代workflow引擎
- Luigi - Spotify的batch框架

---

## 📞 联系方式

**问题反馈**: 创建GitHub Issue
**技术讨论**: 项目讨论区
**紧急问题**: 联系项目维护者

---

**文档版本**: v1.0
**最后更新**: 2026-03-16
