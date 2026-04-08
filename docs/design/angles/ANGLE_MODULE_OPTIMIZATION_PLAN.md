# 角度模块优化方案

**目标**: 将角度计算模块重构为符合Immunex框架标准的模块

**日期**: 2026-03-18

---

## 问题总结

### 当前架构问题
1. ❌ 不符合"6大模块设计标准"（输入明确、输出明确、副作用明确、错误明确、可测试、可调度）
2. ❌ 无法集成到Pipeline架构
3. ❌ 链识别逻辑重复实现
4. ❌ 缺少批处理支持
5. ❌ 输出格式不统一

---

## 优化方案

### Phase 1: 标准化输入/输出 (P0 - 必需)

#### 1.1 创建标准化数据结构

**文件**: `immunex/analysis/angles/angle_data_structures.py` (~150行)

```python
from dataclasses import dataclass, field
from typing import Optional, Dict, List, Any
from pathlib import Path
import numpy as np


@dataclass
class DockingAngleInput:
    """
    标准化角度计算输入

    遵循Immunex模块设计标准 - 输入明确
    """
    # 必需参数
    topology: str
    trajectory: Optional[str] = None

    # 链选择（可选 - 支持自动识别）
    mhc_selection: Optional[str] = None
    tcr_alpha_selection: Optional[str] = None
    tcr_beta_selection: Optional[str] = None

    # 自动链识别开关
    auto_identify_chains: bool = True
    use_anarci: bool = True

    # 分析参数
    stride: int = 1
    output_dir: Optional[str] = None

    def validate(self) -> None:
        """输入验证"""
        # 文件存在性检查
        if not Path(self.topology).exists():
            raise FileNotFoundError(f"Topology not found: {self.topology}")

        if self.trajectory and not Path(self.trajectory).exists():
            raise FileNotFoundError(f"Trajectory not found: {self.trajectory}")

        # 参数范围检查
        if self.stride < 1:
            raise ValueError(f"stride must be >= 1, got {self.stride}")

        # 模式验证
        if not self.auto_identify_chains:
            if not (self.mhc_selection and self.tcr_alpha_selection and self.tcr_beta_selection):
                raise ValueError(
                    "When auto_identify_chains=False, must provide all three selections: "
                    "mhc_selection, tcr_alpha_selection, tcr_beta_selection"
                )


@dataclass
class DockingAngleResult:
    """
    标准化角度计算输出

    遵循Immunex模块设计标准 - 输出明确
    """
    # 执行状态
    success: bool

    # 单帧结果
    crossing_angle: Optional[float] = None
    incident_angle: Optional[float] = None

    # 轨迹结果
    times: Optional[np.ndarray] = None
    crossing_angles: Optional[np.ndarray] = None
    incident_angles: Optional[np.ndarray] = None

    # 统计信息
    statistics: Optional[Dict[str, Any]] = None
    # 示例: {
    #   'crossing_mean': 45.2,
    #   'crossing_std': 3.1,
    #   'incident_mean': 78.5,
    #   'incident_std': 4.2,
    #   'n_frames': 1001
    # }

    # 输出文件
    output_files: List[str] = field(default_factory=list)

    # 元数据
    metadata: Dict[str, Any] = field(default_factory=dict)
    # 示例: {
    #   'module_version': '3.0.0',
    #   'timestamp': '2026-03-18T15:30:00',
    #   'chain_identifications': {...},
    #   'alignment_quality': {...}
    # }

    # 错误信息（如果失败）
    error_message: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        """支持序列化"""
        result = {
            'success': self.success,
            'crossing_angle': self.crossing_angle,
            'incident_angle': self.incident_angle,
            'statistics': self.statistics,
            'output_files': self.output_files,
            'metadata': self.metadata,
            'error_message': self.error_message
        }

        # 数组转换为列表（JSON兼容）
        if self.times is not None:
            result['times'] = self.times.tolist()
        if self.crossing_angles is not None:
            result['crossing_angles'] = self.crossing_angles.tolist()
        if self.incident_angles is not None:
            result['incident_angles'] = self.incident_angles.tolist()

        return result
```

#### 1.2 重构核心分析器

**文件**: `immunex/analysis/angles/docking_angle_analyzer.py` (~300行)

```python
from typing import Optional, Callable
from datetime import datetime
from .angle_data_structures import DockingAngleInput, DockingAngleResult
from .docking_angles_primary import DockingAnglePrimaryAnalyzer


class DockingAngleAnalyzer:
    """
    符合Immunex标准的角度分析器

    遵循6大设计原则：
    1. 输入明确 - DockingAngleInput dataclass
    2. 输出明确 - DockingAngleResult dataclass
    3. 副作用明确 - 记录所有输出文件
    4. 错误明确 - 结构化错误处理
    5. 可测试 - 输入验证 + mock友好
    6. 可调度 - 支持进度回调和取消
    """

    def __init__(self):
        self._cancel_flag = False
        self.progress_callback: Optional[Callable[[float, str], None]] = None

    def set_progress_callback(self, callback: Callable[[float, str], None]):
        """设置进度回调（可调度）"""
        self.progress_callback = callback

    def _report_progress(self, progress: float, message: str):
        """报告进度"""
        if self.progress_callback:
            self.progress_callback(progress, message)

    def cancel(self):
        """取消操作"""
        self._cancel_flag = True

    def is_cancelled(self) -> bool:
        return self._cancel_flag

    def analyze(self, input_params: DockingAngleInput) -> DockingAngleResult:
        """
        执行角度分析

        Parameters
        ----------
        input_params : DockingAngleInput
            标准化输入参数

        Returns
        -------
        DockingAngleResult
            标准化输出结果
        """
        try:
            # 0. 检查取消
            if self.is_cancelled():
                return DockingAngleResult(
                    success=False,
                    error_message="Analysis cancelled by user"
                )

            # 1. 输入验证
            self._report_progress(0.0, "Validating input parameters")
            input_params.validate()

            # 2. 初始化底层分析器
            self._report_progress(0.1, "Initializing analyzer")
            primary_analyzer = DockingAnglePrimaryAnalyzer(
                topology=input_params.topology,
                trajectory=input_params.trajectory,
                auto_identify_chains=input_params.auto_identify_chains,
                use_anarci=input_params.use_anarci
            )

            # 3. 执行分析
            if input_params.trajectory:
                # 轨迹分析
                result = self._analyze_trajectory(
                    primary_analyzer, input_params
                )
            else:
                # 单帧分析
                result = self._analyze_single_frame(
                    primary_analyzer, input_params
                )

            self._report_progress(1.0, "Analysis complete")
            return result

        except Exception as e:
            # 错误明确
            return DockingAngleResult(
                success=False,
                error_message=f"{type(e).__name__}: {str(e)}"
            )

    def _analyze_trajectory(
        self,
        analyzer: DockingAnglePrimaryAnalyzer,
        input_params: DockingAngleInput
    ) -> DockingAngleResult:
        """轨迹分析"""
        self._report_progress(0.2, "Calculating docking angles for trajectory")

        # 调用底层方法
        times, crossing, incident = analyzer.calculate_docking_angles_trajectory(
            mhc_selection=input_params.mhc_selection,
            tcr_alpha_selection=input_params.tcr_alpha_selection,
            tcr_beta_selection=input_params.tcr_beta_selection,
            stride=input_params.stride,
            output_file=None  # 我们自己管理输出
        )

        self._report_progress(0.8, "Calculating statistics")

        # 统计信息
        statistics = {
            'crossing_mean': float(np.mean(crossing)),
            'crossing_std': float(np.std(crossing)),
            'crossing_min': float(np.min(crossing)),
            'crossing_max': float(np.max(crossing)),
            'incident_mean': float(np.mean(incident)),
            'incident_std': float(np.std(incident)),
            'incident_min': float(np.min(incident)),
            'incident_max': float(np.max(incident)),
            'n_frames': len(times)
        }

        # 保存输出文件（副作用明确）
        output_files = []
        if input_params.output_dir:
            output_dir = Path(input_params.output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)

            # CSV文件
            csv_file = output_dir / "docking_angles.csv"
            self._save_csv(csv_file, times, crossing, incident)
            output_files.append(str(csv_file))

        # 元数据
        metadata = {
            'module_version': '3.0.0',
            'timestamp': datetime.now().isoformat(),
            'topology': input_params.topology,
            'trajectory': input_params.trajectory,
            'stride': input_params.stride,
            'auto_identify_chains': input_params.auto_identify_chains
        }

        # 添加链识别信息
        if analyzer.auto_identify_chains and analyzer.identifications:
            metadata['chain_identifications'] = analyzer.identifications

        return DockingAngleResult(
            success=True,
            times=times,
            crossing_angles=crossing,
            incident_angles=incident,
            statistics=statistics,
            output_files=output_files,
            metadata=metadata
        )

    def _analyze_single_frame(
        self,
        analyzer: DockingAnglePrimaryAnalyzer,
        input_params: DockingAngleInput
    ) -> DockingAngleResult:
        """单帧分析"""
        self._report_progress(0.5, "Calculating docking angles for single frame")

        crossing, incident = analyzer.calculate_docking_angles(
            mhc_selection=input_params.mhc_selection,
            tcr_alpha_selection=input_params.tcr_alpha_selection,
            tcr_beta_selection=input_params.tcr_beta_selection
        )

        metadata = {
            'module_version': '3.0.0',
            'timestamp': datetime.now().isoformat(),
            'topology': input_params.topology
        }

        return DockingAngleResult(
            success=True,
            crossing_angle=crossing,
            incident_angle=incident,
            metadata=metadata
        )

    def _save_csv(self, filepath, times, crossing, incident):
        """保存CSV文件"""
        import csv
        with open(filepath, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['Time(ps)', 'Crossing(deg)', 'Incident(deg)'])
            for t, c, i in zip(times, crossing, incident):
                writer.writerow([f'{t:.2f}', f'{c:.4f}', f'{i:.4f}'])
```

---

### Phase 2: Pipeline集成 (P0 - 必需)

#### 2.1 创建Pipeline Node

**文件**: `immunex/pipeline/nodes/docking_angle_node.py` (~150行)

```python
from typing import Optional
from immunex.pipeline.pipeline_node import PipelineNode
from immunex.pipeline.pipeline_context import PipelineContext
from immunex.analysis.angles import DockingAngleAnalyzer, DockingAngleInput


class DockingAngleNode(PipelineNode):
    """
    对接角度分析Pipeline节点

    职责：
    1. 从context获取输入参数
    2. 调用DockingAngleAnalyzer
    3. 将结果写入context
    4. 处理错误
    """

    def __init__(self, stride: int = 1, output_basename: str = "docking_angles"):
        self.stride = stride
        self.output_basename = output_basename
        self.analyzer = DockingAngleAnalyzer()

    def execute(self, context: PipelineContext) -> PipelineContext:
        """执行角度分析"""
        try:
            # 1. 校验输入
            if 'trajectory_processed' not in context.metadata:
                context.add_error("Missing processed trajectory")
                context.should_stop = True
                return context

            # 2. 构建输入参数
            input_params = DockingAngleInput(
                topology=context.topology,
                trajectory=context.metadata['trajectory_processed'],
                mhc_selection=context.selections.get('mhc'),
                tcr_alpha_selection=context.selections.get('tcr_alpha'),
                tcr_beta_selection=context.selections.get('tcr_beta'),
                auto_identify_chains=(
                    context.selections.get('mhc') is None
                ),
                stride=self.stride,
                output_dir=context.get_output_path("angles")
            )

            # 3. 执行分析
            result = self.analyzer.analyze(input_params)

            # 4. 处理结果
            if result.success:
                # 写入context
                context.results['docking_angles'] = {
                    'crossing_angle': result.crossing_angle,
                    'incident_angle': result.incident_angle,
                    'statistics': result.statistics,
                    'output_files': result.output_files
                }

                # 合并metadata
                context.metadata.update({
                    'angle_analysis': result.metadata
                })
            else:
                context.add_error(f"Angle analysis failed: {result.error_message}")
                context.should_stop = True

        except Exception as e:
            context.add_error(f"DockingAngleNode error: {str(e)}")
            context.should_stop = True

        return context
```

#### 2.2 创建标准Pipeline

**文件**: `immunex/pipeline/standard_analysis_pipelines.py` (添加)

```python
class DockingAnglePipeline(Pipeline):
    """对接角度分析Pipeline"""

    def __init__(self, stride: int = 1):
        nodes = [
            ChainIdentificationNode(),  # 复用统一链识别
            DockingAngleNode(stride=stride),
            ReportNode()
        ]
        super().__init__(nodes=nodes)
```

**使用示例**:
```python
from pathlib import Path

from immunex.core import discover_tasks
from immunex.pipeline import DockingAnglePipeline, BatchExecutor

# 发现任务
report = discover_tasks(Path("data/simulations"))

# 批量执行
pipeline = DockingAnglePipeline(stride=10)
executor = BatchExecutor(max_workers=4)
results = executor.execute_pipeline(report, pipeline)
```

---

### Phase 3: 消除链识别重复 (P1 - 推荐)

**策略**: 不在角度模块内部实现链识别，而是：

1. **从PipelineContext获取**（推荐）:
```python
# Pipeline模式
context = PipelineContext(...)
context = ChainIdentificationNode().execute(context)
# context.metadata['chain_mapping'] = {'mhc_alpha': 'A', ...}

context = DockingAngleNode().execute(context)
# 从context.metadata获取链信息
```

2. **使用独立的链识别适配器**（兼容旧代码）:
```python
# 独立使用模式
from immunex.analysis import ChainIdentificationAdapter

adapter = ChainIdentificationAdapter()
identifications = adapter.identify_chains(topology)

analyzer = DockingAngleAnalyzer()
result = analyzer.analyze(DockingAngleInput(
    topology=topology,
    trajectory=trajectory,
    mhc_selection=identifications['mhc_alpha_selection'],
    tcr_alpha_selection=identifications['tcr_alpha_selection'],
    tcr_beta_selection=identifications['tcr_beta_selection'],
    auto_identify_chains=False  # 已手动识别
))
```

**修改建议**:
- 保留`docking_angles_primary.py`的当前API（向后兼容）
- 新代码使用`DockingAngleAnalyzer` + Pipeline模式
- 逐步迁移旧用户代码

---

### Phase 4: 批处理标准化 (P1 - 推荐)

**删除**: 自定义批处理循环

**使用**: `BatchExecutor` + `Pipeline`

**迁移示例**:
```python
# 旧代码 ❌
tasks = [...]
for task in tasks:
    analyzer = DockingAnglePrimaryAnalyzer(
        task['topology'], task['trajectory']
    )
    times, crossing, incident = analyzer.calculate_docking_angles_trajectory(...)

# 新代码 ✓
from immunex.pipeline import DockingAnglePipeline, BatchExecutor

pipeline = DockingAnglePipeline(stride=10)
executor = BatchExecutor(max_workers=4)
results = executor.execute_pipeline(tasks, pipeline)

# 统一的结果汇总
summary = executor.summarize_results(results)
```

---

## 实施计划

### 优先级

**P0 - 必须完成** (核心功能):
1. ✅ 创建`angle_data_structures.py` (标准化输入/输出)
2. ✅ 创建`DockingAngleAnalyzer` (符合6大标准)
3. ✅ 创建`DockingAngleNode` (Pipeline集成)
4. ✅ 添加到`standard_analysis_pipelines.py`
5. ✅ 单元测试

**P1 - 推荐完成** (架构优化):
6. 消除链识别重复（复用`ChainIdentificationAdapter`）
7. 删除自定义批处理代码
8. 文档更新

**P2 - 可选完成** (增强功能):
9. 添加更多角度类型（如Swing/Tilt）
10. 可视化生成（PyMOL脚本、角度分布图）

### 时间估算

- **Phase 1**: 1天（P0）
- **Phase 2**: 0.5天（P0）
- **Phase 3**: 0.5天（P1）
- **Phase 4**: 0.5天（P1）

**总计**: 2.5天

---

## 向后兼容性

**保留旧API**:
```python
# 旧代码仍然工作 ✓
from immunex.analysis.angles import DockingAnglePrimaryAnalyzer

analyzer = DockingAnglePrimaryAnalyzer('md.tpr', 'md.xtc')
crossing, incident = analyzer.calculate_docking_angles(...)
```

**推荐新API**:
```python
# 新代码使用标准化接口 ✓
from immunex.analysis.angles import DockingAngleAnalyzer, DockingAngleInput

analyzer = DockingAngleAnalyzer()
result = analyzer.analyze(DockingAngleInput(
    topology='md.tpr',
    trajectory='md.xtc',
    stride=10
))

print(f"Crossing: {result.statistics['crossing_mean']:.2f}°")
```

---

## 测试策略

### 单元测试

**文件**: `tests/test_docking_angle_analyzer.py`

```python
def test_input_validation():
    """测试输入验证"""
    with pytest.raises(FileNotFoundError):
        DockingAngleInput(topology="nonexistent.pdb").validate()

def test_single_frame_analysis():
    """测试单帧分析"""
    analyzer = DockingAngleAnalyzer()
    result = analyzer.analyze(DockingAngleInput(
        topology="test_data/1ao7.pdb"
    ))
    assert result.success
    assert result.crossing_angle is not None

def test_trajectory_analysis():
    """测试轨迹分析"""
    analyzer = DockingAngleAnalyzer()
    result = analyzer.analyze(DockingAngleInput(
        topology="test_data/md.tpr",
        trajectory="test_data/md.xtc",
        stride=100
    ))
    assert result.success
    assert result.statistics is not None
```

### 集成测试

**文件**: `tests/test_docking_angle_pipeline.py`

```python
def test_pipeline_integration():
    """测试Pipeline集成"""
    context = PipelineContext(
        system_id="1ao7",
        topology="test_data/md.tpr",
        trajectory_raw="test_data/md.xtc"
    )

    pipeline = DockingAnglePipeline(stride=10)
    result = pipeline.execute(context)

    assert not result.has_errors()
    assert 'docking_angles' in result.results
```

---

## 受益总结

### 1. 一致性
- 所有模块遵循相同的设计标准
- 统一的输入/输出格式
- 统一的错误处理

### 2. 可维护性
- 清晰的职责分离
- 减少代码重复
- 易于测试

### 3. 可扩展性
- 易于添加新的角度类型
- 易于集成到不同的workflow
- 支持批处理和并行

### 4. 用户体验
- Pipeline模式：无缝组合多个分析
- 标准化输出：易于后处理
- 进度报告：长时间任务反馈

---

## 参考文献

- `docs/design/angles/ANGLE_MODULE_REDESIGN.md` - 角度模块设计文档
- `docs/specs/DATA_STRUCTURE_STANDARD.md` - 数据结构标准
- `CLAUDE.md` - 模块设计6大原则
