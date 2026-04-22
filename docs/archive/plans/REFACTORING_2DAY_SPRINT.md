# Immunex 2天冲刺重构计划

**开始时间**: 2026-03-16
**截止时间**: 2026-03-18
**总工时**: 16小时 (2天 × 8小时)

---

## 🎯 核心目标

**MVP原则**: 用2天时间交付最小可用的四层架构，证明概念可行

### 必须交付
1. ✅ 核心基础设施 (PipelineContext, PipelineNode基类)
2. ✅ 1个重构的Core Module (RMSDCalculator)
3. ✅ 2-3个Pipeline Nodes (Preprocess, RMSD, Report)
4. ✅ 1个完整的示例Pipeline (StandardTrajectoryPipeline)
5. ✅ BatchExecutor (简化版)
6. ✅ TaskDiscovery工具
7. ✅ 基础文档和使用示例

### 暂不交付 (后续渐进)
- ❌ 全面的单元测试覆盖 (只测核心功能)
- ❌ YAML配置解析 (先用Python API)
- ❌ CLI命令 (先用脚本)
- ❌ 迁移所有现有脚本 (提供1个迁移示例即可)
- ❌ 性能优化 (先保证功能正确)

---

## 📅 Day 1: 基础设施 + 核心模块

### Hour 1-2: 核心数据结构 ⏰ 2小时

**任务**:
- [ ] 创建 `immunex/core/context.py`
  - [ ] `PipelineContext` dataclass (核心字段)
  - [ ] `get_output_path()` 方法
  - [ ] `to_dict()` 序列化
  - [ ] `copy()` 深拷贝

- [ ] 创建 `immunex/core/exceptions.py`
  - [ ] `ImmunexError` 基类
  - [ ] `InputValidationError`
  - [ ] `PipelineError`

**产出**:
```
immunex/core/
├── context.py        # ~150行
└── exceptions.py     # ~50行
```

**验收**: 可以创建PipelineContext实例并序列化

---

### Hour 3-4: PipelineNode基类 + 重构1个模块 ⏰ 2小时

**任务**:
- [ ] 创建 `immunex/core/base_node.py`
  - [ ] `PipelineNode` 抽象基类
  - [ ] `execute()` 抽象方法
  - [ ] `validate_inputs()` 方法

- [ ] 重构 `immunex/analysis/trajectory/rmsd.py`
  - [ ] 创建 `RMSDInput` dataclass
  - [ ] 创建 `RMSDResult` dataclass
  - [ ] 重构 `RMSDCalculator.calculate()` 方法
  - [ ] 输入验证

**产出**:
```
immunex/core/
└── base_node.py      # ~80行

immunex/analysis/trajectory/
└── rmsd.py           # 重构，添加Input/Result类
```

**验收**: RMSDCalculator可以接收RMSDInput并返回RMSDResult

---

### Hour 5-6: 创建Pipeline Nodes ⏰ 2小时

**任务**:
- [ ] 创建 `immunex/pipeline/nodes/`目录

- [ ] 创建 `immunex/pipeline/nodes/rmsd_node.py`
  - [ ] `RMSDNode` 类
  - [ ] `execute()` 实现
  - [ ] 调用 `RMSDCalculator`
  - [ ] 更新 `context.results['rmsd']`

- [ ] 创建 `immunex/pipeline/nodes/preprocess_node.py`
  - [ ] `PreprocessNode` 类
  - [ ] 调用 `PBCProcessor`
  - [ ] 更新 `context.trajectory_processed`

**产出**:
```
immunex/pipeline/
├── __init__.py
└── nodes/
    ├── __init__.py
    ├── preprocess_node.py   # ~100行
    └── rmsd_node.py         # ~80行
```

**验收**: 节点可以接收context，调用底层模块，返回更新后的context

---

### Hour 7-8: Pipeline基类 + 示例Pipeline ⏰ 2小时

**任务**:
- [ ] 创建 `immunex/pipeline/base_pipeline.py`
  - [ ] `Pipeline` 基类
  - [ ] `execute()` 方法（顺序执行节点）
  - [ ] 错误处理
  - [ ] `should_stop` 控制

- [ ] 创建 `immunex/pipeline/standard_pipelines.py`
  - [ ] `StandardTrajectoryPipeline` 类
  - [ ] 定义节点列表
  - [ ] 实现 `execute()`

**产出**:
```
immunex/pipeline/
├── base_pipeline.py        # ~120行
└── standard_pipelines.py   # ~100行
```

**验收**: 可以创建Pipeline实例并执行单个任务

---

## 📅 Day 2: BatchExecutor + 集成测试

### Hour 1-2: TaskDiscovery工具 ⏰ 2小时

**任务**:
- [ ] 创建 `immunex/core/task_discovery.py`
  - [ ] `TaskDiscovery` 类
  - [ ] `discover()` 方法
  - [ ] 默认发现规则 (查找md.xtc + md.tpr)
  - [ ] `_create_context()` 方法

**产出**:
```
immunex/core/
└── task_discovery.py    # ~150行
```

**验收**: 可以扫描目录并返回PipelineContext列表

---

### Hour 3-4: BatchExecutor ⏰ 2小时

**任务**:
- [ ] 创建 `immunex/pipeline/batch_executor.py`
  - [ ] `BatchExecutor` 类
  - [ ] `execute_pipeline()` 方法
  - [ ] 并行执行 (`ProcessPoolExecutor`)
  - [ ] 错误隔离
  - [ ] `summarize_results()` 统计

**产出**:
```
immunex/pipeline/
└── batch_executor.py    # ~200行
```

**验收**: 可以批量执行Pipeline并返回结果汇总

---

### Hour 5-6: 集成测试 + 示例脚本 ⏰ 2小时

**任务**:
- [ ] 创建测试脚本 `tests/test_integration_2day.py`
  - [ ] 测试单任务Pipeline
  - [ ] 测试批量执行
  - [ ] 测试错误处理

- [ ] 创建示例脚本 `examples/quick_start_new_architecture.py`
  - [ ] 单任务示例
  - [ ] 批量处理示例
  - [ ] 详细注释

**产出**:
```
tests/
└── test_integration_2day.py    # ~150行

examples/
└── quick_start_new_architecture.py    # ~100行
```

**验收**: 示例可成功运行，输出正确

---

### Hour 7-8: 文档 + 兼容性包装 ⏰ 2小时

**任务**:
- [ ] 创建 `docs/architecture/NEW_ARCHITECTURE_QUICKSTART.md`
  - [ ] 架构简介
  - [ ] 使用示例
  - [ ] API文档
  - [ ] 迁移指南

- [ ] 创建兼容性包装示例 `scripts/batch_rmsd_new.py`
  - [ ] 使用新架构
  - [ ] 保持原有CLI接口
  - [ ] 添加对比说明

**产出**:
```
docs/
└── NEW_ARCHITECTURE_QUICKSTART.md    # ~500行

scripts/
└── batch_rmsd_new.py    # ~80行 (示例)
```

**验收**: 文档清晰，兼容性脚本可运行

---

## 📦 最终交付清单

### 代码文件 (约1500行新代码)
```
immunex/core/
├── context.py              # PipelineContext
├── exceptions.py           # 异常类
├── base_node.py            # PipelineNode基类
└── task_discovery.py       # TaskDiscovery

immunex/analysis/trajectory/
└── rmsd.py                 # 重构版 (添加Input/Result)

immunex/pipeline/
├── base_pipeline.py        # Pipeline基类
├── standard_pipelines.py   # 示例Pipeline
├── batch_executor.py       # BatchExecutor
└── nodes/
    ├── preprocess_node.py  # PreprocessNode
    └── rmsd_node.py        # RMSDNode

examples/
└── quick_start_new_architecture.py    # 使用示例

scripts/
└── batch_rmsd_new.py       # 兼容性示例

tests/
└── test_integration_2day.py    # 集成测试

docs/
└── NEW_ARCHITECTURE_QUICKSTART.md    # 快速开始文档
```

### 功能验收
- [ ] 可以创建PipelineContext并传递数据
- [ ] 可以创建自定义Pipeline
- [ ] 可以批量执行Pipeline
- [ ] 可以正确处理错误
- [ ] 可以自动发现任务
- [ ] 示例脚本可成功运行
- [ ] 输出结果正确

---

## 🚀 执行策略

### 优先级排序
```
P0 (必须完成):
- PipelineContext
- PipelineNode基类
- RMSDCalculator重构
- 1个Pipeline示例
- BatchExecutor基础功能
- TaskDiscovery
- 使用示例

P1 (尽量完成):
- 错误处理
- 基础测试
- 文档

P2 (可延后):
- 完整测试覆盖
- 性能优化
- 更多示例
```

### 时间管理
- **严格控制每个时间块** - 每2小时完成1个模块
- **快速迭代** - 先实现核心功能，细节后续优化
- **持续验收** - 每个模块完成后立即验证
- **果断取舍** - 如果某个功能耗时过长，果断简化

### 风险控制
**如果进度落后**:
- Hour 5-6: 砍掉PreprocessNode，只保留RMSDNode
- Hour 7-8: 简化Pipeline，只做最基本的顺序执行
- Day 2: 砍掉集成测试，只保留示例脚本
- Day 2: 简化文档，只写核心API

**底线交付** (最小化):
```
必须有:
1. PipelineContext (可传递数据)
2. 1个Node (RMSDNode)
3. 1个Pipeline (3行代码串起来)
4. BatchExecutor (并行执行)
5. 1个示例脚本 (证明可用)

这样至少证明了架构可行
```

---

## 📊 进度追踪

### Day 1
- [ ] Hour 1-2: ⬜ PipelineContext + 异常类
- [ ] Hour 3-4: ⬜ PipelineNode基类 + RMSD重构
- [ ] Hour 5-6: ⬜ 创建Nodes
- [ ] Hour 7-8: ⬜ Pipeline基类 + 示例Pipeline

### Day 2
- [ ] Hour 1-2: ⬜ TaskDiscovery
- [ ] Hour 3-4: ⬜ BatchExecutor
- [ ] Hour 5-6: ⬜ 集成测试 + 示例
- [ ] Hour 7-8: ⬜ 文档 + 兼容性包装

### 完成标准
- [ ] 示例脚本可成功运行
- [ ] 批量处理功能正常
- [ ] 文档完整
- [ ] 代码可读可维护

---

## 🎯 Day 2结束时的效果

**用户可以做什么**:
```python
# 方式1: 使用预定义Pipeline
from immunex.pipeline import StandardTrajectoryPipeline, BatchExecutor
from immunex.core import TaskDiscovery

# 发现任务
discovery = TaskDiscovery()
tasks = discovery.discover("/data/simulations")

# 批量执行
executor = BatchExecutor(max_workers=4)
results = executor.execute_pipeline(tasks, StandardTrajectoryPipeline())

# 查看结果
summary = executor.summarize_results(results)
print(f"成功: {summary['successful']}/{summary['total_tasks']}")

# 方式2: 自定义Pipeline
from immunex.pipeline import Pipeline
from immunex.pipeline.nodes import RMSDNode, PreprocessNode

class MyPipeline(Pipeline):
    nodes = [
        PreprocessNode(),
        RMSDNode(selection="protein")
    ]

# 执行
results = executor.execute_pipeline(tasks, MyPipeline())
```

**架构已经可用**:
- ✅ 数据通过Context传递，不再是字符串路径
- ✅ 可以自由组合Pipeline
- ✅ 批处理逻辑统一，不需要重复写
- ✅ 职责清晰，易于扩展

**后续可以做**:
- 渐进式迁移其他模块
- 添加更多Nodes
- 完善测试
- 添加YAML配置支持
- 添加CLI命令

---

## 💪 开始行动

**立即开始**:
```bash
# 1. 创建目录结构
mkdir -p immunex/core immunex/pipeline/nodes tests examples

# 2. 开始编写第一个文件
# immunex/core/context.py
```

**我会在这2天内持续提供支持**:
- 代码实现
- 问题解答
- 进度检查
- 及时调整

让我们开始吧！🚀
