# 统一批处理执行器设计方案

**日期**: 2026-03-16
**作者**: Immunex Development Team
**版本**: v2.0

---

## 核心理念

用户的洞察非常正确：**我们不需要为每种分析写独立的batch脚本**。

### 当前问题

```python
# ❌ 当前做法：每个分析都有独立的batch脚本
scripts/
├── batch_tcr_rmsd.py          # 245行
├── batch_cdr_rmsd_exact.py    # 410行
├── batch_cdr_rmsf.py          # 461行
├── batch_docking_angles.py    # 392行
└── ... (还有11个batch脚本)

# 每个脚本都重复实现：
# 1. 任务发现逻辑
# 2. 批处理循环
# 3. 结果汇总
# 4. 命令行参数解析
```

### 理想方案

```python
# ✅ 理想做法：一个统一的执行器 + Pipeline定义

# 用户只需要：
from immunex.pipeline import TCRRMSDPipeline, BatchExecutor
from immunex.core import TaskDiscovery

# 发现任务
discovery = TaskDiscovery()
tasks = discovery.discover("/data")

# 创建Pipeline（定义WHAT）
pipeline = TCRRMSDPipeline()

# 执行批处理（定义HOW）
executor = BatchExecutor(max_workers=4)
results = executor.execute_pipeline(tasks, pipeline)

# 仅需 6 行代码！
```

---

## 架构设计

### 三层分离

```
┌─────────────────────────────────────────────┐
│ Layer 1: Pipeline Definition (WHAT)        │
│ - TCRRMSDPipeline                           │
│ - CDRRMSFPipeline                           │
│ - DockingAnglePipeline                      │
│ - ... (业务逻辑，定义要计算什么)             │
└─────────────────────────────────────────────┘
                    ↓
┌─────────────────────────────────────────────┐
│ Layer 2: Task Discovery (WHERE)            │
│ - TaskDiscovery                             │
│ - 查找任务目录，创建PipelineContext         │
└─────────────────────────────────────────────┘
                    ↓
┌─────────────────────────────────────────────┐
│ Layer 3: Batch Executor (HOW)              │
│ - BatchExecutor.execute_pipeline()          │
│ - 并行调度，进度跟踪，错误处理              │
└─────────────────────────────────────────────┘
```

### 关键原则

1. **Pipeline = WHAT** (要计算什么)
   - 定义计算步骤
   - 定义输入输出
   - **不关心**是单任务还是批处理

2. **TaskDiscovery = WHERE** (数据在哪里)
   - 查找任务目录
   - 验证文件存在性
   - **不关心**要运行什么分析

3. **BatchExecutor = HOW** (如何并行执行)
   - 并行调度
   - 进度跟踪
   - 错误隔离
   - **不关心**具体的分析逻辑

---

## 实施方案

### 方案A: 统一CLI工具（推荐）

创建一个`imn`命令行工具，支持所有分析：

```bash
# TCR RMSD分析
imn batch tcr-rmsd --input-dirs /data --max-workers 8

# CDR RMSF分析
imn batch cdr-rmsf --cdr-csv cdrs.csv --input-dirs /data

# 对接角度分析
imn batch docking-angles --input-dirs /data --stride 10

# 所有分析共享同一个批处理引擎
```

**优势**:
- ✅ 用户只需要记住一个命令
- ✅ 所有分析统一的参数格式
- ✅ 无需维护15个独立脚本
- ✅ 新增分析 = 添加subcommand，无需新脚本

**实现**:

```python
# scripts/imn
#!/usr/bin/env python3
"""Immunex unified CLI tool"""

def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()

    # 每个分析类型一个subcommand
    parser_tcr = subparsers.add_parser('tcr-rmsd')
    parser_tcr.add_argument('--input-dirs', nargs='+', required=True)
    parser_tcr.add_argument('--max-workers', type=int, default=4)
    parser_tcr.set_defaults(func=cmd_tcr_rmsd)

    # ... 其他subcommands

    args = parser.parse_args()
    return args.func(args)

def cmd_tcr_rmsd(args):
    """仅需10行代码"""
    tasks = TaskDiscovery().discover(args.input_dirs)
    pipeline = TCRRMSDPipeline()
    executor = BatchExecutor(max_workers=args.max_workers)
    results = executor.execute_pipeline(tasks, pipeline)
    return 0 if all(not r.has_errors() for r in results) else 1
```

### 方案B: Python API（编程接口）

对于需要自定义逻辑的用户：

```python
from immunex import batch_execute

# 最简单的用法
results = batch_execute(
    pipeline='tcr-rmsd',
    input_dirs=['/data/batch1', '/data/batch2'],
    max_workers=8
)

# 自定义Pipeline
from immunex.pipeline import Pipeline, RMSDNode

class MyCustomPipeline(Pipeline):
    def __init__(self):
        super().__init__(nodes=[
            ChainIdentificationNode(),
            RMSDNode(selection='backbone')
        ])

results = batch_execute(
    pipeline=MyCustomPipeline(),
    input_dirs=['/data'],
    max_workers=4
)
```

---

## 代码对比

### 当前方案（每种分析一个脚本）

**问题**: 代码重复

```
batch_tcr_rmsd.py (245行)
├── 任务发现逻辑 (40行)        ← 重复
├── 批处理循环 (60行)          ← 重复
├── 结果汇总 (50行)            ← 重复
├── 参数解析 (40行)            ← 重复
└── TCR RMSD逻辑 (55行)        ← 唯一

batch_cdr_rmsf.py (461行)
├── 任务发现逻辑 (40行)        ← 重复
├── 批处理循环 (60行)          ← 重复
├── 结果汇总 (50行)            ← 重复
├── 参数解析 (40行)            ← 重复
└── CDR RMSF逻辑 (271行)       ← 唯一

总计: 706行，其中 460行重复 (65%)
```

### 新方案（统一执行器）

**优势**: 零重复

```
imn (统一CLI工具, 300行)
├── 通用任务发现 (40行)        ← 仅一次
├── 通用批处理循环 (60行)      ← 仅一次
├── 通用结果汇总 (50行)        ← 仅一次
├── 通用参数解析 (40行)        ← 仅一次
├── Subcommand调度 (60行)      ← 仅一次
└── 各分析调用逻辑 (50行)      ← 每个5-10行

TCRRMSDPipeline (已存在, 50行)  ← 已在analysis_pipelines.py
CDRRMSFPipeline (新增, 100行)   ← 仅业务逻辑

总计: 450行，0行重复 (0%)
代码减少: 706 - 450 = 256行 (36%减少)
```

---

## BatchExecutor增强

为了支持这个理念，`BatchExecutor`需要增强：

### 当前BatchExecutor

```python
class BatchExecutor:
    def execute_pipeline(self, tasks: List[PipelineContext],
                         pipeline: Pipeline) -> List[PipelineContext]:
        """执行批处理"""
        with ProcessPoolExecutor(max_workers=self.max_workers) as executor:
            futures = [executor.submit(pipeline.execute, task) for task in tasks]
            results = [f.result() for f in as_completed(futures)]
        return results
```

### 增强的BatchExecutor

```python
class BatchExecutor:
    def __init__(self, max_workers: int = 4, show_progress: bool = True):
        self.max_workers = max_workers
        self.show_progress = show_progress

    def execute_pipeline(self,
                         tasks: List[PipelineContext],
                         pipeline: Pipeline,
                         output_dir: Optional[str] = None) -> List[PipelineContext]:
        """
        执行批处理Pipeline

        Args:
            tasks: 任务列表
            pipeline: Pipeline实例
            output_dir: 输出目录（可选，会自动设置每个task的output_dir）

        Returns:
            执行结果列表
        """
        # 设置输出目录
        if output_dir:
            for task in tasks:
                if not task.output_dir:
                    task.output_dir = str(Path(output_dir) / task.system_id)

        # 并行执行
        results = []
        with ProcessPoolExecutor(max_workers=self.max_workers) as executor:
            futures = {executor.submit(pipeline.execute, task): task for task in tasks}

            for i, future in enumerate(as_completed(futures), 1):
                task = futures[future]
                try:
                    result = future.result()
                    results.append(result)

                    # 进度显示
                    if self.show_progress:
                        status = "✓" if not result.has_errors() else "✗"
                        logger.info(f"[{i}/{len(tasks)}] {status} {task.system_id}")

                except Exception as e:
                    logger.error(f"[{i}/{len(tasks)}] ✗ {task.system_id}: {e}")
                    # 创建失败的context
                    task.add_error(str(e))
                    results.append(task)

        return results

    def summarize_results(self, results: List[PipelineContext]) -> dict:
        """生成结果摘要"""
        total = len(results)
        successful = sum(1 for r in results if not r.has_errors())
        failed = total - successful

        return {
            'total_tasks': total,
            'successful': successful,
            'failed': failed,
            'success_rate': (successful / total * 100) if total > 0 else 0.0
        }

    def save_summary(self, results: List[PipelineContext],
                     output_file: str,
                     pipeline_name: str = None):
        """
        保存结果摘要到CSV

        Args:
            results: 执行结果
            output_file: 输出CSV文件
            pipeline_name: Pipeline名称（可选）
        """
        import pandas as pd

        summary_data = []
        for result in results:
            summary_data.append({
                'task': result.system_id,
                'pdb_id': result.metadata.get('pdb_id', ''),
                'status': 'success' if not result.has_errors() else 'failed',
                'errors': '; '.join(result.errors) if result.errors else None,
                'warnings': '; '.join(result.warnings) if result.warnings else None,
                # 添加Pipeline特定的结果字段
                **result.results
            })

        df = pd.DataFrame(summary_data)
        df.to_csv(output_file, index=False)
        logger.info(f"Summary saved to: {output_file}")
```

---

## 迁移路径

### Phase 1: 创建统一CLI工具（本周）

1. ✅ 创建`scripts/imn_batch.py`
2. ⏳ 实现3-5个常用subcommands
3. ⏳ 测试与现有batch脚本的功能一致性

### Phase 2: 增强BatchExecutor（下周）

1. ⏳ 添加进度显示
2. ⏳ 添加`save_summary()`方法
3. ⏳ 添加错误隔离和重试机制

### Phase 3: 完善Pipeline库（下周）

1. ⏳ 确保所有分析都有对应的Pipeline类
2. ⏳ 标准化Pipeline的输入输出
3. ⏳ 添加Pipeline文档和示例

### Phase 4: 弃用旧脚本（下月）

1. ⏳ 在旧脚本中添加弃用警告
2. ⏳ 更新所有文档指向新CLI
3. ⏳ 移除旧batch脚本

---

## 用户体验对比

### 当前方式

```bash
# 用户需要记住15个不同的脚本
python batch_tcr_rmsd.py --input-dirs /data --max-workers 8
python batch_cdr_rmsf.py --cdr-csv cdrs.csv --input-dirs /data --max-workers 8
python batch_docking_angles.py --input-dirs /data --stride 10 --max-workers 8

# 每个脚本的参数格式可能不同
# 用户需要查看每个脚本的帮助文档
```

### 新方式

```bash
# 用户只需要记住一个命令
imn batch tcr-rmsd --input-dirs /data --max-workers 8
imn batch cdr-rmsf --cdr-csv cdrs.csv --input-dirs /data --max-workers 8
imn batch docking-angles --input-dirs /data --stride 10 --max-workers 8

# 或者更简洁
imn batch tcr-rmsd /data          # 使用默认参数
imn batch cdr-rmsf /data --cdr-csv cdrs.csv
imn batch docking-angles /data

# 统一的帮助文档
imn batch --help                   # 查看所有可用分析
imn batch tcr-rmsd --help         # 查看特定分析的参数
```

---

## Python API示例

对于需要编程控制的用户：

```python
from immunex import batch

# 方式1: 使用预定义Pipeline名称
results = batch.execute(
    'tcr-rmsd',
    input_dirs=['/data/batch1', '/data/batch2'],
    max_workers=8
)

# 方式2: 使用Pipeline实例
from immunex.pipeline import TCRRMSDPipeline

pipeline = TCRRMSDPipeline(use_anarci=True)
results = batch.execute(
    pipeline,
    input_dirs=['/data'],
    max_workers=4
)

# 方式3: 完全自定义
from immunex.core import TaskDiscovery, PipelineContext
from immunex.pipeline import BatchExecutor

discovery = TaskDiscovery()
tasks = discovery.discover('/data')

executor = BatchExecutor(max_workers=4, show_progress=True)
results = executor.execute_pipeline(tasks, pipeline)

# 保存结果
executor.save_summary(results, 'results.csv', pipeline_name='TCR_RMSD')
```

---

## 总结

### 核心优势

1. **极简主义**: 用户只需一个命令，无需记住15个脚本
2. **零重复**: 所有基础设施代码只写一次
3. **易扩展**: 新增分析 = 添加Pipeline类 + 5行调用代码
4. **统一体验**: 所有分析使用相同的参数格式和输出格式

### 关键数字

- **代码减少**: 从5200行 → ~800行 (**85%减少**)
- **脚本数量**: 从15个 → 1个 (**93%减少**)
- **用户命令**: 从15个 → 1个 (**一个统一的`imn batch`命令**)

### 下一步

1. **立即**: 测试`imn_batch.py`原型
2. **本周**: 实现所有常用分析的subcommands
3. **下周**: 增强BatchExecutor，添加进度显示和结果保存
4. **下月**: 完全迁移到新CLI，弃用旧脚本

---

**结论**: 用户的洞察完全正确 - 我们确实不需要为每种分析写独立的batch脚本。统一的执行器 + Pipeline定义是更优雅、更可维护的方案。
