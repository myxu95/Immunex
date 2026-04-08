# PBC-RMSD Quality Pipeline 使用指南

本指南介绍如何使用Immunex的前处理+质量评估管线。

## 概述

PBC-RMSD管线提供了从原始MD轨迹到质量评估的完整工作流:

1. **PBC校正** - 移除周期性边界条件伪影
2. **Post-PBC验证** - 验证PBC处理质量
3. **RMSD计算** - 计算骨架RMSD
4. **收敛性分析** - 分析轨迹收敛性和稳定性
5. **质量评分** - 分配A/B/C/D等级

## 快速开始

### 单轨迹处理

```bash
python scripts/pbc_rmsd_workflow.py \
    --trajectory /path/to/md.xtc \
    --topology /path/to/md.tpr \
    --output-dir ./processed
```

### 批量处理

```bash
python scripts/batch_quality_pipeline.py \
    --base-dir /data/md_simulations \
    --output-dir /data/processed \
    --max-workers 4 \
    --enable-quality-check
```

## Python API使用

### 单轨迹处理

```python
from immunex.pipeline import PBCRMSDPipeline

# 初始化管线
pipeline = PBCRMSDPipeline(gmx_executable="gmx")

# 处理单个轨迹
results = pipeline.process_single_trajectory(
    trajectory="md.xtc",
    topology="md.tpr",
    output_dir="./processed",
    pbc_method="2step",
    rmsd_selection="backbone",
    run_quality_check=True
)

# 查看结果
print(f"质量等级: {results['overall_grade']}")
print(f"是否合格: {results['is_qualified']}")
print(f"RMSD均值: {results['rmsd_metrics']['mean_rmsd']:.3f} nm")
```

### 批量处理

```python
from immunex.pipeline import PBCRMSDPipeline

# 定义任务列表
tasks = [
    {
        'trajectory': 'task1/md.xtc',
        'topology': 'task1/md.tpr',
        'output_dir': 'task1/processed'
    },
    {
        'trajectory': 'task2/md.xtc',
        'topology': 'task2/md.tpr',
        'output_dir': 'task2/processed'
    }
]

# 批量处理
pipeline = PBCRMSDPipeline()
batch_results = pipeline.batch_process(
    tasks=tasks,
    max_workers=4
)

# 查看汇总
print(f"成功: {batch_results['successful']}")
print(f"失败: {batch_results['failed']}")
print(f"等级分布: {batch_results['summary']}")
```

### 仅质量评估(已有PBC处理的轨迹)

```python
from immunex.pipeline import QualityAssessmentPipeline

# 初始化质量评估管线
pipeline = QualityAssessmentPipeline()

# 运行质量评估
results = pipeline.run_comprehensive_assessment(
    trajectory="md_pbc.xtc",
    topology="md.tpr",
    rmsd_selection="backbone",
    output_dir="./quality_reports"
)

# 生成报告
pipeline.generate_quality_report(
    results=results,
    output_file="quality_report.md",
    format="markdown"
)
```

## 质量评分标准

### Post-PBC验证

- **质心漂移** < 1.0 nm (可配置)
- **Rg标准差/均值** < 0.15 (可配置)
- **帧间RMSD跳跃** < 0.5 nm (可配置)

### RMSD收敛性分级

| 等级 | 均值RMSD | 标准差 | 收敛时间 |
|------|----------|--------|----------|
| A (优秀) | < 0.3 nm | < 0.05 nm | < 20% |
| B (良好) | 0.3-0.5 nm | 0.05-0.1 nm | < 30% |
| C (合格) | 0.5-0.8 nm | 0.1-0.15 nm | < 50% |
| D (不合格) | > 0.8 nm | > 0.15 nm | 未收敛 |

### 综合评分

最终等级取Post-PBC验证和RMSD分析中较差的一个。

## CLI参数详解

### pbc_rmsd_workflow.py

**必需参数:**
- `--trajectory`: 输入轨迹文件
- `--topology`: 拓扑文件
- `--output-dir`: 输出目录

**可选参数:**
- `--pbc-method`: PBC方法(`2step`或`3step`, 默认`2step`)
- `--dt`: 输出轨迹时间步长(ps)
- `--rmsd-selection`: RMSD选择字符串(默认`backbone`)
- `--no-quality-check`: 跳过质量检查
- `--validation-stride`: 验证采样间隔(默认1)
- `--max-com-drift`: 最大质心漂移(nm, 默认1.0)
- `--max-rg-std-ratio`: 最大Rg标准差比(默认0.15)
- `--no-report`: 不生成报告
- `--verbose`: 详细日志

### batch_quality_pipeline.py

**必需参数:**
- `--base-dir`: 包含MD任务的基础目录

**可选参数:**
- `--output-dir`: 覆盖输出目录
- `--max-workers`: 最大并行数(默认4)
- `--pbc-method`: PBC方法(默认`2step`)
- `--enable-quality-check`: 启用质量检查
- `--config`: 配置文件(YAML)
- `--summary-file`: 保存汇总JSON文件
- `--fail-fast`: 遇到失败立即停止

## 配置文件

使用YAML配置文件可以避免重复输入参数:

```yaml
# scripts/configs/default_pbc_quality.yaml

pipeline:
  pbc_method: "2step"
  dt: 100.0
  gmx_executable: "gmx"

quality_check:
  enable_post_pbc_validation: true
  enable_rmsd_convergence: true
  rmsd_selection: "backbone"

  post_pbc_thresholds:
    max_com_drift: 1.0
    max_rg_std_ratio: 0.15
    max_frame_jump_rmsd: 0.5

batch:
  max_workers: 4
  retry_failed: true
```

使用配置文件:

```bash
python scripts/batch_quality_pipeline.py \
    --config scripts/configs/default_pbc_quality.yaml \
    --base-dir /data/md_simulations
```

## 输出文件

### 单轨迹输出

```
output_dir/
├── md_pbc.xtc                      # PBC校正后的轨迹
├── rmsd.csv                        # RMSD时间序列
├── md_quality_metrics.json         # 质量指标(JSON)
├── md_quality_report.md            # 质量报告(Markdown)
└── md_convergence.txt              # 收敛性报告(文本)
```

### 批量处理输出

```
base_dir/
├── task1/
│   └── processed/
│       ├── md_pbc.xtc
│       ├── rmsd.csv
│       └── ...
├── task2/
│   └── processed/
│       └── ...
└── summary.json                    # 批量处理汇总
```

## 常见问题

### Q: 什么时候使用2步法vs 3步法?

**推荐使用2步法** (默认):
- 适用于大多数TCR-pMHC复合物
- 更快,问题更少
- 基于实践经验优化

**使用3步法**:
- 标准的GROMACS工作流
- 需要显式居中和完整化
- 某些特殊体系可能需要

### Q: 如何调整质量评分阈值?

通过配置文件或命令行参数:

```bash
python scripts/pbc_rmsd_workflow.py \
    --trajectory md.xtc \
    --topology md.tpr \
    --output-dir ./processed \
    --max-com-drift 1.5 \
    --max-rg-std-ratio 0.2
```

### Q: 批量处理时如何处理失败的任务?

默认情况下,单个任务失败不影响其他任务:

```python
batch_results = pipeline.batch_process(
    tasks=tasks,
    max_workers=4,
    fail_fast=False  # 继续处理其他任务
)

# 检查失败任务
for task_result in batch_results['task_results']:
    if task_result['status'] == 'failed':
        print(f"失败: {task_result['task_name']}")
        print(f"错误: {task_result['error']}")
```

## 进阶使用

### 自定义Post-PBC验证

```python
from immunex.analysis.quality import PostPBCValidator

# 自定义阈值
validator = PostPBCValidator(
    max_com_drift=1.5,
    max_rg_std_ratio=0.2,
    max_frame_jump_rmsd=0.8
)

# 运行验证
results = validator.comprehensive_validation(
    trajectory="md_pbc.xtc",
    topology="md.tpr",
    selection="protein",
    stride=10  # 每10帧采样一次
)
```

### 自定义RMSD收敛性分析

```python
from immunex.analysis.trajectory import RMSDConvergenceAnalyzer

analyzer = RMSDConvergenceAnalyzer(
    n_blocks=10,              # 分成10个block
    moving_avg_window=200,    # 200帧移动平均
    convergence_threshold=0.03 # 3 Å收敛阈值
)

metrics = analyzer.calculate_convergence_metrics(times, rmsd_values)
grade = analyzer.assign_quality_grade(metrics)
```

## 相关文档

- [模块设计文档](../ARCHITECTURE.md)
- [PBC处理指南](../docs/PBC_PROCESSING_METHODS.md)
- [RMSD计算指南](../docs/RMSD_CALCULATION_GUIDE.md)
