# Preprocess and Quality Guide

本指南描述当前活动主线中的前处理与质量评估流程。旧的大一统前处理 pipeline 和脚本式入口已经下线，当前统一使用：

- `imn preprocess`
- `imn quality`
- `imn batch preprocess`
- `PreprocessOnlyPipeline`
- `PreprocessQualityPipeline`
- `QualityAssessmentPipeline`

## 推荐工作流

### 单体系：前处理 + 预处理质量检查

```bash
imn preprocess \
  --trajectory /path/to/md.xtc \
  --topology /path/to/md.tpr \
  --output-dir ./output/system_a
```

默认会输出：

- `md_processed.xtc`
- `md_processed_converted.pdb`
- `analysis/rmsd/rmsd.xvg`
- `quality` 相关汇总文件

### 单体系：已有预处理轨迹时只做质量评估

```bash
imn quality \
  --trajectory /path/to/md_processed.xtc \
  --topology /path/to/md.tpr \
  --output-dir ./output/system_a_quality
```

### 批量前处理

```bash
imn batch preprocess input/database/parallel2 \
  --workers 4 \
  --output-dir output/batch_preprocess
```

## Python API

### 仅前处理

```python
from immunex.core.context import PipelineContext
from immunex.pipeline import PreprocessOnlyPipeline

context = PipelineContext(
    system_id="1oga_run2",
    topology="md.tpr",
    trajectory_raw="md.xtc",
    output_dir="output/1oga_run2",
)

pipeline = PreprocessOnlyPipeline(method="2step")
pipeline.execute(context)
```

### 前处理 + 质量检查

```python
from immunex.core.context import PipelineContext
from immunex.pipeline import PreprocessQualityPipeline

context = PipelineContext(
    system_id="1oga_run2",
    topology="md.tpr",
    trajectory_raw="md.xtc",
    output_dir="output/1oga_run2",
)

pipeline = PreprocessQualityPipeline(method="2step", rmsd_selection="backbone")
pipeline.execute(context)
```

### 已处理轨迹的质量评估

```python
from immunex.pipeline import QualityAssessmentPipeline

pipeline = QualityAssessmentPipeline()
results = pipeline.run_comprehensive_assessment(
    trajectory="md_processed.xtc",
    topology="md.tpr",
    rmsd_selection="backbone",
    output_dir="output/quality_only",
)
```

## 输出约定

当前主线输出遵循“按用途分层”的结构：

```text
output/system_a/
├── md_processed.xtc
├── md_processed_converted.pdb
├── analysis/
│   └── rmsd/
│       └── rmsd.xvg
└── quality/
    ├── quality_summary.json
    └── convergence_report.txt
```

实际字段和子目录可能随节点能力扩展，但入口与主线命名保持以上约定。

## 说明

- `2step` 仍是当前推荐的默认前处理方法。
- `QualityAssessmentPipeline` 仅负责“已处理轨迹”的质量评估，不负责 PBC 前处理。
- 如果你在旧文档、旧脚本中看到旧版前处理 pipeline，请将其视为历史实现，不再作为当前主线使用。
