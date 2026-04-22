# 角度计算脚本梳理

这份文档用于说明 Immunex 当前“结合角度分析”链路涉及哪些脚本，以及它们分别负责什么。

## 总体结构

当前角度计算不是单一脚本，而是分成 4 层：

1. `analysis/angles`
   负责真正的角度算法、几何定义、输入输出数据结构。
2. `pipeline/nodes`
   负责把角度分析包装成一个可插入 pipeline 的步骤。
3. `pipeline`
   负责把链识别和角度分析编排成完整流程。
4. `cli/commands`
   负责把批量角度分析暴露给 CLI。

## 核心文件

### 1. 算法核心

- [analyzer.py](/home/xumy/work/development/Immunex/immunex/analysis/angles/analyzer.py)

这是当前角度模块最核心的文件，主要负责：

- 读取拓扑和轨迹
- 识别 TCR 参考半胱氨酸
- 识别 MHC groove 区域
- 构建 `TCR axis / groove axis / groove normal`
- 逐帧计算 `Crossing / Incident`
- 输出 `docking_angles.csv`

当前这个文件内部还包含几组关键内部类：

- `_GeometryUtils`
  负责平面拟合、向量夹角、角度规范化。
- `_TCRAxisCalculator`
  负责保守半胱氨酸检测与 `TCR axis` 构建。
- `_MHCGrooveDetector`
  负责 groove 区域识别。
- `_StructureGrooveFallback`
  负责当序列映射不可靠时，直接从结构中寻找 groove 双壁。
- `_SequenceAligner`
  负责旧的参考序列映射逻辑。

### 2. 输入输出结构

- [angle_data_structures.py](/home/xumy/work/development/Immunex/immunex/analysis/angles/angle_data_structures.py)

负责定义：

- `DockingAngleInput`
- `DockingAngleResult`

这里是分析模块和 pipeline/node 之间的正式数据契约。

### 3. 分析模块导出

- [__init__.py](/home/xumy/work/development/Immunex/immunex/analysis/angles/__init__.py)

负责对外导出角度分析模块的公开接口。

## Pipeline 层

### 4. 角度分析节点

- [docking_angle_node.py](/home/xumy/work/development/Immunex/immunex/pipeline/nodes/docking_angle_node.py)

这是 pipeline 里的角度节点，负责：

- 从 `PipelineContext` 取输入
- 构造 `DockingAngleInput`
- 调用 `DockingAngleAnalyzer`
- 把结果写回 `context.results['docking_angles']`
- 把元数据写回 `context.metadata['angle_analysis']`

### 5. 角度分析流程

- [analysis_pipelines.py](/home/xumy/work/development/Immunex/immunex/pipeline/analysis_pipelines.py)

这里定义了：

- `DockingAnglePipeline`

当前流程是：

1. `ChainIdentificationNode`
2. `DockingAngleNode`

也就是说，角度分析默认会先做链识别，再做角度计算。

## CLI 层

### 6. 批量角度分析入口

- [batch.py](/home/xumy/work/development/Immunex/immunex/cli/commands/batch.py)

当前 CLI 注册入口是：

```bash
imn batch angle <base_dir>
```

这个文件负责：

- 注册 `batch angle`
- 发现任务
- 调用 `DockingAnglePipeline`
- 汇总批量结果
- 生成：
  - `docking_angle_variation_summary.csv`
  - `docking_angle_variation_summary.json`
  - `docking_angle_variation_summary.md`
  - 排序图

当前还支持：

- `--stride`
- `--limit`
- `--print-each-frame`
- `--workers`

## 测试文件

### 7. 批量汇总与节点测试

- [test_batch_angle_summary.py](/home/xumy/work/development/Immunex/tests/test_batch_angle_summary.py)

当前覆盖的重点是：

- 角度时间序列汇总
- 批量 summary 产物
- `DockingAngleNode` 输入构造行为
- 链映射到手工 selection 的逻辑

### 8. ANARCI 半胱氨酸检测相关测试

- [test_anarci_cysteine_detection.py](/home/xumy/work/development/Immunex/tests/test_anarci_cysteine_detection.py)

主要覆盖 `TCR` 保守半胱氨酸检测这一层。

## 相关但不是主链的文件

### 9. 旧/旁路线实现

- [complex_angle_analyzer.py](/home/xumy/work/development/Immunex/immunex/analysis/trajectory/complex_angle_analyzer.py)

这个文件和角度分析有关，但**不是当前正式 CLI 与 pipeline 主路径**。当前主路径仍然是 `analysis/angles/analyzer.py`。

## 当前主路径

如果只看“现在真正会被调用的主链”，顺序是：

1. [batch.py](/home/xumy/work/development/Immunex/immunex/cli/commands/batch.py)
2. [analysis_pipelines.py](/home/xumy/work/development/Immunex/immunex/pipeline/analysis_pipelines.py)
3. [docking_angle_node.py](/home/xumy/work/development/Immunex/immunex/pipeline/nodes/docking_angle_node.py)
4. [angle_data_structures.py](/home/xumy/work/development/Immunex/immunex/analysis/angles/angle_data_structures.py)
5. [analyzer.py](/home/xumy/work/development/Immunex/immunex/analysis/angles/analyzer.py)

## 当前最值得继续优化的点

如果后续继续改角度算法，优先关注：

1. [analyzer.py](/home/xumy/work/development/Immunex/immunex/analysis/angles/analyzer.py)
   原因：几何定义和逐帧稳定性都在这里。
2. [docking_angle_node.py](/home/xumy/work/development/Immunex/immunex/pipeline/nodes/docking_angle_node.py)
   原因：这里决定参数如何进入分析器。
3. [batch.py](/home/xumy/work/development/Immunex/immunex/cli/commands/batch.py)
   原因：这里决定批量运行和结果汇总怎么对外暴露。
