# Analysis Platform Rules

## 定位

Immunex 是一个用于构建分析流程的平台，而不是一组松散脚本的集合。

平台设计目标：
- 高复用分析元件可被多个流程共享
- pipeline 负责编排，而不是承载厚重业务实现
- 专用功能应尽量局部化，避免污染通用层

## 分层原则

### analysis 层

`immunex/analysis` 放真正做事的分析元件与语义元件。

适合放入 analysis 的内容：
- 轨迹分析计算器
- 结构分析计算器
- residue-pair 注释器
- 链识别、CDR 检测、语义分类等可复用能力
- 热图、矩阵、统计汇总等可复用结果构建器

判断标准：
- 这段代码是否在回答“如何计算”
- 这段代码是否不依赖具体 pipeline 顺序
- 这段代码是否可能在第二条或第三条流程中复用

如果答案是“是”，优先进入 `analysis`。

### node 层

`immunex/pipeline/nodes` 放流程步骤包装器。

适合放入 node 的内容：
- 对 analysis 元件的调用包装
- context 读写
- 输入校验
- 步骤级输出落盘
- 步骤级异常处理

判断标准：
- 这段代码是否在回答“流程中的一个步骤做什么”
- 这段代码是否需要和 `PipelineContext` 交互

如果答案是“是”，优先进入 `node`。

### pipeline 层

`immunex/pipeline` 放流程编排定义。

适合放入 pipeline 的内容：
- 节点顺序
- 默认参数组合
- 针对某条 workflow 的装配关系

不适合放入 pipeline 的内容：
- 复杂计算逻辑
- 厚重结果加工逻辑
- 细粒度语义判断逻辑
- 大量文件组织细节

判断标准：
- 这段代码是否只是在回答“先做什么，再做什么”

如果答案是“是”，可以进入 `pipeline`。

### support / 临时专用逻辑

对于复用性低、但又不适合直接塞进 pipeline 的逻辑：
- 优先放入专用 node
- 或放入对应领域下的 support/helper 模块

禁止把 pipeline 当作低复用逻辑的垃圾桶。

## 设计决策规则

新增功能时，按下面顺序判断归属：

1. 这是计算能力还是语义能力吗？
   如果是，优先放 `analysis`
2. 这是流程中的一个独立步骤吗？
   如果是，放 `pipeline/nodes`
3. 这只是装配顺序和默认参数吗？
   如果是，放 `pipeline`
4. 这是一次性但仍有内部复杂度的逻辑吗？
   放专用 node 或 support，不直接堆进 pipeline

## 对 pipeline 的硬约束

所有 pipeline 应遵守：
- 尽量保持薄
- 不直接承载复杂算法实现
- 不直接承载厚重语义映射逻辑
- 不直接承载大量结果整理逻辑，除非该逻辑确实只属于这一条流程

如果一个 pipeline 文件同时在做：
- 顺序编排
- 计算实现
- 结果分类
- 可视化绘制

说明分层已经失效，需要拆分。

## 对 contact 类流程的应用

以 contact pipeline 为例：
- residue contact 计算进入 `analysis/trajectory`
- residue 语义注释进入 `analysis/topology`
- contact heatmap 绘制进入 `analysis/topology`
- contact_frequency/contact_annotation/contact_heatmap 分别作为 node
- ContactFrequencyPipeline 只负责把这些 node 串起来

## 默认实践

后续新增分析流程时，默认采用下面模式：
- 先抽 analysis 元件
- 再包装 node
- 最后定义 pipeline

除非逻辑极小且明显一次性，否则不要直接把业务实现写进 pipeline。
