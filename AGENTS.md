# AGENTS.md

## 仓库级执行规则

Immunex 是一个分析流程平台，不是脚本集合。

后续在本仓库中实现功能时，默认遵守以下规则：

1. 高复用计算能力、语义能力、结果构建能力优先进入 `immunex/analysis`
2. 流程步骤包装、context 读写、步骤级落盘优先进入 `immunex/pipeline/nodes`
3. `immunex/pipeline` 只负责流程编排、默认参数和节点装配
4. 不要把低复用逻辑直接堆进 pipeline；优先使用专用 node 或 support 模块承接
5. 如果一段逻辑以后可能被第二条流程复用，就不应直接写死在 pipeline 内
6. pipeline 应保持薄；复杂计算、语义判断、可视化构建应尽量下沉

## 新增功能时的归属判断

按下面顺序判断：

1. 这是“如何计算”的逻辑吗？
   放 `analysis`
2. 这是“流程中的一个步骤”吗？
   放 `pipeline/nodes`
3. 这是“先做什么后做什么”的装配关系吗？
   放 `pipeline`
4. 这是低复用但仍然有内部复杂度的逻辑吗？
   放专用 node 或 support，不直接写进 pipeline 主体

## 当前固定约定

- contact 相关流程采用“analysis 元件 + node + pipeline”的分层模式
- TCR 非 CDR 区域统一命名为 `non_cdr`，不再使用 `framework` 作为用户可见术语
- 结果目录优先按用途分层，如 `tables/`、`summaries/`、`heatmaps/`

## 参考文档

完整规则见：
- `docs/architecture/ANALYSIS_PLATFORM_RULES.md`
