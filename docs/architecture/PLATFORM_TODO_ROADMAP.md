# Immunex 平台收口待办与路线图

**状态**：当前执行清单  
**目标**：把仓库从“重构进行中”收口为“入口清晰、工作流稳定、结果可交付”的分析平台。

## 优先级 P0：入口与产品面统一

### 1. 统一根 README
- 目标：让仓库首页准确反映当前 `imn` CLI、架构分层、核心工作流与文档入口。
- 问题：现有根 README 仍混有旧工作流、旧路径和过时说明。
- 产出：
  - 当前能力总览
  - 核心命令入口
  - 推荐文档索引
  - 单体系与批处理工作流说明

### 2. 固化“默认产品入口”
- 目标：明确用户应优先使用哪些入口。
- 推荐入口：
  - `imn preprocess`
  - `imn quality`
  - `imn rmsd`
  - `imn rmsf`
  - `imn contact`
  - `imn identity`
  - `imn bsa`
  - `imn nma`
  - `imn inter_cluster`
  - `imn batch ...`
  - `imn report interaction`

### 3. 文档索引收口
- 目标：保证根 `README.md` 作为唯一总入口，直接指向稳定文档。
- 动作：
  - 为新增模块补索引
  - 将“现行文档”和“历史文档”继续分离

## 优先级 P1：报告链标准化

### 4. 固化单体系报告契约
- 目标：把 `report interaction` 定义为标准单体系报告入口。
- 需要稳定的 section：
  - Overview
  - Quality
  - Interface
  - Flexibility
  - Occupancy
  - Contact
  - Interactions
  - Cluster
  - Downloads

### 5. 统一 section 裁剪机制
- 目标：支持默认全量、按需裁剪。
- 动作：
  - 保持 `--include-sections`
  - 保持 `--exclude-sections`
  - 明确 section 名称契约

### 6. 报告资源分组规则
- 目标：所有 section 的下载资源、图、表统一收口。
- 重点：
  - `Downloads` 使用目录式分组
  - `Cluster` 使用摘要输出，不暴露冗余中间产物

## 优先级 P2：跨模块解释能力

### 7. 建立 Quality / BSA / RMSF / Occupancy / NMA / Cluster 的交叉解释
- 目标：从“多模块并列展示”升级成“模块间有逻辑关系”。
- 先做：
  - `BSA + interaction richness`
  - `RMSF + persistent contact`
  - `NMA key residues + occupancy hotspots`
  - `Cluster state + interface signature`

### 8. 形成标准摘要规则
- 目标：同一类结果在不同体系下给出一致的自动摘要。
- 例如：
  - dominant interface state
  - largest cluster
  - dominant interaction family
  - key biological identity

## 优先级 P3：批处理与平台化能力

### 9. 批处理汇总页
- 目标：在单体系 HTML 之外，提供 batch index / summary 页面。
- 重点：
  - 每个体系一张卡片
  - 支持从总览页进入单体系报告

### 10. 新模块的 batch 对称化
- 目标：新模块不只支持单体系，还能进入 batch 汇总。
- 当前重点模块：
  - `bsa`
  - `nma`
  - `inter_cluster`

## 优先级 P4：仓库清理与发布准备

### 11. 历史目录与迁移残留清理
- 目标：减少新用户辨识成本。
- 重点区域：
  - 历史 `aftermd/` 残留逻辑
  - 根目录旧文档
  - `development/` 中需要明确标记为实验/归档的内容

### 12. 发布前一致性检查
- 目标：确保 CLI、文档、示例、测试一致。
- 检查项：
  - CLI 帮助文本
  - docs 索引
  - examples 是否仍引用旧路径
  - README 是否与当前命令匹配

## 当前建议执行顺序

1. 根 README 收口
2. docs 索引补齐
3. 报告链 section 与下载契约稳定
4. 模块交叉解释增强
5. batch 汇总页
6. 仓库清理与发布前校验

## 当前已完成的基础

- `analysis -> node -> pipeline -> cli` 主架构已成型
- `contact / typed interactions / occupancy` 已接通
- `identity / bsa / rmsf / nma / inter_cluster` 已进入单体系主链
- `report interaction` 已可生成综合 HTML 报告

## 当前最重要的原则

- 继续保持 `analysis`、`pipeline/nodes`、`pipeline` 的分层
- 新功能先进入标准入口，再考虑 HTML 和 batch
- 默认全量报告，按需裁剪
- 输出优先摘要化，不在最终报告中暴露过多中间产物
