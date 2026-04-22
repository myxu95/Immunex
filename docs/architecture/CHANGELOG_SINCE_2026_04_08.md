# 2026-04-08 之后的软件改进记录

本文档用于汇总 `2026-04-08` 之后 `Immunex` 的软件层改进，覆盖：

- 已提交到 Git 历史的正式修改
- 当前工作区中尚未发布、但已经实现的主要改动

文档目标不是逐文件罗列，而是从产品面、架构面和能力面说明软件在 `2026-04-08` 之后发生了什么变化。

## 1. 追溯范围

- 本仓库最早可稳定追溯的 Git 提交时间：`2025-09-22`
- 本文档关注范围：`2026-04-08` 之后
- 当前工作区状态：存在一批尚未提交的主线改动，因此本文档分为：
  - `已提交改进`
  - `未发布改进（当前工作区）`

## 2. 已提交改进

### 2.1 2026-04-08：建立 Immunex 产品表面

对应提交：
- `f906e9b` `Introduce Immunex product surface and core analysis platform`

本次改进的核心意义是：仓库开始从“脚本集合”转向“产品化分析平台”。

主要变化：
- 建立 `Immunex` 对外产品表面
- 收口根目录入口与主 README
- 明确 `imn` 命令行为主入口
- 强化 `analysis -> node -> pipeline -> cli -> report` 的主线结构
- 让仓库从开发态逐步转向对外发布态

### 2.2 2026-04-08：公开仓库面移除测试目录

对应提交：
- `709d90b` `Remove tests from published repository surface`

主要变化：
- 将 `tests/` 从公开仓库表面移除
- 通过 `.gitignore` 控制测试目录不再误进入发布面
- 将“内部验证”与“外部产品面”区分开

### 2.3 2026-04-09：开发工作区结构清理

对应提交：
- `041c5c4` `Clean development workspace structure`

主要变化：
- 整理 `development/` 下的开发期资源
- 归档历史资产、实验脚本和内部报告
- 清理工作区顶层散落内容
- 让主目录与产品面更干净

### 2.4 2026-04-09：环境定义与实际依赖对齐

对应提交：
- `73684ef` `Align environment files with active dependencies`

主要变化：
- 校正 `environment.yml`
- 校正 `requirements.txt`
- 校正 `setup.py`
- 补齐当前主线真正依赖的包
- 统一 Python 版本要求与实际环境

## 3. 未发布改进（当前工作区）

以下内容已经在当前工作区实现，但截至本文档编写时仍未形成正式 Git 发布节点。

### 3.1 报告系统增强

#### 3.1.1 单体系 HTML 报告结构增强

已实现的方向：
- 报告首页布局优化
- 左侧身份卡与指标卡分层
- 新增 `Cluster` section
- 新增 `RRCS` section
- `Downloads` 区改为目录式分组，而不再是平铺文件列表

#### 3.1.2 section 语义重构

为了减少内容重叠，报告中的若干 section 已重构为更清楚的阅读任务：

- `Occupancy` 重命名并重构为更强调稳定性的页面语义
- `Contact` 重构为全局界面格局层
- `Interactions` 重构为具体残基对解释层

目标是让用户区分：
- 全局界面格局
- 稳定性/持续性
- 具体残基对机制

#### 3.1.3 Cluster 报告接入

已接入的内容包括：
- `Cluster` 菜单入口
- `Cluster` 摘要指标
- 时间演化图
- cluster summary
- cluster feature digest
- 目录化下载项

当前 cluster 报告强调：
- 状态占比
- 代表构象
- 簇特异性 signature
- 界面差异解释

### 3.2 新增 RRCS 模块

当前工作区中已新增 RRCS 主线能力。

核心新增内容：
- 新分析模块：
  - `immunex/analysis/interactions/rrcs.py`
- 新 node：
  - `RRCSNode`
- 新 pipeline：
  - `RRCSPipeline`
- 新 CLI：
  - `imn rrcs`

当前 RRCS 的方法学实现要点：
- 参照 `gmx_RRCS` 项目的定义
- 非氢原子对打分
- `3.23–4.63 Å` 分段线性距离映射
- residue-pair RRCS 为 atom-pair score 求和
- 同链近邻残基自动去除 backbone 邻接影响

当前 RRCS 已经支持：
- 默认 pair scope
- 收紧 pair scope
- 自定义 pair 文件
- batch 批量执行
- compare 模块接入
- 单体系 HTML 报告接入

### 3.3 新增 compare 模块

当前工作区中已实现一个新的双条件对比框架。

主要内容：
- 新分析目录：
  - `immunex/analysis/comparison/`
- 新 node：
  - `SystemComparisonNode`
- 新 CLI：
  - `imn compare systems`
- 新报告构建脚本：
  - `scripts/build_comparison_report_html.py`

设计原则：
- 不让 compare 模块重新跑单体系分析
- 先读取两个已有 case root
- 再构造对比摘要、对比表和对比 HTML

当前 compare 已支持：
- 中性语义对比
- `comparison_mode`
- `comparison_context`
- `Identity / Quality / Interface / Flexibility / Interaction` 五层基础对比
- `RRCS` 对比接入

### 3.4 RRCS 报告呈现优化

为避免 RRCS 结果难以理解，当前工作区中对 RRCS 的展示做了几项改造：

- 报告中不再突出“候选 pair 总数”
- 优先突出“识别到的有信号配对数”
- 抽象的 RRCS 图改成更容易解释的热点界面区域和热点残基对
- 标签从原始链位点编码改为更可读的残基对表达

目标是让用户更直接地回答：
- 哪些区域最活跃
- 哪些残基对最强
- 信号覆盖是否稳定

### 3.5 pipeline/nodes 目录结构重构

当前工作区中，`immunex/pipeline/nodes/` 已经完成从顶层平铺到分组子包的渐进式重构。

当前结构：
- `nodes/trajectory/`
- `nodes/topology/`
- `nodes/interactions/`
- `nodes/interface/`
- `nodes/allostery/`
- `nodes/comparison/`
- `nodes/geometry/`

顶层 `nodes/__init__.py` 的角色已经变成：
- 统一导出层
- 兼容层

这次重构的意义：
- 降低顶层节点目录的平铺复杂度
- 让 node 按职责分组
- 为后续 workflow 化、可视化编排和 node 注册机制打基础

### 3.6 冗余脚本与旧实现清理

当前工作区已经执行了大规模瘦身，核心目标是：
- 移除明显与当前主线重复的旧实现
- 移除废弃脚本、遗留兼容层和游离工具

已删除或下沉的一类内容包括：
- 旧的 comparative analysis 实现
- 旧 trajectory 综合分析器
- 旧 RMSD/距离/RDF/接触数脚本
- `immunex/legacy`
- 一批 SLURM 辅助生成器
- 一批旧 scripts 和独立工具脚本

这次瘦身的原则是：
- 保留主线
- 删除高置信度冗余
- 不把实验性和历史实现继续放在活动代码树里

### 3.7 `examples/` 与 `docs/` 持续清理

当前工作区已经继续做了两类清理：

#### `examples/`
- 将旧示例下沉到归档位置
- 增加活动示例入口说明
- 降低过时示例对新用户的误导

#### `docs/`
- 去掉活动文档中明显过时的旧命令叙述
- 迁移部分迁移期文档到 archive
- 保持活动文档面更贴近当前主线

## 4. 当前改进的整体方向

从 `2026-04-08` 之后的软件演化来看，当前改进可以概括成 4 条主线：

### 4.1 从脚本集合走向平台

表现为：
- 产品入口收口
- 环境定义收口
- report 主线明确
- compare / RRCS / cluster 等能力以平台形式接入

### 4.2 从单模块分析走向解释型结果系统

表现为：
- 报告页结构增强
- cluster / RRCS / compare 加入解释层
- 不再只是给原始数值和表格，而是开始给“如何读”的层次

### 4.3 从平铺实现走向分层和分组

表现为：
- `analysis` 子模块持续扩充
- `pipeline/nodes` 按领域分组
- 旧遗留脚本逐步移除

### 4.4 为未来的自定义 workflow / 可视化编排做准备

`nodes` 分组化和 pipeline 主线收口的一个重要价值是：
- 为后续的 node manifest
- registry
- validator
- workflow builder
- 图形化拖拽编排

打下基础。

## 5. 当前状态说明

截至本文档编写时：

- `2026-04-08` 之后的正式 Git 提交已存在
- 但工作区中还有一批重大改进尚未发布
- 因此当前软件状态应视为：
  - `已发布产品面 + 未发布主线增强`

如果后续要正式发版，建议：

1. 将未发布改进分批提交
2. 建立正式 `CHANGELOG.md`
3. 对外区分：
   - `released`
   - `unreleased`

## 6. 建议的后续动作

为了让 `2026-04-08` 之后的改进真正完成收口，建议下一步做：

1. 提交当前工作区主线改动
2. 建立正式版本日志文件
3. 对 compare / RRCS / report / nodes 重构补充公开文档
4. 决定哪些未发布功能进入下一次正式发布面

