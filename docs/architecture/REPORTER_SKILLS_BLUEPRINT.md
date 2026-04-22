# Reporter Skills Blueprint

## 目标

为 Immunex 后续接入大语言模型分析能力提供四类 skill 分类与固定路径。

这些 skill 负责“解释和建议”，不负责底层分析计算。

## 固定路径

- `skills/reporter-query/SKILL.md`
- `skills/reporter-diagnostic/SKILL.md`
- `skills/reporter-compare/SKILL.md`
- `skills/reporter-design/SKILL.md`

## 分类说明

### 1. Query skills
最轻量。
用于直接回答：
- 哪些残基最关键
- 哪些 hotspot 最强
- 哪个 cluster 占比最高
- 哪个区域 RRCS 最强

### 2. Diagnostic skills
用于单体系诊断。
重点关注：
- 稳定性
- 界面构成
- persistence
- RRCS
- cluster
- NMA mechanical hotspots

### 3. Compare skills
用于双体系比较。
适用场景包括：
- WT vs mutant
- binder vs weak binder
- standard MD vs REST2
- 不同 binding mode

### 4. Design skills
最高层。
用于：
- 候选突变位点推荐
- 候选实验验证位点推荐
- 重点可视化区域推荐

## 推荐实施顺序

1. 先让 `query` 和 `diagnostic` 跑通
2. 再补 `compare`
3. 最后再用 `design` 承接候选推荐

## 后续建议

后面应补：
- `reporter_context` 统一输入契约
- stage 化输出格式
- 证据来源规范
- 高/中/低置信度标签
