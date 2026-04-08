# Contact / Interaction Module Design

## 目标

Immunex 的 `Contact / Interaction Module` 应当服务于一条明确主线：

- 面向 `TCR-pMHC MD trajectory`
- 先得到稳定的 `pair-level dynamic results`
- 再进行 `typed interaction` 判定
- 最后做 `region-level interpretation`

本设计文档用于固定：

- 模块分层
- 统一结果契约
- 与 `analysis / node / pipeline` 的映射关系
- 后续扩展路线

---

## 核心原则

### 原则 1：分析层级和工程层级分离

下面两套分层不能混淆：

- 分析层级
  - 粗粒度接触
  - 类型化相互作用
  - 区域级相互作用
- 工程层级
  - `analysis`
  - `node`
  - `pipeline`

`粗粒度接触 / 类型化相互作用 / 区域级相互作用` 是结果语义分层。  
`analysis / node / pipeline` 是实现分层。

### 原则 2：区域级分析不能重新发明计算器

以下内容都必须基于同一个 `pair-level annotated table` 做过滤和汇总：

- `peptide ↔ CDR3β`
- `peptide ↔ TCR total`
- `pMHC groove ↔ TCRα/β`
- `CDR loop ↔ peptide / MHC`

区域级分析只能：

- filter
- aggregate
- visualize

不能重新写一套单独几何主循环。

### 原则 3：typed interaction 是 pair universe 的子层

`contact` 是最宽松的邻近关系。  
`hydrogen bond / salt bridge / hydrophobic / π–π / cation–π` 是更强约束的 interaction family。

因此系统必须允许：

- 同一对 residue 既有 `contact`
- 又有 `hbond`
- 又有 `salt_bridge`

这些 family 不互斥。

### 原则 4：先统一 schema，再扩 interaction family

在 `typed interaction` 扩展前，必须先固定统一的 `pair schema`。  
否则每种 interaction 都会长成一套独立格式，后续不可维护。

---

## 模块分层

### Layer A：粗粒度接触层

职责：

- 计算 `residue-residue minimum heavy atom distance`
- 计算 `contact occupancy`
- 生成 `contact frequency matrix`

特点：

- 只做几何邻近判定
- 不解释 biological semantics
- 不关心 `CDR / peptide / HLA`

建议归属：

- `immunex/analysis/trajectory/residue_contacts.py`

当前状态：

- 已实现并已接入主 pipeline

### Layer B：类型化相互作用层

职责：

- 在 pair-level 邻近关系上识别 interaction family

第一批 family：

- `hydrogen_bond`
- `salt_bridge`
- `hydrophobic_contact`
- `pi_pi`
- `cation_pi`

特点：

- 依赖几何阈值、原子类型、电荷、芳香环中心等规则
- 输出必须仍然是 pair-level table

建议归属：

- `immunex/analysis/interactions/`

建议文件：

- `hydrogen_bonds.py`
- `salt_bridges.py`
- `hydrophobic_contacts.py`
- `pi_interactions.py`

当前状态：

- `hydrogen_bond` 已接入统一主链
- `salt_bridge` 已接入统一主链
- `hydrophobic_contact` 已接入统一主链
- `pi_pi` 已接入统一主链
- `cation_pi` 已接入统一主链

### Layer C：语义注释层

职责：

- 为 pair-level 结果增加结构语义

至少包括：

- `component_1 / component_2`
- `complex_side_1 / complex_side_2`
- `tcr_chain`
- `tcr_region`
- `tcr_region_detailed`
- `phla_region`
- `interaction_class`

未来必须补充：

- `groove_region`

说明：

如果后续要支持：

- `pMHC groove ↔ TCRα/β`

则 `groove` 必须进入 residue semantics 体系，而不能只在汇总阶段临时推导。

建议归属：

- `immunex/analysis/topology/residue_pair_annotation.py`

当前状态：

- 当前由以下模块共同承担：
  - `complex_residue_semantics.py`
  - `tcr_residue_semantics.py`
  - `contact_annotation.py`

建议后续把 `contact_annotation.py` 升级成更通用的 `residue_pair_annotation`

### Layer D：区域级汇总层

职责：

- 从 annotated pair table 中按规则切专题结果
- 生成汇总表、排名、候选位点和热图

典型专题：

- `peptide ↔ CDR3β`
- `peptide ↔ TCR total`
- `pMHC groove ↔ TCRα/β`
- `CDR loop ↔ peptide / MHC`

建议归属：

- `immunex/analysis/topology/region_interaction_summary.py`

当前状态：

- contact 模块已有部分 region summary 能力
- 但尚未抽象成统一 summary owner

---

## 与 Node / Pipeline 的关系

### analysis 层负责“做什么”

建议最终形成以下能力模块：

- `residue_contacts.py`
- `hydrogen_bonds.py`
- `salt_bridges.py`
- `hydrophobic_contacts.py`
- `pi_interactions.py`
- `residue_pair_annotation.py`
- `region_interaction_summary.py`
- `contact_heatmap.py`

### node 层负责“流程步骤包装”

建议节点分工：

- `ContactFrequencyNode`
  - 调粗粒度接触层
- `InteractionTypingNode`
  - 调类型化 interaction 层
- `ResiduePairAnnotationNode`
  - 调语义注释层
- `RegionInteractionSummaryNode`
  - 调区域级汇总层
- `InteractionHeatmapNode`
  - 调绘图模块

### pipeline 层负责“怎么串”

推荐的 contact / interaction 主流程：

1. `ChainIdentificationNode`
2. `CDRDetectionNode`
3. `ContactFrequencyNode`
4. `ResiduePairAnnotationNode`
5. `RegionInteractionSummaryNode`
6. `InteractionHeatmapNode`

针对 typed interaction 的流程可以派生，例如：

1. `ChainIdentificationNode`
2. `CDRDetectionNode`
3. `InteractionTypingNode`
4. `ResiduePairAnnotationNode`
5. `RegionInteractionSummaryNode`

---

## 统一 Pair Schema

### 必选字段

所有 `contact / typed interaction` 输出都必须至少包含：

- `interaction_family`
- `residue_1`
- `residue_2`
- `resid_1`
- `resid_2`
- `resname_1`
- `resname_2`
- `chain_1`
- `chain_2`
- `component_1`
- `component_2`
- `frames_present`
- `total_frames`
- `occupancy`
- `frequency`

### 接触层补充字段

- `min_distance_angstrom`
- `mean_distance_angstrom`
- `distance_metric`

### typed interaction 补充字段

- `geometry_score`
- `donor_acceptor_info`
- `charge_pair_info`
- `ring_center_distance`
- `ring_plane_angle`
- `criteria_passed`

### 语义注释字段

- `complex_side_1`
- `complex_side_2`
- `phla_region_1`
- `phla_region_2`
- `tcr_chain`
- `tcr_region`
- `tcr_region_detailed`
- `interaction_class`

### 保留扩展字段

- `metadata_json`

当不同 interaction family 的专属几何字段差异太大时，允许额外几何信息进入 `metadata_json`，但不允许省略统一主字段。

---

## 输出结构建议

### 根目录原则

不要让每种 interaction family 复制整棵目录树。

推荐结构：

```text
analysis/
  interactions/
    pair_tables/
    region_views/
    heatmaps/
```

### pair_tables

按 interaction family 放主表：

- `contact_pairs.csv`
- `hydrogen_bond_pairs.csv`
- `salt_bridge_pairs.csv`
- `hydrophobic_pairs.csv`
- `pi_pairs.csv`

### region_views

按专题视角放过滤结果：

- `peptide_CDR3b/`
- `peptide_TCR_total/`
- `groove_TCRa/`
- `groove_TCRb/`
- `CDR_loop_peptide/`
- `CDR_loop_MHC/`

### heatmaps

放 pair matrix 风格图：

- `pep_TCRa_heatmap.png`
- `pep_TCRb_heatmap.png`
- `hla_TCRa_heatmap.png`
- `hla_TCRb_heatmap.png`

---

## 当前大问题

### 问题 1：缺统一 pair schema

这是最危险的问题。  
如果现在继续扩 typed interaction，而不先统一 schema，后面输出会彻底碎裂。

### 问题 2：`groove` 语义尚未进入 residue semantics

目前要支持：

- `pMHC groove ↔ TCRα/β`

但 `groove` 还不是稳定语义字段。  
这必须尽快补到 residue semantics 层。

### 问题 3：区域级分析容易重新造轮子

如果不强制规定“区域级分析只做过滤和汇总”，后续会不可避免地长出一批互相重复的小计算器。

### 问题 4：typed interaction 尚未形成统一入口

当前 `hydrogen bond` 有底层能力，但还没有统一 pair schema、统一 annotation、统一 region summary 主链。

---

## 推荐实施顺序

### P0：先收 schema 和语义层

1. 固定统一 pair schema
2. 将 `contact_annotation` 升级为 `residue_pair_annotation`
3. 把 `groove` 正式加入 residue semantics
4. 抽出 `region_interaction_summary`

### P1：再收 typed interaction 主链

1. `hydrogen_bond`
2. `salt_bridge`
3. `hydrophobic`
4. `pi / cation_pi`

每种 family 都必须走：

- pair table
- residue pair annotation
- region summary

### P2：最后做统一汇总与可视化

1. family-level comparison
2. region-level comparison
3. mutation candidate prioritization
4. cross-family heatmap / ranking

---

## 最终结论

Immunex 的 `Contact / Interaction Module` 应当被实现为：

- 一个统一的 `pair-level dynamics` 框架
- 上面叠加 `typed interaction families`
- 再在统一语义体系上做 `region-level interpretation`

因此：

- 三层分析想法是对的
- 它不与 `node -> pipeline` 冲突
- 真正需要先解决的是：
  - 统一 pair schema
  - groove 语义入模
  - 区域级分析只做过滤与汇总
