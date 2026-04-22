# pHLA-TCR 关键 Interface 聚类 Pipeline 设计

## 1. 目标

本设计文档定义一条新的分析主线：

- `pHLA-TCR interface state clustering`

它的目标不是做“纯 CDR3 几何聚类”，而是识别：

- `CDR3α / CDR3β ↔ peptide`
- `CDR3α / CDR3β ↔ α1-helix / α2-helix`

这些关键结合界面在轨迹中出现的**主要功能状态**。

因此，本模块要回答的问题是：

1. 轨迹里存在几类主要的关键 interface 状态
2. 各状态出现的帧数和时间占比是多少
3. 每类状态的代表帧是什么
4. 不同状态下：
   - `CDR3 ↔ peptide`
   - `CDR3 ↔ MHC groove`
   的接触模式有什么差异

---

## 2. 方法学定位

### 2.1 这不是“纯构象聚类”

本模块不把目标定义为：

- `CDR3 loop conformational clustering`

因为在当前任务里，真正关心的是：

- 功能相关的结合界面状态

也就是说，我们不是只想知道：

- `CDR3` 长得像不像

而是更想知道：

- `CDR3` 和 `peptide / groove` 是如何配对的
- 哪几类 residue-pair 组合构成主要界面模式

### 2.2 这也不是“所有原子联合 RMSD”

不推荐直接把：

- `CDR3α`
- `CDR3β`
- `peptide`

所有原子简单拼接成一个对象，再用单一联合 `RMSD` 做聚类。  
原因：

1. 会把 `peptide` 自身形变和 `CDR3` 的功能状态混在一起
2. 无法突出关键 residue-pair 的切换
3. 不利于后续解释 cluster 差异

### 2.3 推荐方法学定义

推荐将本模块表述为：

- `interface-aware conformational clustering`

其核心思想是：

- 以界面主链几何为主
- 以关键侧链朝向为辅
- 以 `CDR3 ↔ peptide / groove` 相互作用模式为功能约束

---

## 3. 特征表示设计

最终聚类不依赖单一特征，而依赖联合表示。

### 3.1 几何主特征 `D_geometry`

第一版建议包含：

1. `CDR3α backbone heavy-atom RMSD`
2. `CDR3β backbone heavy-atom RMSD`
3. `peptide backbone heavy-atom RMSD`
4. `CDR3 tip ↔ peptide core` 相对距离
5. `CDR3α/β ↔ groove centroid` 相对距离

说明：

- backbone 仍然是界面构象主干
- peptide 可以进入，但不是以“独立对象形变”主导，而是作为界面几何的一部分

### 3.2 侧链特征 `D_sidechain`

第一版建议仅对**关键界面侧链**建模，不对所有侧链等权建模。

关键残基来源：

1. `CDR3α / CDR3β` 残基
2. 高频 `contact / occupancy` hotspot 残基
3. 关键化学类型：
   - `Tyr/Trp/Phe`
   - `Arg/Lys/Asp/Glu`
   - `Ser/Thr/Asn/Gln/His`

推荐特征：

1. `χ1`
2. 对长侧链或芳香侧链再加 `χ2`

角度特征表示建议：

- 用 `sin(χ)` / `cos(χ)` 代替原始角度

这样可以避免 `-180/180` 跳变。

### 3.3 相互作用特征 `D_interaction`

这是本模块与纯构象聚类的关键差别。

第一版建议只使用**简化的二值/频率指纹**，不要一上来引入过多类型。

推荐包含：

1. `CDR3α residue ↔ peptide residue` contact fingerprint
2. `CDR3β residue ↔ peptide residue` contact fingerprint
3. `CDR3α/β ↔ alpha1_helix` contact fingerprint
4. `CDR3α/β ↔ alpha2_helix` contact fingerprint

第一版可以只用：

- `coarse contact`

后续再扩展到：

- `hbond`
- `saltbridge`
- `hydrophobic`
- `pipi`
- `cationpi`

---

## 4. 联合距离定义

推荐把聚类距离写成：

\[
D_{joint} = w_g D_{geometry} + w_s D_{sidechain} + w_i D_{interaction}
\]

第一版推荐权重：

- `w_g = 0.45`
- `w_s = 0.25`
- `w_i = 0.30`

说明：

1. 几何仍然是主导
2. 侧链用于区分功能上不同的局部状态
3. interaction fingerprint 用于把 cluster 拉向真正的界面功能状态

如果后续发现：

- 几何 cluster 很清楚，但功能差异不明显  
则可以增加 `w_i`。

如果后续发现：

- interaction 噪声太大，把 cluster 切得过碎  
则降低 `w_i`。

---

## 5. 对齐策略

这是实施时最关键的技术点之一。

### 5.1 不推荐

不推荐直接对：

- `CDR3 + peptide`

整体做对齐，因为这会把真实界面相对摆动吃掉。

### 5.2 推荐

第一版推荐先对：

- `TCR framework / V-domain core`
或
- `pMHC groove core`

做统一刚体对齐，然后再提取界面特征。

建议优先方案：

1. 用 `TCR framework` 对齐
2. 在这个坐标系下比较：
   - `CDR3`
   - `peptide`
   - `CDR3 ↔ peptide/groove`

理由：

- 这样更能保留 `CDR3` 及其相对界面关系的真实变化
- 又不会被整个复合物整体平移/旋转干扰

---

## 6. 聚类算法选择

### 6.1 第一版推荐：层次聚类

输入：

- 预先计算好的 `D_joint`

推荐：

- `agglomerative hierarchical clustering`

原因：

1. 适合结构状态问题
2. 不需要一开始就硬指定 `k`
3. dendrogram 可解释性强

### 6.2 截断方式

第一版建议：

- 用距离阈值切树

而不是一开始就固定 cluster 数。

例如：

- `distance_cutoff = 0.25`

后续可再根据：

- silhouette
- Davies-Bouldin
- occupancy balance

来辅助选阈值。

### 6.3 第二阶段备选

后续如果要增强，可考虑：

- `k-medoids`
- `HDBSCAN`

但第一版不建议优先。

---

## 7. analysis 层设计

根据仓库分层规则，核心计算应进入 `analysis`。

建议新增目录：

- `immunex/analysis/conformation/`

建议文件：

### 7.1 `interface_representation.py`

职责：

1. 提取对齐所需原子组
2. 提取：
   - backbone 几何特征
   - sidechain 特征
   - interaction fingerprint
3. 生成帧级特征矩阵

建议输出：

- `geometry_feature_frame`
- `sidechain_feature_frame`
- `interaction_feature_frame`

### 7.2 `interface_distance.py`

职责：

1. 计算各类距离矩阵：
   - `D_geometry`
   - `D_sidechain`
   - `D_interaction`
2. 组合成：
   - `D_joint`

### 7.3 `interface_clustering.py`

职责：

1. 对 `D_joint` 执行层次聚类
2. 生成：
   - `cluster_labels`
   - `cluster_sizes`
   - `representative_frames`
   - `cluster_population`

### 7.4 `cluster_characterization.py`

职责：

1. 对每个 cluster 做解释性汇总：
   - top contacts
   - peptide contact signature
   - alpha1/alpha2 contact bias
   - representative residue-pairs

---

## 8. node 层设计

建议新增：

### 8.1 `interface_clustering_node.py`

输入：

- `topology`
- `trajectory`
- `structure_pdb`
- `chain_mapping`
- `cdr_detection`

职责：

1. 调用 `analysis/conformation` 相关元件
2. 落盘 cluster 结果
3. 写入 `context.results["interface_clustering"]`

建议输出路径：

- `analysis/conformation/interface_clustering/`

### 8.2 `interface_cluster_characterization_node.py`

职责：

1. 从 cluster labels 出发做每个 cluster 的功能注释
2. 输出 cluster 级别摘要表与图

如果第一版想简化，也可以先合并到一个 node 中，后续再拆。

---

## 9. pipeline 设计

建议新增：

- `InterfaceClusteringPipeline`

### 9.1 推荐节点顺序

1. `ChainIdentificationNode`
2. `CDRDetectionNode`
3. `InterfaceClusteringNode`
4. `InterfaceClusterCharacterizationNode`

### 9.2 Pipeline 目标

输出一套单体系结果，至少包含：

1. 帧级 cluster 指派
2. cluster 占比
3. representative frame
4. cluster 差异总结

### 9.3 暂不建议

第一版不要直接融入 HTML 总报告，先把功能跑通。

---

## 10. 输出设计

建议输出目录：

- `analysis/conformation/interface_clustering/`

### 10.1 表格

#### `frame_cluster_assignments.csv`

字段建议：

- `frame`
- `time_ps`
- `cluster_id`

#### `cluster_summary.csv`

字段建议：

- `cluster_id`
- `n_frames`
- `fraction`
- `representative_frame`
- `medoid_frame`
- `mean_distance_to_medoid`

#### `cluster_interface_signature.csv`

字段建议：

- `cluster_id`
- `top_peptide_contacts`
- `top_alpha1_contacts`
- `top_alpha2_contacts`
- `dominant_tcr_chain`
- `dominant_cdr_region`

#### `distance_components_summary.csv`

字段建议：

- `frame_i`
- `frame_j`
- `D_geometry`
- `D_sidechain`
- `D_interaction`
- `D_joint`

这个文件不一定全量保存；可以只保存摘要或压缩格式。

### 10.2 矩阵文件

建议保存：

- `joint_distance_matrix.npz`
- `geometry_distance_matrix.npz`
- `interaction_distance_matrix.npz`

### 10.3 图

#### `cluster_dendrogram.png`

展示层次聚类树。

#### `cluster_population.png`

展示各 cluster 占比。

#### `cluster_transition_timeline.png`

展示随时间的 cluster 切换。

#### `cluster_interface_signature.png`

展示不同 cluster 的 peptide/groove 接触差异。

---

## 11. 与现有模块的关系

### 11.1 与 `contact / occupancy`

本模块直接复用：

- residue-pair 语义
- `CDR3 ↔ peptide / α1 / α2` 接触信息

但它的目标不是替代 `occupancy`，而是：

- 在 cluster 层面解释不同功能状态

### 11.2 与 `RMSF`

cluster 可用于解释：

- 某段 `CDR3` 为什么 RMSF 高
- 是因为在两个清晰 cluster 间切换
- 还是因为单一 cluster 内部柔性大

### 11.3 与 `NMA`

NMA 提供：

- 候选关键残基

本模块提供：

- 轨迹里真实出现的主要界面状态

二者互补。

---

## 12. 第一版实施范围建议

为了保证可交付，建议第一版范围控制为：

1. 聚焦：
   - `CDR3α / CDR3β / peptide`
   - 以及 `CDR3 ↔ α1 / α2`
2. 特征：
   - backbone 几何
   - 关键侧链朝向
   - coarse contact fingerprint
3. 聚类：
   - 层次聚类
4. 输出：
   - `frame_cluster_assignments.csv`
   - `cluster_summary.csv`
   - `cluster_population.png`
   - `cluster_dendrogram.png`
   - `cluster_interface_signature.csv`

先不要第一版就做：

- typed interaction 全量加权
- HTML 交互展示
- batch clustering
- 多状态自动比较

---

## 13. 当前推荐实施结论

推荐将该功能正式定义为：

- `pHLA-TCR interface state clustering`

并按下面分层实现：

1. `analysis/conformation`
   - 特征提取
   - 距离构建
   - 聚类
   - cluster 注释
2. `pipeline/nodes`
   - 步骤封装与落盘
3. `pipeline`
   - `InterfaceClusteringPipeline`

一句话总结：

> 这条 pipeline 的核心不是“纯 CDR3 几何聚类”，而是“以几何为主、以侧链和相互作用为辅的关键结合界面状态聚类”。

