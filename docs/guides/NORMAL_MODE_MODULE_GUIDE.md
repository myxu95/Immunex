# Normal Mode 模块技术说明

本文档不是概念简介，而是当前 `Immunex` 中 `normal mode / PRS` 模块的**技术规格说明**。目标是把下面这些问题讲清楚：

- 当前实现到底依赖什么库
- 输入结构如何被转换成弹性网络
- `GNM / ANM / PRS / hinge_score` 在代码里是怎么计算的
- 每个输出字段和图到底是什么意思
- 当前版本与教科书实现相比做了哪些简化
- 结果应该怎样读，哪些地方不能过度解读

对应代码入口：

- 算法实现：
  - [normal_mode.py](/home/xumy/work/development/Immunex/immunex/analysis/allostery/normal_mode.py)
- 节点封装：
  - [normal_mode_node.py](/home/xumy/work/development/Immunex/immunex/pipeline/nodes/normal_mode_node.py)
- pipeline 组装：
  - [analysis_pipelines.py](/home/xumy/work/development/Immunex/immunex/pipeline/analysis_pipelines.py)
- CLI：
  - [nma.py](/home/xumy/work/development/Immunex/immunex/cli/commands/nma.py)
- 最小测试：
  - [test_normal_mode_pipeline.py](/home/xumy/work/development/Immunex/tests/test_normal_mode_pipeline.py)

---

## 1. 模块定位

`normal mode` 模块的定位是：

- 用单个代表结构构建 `Cα` 弹性网络模型
- 计算低频集体运动相关指标
- 输出 residue-level 与 region-level 的关键位点候选
- 为后续与 `Contact / Occupancy / RMSF / BSA` 联合解释提供结构力学线索

它**不是**：

- 逐帧轨迹动力学分析器
- 真实时间尺度上的运动模拟器
- 亲和力变化或自由能差异的直接计算器

换句话说，它回答的是：

- 从当前构象出发，哪些残基更像集体运动中的关键点、铰链点、扰动传播点

而不是：

- 真实 MD 里某个位点一定怎么运动

---

## 2. 当前技术依赖

当前实现没有依赖 `ProDy`，而是完全使用现有 Python 科学计算栈：

- `MDAnalysis`
  - 读取结构
  - 选择 `Cα` 原子
- `numpy`
  - 向量、矩阵与数值运算
- `scipy.spatial.distance.cdist`
  - 构造 residue 间欧氏距离矩阵
- `scipy.linalg.eigh`
  - 对对称实矩阵做特征分解
- `pandas`
  - residue-level 与 region-level 表格构建
- `matplotlib`
  - 输出三张基础图

这意味着当前实现是：

- 轻依赖
- 易部署
- 可读性高

但也意味着：

- 没有直接复用专业 NMA 工具包的高级特性
- 某些指标是最小可用版实现，而不是完整学术软件级实现

---

## 3. 输入、输出与调用链

### 3.1 输入

算法层最核心的输入只有一个：

- `structure_pdb`

但为了生成带语义的结果表，还需要：

- `chain_mapping`
- `cdr_detection`

因此实际调用链是：

1. `ChainIdentificationNode`
2. `CDRDetectionNode`
3. `NormalModeNode`

### 3.2 CLI

当前正式命令：

```bash
imn nma --structure <代表结构PDB> -o <输出目录>
```

示例：

```bash
imn nma \
  --structure output/parallel2_preprocess/1OGA_sd_run2/md_processed_converted.pdb \
  -o output/nma_demo_1OGA
```

参数：

- `--cutoff`
  - ENM 连边 cutoff，单位 Å
  - 默认 `10.0`
- `--low-modes`
  - 用于 `hinge_score` 的低频模式数
  - 默认 `10`
- `--prs-forces`
  - PRS 随机扰动方向数
  - 默认 `8`

### 3.3 输出

默认输出目录：

- `analysis/allostery/`

包含：

- `normal_mode_residue_scores.csv`
- `normal_mode_region_summary.csv`
- `normal_mode_summary.json`
- `mode_mobility_profile.png`
- `hinge_profile.png`
- `prs_ranking.png`

---

## 4. 结构节点是如何选择的

### 4.1 当前节点定义

当前每个 residue 只保留一个节点：

- `Cα`

代码位置：

- [`_select_ca_atoms()`](/home/xumy/work/development/Immunex/immunex/analysis/allostery/normal_mode.py)

逻辑：

1. 从 `chain_mapping` 中读取所有有效链
2. 构造选择串：
   - `chainID A or chainID B ...`
3. 从输入结构里选择：
   - `(<chain_selection>) and name CA`

如果没有任何 `Cα`，直接报错。

### 4.2 当前意味着什么

这一步带来的技术含义是：

- 所有 residue 被等质量、等弹簧常数地简化为一个 `Cα` 节点
- side chain 细节没有进入 ENM
- peptide 与蛋白主体在模型里没有区别对待
- 所有链共享同一套网络构建规则

这是经典 coarse-grained ENM 的最常见简化，也是第一版最合理的实现策略。

---

## 5. 弹性网络是如何构建的

### 5.1 距离矩阵

取出所有 `Cα` 坐标后，先计算 residue 间两两欧氏距离：

\[
D_{ij} = \|r_i - r_j\|
\]

实现：

- `pairwise_distance = cdist(positions, positions)`

### 5.2 接触掩码

当前连边规则是：

- `distance > 0`
- 且 `distance <= cutoff_angstrom`

对应布尔矩阵：

\[
C_{ij} =
\begin{cases}
1, & 0 < D_{ij} \le r_c \\
0, & \text{otherwise}
\end{cases}
\]

实现：

```python
contact_mask = (pairwise_distance > 0.0) & (pairwise_distance <= cutoff_angstrom)
```

当前默认：

- `r_c = 10.0 Å`

### 5.3 弹簧常数

当前所有边都使用统一弹簧常数：

- `gamma = 1.0`

也就是说，当前没有：

- 距离衰减
- residue type 依赖
- 链类型依赖
- 界面特殊加权

这是一种**均匀弹簧网络**。

---

## 6. GNM 的技术实现

### 6.1 Kirchhoff 矩阵

首先构造邻接矩阵：

\[
A_{ij} = C_{ij}
\]

然后构造 Kirchhoff 矩阵：

\[
\Gamma_{ij} =
\begin{cases}
-\gamma A_{ij}, & i \ne j \\
\sum_{k \ne i} \gamma A_{ik}, & i = j
\end{cases}
\]

实现位置：

- [`_build_kirchhoff()`](/home/xumy/work/development/Immunex/immunex/analysis/allostery/normal_mode.py)

代码逻辑：

```python
adjacency = contact_mask.astype(float)
kirchhoff = -gamma * adjacency
np.fill_diagonal(kirchhoff, -kirchhoff.sum(axis=1))
```

### 6.2 特征分解

用 `scipy.linalg.eigh` 对 `kirchhoff` 做特征分解。

因为：

- `eigh` 返回升序特征值
- 零模对应整体平移

所以当前代码保留：

- `eigenvalue > 1e-8`

作为有效正模。

实现：

- [`_positive_modes()`](/home/xumy/work/development/Immunex/immunex/analysis/allostery/normal_mode.py)

阈值：

- `threshold = 1e-8`

### 6.3 协方差矩阵

当前 GNM 协方差矩阵使用所有正模构造：

\[
\mathbf{C}_{GNM} = \sum_{m \in positive} \frac{1}{\lambda_m} u_m u_m^T
\]

实现：

```python
inv_lambda = 1.0 / eigenvalues[positive]
covariance = (eigenvectors[:, positive] * inv_lambda) @ eigenvectors[:, positive].T
```

### 6.4 `gnm_mobility`

当前 `gnm_mobility` 定义为：

\[
M_i^{GNM} = C_{ii}
\]

也就是协方差矩阵对角线。

实现：

```python
mobility = np.diag(covariance)
```

### 6.5 `gnm_modes`

虽然协方差矩阵使用所有正模，但用于保留输出的低频模式只取：

- `positive[:n_low_modes]`

注意：

- `eigh` 输出升序，所以这里拿到的是最小正特征值对应的低频模式

---

## 7. ANM 的技术实现

### 7.1 Hessian 矩阵

ANM 对每一对接触残基 `(i, j)` 构造一个 `3x3` 子块。

定义：

\[
\Delta r_{ij} = r_j - r_i
\]

\[
\hat{u}_{ij} = \frac{\Delta r_{ij}}{\|\Delta r_{ij}\|}
\]

\[
H_{ij}^{(block)} = \gamma \hat{u}_{ij}\hat{u}_{ij}^T
\]

然后：

- 对角块加 `+outer`
- 非对角块加 `-outer`

实现位置：

- [`_build_hessian()`](/home/xumy/work/development/Immunex/immunex/analysis/allostery/normal_mode.py)

当前实现是标准的均匀 ANM block 组装方式。

### 7.2 特征分解

同样使用：

- `eigh(hessian)`

并保留：

- `eigenvalue > 1e-8`

作为有效模式。

注意：

- ANM 理论上有 6 个刚体零模
- 当前代码没有显式写“去 6 个零模”，而是统一用 `> 1e-8` 过滤
- 这是可行的，但对不同数值条件更敏感

### 7.3 协方差矩阵

当前 ANM 协方差也使用所有正模构造：

\[
\mathbf{C}_{ANM} = \sum_{m \in positive} \frac{1}{\lambda_m} v_m v_m^T
\]

### 7.4 `anm_mobility`

对每个 residue，取对应 `3x3` 协方差块的迹：

\[
M_i^{ANM} = \mathrm{tr}(C_{ii}^{3\times3})
\]

实现：

```python
mobility = np.array(
    [
        float(np.trace(covariance[3*i:3*i+3, 3*i:3*i+3]))
        for i in range(n_residues)
    ]
)
```

这表示该 residue 在三维方向上的总方差贡献。

---

## 8. `hinge_score` 的技术实现

### 8.1 当前定义

当前 `hinge_score` 不是教科书里某个固定标准公式，而是一个**启发式铰链近似指标**。

实现步骤：

1. 取 ANM 低频模式：
   - `selected = positive[:min(n_low_modes, positive.size)]`
2. 对每个模式，把特征向量 reshape 成：
   - `(n_residues, 3)`
3. 计算每个 residue 在该 mode 下的向量振幅平方：

\[
a_{i,m} = \|\mathbf{v}_{i,m}\|^2
\]

4. 按特征值倒数加权累加：

\[
A_i = \sum_m \frac{a_{i,m}}{\lambda_m}
\]

5. 做“反向归一化”：

\[
hinge_i = normalize(\max(A) - A_i)
\]

也就是说：

- 振幅越小
- 越像铰链点
- 得分越高

### 8.2 当前局限

这一指标在当前版本里有两个特点：

1. 可用于排序
2. 绝对值区分度偏弱

因此：

- 看 `hinge_score` 的相对前后顺序是有意义的
- 直接解释 `0.99` 和 `0.97` 的绝对差别没有太大意义

---

## 9. PRS 的技术实现

### 9.1 PRS 的思路

PRS 要回答的是：

- 扰动 residue `i` 时，其它 residue 会怎样响应

当前实现完全基于 `ANM covariance`。

### 9.2 扰动方向

当前不是对每个 residue 施加一个固定方向扰动，而是：

1. 用高斯分布生成随机方向
2. 再归一化成单位向量

实现：

```python
rng = np.random.default_rng(20260406)
force_directions = rng.normal(size=(n_force_directions, 3))
force_directions /= np.linalg.norm(force_directions, axis=1, keepdims=True)
```

关键点：

- 随机种子固定为 `20260406`
- 因此同样输入下，结果可复现

### 9.3 单个位点的响应

对于被扰动 residue `j`：

1. 取协方差矩阵中对应的三列块：

\[
B_j = C[:, 3j:3j+3]
\]

2. 对于每个随机单位力方向 `f_k`：

\[
d_k = B_j f_k
\]

3. 将 `d_k` reshape 成 `(n_residues, 3)`，并取每个 residue 位移模长：

\[
r_{i,k} = \|d_{i,k}\|
\]

4. 对所有方向平均，得到这一轮 perturbed residue 对所有 residue 的平均响应：

\[
R_{j,i} = \frac{1}{K}\sum_k r_{i,k}
\]

于是得到一个：

- `response_matrix`
- 维度 `(n_residues, n_residues)`

### 9.4 `prs_effectiveness`

当前定义：

\[
effectiveness_j = mean_i(R_{j,i})
\]

即：

- 第 `j` 个 residue 被扰动时，对整体平均造成多大影响

### 9.5 `prs_sensitivity`

当前定义：

\[
sensitivity_i = mean_j(R_{j,i})
\]

即：

- 第 `i` 个 residue 对所有外部扰动平均有多敏感

### 9.6 当前注意事项

当前 `PRS` 是：

- 基于线性响应近似
- 基于单结构 ANM 协方差
- 基于固定随机种子和有限方向采样

所以它更适合：

- 相对排名
- 候选位点筛选

而不适合解释：

- 精确量化的真实扰动传播幅度

---

## 10. 归一化与综合打分

### 10.1 `_normalize_01`

当前归一化函数：

\[
normalize(x_i)=\frac{x_i-\min(x)}{\max(x)-\min(x)}
\]

如果：

- `max - min < 1e-12`

则返回全零数组。

### 10.2 当前参与归一化的列

在 `_build_residue_frame()` 中，当前会计算：

- `gnm_mobility_norm`
- `anm_mobility_norm`
- `prs_effectiveness_norm`
- `prs_sensitivity_norm`
- `network_degree_norm`

### 10.3 `combined_key_residue_score`

当前综合分数公式是：

\[
score_i =
0.32 \cdot hinge_i +
0.30 \cdot prsEff_i^{norm} +
0.13 \cdot prsSens_i^{norm} +
0.25 \cdot degree_i^{norm}
\]

也就是代码中的：

```python
combined_score = (
    0.32 * hinge_score[idx]
    + 0.30 * prs_eff_norm[idx]
    + 0.13 * prs_sens_norm[idx]
    + 0.25 * network_degree_norm[idx]
)
```

### 10.4 为什么没有把 `gnm/anm_mobility` 直接放进综合分

因为当前版本更想找：

- 可能的关键点
- 可能的铰链点
- 可能的传播点

而不是仅仅找“最能动”的点。

如果把 `mobility` 直接高权重放进去，很容易把高柔性边缘环区排得过高。

---

## 11. 语义注释是如何接入的

算法本身只产生匿名 residue 指标，但 `Immunex` 会在输出阶段自动补语义。

调用：

- `ComplexResidueSemanticAnnotator`

对每个 residue 补这些字段：

- `component`
  - `HLA_alpha`
  - `beta2m`
  - `peptide`
  - `TCR_alpha`
  - `TCR_beta`
- `complex_side`
  - `phla`
  - `tcr`
- `phla_region`
- `mhc_region`
- `mhc_subregion`
  - `alpha1_helix`
  - `alpha2_helix`
  - `non_groove`
  - `peptide`
- `tcr_chain`
  - `alpha`
  - `beta`
- `tcr_region`
  - `CDR1 / CDR2 / CDR3 / non_cdr`
- `tcr_region_detailed`
  - `CDR3_alpha` 等

### 11.1 `region_group` 的派生规则

当前规则：

1. 如果存在 `tcr_region_detailed`，优先使用它
2. 否则如果存在 `mhc_subregion`，使用它
3. 否则如果 `component != "unknown"`，使用 component
4. 否则返回 `unknown`

这个字段是 region-level summary 的主分组键。

---

## 12. 输出文件的逐项解释

### 12.1 `normal_mode_residue_scores.csv`

这是 residue-level 主表。

关键字段解释：

- `chain_id / resid / resname / residue_label`
  - 残基身份
- `component / complex_side`
  - 属于 `TCR` 还是 `pHLA`
- `tcr_region / tcr_region_detailed`
  - 如果是 TCR 残基，属于哪个 CDR
- `mhc_subregion`
  - 如果是 MHC 重链残基，属于 `alpha1_helix / alpha2_helix / non_groove`
- `gnm_mobility / anm_mobility`
  - 未归一化原值
- `*_norm`
  - 0 到 1 归一化值
- `hinge_score`
  - 铰链近似指标
- `prs_effectiveness`
  - 扰动传播输出能力
- `prs_sensitivity`
  - 对外部扰动的敏感性
- `network_degree`
  - ENM 邻接度
- `combined_key_residue_score`
  - 第一轮筛选总分

### 12.2 `normal_mode_region_summary.csv`

这是按 `region_group` 聚合后的区域表。

字段：

- `region_group`
- `n_residues`
- `mean_gnm_mobility`
- `mean_anm_mobility`
- `mean_hinge_score`
- `mean_prs_effectiveness`
- `max_combined_score`

推荐重点看：

- `mean_prs_effectiveness`
- `max_combined_score`

### 12.3 `normal_mode_summary.json`

主要包含：

- `n_residues`
- `network_edges`
- `cutoff_angstrom`
- `n_low_modes`
- `prs_force_directions`
- `mean_*`
- `prs_matrix_shape`
- `top_key_residues`
- `top_hinges`
- `top_regions`

这份文件适合：

- HTML 概览
- 后续批量汇总
- 快速 sanity check

### 12.4 图文件

#### `mode_mobility_profile.png`

- 横轴：按 `chain_id, resid` 排序后的 residue index
- 纵轴：归一化 mobility
- 两条线：
  - `GNM mobility`
  - `ANM mobility`

#### `hinge_profile.png`

- 横轴：residue index
- 纵轴：`hinge_score`
- 填充图展示整体 hinge-like 分布

#### `prs_ranking.png`

- 取 `combined_key_residue_score` 前 `top_n=15`
- 画水平条形图
- 标签格式：
  - `RESIDUE (region_group)`

---

## 13. 结果应该怎样读

推荐顺序：

### 第一步：看区域

先看：

- `normal_mode_region_summary.csv`

判断：

- 是 `CDR3_beta` 更突出
- 还是 `alpha1_helix / alpha2_helix`
- 还是 `peptide`

### 第二步：看残基

再看：

- `normal_mode_residue_scores.csv`

先按：

- `combined_key_residue_score`

降序，再聚焦这些区域：

- `CDR3_alpha`
- `CDR3_beta`
- `alpha1_helix`
- `alpha2_helix`
- `peptide`

### 第三步：做交叉验证

当前最推荐的交叉对象：

- `Contact`
  - 这个 residue 是否真的参与界面
- `Occupancy`
  - 这个 residue 相关 interaction 是否稳定
- `RMSF`
  - 它是稳定控制点还是高柔性噪声位点
- `BSA`
  - 界面是否稳定埋藏

### 一个实用判断模板

如果某个 residue 同时满足：

- `combined_key_residue_score` 高
- `prs_effectiveness` 高
- 位于 `CDR3 / peptide / alpha1 / alpha2`
- 同时也出现在高 occupancy interaction 里

那么它就比“只在 NMA 里高分”的 residue 更值得优先关注。

---

## 14. 当前实现与标准 NMA 软件相比的差异

当前版本是**最小可用实现**，和完整学术软件相比有这些简化：

1. 只使用 `Cα`
2. 统一 `gamma = 1.0`
3. 连边只按单一 cutoff 判断
4. 没有质量加权、B 因子加权或 residue type 权重
5. 没有显式状态比较
6. 没有模式动画
7. `hinge_score` 是启发式定义，不是标准化铰链检测算法
8. `PRS` 基于固定随机种子和有限方向采样
9. 没有把 `contact/occupancy` 先验直接融合进网络

这些简化不会让模块失效，但决定了它更适合：

- 候选筛选
- 结构机制提示

不适合：

- 最终机制定论

---

## 15. 当前已知局限

### 15.1 `hinge_score` 绝对值解释意义弱

当前 `hinge_score` 更适合看排名，不适合把绝对值当成强结论。

### 15.2 `combined_key_residue_score` 是人为加权分

它的价值在于第一轮筛选，而不是严格物理量。

### 15.3 输入结构质量直接影响结果

如果代表结构本身偏离主要构象，NMA 结果也会偏。

### 15.4 没有轨迹平均

当前不是：

- 对 ensemble 平均结构做 NMA

而是：

- 对单个代表结构做 NMA

---

## 16. 当前最推荐的使用方式

推荐流程：

1. 跑 `imn nma`
2. 看 `region_summary`
3. 看 `top residues`
4. 和：
   - `contact`
   - `occupancy`
   - `RMSF`
   - `BSA`
   交叉验证
5. 最后形成候选位点清单

一句话：

> `normal mode` 负责“提出候选”，现有动力学与界面模块负责“验证候选”。

---

## 17. 复现实例

当前真实回归示例：

```bash
imn nma \
  --structure output/parallel2_preprocess/1OGA_sd_run2/md_processed_converted.pdb \
  -o output/nma_demo_1OGA
```

默认参数：

- `cutoff = 10.0`
- `low_modes = 10`
- `prs_forces = 8`

如果要改参数：

```bash
imn nma \
  --structure output/parallel2_preprocess/1OGA_sd_run2/md_processed_converted.pdb \
  -o output/nma_demo_1OGA_alt \
  --cutoff 12.0 \
  --low-modes 12 \
  --prs-forces 16
```

---

## 18. 后续增强方向

如果当前结果表现稳定，最有价值的增强方向是：

1. 用 `contact / occupancy` 做 interface-aware 打分修正
2. 增加 `batch nma`
3. 支持多状态结构比较
4. 引入更稳的 hinge 检测方式
5. 再考虑进入 HTML 报告

当前阶段的重点仍然是：

- residue ranking 是否合理
- region ranking 是否有解释力
- 与现有模块交叉后是否能缩小候选空间

