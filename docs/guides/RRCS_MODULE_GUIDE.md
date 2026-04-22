# RRCS 模块指南

## 模块定位

`RRCS`（Residue-Residue Contact Score）用于量化残基对之间的连续接触强度。它不是二值 contact frequency 的替代，而是一个增强层：

- `contact frequency` 回答“有没有接触、持续多久”
- `RRCS` 回答“接触有多紧、细微重排是否发生”

它尤其适合识别：

- 侧链重排
- 界面轻微收紧/松动
- 二值接触频率相近但打包紧密程度不同的状态

## 方法定义

当前实现参照 `gmx_RRCS` 的公开方法：

1. 只使用非氢原子
2. 对每个原子对距离 `d` 计算分数

\[
s(d)=
\begin{cases}
1, & d \le d_{min} \\
0, & d \ge d_{max} \\
\dfrac{d_{max}-d}{d_{max}-d_{min}}, & d_{min} < d < d_{max}
\end{cases}
\]

默认参数：

- `d_min = 3.23 Å`
- `d_max = 4.63 Å`

单帧 residue-pair RRCS 定义为所有 atom-pair 分数之和：

\[
RRCS_{frame}(i,j)=\sum_{a \in i}\sum_{b \in j} s(d_{ab})
\]

额外规则：

- 对同链且 `|resid_i-resid_j| < 5` 的残基对，自动去掉 backbone 原子（`N/CA/C/O`），避免链内近邻把 RRCS 虚高。

## CLI 用法

### 基本用法

```bash
imn rrcs \
  --structure output/parallel2_preprocess/1OGA_sd_run2/md_processed_converted.pdb \
  --topology input/database/parallel2/1OGA_sd_run2/prod/md.tpr \
  --trajectory output/parallel2_preprocess/1OGA_sd_run2/md_processed.xtc \
  -o output/rrcs_demo_1OGA
```

### 常用参数

- `--stride`：采样步长
- `--radius-min`：最小满分距离，默认 `3.23`
- `--radius-max`：最大零分距离，默认 `4.63`
- `--pair-scope`：默认 pair 范围
- `--pair-file`：自定义 residue pair 文件

### 推荐的 pair scope

- `interface`
- `cdr3_peptide`
- `cdr3_groove`
- `tcr_peptide`
- `tcr_groove`
- `tcr_interface`

如果你只关心关键界面，建议优先从 `cdr3_peptide` 或 `tcr_interface` 开始，避免 pair 数过大。

## 输出文件

输出目录固定为：

- `analysis/interactions/rrcs/`

### `rrcs_timeseries.csv`

逐帧、逐 residue pair 的 RRCS time series。

关键字段：

- `frame`
- `time_ps`
- `chain_id_1`
- `resid_1`
- `residue_label_1`
- `chain_id_2`
- `resid_2`
- `residue_label_2`
- `rrcs`

### `rrcs_pair_summary.csv`

对每个 residue pair 做汇总。

关键字段：

- `mean_rrcs`
- `median_rrcs`
- `max_rrcs`
- `rrcs_nonzero_frames`
- `rrcs_nonzero_fraction`

### `annotated_rrcs_pair_summary.csv`

在 pair summary 基础上补充语义注释，例如：

- `interaction_class`
- `tcr_chain`
- `tcr_region`
- `mhc_subregion`
- `partner_component`

这份文件最适合后续报告和比较模块使用。

### `rrcs_region_summary.csv`

把 RRCS 聚合到区域层，例如：

- `peptide_tcr`
- `hla_tcr`
- `CDR3`
- `non_cdr`
- `alpha1_helix`
- `alpha2_helix`

关键字段：

- `n_pairs`
- `mean_rrcs_sum`
- `mean_rrcs_mean`
- `median_rrcs_mean`
- `max_rrcs_max`
- `rrcs_nonzero_fraction_mean`

### `rrcs_summary.json`

模块摘要，包括：

- 参数
- pair 数量
- 非零 pair 数量
- top residue pairs
- top regions

## 与现有 contact 模块的关系

推荐联合阅读方式：

1. 先看 `contact / occupancy`
   - 哪些 pair 出现
   - 持续性如何
2. 再看 `RRCS`
   - 这些 pair 的接触是否更紧
   - 是否发生 subtle packing shift

典型解读：

- `occupancy` 高、`RRCS` 高：稳定且紧密
- `occupancy` 高、`RRCS` 中等：稳定但相对松散
- `occupancy` 相近、`RRCS` 差异大：界面打包方式改变

## 当前边界

当前版本：

- 已接入 `batch rrcs`
- 已接入 comparison report
- 已接入单体系 interaction HTML 报告

当前尚未做的更进一步能力：

- RRCS typed-interaction 分解
- RRCS 与 NMA/cluster 的自动交叉解释
- RRCS 专用 batch 汇总首页
