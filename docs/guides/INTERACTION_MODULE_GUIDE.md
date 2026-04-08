# Interaction Module Guide

## 模块定位

Immunex 的 `Interaction Module` 面向 `TCR-pMHC MD trajectory`，用于输出三层结果：

1. 粗粒度接触
2. 类型化相互作用
3. 区域级相互作用视图

当前已经可直接使用的 interaction family：

- `contact`
- `hbond`
- `saltbridge`
- `hydrophobic`
- `pipi`
- `cationpi`

---

## 当前主线

### 1. Contact

主流程：

1. 链识别
2. CDR 检测
3. `residue-residue minimum heavy atom distance`
4. pair-level 语义注释
5. 区域级汇总
6. 热图生成

主入口：

- `imn batch contact <base_dir>`

输出根目录示例：

- `analysis/contacts/`

核心结果：

- `residue_contact_frequencies.csv`
- `contact_report.csv`
- `peptide_tcr_contacts.csv`
- `hla_tcr_contacts.csv`
- `groove_tcr_contacts.csv`

---

### 2. Hydrogen Bond

主流程：

1. 链识别
2. CDR 检测
3. 使用 `MDAnalysis HydrogenBondAnalysis` 计算跨界面氢键
4. 聚合成 residue-pair hbond 表
5. pair-level 语义注释
6. 区域级汇总
7. 热图生成

主入口：

- `imn batch hbond <base_dir>`

输出根目录示例：

- `analysis/interactions/hydrogen_bonds/`

核心结果：

- `residue_pair_hbonds.csv`
- `hbond_report.csv`
- `groove_tcr_hbonds.csv`

---

### 3. Salt Bridge

主流程：

1. 链识别
2. CDR 检测
3. 在跨界面 charged atoms 间做距离筛选
4. 聚合成 residue-pair salt bridge 表
5. pair-level 语义注释
6. 区域级汇总
7. 热图生成

主入口：

- `imn batch saltbridge <base_dir>`

输出根目录示例：

- `analysis/interactions/salt_bridges/`

核心结果：

- `residue_pair_salt_bridges.csv`
- `salt_bridge_report.csv`
- `groove_tcr_salt_bridges.csv`

---

### 4. Hydrophobic Contact

主流程：

1. 链识别
2. CDR 检测
3. 在跨界面 sidechain hydrophobic heavy atoms 间做距离筛选
4. 聚合成 residue-pair hydrophobic 表
5. pair-level 语义注释
6. 区域级汇总
7. 热图生成

主入口：

- `imn batch hydrophobic <base_dir>`

输出根目录示例：

- `analysis/interactions/hydrophobic_contacts/`

核心结果：

- `residue_pair_hydrophobic_contacts.csv`
- `hydrophobic_report.csv`
- `groove_tcr_hydrophobic_contacts.csv`

---

### 5. Pi-Pi Interaction

主流程：

1. 链识别
2. CDR 检测
3. 识别跨界面芳香环 residue-pair
4. 基于环心距离与环面夹角筛选 pi-pi
5. pair-level 语义注释
6. 区域级汇总
7. 热图生成

主入口：

- `imn batch pipi <base_dir>`

输出根目录示例：

- `analysis/interactions/pi_interactions/`

核心结果：

- `residue_pair_pi_pi.csv`
- `pi_pi_report.csv`
- `groove_tcr_pi_pi.csv`

---

### 6. Cation-Pi Interaction

主流程：

1. 链识别
2. CDR 检测
3. 识别跨界面芳香环与阳离子 residue-pair
4. 基于环心距离与环法向夹角筛选 cation-pi
5. pair-level 语义注释
6. 区域级汇总
7. 热图生成

主入口：

- `imn batch cationpi <base_dir>`

输出根目录示例：

- `analysis/interactions/cation_pi_interactions/`

核心结果：

- `residue_pair_cation_pi.csv`
- `cation_pi_report.csv`
- `groove_tcr_cation_pi.csv`

---

## 输入要求

所有 typed interaction 批量分析目前都要求任务目录中至少可发现：

- `topology`
- `trajectory`
- `structure`

常见对应文件：

- `md.tpr`
- `processed.xtc`
- `*_converted.pdb`

说明：

- `topology` 主要用于几何/化学判定
- `structure` 主要用于回填链语义和区域语义

---

## 输出结构

### 根目录

每个 interaction family 根目录会输出：

- 原始 residue-pair 表
- 注释后的 report 表
- `peptide_tcr / hla_tcr / groove_tcr` 视图
- `TCRa / TCRb` 热图

### 区域目录

区域目录固定为：

- `cdr1/`
- `cdr2/`
- `cdr3/`
- `non_cdr/`

每个区域目录继续分层：

- `tables/`
- `summaries/`
- `heatmaps/`

其中：

- `tables/` 存筛选后的 pair-level 子表
- `summaries/` 存 residue summary 和 `summary.json`
- `heatmaps/` 存 `TCRa / TCRb` 热图

---

## 统一语义

当前交互模块统一使用以下用户可见语义：

- `peptide`
- `hla`
- `groove`
- `beta2m`
- `TCRa`
- `TCRb`
- `cdr1`
- `cdr2`
- `cdr3`
- `non_cdr`

说明：

- `non_cdr` 是固定术语，不再使用 `framework`
- `groove` 来自 `mhc_region`

---

## CLI 示例

### Contact

```bash
imn batch contact /data/processed_md \
  --workers 4 \
  --stride 100 \
  -o ./output/contact_batch
```

### Hydrogen Bond

```bash
imn batch hbond /data/processed_md \
  --workers 4 \
  --stride 100 \
  --distance-cutoff 3.5 \
  --angle-cutoff 150 \
  -o ./output/hbond_batch
```

### Salt Bridge

```bash
imn batch saltbridge /data/processed_md \
  --workers 4 \
  --stride 100 \
  --distance-cutoff 4.0 \
  -o ./output/saltbridge_batch
```

### Hydrophobic Contact

```bash
imn batch hydrophobic /data/processed_md \
  --workers 4 \
  --stride 100 \
  --distance-cutoff 4.5 \
  -o ./output/hydrophobic_batch
```

### Pi-Pi Interaction

```bash
imn batch pipi /data/processed_md \
  --workers 4 \
  --stride 100 \
  --distance-cutoff 6.5 \
  -o ./output/pipi_batch
```

### Cation-Pi Interaction

```bash
imn batch cationpi /data/processed_md \
  --workers 4 \
  --stride 100 \
  --distance-cutoff 6.0 \
  --normal-angle-cutoff 60 \
  -o ./output/cationpi_batch
```

---

## 适合展示的最小结果集

对每种 interaction family，建议优先展示：

1. `*_report.csv`
2. `groove_tcr_*`
3. `cdr3/tables/peptide_tcr_contacts.csv`
4. `cdr3/heatmaps/`
5. 根目录的 `TCRa / TCRb` 热图

这样可以同时覆盖：

- 全局相互作用
- groove 视图
- CDR3 视图
- TCRa / TCRb 区分

---

## 当前边界

当前已经真正打通：

- `contact`
- `hbond`
- `saltbridge`
- `hydrophobic`
- `pipi`
- `cationpi`

当前仍未实现：

- 更严格的 `pi-pi` 构型分类
- 更严格的 `cation-pi` 电荷态判定

如果继续扩展，应继续复用：

- 统一 pair schema
- 通用 residue-pair annotation
- region summary builder
- 统一 heatmap 生成器
