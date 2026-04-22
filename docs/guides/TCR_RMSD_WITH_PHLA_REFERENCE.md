# 使用pHLA作为参考坐标系计算TCR RMSD

## 核心概念

**重要说明**: PBC处理流程**不需要修改**，保持使用默认设置。定制化发生在**RMSD计算**步骤。

### 标准工作流 vs 定制工作流

| 步骤 | 标准工作流 | 定制工作流（本指南） |
|-----|-----------|-------------------|
| PBC处理 | 使用Backbone对齐 | **相同**：使用Backbone对齐 |
| RMSD计算 | Fit: Backbone, Calc: Backbone | **不同**：Fit: pHLA, Calc: TCR |
| 测量内容 | 整体结构稳定性 | TCR相对于pHLA的运动 |

### 生物学意义

**标准RMSD**:
- 测量整个蛋白复合物相对于初始结构的偏离
- 适用于评估MD模拟整体质量

**pHLA参考系RMSD**:
- pHLA作为"固定"参考坐标系（MHC-peptide相对稳定）
- TCR RMSD测量受体相对于pHLA的运动程度
- 反映TCR识别过程中的动态变化
- 较低RMSD = 稳定结合，较高RMSD = 动态识别

## 完整工作流程

### 第1步：生成索引文件（包含pHLA和TCR组）

```python
from immunex.analysis.topology import (
    IndexGenerator, IndexGenerationInput,
    IndexGenerationMethod, ComponentDefinition
)

# 定义组分
components = [
    ComponentDefinition(
        name="pHLA",
        selection="chainID A B C"  # MHC alpha + beta2m + peptide
    ),
    ComponentDefinition(
        name="pHLA_backbone",
        selection="chainID A B C and (name CA C N O)"
    ),
    ComponentDefinition(
        name="TCR",
        selection="chainID D E"  # TCR alpha + beta
    ),
    ComponentDefinition(
        name="TCR_backbone",
        selection="chainID D E and (name CA C N O)"
    ),
]

# 生成索引
index_gen = IndexGenerator()
result = index_gen.generate(IndexGenerationInput(
    topology="structure.pdb",
    method=IndexGenerationMethod.CUSTOM,
    standardized_pdb="structure.pdb",
    components=components,
    output_file="phla_tcr_index.ndx"
))
```

### 第2步：PBC处理（标准流程，无需修改）

```python
from immunex.analysis.trajectory.pbc import PBCProcessor

pbc_processor = PBCProcessor()

# 使用标准2-step方法
pbc_processor.remove_pbc_2step(
    trajectory="md.xtc",
    topology="md.tpr",
    output="md_pbc.xtc"
)

# 这将使用默认的Backbone进行对齐 - 完全不需要改动!
```

**重要**: PBC处理保持原样，不需要任何特殊参数。

### 第3步：计算TCR RMSD（使用pHLA作为参考）

**方法A：命令行交互**

```bash
gmx rms -f md_pbc.xtc -s md.tpr -n phla_tcr_index.ndx -o tcr_rmsd_phla_ref.xvg
```

当提示选择组时：
1. **For least squares fit**: 选择 `pHLA_backbone`
2. **For RMSD calculation**: 选择 `TCR_backbone`

**方法B：Python脚本（推荐）**

```python
import subprocess

cmd = [
    "gmx", "rms",
    "-f", "md_pbc.xtc",
    "-s", "md.tpr",
    "-n", "phla_tcr_index.ndx",
    "-o", "tcr_rmsd_phla_ref.xvg"
]

# 通过stdin提供组选择
# 第一行: fit组（pHLA_backbone）
# 第二行: calc组（TCR_backbone）
stdin_input = "pHLA_backbone\nTCR_backbone\n"

result = subprocess.run(
    cmd,
    input=stdin_input,
    text=True,
    capture_output=True,
    check=True
)

print(f"RMSD calculation completed: tcr_rmsd_phla_ref.xvg")
```

**方法C：使用RMSDCalculator**

```python
from immunex.analysis.trajectory.rmsd import RMSDCalculator

rmsd_calc = RMSDCalculator(
    topology="md.tpr",
    trajectory="md_pbc.xtc"
)

# 使用自定义组
rmsd_calc.calculate_gromacs_custom_groups(
    reference_group="pHLA_backbone",  # Fit on pHLA
    analysis_group="TCR_backbone",    # Calculate TCR RMSD
    output_file="tcr_rmsd_phla_ref.xvg"
)
```

### 第4步：结果解读

**输出文件**: `tcr_rmsd_phla_ref.xvg`

```
# Time (ps)   RMSD (nm)
0.000         0.000
10.000        0.234
20.000        0.198
...
```

**统计分析**:

```python
import numpy as np
import pandas as pd

# 读取RMSD数据
df = pd.read_csv('tcr_rmsd_phla_ref.xvg', comment=['#', '@'],
                 delim_whitespace=True, names=['Time', 'RMSD'])

print(f"Mean RMSD: {df['RMSD'].mean():.3f} nm")
print(f"Std RMSD: {df['RMSD'].std():.3f} nm")
print(f"Min RMSD: {df['RMSD'].min():.3f} nm")
print(f"Max RMSD: {df['RMSD'].max():.3f} nm")

# 判断稳定性
if df['RMSD'].mean() < 0.2:
    print("Interpretation: Stable TCR-pMHC binding")
elif df['RMSD'].mean() < 0.4:
    print("Interpretation: Moderate TCR dynamics")
else:
    print("Interpretation: Highly dynamic TCR recognition")
```

## 技术细节

### GROMACS gmx rms 的工作原理

```bash
gmx rms -f trajectory -s reference -n index -o output
```

**两次选择提示**:
1. **Group for least squares fit**: 用于结构叠合的参考组
   - 系统会对齐到这个组上
   - 本方法选择: `pHLA_backbone`

2. **Group for RMSD calculation**: 计算RMSD的目标组
   - 在对齐后计算这个组的RMSD
   - 本方法选择: `TCR_backbone`

**计算流程**:
```
每一帧:
  1. 基于pHLA_backbone进行最小二乘叠合
  2. 在叠合后的坐标中计算TCR_backbone的RMSD
  3. 输出 (时间, RMSD)
```

### 索引文件结构

```
[ pHLA ]
1 2 3 4 ... (6119 atoms)

[ pHLA_backbone ]
1 5 9 13 ... (1533 atoms - CA, C, N, O only)

[ TCR ]
6120 6121 6122 ... (4948 atoms)

[ TCR_backbone ]
6120 6124 6128 ... (1294 atoms - CA, C, N, O only)
```

## 对比分析

### 与标准RMSD的对比

```python
# 标准RMSD计算
subprocess.run([
    "gmx", "rms",
    "-f", "md_pbc.xtc",
    "-s", "md.tpr",
    "-o", "rmsd_standard.xvg"
], input="Backbone\nBackbone\n", text=True, check=True)

# pHLA参考系RMSD计算
subprocess.run([
    "gmx", "rms",
    "-f", "md_pbc.xtc",
    "-s", "md.tpr",
    "-n", "phla_tcr_index.ndx",
    "-o", "tcr_rmsd_phla_ref.xvg"
], input="pHLA_backbone\nTCR_backbone\n", text=True, check=True)

# 绘制对比图
import matplotlib.pyplot as plt

df_std = pd.read_csv('rmsd_standard.xvg', ...)
df_tcr = pd.read_csv('tcr_rmsd_phla_ref.xvg', ...)

plt.plot(df_std['Time'], df_std['RMSD'], label='Standard (Backbone)')
plt.plot(df_tcr['Time'], df_tcr['RMSD'], label='TCR vs pHLA')
plt.xlabel('Time (ps)')
plt.ylabel('RMSD (nm)')
plt.legend()
plt.savefig('rmsd_comparison.png', dpi=300)
```

## 使用场景

### 1. TCR工程化设计

评估突变TCR对不同peptide的结合稳定性：

```python
# 比较野生型 vs 突变体
calculate_tcr_rmsd(wt_trajectory, "wt_tcr_rmsd.xvg")
calculate_tcr_rmsd(mut_trajectory, "mut_tcr_rmsd.xvg")

# 比较稳定性
wt_rmsd = analyze_rmsd("wt_tcr_rmsd.xvg")
mut_rmsd = analyze_rmsd("mut_tcr_rmsd.xvg")

if mut_rmsd < wt_rmsd:
    print("Mutation improves binding stability")
```

### 2. Peptide识别特异性研究

研究同一TCR对不同peptide的识别动态：

```python
peptides = ['NLVPMVATV', 'GILGFVFTL', 'LLWNGPMAV']

for peptide in peptides:
    trajectory = f"md_{peptide}.xtc"
    output = f"tcr_rmsd_{peptide}.xvg"
    calculate_tcr_rmsd(trajectory, output)

    rmsd_mean = analyze_rmsd(output)['mean']
    print(f"{peptide}: TCR RMSD = {rmsd_mean:.3f} nm")
```

### 3. 别构效应分析

研究MHC多态性对TCR结合的影响：

```python
mhc_alleles = ['HLA-A*02:01', 'HLA-A*02:02', 'HLA-A*02:03']

for allele in mhc_alleles:
    trajectory = f"md_{allele}.xtc"
    output = f"tcr_rmsd_{allele}.xvg"
    calculate_tcr_rmsd(trajectory, output)
```

## 测试示例

**脚本**: `development/validation/test_tcr_rmsd_with_phla_reference.py`

运行测试：

```bash
python development/validation/test_tcr_rmsd_with_phla_reference.py
```

预期输出：

```
TCR RMSD Calculation with pHLA as Reference Frame
======================================================================

Step 1: Generating index file with pHLA and TCR groups
----------------------------------------------------------------------
Generated index file: development/test_output/tcr_rmsd_phla_ref/index.ndx
  - pHLA: 6119 atoms
  - pHLA_backbone: 1533 atoms
  - TCR: 4948 atoms
  - TCR_backbone: 1294 atoms

Step 2: PBC processing (standard workflow)
----------------------------------------------------------------------
Use normal PBC processing - no changes needed!

Step 3: Calculate TCR RMSD (using pHLA as reference)
----------------------------------------------------------------------
Key command:
gmx rms -f md_pbc.xtc -s md.tpr -n index.ndx -o tcr_rmsd.xvg

When prompted:
  Select group 1 (pHLA_backbone) for least squares fit
  Select group 2 (TCR_backbone) for RMSD calculation

Key point:
  PBC processing unchanged - customization is in RMSD calculation step!
```

## 常见问题

### Q1: 为什么PBC处理不需要修改？

**A**: PBC处理的目的是消除周期性边界伪影，使用整体Backbone对齐即可。定制化的参考坐标系应该在后续分析步骤（RMSD计算）中实现。

### Q2: 能否在PBC步骤就使用pHLA对齐？

**A**: 技术上可以，但**不推荐**。原因：
- PBC处理应该保持标准化，便于不同分析复用
- 定制化分析在后续步骤更灵活
- 分离关注点：PBC处理专注于物理问题，RMSD计算专注于生物学问题

### Q3: 这种方法和标准RMSD有什么区别？

**A**:
- 标准RMSD: 整体结构相对于初始结构的偏离（质量控制）
- 本方法: TCR相对于pHLA的运动（生物学意义）

### Q4: backbone定义是什么？

**A**: CA, C, N, O 四种主链原子。MDAnalysis选择语法: `(name CA C N O)`

### Q5: 如何选择合适的参考组？

**A**:
- 用于TCR-pMHC研究: pHLA_backbone（本指南）
- 用于抗体-抗原研究: antigen_backbone
- 一般原则: 选择相对稳定的结构作为参考

## 参考文献

- **TCR-pMHC结构生物学**: Garcia KC, et al. (2009) Annu Rev Immunol 27:83-117
- **MD模拟研究TCR识别**: Ayres CM, et al. (2020) Front Immunol 11:2235
- **GROMACS RMSD计算**: GROMACS Manual, gmx rms documentation

## 更新记录

**日期**: 2026-03-18
**修改**: 明确PBC处理不需要修改，定制化在RMSD计算步骤
**测试**: 使用1AO7结构验证索引生成
**状态**: 已验证

---

**技术支持**: Immunex Development Team
**文档版本**: 1.0
