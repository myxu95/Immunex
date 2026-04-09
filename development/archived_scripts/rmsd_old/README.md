# Archived RMSD Scripts

This directory contains old RMSD calculation scripts that have been replaced by the unified RMSD interface.

**Archived on**: 2026-01-08

## Archived Scripts

### 1. `calculate_tcr_overall_rmsd.py`
**原功能**: 计算TCR整体RMSD (align to pHLA)

**替代方案**:
```bash
python scripts/calculate_rmsd.py md.tpr md.xtc --align pHLA --calc TCR
```

**Python API**:
```python
from immunex.analysis.trajectory import RMSDInterface
rmsd = RMSDInterface("md.tpr", "md.xtc")
result = rmsd.calculate(align="pHLA", calc="TCR")
```

---

### 2. `batch_cdr3_rmsd.py`
**原功能**: 批量计算CDR3β RMSD (align to TCR)

**替代方案**:
```bash
python scripts/calculate_rmsd.py md.tpr md.xtc \
    --align TCR --calc CDR3_beta \
    --cdr3-beta CASSLGQAYEQYF
```

**批量处理**:
```python
from immunex.analysis.trajectory import RMSDInterface
rmsd = RMSDInterface("md.tpr", "md.xtc")
calculations = [
    {'align': 'TCR', 'calc': 'CDR3_beta',
     'cdr3_sequences': {'beta': 'CASSLGQAYEQYF'}}
]
df = rmsd.batch_calculate(calculations)
```

---

### 3. `generate_tcr_rmsd_index.py`
**原功能**: 生成TCR RMSD的索引文件

**替代方案**:
Index生成现已集成到 `IndexGenerator`，自动处理

```python
from immunex.utils import IndexGenerator
gen = IndexGenerator("md.tpr", auto_standardize=True)
# 自动生成 base_components.ndx（包含TCR和所有固定组分）
```

---

### 4. `regenerate_phla_tcr_index.py`
**原功能**: 重新生成pHLA-TCR索引文件

**替代方案**:
Index生成现已集成到 `IndexGenerator`，自动处理chain标准化

```python
from immunex.utils import IndexGenerator
gen = IndexGenerator("md.tpr", auto_standardize=True)
gen.generate_multi_component_index(['pHLA', 'TCR'])
```

---

### 5. `rmsd_calculator.py`
**原功能**: 通用批量RMSD计算工具

**替代方案**:
新的 `calculate_rmsd.py` 提供更灵活的组分选择

```bash
# 旧方式 (固定group 3)
python scripts/rmsd_calculator.py ./data --fit-group 3 --calc-group 3

# 新方式 (灵活组分选择)
python scripts/calculate_rmsd.py md.tpr md.xtc --align pHLA --calc TCR
```

---

## 新架构优势

### 1. **统一接口**
- 单一入口点：`RMSDInterface` 或 `calculate_rmsd.py`
- 所有组分计算使用相同API

### 2. **混合Index策略**
- 固定组分（HLA, pHLA, peptide, TCR等）使用统一base index
- CDR3组分按需生成，支持不同序列
- 最小化文件数量（1-3个index vs 8+个）

### 3. **自动chain标准化**
- 自动检测并标准化PDB chain顺序
- 基于残基数量排序：C < B < D < E < A
- 支持TPR/GRO自动转换

### 4. **灵活组分选择**
- 支持任意align-calc组合
- 不再局限于固定的group编号
- 语义化命名（pHLA, TCR而非group 20, 21）

## Migration Guide

### 旧脚本 → 新命令对照表

| 旧脚本 | 旧命令 | 新命令 |
|--------|--------|--------|
| `calculate_tcr_overall_rmsd.py` | `python calculate_tcr_overall_rmsd.py` | `python calculate_rmsd.py md.tpr md.xtc --align pHLA --calc TCR` |
| `batch_cdr3_rmsd.py` | `python batch_cdr3_rmsd.py` | `python calculate_rmsd.py md.tpr md.xtc --align TCR --calc CDR3_beta --cdr3-beta SEQ` |
| `rmsd_calculator.py` | `python rmsd_calculator.py ./data -t calpha` | `python calculate_rmsd.py md.tpr md.xtc --align pHLA --calc TCR` |

### 文档参考

- **完整指南**: `docs/RMSD_CALCULATION_GUIDE.md`
- **使用示例**: `examples/rmsd_unified_example.py`
- **API文档**: `immunex/analysis/trajectory/rmsd_interface.py`

---

## 为什么归档？

1. **功能重复**: 多个脚本执行类似任务但接口不统一
2. **维护困难**: 分散的脚本难以维护和更新
3. **Index管理混乱**: 每个脚本独立生成index，导致文件冗余
4. **缺乏灵活性**: 固定的align-calc组合，难以扩展

新的统一接口解决了这些问题，提供更清晰、灵活和高效的RMSD计算方案。

---

**注意**: 这些归档脚本仍然可用，但不再推荐使用。建议迁移到新的统一接口。
