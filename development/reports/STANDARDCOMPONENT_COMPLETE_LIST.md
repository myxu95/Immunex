# IndexManager StandardComponent 完整列表

**更新日期**: 2026-03-17
**版本**: v2.0 (添加CDR loops和Interface组分)

---

## 组分总览

IndexManager现在支持 **15个标准组分**，分为两大类：

### 📊 统计

- **静态组分**: 8个 (基于链ID，始终可用)
- **动态组分**: 7个 (需要额外信息)
  - CDR loops: 6个 (需要序列信息)
  - Interface: 1个 (需要距离计算)

---

## 1. 静态组分 (Static Components)

这些组分基于标准化的链ID，在base index中始终可用。

| 组分名称 | Enum | 链ID | 说明 |
|---------|------|------|------|
| **HLA** | `StandardComponent.HLA` | A, B | MHC复合物（重链+β2m） |
| **HLA_alpha** | `StandardComponent.HLA_ALPHA` | A | MHC重链（肽结合沟） |
| **HLA_beta** | `StandardComponent.HLA_BETA` | B | β2-微球蛋白 |
| **pHLA** | `StandardComponent.PHLA` | A, B, C | pMHC复合物（HLA+肽） |
| **peptide** | `StandardComponent.PEPTIDE` | C | 抗原肽 |
| **TCR** | `StandardComponent.TCR` | D, E | T细胞受体（α+β链） |
| **TCR_alpha** | `StandardComponent.TCR_ALPHA` | D | TCR α链 |
| **TCR_beta** | `StandardComponent.TCR_BETA` | E | TCR β链 |

### 使用示例

```python
from immunex.core import PipelineContext

context = PipelineContext(
    system_id="1ao7",
    topology="md.tpr",
    trajectory_raw="md_pbc.xtc"
)

# 获取静态组分的组号
phla_id = context.index_manager.get_group_id('pHLA')    # "20"
tcr_id = context.index_manager.get_group_id('TCR')      # "21"
```

---

## 2. 动态组分 - CDR Loops

这些组分需要提供CDR序列信息，通过序列匹配生成索引。

| 组分名称 | Enum | 链ID | 说明 |
|---------|------|------|------|
| **CDR1_alpha** | `StandardComponent.CDR1_ALPHA` | D | TCR α链CDR1环 |
| **CDR2_alpha** | `StandardComponent.CDR2_ALPHA` | D | TCR α链CDR2环 |
| **CDR3_alpha** | `StandardComponent.CDR3_ALPHA` | D | TCR α链CDR3环 |
| **CDR1_beta** | `StandardComponent.CDR1_BETA` | E | TCR β链CDR1环 |
| **CDR2_beta** | `StandardComponent.CDR2_BETA` | E | TCR β链CDR2环 |
| **CDR3_beta** | `StandardComponent.CDR3_BETA` | E | TCR β链CDR3环 |

### 使用方法

#### 单个CDR3

```python
# 生成CDR3 beta索引
cdr3_id = context.index_manager.ensure_cdr3_index(
    chain='beta',
    sequence='CASSLGQAYEQYF',
    ca_only=True
)
print(f"CDR3 beta group ID: {cdr3_id}")  # "25"
```

#### 全部CDR loops

```python
# 生成所有CDR loops索引
group_ids = context.index_manager.ensure_cdr_loops_index(
    chain='beta',
    cdr1_seq='SGHNS',
    cdr2_seq='SYNSPL',
    cdr3_seq='CASSLGQAYEQYF',
    ca_only=True
)

# 输出:
# {
#     'CDR1_beta': '25',
#     'CDR2_beta': '26',
#     'CDR3_beta': '27'
# }
```

---

## 3. 动态组分 - Interface

TCR-pMHC接触界面残基，基于距离截断（默认4.5 Å）生成。

| 组分名称 | Enum | 链ID | 说明 |
|---------|------|------|------|
| **INTERFACE** | `StandardComponent.INTERFACE` | A,B,C,D,E | TCR-pMHC接触界面 |

### 使用方法

```python
# 生成Interface索引
interface_id = context.index_manager.ensure_interface_index(
    cutoff=4.5,  # 距离截断（Å）
    ca_only=True
)
print(f"Interface group ID: {interface_id}")  # "28"
```

### Interface定义

- **TCR选择**: `chainID D E and protein`
- **pMHC选择**: `chainID A B C and protein`
- **接触标准**: 任意原子间距离 ≤ cutoff
- **包含**: 所有参与接触的残基（来自TCR和pMHC两侧）

---

## 组分属性

每个StandardComponent都有以下属性：

### 基础属性

```python
comp = StandardComponent.CDR3_BETA

# 链ID列表
comp.chains                    # ['E']

# 标准化名称
comp.component_name            # 'CDR3_beta'
```

### 动态组分标记

```python
# 是否为动态组分
comp.is_dynamic                # True (CDR/Interface)

# 是否需要序列信息
comp.requires_sequence         # True (CDR loops only)

# 是否需要距离计算
comp.requires_distance_calculation  # True (Interface only)
```

---

## 常见分析场景

### Scenario 1: TCR相对pMHC的运动

```python
# RMSD: 对齐pMHC，计算TCR
rmsd_calc.calculate_with_index_manager(
    context=context,
    fit_component='pHLA',
    calc_component='TCR'
)
```

### Scenario 2: CDR3环的灵活性

```python
# RMSF: CDR3 beta环
context.index_manager.ensure_cdr3_index('beta', 'CASSLGQAYEQYF')

# 使用CDR3组进行RMSF分析
rmsf_calc.calculate_with_index_manager(
    context=context,
    component='CDR3_beta'
)
```

### Scenario 3: 界面接触分析

```python
# 生成Interface索引
context.index_manager.ensure_interface_index(cutoff=5.0)

# 分析界面区域的RMSF
rmsf_calc.calculate_with_index_manager(
    context=context,
    component='INTERFACE'
)
```

### Scenario 4: 多CDR loops比较

```python
# 生成所有CDR loops
group_ids = context.index_manager.ensure_cdr_loops_index(
    chain='beta',
    cdr1_seq='...',
    cdr2_seq='...',
    cdr3_seq='...'
)

# 对比不同CDR的RMSD
for cdr_name, group_id in group_ids.items():
    if group_id:
        rmsd_calc.calculate(
            group=group_id,
            output=f"{cdr_name}_rmsd.xvg"
        )
```

---

## 完整列表（按类型）

### 按功能分类

**整体复合物**:
- pHLA (MHC+肽)
- TCR (受体)

**子复合物**:
- HLA (MHC无肽)

**单链**:
- HLA_alpha, HLA_beta (MHC组分)
- peptide (抗原肽)
- TCR_alpha, TCR_beta (TCR组分)

**功能区域**:
- CDR1/2/3_alpha (TCR α链互补决定区)
- CDR1/2/3_beta (TCR β链互补决定区)
- INTERFACE (TCR-pMHC接触面)

---

## API参考

### 组号获取

```python
# 获取组号（字符串）
group_id = index_mgr.get_group_id('pHLA')  # "20"

# 获取完整组信息
group_info = index_mgr.get_group_info('TCR')
# GroupInfo(group_id=21, group_name='TCR', atom_count=800, ...)

# 列出所有可用组分
components = index_mgr.list_available_components()
# ['HLA', 'HLA_alpha', ..., 'pHLA', 'TCR', ...]
```

### CDR索引生成

```python
# 单个CDR3
cdr3_id = index_mgr.ensure_cdr3_index('beta', 'CASSLGQAYEQYF')

# 全部CDR loops
group_ids = index_mgr.ensure_cdr_loops_index(
    chain='beta',
    cdr1_seq='...',
    cdr2_seq='...',
    cdr3_seq='...'
)
```

### Interface索引生成

```python
interface_id = index_mgr.ensure_interface_index(cutoff=4.5)
```

### 组合索引

```python
# 获取包含多个组分的索引文件
combined_index = index_mgr.get_combined_index(
    'pHLA', 'TCR', 'CDR3_beta', 'INTERFACE'
)
```

---

## 注意事项

### CDR序列匹配

- CDR序列必须在蛋白结构中**完全匹配**
- 序列匹配基于单字母氨基酸代码
- 如果序列未找到，返回`None`

### Interface距离截断

- 默认cutoff = 4.5 Å（侧链接触）
- 更大cutoff（如6.0 Å）可包含更广泛的接触区域
- CA-only模式仅选择Cα原子

### Enum值格式

- 静态组分: `['A', 'B']` (列表)
- 动态CDR组分: `(['D'], 'CDR1')` (元组) - 避免Enum别名
- `.chains`属性自动处理两种格式

---

## 扩展性

如需添加更多组分，修改`StandardComponent` Enum：

```python
class StandardComponent(Enum):
    # 现有组分...

    # 新增组分
    BINDING_SITE = (['A', 'C'], 'binding_site')  # 肽结合位点
    TCR_V_DOMAIN = (['D', 'E'], 'v_domain')      # TCR V结构域
```

---

**文档版本**: 2.0
**最后更新**: 2026-03-17
**相关模块**: `immunex/analysis/index_manager.py`
