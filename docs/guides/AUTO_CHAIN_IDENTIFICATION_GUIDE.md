# 自动链识别功能指南

## 概述

自动链识别功能允许用户无需手动指定链信息即可分析TCR-pMHC对接角度。系统会自动：

1. 检测文件格式（PDB vs TPR/GRO）
2. 提取蛋白质链序列
3. 使用ANARCI + 长度启发式识别链类型
4. 生成MDAnalysis选择字符串
5. 计算对接角度

## 主要特性

### ✅ 零配置分析

```python
from immunex.analysis.angles import DockingAnglePrimaryAnalyzer

# 一行代码完成初始化和链识别
analyzer = DockingAnglePrimaryAnalyzer('md.tpr', 'md.xtc')

# 无需手动指定链！
crossing, incident = analyzer.calculate_docking_angles()
```

### ✅ 支持多种文件格式

- **PDB文件**：使用chainID（A, B, C, D, E）
- **TPR/GRO文件**：使用segname（PROA, PROB, PROC, PROD, PROE）

### ✅ 宽松模式验证

系统采用宽松验证模式，仅需：
- **必需**：HLA-α链
- **必需**：至少1个TCR链（α或β）
- **可选**：peptide和β2m（缺失时仅警告）

### ✅ 向后兼容

旧的手动模式仍然完全支持：

```python
# 旧API仍然可用
analyzer = DockingAnglePrimaryAnalyzer('md.tpr', auto_identify_chains=False)
crossing, incident = analyzer.calculate_docking_angles(
    mhc_selection='segname PROA',
    tcr_alpha_selection='segname PROD',
    tcr_beta_selection='segname PROE'
)
```

---

## 技术架构

### 模块组成

```
immunex/utils/
├── tpr_chain_extractor.py          # TPR/GRO链序列提取
├── chain_identification_adapter.py  # 统一PDB/TPR识别接口
└── selection_string_builder.py      # MDAnalysis选择字符串生成
```

### 识别策略

系统使用与 `IntelligentChainIdentifier` 相同的识别策略：

1. **Peptide识别**（确定性，置信度=1.0）
   - 长度 ≤ 20个氨基酸

2. **Beta2-microglobulin识别**（确定性，置信度=1.0）
   - 长度在90-110个氨基酸范围内

3. **TCR链识别**（ANARCI，置信度=0.95）
   - 使用ANARCI识别TCR-α和TCR-β
   - 回退方案：长度启发式（置信度=0.6）

4. **HLA-α识别**（排除法，置信度=0.7-0.9）
   - 剩余的链（去除peptide、β2m、TCR后）
   - 长度验证：230-310 AA

### 数据流

```
用户输入 (topology, trajectory)
    ↓
文件格式检测 (.pdb / .tpr / .gro)
    ↓
ChainIdentificationAdapter.identify_chains()
    ├─> .pdb → IntelligentChainIdentifier (ANARCI + 长度)
    └─> .tpr/.gro → TPRChainSequenceExtractor
        ├─> MDAnalysis提取序列 (按segname)
        ├─> ANARCI识别TCR-α/β
        └─> 长度启发式识别peptide/β2m/HLA-α
    ↓
宽松验证 (validate_identification)
    ├─> 必需: HLA-α ✓
    ├─> 必需: TCR-α or TCR-β (至少1个) ✓
    └─> 可选: peptide, β2m (缺失仅警告)
    ↓
SelectionStringBuilder.build_all_selections()
    └─> 返回: {
            'mhc_selection': 'segname PROA',
            'tcr_alpha_selection': 'segname PROD',
            'tcr_beta_selection': 'segname PROE'
        }
    ↓
DockingAnglePrimaryAnalyzer.calculate_docking_angles()
    └─> 返回: (crossing_angle, incident_angle)
```

---

## API参考

### DockingAnglePrimaryAnalyzer

#### 初始化

```python
DockingAnglePrimaryAnalyzer(
    topology: str,
    trajectory: Optional[str] = None,
    auto_identify_chains: bool = True,
    use_anarci: bool = True
)
```

**参数**：
- `topology`：拓扑文件路径（PDB、TPR、GRO）
- `trajectory`：轨迹文件路径（XTC、TRR），可选
- `auto_identify_chains`：是否自动识别链（默认True）
- `use_anarci`：是否使用ANARCI识别TCR（默认True）

#### 计算对接角度

```python
analyzer.calculate_docking_angles(
    mhc_selection: Optional[str] = None,
    tcr_alpha_selection: Optional[str] = None,
    tcr_beta_selection: Optional[str] = None
) -> Tuple[float, float]
```

**参数**：
- 全部可选（自动模式下）
- 如果提供，则优先使用手动参数

**返回**：
- `crossing_angle`：Crossing角度（度）
- `incident_angle`：Incident角度（度）

#### 轨迹分析

```python
analyzer.calculate_docking_angles_trajectory(
    mhc_selection: Optional[str] = None,
    tcr_alpha_selection: Optional[str] = None,
    tcr_beta_selection: Optional[str] = None,
    stride: int = 1,
    output_file: Optional[str] = None
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]
```

**参数**：
- 选择参数全部可选（自动模式）
- `stride`：帧步长
- `output_file`：CSV输出文件路径

**返回**：
- `times`：时间点（皮秒）
- `crossing_angles`：Crossing角度数组
- `incident_angles`：Incident角度数组

---

## 使用示例

### 示例1：完全自动模式

```python
from immunex.analysis.angles import DockingAnglePrimaryAnalyzer

# PDB文件
analyzer_pdb = DockingAnglePrimaryAnalyzer('1ao7.pdb')
crossing, incident = analyzer_pdb.calculate_docking_angles()
print(f"Crossing: {crossing:.2f}°, Incident: {incident:.2f}°")

# TPR轨迹
analyzer_tpr = DockingAnglePrimaryAnalyzer('md.tpr', 'md.xtc')
times, crossing_traj, incident_traj = analyzer_tpr.calculate_docking_angles_trajectory(
    stride=10,
    output_file='docking_angles.csv'
)
```

### 示例2：检查识别结果

```python
analyzer = DockingAnglePrimaryAnalyzer('md.tpr', 'md.xtc')

# 打印识别摘要
print(analyzer.chain_adapter.get_chain_summary(analyzer.identifications))

# 输出示例:
# Chain Identification Summary:
# ==================================================
#   PROA  : HLA_alpha    | 275 AA | conf=0.90
#   PROB  : beta2m       |  99 AA | conf=1.00
#   PROC  : peptide      |   9 AA | conf=1.00
#   PROD  : TCR_alpha    | 205 AA | conf=0.95 (ANARCI)
#   PROE  : TCR_beta     | 244 AA | conf=0.95 (ANARCI)
# ==================================================

# 访问自动生成的选择字符串
for key, value in analyzer.auto_selections.items():
    print(f"{key}: {value}")
```

### 示例3：混合模式（部分自动 + 部分手动）

```python
# 自动识别所有链
analyzer = DockingAnglePrimaryAnalyzer('md.tpr', 'md.xtc')

# 手动覆盖MHC选择，其他使用自动识别
crossing, incident = analyzer.calculate_docking_angles(
    mhc_selection='segname PROA and resid 50:180'  # 自定义范围
    # tcr_alpha_selection 和 tcr_beta_selection 使用自动值
)
```

### 示例4：批量处理

```python
structures = ['1ao7.pdb', '1bd2.pdb', '1oga.pdb']

results = []
for pdb in structures:
    try:
        analyzer = DockingAnglePrimaryAnalyzer(pdb)
        crossing, incident = analyzer.calculate_docking_angles()
        results.append((pdb, crossing, incident))
        print(f"✓ {pdb}: {crossing:.2f}° / {incident:.2f}°")
    except Exception as e:
        print(f"✗ {pdb}: {e}")
```

---

## 故障排查

### 问题1：链识别失败

**错误信息**：
```
ValueError: Chain identification failed validation. Errors: ERROR: Missing HLA-alpha chain
```

**解决方案**：
1. 检查结构文件是否包含所有必需的链
2. 使用手动模式指定链：
   ```python
   analyzer = DockingAnglePrimaryAnalyzer(topology, auto_identify_chains=False)
   ```

### 问题2：ANARCI识别失败

**警告信息**：
```
WARNING: ANARCI not available: anarci command not found
WARNING: Will use fallback length-based identification
```

**解决方案**：
1. 安装ANARCI：
   ```bash
   pip install anarci
   ```
2. 或禁用ANARCI：
   ```python
   analyzer = DockingAnglePrimaryAnalyzer(topology, use_anarci=False)
   ```

### 问题3：TPR文件无chainID信息

这是正常现象。TPR/GRO文件使用segname而非chainID。系统会自动处理这种差异。

### 问题4：缺失peptide或β2m

**警告信息**：
```
WARNING: Missing peptide chain (optional in graceful mode)
WARNING: Missing beta2-microglobulin chain (optional in graceful mode)
```

这是正常警告。系统仍可继续计算对接角度（仅需HLA-α + 1个TCR链）。

---

## 置信度评分

系统为每个识别的链分配置信度评分（0.0-1.0）：

| 链类型 | 识别方法 | 置信度 |
|--------|----------|--------|
| Peptide | 长度 ≤ 20 AA | 1.0 |
| Beta2m | 长度 90-110 AA | 1.0 |
| TCR-α | ANARCI识别 | 0.95 |
| TCR-β | ANARCI识别 | 0.95 |
| TCR-α | 长度回退（无ANARCI） | 0.6 |
| TCR-β | 长度回退（无ANARCI） | 0.6 |
| HLA-α | 排除法（典型长度250-290） | 0.9 |
| HLA-α | 排除法（扩展长度230-310） | 0.7 |
| HLA-α | 排除法（异常长度） | 0.3 |

---

## 性能考虑

### 计算开销

- **PDB识别**：~0.5秒（ANARCI）
- **TPR识别**：~1秒（序列提取 + ANARCI）
- **无ANARCI**：<0.1秒（仅长度启发式）

### 内存占用

- 轻量级：识别过程仅加载拓扑，不加载轨迹
- 适用于大规模批量处理

---

## 相关文档

- [`docking_angles.md`](docking_angles.md) - 对接角度分析完整文档
- [`INTELLIGENT_CHAIN_STANDARDIZATION_GUIDE.md`](INTELLIGENT_CHAIN_STANDARDIZATION_GUIDE.md) - 链标准化指南
- [`examples/auto_chain_identification_usage.py`](../examples/auto_chain_identification_usage.py) - 完整使用示例

---

## 更新日志

### v2026.03.10

- ✅ 实现TPR/GRO链序列提取 (`tpr_chain_extractor.py`)
- ✅ 统一PDB/TPR识别接口 (`chain_identification_adapter.py`)
- ✅ MDAnalysis选择字符串生成 (`selection_string_builder.py`)
- ✅ 集成到 `DockingAnglePrimaryAnalyzer`
- ✅ 宽松验证模式（允许peptide/β2m缺失）
- ✅ 向后兼容旧API

---

**维护者**：Immunex Development Team
**最后更新**：2026-03-10
