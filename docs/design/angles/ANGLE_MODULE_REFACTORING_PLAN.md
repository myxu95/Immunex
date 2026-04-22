# 角度模块结构优化方案

**当前日期**: 2026-03-18

---

## 当前问题分析

### 文件结构混乱

**当前文件** (9个):
```
immunex/analysis/angles/
├── conserved_cysteine_detector.py    21K  # TCR二硫键检测
├── docking_angles_primary.py         19K  # 旧的主分析器
├── docking_angle_analyzer.py         11K  # 新的标准化分析器 (新增)
├── mhc_sequence_aligner.py           12K  # MHC序列比对
├── mhc_groove_detector.py            7.7K # MHC groove几何
├── angle_data_structures.py          7.4K # 数据结构 (新增)
├── plane_fitting.py                  5.2K # SVD平面拟合
├── vector_angles.py                  3.2K # 向量角度工具
└── __init__.py                       3.1K # 导出
```

**总计**: 2556行代码，分散在9个文件

### 主要问题

1. **功能重复**
   - `docking_angles_primary.py` (19K) - 旧API
   - `docking_angle_analyzer.py` (11K) - 新API
   - 两者功能重叠，造成混乱

2. **过度细分**
   - `plane_fitting.py` (5.2K) - 仅一个SVD函数
   - `vector_angles.py` (3.2K) - 几个工具函数
   - `conserved_cysteine_detector.py` (21K) - TCR特定功能
   - 这些都是角度计算的内部实现，不应暴露

3. **入口不清晰**
   - 用户不知道应该用哪个类
   - `DockingAnglePrimaryAnalyzer` vs `DockingAngleAnalyzer`？
   - 是否需要直接使用`MHCGrooveDetector`？

4. **维护成本高**
   - 修改一个功能需要跨多个文件
   - 不利于理解整体逻辑

---

## 优化方案：3文件精简结构

### 目标结构

```
immunex/analysis/angles/
├── analyzer.py            # 唯一主文件 (~1500行)
│   ├── DockingAngleAnalyzer (公开API)
│   └── 内部辅助类（不导出）:
│       ├── _MHCGrooveDetector
│       ├── _TCRAxisCalculator
│       ├── _SequenceAligner
│       └── _GeometryUtils
├── data_structures.py     # 数据结构 (~200行)
│   ├── DockingAngleInput
│   └── DockingAngleResult
└── __init__.py            # 简洁导出 (~20行)
```

**总计**: 3个文件，~1720行（删除30%冗余代码）

---

## 详细设计

### 1. `analyzer.py` - 唯一主文件

**结构**:
```python
"""
Docking Angle Analyzer - Complete implementation

This module contains all logic for TCR-pMHC docking angle calculation.
Internal helper classes are prefixed with underscore and not exported.

Public API:
    - DockingAngleAnalyzer: Main analyzer class

Author: Immunex Development Team
Date: 2026-03-18
"""

import numpy as np
import MDAnalysis as mda
from typing import Tuple, Optional, Dict
import logging

# ============================================================================
# Internal Helper Classes (Not Exported)
# ============================================================================

class _GeometryUtils:
    """Internal geometry utilities (plane fitting, vector angles, etc.)"""

    @staticmethod
    def fit_plane_svd(coords: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """SVD-based plane fitting"""
        centroid = coords.mean(axis=0)
        centered = coords - centroid
        _, _, vh = np.linalg.svd(centered)
        normal = vh[-1]  # Last singular vector
        return centroid, normal

    @staticmethod
    def angle_between_vectors(v1: np.ndarray, v2: np.ndarray) -> float:
        """Calculate angle between two vectors (degrees)"""
        v1_norm = v1 / np.linalg.norm(v1)
        v2_norm = v2 / np.linalg.norm(v2)
        cos_theta = np.clip(np.dot(v1_norm, v2_norm), -1.0, 1.0)
        return float(np.degrees(np.arccos(cos_theta)))


class _SequenceAligner:
    """Internal MHC sequence alignment"""

    def __init__(self, universe: mda.Universe, mhc_selection: str):
        self.universe = universe
        self.mhc_selection = mhc_selection
        # ... implementation ...

    def get_alpha1_alpha2_residues(self) -> Tuple[list, list]:
        """Align to HLA-A*02:01 and get α1/α2 residue IDs"""
        # ... implementation ...


class _MHCGrooveDetector:
    """Internal MHC groove geometry detector"""

    def __init__(self, universe: mda.Universe, mhc_selection: str):
        self.universe = universe
        self.aligner = _SequenceAligner(universe, mhc_selection)

    def detect_groove_axis(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Calculate groove axis (α1 → α2)"""
        # ... implementation ...

    def detect_groove_plane(self) -> Tuple[np.ndarray, np.ndarray]:
        """Fit groove plane using α1/α2 backbone"""
        # ... implementation ...


class _TCRAxisCalculator:
    """Internal TCR axis calculator using conserved cysteines"""

    def __init__(self, universe: mda.Universe):
        self.universe = universe

    def detect_conserved_cysteines(self, chain_selection: str) -> Dict:
        """Detect Cys23 and Cys104 using ANARCI"""
        # ... implementation ...

    def calculate_tcr_axis(
        self,
        tcr_alpha_selection: str,
        tcr_beta_selection: str
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Calculate TCR axis: COM(Vβ disulfide) - COM(Vα disulfide)"""
        # ... implementation ...


# ============================================================================
# Public API
# ============================================================================

class DockingAngleAnalyzer:
    """
    TCR-pMHC Docking Angle Analyzer

    Calculates Crossing and Incident angles for TCR-pMHC complexes.
    Conforms to Immunex 6 design principles.

    Examples
    --------
    Single frame:
    >>> from immunex.analysis.angles import DockingAngleAnalyzer, DockingAngleInput
    >>> analyzer = DockingAngleAnalyzer()
    >>> result = analyzer.analyze(DockingAngleInput(topology='structure.pdb'))

    Trajectory:
    >>> result = analyzer.analyze(DockingAngleInput(
    ...     topology='md.tpr',
    ...     trajectory='md_pbc.xtc',
    ...     stride=10
    ... ))
    """

    def __init__(self):
        self._cancel_flag = False
        self.progress_callback = None

    def analyze(self, input_params: 'DockingAngleInput') -> 'DockingAngleResult':
        """Main analysis method"""
        # ... implementation using internal helper classes ...

    # ... other methods ...
```

**优势**:
- 所有逻辑在一个文件中，易于理解
- 内部类用`_`前缀标记，不导出
- 清晰的public/private分离

---

### 2. `data_structures.py` - 保持不变

数据结构文件保持独立，因为：
1. 清晰分离数据和逻辑
2. 便于类型检查
3. 文档友好

**无需修改**。

---

### 3. `__init__.py` - 极简导出

```python
"""
Angle Analysis Module

Public API:
    - DockingAngleAnalyzer: Main analyzer
    - DockingAngleInput: Input data structure
    - DockingAngleResult: Result data structure

Internal implementation details are not exported.

Example:
    >>> from immunex.analysis.angles import DockingAngleAnalyzer, DockingAngleInput
    >>> analyzer = DockingAngleAnalyzer()
    >>> result = analyzer.analyze(DockingAngleInput(topology='structure.pdb'))
"""

from .analyzer import DockingAngleAnalyzer
from .data_structures import DockingAngleInput, DockingAngleResult

__all__ = [
    'DockingAngleAnalyzer',
    'DockingAngleInput',
    'DockingAngleResult',
]

__version__ = '4.0.0'
__refactoring_date__ = '2026-03-18'
```

**变化**:
- 仅导出3个类
- 移除所有内部类导出
- 清晰的API边界

---

## 迁移计划

### Step 1: 创建新的`analyzer.py`

合并以下文件内容：
```
conserved_cysteine_detector.py  → _TCRAxisCalculator
docking_angles_primary.py       → (核心逻辑) → DockingAngleAnalyzer
docking_angle_analyzer.py       → DockingAngleAnalyzer (保留新API)
mhc_sequence_aligner.py          → _SequenceAligner
mhc_groove_detector.py           → _MHCGrooveDetector
plane_fitting.py                 → _GeometryUtils.fit_plane_svd()
vector_angles.py                 → _GeometryUtils 静态方法
```

### Step 2: 移除旧文件

删除以下文件：
```bash
rm conserved_cysteine_detector.py
rm docking_angles_primary.py
rm mhc_sequence_aligner.py
rm mhc_groove_detector.py
rm plane_fitting.py
rm vector_angles.py
```

保留：
```
analyzer.py              # 新
data_structures.py       # 保留
__init__.py              # 更新
```

### Step 3: 更新导出

**旧导出** (9个类):
```python
__all__ = [
    'PlaneFitter',
    'MHCSequenceAligner',
    'MHCGrooveDetector',
    'DockingAnglePrimaryAnalyzer',  # 删除
    'ConservedCysteineDetector',    # 删除
    'angle_between_vectors',        # 删除
    'DockingAngleAnalyzer',
    'DockingAngleInput',
    'DockingAngleResult',
]
```

**新导出** (3个类):
```python
__all__ = [
    'DockingAngleAnalyzer',
    'DockingAngleInput',
    'DockingAngleResult',
]
```

### Step 4: 向后兼容（可选）

如果需要保持旧API兼容：
```python
# analyzer.py 末尾
# Backward compatibility aliases (deprecated)
DockingAnglePrimaryAnalyzer = DockingAngleAnalyzer  # Deprecated

# __init__.py
from .analyzer import (
    DockingAngleAnalyzer,
    DockingAnglePrimaryAnalyzer,  # Deprecated, use DockingAngleAnalyzer
)
```

---

## 对比

### Before (9文件)
```
immunex/analysis/angles/
├── conserved_cysteine_detector.py    21K
├── docking_angles_primary.py         19K
├── docking_angle_analyzer.py         11K
├── mhc_sequence_aligner.py           12K
├── mhc_groove_detector.py            7.7K
├── angle_data_structures.py          7.4K
├── plane_fitting.py                  5.2K
├── vector_angles.py                  3.2K
└── __init__.py                       3.1K

总计: 2556行, 9个文件
公开API: 9个类/函数
```

### After (3文件)
```
immunex/analysis/angles/
├── analyzer.py                     ~1500行  # 唯一主文件
├── data_structures.py              ~200行   # 数据结构
└── __init__.py                     ~20行    # 简洁导出

总计: ~1720行, 3个文件 (减少30%)
公开API: 3个类
```

---

## 优势总结

### 1. 清晰性 ✅
- 单一主文件，逻辑集中
- 用户明确知道使用`DockingAngleAnalyzer`
- 内部实现细节隐藏

### 2. 可维护性 ✅
- 修改角度计算逻辑只需编辑一个文件
- 内部类重构不影响外部API
- 减少30%冗余代码

### 3. 易用性 ✅
- 简化的导入: `from immunex.analysis.angles import DockingAngleAnalyzer`
- 不需要理解`MHCGrooveDetector`、`PlaneFitter`等内部实现
- 文档更简洁

### 4. 架构一致性 ✅
- 与RMSD、RMSF等模块结构一致
- 每个模块一个主分析器
- 统一的设计模式

### 5. 性能 ✅
- 减少import开销
- 代码局部性更好（CPU缓存友好）

---

## 风险和缓解

### 风险1: 破坏现有代码

**影响**: 中等
**概率**: 低

**受影响代码**:
```python
# 这些导入会失败 ❌
from immunex.analysis.angles import MHCGrooveDetector
from immunex.analysis.angles import ConservedCysteineDetector
from immunex.analysis.angles import angle_between_vectors
```

**缓解**:
1. 在`__init__.py`中添加弃用警告和别名
2. 保留旧API作为向后兼容层（短期）
3. 提供迁移指南

### 风险2: 单文件过大

**影响**: 低
**概率**: 低

**缓解**:
- 1500行对于Python是合理大小
- 清晰的内部类分隔
- IDE友好（折叠、跳转）

---

## 实施时间

**总时间**: 1天

- Step 1 (创建analyzer.py): 3小时
- Step 2 (删除旧文件): 10分钟
- Step 3 (更新导出): 20分钟
- Step 4 (测试): 2小时
- Step 5 (文档更新): 1小时

---

## 测试策略

### 单元测试
```python
def test_analyzer_basic():
    """测试基本功能"""
    from immunex.analysis.angles import DockingAngleAnalyzer
    analyzer = DockingAngleAnalyzer()
    # ...

def test_internal_classes_not_exported():
    """测试内部类不导出"""
    import immunex.analysis.angles as angles
    assert not hasattr(angles, '_MHCGrooveDetector')
    assert not hasattr(angles, '_TCRAxisCalculator')
```

### 集成测试
```python
def test_pipeline_integration():
    """测试Pipeline集成不受影响"""
    from immunex.pipeline import DockingAnglePipeline
    # ...
```

---

## 推荐

**强烈推荐实施此优化方案**，理由：

1. **大幅简化**: 9文件 → 3文件
2. **清晰API**: 仅3个公开类
3. **易维护**: 逻辑集中
4. **向后兼容**: 可选的兼容层
5. **架构一致**: 符合项目标准

**唯一风险**: 可能破坏直接使用内部类的代码（但这本身就是不推荐的用法）。
