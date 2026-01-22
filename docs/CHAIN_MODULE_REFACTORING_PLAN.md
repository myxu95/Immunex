# Chain标准化模块重构计划

## 一、现状分析

### 现有模块

| 模块 | 大小 | 主要功能 | 依赖 |
|------|------|----------|------|
| `PDBSequenceExtractor` | 11K | 从PDB提取序列 | MDAnalysis |
| `IntelligentChainIdentifier` | 16K | 智能识别链类型 | PDBSequenceExtractor, ANARCIWrapper |
| `PDBChainStandardizer` | 26K | 标准化链ID | IntelligentChainIdentifier (可选) |
| `ShortestChainDetector` | 14K | GRO文件最短链检测 | GROMACS |

### 依赖关系

```
PDBSequenceExtractor (基础)
    ↓
IntelligentChainIdentifier (类型识别)
    ↓ (可选)
PDBChainStandardizer (ID标准化)

ShortestChainDetector (独立，GRO文件)
```

### 识别的问题

1. **功能重叠**
   - PDBChainStandardizer和IntelligentChainIdentifier都有链检测逻辑
   - 长度阈值在多处硬编码

2. **标准不一致**
   ```
   PDBChainStandardizer:
   A=HLA-α, B=β2m, C=peptide, D=TCR-α, E=TCR-β

   IntelligentChainIdentifier:
   标准顺序：A=HLA-α, B=β2m, C=peptide, D=TCR-α, E=TCR-β
   （实际一致，但文档不清晰）
   ```

3. **GRO vs PDB分离**
   - ShortestChainDetector只处理GRO
   - 其他模块只处理PDB
   - 缺乏统一的抽象接口

4. **缺乏统一配置**
   - 长度阈值散布在多个文件
   - pHLA-TCR体系的参数硬编码

5. **文档和测试不足**
   - 缺少使用示例
   - 缺少单元测试
   - 模块间关系不清晰

## 二、改进方案

### 2.1 架构设计原则

**单一职责原则**：
- 序列提取 → PDBSequenceExtractor
- 类型识别 → IntelligentChainIdentifier
- ID标准化 → PDBChainStandardizer
- GRO处理 → ShortestChainDetector

**依赖倒置**：
- 引入抽象接口（Protocol）
- 减少具体类之间的耦合

**配置外部化**：
- 创建ChainConfig类统一管理阈值

### 2.2 新增模块

#### 1. `chain_config.py` - 统一配置

```python
@dataclass
class pHLATCRChainConfig:
    """pHLA-TCR complex chain configuration."""

    # Length thresholds
    peptide_min_length: int = 5
    peptide_max_length: int = 25
    beta2m_min_length: int = 95
    beta2m_max_length: int = 105
    tcr_alpha_min_length: int = 115
    tcr_alpha_max_length: int = 215
    tcr_beta_min_length: int = 209
    tcr_beta_max_length: int = 252
    hla_alpha_min_length: int = 260

    # Standard chain assignment
    standard_chain_order: Dict[str, str] = field(default_factory=lambda: {
        'HLA_alpha': 'A',
        'beta2m': 'B',
        'peptide': 'C',
        'TCR_alpha': 'D',
        'TCR_beta': 'E'
    })

    # ANARCI settings
    use_anarci: bool = True
    anarci_scheme: str = 'imgt'
```

#### 2. `chain_analyzer.py` - 统一分析接口

```python
from typing import Protocol

class ChainAnalyzer(Protocol):
    """Abstract interface for chain analysis."""

    def analyze_chains(self, structure_file: str) -> List[ChainInfo]:
        """Analyze chains from structure file."""
        ...

class PDBChainAnalyzer(ChainAnalyzer):
    """PDB file chain analyzer."""

    def __init__(self, config: ChainConfig):
        self.extractor = PDBSequenceExtractor()
        self.identifier = IntelligentChainIdentifier(config)

    def analyze_chains(self, pdb_file: str) -> List[ChainInfo]:
        # 1. Extract sequences
        sequences = self.extractor.extract_sequences_from_pdb(pdb_file)

        # 2. Identify chain types
        identifications = self.identifier.identify_chains(sequences)

        # 3. Return structured results
        return self._build_chain_info(identifications)

class GROChainAnalyzer(ChainAnalyzer):
    """GRO file chain analyzer."""

    def __init__(self, topology_file: str):
        self.detector = ShortestChainDetector(topology_file)

    def analyze_chains(self, gro_file: str) -> List[ChainInfo]:
        # Implement GRO-specific logic
        ...
```

#### 3. `chain_standardizer_facade.py` - 统一入口

```python
class ChainStandardizerFacade:
    """
    Unified facade for chain standardization operations.

    Provides high-level API for common use cases:
    - PDB standardization
    - GRO shortest chain detection
    - Chain type identification
    """

    def __init__(self, config: Optional[ChainConfig] = None):
        self.config = config or pHLATCRChainConfig()
        self.pdb_analyzer = PDBChainAnalyzer(self.config)
        self.pdb_standardizer = PDBChainStandardizer(self.config)

    def standardize_pdb(
        self,
        input_pdb: str,
        output_pdb: str,
        mode: str = 'intelligent'
    ) -> StandardizationResult:
        """
        Standardize PDB chain IDs.

        Args:
            input_pdb: Input PDB file
            output_pdb: Output standardized PDB file
            mode: 'simple' (length-based) or 'intelligent' (type-based)
        """
        if mode == 'intelligent':
            # Use intelligent identifier
            identifications = self.pdb_analyzer.analyze_chains(input_pdb)
            return self.pdb_standardizer.standardize_by_type(
                input_pdb, output_pdb, identifications
            )
        else:
            # Simple length-based
            return self.pdb_standardizer.standardize_by_length(
                input_pdb, output_pdb
            )

    def detect_shortest_chain(
        self,
        gro_file: str,
        topology_file: str,
        output_dir: str
    ) -> str:
        """Detect shortest chain for PBC centering."""
        detector = ShortestChainDetector(gro_file, topology_file)
        return detector.generate_shortest_chain_index(output_dir)

    def identify_chain_types(self, pdb_file: str) -> Dict[str, ChainIdentification]:
        """Identify chain types without standardization."""
        return self.pdb_analyzer.analyze_chains(pdb_file)
```

### 2.3 模块重构

#### IntelligentChainIdentifier 改进

**当前问题**：
- 长度阈值硬编码
- ANARCI失败后的fallback不够健壮

**改进**：
```python
class IntelligentChainIdentifier:
    def __init__(self, config: ChainConfig, use_anarci: bool = True):
        self.config = config
        self.use_anarci = use_anarci
        self.extractor = PDBSequenceExtractor()

        if use_anarci:
            self.anarci = ANARCIWrapper(
                scheme=config.anarci_scheme,
                allow_fallback=True
            )

    def identify_chains(
        self,
        pdb_file: str
    ) -> List[ChainIdentification]:
        """
        Identify chain types using multi-strategy approach.

        Strategy priority:
        1. High-confidence length-based (peptide, beta2m)
        2. ANARCI for TCR chains
        3. Residual assignment (HLA-alpha)
        4. Fallback to pure length-based
        """
        sequences = self.extractor.extract_sequences_from_pdb(pdb_file)

        # Strategy 1: Identify peptide and beta2m by length
        results = self._identify_by_length(sequences)

        # Strategy 2: Use ANARCI for remaining chains
        if self.use_anarci:
            results = self._refine_with_anarci(results, sequences)

        # Strategy 3: Assign remaining by process of elimination
        results = self._assign_remaining(results)

        return results

    def _identify_by_length(self, sequences: Dict) -> List[ChainIdentification]:
        """Identify high-confidence chains by length."""
        identifications = []

        for chain_id, seq_info in sequences.items():
            length = seq_info['length']

            # Peptide identification (high confidence)
            if self.config.peptide_min_length <= length <= self.config.peptide_max_length:
                identifications.append(ChainIdentification(
                    chain_id=chain_id,
                    length=length,
                    sequence=seq_info['sequence'],
                    chain_type='peptide',
                    confidence=0.95,
                    method='length_based'
                ))

            # Beta2m identification (high confidence)
            elif self.config.beta2m_min_length <= length <= self.config.beta2m_max_length:
                identifications.append(ChainIdentification(
                    chain_id=chain_id,
                    length=length,
                    sequence=seq_info['sequence'],
                    chain_type='beta2m',
                    confidence=0.90,
                    method='length_based'
                ))

            else:
                # Unknown, need further analysis
                identifications.append(ChainIdentification(
                    chain_id=chain_id,
                    length=length,
                    sequence=seq_info['sequence'],
                    chain_type='unknown',
                    confidence=0.0,
                    method='pending'
                ))

        return identifications
```

#### PDBChainStandardizer 改进

**当前问题**：
- 两种模式混在一个类中
- 批处理逻辑和标准化逻辑耦合

**改进**：
```python
class PDBChainStandardizer:
    """
    Standardize PDB chain IDs based on identification results.

    Supports two modes:
    - Length-based: Order by residue count (legacy)
    - Type-based: Order by identified chain types (recommended)
    """

    def __init__(self, config: Optional[ChainConfig] = None):
        self.config = config or pHLATCRChainConfig()

    def standardize_by_type(
        self,
        input_pdb: str,
        output_pdb: str,
        identifications: List[ChainIdentification]
    ) -> StandardizationResult:
        """
        Standardize based on identified chain types.

        Uses config.standard_chain_order for mapping:
        peptide → C, beta2m → B, HLA_alpha → A, etc.
        """
        # Build chain mapping from identifications
        chain_mapping = self._build_type_mapping(identifications)

        # Apply mapping to PDB file
        self._renumber_chains(input_pdb, output_pdb, chain_mapping)

        return StandardizationResult(
            status='OK',
            chain_mapping=chain_mapping,
            method='type_based'
        )

    def standardize_by_length(
        self,
        input_pdb: str,
        output_pdb: str
    ) -> StandardizationResult:
        """Legacy length-based standardization."""
        # Extract chain lengths
        chains = self._analyze_chain_lengths(input_pdb)

        # Sort by length and assign standard IDs
        sorted_chains = sorted(chains, key=lambda x: x.residue_count)
        chain_mapping = {
            chain.chain_id: standard_id
            for chain, standard_id in zip(sorted_chains, ['C', 'B', 'D', 'E', 'A'])
        }

        self._renumber_chains(input_pdb, output_pdb, chain_mapping)

        return StandardizationResult(
            status='OK',
            chain_mapping=chain_mapping,
            method='length_based'
        )
```

### 2.4 新目录结构

```
aftermd/utils/
├── chain/
│   ├── __init__.py
│   ├── config.py                   # ChainConfig, pHLATCRChainConfig
│   ├── types.py                    # ChainInfo, ChainIdentification等数据类
│   ├── sequence_extractor.py       # PDBSequenceExtractor (重命名)
│   ├── chain_identifier.py         # IntelligentChainIdentifier (重命名)
│   ├── chain_standardizer.py      # PDBChainStandardizer (重命名)
│   ├── shortest_chain_detector.py  # ShortestChainDetector (保留)
│   ├── chain_analyzer.py           # 新增: PDBChainAnalyzer, GROChainAnalyzer
│   └── facade.py                   # 新增: ChainStandardizerFacade
```

## 三、实施步骤

### 阶段1: 准备工作（不破坏现有代码）

1. ✅ 分析现有模块功能和依赖
2. 📝 创建改进方案文档
3. ⏭️ 创建新的chain子包目录
4. ⏭️ 创建配置模块（chain/config.py）
5. ⏭️ 创建类型定义模块（chain/types.py）

### 阶段2: 重构核心模块

1. 重构PDBSequenceExtractor
   - 移动到chain/sequence_extractor.py
   - 保持原有API兼容

2. 重构IntelligentChainIdentifier
   - 移动到chain/chain_identifier.py
   - 接受ChainConfig参数
   - 改进多策略识别逻辑
   - 保持原有API兼容（通过wrapper）

3. 重构PDBChainStandardizer
   - 移动到chain/chain_standardizer.py
   - 分离length-based和type-based逻辑
   - 保持原有API兼容

4. 保留ShortestChainDetector
   - 移动到chain/shortest_chain_detector.py
   - 保持原有功能不变

### 阶段3: 新增功能

1. 创建ChainAnalyzer接口
   - chain/chain_analyzer.py
   - PDBChainAnalyzer实现
   - GROChainAnalyzer实现

2. 创建Facade统一入口
   - chain/facade.py
   - 高级API封装

3. 向后兼容层
   - 在原位置保留wrapper
   - 发出DeprecationWarning
   - 指向新位置

### 阶段4: 文档和测试

1. 创建使用文档
   - docs/CHAIN_STANDARDIZATION_GUIDE.md
   - 包含常见使用场景
   - 迁移指南

2. 添加单元测试
   - tests/test_chain_config.py
   - tests/test_chain_identifier.py
   - tests/test_chain_standardizer.py
   - tests/test_chain_facade.py

3. 添加集成测试
   - tests/integration/test_chain_workflow.py
   - 测试完整的标准化流程

### 阶段5: 迁移和清理

1. 更新内部模块引用
2. 更新scripts中的引用
3. 更新文档
4. 删除旧的wrapper（如果确认无人使用）

## 四、兼容性策略

### 向后兼容wrapper示例

```python
# aftermd/utils/intelligent_chain_identifier.py (旧位置)
import warnings
from .chain import IntelligentChainIdentifier as _NewIdentifier
from .chain import ChainIdentification

warnings.warn(
    "Importing from aftermd.utils.intelligent_chain_identifier is deprecated. "
    "Please use: from aftermd.utils.chain import IntelligentChainIdentifier",
    DeprecationWarning,
    stacklevel=2
)

IntelligentChainIdentifier = _NewIdentifier

__all__ = ['IntelligentChainIdentifier', 'ChainIdentification']
```

### 导入路径

**旧路径（保留6个月）**：
```python
from aftermd.utils import IntelligentChainIdentifier
from aftermd.utils import PDBChainStandardizer
```

**新路径（推荐）**：
```python
from aftermd.utils.chain import ChainStandardizerFacade
from aftermd.utils.chain import IntelligentChainIdentifier
from aftermd.utils.chain import pHLATCRChainConfig
```

## 五、使用示例

### 场景1: 简单PDB标准化

```python
from aftermd.utils.chain import ChainStandardizerFacade

facade = ChainStandardizerFacade()
result = facade.standardize_pdb(
    input_pdb="complex.pdb",
    output_pdb="complex_std.pdb",
    mode='intelligent'  # 或 'simple'
)

print(f"Standardization: {result.status}")
print(f"Chain mapping: {result.chain_mapping}")
```

### 场景2: 自定义配置

```python
from aftermd.utils.chain import ChainStandardizerFacade, pHLATCRChainConfig

# Custom configuration for non-standard complexes
config = pHLATCRChainConfig(
    peptide_max_length=30,  # Longer peptide
    use_anarci=False         # Disable ANARCI
)

facade = ChainStandardizerFacade(config)
result = facade.standardize_pdb("complex.pdb", "output.pdb")
```

### 场景3: 仅识别链类型

```python
from aftermd.utils.chain import ChainStandardizerFacade

facade = ChainStandardizerFacade()
identifications = facade.identify_chain_types("complex.pdb")

for ident in identifications:
    print(f"Chain {ident.chain_id}: {ident.chain_type} "
          f"(confidence: {ident.confidence:.2f})")
```

### 场景4: GRO文件最短链检测

```python
from aftermd.utils.chain import ChainStandardizerFacade

facade = ChainStandardizerFacade()
index_file = facade.detect_shortest_chain(
    gro_file="md.gro",
    topology_file="md.tpr",
    output_dir="./"
)

print(f"Index file created: {index_file}")
```

## 六、测试计划

### 单元测试

1. **config.py**
   - 测试默认配置加载
   - 测试自定义配置
   - 测试配置验证

2. **sequence_extractor.py**
   - 测试标准PDB解析
   - 测试修饰残基处理
   - 测试错误PDB处理

3. **chain_identifier.py**
   - 测试长度识别策略
   - 测试ANARCI识别
   - 测试fallback逻辑
   - 测试边界情况（4链、6链）

4. **chain_standardizer.py**
   - 测试type-based标准化
   - 测试length-based标准化
   - 测试链映射正确性

5. **shortest_chain_detector.py**
   - 测试GRO解析
   - 测试最短链检测
   - 测试index文件生成

### 集成测试

1. **完整workflow测试**
   - PDB下载 → 识别 → 标准化 → 验证
   - 多个真实PDB结构测试

2. **边界情况测试**
   - 非标准体系（4链、6链）
   - 缺失链
   - 非常规长度

3. **性能测试**
   - 批量处理100个PDB
   - 内存使用监控

## 七、时间估计

| 阶段 | 任务 | 预计时间 |
|------|------|----------|
| 阶段1 | 准备工作 | 2小时 |
| 阶段2 | 重构核心模块 | 6小时 |
| 阶段3 | 新增功能 | 4小时 |
| 阶段4 | 文档和测试 | 6小时 |
| 阶段5 | 迁移和清理 | 2小时 |
| **总计** | | **20小时** |

## 八、风险和缓解

| 风险 | 影响 | 缓解措施 |
|------|------|----------|
| 破坏现有代码 | 高 | 保留向后兼容wrapper |
| ANARCI依赖问题 | 中 | 提供fallback机制 |
| 性能下降 | 低 | 性能测试和优化 |
| 文档不同步 | 中 | 自动化文档生成 |

## 九、成功标准

1. ✅ 所有现有功能保持正常工作
2. ✅ 新API更简单易用
3. ✅ 单元测试覆盖率 > 80%
4. ✅ 完整的使用文档
5. ✅ 无性能退化
6. ✅ 代码重复减少 > 30%
