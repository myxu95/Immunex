# Immunex 数据文件格式说明

## NPZ 文件格式

NPZ是NumPy的压缩数组存档格式（`.npz`），用于高效存储和加载多维数组数据。

### 1. 变构分析相关 NPZ 文件

#### `correlation_matrix.npz`

**用途**: 存储残基接触动态相关性矩阵

**结构**:
```python
{
    'correlation_matrix': np.ndarray,  # 相关系数矩阵
        # Shape: (N_contacts, N_contacts)
        # Dtype: float64
        # 值范围: [-1.0, 1.0]
        # 示例大小: 41 MB (2294×2294矩阵)

    'contact_labels': np.ndarray,  # 接触对标签（需要pickle）
        # Shape: (N_contacts,)
        # Dtype: object
        # 内容: 字符串标签如 "RES123_A-RES456_B"
}
```

**用途说明**:
- 用于变构分析（allosteric communication）
- 行/列表示不同的残基接触对
- 矩阵元素 (i,j) 表示接触对i和接触对j之间的相关系数
- 正相关：两个接触同时形成/断开（协同运动）
- 负相关：一个接触形成时另一个断开（反协同运动）

**读取示例**:
```python
import numpy as np

data = np.load('correlation_matrix.npz', allow_pickle=True)
corr_matrix = data['correlation_matrix']  # (N, N) 相关矩阵
contact_labels = data['contact_labels']   # (N,) 接触标签

# 找高度正相关的接触对
high_corr_indices = np.where(corr_matrix > 0.7)
print(f"找到 {len(high_corr_indices[0])} 个高度相关的接触对")
```

---

#### `contact_time_series.npz`

**用途**: 存储接触的时间序列二值矩阵

**结构**:
```python
{
    'contact_matrix': np.ndarray,  # 接触时间序列矩阵
        # Shape: (N_contacts, N_frames)
        # Dtype: uint8 (0或1，节省内存)
        # 0 = 接触断开，1 = 接触形成
        # 示例: (2294, 1001) - 2294个接触对，1001帧
        # 大小: 2.2 MB

    'times': np.ndarray,  # 时间戳
        # Shape: (N_frames,)
        # Dtype: float64
        # 单位: ps (皮秒)
        # 示例: [0, 200, 400, ..., 200000]

    'contact_labels': np.ndarray,  # 接触对标签（需要pickle）
        # Shape: (N_contacts,)
        # Dtype: object
}
```

**用途说明**:
- 记录每个残基接触对在轨迹中的形成/断开状态
- 用于计算correlation_matrix（皮尔逊相关系数）
- 用于变构路径分析
- 节省内存：uint8（1字节）而非float64（8字节）

**读取示例**:
```python
import numpy as np

data = np.load('contact_time_series.npz', allow_pickle=True)
contact_matrix = data['contact_matrix']  # (N_contacts, N_frames)
times = data['times']                     # (N_frames,)
labels = data['contact_labels']           # (N_contacts,)

# 计算接触频率
contact_frequencies = contact_matrix.mean(axis=1)  # (N_contacts,)

# 找到最稳定的接触（频率>80%）
stable_contacts = np.where(contact_frequencies > 0.8)[0]
print(f"稳定接触: {labels[stable_contacts]}")
```

---

#### `mosaic_matrix_reordered.npz`

**用途**: MoSAIC聚类分析重排序后的矩阵

**结构**:
```python
{
    'matrix': np.ndarray,  # 重排序后的相关矩阵
        # Shape: (N_contacts, N_contacts)
        # Dtype: float64
        # 按聚类结果重新排列

    'labels': np.ndarray,  # 聚类标签
        # Shape: (N_contacts,)
        # Dtype: int
        # 每个接触对所属的聚类ID

    'reorder_indices': np.ndarray,  # 重排序索引
        # Shape: (N_contacts,)
        # Dtype: int
        # 原始索引到新索引的映射
}
```

**用途说明**:
- 用于可视化变构社区（allosteric communities）
- 聚类后的接触对会在矩阵中聚集成块状结构
- 便于识别功能相关的残基组

---

### 2. 轨迹偏移 NPZ 文件

#### `.{trajectory_name}.xtc_offsets.npz`

**用途**: MDAnalysis缓存文件，用于快速随机访问XTC轨迹帧

**结构**:
```python
{
    'offsets': np.ndarray,  # 帧偏移量（字节位置）
        # Shape: (N_frames,)
        # Dtype: int64
        # 每一帧在XTC文件中的起始字节位置
        # 示例: [0, 560548, 1139268, 1724076, ...]
        # 大小: 8.9 KB (1001帧)

    'size': np.int64,  # XTC文件总大小（字节）
        # Scalar
        # 用于验证文件是否被修改

    'ctime': np.float64,  # 文件创建时间戳
        # Scalar
        # Unix时间戳

    'n_atoms': np.int64,  # 原子总数
        # Scalar
        # 用于验证拓扑匹配
}
```

**用途说明**:
- **性能优化**: XTC是压缩格式，顺序读取快但随机访问慢
- **缓存机制**: 第一次读取时生成offsets，后续访问直接跳转到对应帧
- **加速效果**: 随机访问速度提升10-100倍
- **自动生成**: MDAnalysis在第一次加载轨迹时自动创建
- **缓存失效**: 如果XTC文件被修改（size/ctime变化），会重新生成

**典型使用场景**:
```python
import MDAnalysis as mda

# 第一次加载 - 会生成.xtc_offsets.npz（慢）
u = mda.Universe('md.tpr', 'md_pbc.xtc')

# 随机访问第500帧（快 - 使用cached offsets）
u.trajectory[500]

# 后续加载 - 直接使用cached offsets（快）
u = mda.Universe('md.tpr', 'md_pbc.xtc')
```

---

## LOCK 文件格式

#### `.{trajectory_name}.xtc_offsets.lock`

**用途**: 防止多个进程同时生成offsets.npz文件

**结构**:
```
空文件（0字节）
```

**用途说明**:
- **并发控制**: 多个进程同时读取同一个轨迹时的锁机制
- **流程**:
  1. 进程A发现没有offsets.npz
  2. 进程A创建.lock文件
  3. 进程A生成offsets.npz
  4. 进程A删除.lock文件
  5. 进程B等待.lock消失后使用offsets.npz
- **孤儿锁**: 如果进程崩溃，.lock文件会残留（可以手动删除）

**清理孤儿锁**:
```bash
# 查找所有lock文件
find . -name "*.lock"

# 删除孤儿锁（确认没有进程在使用）
find . -name "*.lock" -delete
```

---

## NPZ vs CSV vs JSON

### 选择建议

| 格式 | 适用场景 | 优势 | 劣势 |
|------|---------|------|------|
| **NPZ** | 大型数值矩阵 | 压缩高效，读取快 | 不可人工查看 |
| **CSV** | 表格数据 | 通用，可用Excel打开 | 大文件加载慢 |
| **JSON** | 元数据，配置 | 结构化，易读 | 不适合数值数组 |

### 示例对比

同样的2294×2294相关矩阵：
- **NPZ压缩**: 41 MB
- **CSV**: ~200 MB
- **JSON**: ~250 MB

结论：**NPZ节省空间75-85%，加载速度快100倍**

---

## 常见操作

### 查看NPZ文件内容

```python
import numpy as np

# 列出所有key
data = np.load('file.npz')
print(f"Keys: {list(data.keys())}")

# 查看每个数组的形状和类型
for key in data.keys():
    print(f"{key}: shape={data[key].shape}, dtype={data[key].dtype}")
```

### 合并多个NPZ文件

```python
import numpy as np

# 读取多个文件
data1 = np.load('task1/correlation_matrix.npz')
data2 = np.load('task2/correlation_matrix.npz')

# 合并（例如取平均）
avg_matrix = (data1['correlation_matrix'] + data2['correlation_matrix']) / 2

# 保存合并结果
np.savez_compressed('merged.npz', correlation_matrix=avg_matrix)
```

### NPZ转CSV（用于可视化）

```python
import numpy as np
import pandas as pd

# 读取NPZ
data = np.load('correlation_matrix.npz', allow_pickle=True)
matrix = data['correlation_matrix']
labels = data['contact_labels']

# 转换为DataFrame
df = pd.DataFrame(matrix, index=labels, columns=labels)

# 保存为CSV（警告：可能很大）
df.to_csv('correlation_matrix.csv')
```

### 检查文件完整性

```python
import numpy as np
from pathlib import Path

def check_npz_integrity(npz_file):
    """检查NPZ文件是否完整可读"""
    try:
        data = np.load(npz_file)
        print(f"✓ {npz_file} 完整")
        print(f"  Keys: {list(data.keys())}")
        return True
    except Exception as e:
        print(f"✗ {npz_file} 损坏: {e}")
        return False

# 批量检查
for npz_file in Path('output/allostery_analysis').rglob('*.npz'):
    check_npz_integrity(npz_file)
```

---

## 数据压缩建议

### 使用压缩保存

```python
import numpy as np

# ✓ 推荐：使用压缩（文件更小）
np.savez_compressed('data.npz', matrix=large_array)

# ✗ 不推荐：不压缩（文件更大）
np.savez('data.npz', matrix=large_array)
```

### 选择合适的数据类型

```python
# 二值数据：使用uint8而非float64
contact_matrix = np.array([[0, 1, 0], [1, 1, 0]], dtype=np.uint8)
# 节省内存: 1 byte vs 8 bytes per element

# 小范围整数：使用int32而非int64
small_ints = np.array([1, 2, 3], dtype=np.int32)
# 节省内存: 4 bytes vs 8 bytes per element
```

---

## 故障排查

### 问题1: "ValueError: Object arrays cannot be loaded when allow_pickle=False"

**原因**: NPZ包含Python对象（如字符串数组）

**解决**:
```python
# 添加allow_pickle=True参数
data = np.load('file.npz', allow_pickle=True)
```

### 问题2: 孤儿锁导致等待

**症状**: MDAnalysis加载轨迹时卡住

**解决**:
```bash
# 删除锁文件
rm .*.xtc_offsets.lock

# 或删除整个缓存（会重新生成）
rm .*.xtc_offsets.*
```

### 问题3: 缓存文件与轨迹不匹配

**症状**: "Trajectory has been modified" 警告

**解决**:
```bash
# 删除旧缓存，重新生成
rm .md_pbc.xtc_offsets.npz
rm .md_pbc.xtc_offsets.lock
```

### 问题4: NPZ文件损坏

**症状**: "Failed to read array" 错误

**解决**:
```python
# 重新生成（如果有原始数据）
import numpy as np
# ... 重新计算 ...
np.savez_compressed('file.npz', matrix=new_matrix)
```

---

## 参考

- NumPy NPZ文档: https://numpy.org/doc/stable/reference/generated/numpy.savez.html
- MDAnalysis缓存机制: https://docs.mdanalysis.org/stable/documentation_pages/coordinates/XTC.html
- Immunex变构分析: `immunex/analysis/allostery/contact_correlation.py`
