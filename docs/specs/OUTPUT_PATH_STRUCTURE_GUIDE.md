# Immunex 输出路径结构详解

## 一、当前输出路径结构

### 1.1 基本结构

```
output_base/                    # 输出根目录
└── {system_id}/               # 任务唯一标识（如：3d39_run2）
    ├── md_pbc.xtc             # PBC校正后的轨迹
    ├── rmsd.csv               # RMSD时间序列数据
    ├── md_quality_metrics.json # 质量指标（机器可读）
    ├── md_quality_report.md   # 质量报告（人类可读）
    └── md_convergence.txt     # 收敛性分析报告
```

### 1.2 实际示例

**单任务测试**：
```
output/test_pbc_single/
└── 3d39_run2/
    ├── md_pbc.xtc              # 614 MB - 校正后的轨迹
    ├── rmsd.csv                # 27 KB - RMSD数据
    ├── md_quality_metrics.json # 43 KB - 详细质量指标
    ├── md_quality_report.md    # 820 bytes - 质量报告
    └── md_convergence.txt      # 783 bytes - 收敛性报告
```

**批处理测试**：
```
output/manual_batch_test/
├── 3d39_run2/
│   ├── md_pbc.xtc
│   ├── rmsd.csv
│   ├── md_quality_metrics.json
│   ├── md_quality_report.md
│   └── md_convergence.txt
└── 3kps_run2/
    ├── md_pbc.xtc
    ├── rmsd.csv
    ├── md_quality_metrics.json
    ├── md_quality_report.md
    └── md_convergence.txt
```

---

## 二、路径组成部分详解

### 2.1 路径层级

```
/home/xumy/work/development/Immunex/output/test_pbc_single/3d39_run2/md_pbc.xtc
 \_____________________________________/ \_______________/ \______/ \________/
              项目根目录                      输出类型      任务ID   文件名
```

#### 层级说明

| 层级 | 含义 | 示例 | 说明 |
|------|------|------|------|
| **项目根** | Immunex项目路径 | `/home/xumy/.../Immunex/` | 固定前缀 |
| **输出类型** | 输出目录分类 | `output/test_pbc_single/` | 区分不同测试/分析类型 |
| **任务ID** | 体系唯一标识 | `3d39_run2/` | 对应输入任务 |
| **文件名** | 具体输出文件 | `md_pbc.xtc` | 分析结果文件 |

### 2.2 输出类型（第2层）

常见的输出目录分类：

```
output/
├── test_pbc_single/          # 单任务PBC测试
├── manual_batch_test/        # 批处理测试
├── allostery_analysis/       # 变构分析
├── cdr_rmsd_exact_analysis/  # CDR RMSD分析
├── cdr3_analysis/            # CDR3专项分析
└── {user_specified}/         # 用户自定义输出
```

**命名规则**：
- 小写字母 + 下划线
- 描述分析类型或批次名称
- 避免使用日期或版本号

### 2.3 任务ID（第3层）

任务目录直接对应输入数据：

```
input/database/parallel2/3d39_run2/  → output/xxx/3d39_run2/
                         ^^^^^^^^^^              ^^^^^^^^^^
                         输入任务ID               输出任务ID（保持一致）
```

**好处**：
- ✅ 输入输出一一对应，易于追溯
- ✅ 支持批处理并行（不同任务目录独立）
- ✅ 避免文件名冲突

### 2.4 文件名（第4层）

遵循命名规范（见下文）。

---

## 三、文件命名规范

### 3.1 当前命名方式（PBC Pipeline）

| 文件名 | 内容 | 格式 | 大小 |
|--------|------|------|------|
| `md_pbc.xtc` | PBC校正后的轨迹 | XTC（二进制） | ~500-700 MB |
| `rmsd.csv` | RMSD时间序列 | CSV（文本） | ~20-30 KB |
| `md_quality_metrics.json` | 详细质量指标 | JSON | ~40-50 KB |
| `md_quality_report.md` | 人类可读报告 | Markdown | <1 KB |
| `md_convergence.txt` | 收敛性分析 | 文本 | <1 KB |

**命名规则**：
- `md_` 前缀：表示MD相关输出
- `_pbc` / `_quality` / `_convergence` 后缀：说明内容类型
- 小写 + 下划线（不用驼峰）

### 3.2 标准命名约定

#### 前缀规则

| 前缀 | 用途 | 示例 |
|------|------|------|
| `md_` | 轨迹相关 | `md_pbc.xtc`, `md_whole.xtc` |
| `rmsd_` | RMSD分析 | `rmsd_protein.xvg`, `rmsd_cdr3.csv` |
| `rmsf_` | RMSF分析 | `rmsf_by_residue.csv` |
| `contact_` | 接触分析 | `contact_frequency.csv` |
| `angle_` | 角度分析 | `angle_twist.xvg` |
| `energy_` | 能量分析 | `energy_quality.json` |

#### 后缀规则

| 后缀 | 含义 | 示例 |
|------|------|------|
| `_pbc` | PBC处理后 | `md_pbc.xtc` |
| `_processed` | 预处理后 | `md_processed.xtc` |
| `_aligned` | 对齐后 | `md_aligned.xtc` |
| `_quality` | 质量相关 | `md_quality_metrics.json` |
| `_convergence` | 收敛性 | `md_convergence.txt` |
| `_summary` | 汇总统计 | `rmsd_summary.json` |
| `_by_residue` | 按残基 | `rmsf_by_residue.csv` |
| `_trajectory` | 时间序列 | `rmsd_trajectory.csv` |

---

## 四、标准vs实际对比

### 4.1 标准结构（规划中）

```
output/
└── 3d39_run2/
    ├── preprocessing/          # 📁 预处理结果
    │   ├── md_pbc.xtc
    │   ├── md_center.xtc
    │   └── preprocessing.log
    │
    ├── quality/               # 📁 质量控制
    │   ├── md_quality_metrics.json
    │   ├── md_quality_report.md
    │   └── md_convergence.txt
    │
    ├── analysis/              # 📁 分析结果
    │   ├── rmsd/
    │   │   ├── rmsd.csv
    │   │   └── rmsd_summary.json
    │   ├── rmsf/
    │   └── angles/
    │
    └── plots/                 # 📁 可视化
        └── rmsd_overview.png
```

### 4.2 当前结构（简化版）

```
output/
└── 3d39_run2/
    ├── md_pbc.xtc              # ← 预处理
    ├── md_quality_metrics.json # ← 质量控制
    ├── md_quality_report.md    # ← 质量控制
    ├── md_convergence.txt      # ← 质量控制
    └── rmsd.csv                # ← 分析结果
```

### 4.3 对比分析

| 方面 | 标准结构 | 当前结构 | 优缺点 |
|------|---------|---------|--------|
| **组织性** | 分类清晰（子目录） | 扁平结构 | 简化：易查找<br>标准：更有序 |
| **扩展性** | 易添加新分析类型 | 文件多了会混乱 | 简化：当前够用<br>标准：长期更好 |
| **可读性** | 需要进入子目录 | 一目了然 | 简化：快速浏览<br>标准：专业 |
| **兼容性** | 符合大型项目规范 | 适合小规模使用 | 简化：学习成本低<br>标准：工业标准 |

---

## 五、路径设计原则

### 5.1 核心原则

#### 1. 输入输出对应
```python
input_path  = "input/database/parallel2/3d39_run2/prod/md.xtc"
output_path = "output/pbc_analysis/3d39_run2/md_pbc.xtc"
                                     ^^^^^^^^^^
                                     保持任务ID一致
```

#### 2. 避免硬编码
```python
# ✗ 错误：硬编码路径
output = "/home/xumy/output/result.xtc"

# ✓ 正确：使用变量和拼接
output = Path(args.output_dir) / task.system_id / "md_pbc.xtc"
```

#### 3. 相对路径 vs 绝对路径

**开发测试**：
```bash
output/test_pbc_single/3d39_run2/      # 相对路径（简洁）
```

**生产环境**：
```bash
/data/projects/immunex_analysis/batch_2026_03/3d39_run2/  # 绝对路径（明确）
```

#### 4. 批次管理

```
output/
├── batch_2026_03_15/          # 按日期批次
│   ├── 3d39_run2/
│   ├── 3kps_run2/
│   └── batch_summary.csv
│
├── rest2_sampling/            # 按实验类型批次
│   ├── 1bd2_run2_sampling/
│   └── ...
│
└── test_runs/                 # 测试批次
    └── ...
```

---

## 六、实际使用场景

### 6.1 单任务指定输出

```python
from immunex.core.context import PipelineContext
from immunex.pipeline import PreprocessQualityPipeline

context = PipelineContext(
    system_id="3d39_run2",
    topology="input/database/parallel2/3d39_run2/prod/md.tpr",
    trajectory_raw="input/database/parallel2/3d39_run2/prod/md.xtc",
    output_dir="output/my_analysis/3d39_run2",
)

pipeline = PreprocessQualityPipeline()
pipeline.execute(context)
```

**生成结构**：
```
output/my_analysis/3d39_run2/
├── md_pbc.xtc
├── rmsd.csv
└── ...
```

### 6.2 批处理自动路径

当前推荐使用 `imn batch ...` 作为标准批处理入口：

```bash
imn batch preprocess input/database/parallel2 \
    --workers 8 \
    --output-dir output/batch_2026_03_17
```

**生成结构**：
```
output/batch_2026_03_17/
├── 3d39_run2/
│   ├── md_processed.xtc
│   ├── md_processed_converted.pdb
│   └── ...
├── 3kps_run2/
│   ├── md_processed.xtc
│   └── ...
└── ...
```

### 6.3 CLI命令输出

```bash
imn report interaction \
    --base-dir output/interaction_case_1OGA \
    --system-id 1OGA_sd_run2 \
    --bsa-root output/bsa_demo_1OGA \
    --rmsf-root output/rmsf_demo_1OGA \
    --identity-root output/identity_demo_1OGA \
    --cluster-root output/interface_cluster_demo_1OGA_s5_v3

# 自动生成：
# output/interaction_case_1OGA/
# └── overview/
#     ├── interaction_report_demo.html
#     ├── bsa/
#     ├── cluster/
#     └── ...
```

---

## 七、路径管理最佳实践

### 7.1 使用PathManager

```python
from immunex.utils import PathManager

# 初始化路径管理器
pm = PathManager(
    base_dir="output/my_project",
    system_id="3d39_run2"
)

# 自动创建标准子目录
pbc_output = pm.get_preprocessing_path("md_pbc.xtc")
# → output/my_project/3d39_run2/preprocessing/md_pbc.xtc

rmsd_output = pm.get_analysis_path("rmsd/rmsd_protein.xvg")
# → output/my_project/3d39_run2/analysis/rmsd/rmsd_protein.xvg
```

### 7.2 避免路径拼接错误

```python
# ✗ 错误：字符串拼接（不跨平台）
output = base_dir + "/" + system_id + "/" + filename

# ✓ 正确：使用pathlib
from pathlib import Path
output = Path(base_dir) / system_id / filename

# ✓ 正确：确保目录存在
output.parent.mkdir(parents=True, exist_ok=True)
```

### 7.3 清理临时文件

```python
from immunex.utils import CleanupManager

cleanup = CleanupManager(output_dir)

# 标记临时文件
cleanup.mark_temporary("md_center.xtc")
cleanup.mark_temporary("md_whole.xtc")

# 保留最终结果
cleanup.mark_permanent("md_pbc.xtc")

# 自动清理
cleanup.cleanup_temporary_files()
```

---

## 八、常见问题

### Q1: 输出目录太多，如何管理？

**解决方案**：使用批次目录分类

```
output/
├── production_runs/          # 生产数据
│   └── batch_2026_03/
├── test_runs/               # 测试数据
│   └── test_pbc_single/
└── archived/                # 归档数据
    └── batch_2025_12/
```

### Q2: 文件名太长/太复杂？

**当前命名**（简洁）：
```
md_pbc.xtc
rmsd.csv
```

**标准命名**（详细）：
```
md_pbc_2step_processed.xtc
rmsd_backbone_protein_trajectory.csv
```

**建议**：
- 小规模项目：用简洁命名
- 大规模项目：用详细命名

### Q3: 如何追溯输入输出关系？

**方法1：保持任务ID一致**
```
input/.../3d39_run2/  →  output/.../3d39_run2/
```

**方法2：保存元数据**
```json
// output/3d39_run2/metadata.json
{
  "system_id": "3d39_run2",
  "input_topology": "input/database/parallel2/3d39_run2/prod/md.tpr",
  "input_trajectory": "input/database/parallel2/3d39_run2/prod/md.xtc",
  "processed_at": "2026-03-17T15:00:00",
  "pipeline": "PreprocessQualityPipeline",
  "version": "1.0.0"
}
```

### Q4: 批处理输出如何汇总？

**自动生成汇总文件**：
```
output/batch_2026_03/
├── 3d39_run2/
├── 3kps_run2/
└── batch_summary.csv         # ← 汇总所有任务的关键指标
```

**batch_summary.csv 示例**：
```csv
system_id,pbc_grade,rmsd_mean,rmsd_converged,output_path
3d39_run2,B,2.621,True,output/batch_2026_03/3d39_run2
3kps_run2,B,3.091,False,output/batch_2026_03/3kps_run2
```

---

## 九、总结

### 当前路径结构特点

✅ **优点**：
- 扁平结构，易于查找
- 文件少，一目了然
- 适合PBC+RMSD单一流程

⚠️ **限制**：
- 不适合复杂多步分析
- 文件多了会混乱
- 缺少分类组织

### 建议

**短期**（当前够用）：
- 保持简单扁平结构
- 文件数量<10个时很高效

**长期**（扩展时）：
- 迁移到标准子目录结构
- 使用PathManager管理路径
- 添加metadata.json追溯

### 快速参考

```python
# 当前输出结构（简化版）
output_dir/
└── {system_id}/
    ├── md_pbc.xtc              # 主输出
    ├── rmsd.csv                # 分析数据
    ├── md_quality_metrics.json # 质量指标
    ├── md_quality_report.md    # 人类可读报告
    └── md_convergence.txt      # 收敛性分析

# 标准输出结构（完整版）
output_dir/
└── {system_id}/
    ├── preprocessing/          # PBC处理
    ├── quality/               # 质量控制
    ├── analysis/              # 分析结果
    ├── plots/                 # 可视化
    └── metadata.json          # 元数据
```

希望这个指南能帮助你理解Immunex的输出路径结构！
