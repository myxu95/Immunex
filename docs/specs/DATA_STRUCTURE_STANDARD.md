# Immunex 数据结构标准规范

**版本**: 1.0
**日期**: 2026-03-16
**状态**: Draft

---

## 一、输入目录结构标准

### 1.1 标准结构（推荐）

```
input_base/
├── {pdb_id}/                    # PDB ID (4字符)
│   ├── {variant}/               # 变体/副本名称（可选）
│   │   ├── md.tpr              # 拓扑文件 (必需)
│   │   ├── md.xtc              # 原始轨迹 (必需)
│   │   ├── md.gro              # 结构文件 (必需)
│   │   ├── md.log              # 日志文件 (可选)
│   │   ├── md.edr              # 能量文件 (可选)
│   │   └── structure.pdb       # PDB结构 (可选)
│   └── ...
└── ...
```

**示例**:
```
/data/simulations/
├── 1ao7/
│   ├── standard/
│   │   ├── md.tpr
│   │   ├── md.xtc
│   │   └── md.gro
│   ├── rest2/
│   │   ├── md.tpr
│   │   ├── md.xtc
│   │   └── md.gro
│   └── replica_1/
│       └── ...
├── 1bd2/
│   └── standard/
│       └── ...
└── 1mi5/
    └── ...
```

### 1.2 简化结构（向后兼容）

```
input_base/
├── {task_name}/                 # 任务名称（可以是PDB ID或任意名称）
│   ├── md.tpr                  # 拓扑文件 (必需)
│   ├── md.xtc                  # 原始轨迹 (必需)
│   ├── md.gro                  # 结构文件 (必需)
│   └── prod/                   # 或在prod子目录中 (可选)
│       ├── md.tpr
│       ├── md.xtc
│       └── md.gro
└── ...
```

**示例**:
```
/data/simulations/
├── 1ao7_run1/
│   ├── md.tpr
│   ├── md.xtc
│   └── md.gro
├── 1bd2_mutant_A/
│   └── prod/
│       ├── md.tpr
│       ├── md.xtc
│       └── md.gro
└── ...
```

### 1.3 文件命名规范

| 文件类型 | 标准命名 | 备选命名 | 说明 |
|---------|---------|---------|------|
| 拓扑文件 | `md.tpr` | `topol.tpr`, `system.tpr` | GROMACS拓扑 |
| 原始轨迹 | `md.xtc` | `traj.xtc`, `md.trr` | 压缩轨迹（首选xtc） |
| 结构文件 | `md.gro` | `conf.gro`, `system.gro` | 最终结构 |
| PDB结构 | `structure.pdb` | `protein.pdb`, `complex.pdb` | 初始PDB |
| 日志文件 | `md.log` | - | GROMACS日志 |
| 能量文件 | `md.edr` | - | 能量轨迹 |

---

## 二、输出目录结构标准

### 2.1 标准输出目录结构

```
output_base/
├── {system_id}/                 # 与输入任务对应的唯一标识
│   ├── preprocessing/          # 预处理结果
│   │   ├── processed.xtc      # PBC处理后的轨迹
│   │   ├── aligned.xtc        # 对齐后的轨迹
│   │   ├── center.xtc         # 中间文件（临时）
│   │   ├── whole.xtc          # 中间文件（临时）
│   │   └── preprocessing.log  # 预处理日志
│   │
│   ├── quality/               # 质量控制结果
│   │   ├── completeness.json  # 完整性检查
│   │   ├── energy_quality.json # 能量质量
│   │   ├── structure_quality.json # 结构质量
│   │   └── quality_report.txt # 综合报告
│   │
│   ├── analysis/              # 分析结果
│   │   ├── rmsd/
│   │   │   ├── rmsd_protein.xvg
│   │   │   ├── rmsd_backbone.xvg
│   │   │   ├── rmsd_cdr3.xvg
│   │   │   └── rmsd_summary.json
│   │   │
│   │   ├── rmsf/
│   │   │   ├── rmsf_protein.xvg
│   │   │   ├── rmsf_by_residue.csv
│   │   │   └── rmsf_summary.json
│   │   │
│   │   ├── angles/
│   │   │   ├── docking_angles.csv
│   │   │   ├── twist_angle.xvg
│   │   │   ├── tilt_angle.xvg
│   │   │   ├── swing_angle.xvg
│   │   │   └── angles_summary.json
│   │   │
│   │   ├── interface/
│   │   │   ├── bsa_trajectory.csv
│   │   │   ├── contact_frequency.csv
│   │   │   ├── contact_atoms.csv
│   │   │   ├── contact_heatmap.png
│   │   │   └── interface_summary.json
│   │   │
│   │   ├── allostery/
│   │   │   ├── correlation_matrix.csv
│   │   │   ├── correlation_heatmap.png
│   │   │   ├── contact_labels.csv
│   │   │   └── allostery_summary.json
│   │   │
│   │   └── hydrogen_bonds/
│   │       ├── hbonds_trajectory.csv
│   │       └── hbonds_summary.json
│   │
│   ├── plots/                 # 可视化结果
│   │   ├── rmsd_overview.png
│   │   ├── rmsf_overview.png
│   │   ├── angles_timeseries.png
│   │   ├── interface_evolution.png
│   │   └── summary_dashboard.png
│   │
│   ├── indices/               # GROMACS索引文件
│   │   ├── protein.ndx
│   │   ├── tcr.ndx
│   │   ├── phla.ndx
│   │   └── interface.ndx
│   │
│   ├── context.json           # PipelineContext序列化
│   └── metadata.json          # 任务元数据
│
└── batch_summary/             # 批次级别汇总
    ├── batch_results.csv
    ├── batch_statistics.json
    ├── failed_tasks.txt
    └── batch_report.html
```

### 2.2 输出文件命名规范

#### 基本规则
1. **小写字母 + 下划线**: `rmsd_protein.xvg`（不使用驼峰或连字符）
2. **描述性名称**: 文件名应清楚表明内容
3. **标准扩展名**: `.xvg`（数据）, `.csv`（表格）, `.json`（结构化数据）, `.png`（图像）

#### 分析类型前缀

| 分析类型 | 前缀 | 示例 |
|---------|-----|------|
| RMSD | `rmsd_` | `rmsd_protein.xvg`, `rmsd_cdr3.xvg` |
| RMSF | `rmsf_` | `rmsf_by_residue.csv` |
| 角度 | `angle_`, `{type}_angle` | `twist_angle.xvg`, `docking_angles.csv` |
| 距离 | `distance_`, `dist_` | `distance_com.xvg` |
| 接触 | `contact_` | `contact_frequency.csv` |
| 氢键 | `hbonds_`, `hbond_` | `hbonds_trajectory.csv` |
| 能量 | `energy_` | `energy_quality.json` |

#### 后缀约定

| 后缀 | 含义 | 示例 |
|-----|------|------|
| `_trajectory` | 时间序列数据 | `rmsd_trajectory.xvg` |
| `_summary` | 汇总统计 | `rmsd_summary.json` |
| `_by_residue` | 按残基 | `rmsf_by_residue.csv` |
| `_heatmap` | 热图 | `contact_heatmap.png` |
| `_overview` | 总览图 | `rmsd_overview.png` |

---

## 三、PipelineContext 路径管理增强

### 3.1 新增方法

```python
class PipelineContext:
    def get_preprocessing_path(self, filename: str) -> str:
        """获取预处理文件路径"""
        return self.get_subdir_path("preprocessing", filename)

    def get_analysis_path(self, analysis_type: str, filename: str) -> str:
        """获取分析结果路径"""
        return self.get_subdir_path(f"analysis/{analysis_type}", filename)

    def get_plot_path(self, filename: str) -> str:
        """获取绘图文件路径"""
        return self.get_subdir_path("plots", filename)

    def get_quality_path(self, filename: str) -> str:
        """获取质量控制文件路径"""
        return self.get_subdir_path("quality", filename)

    def get_index_path(self, filename: str) -> str:
        """获取索引文件路径"""
        return self.get_subdir_path("indices", filename)

    def get_subdir_path(self, subdir: str, filename: str) -> str:
        """通用子目录路径获取"""
        if self.output_dir is None:
            self.output_dir = f"./results/{self.system_id}"

        subdir_path = Path(self.output_dir) / subdir
        subdir_path.mkdir(parents=True, exist_ok=True)

        return str(subdir_path / filename)

    def get_temp_path(self, filename: str) -> str:
        """获取临时文件路径（自动注册到temporary_files）"""
        temp_path = self.get_subdir_path("temp", filename)
        if temp_path not in self.temporary_files:
            self.temporary_files.append(temp_path)
        return temp_path
```

### 3.2 使用示例

```python
# 预处理节点
processed_traj = context.get_preprocessing_path("processed.xtc")

# RMSD节点
rmsd_output = context.get_analysis_path("rmsd", "rmsd_protein.xvg")

# 角度节点
angles_output = context.get_analysis_path("angles", "docking_angles.csv")

# 绘图
plot_file = context.get_plot_path("rmsd_overview.png")

# 临时文件（自动标记为可清理）
temp_index = context.get_temp_path("temp_protein.ndx")
```

---

## 四、DiscoveryReport 契约增强

### 4.1 多层级结构支持

```python
from pathlib import Path
from immunex.core import discover_tasks

report = discover_tasks(
    Path("/data/batch"),
    task_depth=2,
    required_files=["structure", "topology", "trajectory"],
)
```

### 4.2 发现规则优先级

1. **显式 task.yaml / immunex_task.yaml**
2. **nested结构规则** ({pdb_id}/{variant}/)
3. **flat结构规则** ({task_name}/)
4. **prod子目录规则** ({task_name}/prod/)

---

## 五、system_id 命名规范

### 5.1 标准格式

```
{pdb_id}_{variant}_{suffix}
```

**组成部分**:
- `pdb_id`: 4字符PDB ID（小写）
- `variant`: 变体/副本名称（可选）
- `suffix`: 额外后缀（可选）

### 5.2 示例

| 输入目录结构 | system_id |
|------------|-----------|
| `1ao7/standard/` | `1ao7_standard` |
| `1ao7/rest2/` | `1ao7_rest2` |
| `1bd2/replica_1/` | `1bd2_replica_1` |
| `1ao7_run1/` | `1ao7_run1` |
| `complex_A/` | `complex_A` |

### 5.3 命名规则

1. **唯一性**: system_id在批次中必须唯一
2. **可读性**: 使用下划线分隔，避免特殊字符
3. **一致性**: 同一PDB的不同变体使用统一前缀
4. **文件系统安全**: 仅使用字母、数字、下划线

---

## 六、实施计划

### Phase 1: 增强 PipelineContext (1小时)
- [ ] 添加路径管理方法
- [ ] 更新单元测试
- [ ] 更新文档

### Phase 2: 增强 DiscoveryReport 契约 (2小时)
- [ ] 实现多层级结构支持
- [ ] 补充 stage-specific required_files
- [ ] 实现system_id标准化生成
- [ ] 更新单元测试

### Phase 3: 更新 Pipeline Nodes (2小时)
- [ ] PreprocessNode使用新路径方法
- [ ] RMSDNode使用新路径方法
- [ ] 创建示例脚本
- [ ] 测试完整流程

### Phase 4: 文档和示例 (1小时)
- [ ] 更新快速开始指南
- [ ] 创建目录结构迁移指南
- [ ] 添加完整示例

---

## 七、向后兼容性

### 7.1 兼容策略

1. **公开主路径统一**: 默认使用 `discover_tasks(...)` 和 `DiscoveryReport`
2. **阶段化输入要求**: 通过 `required_files` 声明不同 pipeline 的输入契约
3. **legacy代码隔离**: 历史实现仅保留在 `immunex.legacy`

### 7.2 迁移建议

```python
# 旧方式（仍然可用）
output_file = os.path.join(output_dir, task_name, "rmsd.xvg")

# 新方式（推荐）
output_file = context.get_analysis_path("rmsd", "rmsd_protein.xvg")
```

---

## 八、最佳实践

### 8.1 输入数据组织

1. **使用nested结构处理多变体**: `{pdb_id}/{variant}/`
2. **使用描述性variant名称**: `standard`, `rest2`, `mutant_A`
3. **保持文件命名一致**: 统一使用 `md.tpr`, `md.xtc`, `md.gro`

### 8.2 输出数据组织

1. **使用PipelineContext路径方法**: 避免硬编码路径
2. **分类存储结果**: 按分析类型分目录
3. **保存context.json**: 便于结果追溯和重现
4. **及时清理临时文件**: 使用`get_temp_path()`自动管理

### 8.3 批处理

1. **使用 `discover_tasks(...)` 自动发现**: 减少手动配置
2. **统一输出根目录**: 便于结果汇总
3. **记录批次元数据**: 保存批次配置和统计

---

## 九、FAQ

### Q1: 如何处理不符合标准的现有数据？

**A**: 新的发现契约支持显式任务清单，先把不规则目录映射成标准任务，再统一进入后续 pipeline：

```python
from immunex.core import discover_tasks_from_list, discovery_report_to_contexts

task_list = [
    {
        "task_id": "custom_task",
        "structure": "legacy_inputs/my_structure.gro",
        "topology": "legacy_inputs/my_topol.tpr",
        "trajectory_raw": "legacy_inputs/my_traj.xtc",
    }
]

report = discover_tasks_from_list(task_list, source_root=base_dir)
contexts = discovery_report_to_contexts(report)
```

### Q2: 如何在不同项目间复用配置？

**A**: 使用PipelineContext序列化功能：

```python
# 保存配置
context.save("project_config.json")

# 其他项目加载
context = PipelineContext.load("project_config.json")
context.system_id = "new_task"  # 修改必要字段
```

### Q3: 如何实现自定义输出目录结构？

**A**: 覆盖`output_dir`参数，使用`get_subdir_path()`灵活组织：

```python
context.output_dir = "/custom/output/path"
custom_result = context.get_subdir_path("my_analysis/custom", "result.dat")
```

---

**变更历史**:
- 2026-03-16: 初版发布
