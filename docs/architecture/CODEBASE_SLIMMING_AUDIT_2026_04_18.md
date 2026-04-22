# 代码库瘦身审计（2026-04-18）

## 目的

本审计用于识别仓库中仍然存在的冗余脚本、旧实现和边缘工具链，目标是：

1. 明确当前产品主线实际依赖的代码面
2. 找出与现主线重复、已经被替代或长期游离的脚本
3. 为后续“删除 / 归档 / 保留”提供可审阅依据

本次审计**不直接删除文件**，只给出分级结论。

## 审计方法

本次判断基于三类证据：

1. **活动主线引用**
   - 检查 `immunex/`、`scripts/`、`docs/`、`examples/`、根 `README.md`
   - 看文件名或模块是否仍被当前主线显式引用
2. **架构重叠**
   - 判断该文件是否已被当前 `analysis -> node -> pipeline -> cli` 主线替代
3. **时间线辅助**
   - Linux 下稳定的“创建时间”通常不可用
   - 本文使用**文件修改时间（mtime）**作为辅助，不单独作为删除依据

审计时间：
- `2026-04-18 14:41:41 CST`

## 当前主线判断

当前仓库的活动产品面，已经明显集中在以下路径：

- `immunex/analysis/`
- `immunex/pipeline/analysis_pipelines.py`
- `immunex/pipeline/nodes/*`
- `immunex/cli/commands/*`
- `scripts/build_interaction_demo_html.py`
- `scripts/build_comparison_report_html.py`

其中近两周明确新增或强化的主线包括：

- `RRCS`
- `compare`
- `report interaction`
- `interface clustering`
- `nodes/` 子包化重构

因此，凡是不进入上述主线、且又与这些能力重复的文件，都应优先视作瘦身候选。

---

## A. 高置信度冗余候选

这批文件要么已经被当前主线替代，要么明显属于旧实现残留。  
建议优先处理，处理方式可以是：

- 直接删除
- 或下沉到 `docs/archive/` / `development/archived_*` 对应说明中

### 1. `immunex/analysis/trajectory/rmsd_analyzer.py`

判断：
- 当前主线 RMSD 已由：
  - `immunex/analysis/trajectory/rmsd_refactored.py`
  - `immunex/pipeline/nodes/trajectory/rmsd_node.py`
  - `immunex/cli/commands/rmsd.py`
  承接
- 活动树中未发现对 `rmsd_analyzer` 的实际引用

结论：
- **高置信度冗余**

建议：
- 删除或归档

---

### 2. `immunex/analysis/trajectory/rmsd_convergence.py`

判断：
- 目前仅在 `quality_assessment_pipeline.py` 这条旧式质量支线中使用
- 不在当前产品主入口中暴露
- 与当前 `PreprocessQualityPipeline + report` 主线存在能力重叠

结论：
- **高置信度旧支线候选**

建议：
- 连同 `quality_assessment_pipeline.py` 一并评估
- 若确认不再作为单独产品面保留，可整体归档

---

### 3. `immunex/pipeline/quality_assessment_pipeline.py`

判断：
- 这是一个旧风格 pipeline，直接在 pipeline 文件中自行装配 analyzer
- 不符合当前“analysis + node + pipeline”主线模式
- 当前 CLI 并未以此作为主要入口

结论：
- **高置信度旧 pipeline 候选**

建议：
- 若无外部依赖，优先归档或删除

---

### 4. `immunex/pipeline/standard_pipelines.py`

判断：
- 目前仅剩 `PreprocessOnlyPipeline`、`PreprocessQualityPipeline`
- 其中前者价值有限，后者与主分析流水线的角色边界并不清晰
- 当前真正的专题能力已经集中在 `analysis_pipelines.py`

结论：
- **高置信度结构重复候选**

建议：
- 评估是否将仍有价值的两个类迁入 `analysis_pipelines.py`
- 完成迁入后移除本文件

---

### 5. `immunex/pipeline/batch_workflow.py`

判断：
- 仍被：
  - `immunex/pipeline/__init__.py`
  - `immunex/cluster/slurm_generator.py`
  引用
- 但其风格偏旧式包装层，且与 `batch_executor.py`、CLI batch 主线存在边界重叠

结论：
- **高优先级重构候选**
- 不是立即删除项，但已明显偏旧

建议：
- 先把 `discover_md_tasks / process_md_tasks / check_task_status` 向新的 batch 主线收口
- 再移除旧入口

---

### 6. `immunex/pipeline/batch_executor.py`

判断：
- 仍被 `pipeline/__init__.py` 暴露
- 仍承担现 batch 主线内部角色
- 但它和 `batch_workflow.py`、CLI `batch.py` 的边界需要进一步清理

结论：
- **高优先级结构收口候选**
- 不是立即删除项

建议：
- 暂保留，但应纳入 batch 体系收口计划

---

## B. 低引用、边缘工具链候选

这批文件不一定“错”，但它们明显不在当前主分析产品面中心。  
如果目标是让仓库更聚焦，可以考虑下沉或移出默认公开面。

### 1. `immunex/cli/batch_pdb.py`

判断：
- 仍被 `immunex/cli/commands/pdb.py` 引用
- 功能上属于独立 PDB 处理工具链，不属于当前 MD 主分析主线

结论：
- **边缘工具链**

建议：
- 若保留，建议明确定位为“结构预处理辅助工具”
- 若要极限瘦身，可迁出主产品面或后续归档

---

### 2. `immunex/utils/pdb_downloader.py`

判断：
- 仅在 `immunex/utils/__init__.py` 暴露
- 没有进入当前主分析流程
- 更像独立工具

结论：
- **边缘工具链**

建议：
- 不建议直接删
- 但建议从默认产品面降级，或迁入更明确的 `tools/` / `development/` 体系

---

### 3. `scripts/concatenate_multimodel_pdbs.py`

判断：
- 仅在示例中出现
- 与当前主分析主线无直接关系
- 明显属于 standalone utility

结论：
- **边缘工具脚本**

建议：
- 可保留
- 但建议下沉为工具脚本，不作为主产品面的一部分强调

---

## C. 可保留但应重新分层的模块

这批文件不是冗余，但当前归属不够理想，容易造成“仓库看起来过胖”。

### 1. `immunex/analysis/free_energy/fel_calculator.py`
### 2. `immunex/analysis/free_energy/fel_visualizer.py`

判断：
- 形成一条相对独立的小功能支线
- 当前不在主 CLI 产品面中
- 不与主线直接冲突

结论：
- **保留**
- 但应明确是否属于公开产品路线

建议：
- 如果近期不会进入产品面，可考虑降低曝光度

---

### 3. `immunex/analysis/structure/contact_map.py`

判断：
- 是独立结构分析工具，不和主线冲突
- 但目前不进入主报告 / 主 CLI 流程

结论：
- **保留**

建议：
- 作为结构分析辅助能力保留即可，不建议优先清理

---

### 4. `immunex/analysis/quality/*`

重点包括：
- `batch_tracker.py`
- `md_completeness.py`
- `post_pbc_validator.py`
- `quality_reporter.py`
- `convergence_checker.py`
- `structure_validator.py`

判断：
- 这批文件内部是自洽的一套质量子系统
- 但当前产品主线并没有完整围绕它展开
- 和现 `quality` CLI、`PreprocessQualityPipeline`、结果汇总层之间有一定重叠

结论：
- **不建议单文件级删除**
- 应按“质量子系统是否保留为正式能力”整体判断

建议：
- 若保留：需要进一步产品化收口
- 若不保留：应整体归档，而不是零碎删文件

---

## D. 当前不建议动的核心主线

以下模块虽然有些引用数不高，但已经进入当前主线，不应仅凭低引用删除：

- `immunex/analysis/interactions/rrcs.py`
- `immunex/analysis/comparison/*`
- `immunex/analysis/conformation/interface_clustering.py`
- `immunex/analysis/allostery/normal_mode.py`
- `immunex/analysis/topology/biological_identity.py`
- `scripts/build_interaction_demo_html.py`
- `scripts/build_comparison_report_html.py`
- `immunex/pipeline/nodes/*` 子包化后的节点
- `immunex/cli/commands/*` 当前活动命令

原因：
- 这批是当前 2026-04-08 之后产品表面和新主线的核心组成部分

---

## 初步瘦身建议顺序

### 第一批：立即处理

建议优先评估并处理：

1. `immunex/analysis/trajectory/rmsd_analyzer.py`
2. `immunex/pipeline/quality_assessment_pipeline.py`
3. `immunex/pipeline/standard_pipelines.py`
4. `immunex/pipeline/batch_workflow.py`（先重构后移除）

### 第二批：边缘工具链降级

1. `immunex/cli/batch_pdb.py`
2. `immunex/utils/pdb_downloader.py`
3. `scripts/concatenate_multimodel_pdbs.py`

### 第三批：专题支线做整体决策

1. `immunex/analysis/quality/*`
2. `immunex/analysis/free_energy/*`
3. `immunex/analysis/structure/contact_map.py`

---

## 建议的处理原则

后续正式瘦身时，建议遵循：

1. **优先删“被替代的旧实现”**
2. **不要先删“低引用但仍自洽的子系统”**
3. **边缘工具先降级，再决定是否删**
4. **对 batch/quality 这种系统级支线，按整条链处理，不按单文件处理**

---

## 当前审计结论

本次扫描后，仓库里最值得优先处理的冗余来源，不再是明显的 `legacy/` 或旧 SLURM 脚本，而是：

- 旧式 `pipeline` 包装层
- 被新主线替代的 RMSD/质量旧实现
- 不属于当前主产品面的边缘工具链

建议下一步动作：

1. 先人工审阅本文件
2. 选定“第一批立即处理清单”
3. 再做一轮**真正的删除/归档提交**

