# GitHub Legacy代码清理指南

**日期**: 2026-04-22
**问题**: GitHub仓库包含大量历史legacy代码和大文件
**Git仓库大小**: 6.6GB

---

## 问题分析

### 1. Git历史中的大文件

最大的文件（前10）：
- `immunex/data/reference/hla/derived/class_i_extracellular.fasta` (4.3MB)
- `immunex/data/reference/hla/derived/class_i_metadata.csv` (4.0MB)
- `immunex/data/reference/hla/raw/class_i/B_prot.fasta` (1.6MB)
- `immunex/data/reference/hla/raw/class_i/A_prot.fasta` (1.4MB)
- `immunex/data/reference/hla/raw/class_i/C_prot.fasta` (1.3MB)
- `development/angle_analysis/Analysis-scripts/Angles/脚本说明.pdf` (953KB)
- `development/angle_analysis/Analysis-scripts/Angles/脚本说明.docx` (794KB)
- `docs/截屏2026-01-07 13.33.39.png` (767KB)
- `development/angle_analysis/Analysis-scripts/TCRdockingangle.png` (623KB)

### 2. Legacy代码问题

- 历史提交中包含已删除的legacy/目录
- 多次重构导致大量文件被删除又重新创建
- 项目名称从"aftermd"改为"immunex"

---

## 解决方案

### 方案1：保守方案 - 仅清理未来提交 ✅ 已完成

**适用场景**: 不想改变历史，只想确保未来干净

**操作**:
```bash
# 1. 更新.gitignore（已完成）
# 2. 提交当前清理工作
git add .
git commit -m "refactor: Phase 1 cleanup - remove legacy code and temporary files"
git push origin main
```

**优点**:
- 简单安全
- 不破坏历史
- 立即生效

**缺点**:
- Git仓库仍然是6.6GB
- 克隆仓库仍然很慢
- 历史中的legacy代码仍然存在

---

### 方案2：温和方案 - 创建干净的发布分支 ⭐ 推荐

**适用场景**: 想要干净的发布版本，但保留完整历史

**操作**:
```bash
# 1. 创建新的orphan分支（无历史）
git checkout --orphan clean-main

# 2. 添加当前所有文件
git add .

# 3. 创建初始提交
git commit -m "Initial commit: Immunex v1.0.0 - Clean codebase

- Complete TCR-pMHC MD analysis platform
- 15 CLI commands (preprocess, quality, rmsd, rmsf, contact, etc.)
- Interactive HTML reports with AI assistant
- Batch processing support
- Comprehensive documentation
"

# 4. 推送到GitHub（作为新的默认分支）
git push origin clean-main

# 5. 在GitHub上设置clean-main为默认分支
# Settings -> Branches -> Default branch -> clean-main

# 6. 可选：重命名旧分支
git push origin main:main-with-history
git branch -D main
git branch -m clean-main main
git push origin main --force
```

**优点**:
- 新用户克隆的是干净版本（约15MB）
- 历史仍然保留在main-with-history分支
- 可以随时切换回历史分支查看

**缺点**:
- 需要force push（团队需要重新克隆）
- 失去了commit历史（但可以在旧分支查看）

---

### 方案3：激进方案 - 使用git filter-repo清理历史 ⚠️ 慎用

**适用场景**: 想要彻底清理历史，减小仓库大小

**前提条件**:
```bash
# 安装git-filter-repo
pip install git-filter-repo
```

**操作**:
```bash
# 1. 备份仓库
cd /home/xumy/work/development/
cp -r Immunex Immunex_backup

cd Immunex

# 2. 删除大文件和legacy目录
git filter-repo --path immunex/legacy --invert-paths
git filter-repo --path development/angle_analysis/Analysis-scripts/Angles/脚本说明.pdf --invert-paths
git filter-repo --path development/angle_analysis/Analysis-scripts/Angles/脚本说明.docx --invert-paths
git filter-repo --path docs/截屏2026-01-07\ 13.33.39.png --invert-paths

# 3. 强制推送到GitHub
git remote add origin <your-repo-url>
git push origin --force --all
git push origin --force --tags

# 4. 通知团队成员重新克隆
```

**优点**:
- 彻底减小仓库大小（可能减少到几百MB）
- 历史更干净
- 克隆速度更快

**缺点**:
- ⚠️ 破坏性操作，不可逆
- 所有团队成员必须重新克隆
- 可能破坏已有的PR和Issues引用
- 需要协调团队

---

### 方案4：最佳实践 - 组合方案 ⭐⭐ 最推荐

**结合方案2和部分方案3**:

```bash
# 1. 创建干净的发布分支（方案2）
git checkout --orphan release-v1.0

# 2. 只添加核心代码
git add immunex/
git add scripts/
git add docs/
git add examples/
git add *.md
git add *.yml
git add *.yaml
git add setup.py
git add .gitignore

# 3. 排除不必要的文件
git reset development/  # 开发目录不包含在发布版
git reset test/         # 测试数据不包含
git reset input/        # 输入数据不包含
git reset output/       # 输出数据不包含

# 4. 创建发布提交
git commit -m "Release: Immunex v1.0.0

Core Features:
- TCR-pMHC MD trajectory analysis platform
- 15 CLI commands for comprehensive analysis
- Interactive HTML reports with AI assistant
- Batch processing and quality control
- Complete documentation and examples

Repository Size: ~15MB (clean release)
Full history available in 'main' branch
"

# 5. 推送发布分支
git push origin release-v1.0

# 6. 在GitHub上创建Release
# - Tag: v1.0.0
# - Branch: release-v1.0
# - 添加Release Notes
```

**优点**:
- 发布版本非常干净（15MB）
- 保留完整开发历史在main分支
- 用户可以选择克隆release分支（快速）或main分支（完整）
- 符合开源项目最佳实践

**缺点**:
- 需要维护两个分支
- 需要定期同步

---

## 推荐执行步骤

### 立即执行（今天）

1. ✅ **提交当前清理工作**
   ```bash
   git add .
   git commit -m "refactor: Phase 1 cleanup - remove legacy code and temporary files

   - Remove immunex/legacy/ directory
   - Clean up 2189 __pycache__ directories
   - Delete temporary files (.pyc, .pyo, .tmp, .bak)
   - Clean development/workspaces/ (24GB)
   - Update .gitignore for future commits
   - Unify CLI commands (remove redundant reporter.py)
   "
   git push origin main
   ```

2. ✅ **创建干净的发布分支**（方案4）
   ```bash
   # 执行上面"方案4"的步骤
   ```

### 后续执行（本周内）

3. **在GitHub上设置**
   - 创建v1.0.0 Release，指向release-v1.0分支
   - 更新README.md，说明两个分支的用途
   - 添加CONTRIBUTING.md，说明开发者应该使用main分支

4. **更新文档**
   - 在README中添加：
     ```markdown
     ## Installation

     For users (recommended):
     ```bash
     git clone -b release-v1.0 https://github.com/your-org/Immunex.git
     ```

     For developers:
     ```bash
     git clone https://github.com/your-org/Immunex.git
     ```
     ```

---

## 长期维护策略

### 分支策略

- **main**: 开发分支，包含完整历史和开发文件
- **release-v1.x**: 发布分支，只包含核心代码
- **feature/***: 功能分支，从main创建

### 发布流程

1. 在main分支开发和测试
2. 准备发布时，合并到release分支：
   ```bash
   git checkout release-v1.0
   git checkout main -- immunex/ scripts/ docs/ examples/
   git commit -m "Release: v1.1.0"
   git tag v1.1.0
   git push origin release-v1.0 --tags
   ```

### 文件大小控制

- 参考数据文件（>1MB）使用Git LFS
- 图片文件压缩后再提交
- 文档使用Markdown而非PDF/DOCX
- 测试数据不提交到仓库

---

## 总结

**当前状态**: Git仓库6.6GB，包含大量历史legacy代码

**推荐方案**: 方案4（组合方案）
- 创建干净的release-v1.0分支（15MB）
- 保留main分支的完整历史
- 用户克隆release分支，开发者使用main分支

**预期效果**:
- 用户克隆时间：从几分钟减少到几秒
- 仓库大小：从6.6GB减少到15MB（release分支）
- 保留完整历史供开发者参考

**风险**: 低（不破坏现有历史，只是创建新分支）

**执行时间**: 约30分钟
