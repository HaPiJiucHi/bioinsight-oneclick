# 差异分析一键软件

一个面向 Windows 用户的本地差异分析桌面软件。双击 `差异分析软件.exe` 后，软件会启动本地 Shiny 分析界面，支持导入表达矩阵和分组表，一键完成两组差异分析，并导出差异表、显著差异基因、火山图、热图、PCA、WGCNA、PPI 和 GSEA。

## 下载

推荐下载 GitHub Release 里的压缩包：

- [DifferentialAnalysisSoftware-v1.1.0.zip](https://github.com/HaPiJiucHi/differential-analysis-software/releases/download/v1.1.0/DifferentialAnalysisSoftware-v1.1.0.zip)

## 快速开始

1. 下载并解压 `DifferentialAnalysisSoftware-v1.1.0.zip`。
2. 双击 `差异分析软件.exe`。
3. 如果提示缺少依赖，点击启动器窗口里的“检查依赖”。
4. 在浏览器界面导入表达矩阵和分组表。
5. 选择对照组、处理组和分析参数，点击“一键开始差异分析”。

如果电脑没有 R，启动器可以把 R 安装到软件同目录下的 `R` 文件夹，尽量减少对系统环境的依赖。

## 主要功能

- 支持 `.csv`、`.tsv`、`.txt`、`.xlsx`、`.xls` 表达矩阵导入。
- 支持分组文件、样本名关键词、手动粘贴三种分组方式。
- 使用 `limma` 进行两组差异分析。
- 输出完整差异表和显著差异基因表。
- 火山图支持颜色调整和 Top 显著基因名称标注。
- 热图支持行聚类、列聚类开关和颜色调整。
- PCA 图用于检查样本分组趋势。
- WGCNA 用于识别与分组相关的共表达模块，并输出模块 hub genes。
- PPI 用于基于 STRING 互作表构建差异基因互作网络，并输出 hub genes。
- GSEA 用于基于全部基因排序解释 GO 通路整体偏向。

## 示例数据建议

内置示例数据是 20 个样本，Normal 10 个、Disease 10 个。推荐分析主线是：

```text
DEG 差异分析 + GSEA 通路解释为主，WGCNA 做模块层面辅助验证，PPI 做候选 hub gene 筛选。
```

原因：

- 示例数据显著差异基因数量充足，验证结果中上调 1624 个、下调 1179 个。
- GSEA 不依赖硬阈值，能解释整套基因排序的通路偏向，适合放在 DEG 之后解释机制。
- 示例 GSEA 靠前结果集中在有丝分裂姐妹染色单体分离、染色体分离、核染色体分离等细胞周期/染色体分离相关过程，NES 为负，提示这些过程整体偏向 Normal 组一侧。
- WGCNA 显示 turquoise 模块偏向 Disease、blue 模块偏向 Normal，说明样本存在清晰共表达模块结构；但 WGCNA 对样本量敏感，应作为辅助证据。
- PPI 在高置信阈值下给出 TOP2A、PRC1、EZH2 等 hub 候选，适合后续实验验证优先级排序。

## 输入格式

表达矩阵：第一列是基因或探针 ID，其余列是样本表达值。

```csv
feature_id,Sample_1,Sample_2,Sample_3,Sample_4
GeneA,8.1,8.3,10.2,10.0
GeneB,5.0,5.2,4.8,4.7
```

分组表：至少包含样本名和分组名两列。

```csv
sample,group
Sample_1,Control
Sample_2,Control
Sample_3,Treatment
Sample_4,Treatment
```

注释表可选：第一列为表达矩阵里的 ID，第二列为基因 symbol。

```csv
probe_id,symbol
1007_s_at,DDR1
1053_at,RFC2
```

## PPI 输入

软件默认会读取同目录下的 `string_interactions.tsv`。如果需要换成自己的物种或自己的基因集，可以从 STRING 导出 interaction 文件后上传。建议至少包含 `node1`、`node2`、`combined_score` 三列。

## 仓库内容

```text
.
├── app.R                         # Shiny 主程序
├── DifferentialAppLauncher.cs    # Windows exe 启动器源码
├── install_dependencies.R        # R 依赖安装脚本
├── string_interactions.tsv       # 示例 PPI 互作表
├── test_app.R                    # 示例数据基础自测脚本
├── dist/                         # 本地打包输出目录，zip 通过 Release 发布
├── docs/                         # 使用文档
├── scripts/publish_to_github.ps1 # 发布脚本
└── .github/ISSUE_TEMPLATE/       # GitHub issue 模板
```

## 文档

- [安装说明](docs/INSTALL.md)
- [使用说明](docs/USAGE.md)
- [常见问题](docs/FAQ.md)
- [版本发布说明](RELEASE_NOTES.md)

## 本地开发

安装依赖：

```powershell
& "C:\Program Files\R\R-4.5.3\bin\Rscript.exe" .\install_dependencies.R
```

运行自测：

```powershell
& "C:\Program Files\R\R-4.5.3\bin\Rscript.exe" .\test_app.R
```

直接启动 Shiny：

```powershell
& "C:\Program Files\R\R-4.5.3\bin\Rscript.exe" -e "shiny::runApp('.', launch.browser = TRUE, host = '127.0.0.1', port = 3838)"
```

## 许可证

本项目使用 MIT License。见 [LICENSE](LICENSE)。
