# 差异分析一键软件

一个面向 Windows 用户的本地差异分析桌面软件。双击 `差异分析软件.exe` 后，软件会启动本地 Shiny 分析界面，支持导入表达矩阵和分组表，一键完成两组差异分析，并导出差异表、显著差异基因、火山图、热图和 PCA。

## 下载

推荐下载 Release 里的压缩包：

- [DifferentialAnalysisSoftware-v1.0.0.zip](https://github.com/HaPiJiucHi/differential-analysis-software/releases/download/v1.0.0/DifferentialAnalysisSoftware-v1.0.0.zip)

正式下载文件会放在 GitHub Release 附件中：

- `DifferentialAnalysisSoftware-v1.0.0.zip`

## 快速开始

1. 下载并解压 `差异分析软件.zip`。
2. 双击 `差异分析软件.exe`。
3. 如果提示缺少依赖，点击启动器窗口里的“检查依赖”。
4. 在浏览器界面导入表达矩阵和分组表。
5. 选择对照组、处理组和分析参数，点击“一键开始差异分析”。

## 主要功能

- 支持 `.csv`、`.tsv`、`.txt`、`.xlsx`、`.xls` 表达矩阵导入。
- 支持分组文件、样本名关键词、手动粘贴三种分组方式。
- 使用 `limma` 进行两组差异分析。
- 输出完整差异表和显著差异基因表。
- 火山图支持颜色调整和 Top 显著基因名称标注。
- 热图支持行聚类、列聚类开关和颜色调整。
- PCA 图用于检查样本分组趋势。
- 如果电脑没有 R，可通过“检查依赖”自动下载 R 并安装到软件同目录。

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

## 仓库内容

```text
.
├── app.R                         # Shiny 主程序
├── DifferentialAppLauncher.cs    # Windows exe 启动器源码
├── install_dependencies.R        # R 依赖安装脚本
├── test_app.R                    # 示例数据自测脚本
├── dist/                         # 本地打包输出目录；zip 通过 Release 发布
├── docs/                         # 使用文档
├── scripts/publish_to_github.ps1  # 一键发布脚本
└── .github/ISSUE_TEMPLATE/       # GitHub issue 模板
```

## 文档

- [安装说明](docs/INSTALL.md)
- [使用说明](docs/USAGE.md)
- [常见问题](docs/FAQ.md)
- [版本发布说明](RELEASE_NOTES.md)

## 本地开发

需要 Windows 和 R。安装依赖：

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

## 发布到 GitHub

本目录已经是一个标准 Git 仓库。登录 GitHub 后，可运行：

```powershell
.\scripts\publish_to_github.ps1
```

脚本会创建公开仓库 `differential-analysis-software`，推送当前代码，并创建 `v1.0.0` Release，把 `dist/差异分析软件.zip` 作为 Release 附件。

## 许可证

本项目使用 MIT License。见 [LICENSE](LICENSE)。
