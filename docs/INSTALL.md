# 安装说明

## 普通用户

1. 下载 `DifferentialAnalysisSoftware-v1.1.1.zip`。
2. 解压到任意英文或中文目录。
3. 双击 `差异分析软件.exe`。
4. 第一次打开时，如果提示缺少 R 或 R 包，点击“检查依赖”。

“检查依赖”会执行两件事：

- 如果电脑没有 R，会下载 CRAN 的 Windows R 安装包，并安装到软件同目录的 `R` 文件夹。
- 安装或检查 Shiny、limma、ggplot2、pheatmap、WGCNA、igraph、clusterProfiler、org.Hs.eg.db、enrichplot 等分析依赖。

## 已安装 R 的用户

软件会按以下顺序查找 R：

1. 软件目录下的 `R` 文件夹。
2. 系统安装目录，例如 `C:\Program Files\R\R-4.5.3`。
3. PATH 中的 `Rscript.exe`。

## 手动安装依赖

在软件目录打开 PowerShell：

```powershell
& "C:\Program Files\R\R-4.5.3\bin\Rscript.exe" .\install_dependencies.R
```

## 网络要求

只有在缺少 R 或 R 包时才需要联网。R 会从 CRAN 下载，Bioconductor 包会从配置的 Bioconductor 镜像下载。安装完成后，基础差异分析、PPI 默认互作表和 GSEA 注释包都可以在本地运行。
