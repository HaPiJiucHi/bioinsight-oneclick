# BioInsight 一键生信分析平台

BioInsight 是一个本地 Windows 桌面式生信分析平台。双击 `BioInsight 一键生信分析平台.exe` 后，软件会启动本地 Shiny 图形界面，用于导入表达矩阵和分组信息，一键完成 DEG、火山图、热图、PCA、GSEA、WGCNA 和 PPI 分析。

## 启动

优先双击：

```text
BioInsight 一键生信分析平台.exe
```

也可以双击批处理启动器：

```text
启动BioInsight.bat
start_bioinsight_app.bat
```

如果提示缺少 R 包，先点击启动器里的“检查依赖”。如果电脑没有 R，启动器会尝试把 R 安装到软件同目录下的 `R` 文件夹。

## 新手怎么准备数据

最少需要两个文件：

1. 表达矩阵：第一列是基因或探针 ID，其余列是样本表达值。
2. 分组表：至少两列，`sample` 和 `group`。

软件左侧“新手：没有表达矩阵/分组表？”里提供：

- 表达矩阵模板下载。
- 分组表模板下载。
- GEO `series_matrix.txt.gz` 上传和拆分。
- 拆出的 GEO 表达矩阵下载。
- 拆出的 GEO 样本信息表下载。

如果数据来自 GEO，推荐下载 `series_matrix.txt.gz`，上传后点击“拆出 GEO 表达矩阵”。在“数据检查”页选择样本信息中的疾病、处理、组织、分型等字段，就可以一键生成分组表。

## 输入格式

表达矩阵示例：

```csv
feature_id,Normal_1,Normal_2,Disease_1,Disease_2
GeneA,8.1,8.3,10.2,10.0
GeneB,5.0,5.2,4.8,4.7
```

分组表示例：

```csv
sample,group
Normal_1,Normal
Normal_2,Normal
Disease_1,Disease
Disease_2,Disease
```

注释表可选：第一列为表达矩阵里的 ID，第二列为基因 symbol。

## 主要功能

1. 导入表达矩阵和分组表，或直接载入当前示例数据。
2. 从 GEO `series_matrix.txt.gz` 拆出表达矩阵和样本信息。
3. 使用 `limma` 完成两组 DEG 分析。
4. 输出完整差异表、显著差异基因表、火山图、热图和 PCA。
5. GSEA 用全部基因排序解释通路整体偏向，支持 GO BP、MF、CC。
6. WGCNA 用于寻找与分组相关的共表达模块，并输出模块相关性和模块 hub genes。
7. PPI 用于基于 STRING 互作表构建差异基因蛋白互作网络，并输出 hub genes。

## 自测

基础自测：

```powershell
& "C:\Program Files\R\R-4.5.3\bin\Rscript.exe" .\test_app.R
```

当前功能已用示例数据验证通过；根目录存在 GEO 示例文件时，也会验证 GEO series matrix 拆分。
