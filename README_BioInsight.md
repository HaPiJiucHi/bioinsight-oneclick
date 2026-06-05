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

如果提示缺少 R 包，先运行：

```powershell
& "C:\Program Files\R\R-4.5.3\bin\Rscript.exe" .\install_dependencies.R
```

说明：启动器会优先调用软件目录下的 `R` 文件夹；如果没有，再查找系统已安装的 R。点击“检查依赖”时，如果电脑没有 R，软件会从 CRAN 下载 Windows R 安装包，并静默安装到当前软件目录的 `R` 文件夹，然后继续安装 R 包依赖。

## 输入格式

表达矩阵支持 `.csv`、`.tsv`、`.txt`、`.xlsx`、`.xls`。第一列是基因或探针 ID，其余列是样本表达值：

```csv
feature_id,Sample_1,Sample_2,Sample_3,Sample_4
GeneA,8.1,8.3,10.2,10.0
GeneB,5.0,5.2,4.8,4.7
```

分组表至少两列：样本名、分组名。样本名必须与表达矩阵列名一致：

```csv
sample,group
Sample_1,Control
Sample_2,Control
Sample_3,Treatment
Sample_4,Treatment
```

注释表可选：第一列为表达矩阵里的 ID，第二列为基因 symbol。

## 主要功能

1. 导入表达矩阵和分组表，或直接载入当前示例数据。
2. 使用 `limma` 完成两组 DEG 分析。
3. 输出完整差异表、显著差异基因表、火山图、热图和 PCA。
4. 火山图支持标注前 N 个显著基因，并可调整上调、下调、稳定基因颜色。
5. 热图支持行聚类、列聚类开关，并可调整低值、中间值、高值颜色。
6. PCA 支持显示分组椭圆和分组中心点。
7. GSEA 用全部基因排序解释通路整体偏向，支持 GO BP、MF、CC。
8. WGCNA 用于寻找与分组相关的共表达模块，并输出模块相关性和模块 hub genes。
9. PPI 用于基于 STRING 互作表构建差异基因蛋白互作网络，并输出 hub genes。

## 当前示例样本建议

这个示例数据是 20 个样本，Normal 10 个、Disease 10 个。更适合采用：

```text
DEG 差异分析 + GSEA 通路解释为主，WGCNA 做模块层面辅助验证，PPI 做候选 hub gene 筛选。
```

原因：DEG 先明确哪些基因显著变化；GSEA 不依赖单个显著阈值，能解释通路整体轻中度同向变化；WGCNA 在 20 个样本下可作为探索性模块证据；PPI 适合从差异基因里筛选后续验证的 hub gene。

## PPI 输入说明

软件默认会读取同目录下的 `string_interactions.tsv`。如果需要换成自己的物种或自己的基因集，可以从 STRING 导出 interaction 文件后上传。建议至少包含 `node1`、`node2`、`combined_score` 三列。

## 自测

基础自测：

```powershell
& "C:\Program Files\R\R-4.5.3\bin\Rscript.exe" .\test_app.R
```

当前新增功能已用示例数据验证通过：DEG、PPI、GSEA、WGCNA 均可返回结果。


