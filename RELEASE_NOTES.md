# Release Notes

## v1.5.2

### 富集分析增强

- 下调基因、合并分析、GO/KEGG 的气泡图改为显著条目优先；如果显著条目不足，会自动补充排名靠前的趋势条目，避免图上只剩几条。
- 新增 `GO All（BP+MF+CC）` 选项，一次运行 GO 三大类，结果表用 `collection` 区分 BP、MF 和 CC。

### PPI 增强

- PPI 自选基因改为从当前表达矩阵/差异结果中搜索选择，不再要求用户随便手打。
- 支持按 gene symbol 或表达矩阵 feature ID 匹配 PPI 节点；ENSEMBL/ENTREZ ID 会尽量用本地 `org.Hs.eg.db` 转成 SYMBOL。
- 支持公司同款 STRING link 表，`combined_score` 可为 0-1 或 0-1000。
- 新增公司风格 Top 6 community 模块图。

### 箱线图与 GSEA

- 箱线图新增“公司 FPKM 风格：log10(value + 0.001)”尺度。导入 FPKM/TPM 矩阵时可复现公司图里负值和样本分布更齐的效果。
- GSEA running enrichment 图标题新增 NES。
- 曲线图下方新增对应通路的 core enrichment genes、leading-edge genes、pvalue、padj 和 NES 表。
- 结果表中的 P 值、padj、qvalue、FDR 默认使用固定小数显示，避免 `1.29e-08` 这种科学计数法。

## v1.5.1

### 结果一致性修复

- 修复 featureCounts 注释表自动识别，上传带 `gene_name` 的注释文件后，火山图、热图、富集分析和 GSEA 会优先使用真实基因名。
- RNA-seq raw counts 默认过滤总 counts 小于 2 的极低表达基因，默认 DESeq2 结果口径更接近常见公司报告流程。
- 修复 GO/KEGG 富集和 GSEA 因 SYMBOL 识别失败而报错的问题。

### GSEA 与 PPI 增强

- GSEA 结果表新增 `core_enrichment_genes`、`leading_edge_genes`、核心基因数量等字段。
- GSEA 新增经典 running enrichment 曲线图，便于查看通路富集峰值和命中基因分布。
- PPI 支持手动输入候选基因，也支持在线 STRING 查询；节点颜色表示上调、下调或不显著，节点大小反映连接数和差异幅度。
- 箱线图默认显示离群点，便于发现特殊样本。

## v1.5.0

### 流程和解释增强

- 数据类型自动识别、log2 转换和 raw count 方法选择移动到“数据检查”页。
- 每个页签顶部新增通俗解释，说明这一页的目的和能看出什么。
- 页面顺序调整为更接近实际分析流程：数据检查、差异结果、火山图、热图、箱线图、PCA、富集分析、GSEA、WGCNA、PPI。

### 新增箱线图

- 新增“箱线图”页，用于查看每个样本整体表达分布。
- 可用于识别明显偏高、偏低或分布异常的特殊样本。
- 支持设置抽样基因数和是否显示离群点。

### 新增差异基因富集分析

- 在 PCA 后、GSEA 前新增“富集分析”页。
- 支持上调基因、下调基因、上调+下调合并、上调/下调分别分析。
- 支持 GO BP、GO MF、GO CC 和 KEGG。
- 输出富集结果表、气泡图和白话解释。

## v1.4.0

### 更智能的一键分析

- 数据类型默认改为“自动识别并选择方法”。
- 软件会根据表达矩阵数值特征自动判断：
  - 非负整数且数值较大：按 RNA-seq raw counts 处理，默认 `DESeq2`。
  - 小数且范围较大：按 RNA-seq TPM/FPKM/RPKM 处理，自动 log2 判断后走 `limma`。
  - 已 log2 的小数或含负值：按芯片/已标准化表达矩阵走 `limma`。
- 用户仍可手动覆盖数据类型和 raw count 分析方法。

### GSEA 与火山图增强

- GSEA 新增 KEGG，适合解释经典信号通路和代谢通路。
- GO BP/MF/CC 与 KEGG 共用结果表、气泡图和自动解释。
- 火山图新增“标注升高基因”和“标注降低基因”独立开关。
- 升高和降低基因可分别设置 Top N 标注数量。

## v1.3.0

### 数据类型与 RNA-seq counts 支持

- “分析结果”页新增数据类型选择：
  - 芯片 / 已标准化表达矩阵：使用 `limma`。
  - RNA-seq TPM / FPKM / RPKM：自动或手动 log2 后使用 `limma`。
  - RNA-seq raw counts：支持 `DESeq2`、`edgeR`、`limma-voom`。
- raw count 输入会检查是否为非负整数，避免把 TPM/FPKM 误当 counts。
- 结果表统一输出 `logFC`、`P.Value`、`adj.P.Val` 和 `change`，后续火山图、热图、PCA、GSEA、WGCNA、PPI 继续复用。
- 页面新增“傻瓜式”数据类型科普，并根据矩阵数值给出粗略判断。
- 依赖安装脚本新增 `DESeq2` 和 `edgeR`。

## v1.2.1

### 新手导入增强

- 左侧导入区新增表达矩阵模板和分组表模板下载。
- 新增 GEO `series_matrix.txt.gz` 上传入口，可拆出表达矩阵和样本信息表。
- “数据检查”页新增 GEO 样本信息预览。
- 可从 GEO 样本信息中选择一列，一键生成分组表。
- 新增 [新手数据准备](docs/DATA_PREP.md) 文档。
- 删除公开仓库中的个人视频发布脚本和相关链接。

## v1.2.0

### 品牌与发布包更新

- 软件更名为 **BioInsight 一键生信分析平台**，不再只强调差异分析，突出 DEG、GSEA、WGCNA、PPI 的完整工作流。
- 新增原创蜡笔手绘风应用图标，并嵌入 Windows `.exe`。
- Windows 启动器文件名更新为 `BioInsight 一键生信分析平台.exe`。
- 发布压缩包更新为 `BioInsight-OneClick-Bioinformatics-v1.2.0.zip`。
- GitHub 首页和安装说明同步使用新名称。

## v1.1.2

### GitHub 首页更新

- 将 README 改为教程型首页，突出软件解决的问题、安装方式、分析流程和结果展示。
- 增加功能截图、结果展示图和 Mermaid 流程图。
- 增加 `scripts/create_readme_assets.R`，可重新生成 README 展示图。

## v1.1.1

### 界面更新

- 左侧栏只保留数据导入、注释表和分组信息。
- 差异分析、火山图、热图、PCA、GSEA、WGCNA、PPI 参数移动到对应页签顶部。
- GSEA 页签移动到 WGCNA 前面。
- PCA 增加“显示分组椭圆”和“显示分组中心点”选项。
- GSEA 页签下方增加说明，解释为什么在差异分析之后使用 GSEA。

## v1.1.0

### 新增功能

- 增加 WGCNA 共表达模块分析。
- 增加 PPI 网络分析，默认支持同目录 `string_interactions.tsv`。
- 增加 GSEA GO 富集分析，支持 BP、MF、CC。
- 在结果页增加“最终建议”，说明当前样本更适合采用哪条分析主线。
- 更新依赖安装脚本，补充 `WGCNA`、`igraph`、`clusterProfiler`、`org.Hs.eg.db`、`enrichplot`。

## v1.0.0

### 功能

- Windows `.exe` 桌面启动器。
- 本地 Shiny 分析界面。
- 表达矩阵和分组表导入。
- `limma` 两组差异分析。
- 火山图、热图、PCA、GSEA、WGCNA、PPI。
