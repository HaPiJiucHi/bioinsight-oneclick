# Changelog

## v1.5.2

- 富集分析绘图改为“显著项优先，不足显示条数时自动补排名靠前的趋势项”，避免下调或合并分析只显示少数几条。
- 富集分析新增 `GO All（BP+MF+CC）`，可一次合并查看 GO 三大类，并在结果表中保留 `collection` 区分来源。
- PPI 自选基因改为从当前表达矩阵/差异表中搜索选择，支持按基因名或表达矩阵 ID 检索。
- PPI 本地/上传 STRING 表支持 node 使用基因名或 feature ID，并新增公司风格 Top 6 模块网络图模式。
- 箱线图新增显示尺度：分析矩阵、原始输入矩阵、公司 FPKM 风格 `log10(value + 0.001)`、按样本中位数中心化。
- GSEA running enrichment 图标题新增 NES；图下新增对应通路的 `core_enrichment_genes`、`leading_edge_genes`、pvalue 和 padj 表。
- 所有结果表中的 P 值、padj、qvalue、FDR 等默认用固定小数显示，避免科学计数法影响阅读。
- 修复本地路径读取 PPI/表格文件时的路径对象兼容问题。

## v1.5.1

- 修复 featureCounts 注释表自动识别，优先使用 `gene_name`、`symbol`、`gene_symbol` 等列，避免火山图和热图显示数字 ID。
- RNA-seq raw counts 默认过滤总 counts 小于 2 的极低表达基因，使结果口径更接近常见公司 DESeq2 流程。
- 修复富集分析和 GSEA 因 SYMBOL 识别失败导致无法运行的问题。
- GSEA 结果表新增 `core_enrichment_genes`、`leading_edge_genes` 和基因数量统计，并新增 running enrichment 曲线图。
- PPI 支持显著差异基因或手动输入候选基因，支持在线 STRING 查询；节点颜色表示上调、下调或不显著，节点大小结合连接数和差异幅度。
- 箱线图默认显示离群点，便于检查特殊样本。

## v1.5.0

- 将数据类型自动识别和转换设置移动到“数据检查”页。
- 新增“箱线图”页，用于检查样本整体表达分布和特殊样本。
- 在 PCA 之后、GSEA 之前新增“富集分析”页，支持上调、下调、合并、上/下调分别分析。
- 富集分析支持 GO BP、GO MF、GO CC 和 KEGG。
- 每个分析页顶部新增白话解释，说明这一页为什么做、能得到什么。
- 页面顺序调整为：数据检查、分析结果、火山图、热图、箱线图、PCA、富集分析、GSEA、WGCNA、PPI。

## v1.4.0

- 数据类型默认改为“自动识别并选择方法”，运行时自动匹配 raw counts、TPM/FPKM/RPKM 或芯片/标准化表达矩阵路线。
- 自动模式下 RNA-seq raw counts 默认使用 `DESeq2`，TPM/FPKM/RPKM 和芯片/标准化矩阵自动使用合适的 `limma` 流程。
- GSEA 新增 KEGG，和 GO BP/MF/CC 共用同一套排序和解释界面。
- 火山图新增升高基因、降低基因独立标注开关和 Top N 设置。

## v1.3.0

- 增加数据类型选择：芯片/已标准化表达矩阵、RNA-seq TPM/FPKM/RPKM、RNA-seq raw counts。
- RNA-seq raw counts 新增 `DESeq2`、`edgeR`、`limma-voom` 三种差异分析方法。
- 结果摘要新增“数据类型”和“方法”，便于确认本次分析实际使用的模型。
- 分析结果页新增小白科普和矩阵粗略判断，帮助用户区分 counts、TPM/FPKM 和芯片/标准化矩阵。
- 依赖安装脚本新增 `DESeq2` 和 `edgeR`。

## v1.2.1

- 删除公开仓库中的个人视频发布脚本和相关入口。
- 增加新手数据准备入口：表达矩阵模板、分组表模板。
- 增加 GEO `series_matrix.txt.gz` 拆分功能，可自动生成表达矩阵和样本信息表。
- 增加 GEO 样本信息列生成分组表的流程。
- 新增 `docs/DATA_PREP.md`，说明新手如何准备表达矩阵和分组信息。

## v1.2.0

- 软件正式更名为 **BioInsight 一键生信分析平台**，突出 DEG、GSEA、WGCNA、PPI 等完整生信分析能力。
- 新增原创蜡笔手绘风 BioInsight 图标，并嵌入 Windows 启动器 exe。
- 发行包重命名为 `BioInsight-OneClick-Bioinformatics-v1.2.0.zip`。
- GitHub 仓库更名为 `bioinsight-oneclick`，README 下载入口同步更新。

## v1.1.2

- 重做 GitHub 首页为教程型项目风格。
- 增加软件界面截图、火山图、热图、PCA、GSEA、WGCNA、PPI 结果展示图。
- 增加 README 展示图生成脚本。

## v1.1.1

- 调整界面布局：左侧只保留导入、注释和分组。
- 将差异分析、火山图、热图、PCA、GSEA、WGCNA、PPI 参数移动到对应页签顶部。
- 将 GSEA 页签移动到 WGCNA 前面。
- PCA 增加分组椭圆和分组中心点显示选项。
- 在 GSEA 页签下方增加“为什么使用 GSEA”的说明。

## v1.1.0

- 增加 WGCNA、PPI、GSEA 高级分析模块。
- 增加结果页最终建议和样本解释。
- 更新 PPI 默认参数，示例数据默认可生成高置信 PPI 网络。
- 更新依赖安装脚本和说明文档。

## v1.0.0

- 发布 Windows 桌面启动器和本地 Shiny 分析界面。
- 增加 R 自动检测和本地安装逻辑。
- 增加颜色调整、火山图标签、热图聚类选项。
- 增加 GitHub 发布文档和 issue 模板。
