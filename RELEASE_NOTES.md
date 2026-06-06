# Release Notes

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
