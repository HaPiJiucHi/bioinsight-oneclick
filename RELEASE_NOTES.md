# Release Notes

## v1.1.2

### GitHub 首页更新

- 将 README 改为教程型首页，突出软件解决的问题、安装方式、分析流程和结果展示。
- 增加功能截图、结果展示图和 Mermaid 流程图。
- 增加 `scripts/create_readme_assets.R`，可重新生成 README 展示图。
- 增加 `docs/DOUYIN_VIDEO_SCRIPT.md`，包含 60 秒抖音讲解脚本、标题、封面文案和评论区引导。

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

### 示例数据验证

- 差异分析：上调 1624 个，下调 1179 个。
- GSEA：靠前结果主要集中在细胞周期、姐妹染色单体分离、染色体分离相关过程。
- WGCNA：turquoise 模块偏向 Disease，blue 模块偏向 Normal。
- PPI：高置信默认参数下 TOP2A、PRC1、EZH2 等位于 hub 前列。

## v1.0.0

首个公开发布版本。

### 功能

- Windows `.exe` 桌面启动器。
- 本地 Shiny 差异分析界面。
- 表达矩阵和分组表导入。
- `limma` 两组差异分析。
- 完整差异表与显著差异基因表下载。
- 火山图、热图、PCA 图输出。
- 火山图显著基因标注。
- 火山图和热图颜色调整。
- 热图行聚类和列聚类开关。
- 缺少 R 时可通过“检查依赖”下载并安装到软件同目录。
