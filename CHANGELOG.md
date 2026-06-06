# Changelog

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
