# 新手数据准备

BioInsight 需要两个核心输入：表达矩阵和分组信息。注释表是可选的。

## 1. 表达矩阵是什么？

表达矩阵就是“行是基因或探针，列是样本”的表格。第一列放基因 ID 或探针 ID，后面每一列是一个样本的表达值。

```csv
feature_id,Normal_1,Normal_2,Disease_1,Disease_2
GeneA,8.1,8.3,10.2,10.0
GeneB,5.0,5.2,4.8,4.7
```

软件左侧可以下载“表达矩阵模板”。

## 2. 我怎么知道自己是芯片、TPM/FPKM，还是 raw counts？

很多新手真正卡住的是“不知道公司或 GEO 给的文件是什么”。可以按下面判断：

| 你看到的数据特征 | 大概率是什么 | 软件里怎么选 |
|---|---|---|
| 值基本都是整数，例如 `0, 1, 2, 50, 1000` | RNA-seq raw counts | `RNA-seq raw counts`，默认 `DESeq2`，并默认过滤总 counts 小于 2 的极低表达基因 |
| 值有很多小数，文件名有 `TPM`、`FPKM`、`RPKM` | RNA-seq TPM/FPKM/RPKM | `RNA-seq TPM / FPKM / RPKM` |
| 值多在 `0-20` 左右，常见小数，GEO 芯片矩阵或 normalized expression | 芯片/已标准化表达矩阵 | `芯片 / 已标准化表达矩阵` |
| 数据里有负值 | 通常是已标准化或中心化后的表达矩阵 | `芯片 / 已标准化表达矩阵` |
| 文件列名已经是 `logFC`、`P.Value`、`padj` | 差异分析结果表，不是表达矩阵 | 不能直接当表达矩阵导入 |

拿不准时，可以在“数据检查”页保持软件默认的“自动识别并选择方法”。软件会先看数值是不是非负整数、是否有小数、数值范围和是否有负值，然后自动选择 counts、TPM/FPKM 或芯片/标准化矩阵路线。判断错了再手动改。

仍然建议优先看文件名和说明书。`counts` 做正式 RNA-seq 差异分析更合适；`TPM/FPKM` 更适合快速作图、聚类、GSEA 和探索性分析。很多测序公司在 DESeq2 前会过滤总 counts 很低的基因，BioInsight 默认按“基因总 counts 至少 2”处理；如果你的公司说明书用了其他阈值，可以在“数据检查”页调整。

## 3. 分组表是什么？

分组表告诉软件每个样本属于哪一组。最少两列：`sample` 和 `group`。

```csv
sample,group
Normal_1,Normal
Normal_2,Normal
Disease_1,Disease
Disease_2,Disease
```

注意：`sample` 必须和表达矩阵的列名完全一致。

## 4. 如果只有 GEO 数据怎么办？

很多新手拿到的是 GEO 页面，而不是现成的表达矩阵和分组表。推荐步骤：

1. 打开 GEO 数据集页面。
2. 找到 `Download family` 或 `Series Matrix File(s)`。
3. 下载 `series_matrix.txt.gz`。
4. 在 BioInsight 左侧“新手：没有表达矩阵/分组表？”里上传这个文件。
5. 点击“拆出 GEO 表达矩阵”。
6. 在“数据检查”页查看 GEO 样本信息。
7. 如果样本信息里有疾病、处理、组织、分型等字段，选择该列，点击“用该列生成分组表”。

拆分后也可以下载：

- `GEO_expression_matrix_日期.csv`
- `GEO_sample_metadata_日期.csv`

## 5. 找不到分组列怎么办？

GEO 样本信息常见字段包括：

- `Sample_title`
- `Sample_source_name_ch1`
- `Sample_characteristics_ch1`
- `disease`
- `treatment`
- `group`
- `tissue`

如果没有一列刚好能分组，可以用软件左侧“手动粘贴”方式，把样本按自己的研究设计填成 `sample,group`。

## 6. 推荐给新手的顺序

1. 先点击“载入当前示例数据”跑通一次。
2. 下载表达矩阵模板和分组表模板。
3. 如果数据来自 GEO，先试 `series_matrix.txt.gz` 自动拆分。
4. 在“数据检查”页确认表达矩阵预览和分组预览都正常。
5. 在“分析结果”页保留默认“自动识别并选择方法”；如果软件判断不符合你的文件说明，再手动选择芯片/TPM-FPKM/raw counts。
6. 再点击“一键开始分析”。
