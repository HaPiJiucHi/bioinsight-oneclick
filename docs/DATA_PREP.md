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

## 2. 分组表是什么？

分组表告诉软件每个样本属于哪一组。最少两列：`sample` 和 `group`。

```csv
sample,group
Normal_1,Normal
Normal_2,Normal
Disease_1,Disease
Disease_2,Disease
```

注意：`sample` 必须和表达矩阵的列名完全一致。

## 3. 如果只有 GEO 数据怎么办？

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

## 4. 找不到分组列怎么办？

GEO 样本信息常见字段包括：

- `Sample_title`
- `Sample_source_name_ch1`
- `Sample_characteristics_ch1`
- `disease`
- `treatment`
- `group`
- `tissue`

如果没有一列刚好能分组，可以用软件左侧“手动粘贴”方式，把样本按自己的研究设计填成 `sample,group`。

## 5. 推荐给新手的顺序

1. 先点击“载入当前示例数据”跑通一次。
2. 下载表达矩阵模板和分组表模板。
3. 如果数据来自 GEO，先试 `series_matrix.txt.gz` 自动拆分。
4. 在“数据检查”页确认表达矩阵预览和分组预览都正常。
5. 再点击“一键开始分析”。
