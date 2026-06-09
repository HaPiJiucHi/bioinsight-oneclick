options(shiny.maxRequestSize = 500 * 1024^2)

required_packages <- c(
  "shiny", "limma", "ggplot2", "readr", "readxl", "DT", "pheatmap",
  "matrixStats", "ggrepel", "colourpicker", "WGCNA", "igraph",
  "clusterProfiler", "org.Hs.eg.db", "enrichplot", "DESeq2", "edgeR"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "缺少 R 包：", paste(missing_packages, collapse = ", "),
    "\n请先运行 install_dependencies.R 安装依赖。"
  )
}

library(shiny)
library(ggplot2)
suppressPackageStartupMessages(library(WGCNA))

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

trim_text <- function(x) {
  trimws(as.character(x))
}

numeric_clean <- function(x) {
  if (is.numeric(x)) {
    return(x)
  }
  x <- trim_text(x)
  x <- gsub(",", "", x, fixed = TRUE)
  x[x %in% c("", "NA", "NaN", "NULL", "null", "-", "--")] <- NA
  suppressWarnings(as.numeric(x))
}

count_delim <- function(lines, delim) {
  sum(vapply(strsplit(lines, delim, fixed = TRUE), length, integer(1)) - 1)
}

guess_delim <- function(path) {
  lines <- readLines(path, n = 5, warn = FALSE, encoding = "UTF-8")
  if (length(lines) == 0) {
    return("\t")
  }
  counts <- c(
    tab = count_delim(lines, "\t"),
    comma = count_delim(lines, ","),
    semicolon = count_delim(lines, ";")
  )
  c(tab = "\t", comma = ",", semicolon = ";")[[which.max(counts)]]
}

read_any_table <- function(file) {
  if (is.list(file)) {
    path <- file$datapath %||% file
    name <- file$name %||% basename(path)
  } else {
    path <- file
    name <- basename(path)
  }
  ext <- tolower(tools::file_ext(name))

  if (ext %in% c("xlsx", "xls")) {
    out <- readxl::read_excel(path, .name_repair = "unique")
    return(as.data.frame(out, check.names = FALSE))
  }

  delim <- switch(
    ext,
    csv = ",",
    tsv = "\t",
    tab = "\t",
    txt = guess_delim(path),
    guess_delim(path)
  )

  out <- readr::read_delim(
    path,
    delim = delim,
    show_col_types = FALSE,
    progress = FALSE,
    guess_max = 10000,
    trim_ws = TRUE
  )
  as.data.frame(out, check.names = FALSE)
}

expression_template <- function() {
  data.frame(
    feature_id = c("GeneA", "GeneB", "GeneC"),
    Normal_1 = c(8.1, 5.0, 10.4),
    Normal_2 = c(8.3, 5.2, 10.1),
    Disease_1 = c(10.2, 4.8, 12.8),
    Disease_2 = c(10.0, 4.7, 12.5),
    check.names = FALSE
  )
}

group_template <- function() {
  data.frame(
    sample = c("Normal_1", "Normal_2", "Disease_1", "Disease_2"),
    group = c("Normal", "Normal", "Disease", "Disease"),
    stringsAsFactors = FALSE
  )
}

read_geo_lines <- function(path) {
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    gzfile(path, open = "rt", encoding = "UTF-8")
  } else {
    file(path, open = "rt", encoding = "UTF-8")
  }
  on.exit(close(con), add = TRUE)
  readLines(con, warn = FALSE)
}

parse_geo_record <- function(line) {
  out <- read.delim(
    text = line,
    header = FALSE,
    sep = "\t",
    quote = "\"",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  as.character(out[1, ])
}

safe_column_name <- function(x) {
  x <- trim_text(x)
  x <- gsub("[^A-Za-z0-9_]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  ifelse(nzchar(x), x, "metadata")
}

parse_geo_series_matrix <- function(file) {
  path <- if (is.list(file) && !is.null(file$datapath)) file$datapath else file
  lines <- read_geo_lines(path)
  begin <- grep("^!series_matrix_table_begin", lines, ignore.case = TRUE)
  end <- grep("^!series_matrix_table_end", lines, ignore.case = TRUE)
  if (length(begin) == 0 || length(end) == 0 || end[1] <= begin[1] + 1) {
    stop("没有在 GEO Series Matrix 中找到表达矩阵区段。请确认上传的是 series_matrix.txt 或 series_matrix.txt.gz。")
  }

  matrix_text <- paste(lines[(begin[1] + 1):(end[1] - 1)], collapse = "\n")
  expr <- read.delim(
    text = matrix_text,
    header = TRUE,
    sep = "\t",
    quote = "\"",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (ncol(expr) < 3) {
    stop("GEO 表达矩阵样本列太少，无法继续分析。")
  }
  names(expr)[1] <- "feature_id"

  samples <- names(expr)[-1]
  sample_meta <- data.frame(sample = samples, stringsAsFactors = FALSE)
  sample_lines <- grep("^!Sample_", lines, value = TRUE)
  used_names <- names(sample_meta)

  for (line in sample_lines) {
    fields <- parse_geo_record(line)
    if (length(fields) < 2) {
      next
    }
    key <- sub("^!", "", fields[1])
    values <- fields[-1]
    if (length(values) != length(samples)) {
      next
    }
    base_name <- safe_column_name(key)
    col_name <- make.unique(c(used_names, base_name))[length(used_names) + 1]
    used_names <- c(used_names, col_name)
    sample_meta[[col_name]] <- trim_text(values)

    colon_pos <- regexpr(":", values, fixed = TRUE)
    if (all(colon_pos > 0)) {
      labels <- trim_text(substr(values, 1, colon_pos - 1))
      if (length(unique(labels)) == 1 && nzchar(labels[1])) {
        parsed_name <- safe_column_name(labels[1])
        parsed_col <- make.unique(c(used_names, parsed_name))[length(used_names) + 1]
        used_names <- c(used_names, parsed_col)
        sample_meta[[parsed_col]] <- trim_text(substr(values, colon_pos + 1, nchar(values)))
      }
    }
  }

  list(expr = expr, samples = sample_meta)
}

geo_group_candidate_cols <- function(sample_meta) {
  cols <- setdiff(names(sample_meta), "sample")
  if (length(cols) == 0) {
    return(character())
  }
  good <- vapply(cols, function(col) {
    values <- trim_text(sample_meta[[col]])
    values <- values[nzchar(values)]
    n <- length(unique(values))
    n >= 2 && n <= min(10, max(2, nrow(sample_meta) - 1))
  }, logical(1))
  out <- cols[good]
  if (length(out) > 0) out else cols
}

needs_log2_transform <- function(mat) {
  x <- mat[is.finite(mat)]
  if (length(x) == 0 || min(x, na.rm = TRUE) < 0) {
    return(FALSE)
  }
  max(x, na.rm = TRUE) > 50 || stats::quantile(x, 0.99, na.rm = TRUE) > 30
}

prepare_expression <- function(df, id_col, log_mode, normalize_between_arrays) {
  if (is.null(df) || nrow(df) == 0) {
    stop("表达矩阵为空。")
  }
  if (!id_col %in% names(df)) {
    stop("找不到选择的 ID 列。")
  }

  sample_cols <- setdiff(names(df), id_col)
  if (length(sample_cols) < 2) {
    stop("表达矩阵至少需要 2 个样本列。")
  }

  feature_id <- trim_text(df[[id_col]])
  valid_feature <- !is.na(feature_id) & nzchar(feature_id)
  if (!any(valid_feature)) {
    stop("ID 列没有可用的基因/探针 ID。")
  }

  mat_df <- df[valid_feature, sample_cols, drop = FALSE]
  mat_df[] <- lapply(mat_df, numeric_clean)
  mat <- as.matrix(mat_df)
  storage.mode(mat) <- "double"

  good_cols <- colSums(is.finite(mat)) > 0
  mat <- mat[, good_cols, drop = FALSE]
  sample_cols <- sample_cols[good_cols]
  if (ncol(mat) < 2) {
    stop("可识别为数值的样本列少于 2 列。")
  }

  good_rows <- rowSums(is.finite(mat)) >= 2
  mat <- mat[good_rows, , drop = FALSE]
  feature_id <- feature_id[valid_feature][good_rows]
  if (nrow(mat) < 2) {
    stop("可用于分析的基因/探针少于 2 行。")
  }

  rownames(mat) <- make.unique(feature_id)
  colnames(mat) <- sample_cols
  input_mat <- mat

  transformed <- FALSE
  if (log_mode == "always" || (log_mode == "auto" && needs_log2_transform(mat))) {
    if (min(mat, na.rm = TRUE) < 0) {
      stop("数据包含负值，不能执行 log2(x + 1) 转换。")
    }
    mat <- log2(mat + 1)
    transformed <- TRUE
  }

  if (isTRUE(normalize_between_arrays)) {
    mat <- limma::normalizeBetweenArrays(mat, method = "quantile")
  }

  list(
    mat = mat,
    input_mat = input_mat,
    feature_id = feature_id,
    transformed = transformed,
    normalized = isTRUE(normalize_between_arrays)
  )
}

prepare_count_expression <- function(df, id_col, min_total_count = 2) {
  if (is.null(df) || nrow(df) == 0) {
    stop("表达矩阵为空。")
  }
  if (!id_col %in% names(df)) {
    stop("找不到选择的 ID 列。")
  }

  sample_cols <- setdiff(names(df), id_col)
  if (length(sample_cols) < 2) {
    stop("raw count 矩阵至少需要 2 个样本列。")
  }

  feature_id <- trim_text(df[[id_col]])
  valid_feature <- !is.na(feature_id) & nzchar(feature_id)
  if (!any(valid_feature)) {
    stop("ID 列没有可用的基因/探针 ID。")
  }

  mat_df <- df[valid_feature, sample_cols, drop = FALSE]
  mat_df[] <- lapply(mat_df, numeric_clean)
  mat <- as.matrix(mat_df)
  storage.mode(mat) <- "double"

  good_cols <- colSums(is.finite(mat)) > 0
  mat <- mat[, good_cols, drop = FALSE]
  sample_cols <- sample_cols[good_cols]
  if (ncol(mat) < 2) {
    stop("可识别为数值的样本列少于 2 列。")
  }

  if (any(mat < 0, na.rm = TRUE)) {
    stop("raw count 不能包含负值。请确认数据类型是否应选择 TPM/FPKM 或芯片/已标准化表达矩阵。")
  }

  finite_values <- mat[is.finite(mat)]
  if (length(finite_values) == 0) {
    stop("raw count 矩阵没有可用数值。")
  }
  integer_like <- abs(finite_values - round(finite_values)) <= 1e-6
  if (mean(integer_like) < 0.99) {
    stop("raw count 应该基本都是整数。当前矩阵包含较多小数，更像 TPM/FPKM 或标准化表达量。")
  }
  mat <- round(mat)

  min_total_count <- numeric_clean(min_total_count)
  if (!is.finite(min_total_count) || min_total_count < 1) {
    min_total_count <- 1
  }
  min_total_count <- floor(min_total_count)

  good_rows <- rowSums(is.finite(mat)) >= 2 & rowSums(mat, na.rm = TRUE) >= min_total_count
  mat <- mat[good_rows, , drop = FALSE]
  feature_id <- feature_id[valid_feature][good_rows]
  if (nrow(mat) < 2) {
    stop("过滤全 0 或无效行后，可用于分析的基因少于 2 行。")
  }

  rownames(mat) <- make.unique(feature_id)
  colnames(mat) <- sample_cols

  list(
    counts = mat,
    feature_id = feature_id,
    transformed = FALSE,
    normalized = FALSE,
    min_total_count = min_total_count
  )
}

align_two_group_samples <- function(mat, group_df, control_group, treatment_group) {
  group_df <- group_df[group_df$sample %in% colnames(mat), , drop = FALSE]
  group_df <- group_df[match(colnames(mat), group_df$sample), , drop = FALSE]
  groups <- group_df$group

  keep_samples <- !is.na(groups) & groups %in% c(control_group, treatment_group)
  mat <- mat[, keep_samples, drop = FALSE]
  groups <- groups[keep_samples]
  if (ncol(mat) < 4) {
    stop("差异分析至少建议每组 2 个样本；当前可用样本总数少于 4。")
  }
  group_counts <- table(groups)
  if (!all(c(control_group, treatment_group) %in% names(group_counts)) ||
      any(group_counts[c(control_group, treatment_group)] < 2)) {
    stop("对照组和处理组都至少需要 2 个匹配样本。")
  }

  list(mat = mat, groups = groups)
}

add_symbols_and_change <- function(deg, feature_id, annotation, p_column,
                                   p_cutoff, logfc_cutoff) {
  deg$row_id <- rownames(deg)
  if (is.null(names(feature_id))) {
    names(feature_id) <- feature_id
  }
  deg$feature_id <- unname(feature_id[match(deg$row_id, names(feature_id))])
  missing_feature <- is.na(deg$feature_id) | !nzchar(deg$feature_id)
  deg$feature_id[missing_feature] <- deg$row_id[missing_feature]
  deg$symbol <- deg$feature_id

  if (!is.null(annotation)) {
    matched_symbol <- annotation$symbol[match(deg$feature_id, annotation$feature_id)]
    use_symbol <- !is.na(matched_symbol) & nzchar(matched_symbol)
    deg$symbol[use_symbol] <- matched_symbol[use_symbol]
  }

  if (!p_column %in% names(deg)) {
    stop("结果表中找不到显著性列：", p_column)
  }
  p_values <- deg[[p_column]]
  deg$change <- ifelse(
    !is.na(p_values) & p_values < p_cutoff & deg$logFC > logfc_cutoff,
    "up",
    ifelse(
      !is.na(p_values) & p_values < p_cutoff & deg$logFC < -logfc_cutoff,
      "down",
      "stable"
    )
  )
  deg$change <- factor(deg$change, levels = c("down", "stable", "up"))

  ordered_cols <- c(
    "row_id", "feature_id", "symbol", "logFC", "log2FoldChange",
    "AveExpr", "logCPM", "baseMean", "t", "stat", "LR", "F",
    "P.Value", "pvalue", "adj.P.Val", "padj", "FDR", "B", "change"
  )
  deg[, c(ordered_cols[ordered_cols %in% names(deg)],
          setdiff(names(deg), ordered_cols)), drop = FALSE]
}

make_analysis_result <- function(deg, mat, groups, feature_id, control_group,
                                 treatment_group, p_cutoff, logfc_cutoff,
                                 p_column, annotation, transformed, normalized,
                                 data_type, analysis_method, input_mat = NULL) {
  deg <- add_symbols_and_change(
    deg = deg,
    feature_id = feature_id,
    annotation = annotation,
    p_column = p_column,
    p_cutoff = p_cutoff,
    logfc_cutoff = logfc_cutoff
  )

  list(
    deg = deg,
    mat = mat,
    input_mat = input_mat,
    groups = groups,
    control_group = control_group,
    treatment_group = treatment_group,
    p_cutoff = p_cutoff,
    logfc_cutoff = logfc_cutoff,
    p_column = p_column,
    transformed = transformed,
    normalized = normalized,
    data_type = data_type,
    analysis_method = analysis_method
  )
}

parse_manual_groups <- function(text) {
  if (is.null(text) || !nzchar(trimws(text))) {
    return(data.frame(sample = character(), group = character()))
  }

  df <- tryCatch(
    read.csv(
      text = text,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ),
    error = function(e) NULL
  )

  if (is.null(df) || ncol(df) < 2) {
    df <- tryCatch(
      read.delim(
        text = text,
        stringsAsFactors = FALSE,
        check.names = FALSE
      ),
      error = function(e) NULL
    )
  }

  if (is.null(df) || ncol(df) < 2) {
    return(data.frame(sample = character(), group = character()))
  }

  data.frame(
    sample = trim_text(df[[1]]),
    group = trim_text(df[[2]]),
    stringsAsFactors = FALSE
  )
}

standardize_group_table <- function(df, sample_col = NULL, group_col = NULL) {
  if (is.null(df) || nrow(df) == 0) {
    return(data.frame(sample = character(), group = character()))
  }
  sample_col <- sample_col %||% names(df)[1]
  group_col <- group_col %||% names(df)[2]
  if (!sample_col %in% names(df) || !group_col %in% names(df)) {
    return(data.frame(sample = character(), group = character()))
  }
  out <- data.frame(
    sample = trim_text(df[[sample_col]]),
    group = trim_text(df[[group_col]]),
    stringsAsFactors = FALSE
  )
  out <- out[nzchar(out$sample) & nzchar(out$group), , drop = FALSE]
  out[!duplicated(out$sample), , drop = FALSE]
}

read_annotation <- function(file, default_annotation = NULL) {
  if (!is.null(file)) {
    df <- read_any_table(file)
  } else {
    df <- default_annotation
  }

  if (is.null(df) || nrow(df) == 0 || ncol(df) < 2) {
    return(NULL)
  }

  lower_names <- tolower(names(df))
  pick_col <- function(candidates, fallback) {
    idx <- match(tolower(candidates), lower_names)
    idx <- idx[!is.na(idx)]
    if (length(idx) > 0) {
      return(names(df)[idx[1]])
    }
    fallback
  }
  id_col <- pick_col(
    c("feature_id", "gene", "gene_id", "id", "ensembl_gene_id", "probe_id"),
    names(df)[1]
  )
  symbol_col <- pick_col(
    c("symbol", "gene_name", "gene_symbol", "hgnc_symbol", "external_gene_name"),
    names(df)[2]
  )

  out <- data.frame(
    feature_id = trim_text(df[[id_col]]),
    symbol = trim_text(df[[symbol_col]]),
    stringsAsFactors = FALSE
  )
  out <- out[nzchar(out$feature_id) & nzchar(out$symbol), , drop = FALSE]
  out[!duplicated(out$feature_id), , drop = FALSE]
}

impute_rows_for_plot <- function(mat) {
  out <- mat
  for (i in seq_len(nrow(out))) {
    missing <- !is.finite(out[i, ])
    if (any(missing)) {
      fill <- stats::median(out[i, is.finite(out[i, ])], na.rm = TRUE)
      if (!is.finite(fill)) {
        fill <- 0
      }
      out[i, missing] <- fill
    }
  }
  out
}

run_differential_analysis <- function(prepared, group_df, control_group,
                                      treatment_group, p_cutoff, logfc_cutoff,
                                      p_column, annotation) {
  mat <- prepared$mat
  group_df <- group_df[group_df$sample %in% colnames(mat), , drop = FALSE]
  group_df <- group_df[match(colnames(mat), group_df$sample), , drop = FALSE]
  groups <- group_df$group

  keep_samples <- !is.na(groups) & groups %in% c(control_group, treatment_group)
  mat <- mat[, keep_samples, drop = FALSE]
  groups <- groups[keep_samples]
  if (ncol(mat) < 4) {
    stop("差异分析至少建议每组 2 个样本；当前可用样本总数少于 4。")
  }
  group_counts <- table(groups)
  if (!all(c(control_group, treatment_group) %in% names(group_counts)) ||
      any(group_counts[c(control_group, treatment_group)] < 2)) {
    stop("对照组和处理组都至少需要 2 个匹配样本。")
  }

  row_ok <- apply(mat, 1, function(x) {
    counts <- tapply(is.finite(x), groups, sum)
    all(counts[c(control_group, treatment_group)] >= 2)
  })
  mat <- mat[row_ok, , drop = FALSE]
  feature_id <- prepared$feature_id[row_ok]
  input_mat <- prepared$input_mat[rownames(mat), colnames(mat), drop = FALSE]
  if (nrow(mat) < 2) {
    stop("按分组过滤后，可用于分析的基因/探针少于 2 行。")
  }

  group_factor <- factor(groups, levels = c(control_group, treatment_group))
  design <- model.matrix(~0 + group_factor)
  colnames(design) <- c("control", "treatment")

  fit <- limma::lmFit(mat, design)
  contrast <- limma::makeContrasts(treatment - control, levels = design)
  fit <- limma::contrasts.fit(fit, contrast)
  fit <- limma::eBayes(fit)
  deg <- limma::topTable(fit, number = Inf, sort.by = "P")

  deg$row_id <- rownames(deg)
  deg$feature_id <- feature_id[match(deg$row_id, rownames(mat))]
  deg$symbol <- deg$feature_id
  if (!is.null(annotation)) {
    matched_symbol <- annotation$symbol[match(deg$feature_id, annotation$feature_id)]
    use_symbol <- !is.na(matched_symbol) & nzchar(matched_symbol)
    deg$symbol[use_symbol] <- matched_symbol[use_symbol]
  }

  p_values <- deg[[p_column]]
  deg$change <- ifelse(
    !is.na(p_values) & p_values < p_cutoff & deg$logFC > logfc_cutoff,
    "up",
    ifelse(
      !is.na(p_values) & p_values < p_cutoff & deg$logFC < -logfc_cutoff,
      "down",
      "stable"
    )
  )
  deg$change <- factor(deg$change, levels = c("down", "stable", "up"))

  ordered_cols <- c(
    "row_id", "feature_id", "symbol", "logFC", "AveExpr", "t",
    "P.Value", "adj.P.Val", "B", "change"
  )
  deg <- deg[, c(ordered_cols[ordered_cols %in% names(deg)],
                 setdiff(names(deg), ordered_cols)), drop = FALSE]

  list(
    deg = deg,
    mat = mat,
    input_mat = input_mat,
    groups = groups,
    control_group = control_group,
    treatment_group = treatment_group,
    p_cutoff = p_cutoff,
    logfc_cutoff = logfc_cutoff,
    p_column = p_column,
    transformed = prepared$transformed,
    normalized = prepared$normalized,
    data_type = "芯片/已标准化表达矩阵",
    analysis_method = "limma"
  )
}

run_deseq2_analysis <- function(prepared, group_df, control_group,
                                treatment_group, p_cutoff, logfc_cutoff,
                                p_column, annotation) {
  aligned <- align_two_group_samples(
    prepared$counts,
    group_df,
    control_group,
    treatment_group
  )
  count_mat <- aligned$mat
  groups <- aligned$groups
  row_ok <- rowSums(count_mat, na.rm = TRUE) > 0
  count_mat <- count_mat[row_ok, , drop = FALSE]
  feature_id <- prepared$feature_id[match(rownames(count_mat), rownames(prepared$counts))]
  names(feature_id) <- rownames(count_mat)

  group_factor <- factor(groups, levels = c(control_group, treatment_group))
  col_data <- data.frame(group = group_factor, row.names = colnames(count_mat))
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = round(count_mat),
    colData = col_data,
    design = ~ group
  )
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  res <- DESeq2::results(
    dds,
    contrast = c("group", treatment_group, control_group)
  )
  deg <- as.data.frame(res)
  deg$logFC <- deg$log2FoldChange
  deg$P.Value <- deg$pvalue
  deg$adj.P.Val <- deg$padj
  deg$t <- deg$stat

  plot_mat <- log2(DESeq2::counts(dds, normalized = TRUE) + 1)
  deg$AveExpr <- rowMeans(plot_mat, na.rm = TRUE)

  make_analysis_result(
    deg = deg,
    mat = plot_mat,
    groups = groups,
    feature_id = feature_id,
    control_group = control_group,
    treatment_group = treatment_group,
    p_cutoff = p_cutoff,
    logfc_cutoff = logfc_cutoff,
    p_column = p_column,
    annotation = annotation,
    transformed = TRUE,
    normalized = TRUE,
    data_type = "RNA-seq raw counts",
    analysis_method = "DESeq2",
    input_mat = count_mat
  )
}

run_edger_analysis <- function(prepared, group_df, control_group,
                               treatment_group, p_cutoff, logfc_cutoff,
                               p_column, annotation) {
  aligned <- align_two_group_samples(
    prepared$counts,
    group_df,
    control_group,
    treatment_group
  )
  count_mat <- aligned$mat
  groups <- aligned$groups
  group_factor <- factor(groups, levels = c(control_group, treatment_group))

  dge <- edgeR::DGEList(counts = count_mat, group = group_factor)
  keep <- edgeR::filterByExpr(dge, group = group_factor)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  if (nrow(dge) < 2) {
    stop("edgeR 过滤低表达基因后，可用于分析的基因少于 2 行。")
  }
  feature_id <- prepared$feature_id[match(rownames(dge), rownames(prepared$counts))]
  names(feature_id) <- rownames(dge)
  dge <- edgeR::calcNormFactors(dge)

  design <- model.matrix(~0 + group_factor)
  colnames(design) <- c("control", "treatment")
  fit <- edgeR::glmQLFit(dge, design)
  qlf <- edgeR::glmQLFTest(fit, contrast = c(-1, 1))
  deg <- edgeR::topTags(qlf, n = Inf, sort.by = "PValue")$table
  deg$P.Value <- deg$PValue
  deg$adj.P.Val <- deg$FDR
  deg$t <- deg$F

  plot_mat <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
  deg$AveExpr <- rowMeans(plot_mat[rownames(deg), , drop = FALSE], na.rm = TRUE)

  make_analysis_result(
    deg = deg,
    mat = plot_mat,
    groups = groups,
    feature_id = feature_id,
    control_group = control_group,
    treatment_group = treatment_group,
    p_cutoff = p_cutoff,
    logfc_cutoff = logfc_cutoff,
    p_column = p_column,
    annotation = annotation,
    transformed = TRUE,
    normalized = TRUE,
    data_type = "RNA-seq raw counts",
    analysis_method = "edgeR",
    input_mat = count_mat
  )
}

run_voom_analysis <- function(prepared, group_df, control_group,
                              treatment_group, p_cutoff, logfc_cutoff,
                              p_column, annotation) {
  aligned <- align_two_group_samples(
    prepared$counts,
    group_df,
    control_group,
    treatment_group
  )
  count_mat <- aligned$mat
  groups <- aligned$groups
  group_factor <- factor(groups, levels = c(control_group, treatment_group))

  dge <- edgeR::DGEList(counts = count_mat, group = group_factor)
  keep <- edgeR::filterByExpr(dge, group = group_factor)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  if (nrow(dge) < 2) {
    stop("limma-voom 过滤低表达基因后，可用于分析的基因少于 2 行。")
  }
  feature_id <- prepared$feature_id[match(rownames(dge), rownames(prepared$counts))]
  names(feature_id) <- rownames(dge)
  dge <- edgeR::calcNormFactors(dge)

  design <- model.matrix(~0 + group_factor)
  colnames(design) <- c("control", "treatment")
  voom_fit <- limma::voom(dge, design, plot = FALSE)
  fit <- limma::lmFit(voom_fit, design)
  contrast <- limma::makeContrasts(treatment - control, levels = design)
  fit <- limma::contrasts.fit(fit, contrast)
  fit <- limma::eBayes(fit)
  deg <- limma::topTable(fit, number = Inf, sort.by = "P")

  make_analysis_result(
    deg = deg,
    mat = voom_fit$E,
    groups = groups,
    feature_id = feature_id,
    control_group = control_group,
    treatment_group = treatment_group,
    p_cutoff = p_cutoff,
    logfc_cutoff = logfc_cutoff,
    p_column = p_column,
    annotation = annotation,
    transformed = TRUE,
    normalized = TRUE,
    data_type = "RNA-seq raw counts",
    analysis_method = "limma-voom",
    input_mat = count_mat
  )
}

run_analysis_by_data_type <- function(expr_df, id_col, data_type, count_method,
                                      log_mode, normalize_between_arrays,
                                      group_df, control_group, treatment_group,
                                      p_cutoff, logfc_cutoff, p_column,
                                      annotation, count_min_total = 2) {
  auto_guess <- NULL
  requested_data_type <- data_type %||% "auto"
  if (identical(requested_data_type, "auto")) {
    auto_guess <- guess_expression_data_type(expr_df, id_col)
    data_type <- if (!is.null(auto_guess)) auto_guess$recommended else "normalized"
  }

  if (identical(data_type, "rnaseq_counts")) {
    prepared <- prepare_count_expression(expr_df, id_col, min_total_count = count_min_total)
    if (identical(count_method, "edger")) {
      out <- run_edger_analysis(
        prepared, group_df, control_group, treatment_group,
        p_cutoff, logfc_cutoff, p_column, annotation
      )
    } else if (identical(count_method, "voom")) {
      out <- run_voom_analysis(
        prepared, group_df, control_group, treatment_group,
        p_cutoff, logfc_cutoff, p_column, annotation
      )
    } else {
      out <- run_deseq2_analysis(
        prepared, group_df, control_group, treatment_group,
        p_cutoff, logfc_cutoff, p_column, annotation
      )
    }
    if (identical(requested_data_type, "auto")) {
      out$data_type <- paste0("自动识别：", out$data_type)
      out$auto_guess <- auto_guess
    }
    return(out)
  }

  normalize <- isTRUE(normalize_between_arrays) && identical(data_type, "normalized")
  prepared <- prepare_expression(expr_df, id_col, log_mode, normalize)
  result <- run_differential_analysis(
    prepared = prepared,
    group_df = group_df,
    control_group = control_group,
    treatment_group = treatment_group,
    p_cutoff = p_cutoff,
    logfc_cutoff = logfc_cutoff,
    p_column = p_column,
    annotation = annotation
  )

  if (identical(data_type, "rnaseq_tpm_fpkm")) {
    result$data_type <- "RNA-seq TPM/FPKM/RPKM"
    result$analysis_method <- "log2 表达量 + limma"
  } else {
    result$data_type <- "芯片/已标准化表达矩阵"
    result$analysis_method <- "limma"
  }
  if (identical(requested_data_type, "auto")) {
    result$data_type <- paste0("自动识别：", result$data_type)
    result$auto_guess <- auto_guess
  }
  result
}

guess_expression_data_type <- function(df, id_col) {
  if (is.null(df) || is.null(id_col) || !id_col %in% names(df)) {
    return(NULL)
  }
  sample_cols <- setdiff(names(df), id_col)
  if (length(sample_cols) < 2) {
    return(NULL)
  }
  mat_df <- df[, sample_cols, drop = FALSE]
  mat_df[] <- lapply(mat_df, numeric_clean)
  mat <- as.matrix(mat_df)
  storage.mode(mat) <- "double"
  values <- mat[is.finite(mat)]
  if (length(values) == 0) {
    return(NULL)
  }

  integer_fraction <- mean(abs(values - round(values)) <= 1e-6)
  zero_fraction <- mean(values == 0)
  min_value <- min(values, na.rm = TRUE)
  max_value <- max(values, na.rm = TRUE)
  q99 <- stats::quantile(values, 0.99, na.rm = TRUE)

  if (min_value < 0) {
    label <- "更像芯片/已标准化表达矩阵"
    reason <- "数据里有负值，raw count 和 TPM/FPKM 通常不会出现负值。"
    recommended <- "normalized"
  } else if (integer_fraction >= 0.99 && q99 >= 20) {
    label <- "更像 RNA-seq raw counts"
    reason <- "大部分值都是整数，并且高分位数较大，符合原始 reads/counts 的常见特征。"
    recommended <- "rnaseq_counts"
  } else if (integer_fraction < 0.99 && (q99 > 30 || max_value > 100)) {
    label <- "更像 RNA-seq TPM/FPKM/RPKM"
    reason <- "数据包含较多小数且数值范围较大，常见于未 log2 的 TPM/FPKM/RPKM。"
    recommended <- "rnaseq_tpm_fpkm"
  } else {
    label <- "更像芯片/已标准化表达矩阵，或已经 log2 的 TPM/FPKM"
    reason <- "数值整体较小，常见于芯片表达矩阵、GEO normalized matrix 或 log2 后表达量。"
    recommended <- "normalized"
  }

  list(
    label = label,
    reason = reason,
    recommended = recommended,
    integer_fraction = integer_fraction,
    zero_fraction = zero_fraction,
    min_value = min_value,
    max_value = max_value,
    q99 = as.numeric(q99)
  )
}

data_type_guide_ui <- function(guess = NULL) {
  guess_text <- if (is.null(guess)) {
    "导入表达矩阵后，这里会给出一个粗略判断。"
  } else {
    paste0(
      guess$label, "：", guess$reason,
      " 数值范围约为 ", signif(guess$min_value, 4), " 到 ",
      signif(guess$max_value, 4), "。"
    )
  }

  div(
    class = "analysis-note",
    h4("不知道自己的数据类型？"),
    p(strong("软件粗略判断："), guess_text),
    p(strong("默认自动模式："), "点击“一键开始分析”时，软件会按这个判断自动选择转换方式和差异分析方法。判断错了再手动切换数据类型。"),
    tags$ul(
      tags$li(strong("raw counts："), "值基本都是 0、1、2 这种整数，文件名常见 count/counts/readcount。正式 RNA-seq 差异分析优先选这个。"),
      tags$li(strong("TPM/FPKM/RPKM："), "值经常带小数，文件名常见 TPM、FPKM、RPKM。适合快速探索、作图和通路解释。"),
      tags$li(strong("芯片/已标准化表达："), "GEO 芯片矩阵或 normalized expression，值通常已经是 log2 后的小数，有时会出现负值。"),
      tags$li(strong("公司给的差异结果表："), "如果文件已经是 logFC、P 值、padj 那种结果表，它不是表达矩阵，不能直接当作输入矩阵。")
    ),
    p("拿不准时：先看文件名；再看数值是不是整数；最后看说明书或让测序公司确认 count、TPM/FPKM 还是 normalized matrix。")
  )
}

table_for_display <- function(deg) {
  out <- deg
  out$row_id <- NULL
  format_table_for_display(out)
}

format_p_decimal <- function(x, digits = 8) {
  x_num <- suppressWarnings(as.numeric(x))
  floor_value <- 10^-digits
  vapply(x_num, function(value) {
    if (!is.finite(value)) {
      return("")
    }
    if (value > 0 && abs(value) < floor_value) {
      return(paste0("<", formatC(floor_value, format = "f", digits = digits)))
    }
    txt <- formatC(value, format = "f", digits = digits)
    txt <- sub("0+$", "", txt)
    txt <- sub("\\.$", ".0", txt)
    txt
  }, character(1))
}

format_table_for_display <- function(df) {
  out <- df
  p_cols <- intersect(
    c("P.Value", "pvalue", "adj.P.Val", "padj", "p.adjust", "qvalue", "FDR", "p_value"),
    names(out)
  )
  for (col in p_cols) {
    out[[col]] <- format_p_decimal(out[[col]])
  }
  out
}

make_volcano <- function(result, colors = NULL, label_top_n = 0,
                         label_up = TRUE, label_down = TRUE,
                         label_up_n = NULL, label_down_n = NULL) {
  deg <- result$deg
  p_values <- pmax(deg[[result$p_column]], .Machine$double.xmin)
  deg$neg_log10_p <- -log10(p_values)
  colors <- colors %||% c(down = "#2563eb", stable = "#9ca3af", up = "#dc2626")

  plot <- ggplot(deg, aes(x = logFC, y = neg_log10_p, color = change)) +
    geom_point(alpha = 0.65, size = 1.7) +
    geom_vline(
      xintercept = c(-result$logfc_cutoff, result$logfc_cutoff),
      linetype = "dashed",
      linewidth = 0.45,
      color = "#444444"
    ) +
    geom_hline(
      yintercept = -log10(result$p_cutoff),
      linetype = "dashed",
      linewidth = 0.45,
      color = "#444444"
    ) +
    scale_color_manual(
      values = colors,
      drop = FALSE
    ) +
    labs(
      x = paste0("log2FC (", result$treatment_group, " / ",
                 result$control_group, ")"),
      y = paste0("-log10(", result$p_column, ")"),
      color = "Change"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "top"
    )

  label_up_n <- label_up_n %||% label_top_n
  label_down_n <- label_down_n %||% label_top_n

  if (requireNamespace("ggrepel", quietly = TRUE) &&
      ((isTRUE(label_up) && label_up_n > 0) || (isTRUE(label_down) && label_down_n > 0))) {
    up_df <- data.frame()
    down_df <- data.frame()
    if (isTRUE(label_up) && label_up_n > 0) {
      up_df <- deg[deg$change == "up" & is.finite(deg[[result$p_column]]), , drop = FALSE]
      up_df <- up_df[order(up_df[[result$p_column]], -abs(up_df$logFC)), , drop = FALSE]
      up_df <- head(up_df, label_up_n)
    }
    if (isTRUE(label_down) && label_down_n > 0) {
      down_df <- deg[deg$change == "down" & is.finite(deg[[result$p_column]]), , drop = FALSE]
      down_df <- down_df[order(down_df[[result$p_column]], -abs(down_df$logFC)), , drop = FALSE]
      down_df <- head(down_df, label_down_n)
    }
    label_parts <- Filter(function(x) nrow(x) > 0, list(up_df, down_df))
    label_df <- if (length(label_parts) > 0) do.call(rbind, label_parts) else data.frame()
    if (nrow(label_df) > 0) {
      label_df <- label_df[!duplicated(label_df$row_id), , drop = FALSE]
    }
    if (nrow(label_df) > 0) {
      plot <- plot +
        ggrepel::geom_text_repel(
          data = label_df,
          aes(label = symbol, color = change),
          size = 3.2,
          max.overlaps = Inf,
          box.padding = 0.35,
          point.padding = 0.2,
          min.segment.length = 0,
          show.legend = FALSE
        )
    }
  }

  plot
}

make_pca_df <- function(result) {
  mat <- impute_rows_for_plot(result$mat)
  vars <- matrixStats::rowVars(mat)
  keep <- is.finite(vars) & vars > 0
  mat <- mat[keep, , drop = FALSE]
  vars <- vars[keep]
  if (nrow(mat) > 2000) {
    mat <- mat[order(vars, decreasing = TRUE)[seq_len(2000)], , drop = FALSE]
  }
  if (nrow(mat) < 2 || ncol(mat) < 3) {
    return(NULL)
  }

  pc <- stats::prcomp(t(mat), center = TRUE, scale. = FALSE)
  if (ncol(pc$x) < 2) {
    return(NULL)
  }
  explained <- summary(pc)$importance[2, 1:2] * 100
  data.frame(
    sample = rownames(pc$x),
    group = result$groups,
    PC1 = pc$x[, 1],
    PC2 = pc$x[, 2],
    PC1_label = sprintf("PC1 (%.1f%%)", explained[1]),
    PC2_label = sprintf("PC2 (%.1f%%)", explained[2]),
    stringsAsFactors = FALSE
  )
}

make_pca_plot <- function(result, show_ellipse = FALSE, show_centers = FALSE) {
  pca_df <- make_pca_df(result)
  if (is.null(pca_df)) {
    return(NULL)
  }
  plot <- ggplot(pca_df, aes(PC1, PC2, color = group, label = sample))
  if (isTRUE(show_ellipse)) {
    group_counts <- table(pca_df$group)
    ellipse_groups <- names(group_counts[group_counts >= 3])
    if (length(ellipse_groups) > 0) {
      plot <- plot +
        stat_ellipse(
          data = pca_df[pca_df$group %in% ellipse_groups, , drop = FALSE],
          aes(group = group),
          type = "norm",
          level = 0.95,
          linewidth = 0.75,
          alpha = 0.7,
          show.legend = FALSE
        )
    }
  }
  plot <- plot +
    geom_point(size = 3, alpha = 0.85) +
    labs(
      x = pca_df$PC1_label[1],
      y = pca_df$PC2_label[1],
      color = "Group"
    ) +
    theme_minimal(base_size = 13) +
    theme(panel.grid.minor = element_blank(), legend.position = "top")
  if (isTRUE(show_centers)) {
    centers <- aggregate(cbind(PC1, PC2) ~ group, data = pca_df, FUN = mean)
    plot <- plot +
      geom_point(
        data = centers,
        aes(PC1, PC2, color = group),
        inherit.aes = FALSE,
        shape = 4,
        size = 5,
        stroke = 1.6,
        show.legend = FALSE
      )
  }
  plot
}

heatmap_matrix <- function(result, top_n = 50) {
  deg <- result$deg
  selected <- deg$row_id[deg$change != "stable"]
  if (length(selected) < 2) {
    selected <- deg$row_id[order(deg[[result$p_column]])]
  }
  selected <- head(selected, top_n)
  n <- result$mat[selected, , drop = FALSE]
  row_labels <- deg$symbol[match(selected, deg$row_id)]
  row_labels[is.na(row_labels) | !nzchar(row_labels)] <- selected[
    is.na(row_labels) | !nzchar(row_labels)
  ]
  rownames(n) <- make.unique(row_labels)
  n <- impute_rows_for_plot(n)
  row_sd <- matrixStats::rowSds(n)
  n[row_sd > 0 & is.finite(row_sd), , drop = FALSE]
}

boxplot_matrix_for_scale <- function(result, scale = "analysis", pseudocount = 0.001) {
  scale <- scale %||% "analysis"
  if (identical(scale, "input")) {
    mat <- result$input_mat %||% result$mat
    y_label <- "Input expression"
  } else if (identical(scale, "log10_input")) {
    mat <- result$input_mat %||% result$mat
    mat[mat < 0] <- NA_real_
    mat <- log10(mat + pseudocount)
    y_label <- paste0("log10(input value + ", pseudocount, ")")
  } else if (identical(scale, "median_centered")) {
    mat <- result$mat
    mat <- sweep(mat, 2, apply(mat, 2, stats::median, na.rm = TRUE), "-")
    y_label <- "Expression minus sample median"
  } else {
    mat <- result$mat
    y_label <- "Expression"
  }
  list(mat = mat, y_label = y_label)
}

make_boxplot_df <- function(result, max_genes = 10000, scale = "analysis",
                            pseudocount = 0.001) {
  scaled <- boxplot_matrix_for_scale(result, scale = scale, pseudocount = pseudocount)
  mat <- impute_rows_for_plot(scaled$mat)
  keep_rows <- rowSums(is.finite(mat)) >= 2
  mat <- mat[keep_rows, , drop = FALSE]
  if (nrow(mat) == 0 || ncol(mat) == 0) {
    return(data.frame())
  }
  max_genes <- as.integer(max_genes %||% 10000)
  if (is.finite(max_genes) && nrow(mat) > max_genes) {
    set.seed(2026)
    mat <- mat[sort(sample(seq_len(nrow(mat)), max_genes)), , drop = FALSE]
  }
  data.frame(
    sample = rep(colnames(mat), each = nrow(mat)),
    group = rep(result$groups, each = nrow(mat)),
    expression = as.vector(mat),
    y_label = scaled$y_label,
    stringsAsFactors = FALSE
  )
}

make_boxplot_plot <- function(result, max_genes = 10000, show_outliers = FALSE,
                              scale = "analysis", pseudocount = 0.001) {
  df <- make_boxplot_df(
    result,
    max_genes = max_genes,
    scale = scale,
    pseudocount = pseudocount
  )
  validate(need(nrow(df) > 0, "当前数据不足以绘制箱线图。"))
  ggplot(df, aes(x = sample, y = expression, fill = group)) +
    geom_boxplot(
      outlier.shape = if (isTRUE(show_outliers)) 16 else NA,
      outlier.size = 0.5,
      linewidth = 0.35
    ) +
    stat_summary(fun = median, geom = "point", shape = 21, size = 1.8, fill = "#111827", color = "#ffffff") +
    labs(x = "Sample", y = unique(df$y_label)[1], fill = "Group") +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "top",
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
    )
}

expression_with_symbols <- function(result) {
  mat <- result$mat
  deg <- result$deg
  labels <- deg$symbol[match(rownames(mat), deg$row_id)]
  labels[is.na(labels) | !nzchar(labels)] <- rownames(mat)[is.na(labels) | !nzchar(labels)]
  rownames(mat) <- make.unique(labels)
  mat
}

run_wgcna_analysis <- function(result, top_n = 5000, soft_power = 0,
                               min_module_size = 30, merge_cut_height = 0.25) {
  tryCatch(
    WGCNA::allowWGCNAThreads(nThreads = 2),
    error = function(e) {
      tryCatch(WGCNA::disableWGCNAThreads(), error = function(...) NULL)
    }
  )
  mat <- expression_with_symbols(result)
  mat <- impute_rows_for_plot(mat)
  vars <- matrixStats::rowVars(mat)
  keep <- is.finite(vars) & vars > 0
  mat <- mat[keep, , drop = FALSE]
  vars <- vars[keep]
  if (nrow(mat) > top_n) {
    mat <- mat[order(vars, decreasing = TRUE)[seq_len(top_n)], , drop = FALSE]
  }
  if (nrow(mat) < 50 || ncol(mat) < 8) {
    stop("WGCNA 至少建议 50 个基因和 8 个样本。当前过滤后数据不足。")
  }

  dat_expr <- as.data.frame(t(mat), check.names = FALSE)
  gsg <- WGCNA::goodSamplesGenes(dat_expr, verbose = 0)
  if (!gsg$allOK) {
    dat_expr <- dat_expr[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
  }
  groups <- result$groups[match(rownames(dat_expr), colnames(result$mat))]
  trait <- ifelse(groups == result$treatment_group, 1, 0)
  names(trait) <- rownames(dat_expr)

  powers <- c(1:10, seq(12, 20, 2))
  chosen_power <- as.integer(soft_power)
  sft_table <- NULL
  if (!is.finite(chosen_power) || chosen_power <= 0) {
    sft <- WGCNA::pickSoftThreshold(
      dat_expr,
      powerVector = powers,
      networkType = "signed",
      verbose = 0
    )
    sft_table <- as.data.frame(sft$fitIndices)
    fit_col <- "SFT.R.sq"
    candidates <- sft_table$Power[is.finite(sft_table[[fit_col]]) & sft_table[[fit_col]] >= 0.8]
    if (length(candidates) > 0) {
      chosen_power <- candidates[1]
    } else {
      chosen_power <- sft_table$Power[which.max(sft_table[[fit_col]])]
    }
    if (!is.finite(chosen_power) || length(chosen_power) == 0) {
      chosen_power <- 6
    }
  }

  net <- WGCNA::blockwiseModules(
    dat_expr,
    power = chosen_power,
    networkType = "signed",
    TOMType = "signed",
    minModuleSize = min_module_size,
    reassignThreshold = 0,
    mergeCutHeight = merge_cut_height,
    numericLabels = FALSE,
    pamRespectsDendro = FALSE,
    saveTOMs = FALSE,
    verbose = 0
  )

  module_colors <- net$colors
  names(module_colors) <- colnames(dat_expr)
  MEs <- WGCNA::orderMEs(WGCNA::moduleEigengenes(dat_expr, module_colors)$eigengenes)
  module_cor <- as.numeric(stats::cor(MEs, trait, use = "p"))
  module_p <- as.numeric(WGCNA::corPvalueStudent(module_cor, nrow(dat_expr)))
  module_names <- sub("^ME", "", colnames(MEs))

  deg <- result$deg
  gene_table <- data.frame(
    symbol = names(module_colors),
    module = as.character(module_colors),
    stringsAsFactors = FALSE
  )
  gene_table$change <- deg$change[match(gene_table$symbol, deg$symbol)]
  gene_table$logFC <- deg$logFC[match(gene_table$symbol, deg$symbol)]
  gene_table$P.Value <- deg$P.Value[match(gene_table$symbol, deg$symbol)]
  gene_table$change <- as.character(gene_table$change)
  gene_table$change[is.na(gene_table$change)] <- "not_tested"

  kme <- WGCNA::signedKME(dat_expr, MEs, outputColumnName = "kME")
  gene_table$kME <- NA_real_
  for (module in unique(gene_table$module)) {
    col <- paste0("kME", module)
    if (col %in% colnames(kme)) {
      idx <- gene_table$module == module
      gene_table$kME[idx] <- kme[gene_table$symbol[idx], col]
    }
  }

  module_summary <- data.frame(
    module = module_names,
    gene_count = as.integer(table(factor(module_colors, levels = module_names))),
    trait_correlation = module_cor,
    p_value = module_p,
    bias = ifelse(module_cor > 0, result$treatment_group, result$control_group),
    stringsAsFactors = FALSE
  )
  module_summary$direction <- ifelse(
    module_summary$trait_correlation > 0,
    paste0("偏向 ", result$treatment_group),
    paste0("偏向 ", result$control_group)
  )
  module_summary$strength <- cut(
    abs(module_summary$trait_correlation),
    breaks = c(-Inf, 0.3, 0.5, 0.7, Inf),
    labels = c("weak", "moderate", "strong", "very strong")
  )
  module_summary$up_genes <- vapply(module_summary$module, function(module) {
    sum(gene_table$module == module & gene_table$change == "up", na.rm = TRUE)
  }, integer(1))
  module_summary$down_genes <- vapply(module_summary$module, function(module) {
    sum(gene_table$module == module & gene_table$change == "down", na.rm = TRUE)
  }, integer(1))
  module_summary <- module_summary[order(module_summary$p_value, -abs(module_summary$trait_correlation)), ]

  hub_genes <- gene_table[gene_table$module != "grey" & is.finite(gene_table$kME), , drop = FALSE]
  hub_genes <- hub_genes[order(hub_genes$module, -abs(hub_genes$kME)), , drop = FALSE]
  hub_genes <- do.call(rbind, lapply(split(hub_genes, hub_genes$module), head, 20))
  rownames(hub_genes) <- NULL

  list(
    module_summary = module_summary,
    gene_table = gene_table,
    hub_genes = hub_genes,
    MEs = MEs,
    trait = trait,
    power = chosen_power,
    sft_table = sft_table,
    control_group = result$control_group,
    treatment_group = result$treatment_group
  )
}

make_wgcna_plot <- function(wgcna_result) {
  df <- wgcna_result$module_summary
  df$module <- factor(df$module, levels = rev(df$module))
  ggplot(df, aes(x = "Group trait", y = module, fill = trait_correlation)) +
    geom_tile(color = "white", linewidth = 0.7) +
    geom_text(aes(label = paste0("r=", sprintf("%.2f", trait_correlation), "\np=", format_p_decimal(p_value, digits = 4))), size = 3.4) +
    scale_fill_gradient2(low = "#2563eb", mid = "white", high = "#dc2626", limits = c(-1, 1)) +
    labs(
      x = paste0(wgcna_result$treatment_group, " = 1, ", wgcna_result$control_group, " = 0"),
      y = "WGCNA module",
      fill = "correlation"
    ) +
    theme_minimal(base_size = 13) +
    theme(panel.grid = element_blank())
}

wgcna_interpretation <- function(wgcna_result) {
  df <- wgcna_result$module_summary
  df <- df[df$module != "grey", , drop = FALSE]
  if (nrow(df) == 0) {
    return("WGCNA 没有识别到非 grey 模块，当前数据的共表达模块结构较弱或参数过严。")
  }
  top <- df[order(df$p_value, -abs(df$trait_correlation)), ][1, ]
  direction <- if (top$trait_correlation > 0) {
    paste0("更偏向 ", wgcna_result$treatment_group, " 组")
  } else {
    paste0("更偏向 ", wgcna_result$control_group, " 组")
  }
  deg_bias <- if (top$up_genes > top$down_genes) {
    paste0("该模块内上调基因更多，说明它更接近 ", wgcna_result$treatment_group, " 相关表达升高信号。")
  } else if (top$down_genes > top$up_genes) {
    paste0("该模块内下调基因更多，说明它更接近 ", wgcna_result$control_group, " 相关表达升高信号。")
  } else {
    "该模块内上调和下调差异基因数量接近，主要体现共表达结构而不是单一上下调方向。"
  }
  paste0(
    "WGCNA 将样本分组编码为 ", wgcna_result$treatment_group, "=1、",
    wgcna_result$control_group, "=0。当前最相关模块是 ",
    top$module, "，模块-分组相关 r=", sprintf("%.2f", top$trait_correlation),
    "，p=", format_p_decimal(top$p_value, digits = 4), "，", direction, "。", deg_bias,
    " 需要注意，WGCNA 解释的是“成组共表达模块”偏向哪一类样本，不等同于单个基因的差异分析。"
  )
}

read_ppi_table <- function(file = NULL, default_path = "string_interactions.tsv") {
  if (!is.null(file)) {
    df <- read_any_table(file)
  } else if (file.exists(default_path)) {
    df <- read.delim(default_path, check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    return(NULL)
  }
  names(df) <- sub("^#", "", names(df))
  if (!all(c("node1", "node2") %in% names(df))) {
    stop("PPI 文件至少需要 node1 和 node2 两列。")
  }
  if (!"combined_score" %in% names(df)) {
    df$combined_score <- 1
  }
  df$node1 <- trim_text(df$node1)
  df$node2 <- trim_text(df$node2)
  df$combined_score <- numeric_clean(df$combined_score)
  df[!is.na(df$node1) & !is.na(df$node2) & nzchar(df$node1) & nzchar(df$node2), , drop = FALSE]
}

run_ppi_analysis <- function(result, ppi_table, score_cutoff = 0.7, max_genes = 1000) {
  if (is.null(ppi_table) || nrow(ppi_table) == 0) {
    stop("没有可用的 PPI 互作表。请上传 STRING 导出的 interaction 文件，或把 string_interactions.tsv 放在软件目录。")
  }
  deg <- result$deg
  sig <- deg[deg$change != "stable" & !is.na(deg$symbol) & nzchar(deg$symbol), , drop = FALSE]
  sig <- sig[order(sig[[result$p_column]], -abs(sig$logFC)), , drop = FALSE]
  sig <- sig[!duplicated(sig$symbol), , drop = FALSE]
  if (nrow(sig) > max_genes) {
    sig <- sig[seq_len(max_genes), , drop = FALSE]
  }
  genes <- unique(sig$symbol)
  if (length(genes) < 2) {
    stop("显著差异基因少于 2 个，无法构建 PPI 网络。")
  }

  score <- ppi_table$combined_score
  if (max(score, na.rm = TRUE) > 1) {
    score <- score / 1000
  }
  edges <- ppi_table[ppi_table$node1 %in% genes & ppi_table$node2 %in% genes & score >= score_cutoff, , drop = FALSE]
  edges$score <- score[ppi_table$node1 %in% genes & ppi_table$node2 %in% genes & score >= score_cutoff]
  if (nrow(edges) == 0) {
    stop("当前阈值下没有 PPI 边。请降低 combined_score 阈值或增加 PPI 基因数。")
  }

  vertices <- data.frame(name = unique(c(edges$node1, edges$node2)), stringsAsFactors = FALSE)
  vertices$change <- sig$change[match(vertices$name, sig$symbol)]
  vertices$logFC <- sig$logFC[match(vertices$name, sig$symbol)]
  graph <- igraph::graph_from_data_frame(edges[, c("node1", "node2", "score")], directed = FALSE, vertices = vertices)
  graph <- igraph::simplify(graph, remove.multiple = TRUE, remove.loops = TRUE, edge.attr.comb = "max")

  centrality <- data.frame(
    symbol = igraph::V(graph)$name,
    degree = as.numeric(igraph::degree(graph)),
    betweenness = as.numeric(igraph::betweenness(graph, normalized = TRUE)),
    closeness = as.numeric(igraph::closeness(graph, normalized = TRUE)),
    stringsAsFactors = FALSE
  )
  centrality$change <- igraph::V(graph)$change
  centrality$logFC <- igraph::V(graph)$logFC
  centrality <- centrality[order(-centrality$degree, -centrality$betweenness, -abs(centrality$logFC)), ]
  rownames(centrality) <- NULL

  communities <- tryCatch(igraph::cluster_louvain(graph), error = function(e) NULL)
  if (!is.null(communities)) {
    centrality$community <- igraph::membership(communities)[centrality$symbol]
  } else {
    centrality$community <- NA_integer_
  }

  list(
    graph = graph,
    hubs = centrality,
    edges = edges,
    score_cutoff = score_cutoff,
    max_genes = max_genes,
    control_group = result$control_group,
    treatment_group = result$treatment_group
  )
}

parse_gene_list <- function(text) {
  if (is.null(text) || !nzchar(trimws(text))) {
    return(character())
  }
  genes <- unlist(strsplit(text, "[,;[:space:]]+"))
  genes <- trim_text(genes)
  unique(genes[!is.na(genes) & nzchar(genes)])
}

map_feature_ids_to_symbols <- function(ids) {
  ids <- unique(trim_text(ids))
  ids <- ids[!is.na(ids) & nzchar(ids)]
  if (length(ids) == 0) {
    return(data.frame(feature_id = character(), symbol = character()))
  }
  out <- data.frame(feature_id = character(), symbol = character())
  ensembl_ids <- ids[grepl("^ENSG[0-9]+", ids, ignore.case = TRUE)]
  entrez_ids <- ids[grepl("^[0-9]+$", ids)]
  if (length(ensembl_ids) > 0) {
    mapped <- tryCatch(
      suppressMessages(clusterProfiler::bitr(
        ensembl_ids,
        fromType = "ENSEMBL",
        toType = "SYMBOL",
        OrgDb = org.Hs.eg.db::org.Hs.eg.db
      )),
      error = function(e) data.frame()
    )
    if (nrow(mapped) > 0) {
      out <- rbind(out, data.frame(feature_id = mapped$ENSEMBL, symbol = mapped$SYMBOL))
    }
  }
  if (length(entrez_ids) > 0) {
    mapped <- tryCatch(
      suppressMessages(clusterProfiler::bitr(
        entrez_ids,
        fromType = "ENTREZID",
        toType = "SYMBOL",
        OrgDb = org.Hs.eg.db::org.Hs.eg.db
      )),
      error = function(e) data.frame()
    )
    if (nrow(mapped) > 0) {
      out <- rbind(out, data.frame(feature_id = mapped$ENTREZID, symbol = mapped$SYMBOL))
    }
  }
  out <- out[!duplicated(out$feature_id), , drop = FALSE]
  out
}

fill_ppi_symbols <- function(gene_info) {
  same_as_feature <- !is.na(gene_info$symbol) & !is.na(gene_info$feature_id) &
    gene_info$symbol == gene_info$feature_id
  needs_mapping <- is.na(gene_info$symbol) | !nzchar(gene_info$symbol) |
    same_as_feature |
    grepl("^ENSG[0-9]+", gene_info$symbol, ignore.case = TRUE) |
    grepl("^[0-9]+$", gene_info$symbol)
  ids_to_map <- unique(c(gene_info$feature_id[needs_mapping], gene_info$query[needs_mapping]))
  mapped <- map_feature_ids_to_symbols(ids_to_map)
  if (nrow(mapped) > 0) {
    by_feature <- mapped$symbol[match(gene_info$feature_id, mapped$feature_id)]
    by_query <- mapped$symbol[match(gene_info$query, mapped$feature_id)]
    use_feature <- !is.na(by_feature) & nzchar(by_feature)
    use_query <- !use_feature & !is.na(by_query) & nzchar(by_query)
    gene_info$symbol[use_feature] <- by_feature[use_feature]
    gene_info$symbol[use_query] <- by_query[use_query]
  }
  gene_info
}

ppi_gene_choices <- function(result) {
  deg <- result$deg
  deg$symbol <- as.character(deg$symbol)
  deg$feature_id <- as.character(deg$feature_id)
  labels <- paste0(
    deg$symbol,
    " | ID: ", deg$feature_id,
    " | ", deg$change,
    " | logFC=", round(deg$logFC, 3),
    " | padj=", format_p_decimal(deg$adj.P.Val, digits = 6)
  )
  values <- deg$feature_id
  valid <- !is.na(values) & nzchar(values)
  stats::setNames(values[valid], labels[valid])
}

ppi_gene_table <- function(result, gene_source = "significant",
                           custom_gene_text = NULL, selected_genes = NULL,
                           max_genes = 1000) {
  deg <- result$deg
  deg$change <- as.character(deg$change)
  if (identical(gene_source, "selected") || identical(gene_source, "custom")) {
    requested <- if (identical(gene_source, "selected")) {
      unique(trim_text(selected_genes %||% character()))
    } else {
      parse_gene_list(custom_gene_text)
    }
    if (length(requested) < 2) {
      stop("自选 PPI 至少需要选择 2 个当前矩阵里的基因。")
    }
    by_symbol <- deg[match(toupper(requested), toupper(deg$symbol)), , drop = FALSE]
    by_feature <- deg[match(toupper(requested), toupper(deg$feature_id)), , drop = FALSE]
    use_feature <- is.na(by_symbol$symbol) | !nzchar(by_symbol$symbol)
    matched <- by_symbol
    matched[use_feature, ] <- by_feature[use_feature, ]
    out <- data.frame(
      query = requested,
      symbol = matched$symbol,
      feature_id = matched$feature_id,
      change = matched$change,
      logFC = matched$logFC,
      P.Value = matched$P.Value,
      adj.P.Val = matched$adj.P.Val,
      stringsAsFactors = FALSE
    )
    missing <- is.na(out$symbol) | !nzchar(out$symbol)
    out$symbol[missing] <- out$query[missing]
    out$feature_id[missing] <- out$query[missing]
    out$change[is.na(out$change) | !nzchar(out$change)] <- "not_tested"
    out <- fill_ppi_symbols(out)
    return(out[!duplicated(out$symbol), , drop = FALSE])
  }

  sig <- deg[deg$change != "stable" & !is.na(deg$symbol) & nzchar(deg$symbol), , drop = FALSE]
  sig <- sig[order(sig[[result$p_column]], -abs(sig$logFC)), , drop = FALSE]
  sig <- sig[!duplicated(sig$symbol), , drop = FALSE]
  if (nrow(sig) > max_genes) {
    sig <- sig[seq_len(max_genes), , drop = FALSE]
  }
  out <- data.frame(
    query = sig$symbol,
    symbol = sig$symbol,
    feature_id = sig$feature_id,
    change = sig$change,
    logFC = sig$logFC,
    P.Value = sig$P.Value,
    adj.P.Val = sig$adj.P.Val,
    stringsAsFactors = FALSE
  )
  out <- fill_ppi_symbols(out)
  out[!duplicated(out$symbol), , drop = FALSE]
}

ppi_alias_map <- function(gene_info) {
  rows <- do.call(rbind, lapply(seq_len(nrow(gene_info)), function(i) {
    aliases <- unique(trim_text(c(
      gene_info$symbol[i],
      gene_info$feature_id[i],
      gene_info$query[i]
    )))
    aliases <- aliases[!is.na(aliases) & nzchar(aliases)]
    data.frame(
      alias = aliases,
      alias_key = toupper(aliases),
      symbol = gene_info$symbol[i],
      stringsAsFactors = FALSE
    )
  }))
  rows <- rows[!duplicated(rows$alias_key), , drop = FALSE]
  rows
}

normalize_ppi_edges_to_genes <- function(ppi_table, gene_info) {
  aliases <- ppi_alias_map(gene_info)
  node1 <- trim_text(ppi_table$node1)
  node2 <- trim_text(ppi_table$node2)
  ppi_table$node1_original <- node1
  ppi_table$node2_original <- node2
  ppi_table$node1 <- aliases$symbol[match(toupper(node1), aliases$alias_key)]
  ppi_table$node2 <- aliases$symbol[match(toupper(node2), aliases$alias_key)]
  ppi_table
}

fetch_string_interactions <- function(genes, score_cutoff = 0.7) {
  genes <- unique(genes[!is.na(genes) & nzchar(genes)])
  if (length(genes) < 2) {
    return(data.frame())
  }
  if (length(genes) > 400) {
    genes <- genes[seq_len(400)]
  }
  required_score <- max(0, min(1000, round(score_cutoff * 1000)))
  url <- paste0(
    "https://string-db.org/api/tsv/network?identifiers=",
    utils::URLencode(paste(genes, collapse = "\r"), reserved = TRUE),
    "&species=9606&required_score=", required_score,
    "&caller_identity=BioInsight"
  )
  df <- tryCatch(
    read.delim(url, check.names = FALSE, stringsAsFactors = FALSE),
    error = function(e) {
      stop("无法在线读取 STRING 互作。请检查网络，或上传 STRING interaction 文件。")
    }
  )
  if (nrow(df) == 0) {
    return(data.frame())
  }
  node1_col <- if ("preferredName_A" %in% names(df)) "preferredName_A" else "stringId_A"
  node2_col <- if ("preferredName_B" %in% names(df)) "preferredName_B" else "stringId_B"
  score_col <- if ("score" %in% names(df)) "score" else if ("combined_score" %in% names(df)) "combined_score" else NULL
  out <- data.frame(
    node1 = trim_text(df[[node1_col]]),
    node2 = trim_text(df[[node2_col]]),
    combined_score = if (!is.null(score_col)) numeric_clean(df[[score_col]]) else 1,
    stringsAsFactors = FALSE
  )
  out[!is.na(out$node1) & !is.na(out$node2) & nzchar(out$node1) & nzchar(out$node2), , drop = FALSE]
}

run_ppi_analysis <- function(result, ppi_table = NULL, score_cutoff = 0.7,
                             max_genes = 1000, gene_source = "significant",
                             custom_gene_text = NULL, selected_genes = NULL,
                             interaction_source = "local") {
  gene_info <- ppi_gene_table(
    result,
    gene_source = gene_source,
    custom_gene_text = custom_gene_text,
    selected_genes = selected_genes,
    max_genes = max_genes
  )
  genes <- unique(gene_info$symbol[!is.na(gene_info$symbol) & nzchar(gene_info$symbol)])
  if (length(genes) < 2) {
    stop("可用于 PPI 的基因少于 2 个，无法构建网络。")
  }

  if (identical(interaction_source, "string_online")) {
    ppi_table <- fetch_string_interactions(genes, score_cutoff = score_cutoff)
  }
  if (is.null(ppi_table) || nrow(ppi_table) == 0) {
    stop("没有可用的 PPI 互作表。可选择在线 STRING，或上传 STRING interaction 文件。")
  }
  ppi_table <- normalize_ppi_edges_to_genes(ppi_table, gene_info)

  score <- ppi_table$combined_score
  if (length(score) == 0 || all(!is.finite(score))) {
    score <- rep(1, nrow(ppi_table))
  }
  if (max(score, na.rm = TRUE) > 1) {
    score <- score / 1000
  }
  keep_edge <- !is.na(ppi_table$node1) & !is.na(ppi_table$node2) &
    ppi_table$node1 %in% genes & ppi_table$node2 %in% genes &
    ppi_table$node1 != ppi_table$node2 & score >= score_cutoff
  edges <- ppi_table[keep_edge, , drop = FALSE]
  edges$score <- score[keep_edge]
  if (nrow(edges) == 0) {
    stop("当前阈值下没有 PPI 边。请降低 combined_score 阈值、增加基因数，或改用在线 STRING。")
  }

  vertices <- data.frame(name = unique(c(edges$node1, edges$node2)), stringsAsFactors = FALSE)
  vertices$change <- gene_info$change[match(vertices$name, gene_info$symbol)]
  vertices$logFC <- gene_info$logFC[match(vertices$name, gene_info$symbol)]
  vertices$feature_id <- gene_info$feature_id[match(vertices$name, gene_info$symbol)]
  vertices$query <- gene_info$query[match(vertices$name, gene_info$symbol)]
  vertices$change[is.na(vertices$change) | !nzchar(vertices$change)] <- "not_tested"
  graph <- igraph::graph_from_data_frame(edges[, c("node1", "node2", "score")], directed = FALSE, vertices = vertices)
  graph <- igraph::simplify(graph, remove.multiple = TRUE, remove.loops = TRUE, edge.attr.comb = "max")

  centrality <- data.frame(
    symbol = igraph::V(graph)$name,
    degree = as.numeric(igraph::degree(graph)),
    betweenness = as.numeric(igraph::betweenness(graph, normalized = TRUE)),
    closeness = as.numeric(igraph::closeness(graph, normalized = TRUE)),
    stringsAsFactors = FALSE
  )
  centrality$change <- igraph::V(graph)$change
  centrality$logFC <- igraph::V(graph)$logFC
  centrality$abs_logFC <- abs(centrality$logFC)
  centrality$feature_id <- igraph::V(graph)$feature_id
  centrality$query <- igraph::V(graph)$query
  centrality <- centrality[order(-centrality$degree, -centrality$betweenness, -abs(centrality$logFC)), ]
  rownames(centrality) <- NULL

  communities <- tryCatch(igraph::cluster_louvain(graph), error = function(e) NULL)
  if (!is.null(communities)) {
    centrality$community <- igraph::membership(communities)[centrality$symbol]
    igraph::V(graph)$community <- as.integer(igraph::membership(communities)[igraph::V(graph)$name])
  } else {
    centrality$community <- NA_integer_
    igraph::V(graph)$community <- NA_integer_
  }

  list(
    graph = graph,
    hubs = centrality,
    edges = edges,
    gene_source = gene_source,
    interaction_source = interaction_source,
    score_cutoff = score_cutoff,
    max_genes = max_genes,
    control_group = result$control_group,
    treatment_group = result$treatment_group
  )
}

draw_ppi_graph <- function(graph, hubs, label_top_n = 15, main = NULL) {
  top_labels <- head(hubs$symbol, label_top_n)
  vertex_color <- ifelse(
    igraph::V(graph)$change == "up", "#dc2626",
    ifelse(igraph::V(graph)$change == "down", "#2563eb", "#9ca3af")
  )
  vertex_logfc <- abs(as.numeric(igraph::V(graph)$logFC))
  vertex_logfc[!is.finite(vertex_logfc)] <- 0
  vertex_size <- 5 + 2.5 * log1p(igraph::degree(graph)) + 2.2 * pmin(vertex_logfc, 4)
  labels <- ifelse(igraph::V(graph)$name %in% top_labels, igraph::V(graph)$name, NA)
  set.seed(2026)
  plot(
    graph,
    layout = igraph::layout_with_fr(graph),
    vertex.color = vertex_color,
    vertex.size = vertex_size,
    vertex.label = labels,
    vertex.label.cex = 0.75,
    vertex.label.color = "#111827",
    edge.width = 1 + 2 * igraph::E(graph)$score,
    edge.color = "#cbd5e1",
    main = main %||% "PPI network: color = direction, size = degree + |logFC|"
  )
  legend(
    "topleft",
    legend = c("up", "down", "other"),
    col = c("#dc2626", "#2563eb", "#9ca3af"),
    pch = 19,
    bty = "n"
  )
}

make_ppi_plot <- function(ppi_result, label_top_n = 15, mode = "network",
                          top_clusters = 6) {
  graph <- ppi_result$graph
  hubs <- ppi_result$hubs
  if (!identical(mode, "clusters")) {
    return(draw_ppi_graph(graph, hubs, label_top_n = label_top_n))
  }

  communities <- igraph::V(graph)$community
  if (all(is.na(communities))) {
    return(draw_ppi_graph(graph, hubs, label_top_n = label_top_n))
  }
  cluster_sizes <- sort(table(communities), decreasing = TRUE)
  selected <- names(head(cluster_sizes, top_clusters))
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(mfrow = c(2, 3), mar = c(1, 1, 4, 1))
  for (cluster_id in selected) {
    vertices <- igraph::V(graph)[communities == cluster_id]
    subgraph <- igraph::induced_subgraph(graph, vertices)
    sub_hubs <- hubs[hubs$symbol %in% igraph::V(subgraph)$name, , drop = FALSE]
    title <- paste0(
      "community ", cluster_id, "\n",
      "proteins: ", igraph::vcount(subgraph), "\n",
      "interactions: ", igraph::ecount(subgraph)
    )
    draw_ppi_graph(subgraph, sub_hubs, label_top_n = label_top_n, main = title)
  }
}

ppi_interpretation <- function(ppi_result) {
  hubs <- ppi_result$hubs
  top <- head(hubs, 10)
  up_n <- sum(hubs$change == "up", na.rm = TRUE)
  down_n <- sum(hubs$change == "down", na.rm = TRUE)
  hub_names <- paste(top$symbol, collapse = ", ")
  bias <- if (up_n > down_n) {
    paste0("网络节点以上调基因为主，更偏向 ", ppi_result$treatment_group, " 组相关的蛋白互作信号。")
  } else if (down_n > up_n) {
    paste0("网络节点以下调基因为主，更偏向 ", ppi_result$control_group, " 组相关的蛋白互作信号。")
  } else {
    "网络中上调和下调节点接近，说明该 PPI 网络可能同时包含两组方向的互作信号。"
  }
  paste0(
    "PPI 网络包含 ", igraph::vcount(ppi_result$graph), " 个节点、",
    igraph::ecount(ppi_result$graph), " 条边。", bias,
    " 当前 degree/betweenness 排名前列的 hub genes 是：", hub_names,
    "。PPI 解释的是差异基因编码蛋白之间的已知/预测互作关系，hub gene 更适合作为后续验证候选，而不是单独证明因果。"
  )
}

enrichment_symbols_by_direction <- function(result, direction) {
  deg <- result$deg
  deg <- deg[!is.na(deg$symbol) & nzchar(deg$symbol), , drop = FALSE]
  if (identical(direction, "up")) {
    return(unique(deg$symbol[deg$change == "up"]))
  }
  if (identical(direction, "down")) {
    return(unique(deg$symbol[deg$change == "down"]))
  }
  unique(deg$symbol[deg$change != "stable"])
}

map_symbols_to_entrez <- function(symbols) {
  symbols <- unique(symbols[!is.na(symbols) & nzchar(symbols)])
  if (length(symbols) == 0) {
    return(data.frame(SYMBOL = character(), ENTREZID = character()))
  }
  suppressMessages(clusterProfiler::bitr(
    symbols,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db::org.Hs.eg.db
  ))
}

run_single_enrichment <- function(entrez_ids, universe_ids, collection,
                                  p_cutoff, group_label) {
  entrez_ids <- unique(entrez_ids[!is.na(entrez_ids) & nzchar(entrez_ids)])
  universe_ids <- unique(universe_ids[!is.na(universe_ids) & nzchar(universe_ids)])
  if (length(entrez_ids) == 0) {
    return(data.frame())
  }

  if (identical(collection, "KEGG")) {
    enriched <- suppressMessages(clusterProfiler::enrichKEGG(
      gene = entrez_ids,
      universe = universe_ids,
      organism = "hsa",
      keyType = "ncbi-geneid",
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      qvalueCutoff = 1
    ))
    enriched <- tryCatch(
      suppressMessages(clusterProfiler::setReadable(
        enriched,
        OrgDb = org.Hs.eg.db::org.Hs.eg.db,
        keyType = "ENTREZID"
      )),
      error = function(e) enriched
    )
  } else {
    ont <- sub("^GO_", "", collection)
    enriched <- suppressMessages(clusterProfiler::enrichGO(
      gene = entrez_ids,
      universe = universe_ids,
      OrgDb = org.Hs.eg.db::org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = ont,
      pAdjustMethod = "BH",
      pvalueCutoff = 1,
      qvalueCutoff = 1,
      readable = TRUE
    ))
  }

  table <- as.data.frame(enriched)
  if (nrow(table) == 0) {
    return(table)
  }
  table$gene_group <- group_label
  table$collection <- ifelse(identical(collection, "KEGG"), "KEGG", gsub("_", " ", collection))
  table <- table[order(table$p.adjust, table$pvalue), , drop = FALSE]
  table
}

run_enrichment_analysis <- function(result, direction = "separate",
                                    collection = "GO_BP", p_cutoff = 0.05,
                                    min_genes = 5) {
  deg <- result$deg
  universe_symbols <- unique(deg$symbol[!is.na(deg$symbol) & nzchar(deg$symbol)])
  mapping <- map_symbols_to_entrez(universe_symbols)
  if (nrow(mapping) == 0) {
    stop("无法把基因 SYMBOL 转换为 ENTREZID。请确认是否为人类基因 SYMBOL，或提供注释表。")
  }
  universe_ids <- unique(mapping$ENTREZID)

  groups <- switch(
    direction,
    up = list("上调基因" = enrichment_symbols_by_direction(result, "up")),
    down = list("下调基因" = enrichment_symbols_by_direction(result, "down")),
    combined = list("上调+下调合并" = enrichment_symbols_by_direction(result, "combined")),
    separate = list(
      "上调基因" = enrichment_symbols_by_direction(result, "up"),
      "下调基因" = enrichment_symbols_by_direction(result, "down")
    ),
    list("上调+下调合并" = enrichment_symbols_by_direction(result, "combined"))
  )

  collections <- if (identical(collection, "GO_ALL")) {
    c("GO_BP", "GO_MF", "GO_CC")
  } else {
    collection
  }

  tables <- unlist(lapply(names(groups), function(label) {
    symbols <- unique(groups[[label]])
    entrez <- unique(mapping$ENTREZID[match(symbols, mapping$SYMBOL)])
    entrez <- entrez[!is.na(entrez)]
    if (length(entrez) < min_genes) {
      return(list(data.frame()))
    }
    lapply(collections, function(one_collection) {
      run_single_enrichment(
        entrez_ids = entrez,
        universe_ids = universe_ids,
        collection = one_collection,
        p_cutoff = p_cutoff,
        group_label = label
      )
    })
  }), recursive = FALSE)
  table <- do.call(rbind, Filter(function(x) nrow(x) > 0, tables))
  if (is.null(table)) {
    table <- data.frame()
  }
  rownames(table) <- NULL

  list(
    table = table,
    direction = direction,
    collection = collection,
    collection_label = ifelse(
      identical(collection, "KEGG"),
      "KEGG",
      ifelse(identical(collection, "GO_ALL"), "GO BP/MF/CC", gsub("_", " ", collection))
    ),
    p_cutoff = p_cutoff,
    min_genes = min_genes,
    control_group = result$control_group,
    treatment_group = result$treatment_group
  )
}

top_enrichment_rows <- function(enrichment_result, show_n = 15) {
  table <- enrichment_result$table
  if (nrow(table) == 0) {
    return(table)
  }
  table$.sig <- is.finite(table$p.adjust) & table$p.adjust <= enrichment_result$p_cutoff
  parts <- split(table, interaction(table$gene_group, table$collection, drop = TRUE))
  out <- do.call(rbind, lapply(parts, function(x) {
    x <- x[order(!x$.sig, x$p.adjust, x$pvalue), , drop = FALSE]
    head(x, show_n)
  }))
  out$.sig <- NULL
  rownames(out) <- NULL
  out
}

make_enrichment_plot <- function(enrichment_result, show_n = 15) {
  table <- top_enrichment_rows(enrichment_result, show_n)
  validate(need(nrow(table) > 0, "当前阈值和基因数下没有可绘制的富集条目。"))
  table$Description <- factor(table$Description, levels = rev(unique(table$Description)))
  plot <- ggplot(table, aes(x = Count, y = Description, color = p.adjust, size = Count)) +
    geom_point(alpha = 0.9) +
    scale_color_gradient(
      low = "#dc2626",
      high = "#2563eb",
      trans = "reverse",
      labels = format_p_decimal
    ) +
    labs(x = "Gene count", y = NULL, color = "padj", size = "genes") +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank())
  if (length(unique(table$gene_group)) > 1 || length(unique(table$collection)) > 1) {
    plot <- plot + facet_wrap(gene_group ~ collection, scales = "free_y")
  }
  plot
}

enrichment_interpretation <- function(enrichment_result) {
  table <- enrichment_result$table
  if (nrow(table) == 0) {
    return(paste0(
      "当前选择下没有得到富集条目。常见原因是显著差异基因太少、基因 ID 不是人类 SYMBOL，或者 KEGG/GO 能匹配到的基因不足。"
    ))
  }
  sig <- table[is.finite(table$p.adjust) & table$p.adjust <= enrichment_result$p_cutoff, , drop = FALSE]
  use_table <- if (nrow(sig) > 0) sig else head(table, 5)
  parts <- split(use_table, use_table$gene_group)
  txt <- paste(vapply(names(parts), function(group_name) {
    top_terms <- paste(head(parts[[group_name]]$Description, 3), collapse = "；")
    paste0(group_name, "主要富集在：", top_terms)
  }, character(1)), collapse = "。")

  if (nrow(sig) == 0) {
    paste0(
      "当前 padj <= ", enrichment_result$p_cutoff,
      " 下没有显著富集条目，但排名靠前的趋势项可作为探索参考。", txt,
      "。富集分析的意义是把一串基因翻译成更容易理解的功能或通路。"
    )
  } else {
    paste0(
      enrichment_result$collection_label, " 富集分析显示：", txt,
      "。这说明差异基因不是随机分散的，而是集中影响了这些功能或通路。"
    )
  }
}

format_gsea_gene_list <- function(x) {
  x <- trim_text(x)
  x[is.na(x)] <- ""
  out <- vapply(x, function(item) {
    genes <- unlist(strsplit(item, "/", fixed = TRUE))
    genes <- trim_text(genes)
    genes <- genes[!is.na(genes) & nzchar(genes)]
    paste(unique(genes), collapse = ", ")
  }, character(1))
  unname(out)
}

count_gsea_gene_list <- function(x) {
  formatted <- format_gsea_gene_list(x)
  out <- vapply(strsplit(formatted, ", ", fixed = TRUE), function(genes) {
    if (length(genes) == 1 && !nzchar(genes[1])) {
      return(0L)
    }
    length(genes)
  }, integer(1))
  unname(out)
}

run_gsea_analysis <- function(result, ontology = "BP", min_size = 10,
                              max_size = 500, p_cutoff = 0.25) {
  deg <- result$deg
  ranked <- deg[!is.na(deg$symbol) & nzchar(deg$symbol) & is.finite(deg$logFC), , drop = FALSE]
  ranked <- ranked[!duplicated(ranked$symbol), , drop = FALSE]
  if (nrow(ranked) < 100) {
    stop("GSEA 至少建议有 100 个可排序基因。")
  }

  mapping <- suppressMessages(clusterProfiler::bitr(
    ranked$symbol,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db::org.Hs.eg.db
  ))
  ranked$ENTREZID <- mapping$ENTREZID[match(ranked$symbol, mapping$SYMBOL)]
  ranked <- ranked[!is.na(ranked$ENTREZID), , drop = FALSE]
  ranked <- ranked[order(abs(ranked$logFC), decreasing = TRUE), , drop = FALSE]
  ranked <- ranked[!duplicated(ranked$ENTREZID), , drop = FALSE]
  gene_list <- ranked$logFC
  names(gene_list) <- ranked$ENTREZID
  gene_list <- sort(gene_list, decreasing = TRUE)
  if (length(gene_list) < 100) {
    stop("SYMBOL 转 ENTREZID 后可用于 GSEA 的基因少于 100 个。请检查物种或注释。")
  }

  if (identical(ontology, "KEGG")) {
    gsea <- suppressMessages(clusterProfiler::gseKEGG(
      geneList = gene_list,
      organism = "hsa",
      keyType = "ncbi-geneid",
      minGSSize = min_size,
      maxGSSize = max_size,
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      verbose = FALSE
    ))
    gsea <- tryCatch(
      suppressMessages(clusterProfiler::setReadable(
        gsea,
        OrgDb = org.Hs.eg.db::org.Hs.eg.db,
        keyType = "ENTREZID"
      )),
      error = function(e) gsea
    )
  } else {
    gsea <- suppressMessages(clusterProfiler::gseGO(
      geneList = gene_list,
      OrgDb = org.Hs.eg.db::org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = ontology,
      minGSSize = min_size,
      maxGSSize = max_size,
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      verbose = FALSE
    ))
    gsea <- suppressMessages(clusterProfiler::setReadable(
      gsea,
      OrgDb = org.Hs.eg.db::org.Hs.eg.db,
      keyType = "ENTREZID"
    ))
  }
  table <- as.data.frame(gsea)
  if (nrow(table) > 0) {
    if ("core_enrichment" %in% names(table)) {
      table$core_enrichment_genes <- format_gsea_gene_list(table$core_enrichment)
      table$leading_edge_genes <- table$core_enrichment_genes
      table$core_enrichment_gene_count <- count_gsea_gene_list(table$core_enrichment)
    }
    if ("leading_edge" %in% names(table)) {
      table$leading_edge_summary <- table$leading_edge
    }
    table$bias <- ifelse(table$NES > 0, result$treatment_group, result$control_group)
    table$direction <- ifelse(
      table$NES > 0,
      paste0("偏向 ", result$treatment_group),
      paste0("偏向 ", result$control_group)
    )
    priority_cols <- c(
      "ID", "Description", "setSize", "enrichmentScore", "NES",
      "pvalue", "p.adjust", "qvalue", "rank",
      "core_enrichment_gene_count", "core_enrichment_genes",
      "leading_edge_genes", "leading_edge_summary",
      "bias", "direction", "core_enrichment", "leading_edge"
    )
    table <- table[, c(priority_cols[priority_cols %in% names(table)],
                       setdiff(names(table), priority_cols)), drop = FALSE]
    table <- table[order(table$p.adjust, -abs(table$NES)), , drop = FALSE]
  }
  list(
    gsea = gsea,
    table = table,
    gene_list = gene_list,
    ontology = ontology,
    collection_label = ifelse(identical(ontology, "KEGG"), "KEGG", paste0("GO ", ontology)),
    p_cutoff = p_cutoff,
    control_group = result$control_group,
    treatment_group = result$treatment_group
  )
}

make_gsea_plot <- function(gsea_result, show_n = 15) {
  table <- gsea_result$table
  table <- table[is.finite(table$p.adjust) & table$p.adjust <= gsea_result$p_cutoff, , drop = FALSE]
  if (nrow(table) == 0) {
    table <- head(gsea_result$table, show_n)
  } else {
    table <- head(table, show_n)
  }
  validate(need(nrow(table) > 0, "当前 GSEA 没有可绘制的通路。"))
  table$Description <- factor(table$Description, levels = rev(table$Description))
  ggplot(table, aes(x = NES, y = Description, color = p.adjust, size = setSize)) +
    geom_point(alpha = 0.9) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "#6b7280") +
    scale_color_gradient(
      low = "#dc2626",
      high = "#2563eb",
      trans = "reverse",
      labels = format_p_decimal
    ) +
    labs(
      x = paste0("NES > 0: ", gsea_result$treatment_group, "; NES < 0: ", gsea_result$control_group),
      y = NULL,
      color = "padj",
      size = "genes"
    ) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank())
}

selected_gsea_running_rows <- function(gsea_result, show_n = 5) {
  table <- gsea_result$table
  table <- table[is.finite(table$p.adjust) & table$p.adjust <= gsea_result$p_cutoff, , drop = FALSE]
  if (nrow(table) == 0) {
    table <- head(gsea_result$table, show_n)
  } else {
    table <- head(table, show_n)
  }
  table
}

short_gene_text <- function(x, max_genes = 8) {
  genes <- unlist(strsplit(x %||% "", ", ", fixed = TRUE))
  genes <- genes[!is.na(genes) & nzchar(genes)]
  if (length(genes) == 0) {
    return("")
  }
  suffix <- if (length(genes) > max_genes) " ..." else ""
  paste0(paste(head(genes, max_genes), collapse = ", "), suffix)
}

gsea_running_gene_table <- function(gsea_result, show_n = 5) {
  table <- selected_gsea_running_rows(gsea_result, show_n = show_n)
  if (nrow(table) == 0) {
    return(data.frame())
  }
  out <- table[, intersect(
    c(
      "ID", "Description", "NES", "pvalue", "p.adjust",
      "core_enrichment_gene_count", "core_enrichment_genes",
      "leading_edge_genes", "bias", "direction"
    ),
    names(table)
  ), drop = FALSE]
  out
}

make_gsea_running_plot <- function(gsea_result, show_n = 5) {
  validate(need(requireNamespace("enrichplot", quietly = TRUE), "需要安装 enrichplot 才能绘制 GSEA 曲线图。"))
  table <- selected_gsea_running_rows(gsea_result, show_n = show_n)
  validate(need(nrow(table) > 0, "当前 GSEA 没有可绘制的通路。"))
  ids <- table$ID[seq_len(min(show_n, nrow(table)))]
  title_lines <- paste0(
    seq_along(ids), ". ",
    table$Description[seq_along(ids)],
    " | NES=", round(table$NES[seq_along(ids)], 3),
    " | padj=", format_p_decimal(table$p.adjust[seq_along(ids)], digits = 4)
  )
  enrichplot::gseaplot2(
    gsea_result$gsea,
    geneSetID = ids,
    title = paste(c(paste0(gsea_result$collection_label, " running enrichment"), title_lines), collapse = "\n"),
    pvalue_table = FALSE
  )
}

gsea_interpretation <- function(gsea_result) {
  table <- gsea_result$table
  collection <- gsea_result$collection_label %||% "GSEA"
  if (nrow(table) == 0) {
    return(paste0(
      "GSEA 没有返回可解释的 ", collection,
      " 条目。建议检查基因 ID 是否为人类 SYMBOL；如果运行 KEGG，还需要能访问 KEGG 在线数据库。"
    ))
  }
  sig <- table[is.finite(table$p.adjust) & table$p.adjust <= gsea_result$p_cutoff, , drop = FALSE]
  if (nrow(sig) == 0) {
    top <- head(table, 3)
    return(paste0(
      "当前阈值 padj <= ", gsea_result$p_cutoff,
      " 下没有显著 GSEA 条目。排名靠前的趋势项包括：",
      paste(top$Description, collapse = "；"),
      "。这说明通路层面信号可能存在但统计强度不足，建议作为探索性结果。"
    ))
  }
  pos <- sig[sig$NES > 0, , drop = FALSE]
  neg <- sig[sig$NES < 0, , drop = FALSE]
  pos_txt <- if (nrow(pos) > 0) paste(head(pos$Description, 3), collapse = "；") else "无明显条目"
  neg_txt <- if (nrow(neg) > 0) paste(head(neg$Description, 3), collapse = "；") else "无明显条目"
  paste0(
    collection, " GSEA 使用全部基因按 logFC 排序，而不是只看过阈值 DEG。NES > 0 表示通路整体偏向 ",
    gsea_result$treatment_group, "，NES < 0 表示偏向 ", gsea_result$control_group,
    "。当前显著条目中，", gsea_result$treatment_group, " 方向主要包括：",
    pos_txt, "；", gsea_result$control_group, " 方向主要包括：", neg_txt,
    "。它的好处是能捕捉一批基因整体轻微但一致的变化，适合解释疾病/处理组的主要生物过程。"
  )
}

final_analysis_recommendation <- function(result, wgcna_result = NULL,
                                          ppi_result = NULL, gsea_result = NULL) {
  change_counts <- table(result$deg$change)
  up <- as.integer(change_counts["up"] %||% 0)
  down <- as.integer(change_counts["down"] %||% 0)
  n_samples <- ncol(result$mat)
  base <- paste0(
    "最终建议：这个样本最适合采用“DEG 差异分析 + GSEA 通路解释”为主线，",
    "WGCNA 作为模块层面的辅助验证，PPI 用来筛选可后续实验验证的 hub genes。"
  )
  why <- paste0(
    "原因是当前数据为两组比较，共 ", n_samples, " 个样本，DEG 已能给出清晰上下调基因（上调 ",
    up, "，下调 ", down, "）。GSEA 不依赖硬阈值，能解释整套基因排序偏向的生物过程；",
    "WGCNA 对样本量更敏感，适合看模块是否整体偏向疾病/处理组；PPI 则依赖 STRING 等外部互作库，适合作为候选核心基因排序。"
  )
  extra <- character()
  if (!is.null(gsea_result)) {
    extra <- c(extra, gsea_interpretation(gsea_result))
  }
  if (!is.null(wgcna_result)) {
    extra <- c(extra, wgcna_interpretation(wgcna_result))
  }
  if (!is.null(ppi_result)) {
    extra <- c(extra, ppi_interpretation(ppi_result))
  }
  paste(c(base, why, extra), collapse = "\n\n")
}

ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      body { background: #f7f8fa; color: #111827; }
      .container-fluid { max-width: 1440px; }
      .app-title { margin: 18px 0 12px; }
      .well { background: #ffffff; border: 1px solid #e5e7eb; box-shadow: none; }
      .nav-tabs { margin-top: 8px; }
      .summary-strip {
        display: grid;
        grid-template-columns: repeat(auto-fit, minmax(130px, 1fr));
        gap: 10px;
        margin: 10px 0 14px;
      }
      .metric {
        background: #ffffff;
        border: 1px solid #e5e7eb;
        border-radius: 8px;
        padding: 10px 12px;
      }
      .metric-label { color: #6b7280; font-size: 12px; }
      .metric-value {
        color: #111827;
        font-size: 18px;
        font-weight: 650;
        line-height: 1.25;
        overflow-wrap: anywhere;
      }
      .hint { color: #6b7280; font-size: 13px; line-height: 1.45; }
      .helper-card {
        background: #f8fafc;
        border: 1px solid #e5e7eb;
        border-radius: 8px;
        padding: 10px 12px;
        margin: 10px 0 12px;
        line-height: 1.55;
      }
      .helper-card summary {
        cursor: pointer;
        font-weight: 650;
        color: #0f766e;
      }
      .template-actions {
        display: flex;
        flex-wrap: wrap;
        gap: 8px;
        margin: 8px 0 10px;
      }
      .empty-state {
        background: #ffffff;
        border: 1px solid #e5e7eb;
        border-radius: 8px;
        padding: 18px 20px;
        margin-top: 8px;
        line-height: 1.7;
      }
      .preview-section { margin-top: 18px; }
      .preview-section h4 { margin: 0 0 10px; line-height: 1.35; }
      .color-control { margin-top: 8px; }
      .tab-controls {
        background: #ffffff;
        border: 1px solid #e5e7eb;
        border-radius: 8px;
        padding: 12px 14px;
        margin: 12px 0 14px;
      }
      .tab-controls h4 { margin-top: 0; }
      .tab-controls .form-group { margin-bottom: 0; }
      .control-grid {
        display: grid;
        grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
        gap: 12px 16px;
        align-items: end;
      }
      .control-actions {
        display: flex;
        flex-wrap: wrap;
        gap: 10px;
        align-items: center;
        margin-top: 12px;
      }
      .analysis-note {
        background: #ffffff;
        border: 1px solid #e5e7eb;
        border-left: 4px solid #0f766e;
        border-radius: 8px;
        padding: 14px 16px;
        margin: 12px 0 16px;
        line-height: 1.65;
        white-space: pre-line;
      }
      .btn-primary { background-color: #0f766e; border-color: #0f766e; }
      .btn-primary:hover, .btn-primary:focus {
        background-color: #115e59; border-color: #115e59;
      }
      @media (max-width: 900px) {
        .summary-strip { grid-template-columns: repeat(2, minmax(120px, 1fr)); }
      }
    "))
  ),
  div(class = "app-title", titlePanel("BioInsight 一键生信分析平台")),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("1. 导入数据"),
      fileInput(
        "expr_file",
        "表达矩阵（行=基因/探针，列=样本）",
        accept = c(".csv", ".tsv", ".txt", ".xlsx", ".xls")
      ),
      actionButton("load_example", "载入当前示例数据", class = "btn-default"),
      tags$details(
        class = "helper-card",
        tags$summary("新手：没有表达矩阵/分组表？"),
        tags$p("最省事的做法：先下载模板看格式；如果是 GEO 数据，可以上传 series_matrix.txt.gz，软件会自动拆出表达矩阵和样本信息。"),
        div(
          class = "template-actions",
          downloadButton("download_expr_template", "表达矩阵模板"),
          downloadButton("download_group_template", "分组表模板")
        ),
        fileInput(
          "geo_series_file",
          "GEO Series Matrix（.txt 或 .txt.gz，可选）",
          accept = c(".txt", ".gz")
        ),
        actionButton("load_geo_series", "拆出 GEO 表达矩阵", class = "btn-default"),
        div(
          class = "template-actions",
          downloadButton("download_geo_expression", "下载拆出的表达矩阵"),
          downloadButton("download_geo_samples", "下载样本信息表")
        ),
        tags$p(
          class = "hint",
          "提示：GEO 页面通常在 Download family 或 Series Matrix File(s) 里下载 series_matrix.txt.gz。临床/分组信息常在样本 title、source 或 characteristics 字段里。"
        )
      ),
      uiOutput("expr_options"),
      fileInput(
        "annotation_file",
        "注释表（可选：第1列ID，第2列symbol）",
        accept = c(".csv", ".tsv", ".txt", ".xlsx", ".xls")
      ),
      tags$hr(),
      h4("2. 分组信息"),
      radioButtons(
        "group_method",
        "分组来源",
        choices = c(
          "分组文件" = "file",
          "样本名关键词" = "keyword",
          "手动粘贴" = "manual"
        ),
        selected = "file"
      ),
      conditionalPanel(
        "input.group_method == 'file'",
        fileInput(
          "group_file",
          "分组表（至少两列：sample, group）",
          accept = c(".csv", ".tsv", ".txt", ".xlsx", ".xls")
        ),
        uiOutput("group_file_options")
      ),
      conditionalPanel(
        "input.group_method == 'keyword'",
        textInput("control_keywords", "对照组关键词", "control|normal|ctrl"),
        textInput("treatment_keywords", "处理组关键词", "case|disease|treat|tumor")
      ),
      conditionalPanel(
        "input.group_method == 'manual'",
        uiOutput("manual_group_ui")
      ),
      uiOutput("group_select_ui")
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        tabPanel(
          "数据检查",
          br(),
          div(
            class = "analysis-note",
            h4("这一页做什么"),
            "先确认表达矩阵、分组表和数据类型有没有问题。软件会根据数值特征自动判断 raw counts、TPM/FPKM 或芯片/标准化表达矩阵；这一步决定后面差异分析用哪种统计方法。"
          ),
          div(
            class = "tab-controls",
            h4("数据类型与转换"),
            div(
              class = "control-grid",
              selectInput(
                "data_type",
                "数据类型",
                choices = c(
                  "自动识别并选择方法（推荐）" = "auto",
                  "芯片 / 已标准化表达矩阵（limma）" = "normalized",
                  "RNA-seq TPM / FPKM / RPKM（log2 后 limma）" = "rnaseq_tpm_fpkm",
                  "RNA-seq raw counts（DESeq2 / edgeR / voom）" = "rnaseq_counts"
                ),
                selected = "auto"
              ),
              conditionalPanel(
                "input.data_type == 'rnaseq_counts' || input.data_type == 'auto'",
                selectInput(
                  "count_method",
                  "raw count 分析方法",
                  choices = c(
                    "DESeq2（推荐）" = "deseq2",
                    "edgeR" = "edger",
                    "limma-voom" = "voom"
                  ),
                  selected = "deseq2"
                ),
                numericInput(
                  "count_min_total",
                  "raw count 低表达过滤：基因总 counts 至少",
                  value = 2,
                  min = 1,
                  step = 1
                )
              ),
              conditionalPanel(
                "input.data_type == 'normalized' || input.data_type == 'rnaseq_tpm_fpkm'",
                selectInput(
                  "log_mode",
                  "log2 转换",
                  choices = c(
                    "自动判断" = "auto",
                    "不转换" = "none",
                    "强制 log2(x + 1)" = "always"
                  ),
                  selected = "auto"
                )
              ),
              conditionalPanel(
                "input.data_type == 'normalized'",
                checkboxInput("normalize_between_arrays", "样本间分位数标准化", FALSE)
              )
            ),
            uiOutput("data_type_guide")
          ),
          uiOutput("data_check_ui")
        ),
        tabPanel(
          "分析结果",
          div(
            class = "analysis-note",
            h4("这一页做什么"),
            "差异分析会找出两组之间表达明显不同的基因，并给出 log2FC、P 值和校正后 P 值。后面的火山图、热图、富集分析、GSEA、WGCNA 和 PPI 都基于这里的结果。"
          ),
          div(
            class = "tab-controls",
            h4("差异分析参数"),
            div(
              class = "control-grid",
              selectInput(
                "p_column",
                "显著性列",
                choices = c("P.Value", "adj.P.Val"),
                selected = "adj.P.Val"
              ),
              numericInput("logfc_cutoff", "log2FC 阈值", value = 1, min = 0, step = 0.1),
              numericInput("p_cutoff", "P 值阈值", value = 0.05, min = 0, max = 1, step = 0.01)
            ),
            div(
              class = "control-actions",
              actionButton("run_analysis", "一键开始分析", class = "btn-primary")
            )
          ),
          uiOutput("final_recommendation"),
          uiOutput("summary_cards"),
          downloadButton("download_all", "下载完整差异表"),
          downloadButton("download_sig", "下载显著差异基因"),
          br(), br(),
          DT::DTOutput("deg_table")
        ),
        tabPanel(
          "火山图",
          div(
            class = "analysis-note",
            h4("这一页做什么"),
            "火山图把每个基因放到一张图里：横轴看变化大小，纵轴看显著程度。它能快速告诉你哪些基因升高、哪些降低，以及最值得优先关注的基因。"
          ),
          div(
            class = "tab-controls",
            h4("火山图参数"),
            div(
              class = "control-grid",
              checkboxInput("volcano_label_enabled", "标注显著基因", TRUE),
              checkboxInput("volcano_label_up", "标注升高基因", TRUE),
              numericInput("volcano_label_up_n", "升高标注 Top N", value = 10, min = 0, max = 100, step = 1),
              checkboxInput("volcano_label_down", "标注降低基因", TRUE),
              numericInput("volcano_label_down_n", "降低标注 Top N", value = 10, min = 0, max = 100, step = 1),
              colourpicker::colourInput("color_up", "上调颜色", "#dc2626"),
              colourpicker::colourInput("color_stable", "稳定颜色", "#9ca3af"),
              colourpicker::colourInput("color_down", "下调颜色", "#2563eb")
            )
          ),
          downloadButton("download_volcano", "下载火山图 PNG"),
          plotOutput("volcano_plot", height = "620px")
        ),
        tabPanel(
          "热图",
          div(
            class = "analysis-note",
            h4("这一页做什么"),
            "热图用颜色显示差异基因在每个样本里的表达模式。它能看出两组样本是否按表达模式分开，也能发现某些样本是不是和自己组内其他样本不太像。"
          ),
          div(
            class = "tab-controls",
            h4("热图参数"),
            div(
              class = "control-grid",
              numericInput("heatmap_top_n", "热图基因数", value = 50, min = 5, max = 300, step = 5),
              checkboxInput("heatmap_cluster_rows", "行聚类", TRUE),
              checkboxInput("heatmap_cluster_cols", "列聚类", TRUE),
              colourpicker::colourInput("heatmap_low_color", "低值颜色", "#2563eb"),
              colourpicker::colourInput("heatmap_mid_color", "中间颜色", "#ffffff"),
              colourpicker::colourInput("heatmap_high_color", "高值颜色", "#dc2626")
            )
          ),
          downloadButton("download_heatmap", "下载热图 PNG"),
          plotOutput("heatmap_plot", height = "720px")
        ),
        tabPanel(
          "箱线图",
          div(
            class = "analysis-note",
            h4("这一页做什么"),
            "箱线图看每个样本整体表达值的分布是否接近。正常情况下，同一批数据的箱体高度和中位数应该比较一致；如果某个样本明显偏高、偏低或分布很怪，可能是异常样本、批次效应或归一化问题。"
          ),
          div(
            class = "tab-controls",
            h4("箱线图参数"),
            div(
              class = "control-grid",
              numericInput("boxplot_max_genes", "最多抽样基因数", value = 10000, min = 500, max = 50000, step = 500),
              selectInput(
                "boxplot_scale",
                "箱线图尺度",
                choices = c(
                  "分析矩阵（差异分析实际使用）" = "analysis",
                  "原始输入矩阵" = "input",
                  "公司 FPKM 风格：log10(value + 0.001)" = "log10_input",
                  "按样本中位数中心化" = "median_centered"
                ),
                selected = "analysis"
              ),
              checkboxInput("boxplot_show_outliers", "显示离群点", TRUE)
            )
          ),
          downloadButton("download_boxplot", "下载箱线图 PNG"),
          plotOutput("boxplot_plot", height = "620px")
        ),
        tabPanel(
          "PCA",
          div(
            class = "analysis-note",
            h4("这一页做什么"),
            "PCA 看样本整体表达谱是否相似。点离得近说明整体表达更像；如果同组样本聚在一起、不同组分开，说明分组信号比较清楚；如果某个点远离所有样本，要重点检查。"
          ),
          div(
            class = "tab-controls",
            h4("PCA 参数"),
            div(
              class = "control-grid",
              checkboxInput("pca_show_ellipse", "显示分组椭圆", TRUE),
              checkboxInput("pca_show_centers", "显示分组中心点", TRUE)
            )
          ),
          downloadButton("download_pca", "下载 PCA PNG"),
          plotOutput("pca_plot", height = "620px")
        ),
        tabPanel(
          "富集分析",
          div(
            class = "analysis-note",
            h4("这一页做什么"),
            "富集分析只拿已经过阈值的显著差异基因来问：这些基因集中在哪些功能、通路或生物过程里。上调、下调分开做能看方向差异；合并做适合看总体受影响的功能主题。"
          ),
          div(
            class = "tab-controls",
            h4("富集分析参数"),
            div(
              class = "control-grid",
              selectInput(
                "enrich_direction",
                "分析哪些差异基因",
                choices = c(
                  "上调基因" = "up",
                  "下调基因" = "down",
                  "上调+下调合并" = "combined",
                  "上调和下调分别分析" = "separate"
                ),
                selected = "separate"
              ),
              selectInput(
                "enrich_collection",
                "富集类型",
                choices = c(
                  "GO BP" = "GO_BP",
                  "GO MF" = "GO_MF",
                  "GO CC" = "GO_CC",
                  "GO All（BP+MF+CC）" = "GO_ALL",
                  "KEGG" = "KEGG"
                ),
                selected = "GO_BP"
              ),
              numericInput("enrich_p_cutoff", "padj 阈值", value = 0.05, min = 0, max = 1, step = 0.01),
              numericInput("enrich_show_n", "图显示条目数", value = 15, min = 5, max = 50, step = 5),
              numericInput("enrich_min_genes", "最少基因数", value = 5, min = 2, max = 50, step = 1)
            ),
            div(
              class = "control-actions",
              actionButton("run_enrichment", "运行富集分析", class = "btn-default")
            )
          ),
          uiOutput("enrichment_interpretation"),
          downloadButton("download_enrichment_table", "下载富集结果表"),
          plotOutput("enrichment_plot", height = "720px"),
          h4("富集分析结果"),
          DT::DTOutput("enrichment_table")
        ),
        tabPanel(
          "GSEA",
          div(
            class = "analysis-note",
            h4("这一页做什么"),
            "GSEA 不只看显著差异基因，而是把全部基因按变化方向排序，判断某条通路是不是整体偏向某一组。它适合发现一组基因整体轻微但方向一致的变化。"
          ),
          div(
            class = "tab-controls",
            h4("GSEA 参数"),
            div(
              class = "control-grid",
              selectInput(
                "gsea_ontology",
                "基因集类型",
                choices = c("GO BP" = "BP", "GO MF" = "MF", "GO CC" = "CC", "KEGG" = "KEGG"),
                selected = "BP"
              ),
              numericInput("gsea_min_size", "最小基因集", value = 10, min = 5, max = 100, step = 5),
              numericInput("gsea_max_size", "最大基因集", value = 500, min = 100, max = 2000, step = 50),
              numericInput("gsea_p_cutoff", "padj 阈值", value = 0.25, min = 0, max = 1, step = 0.05),
              numericInput("gsea_show_n", "图显示条目数", value = 15, min = 5, max = 50, step = 5),
              numericInput("gsea_curve_n", "GSEA 曲线图通路数", value = 5, min = 1, max = 8, step = 1)
            ),
            div(
              class = "control-actions",
              actionButton("run_gsea", "运行 GSEA", class = "btn-default")
            )
          ),
          uiOutput("gsea_interpretation"),
          downloadButton("download_gsea_table", "下载 GSEA 结果表"),
          plotOutput("gsea_plot", height = "720px"),
          h4("GSEA running enrichment 曲线图"),
          plotOutput("gsea_running_plot", height = "760px"),
          h4("GSEA running enrichment 对应核心基因"),
          DT::DTOutput("gsea_running_gene_table"),
          h4("GSEA 结果"),
          DT::DTOutput("gsea_table"),
          div(
            class = "analysis-note",
            "`core_enrichment_genes` / `leading_edge_genes` 是真正推动该通路富集信号的核心基因，适合后续挑候选基因、做验证或放进 PPI 里继续看互作。"
          ),
          div(
            class = "analysis-note",
            h4("为什么使用 GSEA"),
            "GSEA 使用全部基因按表达变化排序，不只依赖显著差异基因阈值。很多通路不是由单个基因巨大变化驱动，而是一组基因整体轻中度同向变化；GSEA 能捕捉这种整体偏移。GO 更适合解释生物过程、分子功能和细胞组分，KEGG 更适合解释代谢通路和经典信号通路。"
          )
        ),
        tabPanel(
          "WGCNA",
          div(
            class = "analysis-note",
            h4("这一页做什么"),
            "WGCNA 会把表达模式相似的基因分成模块，再看哪个模块更接近分组差异。它不是找单个基因，而是找一群一起变化的基因模块，适合做机制线索。"
          ),
          div(
            class = "tab-controls",
            h4("WGCNA 参数"),
            div(
              class = "control-grid",
              numericInput("wgcna_top_n", "高变基因数", value = 5000, min = 500, max = 20000, step = 500),
              numericInput("wgcna_soft_power", "soft power（0=自动）", value = 0, min = 0, max = 30, step = 1),
              numericInput("wgcna_min_module_size", "最小模块基因数", value = 30, min = 10, max = 200, step = 5),
              numericInput("wgcna_merge_cut_height", "模块合并阈值", value = 0.25, min = 0.05, max = 0.5, step = 0.05)
            ),
            div(
              class = "control-actions",
              actionButton("run_wgcna", "运行 WGCNA", class = "btn-default")
            )
          ),
          uiOutput("wgcna_interpretation"),
          downloadButton("download_wgcna_modules", "下载 WGCNA 模块表"),
          downloadButton("download_wgcna_hubs", "下载 WGCNA hub genes"),
          plotOutput("wgcna_plot", height = "620px"),
          h4("模块-分组相关性"),
          DT::DTOutput("wgcna_module_table"),
          h4("模块 hub genes"),
          DT::DTOutput("wgcna_hub_table")
        ),
        tabPanel(
          "PPI",
          div(
            class = "analysis-note",
            h4("这一页做什么"),
            "PPI 把差异基因对应的蛋白放进互作网络里，看哪些蛋白连接最多、位置最关键。它适合从差异基因里筛选后续实验优先验证的 hub genes。"
          ),
          div(
            class = "tab-controls",
            h4("PPI 参数"),
            div(
              class = "control-grid",
              selectInput(
                "ppi_gene_source",
                "PPI 基因来源",
                choices = c(
                  "显著差异基因（默认）" = "significant",
                  "从当前表达矩阵中选择基因" = "selected"
                ),
                selected = "significant"
              ),
              selectInput(
                "ppi_interaction_source",
                "互作来源",
                choices = c(
                  "在线 STRING 自动获取（推荐，会发送基因名）" = "string_online",
                  "本地/上传 interaction 表" = "local"
                ),
                selected = "string_online"
              )
            ),
            conditionalPanel(
              "input.ppi_gene_source == 'selected'",
              selectizeInput(
                "ppi_selected_genes",
                "选择当前矩阵里的基因（可搜基因名或ID）",
                choices = NULL,
                multiple = TRUE,
                options = list(placeholder = "输入基因名或表达矩阵ID搜索，例如 LEPR 或 ENSG...")
              )
            ),
            fileInput(
              "ppi_file",
              "PPI 互作表（可选，默认读取 string_interactions.tsv）",
              accept = c(".tsv", ".txt", ".csv", ".xlsx", ".xls")
            ),
            div(
              class = "analysis-note",
              "PPI 图中红色表示上调，蓝色表示下调，灰色表示未显著或未在差异表中；节点越大，说明连接数和 |logFC| 综合更高。在线 STRING 会把基因名发送到 STRING 数据库查询互作；不想联网时请上传 STRING interaction 文件。"
            ),
            div(
              class = "control-grid",
              numericInput("ppi_score_cutoff", "combined_score 阈值", value = 0.7, min = 0, max = 1, step = 0.05),
              numericInput("ppi_max_genes", "最多差异基因数", value = 1000, min = 20, max = 2000, step = 20),
              numericInput("ppi_label_top_n", "标注 hub 数", value = 15, min = 0, max = 100, step = 1),
              selectInput(
                "ppi_plot_mode",
                "PPI 图模式",
                choices = c(
                  "整体网络" = "network",
                  "公司风格 Top 6 模块" = "clusters"
                ),
                selected = "network"
              )
            ),
            div(
              class = "control-actions",
              actionButton("run_ppi", "运行 PPI", class = "btn-default")
            )
          ),
          uiOutput("ppi_interpretation"),
          downloadButton("download_ppi_hubs", "下载 PPI hub genes"),
          downloadButton("download_ppi_edges", "下载 PPI edges"),
          plotOutput("ppi_plot", height = "720px"),
          h4("PPI hub genes"),
          DT::DTOutput("ppi_hub_table")
        )
      )
    )
  )
)

server <- function(input, output, session) {
  data_state <- reactiveVal(NULL)
  geo_state <- reactiveVal(NULL)

  observeEvent(input$expr_file, {
    df <- read_any_table(input$expr_file)
    geo_state(NULL)
    data_state(list(
      expr = df,
      groups = NULL,
      annotation = NULL,
      label = input$expr_file$name
    ))
  }, ignoreInit = TRUE)

  observeEvent(input$load_example, {
    if (!file.exists("step2output.Rdata")) {
      showNotification("当前目录没有 step2output.Rdata 示例文件。", type = "error")
      return()
    }
    env <- new.env(parent = emptyenv())
    load("step2output.Rdata", envir = env)
    expr_df <- data.frame(
      feature_id = rownames(env$exp),
      as.data.frame(env$exp, check.names = FALSE),
      check.names = FALSE
    )
    group_df <- data.frame(
      sample = colnames(env$exp),
      group = as.character(env$Group),
      stringsAsFactors = FALSE
    )
    annotation_df <- if (exists("ids", envir = env)) env$ids else NULL
    geo_state(NULL)
    data_state(list(
      expr = expr_df,
      groups = group_df,
      annotation = annotation_df,
      label = "step2output.Rdata"
    ))
    updateRadioButtons(session, "group_method", selected = "file")
  }, ignoreInit = TRUE)

  observeEvent(input$load_geo_series, {
    req(input$geo_series_file)
    parsed <- tryCatch(
      parse_geo_series_matrix(input$geo_series_file),
      error = function(e) {
        showNotification(e$message, type = "error", duration = 8)
        NULL
      }
    )
    if (is.null(parsed)) {
      return()
    }

    geo_state(parsed)
    data_state(list(
      expr = parsed$expr,
      groups = NULL,
      annotation = NULL,
      label = paste0(input$geo_series_file$name, " / GEO")
    ))
    updateRadioButtons(session, "group_method", selected = "manual")
    showNotification(
      paste0(
        "已拆出表达矩阵：", nrow(parsed$expr), " 个特征，",
        ncol(parsed$expr) - 1, " 个样本。请在数据检查页选择样本信息列生成分组，或手动填写分组。"
      ),
      type = "message",
      duration = 8
    )
  }, ignoreInit = TRUE)

  expr_df <- reactive({
    state <- data_state()
    validate(need(!is.null(state), "请先导入表达矩阵，或点击“载入当前示例数据”。"))
    state$expr
  })

  output$expr_options <- renderUI({
    df <- expr_df()
    selectInput(
      "expr_id_col",
      "基因/探针 ID 列",
      choices = names(df),
      selected = names(df)[1]
    )
  })

  sample_names <- reactive({
    df <- expr_df()
    req(input$expr_id_col)
    setdiff(names(df), input$expr_id_col)
  })

  output$data_check_ui <- renderUI({
    state <- data_state()
    geo <- geo_state()
    if (is.null(state)) {
      return(div(
        class = "empty-state",
        h4("请先导入数据"),
        div(
          class = "hint",
          "表达矩阵格式：第一列是基因或探针 ID，其余列是样本表达值。分组表格式：至少包含 sample 和 group 两列，样本名必须与表达矩阵列名一致。新手可以先在左侧下载模板，或上传 GEO series_matrix.txt.gz 自动拆出表达矩阵和样本信息。"
        )
      ))
    }

    tagList(
      div(
        class = "hint",
        "表达矩阵格式：第一列是基因或探针 ID，其余列是样本表达值。分组表格式：至少包含 sample 和 group 两列，样本名必须与表达矩阵列名一致。"
      ),
      if (!is.null(geo)) div(
        class = "preview-section",
        h4("GEO 样本信息"),
        div(
          class = "hint",
          "如果样本信息里已经有疾病、处理、组织、分型等字段，可以选择一列作为分组。生成后仍可在左侧手动修改。"
        ),
        uiOutput("geo_group_builder_ui"),
        DT::DTOutput("geo_sample_preview")
      ),
      div(
        class = "preview-section",
        h4("表达矩阵预览"),
        DT::DTOutput("expr_preview")
      ),
      div(
        class = "preview-section",
        h4("分组预览"),
        DT::DTOutput("group_preview")
      )
    )
  })

  output$geo_group_builder_ui <- renderUI({
    geo <- geo_state()
    req(geo)
    choices <- geo_group_candidate_cols(geo$samples)
    if (length(choices) == 0) {
      return(div(class = "hint", "这个 GEO 文件没有解析到可用的样本信息列，请改用手动分组。"))
    }
    div(
      class = "control-actions",
      selectInput(
        "geo_group_col",
        "用哪一列作为分组",
        choices = choices,
        selected = choices[1],
        width = "260px"
      ),
      actionButton("apply_geo_group", "用该列生成分组表", class = "btn-default")
    )
  })

  observeEvent(input$apply_geo_group, {
    geo <- geo_state()
    req(geo, input$geo_group_col)
    values <- trim_text(geo$samples[[input$geo_group_col]])
    group_df <- data.frame(
      sample = geo$samples$sample,
      group = values,
      stringsAsFactors = FALSE
    )
    group_df <- group_df[nzchar(group_df$group), , drop = FALSE]
    if (length(unique(group_df$group)) < 2) {
      showNotification("这一列不足两个分组，请换一列，或手动填写分组。", type = "warning", duration = 6)
      return()
    }
    state <- data_state()
    data_state(list(
      expr = state$expr,
      groups = group_df,
      annotation = state$annotation,
      label = state$label
    ))
    updateRadioButtons(session, "group_method", selected = "file")
    showNotification("已根据 GEO 样本信息生成分组表。", type = "message", duration = 5)
  }, ignoreInit = TRUE)

  output$geo_sample_preview <- DT::renderDT({
    geo <- geo_state()
    validate(need(!is.null(geo), "请先上传并拆出 GEO Series Matrix。"))
    DT::datatable(
      geo$samples,
      rownames = FALSE,
      filter = "top",
      options = list(pageLength = 10, scrollX = TRUE)
    )
  })

  output$expr_preview <- DT::renderDT({
    state <- data_state()
    if (is.null(state)) {
      return(DT::datatable(
        data.frame(提示 = "请先导入表达矩阵，或点击“载入当前示例数据”。"),
        rownames = FALSE,
        options = list(dom = "t", ordering = FALSE)
      ))
    }
    df <- expr_df()
    DT::datatable(
      head(df, 12),
      rownames = FALSE,
      options = list(pageLength = 12, scrollX = TRUE)
    )
  })

  group_file_df <- reactive({
    req(input$group_file)
    read_any_table(input$group_file)
  })

  output$group_file_options <- renderUI({
    req(input$group_file)
    df <- group_file_df()
    tagList(
      selectInput("group_sample_col", "样本名列", choices = names(df), selected = names(df)[1]),
      selectInput("group_group_col", "分组列", choices = names(df), selected = names(df)[min(2, ncol(df))])
    )
  })

  output$manual_group_ui <- renderUI({
    samples <- sample_names()
    textAreaInput(
      "manual_group_text",
      "手动分组",
      value = paste(c("sample,group", paste0(samples, ",")), collapse = "\n"),
      rows = min(max(length(samples) + 1, 6), 18),
      placeholder = "sample,group\nSample_1,Control\nSample_2,Treatment"
    )
  })

  group_table <- reactive({
    samples <- sample_names()
    state <- data_state()

    if (identical(input$group_method, "file")) {
      if (!is.null(input$group_file)) {
        raw <- group_file_df()
        out <- standardize_group_table(raw, input$group_sample_col, input$group_group_col)
      } else if (!is.null(state$groups)) {
        out <- standardize_group_table(state$groups, "sample", "group")
      } else {
        out <- data.frame(sample = character(), group = character())
      }
    } else if (identical(input$group_method, "keyword")) {
      control_hit <- grepl(input$control_keywords %||% "", samples, ignore.case = TRUE)
      treatment_hit <- grepl(input$treatment_keywords %||% "", samples, ignore.case = TRUE)
      group <- rep(NA_character_, length(samples))
      group[control_hit & !treatment_hit] <- "Control"
      group[treatment_hit & !control_hit] <- "Treatment"
      out <- data.frame(sample = samples, group = group, stringsAsFactors = FALSE)
      out <- out[!is.na(out$group), , drop = FALSE]
    } else {
      out <- parse_manual_groups(input$manual_group_text)
    }

    out <- out[out$sample %in% samples, , drop = FALSE]
    out[match(intersect(samples, out$sample), out$sample), , drop = FALSE]
  })

  output$group_preview <- DT::renderDT({
    state <- data_state()
    if (is.null(state)) {
      return(DT::datatable(
        data.frame(提示 = "请先导入表达矩阵，或点击“载入当前示例数据”。"),
        rownames = FALSE,
        options = list(dom = "t", ordering = FALSE)
      ))
    }
    df <- group_table()
    if (nrow(df) == 0) {
      return(DT::datatable(
        data.frame(提示 = "尚未识别到分组信息，请上传分组表、调整关键词，或手动粘贴 sample,group。"),
        rownames = FALSE,
        options = list(dom = "t", ordering = FALSE)
      ))
    }
    DT::datatable(
      df,
      rownames = FALSE,
      options = list(pageLength = 20, scrollX = TRUE)
    )
  })

  output$group_select_ui <- renderUI({
    df <- group_table()
    choices <- sort(unique(df$group))
    if (length(choices) < 2) {
      return(tags$div(
        class = "hint",
        "当前还没有识别到两个分组。请上传分组表、调整关键词，或手动粘贴 sample,group。"
      ))
    }
    tagList(
      selectInput("control_group", "对照组", choices = choices, selected = choices[1]),
      selectInput("treatment_group", "处理组", choices = choices, selected = choices[min(2, length(choices))])
    )
  })

  annotation_df <- reactive({
    state <- data_state()
    read_annotation(input$annotation_file, state$annotation)
  })

  output$data_type_guide <- renderUI({
    guess <- tryCatch(
      guess_expression_data_type(expr_df(), input$expr_id_col),
      error = function(e) NULL
    )
    data_type_guide_ui(guess)
  })

  analysis_result <- eventReactive(input$run_analysis, {
    validate(need(!identical(input$control_group, input$treatment_group),
                  "对照组和处理组不能相同。"))

    withProgress(message = "正在运行差异分析", value = 0, {
      incProgress(0.2, detail = "整理表达矩阵")
      method_label <- if (identical(input$data_type, "auto")) {
        "自动识别数据类型并选择模型"
      } else if (identical(input$data_type, "rnaseq_counts")) {
        paste0("拟合 ", toupper(input$count_method %||% "deseq2"), " 模型")
      } else {
        "拟合 limma 模型"
      }
      incProgress(0.45, detail = method_label)
      result <- run_analysis_by_data_type(
        expr_df = expr_df(),
        id_col = input$expr_id_col,
        data_type = input$data_type %||% "normalized",
        count_method = input$count_method %||% "deseq2",
        log_mode = input$log_mode %||% "auto",
        normalize_between_arrays = isTRUE(input$normalize_between_arrays),
        group_df = group_table(),
        control_group = input$control_group,
        treatment_group = input$treatment_group,
        p_cutoff = input$p_cutoff,
        logfc_cutoff = input$logfc_cutoff,
        p_column = input$p_column,
        annotation = annotation_df(),
        count_min_total = input$count_min_total %||% 2
      )
      incProgress(0.8, detail = "生成表格和图形")
      result
    })
  })

  observeEvent(analysis_result(), {
    updateSelectizeInput(
      session,
      "ppi_selected_genes",
      choices = ppi_gene_choices(analysis_result()),
      selected = character(),
      server = TRUE
    )
  }, ignoreInit = TRUE)

  volcano_colors <- reactive({
    c(
      down = input$color_down %||% "#2563eb",
      stable = input$color_stable %||% "#9ca3af",
      up = input$color_up %||% "#dc2626"
    )
  })

  volcano_label_settings <- reactive({
    if (!isTRUE(input$volcano_label_enabled)) {
      return(list(label_up = FALSE, label_down = FALSE, label_up_n = 0, label_down_n = 0))
    }
    list(
      label_up = isTRUE(input$volcano_label_up),
      label_down = isTRUE(input$volcano_label_down),
      label_up_n = as.integer(input$volcano_label_up_n %||% 0),
      label_down_n = as.integer(input$volcano_label_down_n %||% 0)
    )
  })

  heatmap_colors <- reactive({
    grDevices::colorRampPalette(c(
      input$heatmap_low_color %||% "#2563eb",
      input$heatmap_mid_color %||% "#ffffff",
      input$heatmap_high_color %||% "#dc2626"
    ))(100)
  })

  wgcna_result <- eventReactive(input$run_wgcna, {
    result <- analysis_result()
    withProgress(message = "正在运行 WGCNA", value = 0, {
      incProgress(0.2, detail = "筛选高变基因")
      out <- run_wgcna_analysis(
        result,
        top_n = input$wgcna_top_n,
        soft_power = input$wgcna_soft_power,
        min_module_size = input$wgcna_min_module_size,
        merge_cut_height = input$wgcna_merge_cut_height
      )
      incProgress(0.9, detail = "生成模块解释")
      out
    })
  })

  ppi_result <- eventReactive(input$run_ppi, {
    result <- analysis_result()
    withProgress(message = "正在运行 PPI", value = 0, {
      incProgress(0.3, detail = "读取互作表")
      ppi_table <- read_ppi_table(input$ppi_file)
      incProgress(0.6, detail = "构建网络")
      run_ppi_analysis(
        result,
        ppi_table,
        score_cutoff = input$ppi_score_cutoff,
        max_genes = input$ppi_max_genes,
        gene_source = input$ppi_gene_source %||% "significant",
        custom_gene_text = input$ppi_custom_genes %||% "",
        selected_genes = input$ppi_selected_genes %||% character(),
        interaction_source = input$ppi_interaction_source %||% "local"
      )
    })
  })

  enrichment_result <- eventReactive(input$run_enrichment, {
    result <- analysis_result()
    withProgress(message = "正在运行富集分析", value = 0, {
      incProgress(0.25, detail = "整理显著差异基因")
      out <- run_enrichment_analysis(
        result,
        direction = input$enrich_direction,
        collection = input$enrich_collection,
        p_cutoff = input$enrich_p_cutoff,
        min_genes = input$enrich_min_genes
      )
      incProgress(0.9, detail = "生成富集解释")
      out
    })
  })

  gsea_result <- eventReactive(input$run_gsea, {
    result <- analysis_result()
    withProgress(message = "正在运行 GSEA", value = 0, {
      incProgress(0.3, detail = "基因 ID 转换")
      out <- run_gsea_analysis(
        result,
        ontology = input$gsea_ontology,
        min_size = input$gsea_min_size,
        max_size = input$gsea_max_size,
        p_cutoff = input$gsea_p_cutoff
      )
      incProgress(0.9, detail = "生成 GSEA 解释")
      out
    })
  })

  output$final_recommendation <- renderUI({
    result <- analysis_result()
    w <- tryCatch(wgcna_result(), error = function(e) NULL)
    p <- tryCatch(ppi_result(), error = function(e) NULL)
    g <- tryCatch(gsea_result(), error = function(e) NULL)
    div(class = "analysis-note", final_analysis_recommendation(result, w, p, g))
  })

  output$summary_cards <- renderUI({
    result <- analysis_result()
    change_counts <- table(result$deg$change)
    up <- as.integer(change_counts["up"] %||% 0)
    down <- as.integer(change_counts["down"] %||% 0)
    stable <- as.integer(change_counts["stable"] %||% 0)
    div(
      class = "summary-strip",
      div(class = "metric", div(class = "metric-label", "数据类型"), div(class = "metric-value", result$data_type %||% "表达矩阵")),
      div(class = "metric", div(class = "metric-label", "方法"), div(class = "metric-value", result$analysis_method %||% "limma")),
      div(class = "metric", div(class = "metric-label", "基因/探针"), div(class = "metric-value", nrow(result$deg))),
      div(class = "metric", div(class = "metric-label", "样本"), div(class = "metric-value", ncol(result$mat))),
      div(class = "metric", div(class = "metric-label", "上调"), div(class = "metric-value", up)),
      div(class = "metric", div(class = "metric-label", "下调"), div(class = "metric-value", down)),
      div(class = "metric", div(class = "metric-label", "稳定"), div(class = "metric-value", stable))
    )
  })

  output$deg_table <- DT::renderDT({
    result <- analysis_result()
    DT::datatable(
      table_for_display(result$deg),
      rownames = FALSE,
      filter = "top",
      extensions = c("Buttons"),
      options = list(
        pageLength = 25,
        scrollX = TRUE,
        dom = "Bfrtip",
        buttons = c("copy", "csv")
      )
    )
  })

  output$volcano_plot <- renderPlot({
    labels <- volcano_label_settings()
    make_volcano(
      analysis_result(),
      colors = volcano_colors(),
      label_up = labels$label_up,
      label_down = labels$label_down,
      label_up_n = labels$label_up_n,
      label_down_n = labels$label_down_n
    )
  })

  output$pca_plot <- renderPlot({
    plot <- make_pca_plot(
      analysis_result(),
      show_ellipse = isTRUE(input$pca_show_ellipse),
      show_centers = isTRUE(input$pca_show_centers)
    )
    validate(need(!is.null(plot), "当前数据不足以绘制 PCA。"))
    plot
  })

  output$heatmap_plot <- renderPlot({
    result <- analysis_result()
    n <- heatmap_matrix(result, input$heatmap_top_n)
    validate(need(nrow(n) >= 2, "当前显著基因或高排序基因不足以绘制热图。"))
    annotation_col <- data.frame(Group = result$groups)
    rownames(annotation_col) <- colnames(n)
    pheatmap::pheatmap(
      n,
      show_colnames = FALSE,
      show_rownames = nrow(n) <= 80,
      cluster_rows = isTRUE(input$heatmap_cluster_rows),
      cluster_cols = isTRUE(input$heatmap_cluster_cols),
      scale = "row",
      annotation_col = annotation_col,
      color = heatmap_colors(),
      breaks = seq(-3, 3, length.out = 101),
      border_color = NA
    )
  })

  output$boxplot_plot <- renderPlot({
    make_boxplot_plot(
      analysis_result(),
      max_genes = input$boxplot_max_genes,
      show_outliers = isTRUE(input$boxplot_show_outliers),
      scale = input$boxplot_scale %||% "analysis"
    )
  })

  output$wgcna_interpretation <- renderUI({
    div(class = "analysis-note", wgcna_interpretation(wgcna_result()))
  })

  output$wgcna_plot <- renderPlot({
    make_wgcna_plot(wgcna_result())
  })

  output$wgcna_module_table <- DT::renderDT({
    DT::datatable(
      format_table_for_display(wgcna_result()$module_summary),
      rownames = FALSE,
      filter = "top",
      options = list(pageLength = 15, scrollX = TRUE)
    )
  })

  output$wgcna_hub_table <- DT::renderDT({
    DT::datatable(
      format_table_for_display(wgcna_result()$hub_genes),
      rownames = FALSE,
      filter = "top",
      options = list(pageLength = 20, scrollX = TRUE)
    )
  })

  output$ppi_interpretation <- renderUI({
    div(class = "analysis-note", ppi_interpretation(ppi_result()))
  })

  output$ppi_plot <- renderPlot({
    make_ppi_plot(
      ppi_result(),
      label_top_n = input$ppi_label_top_n,
      mode = input$ppi_plot_mode %||% "network"
    )
  })

  output$ppi_hub_table <- DT::renderDT({
    DT::datatable(
      format_table_for_display(ppi_result()$hubs),
      rownames = FALSE,
      filter = "top",
      options = list(pageLength = 20, scrollX = TRUE)
    )
  })

  output$enrichment_interpretation <- renderUI({
    div(class = "analysis-note", enrichment_interpretation(enrichment_result()))
  })

  output$enrichment_plot <- renderPlot({
    make_enrichment_plot(enrichment_result(), show_n = input$enrich_show_n)
  })

  output$enrichment_table <- DT::renderDT({
    DT::datatable(
      format_table_for_display(enrichment_result()$table),
      rownames = FALSE,
      filter = "top",
      options = list(pageLength = 20, scrollX = TRUE)
    )
  })

  output$gsea_interpretation <- renderUI({
    div(class = "analysis-note", gsea_interpretation(gsea_result()))
  })

  output$gsea_plot <- renderPlot({
    make_gsea_plot(gsea_result(), show_n = input$gsea_show_n)
  })

  output$gsea_running_plot <- renderPlot({
    make_gsea_running_plot(gsea_result(), show_n = input$gsea_curve_n)
  })

  output$gsea_running_gene_table <- DT::renderDT({
    DT::datatable(
      format_table_for_display(gsea_running_gene_table(gsea_result(), show_n = input$gsea_curve_n)),
      rownames = FALSE,
      filter = "top",
      options = list(pageLength = 10, scrollX = TRUE)
    )
  })

  output$gsea_table <- DT::renderDT({
    DT::datatable(
      format_table_for_display(gsea_result()$table),
      rownames = FALSE,
      filter = "top",
      options = list(pageLength = 20, scrollX = TRUE)
    )
  })

  output$download_expr_template <- downloadHandler(
    filename = function() "expression_matrix_template.csv",
    content = function(file) {
      write.csv(expression_template(), file, row.names = FALSE, fileEncoding = "UTF-8")
    }
  )

  output$download_group_template <- downloadHandler(
    filename = function() "group_table_template.csv",
    content = function(file) {
      write.csv(group_template(), file, row.names = FALSE, fileEncoding = "UTF-8")
    }
  )

  output$download_geo_expression <- downloadHandler(
    filename = function() paste0("GEO_expression_matrix_", Sys.Date(), ".csv"),
    content = function(file) {
      geo <- geo_state()
      validate(need(!is.null(geo), "请先上传并拆出 GEO Series Matrix。"))
      write.csv(geo$expr, file, row.names = FALSE, fileEncoding = "UTF-8")
    }
  )

  output$download_geo_samples <- downloadHandler(
    filename = function() paste0("GEO_sample_metadata_", Sys.Date(), ".csv"),
    content = function(file) {
      geo <- geo_state()
      validate(need(!is.null(geo), "请先上传并拆出 GEO Series Matrix。"))
      write.csv(geo$samples, file, row.names = FALSE, fileEncoding = "UTF-8")
    }
  )

  output$download_all <- downloadHandler(
    filename = function() paste0("DEG_all_", Sys.Date(), ".csv"),
    content = function(file) {
      write.csv(
        table_for_display(analysis_result()$deg),
        file,
        row.names = FALSE,
        fileEncoding = "UTF-8"
      )
    }
  )

  output$download_sig <- downloadHandler(
    filename = function() paste0("DEG_significant_", Sys.Date(), ".csv"),
    content = function(file) {
      deg <- table_for_display(analysis_result()$deg)
      write.csv(
        deg[deg$change != "stable", , drop = FALSE],
        file,
        row.names = FALSE,
        fileEncoding = "UTF-8"
      )
    }
  )

  output$download_volcano <- downloadHandler(
    filename = function() paste0("volcano_", Sys.Date(), ".png"),
    content = function(file) {
      labels <- volcano_label_settings()
      ggplot2::ggsave(
        file,
        make_volcano(
          analysis_result(),
          colors = volcano_colors(),
          label_up = labels$label_up,
          label_down = labels$label_down,
          label_up_n = labels$label_up_n,
          label_down_n = labels$label_down_n
        ),
        width = 8,
        height = 6,
        dpi = 300
      )
    }
  )

  output$download_pca <- downloadHandler(
    filename = function() paste0("pca_", Sys.Date(), ".png"),
    content = function(file) {
      plot <- make_pca_plot(
        analysis_result(),
        show_ellipse = isTRUE(input$pca_show_ellipse),
        show_centers = isTRUE(input$pca_show_centers)
      )
      validate(need(!is.null(plot), "当前数据不足以绘制 PCA。"))
      ggplot2::ggsave(file, plot, width = 8, height = 6, dpi = 300)
    }
  )

  output$download_boxplot <- downloadHandler(
    filename = function() paste0("boxplot_", Sys.Date(), ".png"),
    content = function(file) {
      ggplot2::ggsave(
        file,
        make_boxplot_plot(
          analysis_result(),
          max_genes = input$boxplot_max_genes,
          show_outliers = isTRUE(input$boxplot_show_outliers),
          scale = input$boxplot_scale %||% "analysis"
        ),
        width = 9,
        height = 6,
        dpi = 300
      )
    }
  )

  output$download_heatmap <- downloadHandler(
    filename = function() paste0("heatmap_", Sys.Date(), ".png"),
    content = function(file) {
      result <- analysis_result()
      n <- heatmap_matrix(result, input$heatmap_top_n)
      validate(need(nrow(n) >= 2, "当前显著基因或高排序基因不足以绘制热图。"))
      annotation_col <- data.frame(Group = result$groups)
      rownames(annotation_col) <- colnames(n)
      grDevices::png(file, width = 2200, height = 1800, res = 300)
      on.exit(grDevices::dev.off(), add = TRUE)
      pheatmap::pheatmap(
        n,
        show_colnames = FALSE,
        show_rownames = nrow(n) <= 80,
        cluster_rows = isTRUE(input$heatmap_cluster_rows),
        cluster_cols = isTRUE(input$heatmap_cluster_cols),
        scale = "row",
        annotation_col = annotation_col,
        color = heatmap_colors(),
        breaks = seq(-3, 3, length.out = 101),
        border_color = NA
      )
    }
  )

  output$download_enrichment_table <- downloadHandler(
    filename = function() paste0("ORA_enrichment_", Sys.Date(), ".csv"),
    content = function(file) {
      write.csv(enrichment_result()$table, file, row.names = FALSE, fileEncoding = "UTF-8")
    }
  )

  output$download_wgcna_modules <- downloadHandler(
    filename = function() paste0("WGCNA_modules_", Sys.Date(), ".csv"),
    content = function(file) {
      write.csv(wgcna_result()$module_summary, file, row.names = FALSE, fileEncoding = "UTF-8")
    }
  )

  output$download_wgcna_hubs <- downloadHandler(
    filename = function() paste0("WGCNA_hub_genes_", Sys.Date(), ".csv"),
    content = function(file) {
      write.csv(wgcna_result()$hub_genes, file, row.names = FALSE, fileEncoding = "UTF-8")
    }
  )

  output$download_ppi_hubs <- downloadHandler(
    filename = function() paste0("PPI_hub_genes_", Sys.Date(), ".csv"),
    content = function(file) {
      write.csv(ppi_result()$hubs, file, row.names = FALSE, fileEncoding = "UTF-8")
    }
  )

  output$download_ppi_edges <- downloadHandler(
    filename = function() paste0("PPI_edges_", Sys.Date(), ".csv"),
    content = function(file) {
      write.csv(ppi_result()$edges, file, row.names = FALSE, fileEncoding = "UTF-8")
    }
  )

  output$download_gsea_table <- downloadHandler(
    filename = function() paste0("GSEA_", Sys.Date(), ".csv"),
    content = function(file) {
      write.csv(gsea_result()$table, file, row.names = FALSE, fileEncoding = "UTF-8")
    }
  )
}

shinyApp(ui, server)
