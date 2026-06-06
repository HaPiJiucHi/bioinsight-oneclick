options(shiny.maxRequestSize = 500 * 1024^2)

required_packages <- c(
  "shiny", "limma", "ggplot2", "readr", "readxl", "DT", "pheatmap",
  "matrixStats", "ggrepel", "colourpicker", "WGCNA", "igraph",
  "clusterProfiler", "org.Hs.eg.db", "enrichplot"
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
  path <- file$datapath %||% file
  name <- file$name %||% basename(path)
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
    feature_id = feature_id,
    transformed = transformed,
    normalized = isTRUE(normalize_between_arrays)
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

  out <- data.frame(
    feature_id = trim_text(df[[1]]),
    symbol = trim_text(df[[2]]),
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
    groups = groups,
    control_group = control_group,
    treatment_group = treatment_group,
    p_cutoff = p_cutoff,
    logfc_cutoff = logfc_cutoff,
    p_column = p_column,
    transformed = prepared$transformed,
    normalized = prepared$normalized
  )
}

table_for_display <- function(deg) {
  out <- deg
  out$row_id <- NULL
  out
}

make_volcano <- function(result, colors = NULL, label_top_n = 0) {
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

  if (label_top_n > 0 && requireNamespace("ggrepel", quietly = TRUE)) {
    label_df <- deg[deg$change != "stable" & is.finite(deg[[result$p_column]]), , drop = FALSE]
    label_df <- label_df[order(label_df[[result$p_column]], -abs(label_df$logFC)), , drop = FALSE]
    label_df <- head(label_df, label_top_n)
    if (nrow(label_df) > 0) {
      plot <- plot +
        ggrepel::geom_text_repel(
          data = label_df,
          aes(label = symbol),
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
    geom_text(aes(label = sprintf("r=%.2f\np=%.2g", trait_correlation, p_value)), size = 3.4) +
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
    "，p=", signif(top$p_value, 3), "，", direction, "。", deg_bias,
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

make_ppi_plot <- function(ppi_result, label_top_n = 15) {
  graph <- ppi_result$graph
  hubs <- ppi_result$hubs
  top_labels <- head(hubs$symbol, label_top_n)
  vertex_color <- ifelse(
    igraph::V(graph)$change == "up", "#dc2626",
    ifelse(igraph::V(graph)$change == "down", "#2563eb", "#9ca3af")
  )
  vertex_size <- 5 + 3 * log1p(igraph::degree(graph))
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
    main = "PPI network"
  )
  legend(
    "topleft",
    legend = c("up", "down", "other"),
    col = c("#dc2626", "#2563eb", "#9ca3af"),
    pch = 19,
    bty = "n"
  )
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
  table <- as.data.frame(gsea)
  if (nrow(table) > 0) {
    table$bias <- ifelse(table$NES > 0, result$treatment_group, result$control_group)
    table$direction <- ifelse(
      table$NES > 0,
      paste0("偏向 ", result$treatment_group),
      paste0("偏向 ", result$control_group)
    )
    table <- table[order(table$p.adjust, -abs(table$NES)), , drop = FALSE]
  }
  list(
    gsea = gsea,
    table = table,
    gene_list = gene_list,
    ontology = ontology,
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
    scale_color_gradient(low = "#dc2626", high = "#2563eb", trans = "reverse") +
    labs(
      x = paste0("NES > 0: ", gsea_result$treatment_group, "; NES < 0: ", gsea_result$control_group),
      y = NULL,
      color = "padj",
      size = "genes"
    ) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank())
}

gsea_interpretation <- function(gsea_result) {
  table <- gsea_result$table
  if (nrow(table) == 0) {
    return("GSEA 没有返回可解释的 GO 条目。建议检查基因 ID 是否为人类 SYMBOL，或改用自定义基因集。")
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
    "GSEA 使用全部基因按 logFC 排序，而不是只看过阈值 DEG。NES > 0 表示通路整体偏向 ",
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
        grid-template-columns: repeat(5, minmax(120px, 1fr));
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
      .metric-value { color: #111827; font-size: 22px; font-weight: 650; }
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
          uiOutput("data_check_ui")
        ),
        tabPanel(
          "分析结果",
          div(
            class = "tab-controls",
            h4("差异分析参数"),
            div(
              class = "control-grid",
              selectInput(
                "log_mode",
                "log2 转换",
                choices = c(
                  "自动判断" = "auto",
                  "不转换" = "none",
                  "强制 log2(x + 1)" = "always"
                ),
                selected = "auto"
              ),
              checkboxInput("normalize_between_arrays", "样本间分位数标准化", FALSE),
              numericInput("logfc_cutoff", "log2FC 阈值", value = 1, min = 0, step = 0.1),
              numericInput("p_cutoff", "P 值阈值", value = 0.05, min = 0, max = 1, step = 0.01),
              selectInput(
                "p_column",
                "显著性列",
                choices = c("P.Value", "adj.P.Val"),
                selected = "P.Value"
              )
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
            class = "tab-controls",
            h4("火山图参数"),
            div(
              class = "control-grid",
              checkboxInput("volcano_label_enabled", "标注显著基因", TRUE),
              numericInput("volcano_label_top_n", "标注前 N 个基因", value = 10, min = 0, max = 100, step = 1),
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
          "PCA",
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
          "GSEA",
          div(
            class = "tab-controls",
            h4("GSEA 参数"),
            div(
              class = "control-grid",
              selectInput("gsea_ontology", "GO 类型", choices = c("BP", "MF", "CC"), selected = "BP"),
              numericInput("gsea_min_size", "最小基因集", value = 10, min = 5, max = 100, step = 5),
              numericInput("gsea_max_size", "最大基因集", value = 500, min = 100, max = 2000, step = 50),
              numericInput("gsea_p_cutoff", "padj 阈值", value = 0.25, min = 0, max = 1, step = 0.05),
              numericInput("gsea_show_n", "图显示条目数", value = 15, min = 5, max = 50, step = 5)
            ),
            div(
              class = "control-actions",
              actionButton("run_gsea", "运行 GSEA", class = "btn-default")
            )
          ),
          uiOutput("gsea_interpretation"),
          downloadButton("download_gsea_table", "下载 GSEA 结果表"),
          plotOutput("gsea_plot", height = "720px"),
          h4("GSEA 结果"),
          DT::DTOutput("gsea_table"),
          div(
            class = "analysis-note",
            h4("为什么使用 GSEA"),
            "GSEA 使用全部基因按表达变化排序，不只依赖显著差异基因阈值。很多通路不是由单个基因巨大变化驱动，而是一组基因整体轻中度同向变化；GSEA 能捕捉这种整体偏移，因此适合放在差异分析之后解释样本更偏向哪些生物过程。"
          )
        ),
        tabPanel(
          "WGCNA",
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
            class = "tab-controls",
            h4("PPI 参数"),
            fileInput(
              "ppi_file",
              "PPI 互作表（可选，默认读取 string_interactions.tsv）",
              accept = c(".tsv", ".txt", ".csv", ".xlsx", ".xls")
            ),
            div(
              class = "control-grid",
              numericInput("ppi_score_cutoff", "combined_score 阈值", value = 0.7, min = 0, max = 1, step = 0.05),
              numericInput("ppi_max_genes", "最多差异基因数", value = 1000, min = 20, max = 2000, step = 20),
              numericInput("ppi_label_top_n", "标注 hub 数", value = 15, min = 0, max = 100, step = 1)
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

  analysis_result <- eventReactive(input$run_analysis, {
    validate(need(!identical(input$control_group, input$treatment_group),
                  "对照组和处理组不能相同。"))

    withProgress(message = "正在运行差异分析", value = 0, {
      incProgress(0.2, detail = "整理表达矩阵")
      prepared <- prepare_expression(
        expr_df(),
        input$expr_id_col,
        input$log_mode,
        input$normalize_between_arrays
      )
      incProgress(0.45, detail = "拟合 limma 模型")
      result <- run_differential_analysis(
        prepared = prepared,
        group_df = group_table(),
        control_group = input$control_group,
        treatment_group = input$treatment_group,
        p_cutoff = input$p_cutoff,
        logfc_cutoff = input$logfc_cutoff,
        p_column = input$p_column,
        annotation = annotation_df()
      )
      incProgress(0.8, detail = "生成表格和图形")
      result
    })
  })

  volcano_colors <- reactive({
    c(
      down = input$color_down %||% "#2563eb",
      stable = input$color_stable %||% "#9ca3af",
      up = input$color_up %||% "#dc2626"
    )
  })

  volcano_label_top_n <- reactive({
    if (!isTRUE(input$volcano_label_enabled)) {
      return(0)
    }
    as.integer(input$volcano_label_top_n %||% 0)
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
        max_genes = input$ppi_max_genes
      )
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
    make_volcano(
      analysis_result(),
      colors = volcano_colors(),
      label_top_n = volcano_label_top_n()
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

  output$wgcna_interpretation <- renderUI({
    div(class = "analysis-note", wgcna_interpretation(wgcna_result()))
  })

  output$wgcna_plot <- renderPlot({
    make_wgcna_plot(wgcna_result())
  })

  output$wgcna_module_table <- DT::renderDT({
    DT::datatable(
      wgcna_result()$module_summary,
      rownames = FALSE,
      filter = "top",
      options = list(pageLength = 15, scrollX = TRUE)
    )
  })

  output$wgcna_hub_table <- DT::renderDT({
    DT::datatable(
      wgcna_result()$hub_genes,
      rownames = FALSE,
      filter = "top",
      options = list(pageLength = 20, scrollX = TRUE)
    )
  })

  output$ppi_interpretation <- renderUI({
    div(class = "analysis-note", ppi_interpretation(ppi_result()))
  })

  output$ppi_plot <- renderPlot({
    make_ppi_plot(ppi_result(), label_top_n = input$ppi_label_top_n)
  })

  output$ppi_hub_table <- DT::renderDT({
    DT::datatable(
      ppi_result()$hubs,
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

  output$gsea_table <- DT::renderDT({
    DT::datatable(
      gsea_result()$table,
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
      ggplot2::ggsave(
        file,
        make_volcano(
          analysis_result(),
          colors = volcano_colors(),
          label_top_n = volcano_label_top_n()
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
    filename = function() paste0("GSEA_GO_", Sys.Date(), ".csv"),
    content = function(file) {
      write.csv(gsea_result()$table, file, row.names = FALSE, fileEncoding = "UTF-8")
    }
  )
}

shinyApp(ui, server)
