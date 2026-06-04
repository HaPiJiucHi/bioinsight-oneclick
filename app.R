options(shiny.maxRequestSize = 500 * 1024^2)

required_packages <- c(
  "shiny", "limma", "ggplot2", "readr", "readxl", "DT", "pheatmap",
  "matrixStats", "ggrepel", "colourpicker"
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

make_pca_plot <- function(result) {
  pca_df <- make_pca_df(result)
  if (is.null(pca_df)) {
    return(NULL)
  }
  ggplot(pca_df, aes(PC1, PC2, color = group, label = sample)) +
    geom_point(size = 3, alpha = 0.85) +
    labs(
      x = pca_df$PC1_label[1],
      y = pca_df$PC2_label[1],
      color = "Group"
    ) +
    theme_minimal(base_size = 13) +
    theme(panel.grid.minor = element_blank(), legend.position = "top")
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
      .btn-primary { background-color: #0f766e; border-color: #0f766e; }
      .btn-primary:hover, .btn-primary:focus {
        background-color: #115e59; border-color: #115e59;
      }
      @media (max-width: 900px) {
        .summary-strip { grid-template-columns: repeat(2, minmax(120px, 1fr)); }
      }
    "))
  ),
  div(class = "app-title", titlePanel("差异分析一键软件")),
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
      uiOutput("expr_options"),
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
      uiOutput("group_select_ui"),
      tags$hr(),
      h4("3. 分析参数"),
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
      ),
      checkboxInput("volcano_label_enabled", "火山图标注显著基因", TRUE),
      numericInput("volcano_label_top_n", "标注前 N 个显著基因", value = 10, min = 0, max = 100, step = 1),
      numericInput("heatmap_top_n", "热图基因数", value = 50, min = 5, max = 300, step = 5),
      checkboxInput("heatmap_cluster_rows", "热图行聚类", TRUE),
      checkboxInput("heatmap_cluster_cols", "热图列聚类", TRUE),
      tags$details(
        tags$summary("颜色设置"),
        colourpicker::colourInput("color_up", "上调颜色", "#dc2626"),
        colourpicker::colourInput("color_stable", "稳定颜色", "#9ca3af"),
        colourpicker::colourInput("color_down", "下调颜色", "#2563eb"),
        colourpicker::colourInput("heatmap_low_color", "热图低值颜色", "#2563eb"),
        colourpicker::colourInput("heatmap_mid_color", "热图中间颜色", "#ffffff"),
        colourpicker::colourInput("heatmap_high_color", "热图高值颜色", "#dc2626")
      ),
      fileInput(
        "annotation_file",
        "注释表（可选：第1列ID，第2列symbol）",
        accept = c(".csv", ".tsv", ".txt", ".xlsx", ".xls")
      ),
      actionButton("run_analysis", "一键开始差异分析", class = "btn-primary")
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        tabPanel(
          "数据检查",
          br(),
          div(
            class = "empty-state",
            div(
              class = "hint",
              "表达矩阵格式：第一列是基因或探针 ID，其余列是样本表达值。分组表格式：至少包含 sample 和 group 两列，样本名必须与表达矩阵列名一致。"
            )
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
        ),
        tabPanel(
          "分析结果",
          uiOutput("summary_cards"),
          downloadButton("download_all", "下载完整差异表"),
          downloadButton("download_sig", "下载显著差异基因"),
          br(), br(),
          DT::DTOutput("deg_table")
        ),
        tabPanel(
          "火山图",
          br(),
          downloadButton("download_volcano", "下载火山图 PNG"),
          plotOutput("volcano_plot", height = "620px")
        ),
        tabPanel(
          "热图",
          br(),
          downloadButton("download_heatmap", "下载热图 PNG"),
          plotOutput("heatmap_plot", height = "720px")
        ),
        tabPanel(
          "PCA",
          br(),
          downloadButton("download_pca", "下载 PCA PNG"),
          plotOutput("pca_plot", height = "620px")
        )
      )
    )
  )
)

server <- function(input, output, session) {
  data_state <- reactiveVal(NULL)

  observeEvent(input$expr_file, {
    df <- read_any_table(input$expr_file)
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
    data_state(list(
      expr = expr_df,
      groups = group_df,
      annotation = annotation_df,
      label = "step2output.Rdata"
    ))
    updateRadioButtons(session, "group_method", selected = "file")
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
    if (is.null(state)) {
      return(div(
        class = "empty-state",
        h4("请先导入数据"),
        div(
          class = "hint",
          "表达矩阵格式：第一列是基因或探针 ID，其余列是样本表达值。分组表格式：至少包含 sample 和 group 两列，样本名必须与表达矩阵列名一致。"
        )
      ))
    }

    tagList(
      div(
        class = "hint",
        "表达矩阵格式：第一列是基因或探针 ID，其余列是样本表达值。分组表格式：至少包含 sample 和 group 两列，样本名必须与表达矩阵列名一致。"
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
    plot <- make_pca_plot(analysis_result())
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
      plot <- make_pca_plot(analysis_result())
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
}

shinyApp(ui, server)
