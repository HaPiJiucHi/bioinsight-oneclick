source("app.R")

if (!file.exists("step2output.Rdata")) {
  message("step2output.Rdata not found; skipping bundled example analysis.")
  quit(status = 0)
}

load("step2output.Rdata")

expr_df <- data.frame(
  feature_id = rownames(exp),
  as.data.frame(exp, check.names = FALSE),
  check.names = FALSE
)
group_df <- data.frame(
  sample = colnames(exp),
  group = as.character(Group),
  stringsAsFactors = FALSE
)

anno <- read_annotation(NULL, ids)
smart_anno <- read_annotation(
  NULL,
  data.frame(
    Gene = c("ENSG1", "ENSG2"),
    C1 = c(10, 20),
    C2 = c(12, 18),
    gene_name = c("TP53", "EGFR"),
    check.names = FALSE
  )
)
stopifnot(identical(smart_anno$symbol, c("TP53", "EGFR")))
stopifnot(identical(parse_gene_list("LEPR HSD3BP5, PDE5A\nMYLK;RGS18"), c("LEPR", "HSD3BP5", "PDE5A", "MYLK", "RGS18")))
stopifnot(identical(format_gsea_gene_list(c("TP53/EGFR/MYC", "")), c("TP53, EGFR, MYC", "")))
stopifnot(identical(count_gsea_gene_list(c("TP53/EGFR/MYC", "")), c(3L, 0L)))
prepared <- prepare_expression(expr_df, "feature_id", "none", FALSE)
result <- run_differential_analysis(
  prepared = prepared,
  group_df = group_df,
  control_group = "Normal",
  treatment_group = "Disease",
  p_cutoff = 0.05,
  logfc_cutoff = 1,
  p_column = "P.Value",
  annotation = anno
)

print(dim(result$deg))
print(table(result$deg$change))
print(head(table_for_display(result$deg), 3))
print(dim(heatmap_matrix(result, 50)))
print(!is.null(make_volcano(result)))
print(!is.null(make_volcano(result, label_up = TRUE, label_down = FALSE, label_up_n = 5, label_down_n = 0)))
print(!is.null(make_boxplot_plot(result, max_genes = 1000)))
print(!is.null(make_pca_plot(result)))

enrich <- run_enrichment_analysis(
  result,
  direction = "separate",
  collection = "GO_BP",
  p_cutoff = 0.05,
  min_genes = 5
)
print(dim(enrich$table))
if (nrow(enrich$table) > 0) {
  print(!is.null(make_enrichment_plot(enrich, show_n = 5)))
  print(substr(enrichment_interpretation(enrich), 1, 120))
}

set.seed(2026)
gene_n <- 200
sample_n <- 6
base_mu <- stats::rgamma(gene_n, shape = 2, rate = 0.1)
count_mat <- sapply(seq_len(sample_n), function(i) {
  fold <- if (i > 3) c(rep(4, 20), rep(1, gene_n - 20)) else rep(1, gene_n)
  stats::rnbinom(gene_n, mu = base_mu * fold, size = 10)
})
colnames(count_mat) <- paste0(rep(c("C", "T"), each = 3), seq_len(sample_n))
count_df <- data.frame(
  feature_id = paste0("Gene", seq_len(gene_n)),
  as.data.frame(count_mat, check.names = FALSE),
  check.names = FALSE
)
count_groups <- data.frame(
  sample = colnames(count_mat),
  group = rep(c("C", "T"), each = 3),
  stringsAsFactors = FALSE
)

low_count_df <- data.frame(
  feature_id = c("keep_total_2", "keep_total_3", "drop_total_1", "drop_total_0"),
  C1 = c(1, 1, 1, 0),
  C2 = c(1, 1, 0, 0),
  C3 = c(0, 1, 0, 0),
  T4 = c(0, 0, 0, 0),
  T5 = c(0, 0, 0, 0),
  T6 = c(0, 0, 0, 0),
  check.names = FALSE
)
filtered_counts <- prepare_count_expression(low_count_df, "feature_id")
stopifnot(identical(filtered_counts$feature_id, c("keep_total_2", "keep_total_3")))

for (method in c("deseq2", "edger", "voom")) {
  count_result <- run_analysis_by_data_type(
    expr_df = count_df,
    id_col = "feature_id",
    data_type = "rnaseq_counts",
    count_method = method,
    log_mode = "auto",
    normalize_between_arrays = FALSE,
    group_df = count_groups,
    control_group = "C",
    treatment_group = "T",
    p_cutoff = 0.05,
    logfc_cutoff = 1,
    p_column = "adj.P.Val",
    annotation = NULL
  )
  stopifnot(all(c("logFC", "P.Value", "adj.P.Val", "change") %in% names(count_result$deg)))
  stopifnot(ncol(count_result$mat) == 6)
  print(c(method = method, genes = nrow(count_result$deg), samples = ncol(count_result$mat)))
}

auto_count_result <- run_analysis_by_data_type(
  expr_df = count_df,
  id_col = "feature_id",
  data_type = "auto",
  count_method = "deseq2",
  log_mode = "auto",
  normalize_between_arrays = FALSE,
  group_df = count_groups,
  control_group = "C",
  treatment_group = "T",
  p_cutoff = 0.05,
  logfc_cutoff = 1,
  p_column = "adj.P.Val",
  annotation = NULL
)
stopifnot(grepl("自动识别", auto_count_result$data_type))
stopifnot(identical(auto_count_result$analysis_method, "DESeq2"))
print(c(auto = auto_count_result$data_type, method = auto_count_result$analysis_method))

if (file.exists("GSE7305_series_matrix.txt.gz")) {
  geo <- parse_geo_series_matrix("GSE7305_series_matrix.txt.gz")
  stopifnot(nrow(geo$expr) > 0)
  stopifnot(nrow(geo$samples) == ncol(geo$expr) - 1)
  print(dim(geo$expr))
  print(dim(geo$samples))
}
