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
print(!is.null(make_pca_plot(result)))
