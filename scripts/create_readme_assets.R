args <- commandArgs(FALSE)
script_arg <- args[grepl("^--file=", args)]
if (length(script_arg) > 0) {
  script_path <- sub("^--file=", "", script_arg[1])
  script_dir <- dirname(normalizePath(script_path))
} else {
  script_dir <- getwd()
}
setwd(normalizePath(file.path(script_dir, "..")))
source("app.R")

dir.create("docs/assets", showWarnings = FALSE, recursive = TRUE)

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

ggplot2::ggsave(
  "docs/assets/result-volcano.png",
  make_volcano(result, label_top_n = 10),
  width = 8,
  height = 6,
  dpi = 180
)

ggplot2::ggsave(
  "docs/assets/result-pca.png",
  make_pca_plot(result, show_ellipse = TRUE, show_centers = TRUE),
  width = 8,
  height = 6,
  dpi = 180
)

heatmap_data <- heatmap_matrix(result, top_n = 50)
annotation_col <- data.frame(Group = result$groups)
rownames(annotation_col) <- colnames(heatmap_data)
grDevices::png("docs/assets/result-heatmap.png", width = 1500, height = 1200, res = 180)
pheatmap::pheatmap(
  heatmap_data,
  show_colnames = FALSE,
  show_rownames = nrow(heatmap_data) <= 80,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  scale = "row",
  annotation_col = annotation_col,
  color = grDevices::colorRampPalette(c("#2563eb", "#ffffff", "#dc2626"))(100),
  breaks = seq(-3, 3, length.out = 101),
  border_color = NA
)
grDevices::dev.off()

ppi <- run_ppi_analysis(result, read_ppi_table(), score_cutoff = 0.7, max_genes = 1000)
grDevices::png("docs/assets/result-ppi.png", width = 1500, height = 1200, res = 180)
make_ppi_plot(ppi, label_top_n = 10)
grDevices::dev.off()

gsea <- run_gsea_analysis(result, ontology = "BP", min_size = 10, max_size = 500, p_cutoff = 0.25)
ggplot2::ggsave(
  "docs/assets/result-gsea.png",
  make_gsea_plot(gsea, show_n = 12),
  width = 9,
  height = 7,
  dpi = 180
)

wgcna <- run_wgcna_analysis(
  result,
  top_n = 500,
  soft_power = 6,
  min_module_size = 20,
  merge_cut_height = 0.25
)
ggplot2::ggsave(
  "docs/assets/result-wgcna.png",
  make_wgcna_plot(wgcna),
  width = 8,
  height = 5,
  dpi = 180
)

cat("README assets generated in docs/assets\n")
