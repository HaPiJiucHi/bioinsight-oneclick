options(repos = c(CRAN = "https://mirrors.ustc.edu.cn/CRAN/"))
options(BioC_mirror = "https://mirrors.westlake.edu.cn/bioconductor")

cran_packages <- c(
  "shiny", "ggplot2", "readr", "readxl", "DT", "pheatmap", "matrixStats",
  "ggrepel", "colourpicker"
)
bioc_packages <- c("limma")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", ask = FALSE, update = FALSE)
}

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, ask = FALSE, update = FALSE)
  }
}

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

cat("依赖检查完成。\n")
