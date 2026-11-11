#!/usr/bin/env Rscript

# GSE84465 visualization using the published cell annotations and tSNE
# embedding. This script does not perform dimensionality reduction,
# clustering, marker discovery, or cell-type reannotation.
#
# Usage:
#   Rscript 03_scRNA_GSE84465_plots.R [project_root]

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) args[[1]] else getwd()
project_root <- normalizePath(project_root, mustWork = TRUE)

data_dir <- file.path(project_root, "data")
output_dir <- file.path(project_root, "results", "03_scRNA_GSE84465_plots")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

object_file <- file.path(data_dir, "GSE84465_annotated.rds")
cell_type_column <- "cell_type"
reduction_name <- "tsne"
feature_genes <- c("PPARG", "FARP1", "S100A16")

required_packages <- c(
  "Seurat", "SeuratObject", "ggplot2", "ggsci", "patchwork", "scales"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    "Missing R packages: ", paste(missing_packages, collapse = ", "),
    ". Install them before running this script."
  )
}

if (!file.exists(object_file)) {
  stop(
    "Required Seurat object not found: ", object_file,
    "\nThe object must contain the published GSE84465 cell annotations ",
    "and the original tSNE embedding."
  )
}

scRNA_GBM <- readRDS(object_file)
if (!inherits(scRNA_GBM, "Seurat")) {
  stop(basename(object_file), " is not a Seurat object.")
}
if (!cell_type_column %in% colnames(scRNA_GBM[[]])) {
  stop(
    "Metadata column '", cell_type_column, "' is absent from ",
    basename(object_file), "."
  )
}
if (!reduction_name %in% SeuratObject::Reductions(scRNA_GBM)) {
  stop(
    "Reduction '", reduction_name, "' is absent from ", basename(object_file),
    ". Supply the object containing the original GSE84465 tSNE coordinates; ",
    "this plotting script does not calculate a new embedding."
  )
}

if ("RNA" %in% SeuratObject::Assays(scRNA_GBM)) {
  SeuratObject::DefaultAssay(scRNA_GBM) <- "RNA"
}

missing_genes <- setdiff(feature_genes, rownames(scRNA_GBM))
if (length(missing_genes) > 0) {
  stop(
    "Genes absent from the active assay: ",
    paste(missing_genes, collapse = ", ")
  )
}

# Use the cell-type annotations stored in the supplied object without
# modification.
if (anyNA(scRNA_GBM[[cell_type_column, drop = TRUE]])) {
  stop("Metadata column '", cell_type_column, "' contains missing values.")
}
SeuratObject::Idents(scRNA_GBM) <- cell_type_column

cell_type_counts <- as.data.frame(
  table(scRNA_GBM[[cell_type_column, drop = TRUE]])
)
colnames(cell_type_counts) <- c("Cell_type", "Cells")
write.table(
  cell_type_counts,
  file.path(output_dir, "GSE84465_cell_type_counts.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

save_plot_pair <- function(plot, stem, width, height) {
  ggplot2::ggsave(
    file.path(output_dir, paste0(stem, ".pdf")),
    plot = plot, width = width, height = height, units = "in", bg = "white"
  )
  ggplot2::ggsave(
    file.path(output_dir, paste0(stem, ".tiff")),
    plot = plot, width = width, height = height, units = "in",
    dpi = 600, compression = "lzw", bg = "white"
  )
}

# Fig. 2E: cell-type map followed by PPARG, FARP1, and S100A16 expression.
p_dim <- Seurat::DimPlot(
  scRNA_GBM,
  reduction = reduction_name,
  group.by = cell_type_column,
  pt.size = 1
) +
  ggsci::scale_color_d3("category20c") +
  ggplot2::labs(title = "GSE84465 cell types") +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

p_features <- Seurat::FeaturePlot(
  scRNA_GBM,
  features = feature_genes,
  cols = c("#F0F0F0", "#CB181D"),
  reduction = reduction_name,
  order = TRUE,
  label = TRUE,
  ncol = length(feature_genes),
  pt.size = 1,
  combine = FALSE
)

fig2e <- patchwork::wrap_plots(
  c(list(p_dim), p_features),
  ncol = 4
)
save_plot_pair(
  fig2e,
  "Fig2E_GSE84465_DimPlot_FeaturePlot",
  width = 16,
  height = 4.8
)

# Fig. 2F: PPARG expression summarized by the published cell types.
fig2f <- suppressMessages(
  Seurat::DotPlot(
    scRNA_GBM,
    dot.scale = 5,
    features = "PPARG",
    group.by = cell_type_column,
    col.min = -0.5,
    col.max = 2.5,
    scale = TRUE
  ) +
    ggplot2::scale_color_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = 0,
      limits = c(-0.5, 2.5),
      oob = scales::squish,
      name = "Scaled average\nexpression"
    ) +
    ggplot2::labs(x = NULL, y = NULL, title = "PPARG expression") +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
)

save_plot_pair(
  fig2f,
  "Fig2F_GSE84465_PPARG_DotPlot",
  width = 5.5,
  height = 5.5
)

plot_manifest <- data.frame(
  Figure = c("Fig. 2E", "Fig. 2F"),
  Output_stem = c(
    "Fig2E_GSE84465_DimPlot_FeaturePlot",
    "Fig2F_GSE84465_PPARG_DotPlot"
  ),
  Genes = c(paste(feature_genes, collapse = ", "), "PPARG"),
  Reduction = c(reduction_name, "Not applicable"),
  Annotation = cell_type_column,
  stringsAsFactors = FALSE
)
write.table(
  plot_manifest,
  file.path(output_dir, "plot_manifest.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
message("GSE84465 plots completed. Results: ", output_dir)
