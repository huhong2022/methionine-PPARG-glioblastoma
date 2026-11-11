#!/usr/bin/env Rscript

# U251 NC-versus-MR heatmap and volcano-plot visualization.
#
# Input-to-figure mapping:
#   U251_all_gene_DEG.tsv     -> Supplementary Fig. 13B
#   U251_transporter_DEG.tsv  -> Fig. 7A
#   U251_TRANSPORTERS_TPM.tsv -> Supplementary Fig. 13D
#
# Run from the repository root:
#   Rscript 05_U251_DEG_visualization.R
# or:
#   Rscript 05_U251_DEG_visualization.R /path/to/repository

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) args[[1]] else getwd()
project_root <- normalizePath(project_root, mustWork = TRUE)

data_dir <- file.path(project_root, "data")
gene_list_dir <- file.path(project_root, "gene_lists")
output_dir <- file.path(project_root, "results", "05_U251_DEG_visualization")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

transporter_expression_file <- file.path(data_dir, "U251_TRANSPORTERS_TPM.tsv")
all_gene_deg_file <- file.path(data_dir, "U251_all_gene_DEG.tsv")
transporter_deg_file <- file.path(data_dir, "U251_transporter_DEG.tsv")
matched_transporter_file <- file.path(
  gene_list_dir, "human_transporters_detected_in_U251.txt"
)

control_pattern <- "^N[0-9]+$"
treatment_pattern <- "^T[0-9]+$"
control_label <- "NC"
treatment_label <- "MR"
expected_replicates_per_group <- 4L
expected_detected_transporters <- 253L

log2fc_threshold <- 1
pvalue_threshold <- 0.05
label_genes <- c(
  "SLC25A1", "SLC2A4RG", "ABCD1", "SLC25A11", "KCNN4",
  "SLC30A1", "SLC19A2", "SLC6A9", "SLC7A5", "SLC3A2"
)

required_packages <- c("ggplot2", "ggrepel", "pheatmap", "RColorBrewer")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "))
}

read_tsv <- function(path) {
  if (!file.exists(path)) stop("Required input file not found: ", path)
  read.delim(
    path, header = TRUE, check.names = FALSE, quote = "",
    comment.char = "", stringsAsFactors = FALSE
  )
}

read_transporter_matrix <- function(path) {
  x <- read_tsv(path)
  if (ncol(x) < 3) stop("Transporter matrix must contain genes and >=2 samples.")
  genes <- trimws(as.character(x[[1]]))
  if (anyNA(genes) || any(!nzchar(genes))) stop("Transporter gene names are missing.")
  if (anyDuplicated(genes)) stop("Transporter matrix contains duplicated gene names.")

  matrix_data <- as.matrix(x[, -1, drop = FALSE])
  storage.mode(matrix_data) <- "numeric"
  if (anyNA(matrix_data)) stop("Transporter matrix contains NA/non-numeric values.")
  if (any(matrix_data < 0)) stop("TPM expression must not be negative.")
  rownames(matrix_data) <- genes
  matrix_data
}

assign_groups <- function(sample_names) {
  group <- ifelse(
    grepl(control_pattern, sample_names), control_label,
    ifelse(grepl(treatment_pattern, sample_names), treatment_label, NA_character_)
  )
  if (anyNA(group)) {
    stop(
      "Unable to assign NC/MR groups for samples: ",
      paste(sample_names[is.na(group)], collapse = ", ")
    )
  }
  group <- factor(group, levels = c(control_label, treatment_label))
  observed <- table(group)
  if (
    observed[[control_label]] != expected_replicates_per_group ||
      observed[[treatment_label]] != expected_replicates_per_group
  ) {
    stop(
      "Expected four NC and four MR biological replicates; observed ",
      observed[[control_label]], " and ", observed[[treatment_label]], "."
    )
  }
  group
}

find_column <- function(x, candidates, label, required = TRUE) {
  hit <- candidates[candidates %in% colnames(x)]
  if (length(hit) == 0) {
    if (required) stop("Unable to identify ", label, " column.")
    return(NA_character_)
  }
  hit[[1]]
}

format_volcano_table <- function(x, source_name) {
  gene_col <- find_column(x, c("Gene_Name", "Gene", "gene", "SYMBOL"), "gene")
  fc_col <- find_column(x, c("log2FC", "logFC", "log2FoldChange"), "log2 fold-change")
  p_col <- find_column(x, c("p_value", "P.Value", "pvalue", "p_value_raw"), "raw P-value")
  q_col <- find_column(
    x, c("q_value", "adj.P.Val", "padj", "p.adjust", "FDR"),
    "adjusted P-value", required = FALSE
  )
  group_col <- find_column(
    x, c("group", "Group", "Direction"), "original group", required = FALSE
  )

  out <- data.frame(
    Gene = trimws(as.character(x[[gene_col]])),
    log2FC = suppressWarnings(as.numeric(x[[fc_col]])),
    P_value = suppressWarnings(as.numeric(x[[p_col]])),
    Q_value = if (is.na(q_col)) NA_real_ else suppressWarnings(as.numeric(x[[q_col]])),
    Original_group = if (is.na(group_col)) NA_character_ else as.character(x[[group_col]]),
    stringsAsFactors = FALSE
  )
  invalid <- !nzchar(out$Gene) | !is.finite(out$log2FC) | !is.finite(out$P_value)
  if (any(invalid)) {
    stop(source_name, " contains ", sum(invalid), " rows with invalid gene/FC/P values.")
  }
  if (any(out$P_value < 0 | out$P_value > 1)) {
    stop(source_name, " contains P-values outside [0, 1].")
  }
  duplicate_genes <- unique(out$Gene[duplicated(out$Gene)])
  if (length(duplicate_genes) > 0) {
    warning(
      source_name, " contains duplicated gene symbols: ",
      paste(duplicate_genes, collapse = ", "),
      ". Rows are retained because this script only visualizes the supplied DEG results."
    )
  }

  out$P_value[out$P_value == 0] <- .Machine$double.xmin
  out$Direction <- "Not significant"
  out$Direction[
    out$P_value < pvalue_threshold & out$log2FC >= log2fc_threshold
  ] <- "Up"
  out$Direction[
    out$P_value < pvalue_threshold & out$log2FC <= -log2fc_threshold
  ] <- "Down"
  out$Direction <- factor(
    out$Direction,
    levels = c("Down", "Not significant", "Up")
  )
  out
}

make_volcano <- function(x, title, output_stem, genes_to_label = character(0)) {
  x$Label <- ifelse(x$Gene %in% genes_to_label, x$Gene, NA_character_)

  p <- ggplot2::ggplot(x, ggplot2::aes(x = log2FC, y = -log10(P_value))) +
    ggplot2::geom_point(
      ggplot2::aes(color = Direction), size = 2.3, alpha = 0.8
    ) +
    ggplot2::scale_color_manual(
      values = c(
        "Down" = "#00468B",
        "Not significant" = "grey72",
        "Up" = "#E64B35"
      ),
      drop = FALSE
    ) +
    ggplot2::geom_hline(
      yintercept = -log10(pvalue_threshold), linetype = 2, linewidth = 0.6
    ) +
    ggplot2::geom_vline(
      xintercept = c(-log2fc_threshold, log2fc_threshold),
      linetype = 2, linewidth = 0.6
    ) +
    ggrepel::geom_label_repel(
      ggplot2::aes(label = Label), na.rm = TRUE, size = 3.2,
      min.segment.length = 0, box.padding = 0.35,
      point.padding = 0.4, max.overlaps = Inf, show.legend = FALSE
    ) +
    ggplot2::labs(
      x = "Log2(fold change): MR versus NC",
      y = "-log10(P value)",
      color = NULL,
      title = title
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  ggplot2::ggsave(
    file.path(output_dir, paste0(output_stem, ".pdf")),
    p, width = 6.2, height = 5.2, units = "in", bg = "white"
  )
  ggplot2::ggsave(
    file.path(output_dir, paste0(output_stem, ".tiff")),
    p, width = 6.2, height = 5.2, units = "in", dpi = 600,
    compression = "lzw", bg = "white"
  )
  x$Label <- NULL
  x
}

# Supplementary Fig. 13D: expression heatmap for transporter candidates that
# were represented in the U251 RNA-seq dataset.
transporter_tpm <- read_transporter_matrix(transporter_expression_file)
if (nrow(transporter_tpm) != expected_detected_transporters) {
  stop(
    "Expected ", expected_detected_transporters,
    " transporter candidates with quantified entries; observed ",
    nrow(transporter_tpm), "."
  )
}

if (file.exists(matched_transporter_file)) {
  expected_genes <- unique(trimws(readLines(matched_transporter_file, warn = FALSE)))
  expected_genes <- expected_genes[nzchar(expected_genes)]
  if (!setequal(rownames(transporter_tpm), expected_genes)) {
    stop(
      "The transporter expression genes do not match ",
      basename(matched_transporter_file), "."
    )
  }
}

sample_group <- assign_groups(colnames(transporter_tpm))
sample_table <- data.frame(
  Sample = colnames(transporter_tpm),
  Group = as.character(sample_group),
  Biological_replicate = as.integer(sub("^[NT]", "", colnames(transporter_tpm))),
  stringsAsFactors = FALSE
)
write.table(
  sample_table,
  file.path(output_dir, "U251_sample_groups.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

transporter_log2 <- log2(transporter_tpm + 1)
annotation_col <- data.frame(Group = sample_group)
rownames(annotation_col) <- colnames(transporter_log2)
annotation_colors <- list(Group = c("NC" = "#3E6C8E", "MR" = "#A4473D"))
heatmap_height <- max(8, 0.055 * nrow(transporter_log2))

heatmap_plot <- pheatmap::pheatmap(
  transporter_log2,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  border_color = NA,
  angle_col = 45,
  fontsize_row = 4,
  fontsize_col = 10,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  color = grDevices::colorRampPalette(
    rev(RColorBrewer::brewer.pal(11, "RdBu"))
  )(100),
  main = "Human transporter candidates detected in U251 cells: NC versus MR",
  silent = TRUE
)

grDevices::pdf(
  file.path(output_dir, "SFig13D_transporter_heatmap.pdf"),
  width = 7.5, height = heatmap_height, useDingbats = FALSE
)
grid::grid.newpage()
grid::grid.draw(heatmap_plot$gtable)
grDevices::dev.off()

grDevices::tiff(
  file.path(output_dir, "SFig13D_transporter_heatmap.tiff"),
  width = 7.5, height = heatmap_height, units = "in", res = 600,
  compression = "lzw"
)
grid::grid.newpage()
grid::grid.draw(heatmap_plot$gtable)
grDevices::dev.off()

# Supplementary Fig. 13B: all-gene differential-expression volcano plot.
all_gene_deg <- format_volcano_table(
  read_tsv(all_gene_deg_file), basename(all_gene_deg_file)
)
all_gene_plot_data <- make_volcano(
  all_gene_deg,
  "All genes: MR versus NC",
  "SFig13B_all_gene_volcano"
)
write.table(
  all_gene_plot_data,
  file.path(output_dir, "SFig13B_all_gene_volcano_data.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Fig. 7A: volcano plot restricted to transporter candidates represented in
# the U251 RNA-seq dataset.
transporter_deg <- format_volcano_table(
  read_tsv(transporter_deg_file), basename(transporter_deg_file)
)
if (!setequal(transporter_deg$Gene, rownames(transporter_tpm))) {
  stop(
    basename(transporter_deg_file),
    " does not contain the same detected transporter-candidate set as the ",
    "transporter expression matrix."
  )
}
transporter_plot_data <- make_volcano(
  transporter_deg,
  "Differential expression of detected human transporter candidates",
  "Fig7A_transporter_volcano",
  genes_to_label = label_genes
)
write.table(
  transporter_plot_data,
  file.path(output_dir, "Fig7A_transporter_volcano_data.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

audit_table <- rbind(
  data.frame(
    Figure = "S Fig. 13B",
    Input = basename(all_gene_deg_file),
    Rows = nrow(all_gene_plot_data),
    Down = sum(all_gene_plot_data$Direction == "Down"),
    Not_significant = sum(all_gene_plot_data$Direction == "Not significant"),
    Up = sum(all_gene_plot_data$Direction == "Up")
  ),
  data.frame(
    Figure = "Fig. 7A",
    Input = basename(transporter_deg_file),
    Rows = nrow(transporter_plot_data),
    Down = sum(transporter_plot_data$Direction == "Down"),
    Not_significant = sum(transporter_plot_data$Direction == "Not significant"),
    Up = sum(transporter_plot_data$Direction == "Up")
  )
)
write.table(
  audit_table,
  file.path(output_dir, "volcano_plot_audit.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
message("U251 visualization completed. Results: ", output_dir)

