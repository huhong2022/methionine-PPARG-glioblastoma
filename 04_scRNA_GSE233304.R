#!/usr/bin/env Rscript

# GSE233304 single-cell RNA-seq analysis
# Three matched patients (tumor CD45+ and PBMC CD45+; six libraries).
# The workflow performs global preprocessing, marker-module annotation,
# split UMAP visualization, and patient-level paired PPARG comparisons.

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1L) args[[1]] else getwd()
run_mode <- if (length(args) >= 2L) tolower(args[[2]]) else "preprocess"
if (!run_mode %in% c("preprocess", "reuse")) {
  stop("The second argument must be 'preprocess' or 'reuse'.")
}

project_root <- normalizePath(project_root, winslash = "/", mustWork = FALSE)
data_dir <- file.path(project_root, "data", "GSE233304")
output_dir <- file.path(project_root, "results", "04_scRNA_GSE233304")
object_dir <- file.path(output_dir, "objects")
annotation_dir <- file.path(output_dir, "annotation")
plot_dir <- file.path(output_dir, "plots")
statistics_dir <- file.path(output_dir, "statistics")

for (path in c(output_dir, object_dir, annotation_dir, plot_dir, statistics_dir)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

required_packages <- c(
  "Seurat", "SeuratObject", "harmony", "ggplot2", "ggsci", "patchwork",
  "dplyr", "tidyr", "ggpubr"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]
if (length(missing_packages)) {
  stop("Install the following packages before running this script: ",
       paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(ggsci)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(ggpubr)
})
set.seed(42)
options(future.globals.maxSize = 20 * 1024^3)

# User-adjustable settings ---------------------------------------------------
run_preprocessing <- identical(run_mode, "preprocess")
preprocessed_rds <- file.path(object_dir, "GSE233304_preprocessed.rds")
final_rds <- file.path(object_dir, "GSE233304_final.rds")
flip_umap1 <- TRUE
point_size <- 0.3

# Optional manual corrections after marker review: cluster ID -> label.
# Example: manual_annotation_corrections <- c("18" = "Plasma_cell")
manual_annotation_corrections <- c()

# 1. Preprocess the six libraries --------------------------------------------
if (run_preprocessing) {
  if (!dir.exists(data_dir)) {
    stop("Data directory not found: ", data_dir)
  }
  samples <- data.frame(
    file = c(
      "GSM7421861_S2_filtered_feature_bc_matrix.h5",
      "GSM7421862_S3_filtered_feature_bc_matrix.h5",
      "GSM7421864_S5_filtered_feature_bc_matrix.h5",
      "GSM7421865_S6_filtered_feature_bc_matrix.h5",
      "GSM7421869_S10_filtered_feature_bc_matrix.h5",
      "GSM7421870_S11_filtered_feature_bc_matrix.h5"
    ),
    sample_id = c("S2", "S3", "S5", "S6", "S10", "S11"),
    orig.ident = rep(c("Patient_1", "Patient_2", "Patient_9"), each = 2),
    sample = rep(c("Tumor_CD45", "PBMCs_CD45"), 3),
    stringsAsFactors = FALSE
  )
  paths <- file.path(data_dir, samples$file)
  if (!all(file.exists(paths))) {
    stop("H5 files not found:\n", paste(samples$file[!file.exists(paths)], collapse = "\n"))
  }

  seurat_list <- lapply(seq_len(nrow(samples)), function(i) {
    counts <- Seurat::Read10X_h5(paths[i], use.names = TRUE, unique.features = TRUE)
    if (is.list(counts)) counts <- counts[["Gene Expression"]]
    if (is.null(counts)) stop("Missing Gene Expression matrix in ", samples$file[i])
    obj <- CreateSeuratObject(counts, project = samples$orig.ident[i],
                              min.cells = 3, min.features = 200)
    obj <- RenameCells(obj, add.cell.id = paste(samples$orig.ident[i], samples$sample[i], sep = "_"))
    obj$orig.ident <- samples$orig.ident[i]
    obj$sample <- samples$sample[i]
    obj$sample_id <- samples$sample_id[i]
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
    obj
  })
  names(seurat_list) <- samples$sample_id
  qc <- data.frame(samples[, c("sample_id", "orig.ident", "sample")],
                   n_cells = vapply(seurat_list, ncol, numeric(1)))
  write.table(qc, file.path(statistics_dir, "sample_cell_counts.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  scRNA_CD45 <- merge(seurat_list[[1]], y = seurat_list[-1], project = "GSE233304_CD45")
  rm(seurat_list)
  invisible(gc())
  message("Cells after CreateSeuratObject (no additional QC): ", ncol(scRNA_CD45))

  # Merge Seurat v5 RNA layers before normalization.
  if (inherits(scRNA_CD45[["RNA"]], "Assay5")) {
    layer_names <- Layers(scRNA_CD45[["RNA"]])
    layer_base <- sub("\\..*$", "", layer_names)
    if (sum(layer_base == "counts") > 1 || sum(layer_base == "data") > 1) {
      scRNA_CD45 <- JoinLayers(scRNA_CD45, assay = "RNA")
    }
  }
  scRNA_CD45 <- NormalizeData(scRNA_CD45, assay = "RNA", verbose = FALSE)
  scRNA_CD45 <- SCTransform(scRNA_CD45, assay = "RNA", ncells = 5000,
                            variable.features.n = 3000, vst.flavor = "v2",
                            seed.use = 1448145, verbose = FALSE)
  scRNA_CD45 <- RunPCA(scRNA_CD45, assay = "SCT", npcs = 50,
                       seed.use = 42, verbose = FALSE)
  DefaultAssay(scRNA_CD45) <- "SCT"
  scRNA_CD45 <- harmony::RunHarmony(
    scRNA_CD45, group.by.vars = "orig.ident", reduction.use = "pca",
    dims.use = 1:50, reduction.save = "harmony", project.dim = TRUE,
    verbose = FALSE
  )
  scRNA_CD45 <- FindNeighbors(scRNA_CD45, reduction = "harmony", dims = 1:15,
                              k.param = 20, graph.name = c("SCT_nn", "SCT_snn"),
                              verbose = FALSE)
  scRNA_CD45 <- FindClusters(scRNA_CD45, graph.name = "SCT_snn",
                             resolution = 0.8, algorithm = 1,
                             random.seed = 0, verbose = FALSE)
  scRNA_CD45 <- RunUMAP(scRNA_CD45, reduction = "harmony", dims = 1:10,
                        n.neighbors = 30, min.dist = 0.3, metric = "cosine",
                        umap.method = "uwot", seed.use = 42, verbose = FALSE)

  # Coordinate reflection affects presentation only, not cell membership.
  if (flip_umap1) {
    xy <- Embeddings(scRNA_CD45, "umap")
    xy[, 1] <- -xy[, 1]
    scRNA_CD45[["umap"]] <- SeuratObject::CreateDimReducObject(
      embeddings = xy, key = "UMAP_", assay = "SCT")
  }
  saveRDS(scRNA_CD45, preprocessed_rds)
} else {
  if (!file.exists(preprocessed_rds)) {
    stop("Preprocessed RDS not found: ", preprocessed_rds,
         ". Run with 'preprocess' first.")
  }
  scRNA_CD45 <- readRDS(preprocessed_rds)
}
stopifnot(inherits(scRNA_CD45, "Seurat"))
stopifnot(
  "RNA" %in% names(scRNA_CD45@assays),
  "umap" %in% names(scRNA_CD45@reductions),
  "seurat_clusters" %in% names(scRNA_CD45@meta.data)
)
message("Annotation input: ", ncol(scRNA_CD45), " cells; ",
        length(unique(scRNA_CD45$seurat_clusters)), " global clusters.")

# 2. Marker-module annotation -------------------------------------------------
immune_markers <- list(
  T_cell = c("CD3D", "CD3E", "CD3G", "TRAC", "TRBC1", "TRBC2"),
  CD4_T = c("CD4", "IL7R", "LTB", "MAL", "CCR7", "TCF7", "LEF1"),
  CD8_T = c("CD8A", "CD8B", "CCL5", "GZMK", "NKG7", "CTSW"),
  Treg = c("FOXP3", "IL2RA", "CTLA4", "IKZF2", "TIGIT", "TNFRSF18"),
  Exhausted_T = c("PDCD1", "LAG3", "HAVCR2", "TIGIT", "TOX", "ENTPD1", "CXCL13"),
  Cytotoxic_T_NK = c("NKG7", "GNLY", "PRF1", "GZMB", "CCL5", "CTSW", "FGFBP2"),
  NK = c("NKG7", "GNLY", "KLRD1", "KLRF1", "NCAM1", "XCL1", "XCL2"),
  Cycling = c("MKI67", "TOP2A", "UBE2C", "CDK1", "STMN1", "TUBB"),
  B_cell = c("CD79A", "CD79B", "MS4A1", "CD19", "CD37", "CD74", "BANK1"),
  Plasma_cell = c("JCHAIN", "MZB1", "SDC1", "XBP1", "PRDM1", "DERL3"),
  Microglia = c("P2RY12", "TMEM119", "P2RY13", "GPR34", "SLC2A5", "SALL1", "HEXB", "BIN1"),
  Classical_Monocyte = c("FCN1", "VCAN", "CCR2", "S100A8", "S100A9", "CD14", "SELL"),
  Macrophage_C1QC = c("C1QA", "C1QB", "C1QC", "APOE", "MERTK", "CD163", "MSR1"),
  Macrophage_SPP1 = c("SPP1", "GPNMB", "LGALS3", "APOC1", "FN1", "MMP12", "TREM2"),
  MHCII_Myeloid = c("HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "CD74", "CIITA"),
  cDC1 = c("XCR1", "CLEC9A", "BATF3", "CADM1", "IRF8"),
  cDC2 = c("CD1C", "CLEC10A", "FCER1A", "CD1E"),
  LAMP3_DC = c("LAMP3", "CCR7", "FSCN1", "CCL19", "CCL22", "CD274", "IDO1"),
  pDC = c("CLEC4C", "GZMB", "JCHAIN", "IL3RA", "TCF4", "IRF7"),
  Neutrophil = c("FCGR3B", "CSF3R", "CXCR2", "S100A8", "S100A9", "MMP8", "MMP9", "FPR1"),
  Mast_cell = c("KIT", "CPA3", "TPSAB1", "TPSB2", "MS4A2", "HDC")
)

# Calculate marker-module scores and assign cluster-level candidate labels.
DefaultAssay(scRNA_CD45) <- "RNA"
if (inherits(scRNA_CD45[["RNA"]], "Assay5")) {
  ly <- Layers(scRNA_CD45[["RNA"]])
  if (sum(grepl("^counts(\\.|$)", ly)) > 1 ||
      sum(grepl("^data(\\.|$)", ly)) > 1) {
    scRNA_CD45 <- JoinLayers(scRNA_CD45, assay = "RNA")
  }
  if (!"data" %in% Layers(scRNA_CD45[["RNA"]])) {
    scRNA_CD45 <- NormalizeData(scRNA_CD45, assay = "RNA", verbose = FALSE)
  }
}
score_markers <- lapply(immune_markers, function(g) intersect(g, rownames(scRNA_CD45)))
score_markers <- score_markers[lengths(score_markers) >= 2]
if (length(score_markers) < 2) stop("Too few matched marker panels: check gene symbols.")

existing_score_columns <- grep("^Score_", colnames(scRNA_CD45@meta.data), value = TRUE)
if (length(existing_score_columns)) {
  scRNA_CD45@meta.data[, existing_score_columns] <- NULL
}
for (celltype in names(score_markers)) {
  scRNA_CD45 <- AddModuleScore(
    scRNA_CD45, features = list(score_markers[[celltype]]),
    assay = "RNA", name = paste0("Score_", celltype, "_"), seed = 1234
  )
}
score_cols <- paste0("Score_", names(score_markers), "_1")
stopifnot(all(score_cols %in% colnames(scRNA_CD45@meta.data)))
cluster_id <- as.character(scRNA_CD45$seurat_clusters)
means <- aggregate(scRNA_CD45@meta.data[, score_cols, drop = FALSE],
                   by = list(seurat_clusters = cluster_id), FUN = mean)
means <- means[order(as.numeric(means$seurat_clusters)), , drop = FALSE]
write.table(means, file.path(annotation_dir, "cluster_module_scores.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

winner <- max.col(as.matrix(means[, score_cols, drop = FALSE]), ties.method = "first")
annotation_11 <- data.frame(
  seurat_clusters = as.character(means$seurat_clusters),
  module_top1 = names(score_markers)[winner],
  stringsAsFactors = FALSE
)
write.table(annotation_11,
            file.path(annotation_dir, "cluster_annotation_candidates.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cluster_map <- setNames(annotation_11$module_top1, annotation_11$seurat_clusters)
scRNA_CD45$module_celltype <- unname(cluster_map[cluster_id])
scRNA_CD45$celltype_11 <- as.character(scRNA_CD45$module_celltype)

# Optional cluster-label corrections after marker review.
target_11 <- c("MHCII_Myeloid", "Macrophage_SPP1", "Cytotoxic_T_NK", "T_cell",
               "Cycling", "Treg", "Classical_Monocyte", "Macrophage_C1QC",
               "B_cell", "CD4_T", "Plasma_cell")
if (length(manual_annotation_corrections)) {
  if (is.null(names(manual_annotation_corrections)) ||
      !all(names(manual_annotation_corrections) %in% unique(cluster_id)) ||
      !all(manual_annotation_corrections %in% target_11)) {
    stop("manual_annotation_corrections contains invalid cluster IDs or labels.")
  }
  affected <- cluster_id %in% names(manual_annotation_corrections)
  scRNA_CD45$celltype_11[affected] <- unname(
    manual_annotation_corrections[cluster_id[affected]]
  )
}
saveRDS(scRNA_CD45, file.path(object_dir, "GSE233304_annotation_candidates.rds"))

# 3. Merge related myeloid labels to obtain the final nine cell types --------
map_11_to_9 <- c(
  MHCII_Myeloid = "TAMs", Macrophage_SPP1 = "TAMs", Macrophage_C1QC = "TAMs",
  Classical_Monocyte = "Mono/Mac", Cytotoxic_T_NK = "Cytotoxic_T_NK",
  T_cell = "T_cell", CD4_T = "CD4_T", Treg = "Treg", Cycling = "Cycling",
  B_cell = "B_cell", Plasma_cell = "Plasma_cell"
)
unknown <- setdiff(unique(scRNA_CD45$celltype_11), names(map_11_to_9))
if (length(unknown) || anyNA(scRNA_CD45$celltype_11)) {
  stop("Candidate labels outside the predefined annotation set: ",
       paste(unknown, collapse = ", "),
       ". Review cluster_annotation_candidates.tsv and, if justified, set ",
       "manual_annotation_corrections.")
}
scRNA_CD45$celltype_9_mapped <- unname(map_11_to_9[scRNA_CD45$celltype_11])
scRNA_CD45$celltype_final <- scRNA_CD45$celltype_9_mapped

# Preserve the observed cell-type order in plots.
if (is.factor(scRNA_CD45$celltype_9_mapped)) {
  celltype_order <- levels(droplevels(scRNA_CD45$celltype_9_mapped))
} else {
  celltype_order <- unique(as.character(scRNA_CD45$celltype_9_mapped))
}
celltype_order <- celltype_order[!is.na(celltype_order)]
scRNA_CD45$celltype_9_mapped <- factor(scRNA_CD45$celltype_9_mapped,
                                         levels = celltype_order)
scRNA_CD45$celltype_final <- scRNA_CD45$celltype_9_mapped

# PBMC is shown on the left and tumor on the right.
sample_order <- c("PBMCs_CD45", "Tumor_CD45")
if (!setequal(unique(as.character(scRNA_CD45$sample)), sample_order)) {
  stop("Unexpected sample labels: ", paste(unique(scRNA_CD45$sample), collapse = ", "))
}
scRNA_CD45$sample <- factor(as.character(scRNA_CD45$sample), levels = sample_order)
Idents(scRNA_CD45) <- "celltype_final"

if (length(celltype_order) != 9L) {
  warning("Only ", length(celltype_order), " of the 9 target classes are observed. ",
          "Labels are not fabricated for absent classes.")
}
write.table(data.frame(
  cell = colnames(scRNA_CD45), patient = as.character(scRNA_CD45$orig.ident),
  sample = as.character(scRNA_CD45$sample),
  seurat_clusters = as.character(scRNA_CD45$seurat_clusters),
  celltype_11 = as.character(scRNA_CD45$celltype_11),
  celltype_final = as.character(scRNA_CD45$celltype_final)
), file.path(annotation_dir, "final_cell_annotation.tsv"),
sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(scRNA_CD45, final_rds)

# 4. Split UMAP ---------------------------------------------------------------
if (length(celltype_order) > 20L) stop("category20 supports at most 20 cell types")
celltype_colors <- setNames(ggsci::pal_d3(palette = "category20")(
  length(celltype_order)), celltype_order)

plot_umap <- DimPlot(
  object = scRNA_CD45, reduction = "umap", group.by = "celltype_final",
  split.by = "sample", ncol = 2, pt.size = point_size, label = FALSE,
  raster = FALSE,
  cols = celltype_colors, combine = TRUE
)
plot_umap <- plot_umap + patchwork::plot_layout(guides = "collect")
plot_umap <- plot_umap &
  theme_classic(base_size = 11) &
  guides(color = guide_legend(title = NULL,
                              override.aes = list(size = 3, alpha = 1))) &
  theme(
    strip.text = element_text(size = 11, face = "bold", colour = "black"),
    strip.background = element_blank(),
    axis.title = element_text(size = 11, colour = "black"),
    axis.text = element_text(size = 9, colour = "black"),
    axis.line = element_line(colour = "black", linewidth = 0.5),
    axis.ticks = element_line(colour = "black"),
    legend.position = "right",
    legend.text = element_text(size = 9, colour = "black"),
    legend.key.height = grid::unit(0.4, "cm"),
    legend.key.width = grid::unit(0.4, "cm"),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA)
  )
print(plot_umap)
ggsave(file.path(plot_dir, "GSE233304_split_UMAP.pdf"),
       plot = plot_umap, width = 16, height = 8, units = "in", bg = "white")
ggsave(file.path(plot_dir, "GSE233304_split_UMAP.tiff"),
       plot = plot_umap, width = 16, height = 8, units = "in", dpi = 600,
       compression = "lzw", bg = "white")

# 5. Patient-level paired PPARG analysis -------------------------------------
DefaultAssay(scRNA_CD45) <- "RNA"
if (inherits(scRNA_CD45[["RNA"]], "Assay5")) {
  rna_layers <- Layers(scRNA_CD45[["RNA"]])
  if (sum(grepl("^data(\\.|$)", rna_layers)) > 1L) {
    scRNA_CD45 <- JoinLayers(scRNA_CD45, assay = "RNA")
  }
  if (!"data" %in% Layers(scRNA_CD45[["RNA"]])) {
    scRNA_CD45 <- NormalizeData(scRNA_CD45, assay = "RNA", verbose = FALSE)
  }
}

pparg_cell_data <- tryCatch(
  FetchData(
    scRNA_CD45,
    vars = c("PPARG", "sample", "celltype_final", "orig.ident"),
    layer = "data"
  ),
  error = function(e) {
    FetchData(
      scRNA_CD45,
      vars = c("PPARG", "sample", "celltype_final", "orig.ident"),
      slot = "data"
    )
  }
)
colnames(pparg_cell_data)[colnames(pparg_cell_data) == "orig.ident"] <- "patient"

pparg_cell_data$comparison_group <- dplyr::case_when(
  pparg_cell_data$sample == "PBMCs_CD45" &
    pparg_cell_data$celltype_final == "Mono/Mac" ~ "PBMC Mono/Mac",
  pparg_cell_data$sample == "Tumor_CD45" &
    pparg_cell_data$celltype_final == "Mono/Mac" ~ "Tumor Mono/Mac",
  pparg_cell_data$sample == "Tumor_CD45" &
    pparg_cell_data$celltype_final == "TAMs" ~ "Tumor TAMs",
  TRUE ~ NA_character_
)
pparg_cell_data <- pparg_cell_data |>
  dplyr::filter(!is.na(comparison_group))

required_groups <- c("PBMC Mono/Mac", "Tumor Mono/Mac", "Tumor TAMs")
missing_groups <- setdiff(required_groups, unique(pparg_cell_data$comparison_group))
if (length(missing_groups)) {
  stop("Missing PPARG comparison groups: ", paste(missing_groups, collapse = ", "))
}
pparg_cell_data$comparison_group <- factor(
  pparg_cell_data$comparison_group,
  levels = required_groups
)

# This cell-level plot is descriptive; inference is performed using patients.
plot_pparg_cells <- ggplot(
  pparg_cell_data,
  aes(x = comparison_group, y = PPARG, fill = comparison_group)
) +
  geom_violin(scale = "width", trim = TRUE, linewidth = 0.3) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", linewidth = 0.3) +
  scale_fill_manual(values = c(
    "PBMC Mono/Mac" = "#404040",
    "Tumor Mono/Mac" = "#B65A50",
    "Tumor TAMs" = "#7F1D1D"
  )) +
  labs(x = NULL, y = "PPARG expression (log-normalized)") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1, colour = "black"),
    axis.text.y = element_text(colour = "black"),
    legend.position = "none"
  )

ggsave(file.path(plot_dir, "GSE233304_PPARG_cell_distribution.pdf"),
       plot = plot_pparg_cells, width = 6, height = 5, units = "in", bg = "white")
ggsave(file.path(plot_dir, "GSE233304_PPARG_cell_distribution.tiff"),
       plot = plot_pparg_cells, width = 6, height = 5, units = "in", dpi = 600,
       compression = "lzw", bg = "white")

# Summarize each patient before testing to avoid treating cells as replicates.
pparg_patient_summary <- pparg_cell_data |>
  dplyr::group_by(patient, comparison_group) |>
  dplyr::summarise(
    PPARG_mean_log = mean(PPARG, na.rm = TRUE),
    PPARG_mean = log1p(mean(expm1(PPARG), na.rm = TRUE)),
    PPARG_positive_percent = mean(PPARG > 0, na.rm = TRUE) * 100,
    n_cells = dplyr::n(),
    .groups = "drop"
  )

write.table(
  pparg_patient_summary,
  file.path(statistics_dir, "PPARG_patient_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

p_value_label <- function(p) {
  if (is.na(p)) return("P = NA")
  if (p < 0.001) return("P < 0.001")
  paste0("P = ", formatC(p, digits = 3, format = "f"))
}

build_paired_comparison <- function(data, group_1, group_2, comparison_name) {
  paired_wide <- data |>
    dplyr::filter(comparison_group %in% c(group_1, group_2)) |>
    dplyr::select(patient, comparison_group, PPARG_mean) |>
    tidyr::pivot_wider(names_from = comparison_group, values_from = PPARG_mean) |>
    dplyr::filter(!is.na(.data[[group_1]]), !is.na(.data[[group_2]]))

  if (nrow(paired_wide) < 2L) {
    stop("Fewer than two complete patient pairs for: ", comparison_name)
  }

  test_result <- stats::t.test(
    paired_wide[[group_2]], paired_wide[[group_1]],
    paired = TRUE, alternative = "two.sided"
  )

  paired_long <- paired_wide |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(c(group_1, group_2)),
      names_to = "comparison_group",
      values_to = "PPARG_mean"
    ) |>
    dplyr::mutate(
      comparison = comparison_name,
      comparison_group = factor(comparison_group, levels = c(group_1, group_2))
    )

  statistics <- data.frame(
    comparison = comparison_name,
    group_1 = group_1,
    group_2 = group_2,
    n_pairs = nrow(paired_wide),
    mean_difference_group_2_minus_group_1 = mean(
      paired_wide[[group_2]] - paired_wide[[group_1]]
    ),
    t_statistic = unname(test_result$statistic),
    degrees_of_freedom = unname(test_result$parameter),
    p_value = test_result$p.value,
    stringsAsFactors = FALSE
  )

  annotation <- data.frame(
    comparison = comparison_name,
    group1 = group_1,
    group2 = group_2,
    y.position = max(paired_long$PPARG_mean, na.rm = TRUE) * 1.12 + 0.02,
    p_label = p_value_label(test_result$p.value),
    stringsAsFactors = FALSE
  )

  list(long = paired_long, statistics = statistics, annotation = annotation)
}

comparison_monomac <- build_paired_comparison(
  pparg_patient_summary,
  "PBMC Mono/Mac", "Tumor Mono/Mac",
  "PBMC Mono/Mac vs Tumor Mono/Mac"
)
comparison_tams <- build_paired_comparison(
  pparg_patient_summary,
  "PBMC Mono/Mac", "Tumor TAMs",
  "PBMC Mono/Mac vs Tumor TAMs"
)

paired_plot_data <- dplyr::bind_rows(
  comparison_monomac$long,
  comparison_tams$long
)
paired_statistics <- dplyr::bind_rows(
  comparison_monomac$statistics,
  comparison_tams$statistics
)
paired_annotations <- dplyr::bind_rows(
  comparison_monomac$annotation,
  comparison_tams$annotation
)

write.table(
  paired_plot_data,
  file.path(statistics_dir, "PPARG_paired_plot_data.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  paired_statistics,
  file.path(statistics_dir, "PPARG_paired_test_results.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

plot_pparg_paired <- ggplot(
  paired_plot_data,
  aes(x = comparison_group, y = PPARG_mean, group = patient)
) +
  geom_line(colour = "grey65", linewidth = 0.6) +
  geom_point(aes(colour = comparison_group), size = 3) +
  ggpubr::stat_pvalue_manual(
    paired_annotations,
    label = "p_label",
    xmin = "group1", xmax = "group2", y.position = "y.position",
    tip.length = 0.01, bracket.size = 0.4
  ) +
  facet_wrap(~ comparison, scales = "free_x") +
  scale_colour_manual(values = c(
    "PBMC Mono/Mac" = "#404040",
    "Tumor Mono/Mac" = "#B65A50",
    "Tumor TAMs" = "#7F1D1D"
  )) +
  labs(x = NULL, y = "Patient-level mean PPARG expression") +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 25, hjust = 1, colour = "black"),
    axis.text.y = element_text(colour = "black"),
    legend.position = "none"
  )

ggsave(file.path(plot_dir, "GSE233304_PPARG_paired_patient.pdf"),
       plot = plot_pparg_paired, width = 9, height = 5, units = "in", bg = "white")
ggsave(file.path(plot_dir, "GSE233304_PPARG_paired_patient.tiff"),
       plot = plot_pparg_paired, width = 9, height = 5, units = "in", dpi = 600,
       compression = "lzw", bg = "white")

capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
message("Analysis complete. Final object: ", final_rds)
