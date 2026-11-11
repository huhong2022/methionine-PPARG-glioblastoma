#!/usr/bin/env Rscript

# ssGSEA, ESTIMATE, and PPARG-high versus PPARG-low GSEA.
# Usage:
#   Rscript 02_ssGSEA_GSEA.R all
#   Rscript 02_ssGSEA_GSEA.R ssgsea
#   Rscript 02_ssGSEA_GSEA.R estimate
#   Rscript 02_ssGSEA_GSEA.R gsea
# An optional second argument specifies the repository root.

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
task <- if (length(args) >= 1) tolower(args[[1]]) else "all"
valid_tasks <- c("all", "ssgsea", "estimate", "gsea")
if (!task %in% valid_tasks) {
  stop("First argument must be one of: ", paste(valid_tasks, collapse = ", "))
}
project_root <- if (length(args) >= 2) args[[2]] else getwd()
project_root <- normalizePath(project_root, mustWork = TRUE)

data_dir <- file.path(project_root, "data")
gene_set_dir <- file.path(project_root, "gene_lists")
output_dir <- file.path(project_root, "results", "02_ssGSEA_GSEA")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

expression_file <- file.path(data_dir, "TCGA_GBM_IDHWT_TPM.tsv")
methionine_gmt <- file.path(
  gene_set_dir, "methionine_metabolism_activity_gene_set.gmt"
)
go_gmt <- file.path(gene_set_dir, "c5.go.v7.4.symbols.gmt")
target_gene <- "PPARG"
target_terms <- c(
  "GOBP_NEGATIVE_REGULATION_OF_INFLAMMATORY_RESPONSE",
  "GOBP_NEGATIVE_REGULATION_OF_T_CELL_PROLIFERATION"
)
min_gs_size <- 10
max_gs_size <- 500
random_seed <- 123L

check_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop("Missing R packages: ", paste(missing, collapse = ", "))
  }
}

read_expression <- function(path) {
  if (!file.exists(path)) stop("Expression file not found: ", path)
  x <- read.delim(
    path, header = TRUE, check.names = FALSE, quote = "",
    comment.char = "", stringsAsFactors = FALSE
  )
  if (ncol(x) < 3) stop("Expression table must contain a gene column and >=2 samples.")
  gene_id <- as.character(x[[1]])
  matrix_data <- as.matrix(x[, -1, drop = FALSE])
  storage.mode(matrix_data) <- "numeric"
  if (anyNA(matrix_data)) stop("Expression matrix contains NA or non-numeric values.")
  valid <- !is.na(gene_id) & nzchar(trimws(gene_id))
  matrix_data <- matrix_data[valid, , drop = FALSE]
  rownames(matrix_data) <- gene_id[valid]
  matrix_data <- limma::avereps(matrix_data)
  if (any(matrix_data < 0)) stop("TPM matrix contains negative values.")
  matrix_data
}

read_gmt_simple <- function(path) {
  if (!file.exists(path)) stop("GMT file not found: ", path)
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  sets <- lapply(lines, function(line) {
    fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
    if (length(fields) < 3) stop("Invalid GMT line in ", path)
    unique(fields[3:length(fields)])
  })
  names(sets) <- vapply(lines, function(line) {
    strsplit(line, "\t", fixed = TRUE)[[1]][1]
  }, character(1))
  sets
}

gmt_to_term2gene <- function(path) {
  sets <- read_gmt_simple(path)
  do.call(rbind, lapply(names(sets), function(term) {
    data.frame(term = term, gene = sets[[term]], stringsAsFactors = FALSE)
  }))
}

if (!file.exists(expression_file)) {
  stop(
    "Required file is missing: ", expression_file,
    "\nRename the untransformed gene-by-sample TPM matrix to ",
    "TCGA_GBM_IDHWT_TPM.tsv."
  )
}

check_packages("limma")
expression_tpm <- read_expression(expression_file)
expression_log2 <- log2(expression_tpm + 1)

patient_id <- substr(colnames(expression_log2), 1, 12)
duplicate_patients <- names(which(table(patient_id) > 1))
if (length(duplicate_patients) > 0) {
  duplicate_audit <- data.frame(
    Sample = colnames(expression_log2),
    Patient = patient_id,
    Is_duplicated_patient = patient_id %in% duplicate_patients,
    stringsAsFactors = FALSE
  )
  write.table(
    duplicate_audit,
    file.path(output_dir, "TCGA_duplicate_patient_audit.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  message(
    length(duplicate_patients),
    " patients have multiple tumor samples. All samples are retained for the ",
    "sample-level expression analysis; see TCGA_duplicate_patient_audit.tsv."
  )
}

run_ssgsea <- function() {
  check_packages("GSVA")
  gene_sets <- read_gmt_simple(methionine_gmt)
  if (!"Methionine metabolism" %in% names(gene_sets)) {
    stop("Expected 'Methionine metabolism' set not found in ", methionine_gmt)
  }
  if (length(gene_sets[["Methionine metabolism"]]) != 38) {
    stop("Methionine metabolism gene set must contain exactly 38 unique genes.")
  }
  overlap <- intersect(rownames(expression_log2), gene_sets[["Methionine metabolism"]])
  if (length(overlap) < 10) {
    stop("Only ", length(overlap), " methionine genes overlap the expression matrix.")
  }
  writeLines(
    c(
      paste0("Gene set size: ", length(gene_sets[["Methionine metabolism"]])),
      paste0("Genes found in expression matrix: ", length(overlap)),
      paste0("Missing genes: ", paste(setdiff(gene_sets[["Methionine metabolism"]], overlap), collapse = ", "))
    ),
    file.path(output_dir, "ssGSEA_gene_overlap.txt")
  )

  if ("ssgseaParam" %in% getNamespaceExports("GSVA")) {
    parameter <- GSVA::ssgseaParam(
      exprData = expression_log2,
      geneSets = gene_sets,
      normalize = TRUE
    )
    score <- GSVA::gsva(parameter, verbose = TRUE)
  } else {
    score <- GSVA::gsva(
      expression_log2, gene_sets, method = "ssgsea",
      ssgsea.norm = TRUE, verbose = TRUE
    )
  }
  score_table <- data.frame(
    id = colnames(score),
    t(score),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  colnames(score_table)[2] <- "Methionine metabolism Score"
  write.table(
    score_table,
    file.path(output_dir, "TCGA_GBM_IDHWT_methionine_ssGSEA_scores.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
}

run_estimate <- function() {
  check_packages("estimate")

  # Load the ESTIMATE reference datasets explicitly for consistent execution
  # across package installations.
  estimate_data_env <- new.env(parent = asNamespace("estimate"))
  for (dataset_name in c("common_genes", "SI_geneset")) {
    suppressWarnings(utils::data(
      list = dataset_name,
      package = "estimate",
      envir = estimate_data_env
    ))
    if (!exists(dataset_name, envir = estimate_data_env, inherits = FALSE)) {
      stop(
        "The estimate package is missing its built-in dataset: ", dataset_name,
        ". Reinstall the official estimate package with its data files."
      )
    }
  }
  common_genes <- get("common_genes", envir = estimate_data_env)
  if (!is.data.frame(common_genes) ||
      !all(c("GeneSymbol", "EntrezID") %in% names(common_genes))) {
    stop("estimate::common_genes has an unexpected structure; please reinstall estimate.")
  }

  filterCommonGenes_fixed <- estimate::filterCommonGenes
  environment(filterCommonGenes_fixed) <- estimate_data_env
  estimateScore_fixed <- estimate::estimateScore
  environment(estimateScore_fixed) <- estimate_data_env

  estimate_input <- data.frame(
    GeneSymbol = rownames(expression_log2),
    expression_log2,
    check.names = FALSE
  )
  estimate_input_path <- file.path(output_dir, "ESTIMATE_log2TPM_input.tsv")
  common_genes_path <- file.path(output_dir, "ESTIMATE_common_genes.gct")
  raw_score_path <- file.path(output_dir, "ESTIMATE_raw_scores.gct")
  write.table(
    estimate_input, estimate_input_path, sep = "\t", quote = FALSE,
    row.names = FALSE
  )

  filterCommonGenes_fixed(
    input.f = estimate_input_path,
    output.f = common_genes_path,
    id = "GeneSymbol"
  )
  if (!file.exists(common_genes_path) || file.info(common_genes_path)$size == 0) {
    stop("ESTIMATE common-gene GCT was not created or is empty.")
  }

  estimateScore_fixed(
    input.ds = common_genes_path,
    output.ds = raw_score_path,
    platform = "illumina"
  )
  if (!file.exists(raw_score_path) || file.info(raw_score_path)$size == 0) {
    stop("ESTIMATE score GCT was not created or is empty.")
  }

  scores <- read.delim(
    raw_score_path, skip = 2, header = TRUE, check.names = FALSE,
    stringsAsFactors = FALSE
  )
  if (ncol(scores) < 4) {
    stop("Unexpected ESTIMATE score GCT format: fewer than 4 columns.")
  }
  rownames(scores) <- scores[[1]]
  score_matrix <- t(as.matrix(scores[, 3:ncol(scores), drop = FALSE]))
  storage.mode(score_matrix) <- "numeric"
  if (anyNA(score_matrix) || !all(c("StromalScore", "ImmuneScore", "ESTIMATEScore")
                                    %in% colnames(score_matrix))) {
    stop("ESTIMATE score output is missing expected score columns or has NA values.")
  }
  score_table <- data.frame(
    id = rownames(score_matrix), score_matrix, check.names = FALSE
  )
  write.table(
    score_table,
    file.path(output_dir, "TCGA_GBM_IDHWT_ESTIMATE_scores.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  message("ESTIMATE scores saved for ", nrow(score_table), " samples.")
}

run_gsea <- function() {
  check_packages(c("clusterProfiler", "BiocParallel", "GseaVis", "cowplot", "ggplot2"))
  if (!target_gene %in% rownames(expression_log2)) {
    stop("Target gene not found in expression matrix: ", target_gene)
  }
  if (!file.exists(go_gmt)) {
    stop(
      "MSigDB GO GMT file not found: ", go_gmt,
      "\nDownload c5.go.v7.4.symbols.gmt from MSigDB and place it in gene_lists/."
    )
  }
  BiocParallel::register(BiocParallel::SerialParam(), default = TRUE)

  target_expression <- expression_log2[target_gene, ]
  cutoff <- stats::median(target_expression, na.rm = TRUE)
  group <- factor(
    ifelse(target_expression > cutoff, "High", "Low"),
    levels = c("Low", "High")
  )
  group_table <- data.frame(
    Sample = names(target_expression),
    PPARG_log2_TPM_plus_1 = as.numeric(target_expression),
    Median_cutoff = cutoff,
    Group = as.character(group),
    stringsAsFactors = FALSE
  )
  write.table(
    group_table,
    file.path(output_dir, "PPARG_median_groups.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  if (length(unique(group)) != 2) stop("Median split did not create two groups.")

  design <- stats::model.matrix(~0 + group)
  colnames(design) <- levels(group)
  fit <- limma::lmFit(expression_log2, design)
  contrast <- limma::makeContrasts(High_vs_Low = High - Low, levels = design)
  fit <- limma::eBayes(limma::contrasts.fit(fit, contrast))
  deg <- limma::topTable(fit, coef = "High_vs_Low", number = Inf, sort.by = "P")
  deg$Gene <- rownames(deg)
  deg <- deg[, c("Gene", setdiff(colnames(deg), "Gene"))]
  write.table(
    deg,
    file.path(output_dir, "PPARG_High_vs_Low_limma.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )

  rank_metric <- fit$t[, "High_vs_Low"]
  names(rank_metric) <- rownames(fit$t)
  rank_metric <- rank_metric[
    is.finite(rank_metric) & !is.na(names(rank_metric)) & nzchar(names(rank_metric))
  ]
  rank_metric <- rank_metric[names(rank_metric) != target_gene]
  rank_metric <- sort(rank_metric, decreasing = TRUE)
  rank_table <- data.frame(
    Gene = names(rank_metric),
    Moderated_t = as.numeric(rank_metric),
    stringsAsFactors = FALSE
  )
  write.table(
    rank_table,
    file.path(output_dir, "PPARG_GSEA_rank_metric.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )

  term2gene <- gmt_to_term2gene(go_gmt)
  missing_terms <- setdiff(target_terms, unique(term2gene$term))
  if (length(missing_terms) > 0) {
    stop("Target pathways absent from GMT: ", paste(missing_terms, collapse = ", "))
  }
  overlap_table <- do.call(rbind, lapply(target_terms, function(term) {
    genes <- unique(term2gene$gene[term2gene$term == term])
    data.frame(
      Term = term,
      GMT_genes = length(genes),
      Ranked_genes_overlap = length(intersect(genes, names(rank_metric))),
      stringsAsFactors = FALSE
    )
  }))
  write.table(
    overlap_table,
    file.path(output_dir, "GSEA_target_term_overlap.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )

  set.seed(random_seed)
  gsea_args <- list(
    geneList = rank_metric,
    TERM2GENE = term2gene,
    minGSSize = min_gs_size,
    maxGSSize = max_gs_size,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    verbose = TRUE,
    seed = TRUE
  )
  if ("by" %in% names(formals(clusterProfiler::GSEA))) {
    gsea_args$by <- "fgsea"
  }
  gsea_result <- do.call(clusterProfiler::GSEA, gsea_args)
  result_table <- as.data.frame(gsea_result)
  if (nrow(result_table) == 0) stop("GSEA returned no pathways.")
  write.table(
    result_table,
    file.path(output_dir, "PPARG_GSEA_all_results.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  write.table(
    result_table[!is.na(result_table$p.adjust) & result_table$p.adjust < 0.05, , drop = FALSE],
    file.path(output_dir, "PPARG_GSEA_FDR_less_0.05.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )

  # Preserve the complete GSEA result object for plotting and audit.
  saveRDS(
    gsea_result,
    file = file.path(output_dir, "PPARG_GSEA_result.rds")
  )

  # Publication-style GSEA plotting settings.
  high_color <- "#A4473D"     # PPARG-high: muted dark red
  low_color <- "#3E6C8E"      # PPARG-low: muted blue
  neutral_color <- "#4A4A4A"

  format_pvalue <- function(x) {
    if (length(x) != 1L || is.na(x) || !is.finite(x)) return("NA")
    if (x < 0.001) return("<0.001")
    sprintf("%.3f", x)
  }

  # Ranking statistic is the moderated t of High - Low:
  # positive rank -> PPARG-high, negative rank -> PPARG-low.
  positive_group <- "PPARG-high"
  negative_group <- "PPARG-low"

  plot_audit <- list()

  # Draw one plot for each target GO term.
  for (term in target_terms) {
    row_index <- which(
      (!is.na(result_table$ID) & result_table$ID == term) |
      (!is.na(result_table$Description) & result_table$Description == term)
    )
    if (length(row_index) == 0L) {
      warning("Target term not returned by GSEA: ", term)
      plot_audit[[term]] <- data.frame(
        Term = term, Status = "not_returned", NES = NA_real_,
        P_value = NA_real_, FDR = NA_real_, stringsAsFactors = FALSE
      )
      next
    }

    term_result <- result_table[row_index[1L], , drop = FALSE]
    plot_term_id <- as.character(term_result$ID[[1L]])
    nes_value <- as.numeric(term_result$NES[[1L]])
    p_value <- as.numeric(term_result$pvalue[[1L]])
    fdr_value <- as.numeric(term_result$p.adjust[[1L]])

    if (!is.finite(nes_value)) {
      warning("Target term has no finite NES, skipping plot: ", term)
      plot_audit[[term]] <- data.frame(
        Term = term, Status = "missing_NES", NES = nes_value,
        P_value = p_value, FDR = fdr_value, stringsAsFactors = FALSE
      )
      next
    }

    enriched_group <- if (nes_value > 0) {
      positive_group
    } else if (nes_value < 0) {
      negative_group
    } else {
      "Neither group (NES = 0)"
    }

    curve_color <- if (nes_value > 0) {
      high_color
    } else if (nes_value < 0) {
      low_color
    } else {
      neutral_color
    }

    # Use the returned pathway description as the plot title.
    term_title <- as.character(term_result$Description[[1L]])
    if (is.na(term_title) || !nzchar(term_title)) term_title <- term
    term_title <- sub("^GOBP_", "", term_title)
    term_title <- gsub("_", " ", term_title)
    term_title <- tools::toTitleCase(tolower(term_title))

    direction_label <- paste0("Enriched in ", enriched_group)
    statistics_label <- paste0(
      "NES = ", sprintf("%.2f", nes_value),
      "   |   nominal P ", format_pvalue(p_value),
      "   |   FDR ", format_pvalue(fdr_value)
    )

    # Enrichment curve and ranked-gene barcode.
    p_gsea <- GseaVis::gseaNb(
      object = gsea_result,
      geneSetID = plot_term_id,
      subPlot = 2,
      lineSize = 1.1,
      termWidth = 60,
      htHeight = 0.25,
      addPval = FALSE,
      curveCol = curve_color,
      htCol = c(low_color, high_color),
      rankCol = "grey70",
      segCol = "grey35",
      rmPrefix = TRUE,
      legend.position = "none",
      base_size = 12
    )

    # Pathway title, enrichment direction, and statistics.
    p_header <- cowplot::ggdraw() +
      cowplot::draw_label(
        label = term_title,
        x = 0.5, y = 0.78, hjust = 0.5, vjust = 0.5,
        fontface = "bold", size = 13, color = "#222222"
      ) +
      cowplot::draw_label(
        label = direction_label,
        x = 0.5, y = 0.39, hjust = 0.5, vjust = 0.5,
        fontface = "bold", size = 11, color = curve_color
      ) +
      cowplot::draw_label(
        label = statistics_label,
        x = 0.5, y = 0.08, hjust = 0.5, vjust = 0.5,
        size = 9.5, color = "#333333"
      )

    # PPARG group labels and ranking direction.
    p_footer <- cowplot::ggdraw() +
      cowplot::draw_label(
        label = paste0(positive_group, "\nPositive rank"),
        x = 0.02, y = 0.52, hjust = 0, vjust = 0.5,
        fontface = "bold", size = 10, color = high_color
      ) +
      cowplot::draw_label(
        label = "Ranked gene list",
        x = 0.5, y = 0.52, hjust = 0.5, vjust = 0.5,
        size = 9, color = "#555555"
      ) +
      cowplot::draw_label(
        label = paste0(negative_group, "\nNegative rank"),
        x = 0.98, y = 0.52, hjust = 1, vjust = 0.5,
        fontface = "bold", size = 10, color = low_color
      )

    # Assemble the final figure.
    gsea_final <- cowplot::plot_grid(
      p_header, p_gsea, p_footer,
      ncol = 1,
      rel_heights = c(0.18, 1, 0.11),
      align = "v"
    )

    print(gsea_final)

    # Save vector PDF and 600-dpi TIFF outputs.
    safe_term <- sub("^GOBP_", "", term)
    output_prefix <- file.path(output_dir, paste0("SFig6_GSEA_", safe_term))

    ggplot2::ggsave(
      filename = paste0(output_prefix, ".pdf"),
      plot = gsea_final,
      width = 8, height = 6.5, units = "in", bg = "white"
    )
    ggplot2::ggsave(
      filename = paste0(output_prefix, ".tiff"),
      plot = gsea_final,
      width = 8, height = 6.5, units = "in",
      dpi = 600, compression = "lzw", bg = "white"
    )

    plot_audit[[term]] <- data.frame(
      Term = term, Status = "saved", NES = nes_value,
      P_value = p_value, FDR = fdr_value, stringsAsFactors = FALSE
    )
    message("GSEA plot saved: ", output_prefix, ".pdf / .tiff")
  }

  if (length(plot_audit) > 0L) {
    write.table(
      do.call(rbind, plot_audit),
      file.path(output_dir, "GSEA_target_plot_status.tsv"),
      sep = "\t", quote = FALSE, row.names = FALSE
    )
  }

}

if (task %in% c("all", "ssgsea")) run_ssgsea()
if (task %in% c("all", "estimate")) run_estimate()
if (task %in% c("all", "gsea")) run_gsea()

writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
message("Requested analysis completed. Results: ", output_dir)

