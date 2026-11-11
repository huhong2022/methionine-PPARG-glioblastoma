#!/usr/bin/env Rscript
# ============================================================================
# Bulk RNA-seq analyses for Supplementary Figures 1A-C, 6C, and 11A.
# Usage:
#   Rscript 01_bulk_RNA_analysis.R [project_root]
# ============================================================================
options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) args[[1]] else getwd()
project_root <- normalizePath(project_root, mustWork = TRUE)

data_dir <- file.path(project_root, "data")
output_dir <- file.path(project_root, "results", "01_bulk_RNA_analysis")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

required_packages <- c("ggplot2", "ComplexHeatmap", "circlize")
missing_packages <- required_packages[!vapply(required_packages,
  requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "),
       ". Install ggplot2/circlize from CRAN and ComplexHeatmap from Bioconductor.")
}

# Required cohort-level inputs.
tcga_score_file <- file.path(data_dir, "TCGA_methionine_scores.tsv")
tcga_clinical_file <- file.path(data_dir, "TCGA_clinical.tsv")
tcga_all_file <- file.path(data_dir, "TCGA_all_glioma_scores.tsv")
tcga_idhwt_file <- file.path(data_dir, "TCGA_GBM_IDHWT_scores.tsv")
tcga_heatmap_file <- file.path(data_dir, "TCGA_GBM_IDHWT_immunosuppression.tsv")
cgga_score_file <- file.path(data_dir, "CGGA_methionine_scores.tsv")
cgga_gbm_file <- file.path(data_dir, "CGGA_GBM_PPARG_CD274.tsv")
gravendeel_file <- file.path(data_dir, "Gravendeel_PPARG_CD274.tsv")

# Optional inputs for additional S Fig. 11A cohort subsets. Their availability
# is recorded in the output status table.
cgga_all_file <- file.path(data_dir, "CGGA_all_glioma_PPARG_CD274.tsv")
gravendeel_clinical_file <- file.path(data_dir, "Gravendeel_clinical.tsv")

score_column <- "Methionine metabolism Score"
violin_levels <- c("Normal", "WHO II", "WHO III", "GBM-IDH-WT", "WHO 4 IDH-mut")

read_tsv <- function(path) {
  if (!file.exists(path)) stop("Input file not found: ", path)
  if (is.na(file.info(path)$size) || file.info(path)$size == 0L) {
    stop("Input file is empty: ", path)
  }
  x <- read.delim(path, header = TRUE, check.names = FALSE, quote = "",
                  comment.char = "", fileEncoding = "UTF-8",
                  stringsAsFactors = FALSE)
  # Duplicate headers are handled by the CGGA parser using column position and
  # sample identifiers.
  colnames(x) <- trimws(sub("^\ufeff", "", colnames(x)))
  x
}
require_columns <- function(x, cols, source) {
  missing <- setdiff(cols, colnames(x))
  if (length(missing)) stop(source, " missing columns: ", paste(missing, collapse = ", "))
}
check_unique_ids <- function(ids, source) {
  ids <- as.character(ids)
  if (anyNA(ids) || any(!nzchar(trimws(ids)))) {
    stop(source, " contains missing IDs")
  }
  if (anyDuplicated(ids)) stop(source, " contains duplicated IDs")
}
as_numeric_checked <- function(x, label) {
  z <- suppressWarnings(as.numeric(as.character(x)))
  bad <- is.na(z) & !is.na(x) & nzchar(trimws(as.character(x)))
  if (any(bad)) stop(label, " contains non-numeric values")
  z
}
write_tsv <- function(x, filename) {
  write.table(x, file.path(output_dir, filename), sep = "\t", quote = FALSE,
              row.names = FALSE, na = "NA")
}
check_expected_n <- function(data, n, label) {
  if (nrow(data) != n) warning(label, ": ", nrow(data), " rows; expected ", n)
}

# Column name lookup ignores spaces, punctuation and case. Avoid ambiguous
# mappings: a duplicate clinical/grade field must be handled explicitly.
canonical <- function(x) tolower(gsub("[^A-Za-z0-9]", "", x))
find_column <- function(x, aliases, mandatory = FALSE, label = "column") {
  idx <- which(canonical(colnames(x)) %in% canonical(aliases))
  if (mandatory && !length(idx)) {
    stop("Cannot find ", label, ". Available columns: ",
         paste(colnames(x), collapse = " | "))
  }
  if (!length(idx)) return(NA_integer_)
  idx[[1]]
}
idh_aliases <- c("IDH_mutation_status", "IDH_status", "IDH mutation status",
                 "IDH", "IDH1_status", "IDH1 mutation status", "IDH1/2_status",
                 "IDH1_2_status", "IDH_Status", "IDH1")
grade_aliases <- c("Grade", "WHO_grade", "WHO grade", "WHO_Grade", "Tumor_Grade")

# Deliberately do not guess that unlabelled 0/1 denotes WT/mutant.
normalize_idh <- function(x) {
  z <- tolower(trimws(as.character(x)))
  z <- gsub("[[:space:]_./-]+", "", z)
  ans <- rep(NA_character_, length(z))
  wt <- z %in% c("wt", "wildtype", "wild", "idhwt", "idhwildtype",
                 "idh1wildtype", "idh2wildtype", "idh1wt", "idh2wt", "nonmutant", "unmutated",
                 "negative", "idhnegative") |
    grepl("wildtype|wildtypeidh|nonmutant|unmutated", z)
  mut <- z %in% c("mut", "mutant", "mutated", "idhmut", "idhmutant",
                  "idhmutated", "idh1mutant", "idh2mutant", "idh1mut", "idh2mut", "positive") |
    grepl("mutant|idh[12]?mutated|idh[12]?mutation", z)
  ans[wt] <- "WT"
  ans[mut & !wt] <- "mut"
  ans
}
normalize_grade <- function(x) {
  raw <- toupper(trimws(as.character(x)))
  raw <- gsub("[_-]", " ", raw)
  raw <- gsub("[[:space:]]+", " ", raw)
  out <- rep(NA_character_, length(raw))
  out[raw %in% c("NORMAL", "CONTROL", "NON TUMOR", "NON-TUMOR")] <- "Normal"
  is_2 <- grepl("^(II|2|WHO II|WHO 2|GRADE II|GRADE 2|WHO GRADE II|WHO GRADE 2)$", raw)
  is_3 <- grepl("^(III|3|WHO III|WHO 3|GRADE III|GRADE 3|WHO GRADE III|WHO GRADE 3)$", raw)
  is_4 <- grepl("^(IV|4|GBM|GLIOBLASTOMA|WHO IV|WHO 4|GRADE IV|GRADE 4|WHO GRADE IV|WHO GRADE 4)$", raw) |
    grepl("GBM|GLIOBLASTOMA|WHO 4|WHO IV|GRADE 4|GRADE IV", raw)
  out[is_2] <- "WHO II"
  out[is_3] <- "WHO III"
  out[is_4] <- "WHO 4"
  out
}
idh_from_grade_text <- function(x) {
  x <- tolower(trimws(as.character(x)))
  out <- rep(NA_character_, length(x))
  out[grepl("idh[ -]?(wt|wild[ -]?type)", x)] <- "WT"
  out[grepl("idh[ -]?(mut|mutant|mutated)", x)] <- "mut"
  out
}
clinical_group <- function(grade_text, idh_text) {
  g <- normalize_grade(grade_text)
  idh <- normalize_idh(idh_text)
  inferred <- idh_from_grade_text(grade_text)
  idh[is.na(idh)] <- inferred[is.na(idh)]
  out <- rep(NA_character_, length(g))
  out[g == "Normal" & !is.na(g)] <- "Normal"
  out[g == "WHO II" & !is.na(g)] <- "WHO II"
  out[g == "WHO III" & !is.na(g)] <- "WHO III"
  out[which(g == "WHO 4" & idh == "WT")] <- "GBM-IDH-WT"
  out[which(g == "WHO 4" & idh == "mut")] <- "WHO 4 IDH-mut"
  out
}

normalize_score_columns <- function(data, source) {
  if (!(score_column %in% colnames(data))) {
    pos <- find_column(data, c("Methionine_score", "methionine_score"),
                       mandatory = TRUE, label = paste0(source, " score"))
    colnames(data)[pos] <- score_column
  }
  id_pos <- find_column(data, c("id", "Sample", "sample_id"),
                        mandatory = TRUE, label = paste0(source, " sample ID"))
  colnames(data)[id_pos] <- "id"
  data
}

# CGGA score and clinical fields may be stored side by side. IDH metadata is
# joined by CGGA_ID rather than row position.
clean_cgga_score <- function(raw) {
  cn <- colnames(raw)
  id_pos <- which(canonical(cn) %in% canonical(c("CGGA_ID", "id")))
  sc_pos <- which(canonical(cn) %in% canonical(c(score_column,
                                                 "methionine_score", "Methionine_score")))
  grade_pos <- which(canonical(cn) %in% canonical(grade_aliases))
  idh_pos <- which(canonical(cn) %in% canonical(idh_aliases))
  if (!length(id_pos) || !length(sc_pos) || !length(grade_pos)) {
    stop("CGGA score input requires CGGA_ID, methionine_score and Grade")
  }
  base <- data.frame(
    id = as.character(raw[[id_pos[1]]]),
    Score = as_numeric_checked(raw[[sc_pos[1]]], "CGGA methionine score"),
    Grade = as.character(raw[[grade_pos[1]]]), stringsAsFactors = FALSE
  )
  check_unique_ids(base$id, "CGGA methionine score IDs")
  if (!length(idh_pos)) {
    stop("CGGA score file has no IDH mutation status column. ",
         "Include IDH_mutation_status (with its CGGA_ID) to split WHO grade 4.")
  }
  if (length(id_pos) >= 2L) {
    # The clinical half begins at the second CGGA_ID. Locate its IDH column.
    clinical_id_pos <- id_pos[2]
    if (idh_pos[1] < clinical_id_pos) {
      stop("CGGA IDH field is before the second CGGA_ID; check source layout")
    }
    clinical <- data.frame(
      id = as.character(raw[[clinical_id_pos]]),
      IDH = as.character(raw[[idh_pos[1]]]), stringsAsFactors = FALSE
    )
    ok <- !is.na(clinical$id) & nzchar(trimws(clinical$id))
    clinical <- clinical[ok, , drop = FALSE]
    check_unique_ids(clinical$id, "CGGA clinical IDs")
    index <- match(base$id, clinical$id)
    base$IDH <- clinical$IDH[index]
    # Check grade consistency across the halves when both columns are present.
    grade_clinical_pos <- grade_pos[grade_pos > clinical_id_pos]
    if (length(grade_clinical_pos)) {
      clinical_grade <- as.character(raw[[grade_clinical_pos[1]]])[ok]
      g1 <- normalize_grade(base$Grade)
      g2 <- normalize_grade(clinical_grade[index])
      mismatches <- !is.na(g1) & !is.na(g2) & g1 != g2
      if (any(mismatches)) warning("CGGA grade disagrees between score and clinical tables for ",
                                   sum(mismatches), " matched IDs")
    }
    if (anyNA(index)) warning("CGGA: ", sum(is.na(index)),
                              " score IDs have no matching clinical ID")
  } else {
    # Only one ID column: assumed that Grade, IDH and score belong to that ID.
    base$IDH <- as.character(raw[[idh_pos[1]]])
  }
  colnames(base)[colnames(base) == "Score"] <- score_column
  base
}

clean_tcga_clinical <- function(raw) {
  id_pos <- find_column(raw, c("Id", "Patient", "Patient_ID", "patient_id", "id"),
                        mandatory = TRUE, label = "TCGA clinical patient ID")
  grade_pos <- find_column(raw, grade_aliases, mandatory = TRUE,
                           label = "TCGA clinical Grade")
  idh_pos <- find_column(raw, idh_aliases)
  grade_text <- as.character(raw[[grade_pos]])
  if (is.na(idh_pos) && any(normalize_grade(grade_text) == "WHO 4", na.rm = TRUE) &&
      !any(!is.na(idh_from_grade_text(grade_text)))) {
    stop("TCGA clinical file lacks IDH status. Add column IDH_mutation_status ",
         "(e.g. Wildtype/Mutant) for WHO grade 4 classification.")
  }
  clinical <- data.frame(
    Patient = substr(as.character(raw[[id_pos]]), 1, 12),
    Grade = grade_text,
    IDH = if (!is.na(idh_pos)) as.character(raw[[idh_pos]]) else
      rep(NA_character_, nrow(raw)), stringsAsFactors = FALSE
  )
  check_unique_ids(clinical$Patient, "TCGA clinical patient IDs")
  clinical
}

prepare_tcga_violin <- function(scores, clinical) {
  require_columns(scores, c("id", score_column), "TCGA score file")
  check_unique_ids(scores$id, "TCGA score IDs")
  ids <- as.character(scores$id)
  is_tcga <- grepl("^TCGA-", ids, ignore.case = TRUE)
  is_gtex <- grepl("^GTEX[-_]", ids, ignore.case = TRUE)
  sample_type <- ifelse(is_tcga, substr(ids, 14, 15), NA_character_)
  is_normal <- (is_tcga & sample_type == "11") | is_gtex
  is_tumor <- is_tcga & sample_type %in% c("01", "02")
  idx <- match(substr(ids, 1, 12), clinical$Patient)
  grade <- clinical$Grade[idx]
  idh <- clinical$IDH[idx]
  group <- rep(NA_character_, length(ids))
  group[is_normal] <- "Normal"
  group[is_tumor] <- clinical_group(grade[is_tumor], idh[is_tumor])
  include <- !is.na(group) & is.finite(as_numeric_checked(scores[[score_column]], "TCGA score"))
  audit <- data.frame(id = ids, Sample_type = sample_type, Grade = grade,
                      IDH = idh, Group = group, Included = include,
                      stringsAsFactors = FALSE)
  write_tsv(audit, "TCGA_grade_assignment_audit.tsv")
  if (any(is_tumor & normalize_grade(grade) == "WHO 4" & is.na(group), na.rm = TRUE)) {
    warning("TCGA: unclassified WHO 4 samples excluded due to unknown IDH status")
  }
  data.frame(Cohort = "TCGA", Group = factor(group[include], levels = violin_levels),
             Score = as_numeric_checked(scores[[score_column]][include], "TCGA score"))
}
prepare_cgga_violin <- function(scores) {
  require_columns(scores, c("id", score_column, "Grade", "IDH"), "CGGA score file")
  groups <- clinical_group(scores$Grade, scores$IDH)
  numeric_score <- as_numeric_checked(scores[[score_column]], "CGGA score")
  include <- !is.na(groups) & is.finite(numeric_score)
  audit <- data.frame(id = scores$id, Grade = scores$Grade, IDH = scores$IDH,
                      Group = groups, Included = include, stringsAsFactors = FALSE)
  write_tsv(audit, "CGGA_grade_assignment_audit.tsv")
  if (any(normalize_grade(scores$Grade) == "WHO 4" & is.na(groups), na.rm = TRUE)) {
    warning("CGGA: unclassified WHO 4 samples excluded due to unknown IDH status")
  }
  data.frame(Cohort = "CGGA", Group = factor(groups[include], levels = violin_levels),
             Score = numeric_score[include])
}

save_plot_pair <- function(plot, stem, width = 5.5, height = 4.7) {
  ggplot2::ggsave(file.path(output_dir, paste0(stem, ".pdf")), plot,
                  width = width, height = height, units = "in", bg = "white")
  ggplot2::ggsave(file.path(output_dir, paste0(stem, ".tiff")), plot,
                  width = width, height = height, units = "in", dpi = 600,
                  compression = "lzw", bg = "white")
}
plot_violin <- function(data, title, stem) {
  if (nrow(data) == 0L) stop("No samples available for violin plot: ", title)
  p <- ggplot2::ggplot(data, ggplot2::aes(x = Group, y = Score, fill = Group)) +
    ggplot2::geom_violin(trim = FALSE, alpha = 0.75, color = "grey25") +
    ggplot2::geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
    ggplot2::labs(x = NULL, y = "Methionine metabolism score", title = title) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(legend.position = "none",
       axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
       plot.title = ggplot2::element_text(hjust = 0.5))
  save_plot_pair(p, stem, width = 7.8, height = 5.0)
}

format_p <- function(p) if (!is.finite(p)) "NA" else if (p < 0.001) "< 0.001" else
  paste0("= ", formatC(p, format = "f", digits = 3))
correlation_result <- function(data, x, y, cohort) {
  require_columns(data, c(x, y), cohort)
  a <- as_numeric_checked(data[[x]], paste0(cohort, ":", x))
  b <- as_numeric_checked(data[[y]], paste0(cohort, ":", y))
  keep <- is.finite(a) & is.finite(b)
  if (sum(keep) < 3L) stop(cohort, ": fewer than 3 complete pairs")
  if (stats::sd(a[keep]) == 0 || stats::sd(b[keep]) == 0) {
    stop(cohort, ": zero variance; Pearson correlation undefined")
  }
  ct <- stats::cor.test(a[keep], b[keep], method = "pearson")
  data.frame(Cohort = cohort, X = x, Y = y, Method = "pearson",
             N = sum(keep), R = unname(ct$estimate), P_value = ct$p.value)
}
scatter_plot <- function(data, x, y, title, stem) {
  tab <- correlation_result(data, x, y, title)
  df <- data.frame(x = as_numeric_checked(data[[x]], x),
                   y = as_numeric_checked(data[[y]], y))
  df <- df[is.finite(df$x) & is.finite(df$y), , drop = FALSE]
  lab <- sprintf("Pearson R = %.3f\nP %s\nn = %d",
                 tab$R, format_p(tab$P_value), tab$N)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_point(size = 2.4, alpha = 0.65, color = "#3E6C8E") +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
                         color = "#A4473D", fill = "#D9A7A1", linewidth = 0.8) +
    ggplot2::annotate("text", x = -Inf, y = Inf, label = lab,
                      hjust = -0.05, vjust = 1.1, size = 4) +
    ggplot2::labs(x = x, y = y, title = title) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  save_plot_pair(p, stem)
  tab
}

# S Fig. 1B: immune-gene expression and Pearson correlations.
make_immunosuppression_heatmap <- function(data) {
  require_columns(data, c("id", score_column), basename(tcga_heatmap_file))
  check_unique_ids(data$id, "TCGA heatmap IDs")
  genes <- setdiff(colnames(data), c("id", score_column))
  if (!length(genes) || anyDuplicated(genes)) {
    stop("Heatmap requires nonduplicated immune gene expression columns")
  }
  score <- as_numeric_checked(data[[score_column]], "Heatmap methionine score")
  if (length(score) < 3L || any(!is.finite(score)) || stats::sd(score) == 0) {
    stop("Heatmap score must be finite, nonconstant and have >=3 samples")
  }
  o <- order(score)
  score <- score[o]
  sorted <- data[o, , drop = FALSE]
  mat <- vapply(sorted[, genes, drop = FALSE], as_numeric_checked,
                numeric(nrow(sorted)), label = "immune expression")
  if (is.null(dim(mat))) mat <- matrix(mat, ncol = 1L)
  colnames(mat) <- genes
  rownames(mat) <- as.character(sorted$id)
  if (any(!is.finite(mat))) stop("Heatmap expression contains NA/Inf values")

  cor_table <- do.call(rbind, lapply(seq_along(genes), function(i) {
    v <- mat[, i]
    if (stats::sd(v) == 0) {
      return(data.frame(Gene = genes[i], R = NA_real_, P_value = NA_real_))
    }
    ct <- suppressWarnings(stats::cor.test(score, v, method = "pearson"))
    data.frame(Gene = genes[i], R = unname(ct$estimate), P_value = ct$p.value)
  }))
  cor_table$FDR <- stats::p.adjust(cor_table$P_value, method = "BH")
  cor_table$R_label <- ifelse(!is.na(cor_table$P_value) & cor_table$P_value < 0.05,
                               paste0("R=", signif(cor_table$R, 2)), "")
  write_tsv(cor_table, "SFig1B_immune_gene_correlations.tsv")

  expr <- t(mat)
  keep <- apply(expr, 1, stats::sd) > 0
  if (!any(keep)) stop("No variable immune genes for heatmap")
  if (any(!keep)) warning("Constant immune genes omitted from heatmap: ",
                          paste(rownames(expr)[!keep], collapse = ", "))
  expr <- expr[keep, , drop = FALSE]
  expr <- t(scale(t(expr)))
  mapped <- cor_table[match(rownames(expr), cor_table$Gene), , drop = FALSE]
  label_colors <- ifelse(mapped$R > 0, "firebrick3", "navy")
  label_colors[mapped$R_label == ""] <- "black"
  top_anno <- ComplexHeatmap::HeatmapAnnotation(
    Score = ComplexHeatmap::anno_barplot(score, gp = grid::gpar(fill = "grey40")),
    annotation_name_side = "left")
  left_anno <- ComplexHeatmap::rowAnnotation(
    Rvalue = ComplexHeatmap::anno_text(mapped$R_label, just = "left",
      gp = grid::gpar(fontsize = 10, col = label_colors)))
  ht <- ComplexHeatmap::Heatmap(expr, name = "Expression (z-score)",
    col = circlize::colorRamp2(c(-2.5, 0, 2.5),
                               c("navy", "#F0F4F8", "firebrick3")),
    cluster_rows = nrow(expr) > 1L, cluster_columns = FALSE,
    show_row_names = TRUE, show_column_names = FALSE,
    top_annotation = top_anno, left_annotation = left_anno)
  plot_height <- max(7, 0.18 * nrow(expr))
  save_ht <- function(ext) {
    out <- file.path(output_dir, paste0("SFig1B_immunosuppression_heatmap.", ext))
    if (ext == "pdf") {
      grDevices::pdf(out, width = 11, height = plot_height, useDingbats = FALSE)
    } else {
      grDevices::tiff(out, width = 11, height = plot_height, units = "in",
                      res = 600, compression = "lzw", bg = "white")
    }
    on.exit(grDevices::dev.off(), add = TRUE)
    ComplexHeatmap::draw(ht, heatmap_legend_side = "right",
                         annotation_legend_side = "right")
  }
  save_ht("pdf")
  save_ht("tiff")
  invisible(cor_table)
}

# Exact sample-ID join first, then an unambiguous patient-ID fallback only.
# This supplies ImmuneScore for the IDH-WT cohort without recomputing any score.
join_idhwt_immune <- function(idhwt, all_tcga) {
  check_unique_ids(idhwt$id, "TCGA IDH-WT IDs")
  check_unique_ids(all_tcga$id, "TCGA all-glioma IDs")
  require_columns(all_tcga, c("id", "ImmuneScore"), "TCGA all-glioma input")
  key <- match(idhwt$id, all_tcga$id)
  missing <- is.na(key)
  if (any(missing)) {
    p_all <- substr(as.character(all_tcga$id), 1, 12)
    p_wt <- substr(as.character(idhwt$id[missing]), 1, 12)
    # Only one all-glioma entry per patient permits a safe fallback.
    valid <- !duplicated(p_all) & !duplicated(p_all, fromLast = TRUE)
    candidates <- match(p_wt, p_all)
    candidates[is.na(candidates) | !valid[replace(candidates, is.na(candidates), 1L)]] <- NA_integer_
    key[missing] <- candidates
  }
  merged <- idhwt
  merged$ImmuneScore <- all_tcga$ImmuneScore[key]
  matched <- sum(!is.na(key))
  message("SFig1C TCGA IDH-WT: ImmuneScore matched for ", matched,
          "/", nrow(idhwt), " sample IDs")
  write_tsv(data.frame(id = idhwt$id, Matched_ImmuneScore = !is.na(key),
                       ImmuneScore = merged$ImmuneScore),
            "SFig1C_IDHWT_ImmuneScore_match_audit.tsv")
  merged
}

# For SFig11, Sample or id can identify rows. Clinical/Histology columns are
# retained if present. An external Gravendeel clinical table can be joined.
standardize_correlation_input <- function(raw, name) {
  id_idx <- find_column(raw, c("Sample", "id", "CGGA_ID", "sample_id"),
                        mandatory = TRUE, label = paste0(name, " sample ID"))
  expr_idx <- c(find_column(raw, "PPARG", mandatory = TRUE),
                find_column(raw, "CD274", mandatory = TRUE))
  out <- data.frame(Sample = as.character(raw[[id_idx]]),
                    PPARG = as_numeric_checked(raw[[expr_idx[1]]], "PPARG"),
                    CD274 = as_numeric_checked(raw[[expr_idx[2]]], "CD274"))
  check_unique_ids(out$Sample, name)
  for (col in c("Grade", "Histology", "IDH")) {
    aliases <- switch(col, Grade = grade_aliases,
                      Histology = c("Histology", "Histological_type", "Diagnosis"),
                      IDH = idh_aliases)
    pos <- find_column(raw, aliases)
    if (!is.na(pos)) out[[col]] <- as.character(raw[[pos]])
  }
  out
}
add_gravendeel_clinical <- function(expr, clinical) {
  id_idx <- find_column(clinical, c("Sample", "id", "sample_id"), mandatory = TRUE,
                        label = "Gravendeel clinical sample ID")
  ids <- as.character(clinical[[id_idx]])
  check_unique_ids(ids, "Gravendeel clinical IDs")
  idx <- match(expr$Sample, ids)
  message("Gravendeel: clinical metadata matched for ", sum(!is.na(idx)),
          "/", nrow(expr), " samples")
  for (col in c("Grade", "Histology", "IDH")) {
    aliases <- switch(col, Grade = grade_aliases,
                      Histology = c("Histology", "Histological_type", "Diagnosis"),
                      IDH = idh_aliases)
    pos <- find_column(clinical, aliases)
    if (!is.na(pos) && !(col %in% colnames(expr))) {
      expr[[col]] <- as.character(clinical[[pos]])[idx]
    }
  }
  expr
}
# All-glioma cohort should not include normal samples when pathology exists.
filter_all_glioma <- function(data) {
  if (!"Grade" %in% colnames(data) && !"Histology" %in% colnames(data)) {
    return(data)  # Input must already be filtered to glioma (see notes).
  }
  normal <- rep(FALSE, nrow(data))
  if ("Grade" %in% colnames(data)) {
    normal <- normal | (!is.na(normalize_grade(data$Grade)) &
                          normalize_grade(data$Grade) == "Normal")
  }
  if ("Histology" %in% colnames(data)) {
    normal <- normal | grepl("normal|control|non.?tumou?r|healthy",
                             as.character(data$Histology), ignore.case = TRUE)
  }
  data[!normal, , drop = FALSE]
}
# Conservative GBM subset: prefer explicit histological GBM; when grade 4 and
# IDH are available, include only IDH-WT. Never call IDH-mut grade 4 'GBM'.
filter_gbm <- function(data, label) {
  if (!"Grade" %in% colnames(data) && !"Histology" %in% colnames(data)) {
    return(list(data = NULL, reason = paste0(label, ": Grade/Histology is missing")))
  }
  grade4 <- if ("Grade" %in% colnames(data))
    normalize_grade(data$Grade) == "WHO 4" else rep(FALSE, nrow(data))
  hist_gbm <- if ("Histology" %in% colnames(data))
    grepl("GBM|glioblastoma", as.character(data$Histology), ignore.case = TRUE) else
      rep(FALSE, nrow(data))
  grade4[is.na(grade4)] <- FALSE
  hist_gbm[is.na(hist_gbm)] <- FALSE
  candidates <- grade4 | hist_gbm
  if ("IDH" %in% colnames(data)) {
    status <- normalize_idh(data$IDH)
    include <- candidates & status == "WT"
    include[is.na(include)] <- FALSE
    if (any(candidates & is.na(status))) {
      warning(label, ": grade-4 samples with unknown IDH excluded from GBM subset")
    }
  } else {
    # Without IDH status, histological GBM is allowed but WHO-4-only is not.
    include <- hist_gbm
    if (any(grade4 & !hist_gbm)) {
      warning(label, ": WHO 4 without IDH/GBM histology excluded")
    }
  }
  list(data = data[include, , drop = FALSE], reason = "")
}

# ---------------------------------------------------------------------------
# Read required files and normalize score headers.
# ---------------------------------------------------------------------------
required_files <- c(tcga_score_file, tcga_clinical_file, tcga_all_file,
                    tcga_idhwt_file, tcga_heatmap_file, cgga_score_file,
                    cgga_gbm_file, gravendeel_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) stop("Missing required files:\n", paste(missing_files, collapse = "\n"))

tcga_score <- normalize_score_columns(read_tsv(tcga_score_file), "TCGA score")
tcga_clinical <- clean_tcga_clinical(read_tsv(tcga_clinical_file))
tcga_all <- normalize_score_columns(read_tsv(tcga_all_file), "TCGA all-glioma")
tcga_idhwt <- normalize_score_columns(read_tsv(tcga_idhwt_file), "TCGA IDH-WT")
tcga_heatmap <- normalize_score_columns(read_tsv(tcga_heatmap_file), "TCGA heatmap")
cgga_score <- clean_cgga_score(read_tsv(cgga_score_file))
cgga_gbm <- standardize_correlation_input(read_tsv(cgga_gbm_file), "CGGA GBM")
gravendeel <- standardize_correlation_input(read_tsv(gravendeel_file), "Gravendeel")

for (x in list(list(tcga_all, 701, "TCGA glioma"),
               list(tcga_idhwt, 149, "TCGA IDH-WT GBM"),
               list(tcga_heatmap, 149, "TCGA immune heatmap"),
               list(cgga_gbm, 139, "CGGA GBM"),
               list(gravendeel, 276, "Gravendeel"))) {
  check_expected_n(x[[1]], x[[2]], x[[3]])
}
require_columns(tcga_all, c("id", score_column, "ImmuneScore", "PPARG", "CD274"),
                "TCGA all-glioma")
require_columns(tcga_idhwt, c("id", score_column, "PPARG", "CD274"),
                "TCGA IDH-WT GBM")
check_unique_ids(tcga_all$id, "TCGA all-glioma")
check_unique_ids(tcga_idhwt$id, "TCGA IDH-WT GBM")
if (!setequal(tcga_idhwt$id, tcga_heatmap$id)) {
  stop("TCGA IDH-WT score and immunosuppression heatmap files have different sample IDs")
}
if ("Histology" %in% colnames(cgga_gbm)) {
  good <- grepl("GBM|glioblastoma", as.character(cgga_gbm$Histology), ignore.case = TRUE)
  if (any(!good)) stop("CGGA_GBM_PPARG_CD274.tsv contains non-GBM/missing Histology")
}

# S Fig. 1A: cohort-specific violin plots with WHO grade-4 IDH stratification.
tcga_violin <- prepare_tcga_violin(tcga_score, tcga_clinical)
cgga_violin <- prepare_cgga_violin(cgga_score)
plot_violin(tcga_violin, "TCGA / GTEx", "SFig1A_methionine_score_violin_TCGA")
plot_violin(cgga_violin, "CGGA", "SFig1A_methionine_score_violin_CGGA")
violin_data <- rbind(tcga_violin, cgga_violin)
write_tsv(violin_data, "SFig1A_violin_plot_data.tsv")
counts <- as.data.frame(table(violin_data$Cohort, violin_data$Group))
colnames(counts) <- c("Cohort", "Group", "N")
write_tsv(counts, "SFig1A_group_counts.tsv")
print(counts)

# S Fig. 1B: immunosuppression-associated gene heatmap.
make_immunosuppression_heatmap(tcga_heatmap)

# S Figs. 1C and 6C: all-glioma and IDH-wildtype GBM correlations.
correlation_tables <- list()
correlation_tables[["SFig1C_TCGA_all"]] <- scatter_plot(
  tcga_all, score_column, "ImmuneScore", "TCGA glioma",
  "SFig1C_methionine_score_vs_ImmuneScore_TCGA_all_glioma")
tcga_idhwt_immune <- join_idhwt_immune(tcga_idhwt, tcga_all)
correlation_tables[["SFig1C_TCGA_IDHWT"]] <- scatter_plot(
  tcga_idhwt_immune, score_column, "ImmuneScore", "TCGA IDH-WT GBM",
  "SFig1C_methionine_score_vs_ImmuneScore_TCGA_IDHWT")
correlation_tables[["SFig6C_TCGA_all"]] <- scatter_plot(
  tcga_all, score_column, "PPARG", "TCGA glioma",
  "SFig6C_PPARG_vs_methionine_score_TCGA_all_glioma")
correlation_tables[["SFig6C_TCGA_IDHWT"]] <- scatter_plot(
  tcga_idhwt, score_column, "PPARG", "TCGA IDH-WT GBM",
  "SFig6C_PPARG_vs_methionine_score_TCGA_IDHWT")

# S Fig. 11A: PPARG-CD274 correlations across cohorts and disease subsets.
# Plot availability is recorded in the status table.
sfig11_status <- list()
run_sfig11 <- function(key, data, label, stem, missing_reason = "") {
  if (is.null(data)) {
    sfig11_status[[key]] <<- data.frame(Plot = key, Status = "SKIPPED",
      N = 0L, Reason = missing_reason)
    warning("SFig11A ", key, " SKIPPED: ", missing_reason)
    return(invisible(NULL))
  }
  result <- tryCatch(scatter_plot(data, "PPARG", "CD274", label, stem),
                     error = function(e) e)
  if (inherits(result, "error")) {
    sfig11_status[[key]] <<- data.frame(Plot = key, Status = "SKIPPED",
      N = nrow(data), Reason = conditionMessage(result))
    warning("SFig11A ", key, " SKIPPED: ", conditionMessage(result))
  } else {
    correlation_tables[[paste0("SFig11A_", key)]] <<- result
    sfig11_status[[key]] <<- data.frame(Plot = key, Status = "CREATED",
      N = result$N, Reason = "")
  }
  invisible(NULL)
}
run_sfig11("TCGA_all_glioma", tcga_all, "TCGA glioma",
           "SFig11A_PPARG_vs_CD274_TCGA_all_glioma")
run_sfig11("TCGA_GBM_IDHWT", tcga_idhwt, "TCGA GBM IDH-WT",
           "SFig11A_PPARG_vs_CD274_TCGA_GBM_IDHWT")

if (file.exists(cgga_all_file)) {
  cgga_all <- standardize_correlation_input(read_tsv(cgga_all_file), "CGGA glioma")
  cgga_all <- filter_all_glioma(cgga_all)
  run_sfig11("CGGA_all_glioma", cgga_all, "CGGA glioma",
             "SFig11A_PPARG_vs_CD274_CGGA_all_glioma")
} else {
  run_sfig11("CGGA_all_glioma", NULL, "CGGA glioma", "",
    "Need data/CGGA_all_glioma_PPARG_CD274.tsv: Sample, PPARG, CD274; ideally Grade/Histology")
}
run_sfig11("CGGA_GBM", cgga_gbm, "CGGA GBM",
           "SFig11A_PPARG_vs_CD274_CGGA_GBM")

if (file.exists(gravendeel_clinical_file)) {
  gravendeel <- add_gravendeel_clinical(gravendeel,
                                      read_tsv(gravendeel_clinical_file))
}
# When Grade/Histology is unavailable, the Gravendeel expression table is
# expected to contain glioma samples only.
run_sfig11("Gravendeel_all_glioma", filter_all_glioma(gravendeel),
           "Gravendeel glioma", "SFig11A_PPARG_vs_CD274_Gravendeel_all_glioma")
grav_gbm <- filter_gbm(gravendeel, "Gravendeel")
run_sfig11("Gravendeel_GBM", grav_gbm$data, "Gravendeel GBM",
           "SFig11A_PPARG_vs_CD274_Gravendeel_GBM", grav_gbm$reason)

# Save audit/status and all successfully calculated statistics.
sfig11_report <- do.call(rbind, sfig11_status)
write_tsv(sfig11_report, "SFig11A_plot_status.tsv")
write_tsv(do.call(rbind, correlation_tables), "bulk_correlation_statistics.tsv")
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
message("Finished. Results in: ", output_dir)
message("SFig11A plots created: ", sum(sfig11_report$Status == "CREATED"), "/6")
