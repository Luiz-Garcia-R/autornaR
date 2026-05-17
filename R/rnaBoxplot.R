#' Gene-level comparison and boxplot visualization using limma
#'
#' @description
#' Performs differential expression analysis for a single gene using a linear
#' model fitted with \pkg{limma}, optionally including batch effects.
#' The function computes log2 fold-change (log2FC), 95% confidence interval,
#' p-value, and Cohen's d, and generates a combined violin + boxplot
#' visualization with individual sample points.
#'
#' Gene identifiers can be provided as Ensembl (with or without version),
#' ENTREZID, or gene symbol. Annotation is performed automatically using
#' organism-specific \pkg{org.*.eg.db} databases.
#'
#' @param project \code{rna_project} object created by \code{rna.project()}.
#' @param gene Character. Single gene identifier (Ensembl ID, ENTREZID, or SYMBOL).
#' @param group_col Character. metadata column defining experimental groups.
#' @param contrast Specifies the comparison. Can be:
#'   \itemize{
#'     \item A string: \code{"test_vs_reference"} matching a stored comparison
#'     \item A character vector of length 2: \code{c("test", "reference")}
#'     \item \code{"last"} to use the most recent comparison in the project
#'     \item \code{NULL}: automatically inferred if only two groups are present
#'   }
#' @param batch_col Character. Optional metadata column defining batch effects
#'   to include in the design matrix.
#' @param style Character. Whether use \code{boxplot} or \code{violin} for data
#'   visualization.
#' @param show_signif Logical. Whether to display significance annotation
#'   (asterisks based on p-value) on the plot.
#' @param seed Optional. Set a seed for reproducibility for point \code{jitter}.
#' @param save Logical. Whether to store results in
#'   \code{rna_project$analyses$boxplot}.
#' @param verbose Logical. Whether to print informative messages.
#'
#' @return
#' Invisibly returns a list containing:
#' \itemize{
#'   \item \code{gene}: Gene identifier used in the expression matrix.
#'   \item \code{gene_label}: Gene symbol (if available).
#'   \item \code{contrast}: Comparison vector (\code{test}, \code{reference}).
#'   \item \code{groups}: Named list with \code{test} and \code{reference}.
#'   \item \code{statistics}: List with:
#'     \itemize{
#'       \item \code{logFC}: log2 fold-change (test − reference)
#'       \item \code{CI_low}, \code{CI_high}: 95% confidence interval
#'       \item \code{p_value}: p-value from limma
#'       \item \code{cohen_d}: Effect size
#'       \item \code{scale}: Expression scale used
#'     }
#'   \item \code{design}: Model formula used
#'   \item \code{method}: Statistical method ("limma")
#' }
#'
#' This function is typically used after \code{rna.normalize()} and before
#' downstream analyses such as \code{rna.compare()}, or clustering.
#'
#' @details
#' The contrast is always defined as \strong{test − reference}, which determines
#' the sign of the log2 fold-change. However, for visualization purposes,
#' groups are displayed in the order \strong{reference → test} (left to right)
#' to improve interpretability.
#'
#' The plot combines:
#' \itemize{
#'   \item Violin plot (distribution)
#'   \item Boxplot (median and IQR)
#'   \item Jittered points (individual samples)
#'   \item Optional significance annotation
#' }
#'
#' Results can be stored in the active \code{rna_project} for reproducibility.
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' rna.boxplot(project = my_project,
#'             gene = "TP53")
#'
#' # Using stored comparison
#' rna.boxplot(my_project,
#'             gene = "TP53",
#'             contrast = "treated_vs_control")
#'
#' # Changing style
#' rna.boxplot(my_project,
#'             gene = "TP53",
#'             style = "violin"
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_violin geom_boxplot geom_point
#'   annotate labs theme_minimal theme element_text
#' @importFrom dplyr left_join
#' @importFrom stats model.matrix qt
#'
#' @export

rna.boxplot <- function(project,
                        gene,
                        group_col = "Group",
                        contrast = NULL,
                        batch_col = NULL,
                        style = c("boxplot", "violin"),
                        show_signif = TRUE,
                        seed = NULL,
                        save = TRUE,
                        verbose = TRUE
) {

  style <- match.arg(style)

  # ---------------------------
  # 1) Required packages
  # ---------------------------
  .check_dependencies(c("ggplot2","dplyr"))
  .check_dependencies(c("limma","AnnotationDbi", bioc = TRUE))

  # ---------------------------
  # 2) Get active project
  # ---------------------------
  proj <- project

  expr_mat <- as.matrix(.get_expr(proj))
  metadata <- .get_meta(proj)
  norm_method <- .get_norm_method(proj)

  gene_map <- .get_gene_annotation(proj)
  gene_map <- .align_gene_annotation(gene_map, expr_mat)

  # ---------------------------
  # 4) Resolve gene ID
  # ---------------------------
  gene_ids <- gene_map$gene_id
  gene_symbols <- gene_map$symbol

  valid_idx <- which(
    gene_map$gene_id %in% gene |
      gene_map$symbol %in% gene
  )

  if (length(valid_idx) == 0) {
    stop("Gene not found in expression matrix or annotation.")
  }

  gene_use <- gene_ids[valid_idx[1]]

  gene_label <- gene_symbols[valid_idx[1]]
  if (is.na(gene_label) || gene_label == "") {
    gene_label <- gene_use
  }

  # ---------------------------
  # 4.1) Build design + contrast
  # ---------------------------

  if (!(group_col %in% colnames(metadata))) {
    stop(paste0("Column '", group_col, "' not found in metadata."))
  }

  # group factor
  metadata[[group_col]] <- as.factor(metadata[[group_col]])

  # ---------------------------
  # 4.2) Resolve contrast
  # ---------------------------
  if (is.null(contrast)) {

    groups <- levels(metadata[[group_col]])

    if (length(groups) != 2) {
      stop("contrast must be specified when there are more than 2 groups.")
    }

    contrast_vec <- c(groups[2], groups[1])  # test vs reference

  } else if (is.character(contrast) && length(contrast) == 1) {

    if (contrast == "last") {
      stop("'last' not implemented yet.")
    }

    contrast_vec <- strsplit(contrast, "_vs_")[[1]]

  } else if (is.character(contrast) && length(contrast) == 2) {

    contrast_vec <- contrast

  } else {
    stop("Invalid contrast format.")
  }

  # Ensure factor order
  metadata[[group_col]] <- factor(metadata[[group_col]],
                              levels = c(contrast_vec[2], contrast_vec[1]))

  # ---------------------------
  # 4.3) Design formula
  # ---------------------------
  if (!is.null(batch_col)) {

    if (!(batch_col %in% colnames(metadata))) {
      stop(paste0("Column '", batch_col, "' not found in metadata."))
    }

    design_formula <- as.formula(
      paste0("~ ", batch_col, " + ", group_col)
    )

  } else {
    design_formula <- as.formula(
      paste0("~ 0 + ", group_col)
    )
  }

  design <- model.matrix(design_formula, data = metadata)

  # ---------------------------
  # 4.4) Contrast matrix
  # ---------------------------
  contrast_name <- paste0(
    group_col, contrast_vec[1], " - ", group_col, contrast_vec[2]
  )

  contrast_matrix <- limma::makeContrasts(
    contrasts = contrast_name,
    levels = design
  )

  # ---------------------------
  # 5) Get limma result (cache or compute)
  # ---------------------------
  limma_id <- paste(gene_use, group_col, sep = "_")
  is_log <- norm_method %in% c("log2", "rlog", "vst")

  expr_mat_use <- if (is_log) expr_mat else log2(expr_mat + 1)

  if (save && !is.null(proj$analyses$gene_tests) &&
      limma_id %in% names(proj$analyses$gene_tests)) {

    if (verbose)
      message("[rna.boxplot] Limma result already saved. Using cached result.")

    cached <- proj$analyses$gene_tests[[limma_id]]

    gene_res <- cached$limma_res

    if (is.null(cached$scale) || cached$scale != "log2") {

      ci_low  <- as.numeric(cached$ci_low)
      ci_high <- as.numeric(cached$ci_high)

    } else {

      if (verbose)
        message("[rna.boxplot] Old cache detected. Recomputing limma fit for CI.")

      fit <- limma::lmFit(expr_mat_use, design)
      fit <- limma::contrasts.fit(fit, contrast_matrix)
      fit <- limma::eBayes(fit)

      tt <- limma::topTable(
        fit,
        coef = 1,
        number = Inf,
        sort.by = "none",
        confint = TRUE
      )

      gene_res <- tt[gene_use, ]

      ci_low  <- gene_res$CI.L
      ci_high <- gene_res$CI.R
    }

  } else {

    if (verbose)
      message("[rna.boxplot] No cached result. Running limma fit.")

    if (!is_log) {
      expr_mat_use <- log2(expr_mat + 1)
    } else {
      expr_mat_use <- expr_mat
    }

    fit <- limma::lmFit(expr_mat_use, design)
    fit <- limma::contrasts.fit(fit, contrast_matrix)
    fit <- limma::eBayes(fit)

    tt <- limma::topTable(
      fit,
      coef = 1,
      number = Inf,
      sort.by = "none",
      confint = TRUE
    )

    gene_res <- tt[gene_use, ]

    ci_low  <- gene_res$CI.L
    ci_high <- gene_res$CI.R
  }

  # ---------------------------
  # 6) Effect size (Cohen's d) and IC95
  # ---------------------------
  value_vec <- expr_mat_use[gene_use, ]

  df_sub <- data.frame(
    Sample = colnames(expr_mat),
    Value  = as.numeric(value_vec)
  )

  df_sub <- dplyr::left_join(df_sub, metadata, by = "Sample")
  df_sub <- df_sub[df_sub[[group_col]] %in% contrast_vec, ]

  groups_split <- split(df_sub$Value, df_sub[[group_col]])

  if (any(!contrast_vec %in% names(groups_split))) {
    stop("One or more groups in contrast have no samples after filtering.")
  }

  x <- groups_split[[contrast_vec[1]]]
  y_group <- groups_split[[contrast_vec[2]]]

  nx <- sum(!is.na(x)); ny <- sum(!is.na(y_group))
  mean_diff <- mean(x, na.rm = TRUE) - mean(y_group, na.rm = TRUE)
  sd_pooled <- sqrt(((nx - 1) * sd(x, na.rm = TRUE)^2 + (ny - 1) * sd(y_group, na.rm = TRUE)^2)/(nx + ny - 2))
  cohen_d <- mean_diff / sd_pooled
  se_diff <- sd_pooled * sqrt(1/nx + 1/ny)
  ic95 <- mean_diff + c(-1,1) * qt(0.975, df = nx+ny-2) * se_diff

  # ---------------------------
  # 7) P-value annotation
  # ---------------------------
  pval <- gene_res$P.Value
  signif_label <- if (pval < 0.001) "***"
  else if (pval < 0.01) "**"
  else if (pval < 0.05) "*"
  else "ns"

  p_label <- if (pval < 0.001) {
    "p < 0.001"
  } else {
    paste0("p = ", formatC(pval, format = "f", digits = 3))
  }

  # ---------------------------
  # 8) Plot
  # ---------------------------
  df_long <- data.frame(
    Sample = colnames(expr_mat),
    Value  = as.numeric(value_vec)
  )

  df_long <- dplyr::left_join(df_long, metadata, by = "Sample")
  df_long <- df_long[df_long[[group_col]] %in% contrast_vec, ]
  df_long[[group_col]] <- droplevels(df_long[[group_col]])


  # Boxplot Plot
  if (style == "boxplot") {

  p <- ggplot2::ggplot(
    df_long,
    ggplot2::aes(x = .data[[group_col]],
                 y = Value, fill = .data[[group_col]])) +

    ggplot2::geom_boxplot(
      alpha = 0.75,
      outlier.shape = NA,
      width = 0.7,
      linewidth = 0.7) +

    ggplot2::geom_point(
      ggplot2::aes(
        color = .data[[group_col]]),
      position = ggplot2::position_jitter(width = 0.1),
      shape = 21,
      color = "black",
      size = 1.2,
      alpha = 0.2) +

    ggplot2::labs(
      title = gene_label,
      subtitle = paste0("log2FC = ", round(gene_res$logFC, 3), " | ", p_label),
      x = "",
      y = "Normalized expression") +

    ggplot2::theme_bw(base_size = 12) +

    ggplot2::theme(
      legend.position = "none",
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 12))
    )

  if (show_signif) {
    p <- p + ggplot2::annotation_custom(
      grid::textGrob(
        label = signif_label,
        x = 0.5,
        y = 0.92,
        gp = grid::gpar(col = "grey20", fontsize = 16)
      )
    )
  }

    # Violin Plot
  } else if (style == "violin") {

    p <- ggplot2::ggplot(
      df_long,
      ggplot2::aes(x = .data[[group_col]], y = Value, fill = .data[[group_col]])
    ) +
      ggplot2::geom_violin(trim = FALSE, alpha = 0.5, color = NA) +
      ggplot2::geom_boxplot(width = 0.2, outlier.shape = NA, linewidth = 0.4) +
      ggplot2::geom_point(
        ggplot2::aes(color = .data[[group_col]]),
        position = ggplot2::position_jitter(width = 0.1),
        shape = 21,
        color = "black",
        size = 1.2,
        alpha = 0.2
      ) +
      ggplot2::labs(
        title = gene_label,
        subtitle = paste0(
          "log2FC = ", round(gene_res$logFC, 3),
          " | ", p_label
        ),
        x = "",
        y = "Normalized expression"
      ) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(
        legend.position = "none",
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 12))
      )

    if (show_signif) {
      p <- p + ggplot2::annotation_custom(
        grid::textGrob(
          label = signif_label,
          x = 0.5,
          y = 0.92,
          gp = grid::gpar(col = "grey20", fontsize = 16)
          )
        )
      }
    }

  # ---------------------------
  # 9) Output
  # ---------------------------
  rng_state <- if (!is.null(seed)) .Random.seed else NULL

  params <- list(
    timestamp = Sys.time(),
    gene = gene_use,
    gene_label = gene_label,
    contrast = contrast_vec,
    groups = list(
      reference = contrast_vec[2],
      test = contrast_vec[1]
    ),
    seed = seed,
    rng_state = rng_state
  )

  obj <- list(
    params = params,
    statistics = list(
      logFC = gene_res$logFC,
      CI_low = ci_low,
      CI_high = ci_high,
      p_value = pval,
      cohen_d = cohen_d,
      scale = ifelse(is_log, "log2", "log2_transformed")
    ),
    design = deparse(design_formula),
    method = "limma"
  )

  comparison_id <- paste(
    gene_use,
    contrast_vec[1],
    "vs",
    contrast_vec[2],
    sep = "_"
  )

  # ---------------------------
  # 10) Attach to project
  # ---------------------------
  if (save) {

    if (!is.null(proj$analyses$boxplot[[comparison_id]])) {
      warning("Overwriting existing boxplot with same ID.")
    }

    proj <- .attach_to_project(
      proj,
      obj,
      slot = "analyses",
      subtype = "boxplot",
      prefix = "boxplot",
      id = comparison_id,
      log = list(
        gene = gene_use,
        contrast = contrast_vec,
        group_col = group_col,
        design = deparse(design_formula),
        method = "limma",
        signature = list(
          gene = gene_use,
          contrast = contrast_vec,
          group_col = group_col,
          design = deparse(design_formula)
        )
      )
    )
  }

  print(p)

  # ---------------------------
  # 11) Return
  # ---------------------------
  .print_header("RNA Boxplot results")

  .print_block("Results Summary", function() {
    cat("Gene:               ", gene_use, "\n")
    cat("Gene label:         ", gene_label, "\n")
    cat("Contrast:           ", contrast_vec, "\n")
    cat("CI high:            ", ci_high, "\n")
    cat("CI low:             ", ci_low, "\n")
    cat("Cohen d:            ", cohen_d, "\n")
    cat("p-value:            ", pval, "\n")
  })

  return(invisible(proj))
}
