#' Correlation analysis between RNA-seq samples or groups
#'
#' Computes correlation between two RNA-seq samples or experimental groups
#' using Pearson, Spearman, or Kendall methods. When \code{method = "auto"},
#' the function automatically selects the most appropriate method based on
#' normality, presence of ties, sample size, and outlier fraction.
#'
#' The function requires an active \code{rna_project} object containing
#' normalized expression data and metadata.
#'
#' When mode = "gene", genes can be provided as Ensembl IDs or gene symbols.
#' The function automatically resolves identifiers and displays gene symbols
#' when available.
#'
#'
#' @param project \code{rna_project} object created by \code{rna.project()}.
#' @param method Character. Correlation method to use.
#'   One of:
#'   \itemize{
#'     \item `"auto"` (default): automatic method selection
#'     \item `"pearson"`
#'     \item `"spearman"`
#'     \item `"kendall"`
#'   }
#' @param mode Character. Correlation structure:
#'   \itemize{
#'     \item `"group"`: correlation between experimental groups
#'     \item `"sample"`: correlation between individual samples
#'     \item `"gene"`: correlation between selected genes
#'   }
#' @param group Optional. Subsets the data to a specific experimental group when mode = "gene".
#' @param genes Character vector of length 2 with both gene symbol or ENSG.
#' @param compute_mi Logical. If \code{TRUE}, computes mutual information (MI)
#'   between variables to capture non-linear dependencies not detected by
#'   correlation coefficients.
#' @param mi_method Character. Method used to estimate mutual information.
#'   One of:
#'   \itemize{
#'     \item `"discrete"`: discretization-based estimator
#'     \item `"knn"`: k-nearest neighbors estimator (continuous data, default)
#'   }
#' @param sample_x Character. Name of first sample (required when \code{mode = "sample"}).
#' @param sample_y Character. Name of second sample (required when \code{mode = "sample"}).
#' @param save Logical. Whether to store results in the active \code{rna_project}.
#' @param verbose Logical. If \code{TRUE} (default), prints selected method (when auto)
#'   and displays the correlation plot.
#'
#' @details
#' Automatic method selection follows these rules:
#'
#' \itemize{
#'   \item Extreme outlier fraction (> 10%) -> Kendall
#'   \item Moderate outliers (> 3%) or presence of ties -> Spearman
#'   \item Large sample size (> 5000 observations) -> Pearson
#'   \item Approximate normality (Shapiro–Wilk) -> Pearson
#'   \item Otherwise -> Spearman
#' }
#'
#' Outliers are detected using the IQR rule (1.5 × IQR).
#'
#' The smoothing line in the plot depends on the correlation method:
#' linear regression is used for Pearson, while loess smoothing is used
#' for rank-based methods (Spearman/Kendall) to capture non-linear trends.
#'
#' In addition to correlation, the function can compute Mutual Information (MI),
#' which captures both linear and non-linear dependencies between variables.
#'
#' When MI is enabled, the function also provides:
#' \itemize{
#'   \item `non_linear_signal` — difference between MI and absolute correlation,
#'         indicating potential non-linear structure
#'   \item `relation_type` — heuristic classification of the relationship:
#'     \itemize{
#'       \item "strong linear relationship"
#'       \item "monotonic (possibly non-linear)"
#'       \item "non-linear complex relationship"
#'       \item "weak or ambiguous relationship"
#'       \item "no clear association"
#'     }
#' }
#'
#' Note: Mutual information values are not directly comparable to correlation
#' coefficients and depend on the estimation method.
#'
#' These metrics are exploratory and should be interpreted within the
#' context of the data and experimental design.
#'
#' The function stores results inside the active project under:
#' \code{rna_project$analyses$correlation}.
#'
#' @return
#' An object of class \code{"rna_correlation"} containing:
#' \itemize{
#'   \item `method` — Method used for correlation
#'   \item `mode` — Correlation mode ("group", "sample" or "gene")
#'   \item `estimate` — Correlation coefficient
#'   \item `p.value` — P-value from \code{cor.test}
#'   \item `mi`, `mi_norm`, `mi_method` — Mutual information metrics (if computed)
#'   \item `relation_type` — Heuristic classification of dependency structure
#'   \item `conf.int` — Confidence interval
#'   \item `interpretation` — Strength classification
#'   \item `diagnostics` — Method selection diagnostics
#'   \item `structure` — Structural information about compared entities
#'   \item `plot` — ggplot object
#' }
#'
#' @examples
#' \dontrun{
#' # Group correlation (automatic method selection)
#' rna.corr(project = my_project,
#'          mode = "group")
#'
#' # Sample-to-sample correlation
#' rna.corr(my_project,
#'          mode = "sample",
#'          sample_x = "Control_1",
#'          sample_y = "Treatment_2"
#' )
#'
#' # Gene-to-gene correlation
#' rna.corr(my_project,
#'          mode = "gene",
#'          genes = c("gene1", "gene2"))
#' }
#'
#' @importFrom stats cor.test shapiro.test quantile
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth theme_minimal labs
#'
#' @export

rna.corr <- function(project,
                     method = c("auto", "pearson", "spearman", "kendall"),
                     mode = c("sample", "group", "gene"),
                     group = NULL,
                     genes = NULL,
                     compute_mi = TRUE,
                     mi_method = c("discrete", "knn"),
                     sample_x = NULL,
                     sample_y = NULL,
                     save = TRUE,
                     verbose = TRUE
) {

  method <- match.arg(method)
  mode <- match.arg(mode)
  mi_method <- match.arg(mi_method)

  # ---------------------------
  # 0) Basic checks
  # ---------------------------
  pkgs <- c("FNN", "ggplot2")

  .check_dependencies(pkgs)

  # ---------------------------
  # 1) Get active project
  # ---------------------------
  proj <- project

  expr_mat <- as.matrix(.get_expr(proj))
  metadata <- .get_meta(proj)

  # ---------------------------
  # 2) Compute gene correlation
  # ---------------------------
  if (mode == "gene") {

    if (is.null(genes) || length(genes) != 2) {
      stop("Provide exactly two genes.")
    }

    # mapping
    gene_map <- .get_gene_annotation(proj, expr_mat)

    gene_ids <- gene_map$gene_id
    gene_symbols <- gene_map$symbol

    gene_lookup <- setNames(gene_ids, gene_symbols)
    gene_lookup <- gene_lookup[!duplicated(names(gene_lookup))]

    if (any(duplicated(gene_symbols))) {
      warning("Duplicated gene symbols detected. Using first occurrence.")
    }

    resolve_gene <- function(g) {
      if (g %in% gene_symbols) {
        id <- gene_lookup[g]
        label <- g
      } else if (g %in% gene_ids) {
        id <- g

        # Find Symbol
        idx <- match(g, gene_ids)
        label <- gene_symbols[idx]

        if (is.na(label) || label == "") {
          label <- g
        }
      } else {
        stop("Gene not found: ", g)
      }

      return(list(id = id, label = label))
    }

    g1 <- resolve_gene(genes[1])
    g2 <- resolve_gene(genes[2])

    clean_ids <- function(x) sub("\\..*", "", x)
    expr_ids <- clean_ids(rownames(expr_mat))

    if (!all(clean_ids(c(g1$id, g2$id)) %in% expr_ids)) {
      stop("Genes not found in expression matrix.")
    }

    # subset by group
    if (!is.null(group) && !group %in% metadata$Group) {
      stop("Group not found: ", group)
    }

    if (!is.null(group)) {
      idx <- metadata$Group == group

      if (sum(idx) < 3) {
        stop("Selected group has fewer than 3 samples.")
      }

      expr_sub <- expr_mat[, idx, drop = FALSE]
    } else {
      expr_sub <- expr_mat
    }

    match_idx <- match(clean_ids(g1$id), expr_ids)
    x <- expr_mat[match_idx, ]

    match_idx <- match(clean_ids(g2$id), expr_ids)
    y <- expr_mat[match_idx, ]

    labels <- c(g1$label, g2$label)

    structure_info <- list(
      level = "gene",
      genes = genes,
      group = group
    )

  } else {

    data_obj <- .prepare_corr_data(
      expr_mat,
      metadata,
      mode,
      sample_x = sample_x,
      sample_y = sample_y
    )

    x <- data_obj$x
    y <- data_obj$y
    labels <- data_obj$labels
    structure_info <- data_obj$structure
  }

  structure_info$n_obs <- length(x)

  if (length(x) < 3) {
    stop("Not enough samples to compute correlation.")
  }

  data_df <- data.frame(x = x, y = y)

  # ---------------------------
  # 5) Method diagnostics
  # ---------------------------
  diag_obj <- .run_corr_diagnostics(x, y, method)
  method_used <- diag_obj$method_used
  diagnostics <- diag_obj$diagnostics

  if (verbose && method == "auto")
    message("Selected method: ", method_used)

  # ---------------------------
  # 6) Statistical test
  # ---------------------------
  test <- suppressWarnings(
    stats::cor.test(x, y, method = method_used)
  )

  estimate <- round(unname(test$estimate), 3)

  r_abs <- abs(estimate)

  interpretation <- if (r_abs < 0.80) {
    "weak similarity"
  } else if (r_abs < 0.90) {
    "moderate similarity"
  } else if (r_abs < 0.95) {
    "strong similarity"
  } else if (r_abs < 0.98) {
    "very strong similarity"
  } else {
    "near-identical profiles"
  }

  if (mode == "gene") {
    interpretation <- paste0(interpretation, " (gene-level)")
  }

  # ---------------------------
  # 7) Compute MI
  # ---------------------------
  mi_value <- NULL
  mi_norm  <- NULL

  if (compute_mi) {

    mi_value <- .compute_mi(x, y, method = mi_method)

    if (mi_method == "discrete") {
      mi_norm <- .normalize_mi(mi_value, x, y, method = "discrete")
    }
  }

  non_linear_signal <- if (!is.null(mi_norm)) {
    mi_norm - abs(estimate)
  } else if (!is.null(mi_value)) {
    mi_value - abs(estimate)
  } else {
    NULL
  }

  mi_signal <- if (!is.null(mi_norm)) mi_norm else mi_value

  relation_type <- .classify_relationship(
    r = estimate,
    mi_signal = mi_signal,
    non_linear_signal = non_linear_signal
  )

  # ---------------------------
  # 8) Plot
  # ---------------------------
  smooth_method <- if (method_used == "pearson") "lm" else "loess"

  line_col <- "#2F5D8A"
  ci_fill  <- "#2F5D8A"
  point_col <- "grey60"

  if (mode == "gene") {

    labs_obj <- ggplot2::labs(
      title = paste0(mode, " correlation (", method_used, ")"),
      subtitle = paste0(
        labels[1], " vs ", labels[2],
        if (!is.null(group)) paste0(" (group: ", group, ")") else "",
        "\n",
        "r = ", estimate,
        " | p = ", signif(test$p.value, 3)
      ),
      x = labels[1],
      y = labels[2]
    )

  } else {

    labs_obj <- ggplot2::labs(
      title = paste0(mode, " correlation (", method_used, ")"),
      subtitle = paste0(
        "r = ", estimate,
        " | p = ", signif(test$p.value, 3),
        " | ", interpretation
      ),
      x = labels[1],
      y = labels[2]
    )
  }

  g <- ggplot2::ggplot(data_df, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_point(alpha = 0.25, size = 2.3, color = point_col) +
    ggplot2::geom_smooth(
      method = smooth_method,
      formula = if (smooth_method == "lm") (y ~ x) else NULL,
      span = if (smooth_method == "loess") 0.75 else NULL,
      se = TRUE,
      level = 0.95,
      color = line_col,
      fill = ci_fill,
      alpha = 0.25,
      linewidth = 1.1
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    labs_obj

  if (verbose) print(g)

  # ---------------------------
  # 9) Output
  # ---------------------------
  obj <- list(
    timestamp = Sys.time(),
    method = method_used,
    mode = mode,
    estimate = estimate,
    p.value = test$p.value,
    structure = structure_info,
    mi = mi_value,
    mi_norm = mi_norm,
    mi_method = mi_method,
    relation_type = relation_type,
    conf.int = test$conf.int,
    interpretation = interpretation,
    diagnostics = diagnostics,
    plot = g
  )

  class(obj) <- "rna_correlation"

  # ---------------------------
  # 10) Attach to project
  # ---------------------------
  if (save) {

    proj <- .attach_to_project(
      proj,
      obj,
      slot = "analyses",
      subtype = "correlation",
      prefix = "corr",
      log = list(
        mode = mode,
        compared = paste(labels, collapse = " vs "),
        method = method_used,
        n = length(x),
        r = unname(estimate),
        p = signif(test$p.value, 3)
      )
    )
  }

  # ---------------------------
  # 11) Return
  # ---------------------------
  if (verbose) {

    .print_header("RNA Correlation")

    .print_block("Summary", function() {
      cat("Mode:              ", mode, "\n")
      cat("Compared:          ", paste(labels, collapse = " vs "), "\n")
      cat("Method used:       ", method_used, "\n")
      cat("Correlation (r):   ", estimate, "\n")
      cat("p-value:           ", signif(test$p.value, 3), "\n")
      cat("Observations (n):  ", length(x), "\n")

      if (!is.null(mi_value)) {
        cat("Mutual Information:", round(mi_value, 3), "\n")
      }

      if (!is.null(mi_norm)) {
        cat("Normalized MI:     ", round(mi_norm, 3), "\n")
      }

      if (!is.null(non_linear_signal)) {

        label <- if (non_linear_signal >= 0) {
          "Non-linear signal: "
        } else {
          "Linear dominance: "
        }

        cat(label, round(non_linear_signal, 3), "\n")
      }

      if (!is.null(relation_type)) {
        cat("Relationship type: ", relation_type, "\n")
      }
      cat("Interpretation:    ", interpretation, "\n")
    })

    if (verbose && method == "auto") {
      .print_block("Diagnostics", function() {
        cat("Outlier fraction:  ",
            round(diagnostics$outlier_fraction, 4), "\n")
        cat("Ties detected:     ",
            diagnostics$ties, "\n")
        cat("Normal X:          ",
            diagnostics$normal_x, "\n")
        cat("Normal Y:          ",
            diagnostics$normal_y, "\n")
      })
    }
  }

  return(invisible(proj))
}
