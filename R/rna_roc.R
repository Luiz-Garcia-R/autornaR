#' ROC analysis for gene signatures
#'
#' @description
#' Computes Receiver Operating Characteristic (ROC) curves and Area Under the Curve (AUC)
#' for gene-level or multi-gene signatures using normalized RNA-seq data stored in
#' the active \code{rna_project}.
#'
#' This function provides a flexible framework to evaluate the discriminative
#' performance of individual genes or gene sets in binary classification problems.
#'
#' @param project \code{rna_project} object created by \code{rna.project()}.
#' @param genes Character vector of gene symbols or Ensembl IDs, or \code{"all"} to use all genes.
#' @param group_col Character. Column name in metadata defining binary groups.
#' @param score_method Character. One of \code{"auto"}, \code{"pca"}, or \code{"mean"}.
#' @param show_gene_names Logical. Wheter to display gene names in subtitle.
#' @param plot Logical. Whether to display the ROC curve.
#' @param save Logical. Whether to store results in the active \code{rna_project}.
#'
#' @return
#' An object of class \code{"roc_result"} containing:
#' \describe{
#'   \item{roc_object}{ROC object from \pkg{pROC}.}
#'   \item{auc}{Numeric AUC value.}
#'   \item{auc_ci}{Confidence interval for AUC.}
#'   \item{genes}{Genes used in the analysis.}
#'   \item{method}{Scoring strategy used (\code{"single_gene"}, \code{"logistic"}, or \code{"signature"}).}
#'   \item{score_method}{Method used for multi-gene scoring (\code{"pca"} or \code{"mean"}).}
#'   \item{best_threshold}{Optimal classification threshold (Youden index).}
#'   \item{best_sens}{Sensitivity at optimal threshold.}
#'   \item{best_spec}{Specificity at optimal threshold.}
#' }
#'
#' @details
#' The function adapts the prediction strategy based on the number of genes provided:
#'
#' \strong{Single gene}
#' Expression values are directly used as predictor scores.
#'
#' \strong{Small gene sets (<= 5 genes)}
#' A logistic regression model is fitted using gene expression as predictors.
#'
#' \strong{Larger gene sets (> 5 genes)}
#' A composite signature score is computed:
#' \itemize{
#'   \item \code{pca}: first principal component (default in \code{"auto"} mode)
#'   \item \code{mean}: average expression across genes
#' }
#'
#' The ROC curve is computed using \pkg{pROC}, and the optimal threshold is
#' determined using the Youden index.
#'
#' Gene identifiers can be provided as Ensembl IDs or gene symbols. When possible,
#' gene symbols are retrieved automatically using \pkg{biomaRt}.
#'
#' @section Interpretation guide:
#' - AUC ~ 0.5 indicates no discriminative ability
#' - AUC 0.7–0.8 indicates moderate performance
#' - AUC 0.8–0.9 indicates strong performance
#' - AUC > 0.9 may indicate very strong signal or potential overfitting
#'
#' For multi-gene signatures, interpretation depends on how the score is constructed.
#' PCA-based scores capture dominant variance structure, while mean-based scores
#' assume consistent directionality across genes.
#'
#' @section Important considerations:
#' \strong{Binary outcome required}
#' The grouping variable must contain exactly two levels.
#'
#' \strong{Unsupervised signatures}
#' PCA-based scoring does not use class labels and may not yield optimal classification.
#'
#' \strong{Large gene sets}
#' Using \code{"all"} genes may produce unstable or uninterpretable results,
#' especially without prior feature selection.
#'
#' \strong{Data leakage}
#' When evaluating predictive performance, avoid using genes selected
#' from the same dataset without proper validation.
#'
#' @section Pipeline context:
#' This function is typically used after \code{rna.normalize()} and optionally after
#' feature selection steps such as \code{rna.compare()}. It complements visualization
#' tools like \code{rna.heatmap()} by providing quantitative classification metrics.
#'
#' @section Side effects:
#' If \code{save = TRUE}, results are stored in the \code{rna_project} under the
#' \code{analyses$roc} slot.
#'
#' @examples
#' \dontrun{
#' # Single gene ROC
#' rna.roc(project = my_project,
#'         genes = "CD8A"
#' )
#'
#' # Multi-gene logistic model
#' rna.roc(my_project,
#'         genes = c("CD3D", "CD8A", "GZMB")
#' )
#'
#' # Signature using PCA score
#' rna.roc(my_project,
#'         genes = c("Gene1", "Gene2", "Gene3", "Gene4", "Gene5", "Gene6"),
#'         score_method = "pca"
#' )
#' }
#'
#' @importFrom pROC roc auc
#' @export

rna.roc <- function(project,
                    genes = NULL,
                    group_col = "Group",
                    score_method = c("auto", "pca", "mean"),
                    show_gene_names = TRUE,
                    plot = TRUE,
                    save = TRUE) {

  score_method <- match.arg(score_method)

  # ---------------------------
  # 0) Basic checks
  # ---------------------------
  pkgs <- c("pROC", "ggplot2", "dplyr")

  .check_dependencies(pkgs)

  # ---------------------------
  # 1) Get active project
  # ---------------------------
  proj <- project

  expr_mat <- as.matrix(.get_expr(proj))
  metadata <- .get_meta(proj)
  organism <- .get_organism(proj)

  # ---------------------------
  # 2) Validate input
  # ---------------------------
  if (is.null(metadata)) {
    stop("No normalized data found. Run rna.normalize() first.")
  }

  if (!is.null(organism)) {
    message("[rna.enrich] Using organism from project: ", organism)
  }

  # ---------------------------
  # 3) Gene annotation
  # ---------------------------
  gene_annotation <- .get_gene_annotation(proj)
  gene_map <- .align_gene_annotation(gene_annotation, expr_mat)

  gene_ids <- gene_map$gene_id
  gene_symbols <- gene_map$symbol

  # ---------------------------
  # 4) Select genes (FIXED)
  # ---------------------------
  if (is.null(genes)) {
    stop("Please provide at least one gene or 'all'.")
  }

  rownames_expr <- rownames(expr_mat)
  rownames_clean <- sub("\\..*$", "", rownames_expr)

  gene_map$gene_id <- sub("\\..*$", "", gene_map$gene_id)

  if (length(genes) == 1 && tolower(genes) == "all") {

    genes_use <- rownames_expr
    gene_labels <- genes_use

    message("[rna.roc] Using all genes.")

  } else {

    # match against annotation
    match_idx <- match(genes, gene_map$gene_id)

    # fallback para SYMBOL
    if (any(is.na(match_idx))) {
      symbol_match <- match(genes, gene_map$symbol)
      match_idx[is.na(match_idx)] <- symbol_match[is.na(match_idx)]
    }

    if (all(is.na(match_idx))) {
      stop(
        "No matching genes found in expression matrix.\n",
        "Examples: ", paste(head(genes, 5), collapse = ", ")
      )
    }

    # Valid index
    valid_idx <- match_idx[!is.na(match_idx)]

    # Using matrix position
    genes_use <- rownames_expr[valid_idx]
    gene_labels <- gene_map$symbol[valid_idx]
    gene_labels[is.na(gene_labels) | gene_labels == ""] <- genes_use
  }

  safe_gene_labels <- make.names(gene_labels)

  if (length(genes_use) > 1000) {
    warning("[rna.roc] Large gene set detected. Results may be unstable.")
  }

  # ---------------------------
  # 5) Prepare data
  # ---------------------------
  expr_mat <- expr_mat[genes_use, , drop = FALSE]

  df <- as.data.frame(t(expr_mat))
  colnames(df) <- safe_gene_labels

  if (!group_col %in% colnames(metadata)) {
    stop(
      "Column '", group_col, "' not found in metadata.\n",
      "Available columns: ", paste(colnames(metadata), collapse = ", ")
    )
  }

  df[[group_col]] <- droplevels(as.factor(metadata[[group_col]]))

  if (length(levels(df[[group_col]])) != 2) {
    stop("ROC requires exactly 2 groups.")
  }

  df[[group_col]] <- stats::relevel(df[[group_col]], ref = levels(df[[group_col]])[1])

  # ---------------------------
  # 6) Method
  # ---------------------------
  n_genes <- length(genes_use)

  method <- if (n_genes == 1) {
    "single_gene"
  } else if (n_genes <= 5) {
    "logistic"
  } else {
    "signature"
  }

  # ---------------------------
  # 7) Prediction
  # ---------------------------
  if (method == "single_gene") {

    pred <- df[[safe_gene_labels[1]]]

  } else if (method == "logistic") {

    formula_str <- paste(group_col, "~", paste(safe_gene_labels, collapse = " + "))

    model <- stats::glm(
      as.formula(formula_str),
      data = df,
      family = stats::binomial
    )

    pred <- stats::predict(model, type = "response")

  } else if (method == "signature") {

    if (score_method == "auto") {
      score_method_use <- "pca"
    } else {
      score_method_use <- score_method
    }

    if (score_method_use == "mean") {

      pred <- rowMeans(df[, safe_gene_labels, drop = FALSE])

    } else if (score_method_use == "pca") {

      pca <- stats::prcomp(df[, safe_gene_labels, drop = FALSE], scale. = TRUE)
      pred <- pca$x[, 1]

      df[, safe_gene_labels, drop = FALSE]
      pred <- pca$x[, 1]

    }

    message("[rna.roc] Using signature score: ", score_method_use)
  }

  # ---------------------------
  # 8) ROC
  # ---------------------------

  valid_idx <- which(!is.na(pred) & !is.na(df[[group_col]]))

  pred_clean <- pred[valid_idx]
  group_clean <- df[[group_col]][valid_idx]

  if (length(pred_clean) == 0) {
    stop("No valid data available after removing NA values.")
  }

  if (length(unique(pred_clean)) < 2) {
    stop("Predictor has no variation. ROC cannot be computed.")
  }

  if (length(unique(group_clean)) != 2) {
    stop("ROC requires exactly 2 groups after filtering.")
  }

  roc_obj <- pROC::roc(
    group_clean,
    pred_clean,
    quiet = TRUE
  )

  auc_val <- pROC::auc(roc_obj)
  ci_auc <- pROC::ci.auc(roc_obj)

  coords <- suppressWarnings(
    pROC::coords(
      roc_obj,
      x = "best",
      best.method = "youden",
      ret = c("threshold", "sensitivity", "specificity"),
      transpose = FALSE
    )
  )

  boot <- .bootstrap_auc(
    df = df,
    group_col = group_col,
    safe_gene_labels = safe_gene_labels,
    method = method,
    score_method = ifelse(method == "signature", score_method_use, NA),
    n_boot = 1000
  )

  # ---------------------------
  # 9) Plot
  # ---------------------------
  if (plot) {

    plot(
      1 - roc_obj$specificities,
      roc_obj$sensitivities,
      type = "l",
      col = "#1f77b4",
      lwd = 3,
      main = paste("ROC -", toupper(method)),
      xlab = "False Positive Rate",
      ylab = "True Positive Rate",
      xlim = c(0,1),
      ylim = c(0,1),
      bty = "n"
    )

    abline(0, 1, lty = 2, col = "gray")

    text(
      0.8, 0.1,
      paste("AUC =", round(auc_val, 3)),
      col = "#1f77b4",
      font = 2
    )

    if (show_gene_names && length(genes_use) <= 5) {
      mtext(
        paste("Genes:", paste(gene_labels, collapse = " + ")),
        side = 3
      )
    }

    grid()
  }

  # ---------------------------
  # 10) Output
  # ---------------------------
  params = list(
    timestamp = Sys.time(),
    n_genes = length(gene_labels),
    n_samples = length(roc_obj$response),
    genes = gene_labels,
    method = method,
    score_method = ifelse(method == "signature", score_method_use, NA),
    group_col = group_col
  )

  obj <- list(
    params = params,
    statistics = list(
      auc = as.numeric(auc_val),
      auc_ci = as.numeric(ci_auc),
      auc_boot_ci = boot$ci,
      best_threshold = coords["threshold"],
      best_sens = coords["sensitivity"],
      best_spec = coords["specificity"]
    ),
    roc = roc_obj
  )

  class(obj) <- "roc_result"

  # ---------------------------
  # 11) Attach to project
  # ---------------------------
  if (save) {

  proj <- .attach_to_project(
    proj,
    obj,
    slot = "analyses",
    subtype = "roc",
    prefix = "roc",
    log = list(
      score_method = ifelse(method == "signature", score_method_use, NA),
      genes = gene_labels,
      group_col = group_col,
      auc = as.numeric(auc_val),
      direction = roc_obj$direction
    )
  )
}

  if (plot) print(obj)
  return(invisible(proj))

}

# ---------------------------
# 12) Print S3
# ---------------------------
#' @export
print.roc_result <- function(x, ...) {

  .print_header("ROC Analysis Result")

  # Extractions
  genes <- x$params$genes
  auc <- x$statistics$auc
  auc_ci <- x$statistics$auc_ci
  auc_boot_ci <- x$statistics$auc_boot_ci
  group_col <- x$params$group_col
  method <- x$params$method
  score_method <- x$params$score_method

  .print_block("Overview", function() {

    cat("Group column: ", group_col, "\n", sep = "")
    cat("Number of genes: ", length(genes), "\n", sep = "")

    n_show <- 5

    if (length(genes) <= n_show) {
      cat("Genes used: ", paste(genes, collapse = ", "), "\n", sep = "")
    } else {
      cat(
        "Genes used (first ", n_show, " of ", length(genes), "): ",
        paste(head(genes, n_show), collapse = ", "),
        " ...\n", sep = ""
      )
    }

    if (!is.null(score_method) && !is.na(score_method)) {
      cat("Score method: ", score_method, "\n", sep = "")
    }

    cat("Method: ", method, "\n", sep = "")

    fmt_ci <- function(ci) {
      paste0(round(ci[1], 3), " - ", round(ci[length(ci)], 3))
    }

    cat(
      "AUC (Delong): ", round(auc, 3),
      " [", fmt_ci(auc_ci), "]\n",
      sep = ""
    )

    if (!is.null(auc_boot_ci)) {
      cat("Bootstrap (stability): [", fmt_ci(auc_boot_ci), "]\n",
        sep = ""
      )
    }
    if (abs(diff(auc_ci[c(1,3)]) - diff(auc_boot_ci)) > 0.05) {
        message("[rna.roc] Notice: bootstrap and DeLong CI differ - model may be unstable.")
      }

  })

  invisible(x)
}

#' @export
summary.roc_result <- function(object, ...) {
  return(object$roc_object)
}
