#' Perform dimensionality reduction (PCA, UMAP, t-SNE) on normalized RNA-seq data
#'
#' This function applies multiple dimensionality reduction techniques — PCA,
#' UMAP, and optionally t-SNE — to visualize sample relationships in
#' normalized RNA-seq expression data. It automatically adjusts parameters for
#' small datasets and returns both results and plots.
#'
#' @param norm_data A `normalized_data` object containing an expression matrix
#'   (`expr_matrix`) and sample metadata (`metadata`).
#' @param group_col Character. Name of the column in `metadata` defining sample
#'   groups for coloring in plots (default: `"Group"`).
#' @param n_neighbors Integer. Number of neighbors for UMAP. If `NULL`,
#'   it is set automatically based on the number of samples (default: `NULL`).
#' @param tsne_perplexity Numeric. Perplexity parameter for t-SNE (default: `30`).
#' @param ncomp_pca Integer. Number of PCA components to compute (default: `2`).
#' @param verbose Logical. If `TRUE`, prints progress and displays plots
#'   (default: `TRUE`).
#'
#' @return A list containing:
#' \item{expr_matrix}{The expression matrix used for dimensionality reduction.}
#' \item{metadata}{The metadata used for annotation.}
#' \item{res_dimred}{A list of PCA, UMAP, and (if applicable) t-SNE results.}
#' \item{plots}{A list of `ggplot` objects for PCA, UMAP, and t-SNE visualizations.}
#'
#' @details
#' Dimensionality reduction techniques such as PCA, UMAP, and t-SNE are commonly
#' used to visualize relationships among RNA-seq samples. This function performs:
#' \itemize{
#'   \item \strong{PCA:} Principal Component Analysis using `stats::prcomp`.
#'   \item \strong{UMAP:} Uniform Manifold Approximation and Projection via
#'   the `umap` package.
#'   \item \strong{t-SNE:} t-distributed Stochastic Neighbor Embedding via
#'   `Rtsne` (skipped for <6 samples).
#' }
#'
#' @examples
#' \dontrun{
#' dim_res <- rna.dimred(normalized_data)
#' dim_res$plots$PCA
#' dim_res$plots$UMAP
#' }
#'
#' @export
rna.dimred <- function(
    norm_data,
    group_col = "Group",
    n_neighbors = NULL,
    tsne_perplexity = 30,
    ncomp_pca = 2,
    verbose = TRUE
) {
  # --- Required packages ---
  required_pkgs <- c("ggplot2", "umap", "Rtsne", "dplyr")
  for (p in required_pkgs) {
    if (!requireNamespace(p, quietly = TRUE))
      stop(sprintf("Package '%s' is required but not installed.", p))
  }

  # --- Validate input ---
  if (!inherits(norm_data, "normalized_data"))
    stop("'norm_data' must be of class 'normalized_data'.")

  if (!group_col %in% colnames(norm_data$metadata))
    stop(sprintf("Column '%s' not found in metadata.", group_col))

  expr_mat <- as.matrix(norm_data$expr_matrix)
  metadata <- norm_data$metadata
  samples <- metadata$Sample
  groups <- as.factor(metadata[[group_col]])

  # --- Ensure sample order consistency ---
  if (!all(colnames(expr_mat) == samples))
    expr_mat <- expr_mat[, samples, drop = FALSE]

  # --- Transpose matrix for reduction (samples as rows) ---
  expr_mat_t <- t(expr_mat)
  n_samples <- nrow(expr_mat_t)

  # --- PCA ---
  if (verbose) message("[rna_dimred] Running PCA...")
  pca_res <- stats::prcomp(expr_mat_t, center = TRUE, scale. = TRUE)
  var_explained <- pca_res$sdev^2 / sum(pca_res$sdev^2)

  df_pca <- data.frame(
    PC1 = pca_res$x[, 1],
    PC2 = pca_res$x[, 2],
    Group = groups,
    Sample = samples
  )

  p_pca <- ggplot2::ggplot(
    df_pca,
    ggplot2::aes(PC1, PC2, color = Group, label = Sample)
  ) +
    ggplot2::geom_point(size = 4, alpha = 0.8) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "PCA - Dimensionality Reduction",
      x = sprintf("PC1 (%.1f%%)", 100 * var_explained[1]),
      y = sprintf("PC2 (%.1f%%)", 100 * var_explained[2])
    )

  # --- UMAP ---
  if (is.null(n_neighbors) || n_neighbors >= n_samples)
    n_neighbors <- max(2, floor(n_samples / 2))

  if (verbose)
    message(sprintf("[rna_dimred] UMAP n_neighbors set to %d", n_neighbors))

  umap_res <- umap::umap(expr_mat_t, n_neighbors = n_neighbors, init = "random")

  df_umap <- data.frame(
    UMAP1 = umap_res$layout[, 1],
    UMAP2 = umap_res$layout[, 2],
    Group = groups,
    Sample = samples
  )

  p_umap <- ggplot2::ggplot(
    df_umap,
    ggplot2::aes(UMAP1, UMAP2, color = Group, label = Sample)
  ) +
    ggplot2::geom_point(size = 4, alpha = 0.8) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "UMAP - Dimensionality Reduction")

  # --- t-SNE (only if sufficient samples) ---
  run_tsne <- n_samples > 5
  if (run_tsne) {
    tsne_perplexity <- min(tsne_perplexity, floor((n_samples - 1) / 3))

    if (verbose)
      message(sprintf("[rna_dimred] Running t-SNE (perplexity = %d)...",
                      tsne_perplexity))

    tsne_res <- Rtsne::Rtsne(
      expr_mat_t,
      perplexity = tsne_perplexity,
      check_duplicates = FALSE,
      verbose = FALSE
    )

    df_tsne <- data.frame(
      tSNE1 = tsne_res$Y[, 1],
      tSNE2 = tsne_res$Y[, 2],
      Group = groups,
      Sample = samples
    )

    p_tsne <- ggplot2::ggplot(
      df_tsne,
      ggplot2::aes(tSNE1, tSNE2, color = Group, label = Sample)
    ) +
      ggplot2::geom_point(size = 4, alpha = 0.8) +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = sprintf("t-SNE - Perplexity %d", tsne_perplexity))
  } else {
    if (verbose)
      message("[rna_dimred] t-SNE skipped: insufficient samples (<6) for reliable embedding.")
    tsne_res <- NULL
    p_tsne <- NULL
  }

  # --- Print plots if verbose ---
  if (verbose) {
    print(p_pca)
    print(p_umap)
    if (!is.null(p_tsne)) print(p_tsne)
  }

  # --- Return structured output ---
  invisible(list(
    expr_matrix = expr_mat,
    metadata = metadata,
    res_dimred = list(
      pca_res = pca_res,
      umap_res = umap_res,
      tsne_res = tsne_res
    ),
    plots = list(
      PCA = p_pca,
      UMAP = p_umap,
      tSNE = p_tsne
    )
  ))
}
