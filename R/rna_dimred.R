#' Dimensionality reduction for RNA-seq data (PCA, UMAP, t-SNE)
#'
#' @description
#' Performs dimensionality reduction on normalized RNA-seq expression data
#' stored in the active \code{rna_project}, generating low-dimensional
#' representations of samples using Principal Component Analysis (PCA),
#' Uniform Manifold Approximation and Projection (UMAP), and optionally
#' t-distributed Stochastic Neighbor Embedding (t-SNE).
#'
#' The function automatically adapts parameters based on dataset size and
#' structure, and returns both numerical results and publication-ready plots.
#'
#' @param project \code{rna_project} object created by \code{rna.project()}.
#' @param group_col Character. Metadata column used to define sample groups
#'   for coloring plots (default: \code{"Group"}).
#' @param run_tsne Logical. Whether to compute t-SNE embeddings (default: \code{FALSE}).
#' @param add_ellipse Logical. If \code{TRUE}, adds confidence ellipses to PCA plot
#'   when group sizes allow (default: \code{FALSE}).
#' @param ellipse_level Numeric. Confidence level for PCA ellipses
#'   (default: \code{0.95}).
#' @param colors Character. Optional point colors.
#' @param ncomp_pca Integer. Maximum number of PCA components to compute.
#'   If \code{NULL}, uses the maximum possible (default).
#' @param n_neighbors Integer. Number of neighbors for UMAP. If \code{NULL},
#'   automatically set based on sample size.
#' @param tsne_perplexity Numeric. Initial perplexity for t-SNE
#'   (automatically adjusted if too large; default: \code{30}).
#' @param use_pca_preprocessing Logical. Whether to use PCA-reduced space as input
#'   for UMAP and t-SNE (default: \code{TRUE}).
#' @param pca_method Character. Strategy for selecting number of PCs:
#'   \itemize{
#'     \item \code{"auto"}: selects PCs based on cumulative variance threshold
#'     \item \code{"variance"}: selects PCs reaching \code{pca_var_threshold}
#'     \item \code{"fixed"}: uses \code{ncomp_pca} directly
#'   }
#' @param pca_var_threshold Numeric. Target cumulative variance explained when
#'   using adaptive PCA selection (default: \code{0.8}).
#' @param pca_max_dims Integer. Maximum number of PCs allowed (default: \code{50}).
#' @param pca_min_dims Integer. Minimum number of PCs enforced (default: \code{10}).
#' @param point_size Numeric. Size of points in plots (default: \code{2}).
#' @param seed Optional integer used to set the random seed for reproducibility.
#' If `NULL`, the existing random number generator state is used.
#' @param save Logical. Whether to store results in the active
#'   \code{rna_project} (default: \code{TRUE}).
#' @param verbose Logical. If \code{TRUE}, prints progress messages and displays plots
#'   (default: \code{TRUE}).
#'
#' @details
#' The function expects normalized expression data available in
#' \code{rna_project$data$normalized_data}, including:
#' \itemize{
#'   \item \code{expr_matrix}: gene expression matrix (genes x samples)
#'   \item \code{metadata}: sample metadata
#' }
#'
#' Internally, the expression matrix is transposed so that samples are treated
#' as observations.
#'
#' \strong{Adaptive behavior:}
#' \itemize{
#'   \item UMAP neighbors are set to \code{floor(n_samples / 2)} when unspecified
#'   \item t-SNE perplexity is capped at \code{floor((n_samples - 1)/3)}
#'   \item t-SNE is skipped when \code{n_samples < 6}
#'   \item PCA dimensionality is automatically selected when
#'   \code{use_pca_preprocessing = TRUE}
#' }
#'
#' \strong{Interpretation notes:}
#' \itemize{
#'   \item PCA captures global variance structure and is deterministic
#'   \item UMAP emphasizes local structure and clustering
#'   \item t-SNE is stochastic and sensitive to parameter choice
#' }
#'
#' @return
#' An object of class \code{"dimred_result"} containing:
#' \describe{
#'   \item{res_dimred}{List with dimensionality reduction results:
#'     \itemize{
#'       \item \code{summary}: dataset overview (samples, genes, grouping, seed and rng state)
#'       \item \code{preprocessing}: PCA preprocessing parameters
#'       \item \code{PCA}: PCA model, variance explained, loadings, coordinates
#'       \item \code{UMAP}: UMAP model, coordinates, and parameters
#'       \item \code{tSNE}: t-SNE model and coordinates (if computed)
#'     }
#'   }
#'   \item{plots}{List of \code{ggplot2} objects:
#'     \code{PCA}, \code{UMAP}, and optionally \code{tSNE}}
#' }
#'
#' @section Stored output:
#' When \code{save = TRUE}, results are stored in:
#' \code{rna_project$analyses$dimred}, enabling reproducibility and reuse.
#'
#' @examples
#' \dontrun{
#' # Run default dimensionality reduction
#' rna.dimred(project = my_project)
#'
#' # Run with t-SNE and ellipses
#' rna.dimred(my_project,
#'            run_tsne = TRUE,
#'            add_ellipse = TRUE)
#' }
#'
#' @importFrom stats prcomp dist
#' @importFrom ggplot2 ggplot aes geom_point labs theme_minimal
#' @importFrom ggplot2 stat_ellipse scale_fill_discrete
#'
#' @export

rna.dimred <- function(project,
                       group_col = "Group",
                       run_tsne = FALSE,
                       add_ellipse = FALSE,
                       ellipse_level = 0.95,
                       colors = NULL,
                       ncomp_pca = NULL,
                       n_neighbors = NULL,
                       tsne_perplexity = 30,
                       use_pca_preprocessing = TRUE,
                       pca_method = "auto",
                       pca_var_threshold = 0.8,
                       pca_max_dims = 50,
                       pca_min_dims = 10,
                       point_size = 2,
                       seed = NULL,
                       save = TRUE,
                       verbose = TRUE
) {

  # ===========================================================================
  # 0) Basic checks
  # ===========================================================================
  .check_dependencies(c("dplyr", "ggplot2", "umap"))

  if (run_tsne) {
    .check_dependencies("Rtsne")
  }

  # --- Set seed ---
  old_seed <- .set_seed(seed)
  on.exit(.reset_seed(old_seed), add = TRUE)

  # ===========================================================================
  # 1) Get active project
  # ===========================================================================
  proj <- project

  expr_mat <- .get_expr(proj)
  metadata <- .get_meta(proj)

  samples <- metadata$Sample
  groups <- as.factor(metadata[[group_col]])

  # ===========================================================================
  # 2) Validate input
  # ===========================================================================
  if (is.null(expr_mat)) {
    stop("No expression matrix available.")
  }

  if (is.null(metadata)) {
    stop("No metadata available.")
  }

  if (!"Sample" %in% colnames(metadata)) {
    stop("metadata must contain a 'Sample' column.")
  }

  if (!group_col %in% colnames(metadata)) {
    stop(
      "Column '", group_col, "' not found in metadata.\n",
      "Available columns: ", paste(colnames(metadata), collapse = ", ")
    )
  }

  # ===========================================================================
  # 3) Ensure sample order consistency
  # ===========================================================================
  if (!all(colnames(expr_mat) == samples))
    expr_mat <- expr_mat[, samples, drop = FALSE]

  if (nrow(expr_mat) < 50) {
    warning("Low number of genes for dimensionality reduction.")
  }

  # --- Transpose matrix for reduction ---
  expr_mat_t <- t(expr_mat)
  n_samples <- nrow(expr_mat_t)
  dist_mat <- dist(expr_mat_t)

  # --- Build auxiliar object ---
  ellipse_layer <- NULL

  if (add_ellipse) {
    group_sizes <- table(groups)

    if (all(group_sizes >= 3)) {
      ellipse_layer <- ggplot2::stat_ellipse(
        ggplot2::aes(fill = Group),
        geom = "polygon",
        alpha = 0.15,
        level = ellipse_level,
        color = NA
      )
    } else if (verbose) {
      message("[rna_dimred] Ellipse skipped: at least one group has < 3 samples.")
    }
  }

  # ===========================================================================
  # 4) PCA
  # ===========================================================================
  if (verbose) message("[rna_dimred] Running PCA...")

  if (is.null(ncomp_pca)) {
    ncomp_pca <- min(n_samples - 1, ncol(expr_mat_t))
  }

  pca_res <- stats::prcomp(expr_mat_t,
                           center = TRUE,
                           scale. = TRUE)
  var_explained <- pca_res$sdev^2 / sum(pca_res$sdev^2)
  cum_var <- cumsum(var_explained)

  loadings <- pca_res$rotation

  scree_df <- data.frame(
    PC = seq_along(var_explained),
    variance = var_explained,
    cumulative = cum_var
  )

  top_loadings <- do.call(
    rbind,
    lapply(1:min(5, ncol(loadings)), function(i) {
      pc_load <- loadings[, i]
      ord <- order(abs(pc_load), decreasing = TRUE)

      top_n <- min(10, length(ord))

      data.frame(
        gene = rownames(loadings)[ord][1:top_n],
        loading = pc_load[ord][1:top_n],
        PC = paste0("PC", i)
      )
    })
  )

  # ===========================================================================
  # 5) PCA preprocessing selection
  # ===========================================================================
  if (use_pca_preprocessing) {

    if (pca_method == "auto") {
      n_pcs <- which(cum_var >= pca_var_threshold)[1]
      n_pcs <- min(n_pcs, pca_max_dims)
      n_pcs <- max(n_pcs, pca_min_dims)

    } else if (pca_method == "variance") {
      n_pcs <- which(cum_var >= pca_var_threshold)[1]

    } else if (pca_method == "fixed") {
      n_pcs <- ncomp_pca
    }

    if (verbose) {
      message(sprintf("[rna_dimred] PCA preprocessing: using %d PCs (%s mode)",
                      n_pcs, pca_method))
    }

    input_mat <- pca_res$x[, 1:n_pcs, drop = FALSE]

  } else {
    input_mat <- expr_mat_t
    n_pcs <- NA
  }

  df_pca <- data.frame(
    PC1 = pca_res$x[, 1],
    PC2 = pca_res$x[, 2],
    Group = groups,
    Sample = samples
  )

  p_pca <- ggplot2::ggplot(
    df_pca,
    ggplot2::aes(PC1, PC2, color = Group)
  ) +
    ellipse_layer +
    ggplot2::geom_point(size = point_size, alpha = 0.8) +
    ggplot2::scale_fill_discrete() +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(
      title = "PCA",
      x = sprintf("PC1 (%.1f%%)", 100 * var_explained[1]),
      y = sprintf("PC2 (%.1f%%)", 100 * var_explained[2])
    )

  if (!is.null(colors)) {
    p_pca <- p_pca +
      ggplot2::scale_color_manual(values = colors)
  }

  # ===========================================================================
  # 6) UMAP
  # ===========================================================================
  if (is.null(n_neighbors) || n_neighbors >= n_samples)
    n_neighbors <- max(2, floor(n_samples / 2))

  if (verbose)
    message(sprintf("[rna_dimred] UMAP n_neighbors set to %d", n_neighbors))

  umap_res <- umap::umap(input_mat, n_neighbors = n_neighbors, init = "random")

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
    ggplot2::geom_point(size = point_size, alpha = 0.8) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(title = "UMAP")

  if (!is.null(colors)) {
    p_umap <- p_umap +
      ggplot2::scale_color_manual(values = colors)
  }

  # ===========================================================================
  # 7) t-SNE
  # ===========================================================================
  tsne_res <- NULL
  df_tsne <- NULL
  p_tsne <- NULL

  if (run_tsne && n_samples <= 5) {
    if (verbose) {
      message("[rna_dimred] t-SNE skipped: insufficient samples (<6).")
    }
    run_tsne <- FALSE
  }

  if (run_tsne && n_samples > 5) {
    tsne_perplexity <- min(tsne_perplexity, floor((n_samples - 1) / 3))

    if (verbose) {
      message(sprintf("[rna_dimred] Running t-SNE (perplexity = %d)...",
                      tsne_perplexity))
    }

    tsne_res <- Rtsne::Rtsne(
      input_mat,
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
      ggplot2::geom_point(size = point_size, alpha = 0.8) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::labs(title = sprintf("t-SNE", tsne_perplexity))

    if (!is.null(colors)) {
      p_tsne <- p_tsne +
        ggplot2::scale_color_manual(values = colors)
    }

  }

  # ===========================================================================
  # 8) Print plots
  # ===========================================================================
  if (verbose) {
    print(p_pca)
    print(p_umap)
    if (!is.null(p_tsne)) print(p_tsne)
  }

  # ===========================================================================
  # 9) Output
  # ===========================================================================
  rng_state <- if (!is.null(seed)) .Random.seed else NULL

  obj <- list(
    summary = list(
      timestamp = Sys.time(),
      n_samples = n_samples,
      n_genes = nrow(expr_mat),
      group_col = group_col,
      seed = seed,
      rng_state = rng_state
    ),
    preprocessing = list(
      pca_preprocessing = use_pca_preprocessing,
      pca_method = pca_method,
      n_pcs_used = n_pcs
    ),
    PCA = list(
      model = pca_res,
      variance_explained = var_explained,
      scree = scree_df,
      loadings = loadings,
      top_loadings = top_loadings,
      coordinates = df_pca
    ),
    UMAP = list(
      model = umap_res,
      coordinates = df_umap,
      n_neighbors = n_neighbors
    ),
    tSNE = if (run_tsne) list(
      model = tsne_res,
      coordinates = df_tsne,
      perplexity = tsne_perplexity
    ) else NULL,
    plots = list(
      PCA = p_pca,
      UMAP = p_umap,
      tSNE = p_tsne
    )
  )

  class(obj) <- "dimred_result"

  # ===========================================================================
  # 10) Attach to project
  # ===========================================================================
  if (save) {

    proj <- .attach_to_project(
      proj,
      obj,
      slot = "analyses",
      subtype = "dimred",
      prefix = "dimred",
      log = list(
        n_samples = n_samples,
        method = paste0("PCA + UMAP", ifelse(run_tsne, " + tSNE", "")),
        group_col = group_col,
        pcs = n_pcs,
        seed = seed
      )
    )
  }

  # ===========================================================================
  # 11) Return
  # ===========================================================================
  .print_header("RNA Dimensionality Reduction")

  .print_block("Summary", function() {
    cat("Samples:            ", n_samples, "\n")
    cat("Genes used:         ", nrow(expr_mat), "\n")
    cat("Group column:       ", group_col, "\n")
    cat("PCA components:     ", ncomp_pca, "\n")
    cat("PCA preprocessing:  ", ifelse(use_pca_preprocessing, "Yes", "No"), "\n")
    if (use_pca_preprocessing) {
      cat("PCs used (UMAP/t-SNE): ", n_pcs, "(", pca_method, ")\n")
    }
    cat("Cumulative variance (PC1+PC2): ", round(100 * cum_var[2], 2), "%\n")
    cat("UMAP neighbors:     ", n_neighbors, "\n")
    cat("t-SNE performed:    ", ifelse(run_tsne, "Yes", "No"), "\n")
    cat("Seed: ", seed, "\n")
  })

  .print_block("PCA Variance Explained", function() {
    n_show <- min(5, length(var_explained))
    print(round(100 * var_explained[1:n_show], 2))
  })

  return(invisible(proj))
}
