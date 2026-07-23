#' PCA Biplot with Feature Selection
#'
#' Generates a PCA biplot using previously computed dimensionality reduction
#' results, highlighting genes based on their contribution to principal components
#' and/or association with a grouping variable.
#'
#' The function supports three modes for gene selection:
#' \itemize{
#'   \item \strong{loading}: selects genes with highest contribution to selected PCs
#'   \item \strong{anova}: selects genes most associated with the grouping variable
#'   \item \strong{hybrid}: combines loading magnitude and ANOVA significance
#' }
#'
#' Gene selection can be automatic (top-ranked) or user-defined.
#'
#' @param project \code{rna_project} object created by \code{rna.project()}.
#' @param use_dimred Optional \code{dimred_object} to use in \code{rna.biplot}.
#' @param group_col Character. Column name in metadata used for grouping samples.
#' @param style Character. Method used for gene selection.
#'   One of \code{"loading"}, \code{"anova"}, or \code{"hybrid"}.
#' @param colors Character. Optional point colors.
#' @param n_genes Integer. Number of top genes to display when
#'   \code{biplot_genes} is not provided.
#' @param genes Character vector. Optional vector of gene names to include.
#'   Overrides automatic selection.
#' @param point_size Numeric. Size of sample points in the plot.
#' @param save Logical. Whether to store results in the active \code{rna_project}.
#' @param verbose Logical. Whether to print the plot.
#'
#' @return A list of class \code{biplot_result} containing:
#' \itemize{
#'   \item \code{biplot}: list with selected genes, coordinates, scores, and PCs used.
#'   \item \code{plot}: ggplot object of the biplot.
#' }
#'
#' @details
#' The biplot is constructed from a previously computed PCA stored in
#' \code{rna_project}. Gene vectors are scaled to match the PCA coordinate space.
#'
#' In \code{hybrid} mode, both loading magnitude and ANOVA-derived significance
#' are normalized to the range \code{[0,1]} before combination.
#'
#' @examples
#' \dontrun{
#' # Standard PCA biplot
#' rna.biplot(project = my_project)
#'
#' # PCA biplot with a specific gene set
#' rna.biplot(project = my_project,
#'            n_genes = 2,
#'            genes = c("IL6", "HSP70"))
#' }
#'
#' @importFrom utils tail
#'
#' @export

rna.biplot <- function(project,
                       use_dimred = NULL,
                       group_col = "Group",
                       style = c("loading", "anova", "hybrid"),
                       colors = NULL,
                       n_genes = 10,
                       genes = NULL,
                       point_size = 2,
                       save = TRUE,
                       verbose = TRUE
) {

  # ===========================================================================
  # 0) Basic checks
  # ===========================================================================
  .check_dependencies(c("dplyr", "ggplot2", "ggrepel"))

  # ===========================================================================
  # 1) Get active project
  # ===========================================================================
  proj <- project

  expr_mat <- as.matrix(.get_expr(proj))
  metadata <- .get_meta(proj)
  dimred <- .get_dimred(proj, id = use_dimred)

  # ===========================================================================
  # 2) Validate input
  # ===========================================================================
  pca_res <- dimred$PCA$model
  loadings <- dimred$PCA$loadings
  df_pca <- dimred$PCA$coordinates

  common <- intersect(colnames(expr_mat), metadata$Sample)

  expr_mat <- expr_mat[, common, drop = FALSE]
  metadata <- metadata[match(common, metadata$Sample), , drop = FALSE]

  groups <- as.factor(metadata[[group_col]])

  # ===========================================================================
  # 3) Gene mapping
  # ===========================================================================
  gene_map <- .get_gene_annotation(proj)

  gene_ids <- gene_map$gene_id
  gene_symbols <- gene_map$symbol

  # fallback
  gene_labels <- ifelse(
    is.na(gene_symbols) | gene_symbols == "",
    gene_ids,
    gene_symbols
  )

  names(gene_labels) <- gene_ids

  # ===========================================================================
  # 4) Computate loadings
  # ===========================================================================
  style <- match.arg(style)
  p_biplot <- NULL

  pc_x = 1
  pc_y = 2

  # --- gene selection ---
  if (!is.null(genes)) {
    valid_idx <- which(
      gene_ids %in% genes |
        gene_symbols %in% genes
    )

    if (length(valid_idx) == 0) {
      warning("None of the provided genes were found.")
      selected_genes <- character(0)
    } else {
      selected_genes <- gene_ids[valid_idx]
    }

  } else {

    all_genes <- gene_ids[gene_ids %in% rownames(loadings)]

    # --- 1) Loadings ---
    pc_load <- loadings[, pc_x]
    loading_score <- sqrt(loadings[, pc_x]^2 + loadings[, pc_y]^2)

    # --- 2) ANOVA ---
    if (style %in% c("anova", "hybrid")) {

      pvals <- apply(expr_mat[all_genes, , drop = FALSE], 1, function(g) {
        tryCatch({
          summary(stats::aov(g ~ groups))[[1]][["Pr(>F)"]][1]
        }, error = function(e) NA)
      })

      # Avoiding p=0
      pvals[pvals == 0] <- 1e-300
      anova_score <- -log10(pvals)

    } else {
      anova_score <- rep(1, length(all_genes))
    }

    # --- 3) final Score ---
    if (style == "loading") {
      score <- loading_score

    } else if (style == "anova") {
      score <- anova_score

    } else if (style == "hybrid") {
      scale01 <- function(x) {
        rng <- range(x, na.rm = TRUE)
        if (diff(rng) == 0) return(rep(0, length(x)))
        (x - rng[1]) / diff(rng)
      }

      loading_norm <- scale01(loading_score)
      anova_norm <- scale01(anova_score)

      score <- loading_norm * anova_norm
    }

    # --- 4) Ordering ---
    ord <- order(score, decreasing = TRUE)
    n_use <- min(n_genes, length(ord))
    selected_genes <- all_genes[ord][1:n_use]
  }
  # --- validation ---
  if (length(selected_genes) == 0) {
    warning("No valid genes found for biplot.")
  } else {

    load_df <- data.frame(
      gene_id = selected_genes,
      gene_label = gene_labels[selected_genes],
      PC1 = loadings[selected_genes, pc_x] * pca_res$sdev[pc_x],
      PC2 = loadings[selected_genes, pc_y] * pca_res$sdev[pc_y]
    )

    arrow_scale <- 0.8 * max(abs(df_pca$PC1), abs(df_pca$PC2)) /
      max(abs(load_df$PC1), abs(load_df$PC2))

    load_df$PC1 <- load_df$PC1 * arrow_scale
    load_df$PC2 <- load_df$PC2 * arrow_scale

    df_pca$Group <- groups

  # ===========================================================================
  # 5) Compute Variation
  # ===========================================================================
  var_exp <- (pca_res$sdev^2) / sum(pca_res$sdev^2)

  pc1_var <- round(var_exp[pc_x] * 100, 1)
  pc2_var <- round(var_exp[pc_y] * 100, 1)

  # ===========================================================================
  # 6) Plot
  # ===========================================================================
    p_biplot <- ggplot2::ggplot(df_pca,
                                ggplot2::aes(PC1, PC2, color = Group)) +
      ggplot2::geom_point(size = point_size, alpha = 0.8) +

      # Arrows
      ggplot2::geom_segment(
        data = load_df,
        ggplot2::aes(x = 0, y = 0, xend = PC1, yend = PC2),
        arrow = ggplot2::arrow(length = grid::unit(0.25, "cm")),
        color = "grey30",
        linewidth = 0.6,
        alpha = 0.8,
        inherit.aes = FALSE
      ) +

      # Invisible anchor point
      ggplot2::geom_point(
        data = load_df,
        ggplot2::aes(x = PC1, y = PC2),
        alpha = 0,
        inherit.aes = FALSE
      ) +

      # Smart labels
      ggrepel::geom_text_repel(
        data = load_df,
        ggplot2::aes(x = PC1, y = PC2, label = .data$gene_label),
        size = 3,
        segment.color = "grey50",
        max.overlaps = Inf,
        box.padding = 0.3,
        point.padding = 0.2,
        inherit.aes = FALSE
      ) +

      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = paste0("PCA Biplot (", style, ")"),
        x = paste0("PC", pc_x, " (", pc1_var, "%)"),
        y = paste0("PC", pc_y, " (", pc2_var, "%)")
      )

  if (!is.null(colors)) {
    p_biplot <- p_biplot +
      ggplot2::scale_color_manual(values = colors)
  }

}

  # ===========================================================================
  # 7) Print plots if verbose
  # ===========================================================================
  if (verbose) {
    print(p_biplot)
  }

  # ===========================================================================
  # 8) Output
  # ===========================================================================
  obj <- NULL

  obj <- list(
    params = list(
      method = style,
      n_genes = length(selected_genes),
      group_col = group_col,
      dimred_used = use_dimred,
      timestamp = Sys.time()
    ),

    genes = data.frame(
      gene_id = selected_genes,
      gene_label = gene_labels[selected_genes],
      score = if (exists("score")) score[selected_genes] else NA,
      loading_PC1 = loadings[selected_genes, pc_x],
      loading_PC2 = loadings[selected_genes, pc_y],
      coord_PC1 = load_df$PC1,
      coord_PC2 = load_df$PC2,
      row.names = NULL
    ),

    pca = list(
      pcs = c(pc_x, pc_y),
      var_explained = var_exp,
      sdev = pca_res$sdev
    ),

    samples = df_pca,

    plot = p_biplot
  )

  class(obj) <- "biplot_result"

  # ===========================================================================
  # 9) Attach to project
  # ===========================================================================
  if (save) {

    proj <- .attach_to_project(
      proj,
      obj,
      slot = "analyses",
      subtype = "biplot",
      prefix = "biplot",
      log = list(
        dimred_used = use_dimred,
        n_genes_selected = length(selected_genes),
        mode = style,
        group_col = group_col,
        pcs = c(pc_x, pc_y),
        used_custom_genes = !is.null(genes)
      )
    )
  }

  # ===========================================================================
  # 10) Return
  # ===========================================================================
    has_score <- exists("score")

    .print_header("PCA biplot results")

    .print_block("Biplot Summary", function() {
      cat("Mode:               ", style, "\n")
      cat("Genes selected:     ", length(selected_genes), "\n")
      cat("PCs used:           ", "PC", pc_x, " vs PC", pc_y, "\n")
    })

    if (exists("score")) {
      .print_block("Top contributing genes", function() {
        df <- data.frame(
          gene = gene_labels[selected_genes],
          gene_id = selected_genes,
          score = round(score[selected_genes], 3),
          row.names = NULL
        )
        print(df)
      })
    }

    return(invisible(proj))
}
