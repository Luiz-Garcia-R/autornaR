#' Expression heatmap visualization
#'
#' This function generates a clustered heatmap from normalized RNA-seq data,
#' optionally using the top differentially expressed genes (DEGs), user-specified
#' gene sets, or genes with the highest expression variance.
#'
#' @param norm_data A `normalized_data` object containing `expr_matrix` and `metadata`.
#' @param comp_data An optional result object from [rna_compare()] containing
#'   differential expression results (`$res`). If provided, the top DEGs will
#'   be selected automatically.
#' @param group_col Character. Column name in metadata representing the sample groups
#'   (default: `"Group"`).
#' @param top_n Integer. Number of top genes to display when using `comp_data` or
#'   `use_variance = TRUE` (default: `100`).
#' @param genes Optional character vector of gene identifiers to include manually.
#' @param use_variance Logical. If `TRUE`, selects the most variable genes across samples
#'   (default: `FALSE`).
#' @param zscore Logical. Whether to z-score each gene (row) prior to plotting
#'   (default: `TRUE`).
#' @param cluster_rows Logical. Whether to cluster genes (rows) in the heatmap
#'   (default: `TRUE`).
#' @param cluster_cols Logical. Whether to cluster samples (columns) in the heatmap
#'   (default: `TRUE`).
#' @param show_rownames Logical. Whether to show gene names (default: `FALSE`).
#' @param show_colnames Logical. Whether to show sample names (default: `TRUE`).
#' @param palette Color palette used for the heatmap (default:
#'   `grDevices::colorRampPalette(c("navy", "white", "firebrick3"))(100)`).
#' @param results Logical. Whether to return a list containing the heatmap data
#'   (default: `FALSE`).
#' @param assign_result Logical. Whether to assign the output object to the environment
#'   (default: `FALSE`).
#' @param assign_name Character. Name of the object to assign (default: `"heatmap_data"`).
#' @param envir Environment in which to assign the output (default: `parent.frame()`).
#'
#' @return If `results = TRUE`, returns a list containing:
#' \item{top_genes}{Character vector of genes used in the heatmap.}
#' \item{expr_z_score}{Z-scored expression matrix (if applicable).}
#' \item{annotation_col}{Data frame with sample group annotations.}
#'
#' Otherwise, the function returns invisibly after plotting the heatmap.
#'
#' @details
#' The heatmap is created using the \pkg{pheatmap} package. If `comp_data` is provided,
#' the top `top_n` genes are selected based on the adjusted p-value (`padj`) or
#' raw p-value (`pvalue`) if `padj` is not available.
#' Alternatively, users can provide a list of specific genes or select the most
#' variable genes via `use_variance = TRUE`.
#'
#' @examples
#' \dontrun{
#' rna.heatmap(norm_data, comp_data = comp_result, top_n = 50)
#' rna.heatmap(norm_data, use_variance = TRUE, top_n = 100)
#' rna.heatmap(norm_data, genes = c("GeneA", "GeneB", "GeneC"))
#' }
#'
#' @export
rna.heatmap <- function(
    norm_data,
    comp_data = NULL,
    group_col = "Group",
    top_n = 100,
    genes = NULL,
    use_variance = FALSE,
    zscore = TRUE,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = FALSE,
    show_colnames = TRUE,
    palette = grDevices::colorRampPalette(c("navy", "white", "firebrick3"))(100),
    results = FALSE,
    assign_result = FALSE,
    assign_name = "heatmap_data",
    envir = parent.frame()
) {
  # --- Package check ---
  if (!requireNamespace("pheatmap", quietly = TRUE))
    stop("Please install the 'pheatmap' package.")

  # --- Check normalized data structure ---
  if (!inherits(norm_data, "normalized_data"))
    stop("`norm_data` must be of class 'normalized_data'.")

  expr_mat <- as.matrix(norm_data$expr_matrix)

  # --- Gene selection ---
  if (!is.null(comp_data) && "res" %in% names(comp_data)) {
    score_col <- NULL
    if ("padj" %in% colnames(comp_data$res) &&
        sum(!is.na(comp_data$res$padj)) > 0) {
      score_col <- "padj"
    } else if ("pvalue" %in% colnames(comp_data$res)) {
      score_col <- "pvalue"
      message("[rna_heatmap] 'padj' not available; using 'pvalue' instead.")
    } else {
      stop("`comp_data$res` must contain either 'padj' or 'pvalue'.")
    }
    genes_sel <- rownames(
      head(comp_data$res[order(comp_data$res[[score_col]]), , drop = FALSE], top_n)
    )
  } else if (!is.null(genes)) {
    genes_sel <- genes
  } else if (use_variance) {
    var_genes <- apply(expr_mat, 1, var, na.rm = TRUE)
    genes_sel <- names(sort(var_genes, decreasing = TRUE))[1:min(top_n, length(var_genes))]
  } else {
    stop("Provide `comp_data`, a list of `genes`, or set `use_variance = TRUE`.")
  }

  genes_sel <- genes_sel[genes_sel %in% rownames(expr_mat)]
  mat_top <- expr_mat[genes_sel, , drop = FALSE]

  # --- Apply z-score normalization ---
  if (zscore)
    mat_top <- t(scale(t(mat_top)))

  # --- Prepare annotation for columns ---
  metadata <- norm_data$metadata
  if (!group_col %in% colnames(metadata))
    stop(sprintf("Column '%s' not found in metadata.", group_col))

  rownames(metadata) <- metadata$Sample
  ann_col <- data.frame(Group = as.factor(metadata[colnames(mat_top), group_col]))
  rownames(ann_col) <- colnames(mat_top)

  # --- Plot heatmap ---
  pheatmap::pheatmap(
    mat_top,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    show_rownames = show_rownames,
    show_colnames = show_colnames,
    annotation_col = ann_col,
    color = palette,
    border_color = NA,
    main = if (!is.null(comp_data) && "res" %in% names(comp_data)) {
      sprintf("Top %d DEGs", length(genes_sel))
    } else {
      "Gene Expression Heatmap"
    }
  )

  # --- Prepare results ---
  result <- list(
    top_genes = genes_sel,
    expr_z_score = mat_top,
    annotation_col = ann_col
  )

  if (assign_result)
    assign(assign_name, result, envir = envir)

  if (results)
    return(result)
  else
    invisible(NULL)
}
