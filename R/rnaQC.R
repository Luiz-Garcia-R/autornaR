#' RNA-seq Quality Control (QC)
#'
#' @description
#' Perform basic quality control on normalized RNA-seq data.
#' Generates library size barplots, expression boxplots, density plots, 
#' and sample correlation heatmaps. Highlights potential outlier samples.
#'
#' @param normalized_data An object returned by [rna.normalize()], containing `expr_matrix` and `metadata`.
#' @param group_col Character; column in `metadata` indicating group assignment (default `"Group"`).
#' @param cluster_rows Logical; whether to cluster rows in the correlation heatmap (default `FALSE`).
#' @param cluster_cols Logical; whether to cluster columns in the correlation heatmap (default `FALSE`).
#' @param angle_col Numeric; rotation angle for x-axis labels in plots (default `45`).
#' @param display_corr_values Logical; if `TRUE`, displays correlation coefficients on the heatmap (default `FALSE`).
#' @param highlight_outliers Logical; if `TRUE`, identifies and marks library size outliers (default `TRUE`).
#'
#' @return
#' Returns an object of class `"rnaQC"` containing:
#' \describe{
#'   \item{data_ready}{Expression matrix used for QC.}
#'   \item{corr_matrix_samples}{Correlation matrix of samples.}
#'   \item{plot_density}{Density plot object.}
#'   \item{plot_expr_boxplot}{Expression boxplot object.}
#'   \item{plot_lib_boxplot}{Library size boxplot object.}
#'   \item{summary}{List with library size, genes detected, outliers, and percent genes detected per group.}
#' }
#' 
#' @importFrom ggplot2 ggplot aes geom_bar geom_text geom_boxplot geom_density theme_minimal labs element_text theme
#' @importFrom tidyr pivot_longer
#' @importFrom pheatmap pheatmap
#' @export

rna.qc <- function(normalized_data,
                   group_col = "Group",
                   cluster_rows = FALSE, 
                   cluster_cols = FALSE,
                   angle_col = 45,
                   display_corr_values = FALSE,
                   highlight_outliers = TRUE) {
  
  # --- Dependencies ---
  required_pkgs <- c("ggplot2", "pheatmap", "tidyr")
  missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs)) stop("Please install packages: ", paste(missing_pkgs, collapse = ", "))
  
  # --- Input validation ---
  if (!is.list(normalized_data) || is.null(normalized_data$expr_matrix)) {
    stop("Input must be a list returned by rna.normalize() with an 'expr_matrix'.")
  }
  if (is.null(normalized_data$metadata)) {
    stop("Input object must contain a 'metadata' slot.")
  }
  
  data_mat <- as.matrix(normalized_data$expr_matrix)
  metadata <- normalized_data$metadata
  if (!"Sample" %in% colnames(metadata)) stop("metadata must contain a 'Sample' column.")
  if (!group_col %in% colnames(metadata)) stop(paste0("metadata must contain a column named '", group_col, "'."))
  
  # --- Align samples ---
  common_samples <- intersect(colnames(data_mat), metadata$Sample)
  if (length(common_samples) == 0) stop("No matching samples between metadata and expression matrix.")
  metadata_sub <- metadata[match(common_samples, metadata$Sample), , drop = FALSE]
  data_mat <- data_mat[, common_samples, drop = FALSE]
  
  # --- Library sizes and genes detected ---
  lib_sizes <- colSums(data_mat)
  genes_detected <- colSums(data_mat > 0)
  
  # --- Detect library size outliers ---
  if (highlight_outliers) {
    q <- quantile(lib_sizes, probs = c(0.25, 0.75))
    iqr <- diff(q)
    outliers_lib <- which(lib_sizes < (q[1] - 1.5 * iqr) | lib_sizes > (q[2] + 1.5 * iqr))
  } else {
    outliers_lib <- integer(0)
  }
  
  # --- Percent of genes detected per group ---
  percent_genes_detected <- sapply(unique(metadata_sub[[group_col]]), function(grp) {
    samples_in_group <- metadata_sub$Sample[metadata_sub[[group_col]] == grp]
    mean(colSums(data_mat[, samples_in_group, drop = FALSE] > 0) / nrow(data_mat)) * 100
  })
  
  # --- Prepare long format for plots ---
  df_long <- as.data.frame(t(data_mat))
  df_long$Sample <- rownames(df_long)
  df_long <- tidyr::pivot_longer(df_long, cols = -Sample, names_to = "Gene", values_to = "Expression")
  df_long <- merge(df_long, metadata_sub, by = "Sample", all.x = TRUE, sort = FALSE)
  
  # --- Library size barplot ---
  df_box <- data.frame(Sample = colnames(data_mat), Library_size = lib_sizes,
                       Genes_detected = genes_detected, Group = metadata_sub[[group_col]])
  df_box$Outlier <- df_box$Sample %in% names(outliers_lib)
  p_boxplot <- ggplot2::ggplot(df_box, ggplot2::aes(x = Sample, y = Library_size, fill = Group)) +
    ggplot2::geom_bar(stat = "identity", alpha = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = ifelse(Outlier, "*", "")), vjust = -0.5, color = "red", size = 5) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Library size per sample", x = "", y = "Library size") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = angle_col, hjust = 1),
                   legend.position = "bottom")
  print(p_boxplot)
  
  # --- Expression boxplot ---
  p_expr_box <- ggplot2::ggplot(df_long, ggplot2::aes(x = Sample, y = Expression, fill = Sample)) +
    ggplot2::geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Expression distribution per sample", x = "", y = "Expression") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = angle_col, hjust = 1),
                   legend.position = "none")
  print(p_expr_box)
  
  # --- Density plot per group ---
  p_density <- ggplot2::ggplot(df_long, ggplot2::aes(x = Expression,
                                                     color = .data[[group_col]],
                                                     fill = .data[[group_col]])) +
    ggplot2::geom_density(alpha = 0.3, na.rm = TRUE) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Expression density per group", x = "Expression", y = "Density",
                  color = group_col, fill = group_col) +
    ggplot2::theme(legend.position = "bottom")
  print(p_density)
  
  # --- Sample correlation heatmap ---
  corr_mat <- cor(data_mat, use = "pairwise.complete.obs", method = "pearson")
  corr_text <- if (display_corr_values) {
    matrix(paste0("r = ", formatC(round(corr_mat, 2), format = "f", digits = 2)),
           nrow = nrow(corr_mat), ncol = ncol(corr_mat), dimnames = dimnames(corr_mat))
  } else FALSE
  
  pheatmap::pheatmap(corr_mat,
                     main = "Sample correlation heatmap",
                     silent = FALSE,
                     color = grDevices::colorRampPalette(c("navy", "white", "firebrick3"))(50),
                     border_color = NA,
                     clustering_method = "complete",
                     clustering_distance_rows = "correlation",
                     clustering_distance_cols = "correlation",
                     fontsize = 10,
                     legend = TRUE,
                     show_rownames = TRUE,
                     show_colnames = TRUE,
                     treeheight_row = ifelse(cluster_rows, 30, 0),
                     treeheight_col = ifelse(cluster_cols, 30, 0),
                     cluster_rows = cluster_rows,
                     cluster_cols = cluster_cols,
                     display_numbers = corr_text,
                     number_color = "black",
                     angle_col = angle_col)
  
  # --- Return QC object ---
  obj <- list(
    data_ready = data_mat,
    corr_matrix_samples = corr_mat,
    plot_density = p_density,
    plot_expr_boxplot = p_expr_box,
    plot_lib_boxplot = p_boxplot,
    summary = list(
      n_genes = nrow(data_mat),
      n_samples = ncol(data_mat),
      percent_genes_detected = percent_genes_detected,
      library_size = lib_sizes,
      genes_detected = genes_detected,
      outliers_lib = names(outliers_lib)
    )
  )
  class(obj) <- "rnaQC"
  invisible(obj)
}
