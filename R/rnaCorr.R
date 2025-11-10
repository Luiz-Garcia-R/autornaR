#' Compute and visualize sample correlations from normalized RNA-seq data
#'
#' This function calculates sample-to-sample correlation matrices and optionally
#' produces a correlation heatmap and a group-wise scatterplot of mean expression
#' values. The correlation method can be automatically determined based on
#' data normality or manually specified.
#'
#' @param normalized_data A `normalized_data` object containing an expression
#'   matrix (`expr_matrix`) and metadata (`metadata`).
#' @param method Character. Correlation method: `"pearson"`, `"spearman"`,
#'   `"kendall"`, or `"auto"` for automatic selection based on normality tests
#'   (default: `"auto"`).
#' @param plot Logical. If `TRUE`, plots a correlation heatmap and, when
#'   applicable, a scatterplot (default: `TRUE`).
#' @param verbose Logical. If `TRUE`, prints progress and method selection
#'   messages (default: `TRUE`).
#'
#' @return A numeric correlation matrix.
#'
#' @details
#' The function first checks whether all samples show normal distributions using
#' the Shapiro–Wilk test and whether there are ties. If all samples appear
#' normally distributed and no ties are present, Pearson correlation is used;
#' otherwise, Spearman correlation is applied.
#'
#' If `plot = TRUE`, two plots are produced:
#' \itemize{
#'   \item A heatmap showing all pairwise sample correlations.
#'   \item A scatterplot comparing group-wise mean expression values if two
#'   groups are available.
#' }
#'
#' @examples
#' \dontrun{
#' corr_mat <- rna.corr(normalized_data)
#' }
#'
#' @export
rna.corr <- function(
    normalized_data,
    method = "auto",
    plot = TRUE,
    verbose = TRUE
) {
  # --- Required packages ---
  required_pkgs <- c("dplyr", "ggplot2", "pheatmap", "tidyr")
  for (pkg in required_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop(sprintf("Package '%s' is required but not installed.", pkg))
  }

  # --- Validate input ---
  if (!inherits(normalized_data, "normalized_data"))
    stop("'normalized_data' must be of class 'normalized_data'.")

  if (!"metadata" %in% names(normalized_data))
    stop("'normalized_data' must contain a 'metadata' element.")

  expr_mat <- as.matrix(normalized_data$expr_matrix)
  gene_ids <- rownames(expr_mat)
  metadata <- normalized_data$metadata

  if (!all(c("Sample", "Group") %in% colnames(metadata)))
    stop("Metadata must contain 'Sample' and 'Group' columns.")

  # --- Correlation method selection ---
  if (method == "auto") {
    normality <- apply(expr_mat, 2, function(x)
      tryCatch(stats::shapiro.test(x)$p.value > 0.05, error = function(e) FALSE))
    has_ties <- any(apply(expr_mat, 2, function(x) anyDuplicated(x) > 0))
    method_used <- if (all(normality) && !has_ties) "pearson" else "spearman"

    if (verbose)
      message(sprintf("[rna_corr] Automatically selected correlation method: %s",
                      method_used))
  } else {
    method_used <- method
  }

  # --- Compute correlation matrix ---
  corr_mat <- stats::cor(
    expr_mat,
    method = method_used,
    use = "pairwise.complete.obs"
  )
  corr_text <- matrix(
    paste0("r = ", round(corr_mat, 2)),
    nrow = nrow(corr_mat),
    ncol = ncol(corr_mat),
    dimnames = dimnames(corr_mat)
  )

  # --- Heatmap visualization ---
  if (plot && ncol(expr_mat) > 1) {
    if (verbose) message("[rna_corr] Generating correlation heatmap...")
    pheatmap::pheatmap(
      corr_mat,
      display_numbers = corr_text,
      number_color = "black",
      main = sprintf("Correlation Heatmap (%s)", method_used),
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      breaks = seq(-1, 1, length.out = 101)
    )
  }

  # --- Group-wise scatterplot ---
  df_long <- tidyr::pivot_longer(
    data.frame(Gene = gene_ids, expr_mat, check.names = FALSE),
    -Gene,
    names_to = "Sample",
    values_to = "Expression"
  )

  df_long <- dplyr::left_join(df_long, metadata, by = "Sample")

  corr_df_scatter <- df_long |>
    dplyr::group_by(Gene, Group) |>
    dplyr::summarise(Expression = mean(Expression, na.rm = TRUE), .groups = "drop") |>
    tidyr::pivot_wider(names_from = Group, values_from = Expression) |>
    as.data.frame()

  corr_df_scatter <- corr_df_scatter[, -1, drop = FALSE]  # remove 'Gene' column

  if (plot && ncol(corr_df_scatter) == 2) {
    if (verbose) message("[rna_corr] Generating group-level scatterplot...")

    cols <- colnames(corr_df_scatter)
    test <- stats::cor.test(
      corr_df_scatter[[1]],
      corr_df_scatter[[2]],
      method = method_used
    )

    r_val <- round(test$estimate, 3)
    p_val <- if (test$p.value < 0.001) "<0.001" else signif(test$p.value, 3)
    strength <- dplyr::case_when(
      abs(r_val) < 0.3 ~ "very weak or none",
      abs(r_val) < 0.5 ~ "weak",
      abs(r_val) < 0.7 ~ "moderate",
      abs(r_val) < 0.9 ~ "strong",
      TRUE ~ "very strong"
    )

    g <- ggplot2::ggplot(
      corr_df_scatter,
      ggplot2::aes(x = .data[[cols[1]]], y = .data[[cols[2]]])
    ) +
      ggplot2::geom_point(alpha = 0.6, color = "steelblue") +
      ggplot2::geom_smooth(
        method = ifelse(method_used == "pearson", "lm", "loess"),
        se = FALSE,
        color = "red",
        linetype = "dashed"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = sprintf("Correlation Scatterplot (%s)", method_used),
        subtitle = sprintf("r = %s | p = %s | %s", r_val, p_val, strength),
        x = cols[1],
        y = cols[2]
      )

    print(g)
  }

  invisible(corr_mat)
}
