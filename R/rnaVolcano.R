#' Generate a Volcano Plot of Differential Expression
#'
#' @description
#' Performs a gene-wise t-test between two groups from a `normalized_data` object
#' and produces a volcano plot with optional labeling.
#'
#' @param normalized_data A `normalized_data` object from `rna.normalize()`.
#' @param group_col Character; column in metadata defining groups (default `"Group"`).
#' @param padj_threshold Numeric; adjusted p-value threshold for significance (default `0.05`).
#' @param log2fc_threshold Numeric; absolute log2 fold-change threshold (default `1`).
#' @param results Logical; if TRUE, returns differentially expressed genes and full results (default `FALSE`).
#' @param identify Logical; if TRUE, labels genes on the volcano plot (default `FALSE`).
#'
#' @return Invisibly returns a list with up-regulated genes, down-regulated genes, and full results.
#' @export

rna.volcano <- function(normalized_data,
                        group_col = "Group",
                        padj_threshold = 0.05,
                        log2fc_threshold = 1,
                        results = FALSE,
                        identify = FALSE) {

  # --- Check input ---
  if (!inherits(normalized_data, "normalized_data")) stop("Input must be of class 'normalized_data'.")
  if (!"metadata" %in% names(normalized_data)) stop("'normalized_data' must contain a 'metadata' data.frame.")

  expr_mat <- as.matrix(normalized_data$expr_matrix)
  metadata <- normalized_data$metadata

  if (!"Sample" %in% colnames(metadata)) stop("'metadata' must contain a 'Sample' column.")
  if (!group_col %in% colnames(metadata)) stop(sprintf("Group column '%s' not found in metadata.", group_col))

  groups <- unique(metadata[[group_col]])
  if (length(groups) != 2) stop("Currently only supports exactly 2 groups.")

  g1_samples <- metadata$Sample[metadata[[group_col]] == groups[1]]
  g2_samples <- metadata$Sample[metadata[[group_col]] == groups[2]]

  missing_samples <- setdiff(c(g1_samples, g2_samples), colnames(expr_mat))
  if (length(missing_samples) > 0) stop("Samples not found in normalized_data: ", paste(missing_samples, collapse = ", "))

  # --- Subset expression matrix ---
  expr_sub <- expr_mat[, c(g1_samples, g2_samples), drop = FALSE]
  gene_ids <- rownames(expr_sub)

  # --- Compute t-tests and log2 fold-changes ---
  res_df <- data.frame(
    Gene = gene_ids,
    log2FoldChange = NA_real_,
    pvalue = NA_real_
  )

  for (i in seq_len(nrow(expr_sub))) {
    vals1 <- as.numeric(expr_sub[i, g1_samples])
    vals2 <- as.numeric(expr_sub[i, g2_samples])
    ttest <- try(stats::t.test(vals2, vals1), silent = TRUE)

    if (inherits(ttest, "try-error")) {
      res_df$pvalue[i] <- NA
      res_df$log2FoldChange[i] <- NA
    } else {
      res_df$pvalue[i] <- ttest$p.value
      res_df$log2FoldChange[i] <- mean(vals2, na.rm = TRUE) - mean(vals1, na.rm = TRUE)
    }
  }

  res_df$padj <- stats::p.adjust(res_df$pvalue, method = "BH")
  res_df <- dplyr::filter(res_df, !is.na(pvalue))

  # --- Identify DE genes ---
  up <- res_df$Gene[res_df$padj < padj_threshold & res_df$log2FoldChange > log2fc_threshold]
  down <- res_df$Gene[res_df$padj < padj_threshold & res_df$log2FoldChange < -log2fc_threshold]

  cat("=== Volcano summary ===\n")
  cat("Groups:", groups[1], "vs", groups[2], "\n")
  cat("Up-regulated (", groups[2], "):", length(up), "\n")
  cat("Down-regulated (", groups[1], "):", length(down), "\n")
  cat("padj threshold:", padj_threshold, "| log2FC threshold:", log2fc_threshold, "\n")
  cat("============================\n")
  if (!results) {
    cat("Tip: To retrieve DE genes, use 'results = TRUE'.\n")
    cat("Use 'identify = TRUE' to label genes on the volcano plot.\n\n")
  }

  # --- Prepare plot ---
  res_df$group_color <- factor(
    ifelse(res_df$padj < padj_threshold & res_df$log2FoldChange > log2fc_threshold, "Up",
           ifelse(res_df$padj < padj_threshold & res_df$log2FoldChange < -log2fc_threshold, "Down", "NS"))
  )

  p <- ggplot2::ggplot(res_df, ggplot2::aes(x = log2FoldChange, y = -log10(pvalue), color = group_color)) +
    ggplot2::geom_point(size = 2, na.rm = TRUE) +
    ggplot2::scale_color_manual(values = c("Up" = "#ff3333", "Down" = "#006699", "NS" = "darkgrey")) +
    ggplot2::labs(color = "Regulation",
                  x = "log2 Fold Change",
                  y = "-log10(p-value)",
                  title = paste0(groups[2], " vs ", groups[1])) +
    ggplot2::theme_minimal(base_size = 16) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1),
                   plot.background = ggplot2::element_rect(fill = "white", color = NA),
                   panel.background = ggplot2::element_rect(fill = "white", color = NA))

  if (identify) {
    if (!requireNamespace("ggrepel", quietly = TRUE)) stop("Package 'ggrepel' needed for identify = TRUE.")
    p <- p + ggrepel::geom_text_repel(ggplot2::aes(label = Gene), size = 3, max.overlaps = 50, na.rm = TRUE)
  }

  print(p)

  # --- Return ---
  volcano_degs <- list(up = up, down = down, full_results = res_df)
  if (results) return(volcano_degs)
  invisible(volcano_degs)
}
