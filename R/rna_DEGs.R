#' Extract and visualize top differentially expressed genes (DEGs)
#'
#' @description
#' Identifies the top up- and down-regulated genes from a differential
#' expression analysis performed with \code{rna.compare()}, and optionally
#' generates a bar plot of log2 fold-changes using gene symbols.
#'
#' Gene identifiers are automatically mapped to symbols using \pkg{biomaRt},
#' based on the organism defined in the active \code{rna_project}.
#'
#' @param project \code{rna_project} object created by \code{rna.project()}.
#' @param contrast Identifier of the comparison to use. Can be:
#'   \itemize{
#'     \item A character string matching a stored comparison name
#'     \item A numeric index (position in stored comparisons)
#'     \item \code{NULL} to use the most recent comparison
#'   }
#' @param top_n Integer. Number of genes to select from each tail of the
#'   log2 fold-change distribution (default: 10).
#' @param padj_cutoff Numeric. Adjusted p-value threshold used to define
#'   significant genes (default: 0.05).
#' @param plot Logical. If \code{TRUE}, produces a bar plot of top genes.
#' @param save Logical. Whether to store results in the active
#'   \code{rna_project}.
#' @param verbose Logical. Whether to print informative messages.
#'
#' @details
#' The function ranks genes based on log2 fold-change and extracts the top
#' \code{top_n} most upregulated and downregulated genes.
#'
#' \strong{Contrast interpretation:}
#' The log2 fold-change is defined as \strong{test − reference}, inherited
#' from \code{rna.compare()}. Therefore:
#' \itemize{
#'   \item Positive values indicate higher expression in the test group
#'   \item Negative values indicate higher expression in the reference group
#' }
#'
#' Gene identifiers are mapped to symbols using \pkg{biomaRt}. If mapping fails,
#' original gene IDs are used as fallback.
#'
#' The bar plot displays:
#' \itemize{
#'   \item Top upregulated genes (positive log2FC)
#'   \item Top downregulated genes (negative log2FC)
#'   \item Color-coded regulation direction
#' }
#'
#' @return
#' An object of class \code{"deg_result"} containing:
#' \itemize{
#'   \item \code{comparison}: Comparison identifier
#'   \item \code{contrast}: Vector with \code{test} and \code{reference}
#'   \item \code{deg_table}: Full differential expression table
#'   \item \code{top_genes}: Data frame of selected top genes
#'   \item \code{plot}: \code{ggplot} object (or \code{NULL})
#'   \item \code{summary}: Summary statistics including:
#'     \itemize{
#'       \item Total genes
#'       \item Number of up- and down-regulated genes
#'       \item Number of significant genes
#'       \item Distribution metrics (mean, median, extremes)
#'     }
#'   \item \code{params}: Parameters used (e.g., \code{top_n}, \code{timestamp})
#' }
#'
#' @section Output summary:
#' A summary is printed to the console including:
#' \itemize{
#'   \item Total genes analyzed
#'   \item Number of up- and down-regulated genes
#'   \item Directional bias (percentage of upregulated genes)
#'   \item Number of significant genes (if available)
#'   \item Table of top genes
#' }
#'
#' @examples
#' \dontrun{
#' # Use last comparison
#' rna.degs(project = my_project)
#'
#' # Specify comparison by name
#' rna.degs(my_project,
#'          contrast = "treated_vs_control")
#'
#' # Select more genes
#' rna.degs(my_project,
#'          top_n = 20)
#'
#' # Disable plotting
#' rna.degs(my_project,
#'          plot = FALSE)
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_col labs theme_minimal scale_fill_manual
#'
#' @export

rna.degs <- function(project,
                     contrast = NULL,
                     top_n = 10,
                     padj_cutoff = 0.05,
                     plot = TRUE,
                     save = TRUE,
                     verbose = TRUE
) {

  # ===========================================================================
  # 0) Basic checks
  # ===========================================================================
  .check_dependencies("ggplot2")
  .check_dependencies("biomaRt", bioc = TRUE)

  # ===========================================================================
  # 1) Get active project
  # ===========================================================================
  proj <- project

  comps <- .get_comp(proj)
  gene_map <- .get_gene_annotation(proj)

  # ===========================================================================
  # 2) Validate input
  # ===========================================================================
  contrast <- .get_last_or_selected(
    comps,
    contrast,
    what = "comparison"
  )

  comp_obj <- .get_comp_obj(proj, contrast)

  df <- comp_obj$res

  contrast <- c(
    comp_obj$groups$test,
    comp_obj$groups$reference
  )

  if (!"log2FoldChange" %in% colnames(df)) {
    stop("Comparison result does not contain 'log2FoldChange'.")
  }


  # --- Gene ID sanity check ---
  gene_ids <- rownames(df)

  if (any(grepl("\\.", gene_ids))) {
    stop(
      "[rna.degs] Gene IDs contain version suffixes (e.g. ENSG...\\.6).\n",
      "This function expects cleaned Ensembl IDs.\n",
      "Please run rna.normalize(clean_gene_versions = TRUE) upstream."
    )
  }

  # --- Retrieve annotation ---
  df$Symbol <- gene_map$symbol[match(rownames(df), gene_map$gene_id)]

  df$Symbol[is.na(df$Symbol)] <- rownames(df)[is.na(df$Symbol)]

  # ===========================================================================
  # 3) Select top up/down genes
  # ===========================================================================
  df_sorted <- df[order(df$log2FoldChange), ]
  top_down <- utils::head(df_sorted, top_n)
  top_up   <- utils::tail(df_sorted, top_n)
  top_genes <- rbind(top_down, top_up)

  top_genes$Regulation <- ifelse(top_genes$log2FoldChange > 0, "Upregulated", "Downregulated")

  if ("padj" %in% colnames(df)) {
    sig <- df[df$padj < padj_cutoff & !is.na(df$padj), ]
    n_significant <- nrow(sig)
  } else {
    n_significant <- NA_integer_
  }

  n_up = sum(df$log2FoldChange > 0, na.rm = TRUE)
  n_down = sum(df$log2FoldChange < 0, na.rm = TRUE)

  # ===========================================================================
  # 4) Plot
  # ===========================================================================
  p <- NULL
  if (plot) {
    if (verbose) message(sprintf("[rna.degs] Plotting top %d up- and down-regulated genes...", top_n))
    p <- ggplot2::ggplot(top_genes, ggplot2::aes(
      y = reorder(Symbol, log2FoldChange),
      x = log2FoldChange,
      fill = Regulation
    )) +
      ggplot2::geom_col(position = "dodge") +
      ggplot2::labs(
        y = "",
        x = "log2 fold-change",
        title = sprintf("Top %d Up- and Down-regulated Genes", top_n)
      ) +
      ggplot2::scale_fill_manual(values = c("Upregulated" = "red", "Downregulated" = "blue")) +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "bottom")

    print(p)
  }

  # ===========================================================================
  # 5) Output
  # ===========================================================================
  params <- list(
    top_n = top_n,
    comparison = contrast,
    timestamp = Sys.time()
  )

  summary_metrics <- list(
    total_genes = nrow(df),
    n_up = n_up,
    n_down = n_down,
    mean_log2fc = mean(df$log2FoldChange, na.rm = TRUE),
    median_log2fc = median(df$log2FoldChange, na.rm = TRUE),
    max_log2fc = max(df$log2FoldChange, na.rm = TRUE),
    min_log2fc = min(df$log2FoldChange, na.rm = TRUE),
    n_significant = n_significant,
    padj_cutoff = if ("padj" %in% colnames(df)) 0.05 else NA_real_
  )

  obj <- list(
    timestamp = Sys.time(),
    comparison = contrast,
    contrast = contrast,
    params = params,
    deg_table = df,
    top_genes = top_genes,
    plot = p
  )

  obj$summary <- summary_metrics

  class(obj) <- "deg_result"

  # ===========================================================================
  # 6) Attach to project
  # ===========================================================================
  if (save) {

    proj <- .attach_to_project(
      proj,
      obj,
      slot = "analyses",
      subtype = "degs",
      prefix = "degs",
      log = list(
        comparison = contrast,
        contrast = contrast,
        n_up = n_up,
        n_down = n_down,
        n_significant = n_significant
      )
    )
  }

  # ===========================================================================
  # 7) Return
  # ===========================================================================

  skewness <- mean(df$log2FoldChange > 0, na.rm = TRUE)

  .print_header("Differential Expression Summary")

  .print_block("Overview", function() {
    cat("Comparison:        ", paste(contrast, collapse = " vs "), "\n")
    cat("Total genes:       ", nrow(df), "\n")
    cat("Upregulated:       ", n_up, "\n")
    cat("Downregulated:     ", n_down, "\n")
    cat("Directional bias:  ", round(skewness * 100, 2), "% genes upregulated\n")
  })

  if ("padj" %in% colnames(df)) {

    .print_block("Significant genes", function() {
      cat("padj cutoff:       ", padj_cutoff, "\n")
      cat("Significant:       ", nrow(sig), "\n")
      cat("Up (sig):          ", sum(sig$log2FoldChange > 0), "\n")
      cat("Down (sig):        ", sum(sig$log2FoldChange < 0), "\n")
    })
  }

  .print_block(paste0("Top ", top_n, " genes (by log2FC)"), function() {
    print(top_genes[, c("Symbol", "log2FoldChange", "Regulation")])
  })

  return(invisible(proj))

}
