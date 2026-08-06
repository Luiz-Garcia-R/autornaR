#' Volcano plot for differential expression results
#'
#' @description
#' Generates a volcano plot from differential expression results stored
#' in the active \code{rna_project}. The function visualizes the relationship
#' between effect size (log2 fold-change) and statistical significance
#' (-log10 p-value), highlighting significantly up- and down-regulated genes.
#'
#' @details
#' This function retrieves results produced by \code{rna.compare()} and does not
#' perform differential expression analysis itself.
#'
#' Genes are classified into three categories:
#' \itemize{
#'   \item Upregulated genes (significant and above log2FC threshold)
#'   \item Downregulated genes (significant and below negative log2FC threshold)
#'   \item Non-significant genes
#' }
#'
#' Significance is defined using the adjusted p-value threshold (\code{padj_threshold}),
#' while biological relevance is controlled via the log2 fold-change threshold
#' (\code{log2fc_threshold}).
#'
#' The function supports optional labeling of genes directly on the plot using
#' \code{ggrepel}, which is useful for highlighting key markers.
#'
#' If no comparison is specified, the most recent comparison stored in the
#' \code{rna_project} is used automatically.
#'
#' @param project \code{rna_project} object created by \code{rna.project()}.
#' @param comparison Character. ID of the comparison to visualize.
#'   If \code{NULL}, the most recent comparison is used.
#' @param genes Optional character vector of genes to label on the plot.
#'   Gene identifiers can be provided either as gene symbols (e.g. \code{"TP53"})
#'   or as gene IDs matching the rownames of the differential expression results
#'   (e.g. Ensembl IDs). Matching is performed automatically using the project
#'   annotation, and both representations are supported transparently.
#'   When provided, this argument takes precedence over \code{top_genes}.
#' @param top_genes Optional integer specifying the number of genes to label
#'   automatically based on statistical significance and effect size.
#'   Genes are ranked primarily by adjusted p-value (\code{padj}) and secondarily
#'   by absolute log2 fold-change, prioritizing genes that are both highly
#'   significant and strongly regulated. Ignored if \code{genes} is provided.
#' @param padj_threshold Numeric. Adjusted p-value cutoff for statistical
#'   significance (default \code{0.05}).
#' @param log2fc_threshold Numeric. Absolute log2 fold-change cutoff for
#'   biological relevance (default \code{0.5}).
#' @param colors Character. Optional points colors.
#' @param save Logical. Whether to store results in the active \code{rna_project}
#'   (default \code{TRUE}).
#' @param verbose Logical. If \code{TRUE}, prints the plot and summary to the console
#'   (default \code{TRUE}).
#'
#' @return
#' An object of class \code{"rna_volcano"} containing:
#' \describe{
#'   \item{plot}{The ggplot object representing the volcano plot.}
#'   \item{thresholds}{List with \code{padj} and \code{log2fc} thresholds used.}
#'   \item{DEGs}{List with vectors of up- and down-regulated genes.}
#'   \item{comparison}{Comparison ID used to generate the plot.}
#' }
#'
#' @section \strong{Visualization features}:
#' \itemize{
#'   \item Color-coded points (Up, Down, Non-significant)
#'   \item Horizontal line indicating adjusted p-value threshold
#'   \item Vertical lines indicating log2 fold-change thresholds
#'   \item Optional gene labeling for interpretation
#' }
#'
#' @section \strong{Interpretation guidelines}:
#' \itemize{
#'   \item Genes in the upper left and upper right quadrants are the most significant
#'   and biologically relevant.
#'   \item Symmetry around zero log2FC may indicate balanced regulation across conditions.
#'   \item Large numbers of significant genes may reflect strong biological signal or
#'   potential confounding factors (e.g., batch effects).
#' }
#'
#' @section \strong{Side effects}:
#' If \code{save = TRUE}, the result is stored in the active \code{rna_project}
#' under \code{analyses$volcano}, along with metadata describing thresholds
#' and comparison used.
#'
#' @examples
#' \dontrun{
#'
#' # Generate basic volcano plot
#' rna.volcano(project = my_project)
#'
#' # Specific gene labeling
#' rna.volcano(my_project,
#'             genes = c('gene1', 'gene2', 'gene3')
#' )
#'
#' # Top gene labeling
#' rna.volcano(my_project,
#'             top_genes = 10
#' )
#'
#' # Custom thresholds
#' rna.volcano(my_project,
#'             padj_threshold = 0.01,
#'             log2fc_threshold = 1
#' )
#'
#' }
#'
#' @importFrom grDevices dev.flush
#'
#' @seealso
#' \code{\link{rna.compare}}, \code{\link{rna.sets}}
#'
#' @export

rna.volcano <- function(project,
                        comparison = NULL,
                        genes = NULL,
                        top_genes = NULL,
                        padj_threshold = 0.05,
                        log2fc_threshold = 0.5,
                        colors = NULL,
                        save = TRUE,
                        verbose = TRUE) {

  # ===========================================================================
  # 0) Basic checks
  # ===========================================================================
  pkgs <- c("ggrepel", "ggplot2")

  .check_dependencies(pkgs)

  # ===========================================================================
  # 1) Get active project
  # ===========================================================================
  proj <- project

  comp_obj <- .get_comp_obj(proj, comparison)

  # ===========================================================================
  # 2) Validate input
  # ===========================================================================
  if (!is.null(genes) && !is.null(top_genes)) {
    stop("Use either 'genes' or 'top_genes', not both.")
  }

  if (is.null(proj$analyses$comparison) ||
      length(proj$analyses$comparison) == 0) {
    stop("No comparisons found. Run rna.compare() first.")
  }

  res_df <- comp_obj$res
  groups <- c(comp_obj$groups$test, comp_obj$groups$reference)

  res_df$Gene <- rownames(res_df)
  res_df$highlight <- FALSE
  label_df <- NULL

  required_cols <- c("log2FoldChange", "pvalue", "padj")

  if (!all(required_cols %in% colnames(res_df))) {
    stop("Comparison results missing required columns: ",
         paste(setdiff(required_cols, colnames(res_df)), collapse = ", "))
  }

  # Define comparison
  if (is.null(comparison)) {
    comparison_id <- proj$analyses$comparison$last
  } else {
    comparison_id <- comparison
  }

  # ===========================================================================
  # 3) Gene annotation
  # ===========================================================================
  if (is.null(rownames(res_df))) {
    rownames(res_df) <- res_df$Gene
  }

  gene_annotation <- .get_gene_annotation(proj)
  gene_map <- .align_gene_annotation(gene_annotation, res_df)
  res_df$gene_symbol <- gene_map$symbol

  res_df$gene_symbol <- ifelse(
    is.na(res_df$gene_symbol) | res_df$gene_symbol == "",
    rownames(res_df),
    res_df$gene_symbol
  )

  sig <- res_df[!is.na(res_df$padj) & res_df$padj < padj_threshold, ]

  is_up <- !is.na(res_df$padj) &
    res_df$padj < padj_threshold &
    res_df$log2FoldChange > log2fc_threshold

  is_down <- !is.na(res_df$padj) &
    res_df$padj < padj_threshold &
    res_df$log2FoldChange < -log2fc_threshold

  up <- res_df$Gene[is_up]
  down <- res_df$Gene[is_down]

  res_df$group_color <- factor(
    ifelse(is_up, "Up",
           ifelse(is_down, "Down", "NS"))
  )

  # ===========================================================================
  # 4) Label selection
  # ===========================================================================
    # gene symbol -> ENSEMBL
    if (!is.null(genes)) {

      gene_annotation <- .get_gene_annotation(proj)
      gene_map <- .align_gene_annotation(gene_annotation, res_df)

      genes_ensembl <- gene_map$gene_id[
        match(genes, gene_map$symbol)
      ]

      genes_ensembl <- genes_ensembl[!is.na(genes_ensembl)]

      genes_all <- unique(c(genes, genes_ensembl))

      label_df <- res_df[
        res_df$Gene %in% genes_all |
          res_df$gene_symbol %in% genes,
      ]

    }

    # Select top genes
    if (is.null(genes) && !is.null(top_genes)) {

      label_df <- res_df[
        order(res_df$padj, -abs(res_df$log2FoldChange)),
      ][1:min(top_genes, nrow(res_df)), ]

    }

    res_df$highlight <- res_df$Gene %in% label_df$Gene

  # ===========================================================================
  # 5) Plot
  # ===========================================================================
  p <- ggplot2::ggplot(res_df, ggplot2::aes(x = log2FoldChange, y = -log10(padj + 1e-300), color = group_color)) +

      ggplot2::geom_point(
        na.rm = TRUE,
        ggplot2::aes(
          size = .data$highlight,
          alpha = .data$highlight
        )
      ) +

      ggplot2::scale_alpha_manual(values = c(0.4, 1), guide = "none") +

      ggplot2::scale_size_identity() +

      ggplot2::labs(
        color = "Regulation",
        x = "log2 Fold Change",
        y = "-log10(padj)",
        title = paste0(groups[2], " vs ", groups[1])) +

      ggplot2::geom_vline(
        xintercept = c(-log2fc_threshold, log2fc_threshold),
        linetype = "dashed") +

      ggplot2::geom_hline(
        yintercept = -log10(padj_threshold),
        linetype = "dashed") +

      ggplot2::theme_minimal(
        base_size = 12) +

      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5),
        panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.6),
        plot.background = ggplot2::element_rect(fill = "white", color = NA),
        panel.background = ggplot2::element_rect(fill = "white", color = NA))

  # Color
  if (is.null(colors)) {
      p <- p +
        ggplot2::scale_color_manual(
          values = c("Up" = "#ff3333",
                     "Down" = "#006699",
                     "NS" = "darkgrey"))
  } else {
    p <- p +
      ggplot2::scale_color_manual(values = colors)
    }


  if (!is.null(label_df)) {

    p <- p +
      ggrepel::geom_text_repel(
        data = label_df,
        ggplot2::aes(label = .data$gene_symbol),
        size = 3,
        max.overlaps = Inf,
        na.rm = TRUE
      )
  }

  if (verbose) {
    print(p)
    invisible(dev.flush())
  }

  # ===========================================================================
  # 6) Output
  # ===========================================================================
    selected_genes <- NULL

    if (!is.null(label_df)) {
      selected_genes <- list(
        gene_ids = label_df$Gene,
        gene_symbols = label_df$gene_symbol
      )
    }

    obj <- list(
      timestamp = Sys.time(),
      plot = p,
      thresholds = list(padj = padj_threshold, log2fc = log2fc_threshold),
      DEGs = list(up = up, down = down),
      selected_genes = selected_genes,
      comparison = comparison_id
    )

  class(obj) <- "rna_volcano"

  # ===========================================================================
  # 7) Attach to project
  # ===========================================================================
  if (save) {

    proj <- .attach_to_project(
      proj,
      obj,
      slot = "analyses",
      subtype = "volcano",
      prefix = "volc",
      id = comparison_id,
      log = list(
        comparison = comparison_id,
        padj = padj_threshold,
        log2fc = log2fc_threshold
      )
    )
  }

  # ===========================================================================
  # 8) Return
  # ===========================================================================
  if (verbose) {

    .print_header("RNA Volcano Plot")

    .print_block("Summary", function() {
      cat("Comparison:   ", comparison_id, "\n")
      cat("Groups:       ", groups[1], "vs", groups[2], "\n")
      cat("Thresholds:    padj <", padj_threshold,
          "| log2FC >", log2fc_threshold, "\n")
      cat("DEGs:          Up =", length(up),
          "| Down =", length(down), "\n")
    })

    if (!is.null(selected_genes)) {
      .print_block("Highlighted genes", function() {
        n <- length(selected_genes$gene_symbols)
        cat("n =", n, "\n")

        cat(paste(head(selected_genes$gene_symbols, 8), collapse = ", "), "\n")

        if (n > 8) cat("...\n")
      })
    }
  }

  return(invisible(proj))

}

