#' Functional enrichment analysis (GO) for differential expression results
#'
#' This function performs Gene Ontology (GO) enrichment analysis using the
#' \pkg{clusterProfiler} package for up- and down-regulated genes separately.
#' It supports human, mouse and zebrafish datasets and automatically creates dot plots
#' for significant terms.
#'
#' @param comp_data A result object from [rna_compare()] or a list containing
#'   differential expression results with a column `log2FoldChange`.
#' @param species Character. Either `"human"`, `"mouse"` or `"zebrafish"` (default: `"mouse"`).
#' @param padj_cutoff Numeric. Adjusted p-value threshold for significance
#'   (default: `0.05`).
#' @param log2fc_cutoff Numeric. Absolute log2 fold-change threshold for gene
#'   selection (default: `1`).
#' @param ont Character. Ontology type passed to
#'   [clusterProfiler::enrichGO()], one of `"BP"`, `"MF"`, or `"CC"`
#'   (default: `"BP"`).
#' @param plot Logical. Whether to display enrichment dot plots (default: `TRUE`).
#' @param top_terms Integer. Number of top GO terms to display in plots (default: `10`).
#'
#' @return A list containing:
#' \item{up}{A `clusterProfiler` enrichment result for up-regulated genes.}
#' \item{down}{A `clusterProfiler` enrichment result for down-regulated genes.}
#'
#' @details
#' This function filters genes based on adjusted p-values and log2 fold-change
#' thresholds, then performs GO enrichment separately for up- and
#' down-regulated gene sets using Entrez ID mapping.
#'
#' It supports human (`org.Hs.eg.db`), mouse (`org.Mm.eg.db`) and zebra fish (`org.Dr.eg.db`)
#' annotation databases.
#'
#' @examples
#' \dontrun{
#' enrich_res <- rna.enrich(comp_data, species = "human")
#' enrich_res$up
#' enrich_res$down
#' }
#'
#' @export
rna.enrich <- function(
    comp_data,
    species = "mouse",
    padj_cutoff = 0.05,
    log2fc_cutoff = 1,
    ont = "BP",
    plot = TRUE,
    top_terms = 10
) {
  # --- Required packages ---
  required_pkgs <- c("clusterProfiler", "enrichplot", "ggplot2")
  if (species == "human") required_pkgs <- c(required_pkgs, "org.Hs.eg.db")
  if (species == "mouse") required_pkgs <- c(required_pkgs, "org.Mm.eg.db")
  if (species == "zebrafish") required_pkgs <- c(required_pkgs, "org.Dr.eg.db")

  missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs)) {
    stop(sprintf("Please install the following packages: %s",
                 paste(missing_pkgs, collapse = ", ")))
  }

  # --- Select organism database ---
  OrgDb <- switch(
    species,
    human = org.Hs.eg.db::org.Hs.eg.db,
    mouse = org.Mm.eg.db::org.Mm.eg.db,
    zebrafish = org.Dr.eg.db::org.Dr.eg.db,
    stop("`species` must be either 'human', 'mouse' or 'zebra fish'.")
  )

  # --- Extract results ---
  if (inherits(comp_data, "rnaCompare")) {
    res <- comp_data$res
  } else if (is.list(comp_data) && "comparisons" %in% names(comp_data) &&
             "results" %in% names(comp_data$comparisons)) {
    res <- comp_data$comparisons$results
  } else {
    stop("Input must be an rnaCompare result or contain $comparisons$results.")
  }

  res <- as.data.frame(res)

  # --- Filter up- and down-regulated genes ---
  up_genes <- res[!is.na(res$padj) &
                    res$padj < padj_cutoff &
                    res$log2FoldChange > log2fc_cutoff, ]

  down_genes <- res[!is.na(res$padj) &
                      res$padj < padj_cutoff &
                      res$log2FoldChange < -log2fc_cutoff, ]

  if (nrow(up_genes) == 0 && nrow(down_genes) == 0) {
    message("[rna_enrich] No significant genes found with the given thresholds.")
    return(invisible(NULL))
  }

  # --- Helper function for enrichment + plotting ---
  enrich_and_plot <- function(genes_df, label) {
    if (nrow(genes_df) == 0) {
      message(sprintf("[rna_enrich] No significant %s genes found.", label))
      return(NULL)
    }

    entrez_ids <- clusterProfiler::bitr(
      rownames(genes_df),
      fromType = "ENSEMBL",
      toType = "ENTREZID",
      OrgDb = OrgDb
    )

    if (nrow(entrez_ids) == 0) {
      message(sprintf("[rna_enrich] No %s genes mapped to ENTREZID.", label))
      return(NULL)
    }

    ego <- clusterProfiler::enrichGO(
      gene = entrez_ids$ENTREZID,
      OrgDb = OrgDb,
      keyType = "ENTREZID",
      ont = ont,
      pAdjustMethod = "BH",
      qvalueCutoff = 0.05,
      readable = TRUE
    )

    if (plot && nrow(ego@result) > 0) {
      p <- enrichplot::dotplot(ego, showCategory = top_terms) +
        ggplot2::ggtitle(sprintf("GO Enrichment (%s genes)", label)) +
        ggplot2::theme_minimal()
      print(p)
    }

    invisible(ego)
  }

  # --- Run enrichment for both directions ---
  ego_up <- enrich_and_plot(up_genes, "up-regulated")
  ego_down <- enrich_and_plot(down_genes, "down-regulated")

  invisible(list(up = ego_up, down = ego_down))
}
