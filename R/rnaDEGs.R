#' Identify and visualize top up- and down-regulated genes with gene symbols
#'
#' This function extracts the top up- and down-regulated genes from a
#' differential expression results table (typically from DESeq2 via rna.compare)
#' and optionally produces a bar plot showing their log2 fold changes using gene symbols.
#'
#' @param comp_data A data frame with DE results (from rna.compare), or a list with `res`.
#'   Must include a column named `log2FoldChange`.
#' @param top_n Integer. Number of genes from each end of fold-change distribution (default: 10).
#' @param plot Logical. If TRUE, produces a bar plot of top genes (default: TRUE).
#' @param verbose Logical. If TRUE, prints informative messages (default: TRUE).
#'
#' @return A list with:
#' \item{top_genes}{Data frame of top up- and down-regulated genes.}
#' \item{plot}{A `ggplot` object (or NULL if plot = FALSE).}
#'
#' @export

rna.degs <- function(
    comp_data,
    top_n = 10,
    plot = TRUE,
    verbose = TRUE
) {

  # --- Dependencies ---
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.")
  if (!requireNamespace("biomaRt", quietly = TRUE)) stop("Package 'biomaRt' is required.")

  # --- Extract DE table ---
  if (is.list(comp_data) && "res" %in% names(comp_data)) {
    df <- comp_data$res
  } else if (is.data.frame(comp_data)) {
    df <- comp_data
  } else stop("`comp_data` must be a data frame or list with 'res' element.")

  if (!"log2FoldChange" %in% colnames(df))
    stop("Input must contain 'log2FoldChange' column.")

  # --- Determine organism from rownames prefix ---
  gene_ids <- rownames(df)
  organism <- if (all(grepl("^ENSG", gene_ids))) "human" else
    if (all(grepl("^ENSMUSG", gene_ids))) "mouse" else
      if (all(grepl("^ENSDARG", gene_ids))) "zebrafish" else "human"
  if (verbose) message("[rna.degs] Organism inferred as: ", organism)

  # --- Attempt to retrieve gene symbols via BioMart ---
  ensembl_dataset <- switch(
    organism,
    "human"     = "hsapiens_gene_ensembl",
    "mouse"     = "mmusculus_gene_ensembl",
    "zebrafish" = "drerio_gene_ensembl"
  )

  ensembl <- tryCatch(
    biomaRt::useEnsembl(biomart = "genes", dataset = ensembl_dataset, mirror = "www"),
    error = function(e) biomaRt::useEnsembl(biomart = "genes", dataset = ensembl_dataset, mirror = "uswest")
  )

  symbol_attr <- switch(
    organism,
    "human"     = "hgnc_symbol",
    "mouse"     = "mgi_symbol",
    "zebrafish" = "zfin_id_symbol"
  )

  # Remove version if present
  gene_ids_clean <- sub("\\..*$", "", gene_ids)

  gene_map <- tryCatch({
    biomaRt::getBM(
      attributes = c("ensembl_gene_id", symbol_attr),
      filters = "ensembl_gene_id",
      values = gene_ids_clean,
      mart = ensembl
    )
  }, error = function(e) NULL)

  gene_map$ensembl_gene_id <- as.character(gene_map$ensembl_gene_id)
  gene_map$symbol <- if (symbol_attr %in% colnames(gene_map)) gene_map[[symbol_attr]] else NA_character_

  # Map symbols back to DE table
  df$Symbol <- gene_map$symbol[match(gene_ids_clean, gene_map$ensembl_gene_id)]
  df$Symbol[is.na(df$Symbol)] <- rownames(df)[is.na(df$Symbol)]

  # --- Select top up/down genes ---
  df_sorted <- df[order(df$log2FoldChange), ]
  top_down <- utils::head(df_sorted, top_n)
  top_up   <- utils::tail(df_sorted, top_n)
  top_genes <- rbind(top_down, top_up)

  top_genes$Regulation <- ifelse(top_genes$log2FoldChange > 0, "Upregulated", "Downregulated")

  # --- Plot ---
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

  invisible(list(top_genes = top_genes, plot = p))
}
