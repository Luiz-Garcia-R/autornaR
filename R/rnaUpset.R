#' Generate an UpSet plot of expressed genes across groups
#'
#' @description
#' Produces an UpSet plot showing the overlap of expressed genes across experimental groups.
#' Expression is considered present if normalized counts exceed the specified cutoff.
#'
#' @param comp_data A `normalized_data` object (with DESeq2 comparisons), an `rnaCompare` object, or a `DESeqDataSet`.
#' @param group_col Character; the metadata column defining groups (default `"Group"`).
#' @param cutoff Numeric; expression threshold to consider a gene as expressed (default `5`).
#' @param results Logical; if TRUE, returns lists of genes by group, unique genes, and shared genes (default `FALSE`).
#'
#' @return Invisibly returns gene lists if `results = TRUE`; otherwise prints the UpSet plot.
#'
#' @importFrom ComplexUpset upset
#' @importFrom DESeq2 counts
#' @importFrom SummarizedExperiment colData
#' @importFrom ggplot2 labs
#' @export

rna.upset <- function(comp_data, group_col = "Group", cutoff = 5, results = FALSE) {

  # --- Check required packages ---
  required_pkgs <- c("ComplexUpset", "ggplot2", "dplyr", "tidyr")
  missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs)) stop("Please install packages: ", paste(missing_pkgs, collapse = ", "))

  # --- Extract DESeqDataSet ---
  if (inherits(comp_data, "normalized_data")) {
    if (!"comparisons" %in% names(comp_data) || !"dds" %in% names(comp_data$comparisons)) {
      stop("normalized_data does not contain DESeq2 results. Run rna.compare() first.")
    }
    dds <- comp_data$comparisons$dds
  } else if (inherits(comp_data, "rnaCompare")) {
    dds <- comp_data$dds
  } else if (inherits(comp_data, "DESeqDataSet")) {
    dds <- comp_data
  } else {
    stop("Input must be a normalized_data with comparisons, rnaCompare, or DESeqDataSet object.")
  }

  # --- Validate group column ---
  coldata <- SummarizedExperiment::colData(dds)
  if (!group_col %in% colnames(coldata)) {
    stop(paste0("Column '", group_col, "' not found in DESeqDataSet colData."))
  }

  # --- Extract samples and groups ---
  samples <- colnames(dds)
  groups <- coldata[[group_col]]
  unique_groups <- unique(groups)

  # --- Build binary presence matrix ---
  expr_mat <- DESeq2::counts(dds, normalized = TRUE)
  presence_list <- lapply(unique_groups, function(g) {
    samp_g <- samples[groups == g]
    rowMeans(expr_mat[, samp_g, drop = FALSE]) > cutoff
  })

  presence_df <- as.data.frame(presence_list)
  colnames(presence_df) <- unique_groups
  presence_df$gene <- rownames(expr_mat)

  # --- Create UpSet plot ---
  upset_plot <- ComplexUpset::upset(
    presence_df,
    intersect = colnames(presence_df)[1:length(unique_groups)],
    name = "Gene Presence",
    width_ratio = 0.1
  ) + ggplot2::labs(title = paste("UpSet Plot - Genes expressed (cutoff >", cutoff, ")"))

  print(upset_plot)

  # --- Return gene lists if requested ---
  if (results) {
    genes_by_group <- lapply(unique_groups, function(g) rownames(expr_mat)[presence_df[[g]]])
    names(genes_by_group) <- unique_groups

    genes_unique <- lapply(unique_groups, function(g) {
      setdiff(genes_by_group[[g]], unlist(genes_by_group[unique_groups != g]))
    })

    genes_shared <- Reduce(intersect, genes_by_group)

    return(invisible(list(
      genes_by_group = genes_by_group,
      genes_unique = genes_unique,
      genes_shared = genes_shared
    )))
  }

  invisible(NULL)
}
