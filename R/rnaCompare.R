#' Compare RNA-seq groups using DESeq2
#'
#' Performs differential expression analysis between two groups
#' in a normalized RNA-seq dataset using the DESeq2 package.
#'
#' @param normalized_data A `normalized_data` object containing an expression
#'   matrix (`expr_matrix`) and associated metadata (`metadata`).
#' @param group_col Character. Column name in `metadata` defining the comparison
#'   groups (default: `"Group"`).
#' @param method Character. Currently only `"deseq2"` is supported.
#' @param batch_col Optional character. Column name specifying batch information
#'   (default: `"Batch"`). If not found, batch correction is skipped.
#' @param clean_ids Logical. Whether to clean Ensembl gene IDs by removing
#'   version suffixes (default: `TRUE`).
#' @param assign_result Logical. If `TRUE`, assigns the result object to the
#'   parent environment (default: `TRUE`).
#' @param assign_name Character. Name to assign the output object if
#'   `assign_result = TRUE` (default: `"comp_data"`).
#' @param plot Logical. If `TRUE`, shows basic DESeq2 diagnostic plots
#'   (default: `TRUE`).
#' @param verbose Logical. If `TRUE`, prints progress messages
#'   (default: `TRUE`).
#' @param envir Environment. Environment where the output will be assigned
#'   if `assign_result = TRUE` (default: `parent.frame()`).
#'
#' @return An object of class `"rnaCompare"` containing:
#'   \item{dds}{The DESeqDataSet object.}
#'   \item{res}{A data frame with DESeq2 results.}
#'   \item{contrast}{The contrast used for comparison.}
#'   \item{plots}{A list of plotting functions.}
#'
#' @importFrom DESeq2 DESeq estimateSizeFactors estimateDispersionsGeneEst dispersions "dispersions<-"
#' @importFrom SummarizedExperiment mcols
#'
#' @examples
#' \dontrun{
#' result <- rna.compare(normalized_data)
#' print(result)
#' result$plots$MA()
#' }
#'
#' @export

rna.compare <- function(
    normalized_data,
    group_col = "Group",
    method = "deseq2",
    batch_col = "Batch",
    clean_ids = TRUE,
    assign_result = TRUE,
    assign_name = "comp_data",
    plot = TRUE,
    verbose = TRUE,
    envir = parent.frame()
) {
  # --- Input validation ---
  if (!inherits(normalized_data, "normalized_data"))
    stop("Input must be a 'normalized_data' object.")

  if (!group_col %in% colnames(normalized_data$metadata))
    stop(sprintf("Column '%s' not found in metadata.", group_col))

  if (verbose) message("[rna_compare] Preparing input matrices...")

  counts <- round(as.matrix(normalized_data$expr_matrix))
  metadata <- normalized_data$metadata
  rownames(metadata) <- metadata$Sample
  metadata[[group_col]] <- factor(metadata[[group_col]])

  # --- Optional ID cleaning ---
  if (clean_ids) {
    old_rownames <- rownames(counts)
    new_rownames <- gsub("\\.\\d+$", "", old_rownames)

    if (any(duplicated(new_rownames))) {
      if (verbose)
        message("[rna_compare] Duplicate gene IDs after cleaning; keeping first occurrence.")
      keep <- !duplicated(new_rownames)
      counts <- counts[keep, ]
      rownames(counts) <- new_rownames[keep]
    } else {
      rownames(counts) <- new_rownames
    }

    if (verbose)
      message("[rna_compare] Ensembl IDs cleaned (version suffixes removed).")
  }

  # --- Design formula ---
  design_formula <- if (!is.null(batch_col) && batch_col %in% colnames(metadata)) {
    as.formula(paste("~", batch_col, "+", group_col))
  } else {
    as.formula(paste("~", group_col))
  }

  # --- DESeq2 object creation ---
  if (verbose) message("[rna_compare] Creating DESeq2 dataset...")
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts,
    colData = metadata,
    design = design_formula
  )

  # --- Run DESeq2 with fallback ---
  if (verbose) message("[rna_compare] Running DESeq2...")

  dds_final <- tryCatch({
    DESeq2::DESeq(dds)
  }, error = function(e) {
    if (verbose)
      message("[rna_compare] Parametric fit failed: using gene-wise dispersions.")

    dds <- DESeq2::estimateSizeFactors(dds)
    dds <- DESeq2::estimateDispersionsGeneEst(dds)

    # Set gene-wise dispersions explicitly
    SummarizedExperiment::mcols(dds)
    dispersions(dds) <- mcols(dds)$dispGeneEst

    dds
  })

  # --- Comparison results ---
  groups <- levels(metadata[[group_col]])
  if (length(groups) != 2)
    stop("Currently only two-group comparisons are supported.")

  contrast_use <- c(group_col, groups[2], groups[1])
  res <- DESeq2::results(dds_final, contrast = contrast_use)
  res_df <- as.data.frame(res[order(res$pvalue), ])

  # --- Plot functions ---
  ma_plot <- function() DESeq2::plotMA(
    res,
    ylim = c(-5, 5),
    main = sprintf("MA Plot: %s vs %s", groups[2], groups[1])
  )
  disp_plot <- function() DESeq2::plotDispEsts(dds_final)

  if (plot) {
    suppressMessages(suppressWarnings(ma_plot()))
    suppressMessages(suppressWarnings(disp_plot()))
  }

  # --- Output ---
  out <- list(
    dds = dds_final,
    res = res_df,
    contrast = contrast_use,
    plots = list(MA = ma_plot, Dispersion = disp_plot)
  )
  class(out) <- "rnaCompare"

  if (assign_result) assign(assign_name, out, envir = envir)

  if (verbose)
    message("[rna_compare] Analysis complete. Use `$res` to access results.")

  invisible(out)
}

#' @export
print.rnaCompare <- function(x, ...) {
  cat("Object of class 'rnaCompare'\n",
      "--------------------------------\n",
      sprintf("Comparison: %s vs %s\n", x$contrast[2], x$contrast[3]),
      sprintf("Number of genes: %d\n", nrow(x$res)),
      "--------------------------------\n", sep = "")
  invisible(x)
}
