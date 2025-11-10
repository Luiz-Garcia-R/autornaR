#' Normalize RNA-seq count data
#'
#' @description
#' Normalizes raw RNA-seq counts from an `imp_data` object created by [rna.import()].
#' Supports multiple normalization methods, low-expression filtering, and optional
#' outlier removal at both sample and gene levels.
#'
#' @param imp_data An object of class `"imp_data"` as returned by [rna.import()].
#' @param method Normalization method. One of:
#'   `"none"`, `"log2"`, `"cpm"`, `"tpm"`, `"rlog"`, `"vst"`, `"quantile"`, or `"upper-quartile"`.
#' @param filter_low Logical; if `TRUE`, removes genes with low expression (default `TRUE`).
#' @param min_expr Numeric; minimum count threshold (default `10`).
#' @param filter_min_prop Numeric; minimum proportion of samples expressing at least
#'   `min_expr` counts (default `0.1`).
#' @param remove_outlier_samples Logical; if `TRUE`, detects and removes outlier samples (default `TRUE`).
#' @param remove_outlier_genes Logical; if `TRUE`, detects and removes outlier genes (default `TRUE`).
#' @param outlier_method Character; method for outlier detection, either `"iqr"` or `"zscore"`.
#' @param assign_result Logical; if `TRUE`, assigns the normalized object in the calling environment (default `TRUE`).
#' @param assign_name Character; name to assign the normalized object (default `"normalized_data"`).
#' @param envir Environment where to assign the object.
#' @param help Logical; if `TRUE`, prints short help information and exits.
#'
#' @details
#' The function performs a sequence of cleaning and normalization steps:
#' 1. Optionally filters out genes expressed in fewer than `filter_min_prop`
#'    of samples at a level below `min_expr`.
#' 2. Detects and optionally removes outlier samples or genes using either
#'    IQR or z-score thresholds.
#' 3. Applies the specified normalization method.
#' 4. Returns a clean, normalized expression matrix ready for downstream analyses.
#'
#' TPM normalization requires a `gene_length` column in `imp_data$annotation`.
#'
#' @return
#' An object of class `"normalized_data"` containing:
#' \describe{
#'   \item{expr_matrix}{Normalized expression matrix.}
#'   \item{metadata}{Sample metadata, realigned with the matrix.}
#'   \item{method}{Normalization method used.}
#'   \item{removed_genes}{Number of genes removed.}
#'   \item{removed_samples}{Number of samples removed.}
#'   \item{QC_metrics}{Basic quality control metrics.}
#' }
#'
#' @importFrom stats quantile sd
#' @export

rna.normalize <- function(imp_data,
                          method = c("none", "log2", "cpm", "tpm", "rlog", "vst",
                                     "quantile", "upper-quartile"),
                          filter_low = TRUE,
                          min_expr = 10,
                          filter_min_prop = 0.1,
                          remove_outlier_samples = TRUE,
                          remove_outlier_genes = TRUE,
                          outlier_method = c("iqr", "zscore"),
                          assign_result = TRUE,
                          assign_name = "normalized_data",
                          envir = parent.frame(),
                          help = FALSE) {

  # --- Help shortcut ---
  if (help || missing(imp_data)) {
    message("rna.normalize(): normalize RNA-seq count data from an 'imp_data' object.\n",
            "Use rna.import() first to create the input object.")
    return(invisible(NULL))
  }

  # --- Validation ---
  if (!inherits(imp_data, "imp_data")) {
    stop("Input must be an object of class 'imp_data'.\n",
         "Use rna.import() first, e.g.: imp <- rna.import(df, metadata, format = 'hisat2')")
  }

  method <- match.arg(method)
  outlier_method <- match.arg(outlier_method)

  data_df <- as.data.frame(imp_data$data)
  gene_col <- grep("^geneid$", colnames(data_df), ignore.case = TRUE, value = TRUE)
  if (length(gene_col) == 0)
    stop("Gene identifier column not found (expected: Geneid, GeneID, or geneID).")

  gene_ids <- data_df[[gene_col]]
  counts <- data_df[, setdiff(colnames(data_df), gene_col), drop = FALSE]
  rownames(counts) <- gene_ids
  metadata <- as.data.frame(imp_data$metadata)

  # --- 1. Filter low-expressed genes ---
  removed_genes <- 0
  if (filter_low) {
    prop_expr <- rowMeans(counts >= min_expr, na.rm = TRUE)
    keep <- prop_expr >= filter_min_prop
    removed_genes <- sum(!keep)
    counts <- counts[keep, , drop = FALSE]
    if (removed_genes > 0) {
      message(removed_genes, " gene(s) removed (min count >= ", min_expr,
              " in >= ", filter_min_prop * 100, "% of samples).")
    }
  }

  # --- 2. Outlier detection ---
  outliers_s <- integer(0)
  if (remove_outlier_samples) {
    lib_sizes <- colSums(counts)
    if (outlier_method == "iqr") {
      q <- quantile(lib_sizes, c(0.25, 0.75))
      iqr <- diff(q)
      outliers_s <- which(lib_sizes < (q[1] - 1.5 * iqr) | lib_sizes > (q[2] + 1.5 * iqr))
    } else {
      z <- (lib_sizes - mean(lib_sizes)) / sd(lib_sizes)
      outliers_s <- which(abs(z) > 3)
    }
    if (length(outliers_s) > 0) {
      removed_samples <- colnames(counts)[outliers_s]
      counts <- counts[, -outliers_s, drop = FALSE]
      metadata <- metadata[metadata$Sample %in% colnames(counts), , drop = FALSE]
      message(length(outliers_s), " outlier sample(s) removed: ",
              paste(removed_samples, collapse = ", "))
    }
  }

  if (remove_outlier_genes && nrow(counts) > 2) {
    gene_sds <- apply(counts, 1, sd)
    if (outlier_method == "iqr") {
      q <- quantile(gene_sds, c(0.25, 0.75))
      iqr <- diff(q)
      outliers_g <- which(gene_sds < (q[1] - 1.5 * iqr) | gene_sds > (q[2] + 1.5 * iqr))
    } else {
      z <- (gene_sds - mean(gene_sds)) / sd(gene_sds)
      outliers_g <- which(abs(z) > 3)
    }
    if (length(outliers_g) > 0) {
      counts <- counts[-outliers_g, , drop = FALSE]
      removed_genes <- removed_genes + length(outliers_g)
      message(length(outliers_g), " outlier gene(s) removed.")
    }
  }

  # --- 3. Normalization ---
  norm_counts <- counts
  if (method == "log2") {
    norm_counts <- log2(norm_counts + 1)
  } else if (method == "cpm") {
    lib_sizes <- colSums(counts)
    norm_counts <- t(t(counts) / lib_sizes * 1e6)
  } else if (method == "tpm") {
    if (is.null(imp_data$annotation) ||
        !"gene_length" %in% colnames(imp_data$annotation)) {
      stop("TPM normalization requires a 'gene_length' column in imp_data$annotation.")
    }
    rpk <- counts / (imp_data$annotation$gene_length / 1000)
    norm_counts <- t(t(rpk) / colSums(rpk) * 1e6)
  } else if (method == "quantile") {
    if (!requireNamespace("limma", quietly = TRUE))
      stop("Please install 'limma' to use quantile normalization.")
    norm_counts <- limma::normalizeQuantiles(as.matrix(counts))
  } else if (method == "upper-quartile") {
    uq <- apply(counts, 2, function(x) quantile(x[x > 0], 0.75))
    norm_counts <- t(t(counts) / uq * median(uq))
  } else if (method %in% c("rlog", "vst")) {
    if (!requireNamespace("DESeq2", quietly = TRUE))
      stop("DESeq2 not installed. Install with: BiocManager::install('DESeq2')")

    counts <- round(counts)

    # Create object DESeq2
    dds <- DESeq2::DESeqDataSetFromMatrix(
      countData = counts,
      colData = metadata,
      design = ~ 1
    )

    if (method == "vst") {
      norm_counts <- DESeq2::vst(dds, blind = TRUE)
      norm_counts <- SummarizedExperiment::assay(norm_counts)
    } else {
      norm_counts <- DESeq2::rlog(dds, blind = TRUE)
      norm_counts <- SummarizedExperiment::assay(norm_counts)
    }
  }

  # --- 4. Quality metrics ---
  QC_metrics <- data.frame(
    Sample = colnames(norm_counts),
    Library_size = colSums(norm_counts),
    Genes_detected = colSums(norm_counts > 0),
    row.names = NULL
  )

  # --- 5. Metadata alignment ---
  metadata <- metadata[match(colnames(norm_counts), metadata$Sample), , drop = FALSE]
  if (!all(metadata$Sample == colnames(norm_counts))) {
    warning("Metadata and expression matrix were realigned to match sample order.")
  }

  # --- 6. Final object ---
  obj <- list(
    expr_matrix = norm_counts,
    metadata = metadata,
    method = method,
    removed_genes = removed_genes,
    removed_samples = length(outliers_s),
    QC_metrics = QC_metrics
  )
  class(obj) <- "normalized_data"

  if (assign_result) {
    assign(assign_name, obj, envir = envir)
    print(obj)
    invisible(obj)
  } else {
    return(obj)
  }
}

#' @export
print.normalized_data <- function(x, ...) {
  message("================================")
  message("Object of class 'normalized_data'")
  message("================================")
  message("Samples: ", ncol(x$expr_matrix))
  message("Genes:   ", nrow(x$expr_matrix))
  message("Method:  ", x$method)
  message("Removed genes: ", x$removed_genes)
  message("Removed samples: ", x$removed_samples)
  message("================================")
  invisible(x)
}
