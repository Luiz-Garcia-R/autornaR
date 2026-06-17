#' Normalize RNA-seq count data with integrated QC and filtering
#'
#' @description
#' Performs normalization of raw RNA-seq count data stored in the active
#' \code{rna_project} object (created with \code{rna.project()} and populated
#' via \code{rna.import()}). The function integrates gene filtering,
#' outlier detection, and multiple normalization strategies into a single,
#' reproducible workflow.
#'
#' The resulting normalized dataset is stored internally and used as the
#' standardized input for downstream analyses (e.g., differential expression,
#' clustering, and pathway analysis).
#'
#' @param project \code{rna_project} object created by \code{rna.project()}.
#' @param method Character. Normalization method to apply:
#' \itemize{
#'   \item \code{"none"}: No normalization (raw counts retained)
#'   \item \code{"log2"}: Log2 transformation (\code{log2(count + 1)})
#'   \item \code{"cpm"}: Counts per million (library size normalization)
#'   \item \code{"tpm"}: Transcripts per million (length-normalized expression)
#'   \item \code{"rlog"}: Regularized log transformation (DESeq2)
#'   \item \code{"vst"}: Variance stabilizing transformation (DESeq2)
#'   \item \code{"quantile"}: Quantile normalization (limma)
#'   \item \code{"upper-quartile"}: Upper-quartile normalization
#' }
#'
#' @param filter_low Logical. If \code{TRUE}, removes lowly expressed genes
#'   prior to normalization (default: \code{TRUE}).
#' @param min_expr Numeric. Minimum count threshold used to define expression
#'   (default: \code{10}).
#' @param filter_min_prop Numeric. Minimum proportion of samples in which a gene
#'   must be expressed (above \code{min_expr}) to be retained (default: \code{0.1}).
#' @param remove_outlier_samples Logical. If \code{TRUE}, detects and removes
#'   samples with abnormal library sizes (default: \code{TRUE}).
#' @param remove_outlier_genes Logical. If \code{TRUE}, removes genes with extreme
#'   variance profiles (default: \code{TRUE}).
#' @param outlier_method Character. Method for outlier detection:
#'   \code{"iqr"} (robust, default) or \code{"zscore"}.
#' @param clean_gene_versions Logical. If \code{TRUE}, removes version suffixes
#'   from gene identifiers (e.g., \code{ENSG000001.5 -> ENSG000001}) when using
#'   Ensembl IDs (default: \code{TRUE}).
#' @param gene_id_sep Character. Regular expression defining the separator used
#'   for gene version removal (default: \code{"\\."}).
#' @param save Logical. Whether to store the normalized data inside the active
#'   \code{rna_project} (default: \code{TRUE}).
#'
#' @return
#' Invisibly returns an object of class \code{"normalized_data"} containing:
#' \describe{
#'   \item{expr_matrix}{Normalized expression matrix (genes × samples).}
#'   \item{metadata}{Aligned sample metadata.}
#'   \item{method}{Normalization method used.}
#'   \item{removed_genes}{Number of genes removed during filtering/outlier steps.}
#'   \item{removed_samples}{Number of samples removed as outliers.}
#'   \item{QC_metrics}{Data frame with per-sample QC metrics (library size and detected genes).}
#'   \item{gene_id_version_cleaned}{Logical indicating whether gene IDs were cleaned.}
#'   \item{gene_id_type}{Original gene identifier type (e.g., ENSEMBL, SYMBOL).}
#' }
#'
#' @details
#' The normalization workflow is structured into sequential steps:
#'
#' \strong{1. Gene filtering} \cr
#' Genes with insufficient expression are removed based on the proportion
#' of samples exceeding \code{min_expr}. This reduces noise and improves
#' statistical stability in downstream analyses.
#'
#' \strong{2. Outlier detection} \cr
#' Outlier samples are identified based on library size distributions, while
#' outlier genes are detected using variance-based criteria. Two strategies
#' are available:
#' \itemize{
#'   \item \code{"iqr"}: Robust, recommended for most datasets
#'   \item \code{"zscore"}: Sensitive to extreme deviations
#' }
#'
#' \strong{3. Normalization} \cr
#' The selected normalization method is applied. Methods differ in their
#' assumptions:
#' \itemize{
#'   \item Count-based scaling (\code{cpm}, \code{upper-quartile})
#'   \item Length-aware normalization (\code{tpm})
#'   \item Variance stabilization (\code{vst}, \code{rlog})
#'   \item Distribution alignment (\code{quantile})
#' }
#'
#' \strong{4. Quality control metrics} \cr
#' Basic QC metrics are computed, including library size and number of
#' detected genes per sample, enabling rapid inspection of data quality.
#'
#' \strong{5. Gene ID harmonization} \cr
#' When Ensembl IDs include version suffixes, they are optionally removed
#' to ensure compatibility with annotation databases and downstream tools.
#' If duplicates arise after cleaning, counts are aggregated by gene.
#'
#' \strong{TPM requirement} \cr
#' TPM normalization requires a \code{gene_length} column in
#' \code{rna_project$input$imp_data$annotation}.
#'
#' @section Filtering strategies and use cases:
#'
#' Filtering parameters can strongly influence downstream analyses.
#' The default configuration is optimized for general RNA-seq workflows,
#' particularly differential expression analysis.
#'
#' However, more permissive filtering may be required when the goal is to
#' preserve genes expressed in rare or specific cell populations.
#'
#' \strong{Standard analysis (default)}
#' \itemize{
#'   \item Removes lowly expressed genes to reduce noise
#'   \item Recommended for differential expression and global analyses
#'   \item Default parameters:
#'     \itemize{
#'       \item \code{min_expr = 10}
#'       \item \code{filter_min_prop = 0.1}
#'       \item \code{remove_outlier_genes = TRUE}
#'     }
#' }
#'
#' \strong{Deconvolution / gene set scoring}
#' \itemize{
#'   \item Preserves low-abundance and cell-type-specific genes
#'   \item Recommended for immune profiling, cell signatures, and scoring methods
#'   \item Suggested parameters:
#'     \itemize{
#'       \item \code{min_expr = 1}
#'       \item \code{filter_min_prop = 0.05}
#'       \item \code{remove_outlier_genes = FALSE}
#'     }
#' }
#'
#' \strong{Important:}
#' More stringent filtering may remove biologically relevant markers,
#' especially for rare cell populations.
#'
#' @seealso
#' \code{\link{rna.import}}, \code{\link{rna.project}},  \code{\link{rna.qc}}
#'
#' @examples
#' \dontrun{
#'
#' Example 1: basic normalization using default settings (suitable for general analysis)
#' rna.normalize(project = my_project)
#'
#' # Example 2: log2 transformation only (quick exploratory analysis)
#' rna.normalize(my_project,
#'               method = "log2")
#'
#' # Example 3: Using relaxed filtering thresholds for deconvolution and gene set scoring
#' proj <- rna.normalize(my_project,
#'                       method = "vst",   # or "log2"
#'                       filter_low = TRUE,
#'                       min_expr = 1,
#'                       filter_min_prop = 0.05,
#'                       remove_outlier_genes = FALSE
#' )
#' }
#'
#'
#'
#'
#' @importFrom stats quantile sd
#'
#' @export

rna.normalize <- function(project,
                          method = c("none", "log2", "cpm", "tpm", "rlog", "vst", "quantile", "upper-quartile"),
                          filter_low = TRUE,
                          min_expr = 10,
                          filter_min_prop = 0.1,
                          remove_outlier_samples = TRUE,
                          remove_outlier_genes = TRUE,
                          outlier_method = c("iqr", "zscore"),
                          clean_gene_versions = TRUE,
                          gene_id_sep = "\\.",
                          save = TRUE
) {

  # ---------------------------
  # 1) Get active project
  # ---------------------------
  proj <- project

  # ---------------------------
  # 2) Validate input
  # ---------------------------
  if (is.null(proj$input$imp_data)) {
    stop("No imported data found. Run rna.import() first.")
  }

  # --- Access last imp_data ---
  imp_container <- proj$input$imp_data

  if (is.null(imp_container$last)) {
    stop("No active import found in 'imp_data$last'")
  }

  imp_data <- imp_container[[imp_container$last]]

  method <- match.arg(method)
  outlier_method <- match.arg(outlier_method)

  data_df <- as.data.frame(imp_data$data)
  gene_col <- grep("^geneid$", colnames(data_df), ignore.case = TRUE, value = TRUE)
  if (length(gene_col) == 0)
    stop("Gene identifier column not found (expected: Geneid, GeneID, or geneID).")

  gene_ids <- data_df[[gene_col]]
  counts <- data_df[, setdiff(colnames(data_df), gene_col), drop = FALSE]
  counts <- data.frame(
    lapply(counts, function(col) as.numeric(col)),
    check.names = FALSE
  )
  rownames(counts) <- gene_ids

  if (!is.null(proj$data) && !is.null(proj$data$normalized_data)) {
    warning("Overwriting existing normalized data.")
  }

  # --- Gene ID version cleanup ---
  if (clean_gene_versions && imp_data$gene_id_type == "ENSEMBL") {
    gene_ids_clean <- sub(paste0(gene_id_sep, ".*$"), "", rownames(counts))

    if (anyDuplicated(gene_ids_clean)) {
      message("Gene ID version removal produced duplicated IDs. Aggregating counts by gene.")

      # --- Aggregation by sum ---
      counts <- rowsum(counts, group = gene_ids_clean)
    } else {
      rownames(counts) <- gene_ids_clean
    }
  }

  metadata <- as.data.frame(imp_data$metadata)
  metadata <- metadata[match(colnames(counts), metadata$Sample), , drop = FALSE]

  # ---------------------------
  # 3) Filter low-expressed genes
  # ---------------------------
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

  # ---------------------------
  # 4) Outlier detection
  # ---------------------------
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

  # ---------------------------
  # 5) Normalization
  # ---------------------------
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

    ann <- imp_data$annotation
    ann <- ann[match(rownames(counts), ann$gene_id), ]
    if (any(is.na(ann$gene_length))) {
      stop("Gene length missing for some genes.")
    }
    rpk <- counts / (ann$gene_length / 1000)
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

    # --- Create object DESeq2 ---
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

  # ---------------------------
  # 6) Quality metrics
  # ---------------------------
  QC_metrics <- data.frame(
    Sample = colnames(norm_counts),
    Library_size = colSums(norm_counts),
    Genes_detected = colSums(norm_counts > 0),
    row.names = NULL
  )

  # ---------------------------
  # 7) Metadata alignment
  # ---------------------------
  metadata <- metadata[match(colnames(norm_counts), metadata$Sample), , drop = FALSE]
  if (!all(metadata$Sample == colnames(norm_counts))) {
    warning("Metadata and expression matrix were realigned to match sample order.")
  }

  # ---------------------------
  # 8) Final object
  # ---------------------------
  obj <- list(
    timestamp = Sys.time(),
    expr_matrix = norm_counts,
    metadata = metadata,
    method = method,
    removed_genes = removed_genes,
    removed_samples = length(outliers_s),
    QC_metrics = QC_metrics,
    gene_id_version_cleaned = clean_gene_versions,
    gene_id_type = imp_data$gene_id_type
  )
  class(obj) <- "normalized_data"

  # ---------------------------
  # 9) Attach to project
  # ---------------------------
  if (save) {

  proj <- .attach_to_project(
    proj,
    obj,
    slot = "data",
    subtype = "normalized_data",
    prefix = "norm",
    log = list(
      method = method,
      n_genes = nrow(obj$expr_matrix),
      n_samples = ncol(obj$expr_matrix),
      removed_genes = removed_genes,
      removed_samples = length(outliers_s)
    )
  )
}

  # ---------------------------
  # 10) Return
  # ---------------------------
  .print_header("RNA Normalization")

  .print_block("Summary", function() {
    cat("Method:            ", method, "\n")
    cat("Samples:           ", ncol(obj$expr_matrix), "\n")
    cat("Genes:             ", nrow(obj$expr_matrix), "\n")
    cat("Removed genes:     ", removed_genes, "\n")
    cat("Removed samples:   ", length(outliers_s), "\n")
  })

  .print_block("QC Metrics (first 5 samples)", function() {
    print(utils::head(obj$QC_metrics, 5))
  })

  return(invisible(proj))
}
