#' Save a lightweight RNA project for meta-analysis
#'
#' @description
#' Extracts the essential components from a \code{rna_project} object and saves
#' a simplified, memory-efficient version as an \code{.rds} file. The function
#' removes large intermediate objects (e.g., expression matrices, model objects,
#' and plots) while preserving key results required for downstream analyses and
#' cross-study comparisons.
#'
#' This function is designed as a bridge between individual RNA-seq analyses
#' and meta-analysis workflows (e.g., integration via \code{metaOmicsR}).
#'
#' @param project A \code{rna_project} object generated through the standard
#'   workflow (\code{rna.project()}, \code{rna.import()}, \code{rna.normalize()},
#'   \code{rna.compare()}, etc.).
#' @param file Character. Output file path where the \code{.rds} object will be saved.
#' @param compress Character. Compression method passed to \code{saveRDS()}
#'   (default: \code{"xz"}).
#' @param verbose Logical. If \code{TRUE}, prints a message upon successful save
#'   (default: \code{TRUE}).
#'
#' @return
#' Invisibly returns a cleaned \code{rna_project} object containing:
#' \describe{
#'   \item{study_info}{General study metadata (organism, platform, sample size,
#'   gene identifier type, and contrast definition).}
#'   \item{gene_info}{Gene annotation table (gene_id, symbol, entrez).}
#'   \item{de_results}{Differential expression results with standardized column names.}
#'   \item{enrichment}{Pathway enrichment results (e.g., GSEA output).}
#'   \item{qc}{Quality control summary extracted from project logs.}
#'   \item{version}{Version of the originating \code{rna_project}.}
#' }
#'
#' @details
#' The function performs the following steps:
#'
#' \strong{1. Validation} \cr
#' Ensures the input object is a valid \code{rna_project}.
#'
#' \strong{2. Extraction of latest results} \cr
#' For each component (input, normalization, comparison, enrichment, QC),
#' only the most recent analysis (identified by the \code{"last"} pointer)
#' is retained to avoid redundancy and reduce file size.
#'
#' \strong{3. Contrast parsing} \cr
#' The comparison contrast is parsed to extract:
#' \itemize{
#'   \item \code{condition_tested} (e.g., treated, disease)
#'   \item \code{reference_condition} (e.g., control, baseline)
#' }
#' This information is critical for consistent interpretation of effect
#' directions across studies in meta-analysis.
#'
#' \strong{4. Object simplification} \cr
#' Heavy elements such as expression matrices, model objects, and intermediate
#' structures are removed. Only summary-level data required for integration
#' are preserved.
#'
#' \strong{5. Serialization} \cr
#' The cleaned object is saved as a compressed \code{.rds} file, enabling
#' efficient storage and fast reloading in R environments.
#'
#' @examples
#' \dontrun{
#' # Save a processed project for meta-analysis
#' rna.save(project = my_project,
#'          file = "my_project.rds")
#'
#' # Save with alternative compression
#' rna.save(my_project,
#'          file = "my_project.rds",
#'          compress = "gzip")
#' }
#'
#' @seealso
#' \code{\link{rna.project}}, \code{\link{rna.import}},
#' \code{\link{rna.compare}}, \code{\link{rna.gsea}}
#'
#' @export

rna.save <- function(project,
                     file,
                     compress = "xz",
                     verbose = TRUE
) {

  # ===========================================================================
  # 1) Data validation
  # ===========================================================================

  if (!inherits(project, "rna_project")) {
    stop("Input must be a 'rna_project' object")
  }

  proj <- project

  get_last <- function(x) {
    if (is.null(x)) return(NULL)

    if (is.list(x) && "last" %in% names(x)) {
      id <- x$last

      if (!id %in% names(x)) {
        warning("Invalid 'last' reference detected. Returning NULL.")
        return(NULL)
      }

      return(x[[id]])
    }

    return(x)
  }

  # ===========================================================================
  # 2) Extract contrast info
  # ===========================================================================
  comp <- get_last(proj$analyses$comparison)

  condition_tested <- NA
  reference_condition <- NA

  if (!is.null(comp$contrast)) {

    contrast_raw <- comp$contrast

    contrast_clean <- gsub("^Group", "", contrast_raw)

    split_pattern <- if (grepl(" - ", contrast_clean)) {
      " - "
    } else if (grepl("_vs_", contrast_clean)) {
      "_vs_"
    } else {
      NULL
    }

    if (!is.null(split_pattern)) {
      parts <- strsplit(contrast_clean, split_pattern)[[1]]

      if (length(parts) == 2) {
        condition_tested    <- parts[1]
        reference_condition <- parts[2]
      }
    }
  }

  # ===========================================================================
  # 3) study_info
  # ===========================================================================
  imp <- get_last(proj$input$imp_data)

  study_info <- list(
    project_name = proj$project_info$name,
    organism     = imp$organism,
    n_genes      = imp$n_genes,
    n_samples    = imp$n_samples,
    platform     = imp$detected_format,
    gene_id_type = imp$gene_id_type,
    condition_tested    = condition_tested,
    reference_condition = reference_condition
  )

  # ===========================================================================
  # 4) gene_info
  # ===========================================================================
  de_results <- comp$res

  gene_info <- imp$gene_annotation

  # ===========================================================================
  # 5) DE results
  # ===========================================================================
  de_results$gene_id <- rownames(de_results)
  rownames(de_results) <- NULL

  colnames(de_results)[colnames(de_results) == "log2FoldChange"] <- "logFC"
  colnames(de_results)[colnames(de_results) == "stat"] <- "statistic"

  # ===========================================================================
  # 6) GSEA (enrichment)
  # ===========================================================================
  gsea <- get_last(proj$analyses$gsea)

  enrichment <- NULL
  if (!is.null(gsea)) {
    enrichment <- gsea$gsea_full
  }

  enrichment$leadingEdge <- NULL
  colnames(enrichment)[colnames(enrichment) == "pval"] <- "pvalue"

  # ===========================================================================
  # 7) QC - Logs
  # ===========================================================================

  qc <- get_last(proj$logs$QC)

  # ===========================================================================
  # 8) Build clean object
  # ===========================================================================

  clean <- list(
    study_info = study_info,
    gene_info  = gene_info,
    de_results = de_results,
    enrichment = enrichment,
    qc         = qc,
    version    = proj$version
  )

  class(clean) <- "rna_project"

  # ===========================================================================
  # 9) Save
  # ===========================================================================

  saveRDS(clean, file = file, compress = compress)

  if (verbose) {
    message("rna_project saved: ", file)
  }

  invisible(clean)
}
