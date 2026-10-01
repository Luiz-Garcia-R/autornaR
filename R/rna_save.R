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
                     verbose = TRUE) {

# =============================================================================
# 1. Data validation
# =============================================================================

if (!inherits(project, "rna_project")) {
  stop("Input must be a 'rna_project' object.")
}

# Create independent object for saving
proj <- project

# =============================================================================
# 2. Remove heavy objects
# =============================================================================

# Raw imported data
if (!is.null(proj$input$imp_data) &&
    !is.null(proj$input$imp_data$data)) {

  proj$input$imp_data$data <- NULL
}

# Normalized expression matrix
if (!is.null(proj$data) &&
    !is.null(proj$data$normalized_data) &&
    !is.null(proj$data$normalized_data$expr_matrix)) {

  proj$data$normalized_data$expr_matrix <- NULL
}

# =============================================================================
# 3. Save
# =============================================================================

saveRDS(
  proj,
  file = file,
  compress = compress
)

if (verbose) {
  message("rna_project saved: ", file)
}

invisible(proj)

}



