#' Initialize an RNA-seq analysis project container
#'
#' @description
#' Creates a structured \code{rna_project} object that serves as the central
#' container for all RNA-seq analyses within the \pkg{autornaR} workflow.
#'
#' The project object is designed to store raw inputs, processed data,
#' analysis results, and execution logs in a consistent and reproducible
#' structure, enabling seamless transitions between analysis steps.
#'
#' Once created, the project is automatically assigned to the global
#' environment as \code{rna_project}, becoming the active context for all
#' subsequent functions (e.g., \code{rna.import()}, \code{rna.normalize()},
#' \code{rna.compare()}).
#'
#' @param project_name Character. Optional project identifier used for
#'   labeling outputs and logs.
#'
#' @return
#' Invisibly returns an object of class \code{"rna_project"}, containing:
#' \describe{
#'   \item{project_info}{Metadata about the project (name, creation time, preprocessing steps, list of downstream analyses).}
#'   \item{input}{Raw imported data (e.g., count matrices, annotations).}
#'   \item{data}{Processed data objects (e.g., normalized expression matrices).}
#'   \item{logs}{Execution logs for reproducibility and tracking.}
#'   \item{version}{Package version used to create the project.}
#' }
#'
#' @details
#' The \code{rna_project} object acts as a lightweight state manager for the
#' analysis workflow. Instead of passing data manually between functions,
#' each step reads from and writes to this centralized structure.
#'
#' \strong{Design principles:}
#' \itemize{
#'   \item \emph{Reproducibility}: All steps are stored and traceable via logs.
#'   \item \emph{Modularity}: Each function updates a specific part of the project.
#'   \item \emph{State awareness}: The current analysis state is always explicit.
#'   \item \emph{Safety}: Prevents loss of intermediate results across steps.
#' }
#'
#' \strong{Typical workflow:}
#' \enumerate{
#'   \item Initialize project with \code{rna.project()}
#'   \item Import data using \code{rna.import()}
#'   \item Normalize expression with \code{rna.normalize()}
#'   \item Perform analyses (e.g., \code{rna.compare()}, \code{rna.gsea()})
#' }
#'
#' Calling \code{rna.project()} again will overwrite the current active project.
#'
#' @examples
#' \dontrun{
#'
#' # Initialize a new project
#' my_project <- rna.project(project_name = "my_project")
#'
#' # Continue workflow updating project
#' my_project <- rna.import(...)
#' my_project <- rna.normalize(...)
#' my_project <- rna.compare(...)
#'
#' }
#'
#' @seealso
#' \code{\link{rna.import}},
#' \code{\link{rna.normalize}},
#' \code{\link{rna.qc}},
#' \code{\link{rna.compare}}
#'
#' @export

rna.project <- function(project_name = NULL) {

  obj <- list(
    project_info = list(
      name = if (is.null(project_name)) "rna_project" else project_name,
      created = Sys.time(),
      preprocessing = NULL,
      analyses = list()
    ),
    input = list(
      imp_data = NULL
    ),
    data = NULL,
    logs = list(),
    version = as.character(utils::packageVersion("autornaR"))
  )

  class(obj) <- "rna_project"

  message("New RNA project initialized.")

  return(invisible(obj))
}
