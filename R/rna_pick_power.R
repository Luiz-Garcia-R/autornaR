#' Diagnose soft-thresholding power for WGCNA
#'
#' @description
#' Evaluates candidate soft-thresholding powers for weighted gene co-expression
#' network analysis (WGCNA) using normalized expression data stored in an
#' \code{rna_project} object. The function filters genes with zero or invalid
#' variance, optionally retains the most variable genes, performs basic
#' expression-data quality control, and evaluates scale-free topology fit
#' across the candidate powers.
#'
#' Diagnostic plots and a summary of the evaluated powers are provided to
#' support the selection of a suitable soft-thresholding power before network
#' construction.
#'
#' @param project An \code{rna_project} object containing normalized
#'   gene expression data and sample metadata.
#'
#' @param genes_n Integer or \code{NULL}. Maximum number of genes to retain,
#'   selected by decreasing expression variance across samples. Genes with
#'   non-finite or zero variance are excluded before selection. If
#'   \code{NULL}, all eligible genes are retained. Using a subset of the most
#'   variable genes can reduce computation time during power evaluation.
#'   Defaults to \code{10000}.
#'
#' @param powers Numeric vector of candidate soft-thresholding powers to
#'   evaluate. Values must be finite, positive, and unique. Defaults to
#'   \code{1:20}.
#'
#' @param network_type Character string specifying the network type used
#'   by WGCNA. Supported values are \code{"signed"}, \code{"unsigned"}, and
#'   \code{"signed hybrid"}. A signed network preserves the distinction
#'   between positive and negative correlations; an unsigned network uses
#'   the magnitude of correlations; and a signed hybrid network excludes
#'   negative correlations from the adjacency calculation. Defaults to
#'   \code{"signed"}.
#'
#' @param fit_threshold Numeric value between 0 and 1 defining the minimum
#'   signed scale-free topology fit index used to flag candidate powers.
#'   The signed fit index is calculated from the scale-free topology fit
#'   and the sign of its regression slope. Powers meeting or exceeding this
#'   threshold are marked in the \code{meets_threshold} column of the
#'   returned diagnostic table. This threshold does not automatically select
#'   a final power or guarantee that the resulting network will be suitable.
#'   Defaults to \code{0.90}.
#'
#' @param plot Logical. If \code{TRUE}, displays diagnostic plots showing
#'   the signed scale-free topology fit across candidate powers and the
#'   corresponding mean network connectivity. The horizontal reference line
#'   in the fit plot corresponds to \code{fit_threshold}. Defaults to
#'   \code{TRUE}.
#'
#' @param save Logical. If \code{TRUE}, attaches the diagnostic result to
#'   the project under the \code{analyses} slot with subtype
#'   \code{"wgcna_power"} and invisibly returns the updated project.
#'   If \code{FALSE}, does not attach the result and instead invisibly
#'   returns an object of class \code{power_result}. Defaults to
#'   \code{TRUE}.
#'
#' @section Input data:
#' The expression matrix must be numeric, with genes in rows and samples
#' in columns. Gene identifiers and sample identifiers must be unique.
#' The metadata must identify samples through a \code{Sample} column or
#' through its row names. Sample identifiers are matched to the expression
#' matrix column names; metadata must contain an entry for every expression
#' sample.
#'
#' Expression values must be finite. Genes with non-finite or zero variance
#' are excluded before the most variable genes are selected. WGCNA quality
#' control is then applied to identify and remove additional problematic
#' samples or genes.
#'
#' @section Interpretation:
#' The scale-free topology fit and mean connectivity plots provide
#' complementary diagnostics. A higher scale-free topology fit may support
#' a candidate power, but increasing the power can also reduce network
#' connectivity. Therefore, power selection should consider both diagnostics
#' and the biological objectives of the analysis.
#'
#' The \code{fit_threshold} argument defines a diagnostic criterion rather
#' than a universal cutoff. If no candidate power reaches the specified
#' threshold, the function reports this result without automatically
#' relaxing the criterion. Alternative thresholds may be explored when
#' justified by the dataset and the intended analysis.
#'
#' @section Scope of the analysis:
#' This function evaluates candidate powers using the samples supplied in
#' the project as a whole. It does not test differential co-expression
#' between experimental groups and does not construct the final network
#' or identify co-expression modules. Group-specific or consensus network
#' analyses require additional steps.
#'
#' @return
#' If \code{save = TRUE}, invisibly returns the updated
#' \code{rna_project} object containing the diagnostic result.
#'
#' If \code{save = FALSE}, invisibly returns an object of class
#' \code{power_result}, containing the fit diagnostics, candidate powers,
#'   network type, fit threshold, numbers of genes and samples used,
#'   identifiers of retained genes and samples, identifiers removed during
#'   quality control, and the analysis timestamp.
#'
#' @seealso
#' \code{\link[WGCNA]{pickSoftThreshold}},
#' \code{\link[WGCNA]{goodSamplesGenes}}
#'
#' @references
#' Langfelder P, Horvath S. (2008).
#' WGCNA: an R package for weighted correlation network analysis.
#' \emph{BMC Bioinformatics}, 9, 559.
#' \doi{10.1186/1471-2105-9-559}
#'
#' @export

rna_pick_power <- function(project,
                           genes_n = 10000,
                           powers = 1:20,
                           network_type = c("signed", "unsigned", "signed hybrid"),
                           fit_threshold = 0.90,
                           plot = TRUE,
                           save = TRUE) {

# =============================================================================
# 0. Basic checks
# =============================================================================

.check_dependencies("WGCNA")

# genes_n
if (!is.null(genes_n)) {
  if (length(genes_n) != 1L ||
      !is.numeric(genes_n) ||
      !is.finite(genes_n) ||
      genes_n < 1 ||
      genes_n != floor(genes_n)) {
    stop("`genes_n` must be a positive integer or NULL.")
  }

  genes_n <- as.integer(genes_n)
}

# powers
if (!is.numeric(powers) ||
    length(powers) == 0 ||
    any(!is.finite(powers)) ||
    any(powers <= 0) ||
    anyDuplicated(powers)) {
  stop("`powers` must contain unique, finite positive values.")
}

# network_type
network_type <- match.arg(network_type)

# fit_threshold
if (!is.numeric(fit_threshold) ||
    length(fit_threshold) != 1L ||
    !is.finite(fit_threshold) ||
    fit_threshold < 0 ||
    fit_threshold > 1) {
  stop("`fit_threshold` must be a finite number between 0 and 1.")
}

# plot
if (!is.logical(plot) || length(plot) != 1L || is.na(plot)) {
  stop("`plot` must be TRUE or FALSE.")
}

# save
if (!is.logical(save) || length(save) != 1L || is.na(save)) {
  stop("`save` must be TRUE or FALSE.")
}

# =============================================================================
# 1. Get active project
# =============================================================================

proj <- project

expr_mat <- as.matrix(.get_expr(proj))
metadata <- .get_meta(proj)

# =============================================================================
# 2. Prepare WGCNA data
# =============================================================================

prepared <- .prepare_wgcna_data(
  expr_mat = expr_mat,
  metadata = metadata,
  genes_n = genes_n
)

datExpr <- prepared$datExpr

# =============================================================================
# 2. Evaluate soft-thresholding powers
# =============================================================================

sft <- WGCNA::pickSoftThreshold(
  datExpr,
  powerVector = powers,
  networkType = network_type,
  verbose = 0
)

fit <- as.data.frame(sft$fitIndices)

if (nrow(fit) == 0) {
  stop("WGCNA did not return soft-thresholding diagnostics.")
}

fit$signed_R2 <- -sign(fit$slope) * fit$SFT.R.sq

fit$meets_threshold <- (
  is.finite(fit$signed_R2) &
    fit$signed_R2 >= fit_threshold
)

# =============================================================================
# 3. Diagnostic plots
# =============================================================================

if (plot) {

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)

  graphics::par(mfrow = c(1, 2))

  # Scale-free topology fit
  graphics::plot(
    fit$Power,
    fit$signed_R2,
    type = "n",
    xlab = "Soft-thresholding power",
    ylab = "Signed scale-free topology fit (R^2)",
    main = "Scale-free topology fit"
  )

  graphics::text(
    fit$Power,
    fit$signed_R2,
    labels = fit$Power,
    col = "red"
  )

  graphics::abline(
    h = fit_threshold,
    col = "red",
    lty = 2
  )

  # Mean connectivity
  graphics::plot(
    fit$Power,
    fit$mean.k.,
    type = "b",
    pch = 19,
    xlab = "Soft-thresholding power",
    ylab = "Mean connectivity",
    main = "Mean network connectivity"
  )
}

# =============================================================================
# 4. Store result
# =============================================================================

obj <- list(
  fit_indices = fit,
  powers = powers,
  network_type = network_type,
  fit_threshold = fit_threshold,
  genes_n = genes_n,

  n_samples = prepared$n_samples,

  genes_used = prepared$genes_used,
  samples_used = prepared$samples_used,

  genes_invalid_variance = prepared$genes_invalid_variance,
  genes_excluded_selection = prepared$genes_excluded_selection,
  removed_genes = prepared$genes_removed_qc,
  removed_samples = prepared$samples_removed_qc,

  timestamp = Sys.time()
)

class(obj) <- "power_result"

# =============================================================================
# 5. Console summary
# =============================================================================

cat("\nWGCNA soft-thresholding diagnostics\n")
cat("-----------------------------------\n")
cat("Network type: ", network_type, "\n", sep = "")
cat("Fit threshold: ", fit_threshold, "\n", sep = "")
cat("Samples used: ", obj$n_samples, "\n", sep = "")
cat("Genes used: ", obj$n_genes, "\n", sep = "")

passing <- fit$Power[fit$meets_threshold]

if (length(passing) > 0) {
  cat(
    "Fit threshold (signed R^2): ", fit_threshold, "\n",
    sep = ""
  )
  cat(
    "Powers meeting the threshold: ",
    paste(passing, collapse = ", "),
    "\n",
    sep = ""
  )
} else {
  cat(
    "No tested power reached signed R^2 >= ",
    fit_threshold,
    ".\n",
    sep = ""
  )
}

print(fit, row.names = FALSE)

# =============================================================================
# 6. Attach to project
# =============================================================================

if (save) {

  proj <- .attach_to_project(
    proj,
    obj,
    slot = "analyses",
    subtype = "wgcna_power",
    log = list(
      network_type = network_type,
      fit_threshold = fit_threshold,
      genes_n = genes_n,
      n_genes = obj$n_genes,
      n_samples = obj$n_samples
    )
  )

}

if (save) {
  return(invisible(proj))
}

return(invisible(obj))
}
