# =============================================================================
# RNA GSVA
# =============================================================================

#' Gene Set Variation Analysis (GSVA)
#'
#' @description
#'   Estimates pathway activity scores for individual samples using gene
#'   expression data and predefined gene sets. Gene sets can be obtained from
#'   a previous \code{rna_gsea()} analysis, MSigDB collections, or
#'   user-supplied pathways.
#'
#'   Unlike GSEA, which evaluates pathway enrichment between experimental
#'   conditions using a ranked gene list, GSVA-based methods estimate
#'   pathway-level variation across individual samples, allowing downstream
#'   comparisons and longitudinal analyses.
#'
#' @param project An \code{rna_project} object containing normalized gene
#'   expression data.
#'
#' @param source Character specifying the source of gene sets:
#'   \code{"gsea"} to reuse pathways from a previous GSEA analysis,
#'   \code{"msigdb"} to retrieve gene sets from MSigDB, or
#'   \code{"custom"} to use user-supplied pathways.
#'
#' @param gsea Optional object of class \code{gsea_result}, required when
#'   \code{source = "gsea"}.
#'
#' @param pathways Optional named list of gene sets, required when
#'   \code{source = "custom"}. Each element should contain gene identifiers
#'   corresponding to the expression matrix.
#'
#' @param pathway_scope Character specifying which pathways to retrieve from
#'   a previous GSEA result: \code{"top"}, \code{"significant"}, or
#'   \code{"all"}. Only applicable when \code{source = "gsea"}.
#'
#' @param geneset_collection Character specifying the MSigDB collection
#'   when \code{source = "msigdb"} (e.g., \code{"H"}, \code{"C2"},
#'   \code{"C5"}).
#'
#' @param geneset_subcollection Optional character specifying an MSigDB
#'   subcollection (e.g., \code{"CP:REACTOME"}, \code{"GO:BP"}).
#'
#' @param method Character specifying the pathway scoring method:
#'   \code{"gsva"}, \code{"ssgsea"}, \code{"zscore"}, or \code{"plage"}.
#'
#' @param kcdf Character specifying the kernel distribution used by GSVA
#'   to estimate expression distributions.
#'
#' @param min_size Integer specifying the minimum number of genes required
#'   for a gene set to be included in the analysis.
#'
#' @param max_size Integer specifying the maximum number of genes allowed
#'   in a gene set.
#'
#' @param save Logical. If \code{TRUE}, stores the GSVA result in the
#'   \code{rna_project} object.
#'
#' @section Gene set sources:
#'   Gene sets can be obtained from three different sources:
#'   \itemize{
#'     \item \strong{GSEA:} reuses pathways from a previous \code{rna_gsea()}
#'       analysis. The \code{pathway_scope} argument determines whether the
#'       top-ranked, statistically significant, or all available pathways are
#'       considered.
#'
#'     \item \strong{MSigDB:} retrieves gene sets directly from Molecular
#'       Signatures Database collections, allowing pathway activity estimation
#'       independently of a previous GSEA analysis.
#'
#'     \item \strong{Custom:} accepts user-defined gene sets, allowing the
#'       evaluation of specific biological signatures or experimental gene
#'       lists.
#'   }
#'
#' @section Scoring methods:
#'   Four methods are available for estimating pathway activity:
#'   \itemize{
#'     \item \strong{GSVA:} estimates variation in pathway activity across
#'       samples using a non-parametric approach based on expression
#'       distributions.
#'
#'     \item \strong{ssGSEA:} calculates enrichment scores for each gene set
#'       within individual samples based on gene expression rankings.
#'
#'     \item \strong{z-score:} combines standardized expression values of genes
#'       within each gene set to estimate pathway-level scores.
#'
#'     \item \strong{PLAGE:} uses singular value decomposition (SVD) to derive
#'       pathway activity scores from the expression profiles of genes within
#'       each set.
#'   }
#'
#'   These methods differ in their scoring assumptions and numerical properties.
#'   Scores obtained using different methods should not be assumed to be
#'   directly comparable.
#'
#' @section Interpretation:
#'   The resulting scores represent relative pathway-level expression patterns
#'   across samples and can be used in downstream analyses, including
#'   differential comparisons, clustering, correlation analysis, and
#'   longitudinal modeling.
#'
#'   Higher scores generally indicate stronger representation of the
#'   corresponding gene expression signature, but do not necessarily demonstrate
#'   increased biochemical activity or functional activation of the pathway.
#'   For PLAGE, score direction is arbitrary and should be interpreted with
#'   respect to the underlying expression patterns.
#'
#'   Results depend on the selected gene sets, gene identifier compatibility,
#'   expression data, and scoring method. Gene sets containing too few or too
#'   many genes after matching to the expression matrix are excluded according
#'   to \code{min_size} and \code{max_size}.
#'
#' @section Pipeline context:
#'   GSVA is typically performed after expression normalization and before
#'   downstream statistical analyses. Unlike GSEA, it does not require a
#'   previous differential expression comparison and can therefore be used
#'   to investigate pathway-level variation across all samples in a project.
#'
#' @return
#'   An object of class \code{gsva_result} containing the pathway activity
#'   scores and associated analysis information.
#'
#' @seealso
#'   \code{\link{rna_gsea}},
#'   \code{\link{rna_compare}},
#'   \code{\link[GSVA]{gsva}}
#'
#' @export



rna_gsva <- function(
    project,
    source = c("gsea", "msigdb", "custom"),
    pathways = NULL,
    pathway_scope = c("top", "significant", "all"),
    geneset_collection = "H",
    geneset_subcollection = NULL,
    method = c("gsva", "ssgsea", "zscore", "plage"),
    kcdf = "Gaussian",
    min_size = 10,
    max_size = 500,
    save = TRUE
) {

source <- match.arg(source)
pathway_scope <- match.arg(pathway_scope)
method <- match.arg(method)

pathway_metadata <- NULL

# =============================================================================
# 0. Basic checks
# =============================================================================

.check_dependencies("GSVA", bioc = TRUE)

# =============================================================================
# 1. Get active project
# =============================================================================

proj <- project

expr <- as.matrix(.get_expr(proj))
organism <- .get_organism(proj)
gsea_obj <- .get_gsea(proj)

# Harmonize gene identifiers
expr <- .convert_expr_to_symbols(
  expr = expr,
  proj = proj
)

# =============================================================================
# 2. Obtain pathways
# =============================================================================

if (source == "gsea") {

  if (is.null(gsea_obj)) {
    stop(
      "No GSEA result found in project. Run rna_gsea() first."
    )
  }

  if (pathway_scope == "all") {

    pathways <- gsea_obj$pathways
    pathway_metadata <- gsea_obj$gsea_full

  } else if (pathway_scope == "significant") {

    sig_names <- gsea_obj$gsea_full$pathway[
      gsea_obj$gsea_full$padj <
        gsea_obj$params$padj_gsea
    ]

    pathways <- gsea_obj$pathways[sig_names]

    pathway_metadata <- gsea_obj$gsea_full[
      gsea_obj$gsea_full$padj <
        gsea_obj$params$padj_gsea,
    ]

  } else if (pathway_scope == "top") {

    top_names <- gsea_obj$gsea_top$pathway
    pathways <- gsea_obj$pathways[top_names]
    pathway_metadata <- gsea_obj$gsea_top
  }
}

# ===============================================
# 2.1. msigdb
# ===============================================

if (source == "msigdb") {

  .check_dependencies("msigdbr")

  species_name <- switch(
    organism,
    human = "Homo sapiens",
    mouse = "Mus musculus",
    zebrafish = "Danio rerio"
  )

  if (is.null(geneset_subcollection)) {

    msig <- msigdbr::msigdbr(
      species = species_name,
      collection = geneset_collection
    )

  } else {

    msig <- msigdbr::msigdbr(
      species = species_name,
      collection = geneset_collection,
      subcollection = geneset_subcollection
    )

  }

  pathways <- split(
    msig$gene_symbol,
    msig$gs_name
  )

  pathway_metadata <- data.frame(
    pathway = names(pathways),
    source = "MSigDB"
  )
}

# ===============================================
# 2.2. Custom
# ===============================================

if (source == "custom") {

  if (is.null(pathways)) {
    stop(
      "Provide 'pathways' when source = 'custom'."
    )
  }
}

# =============================================================================
# 3. Filter pathways
# =============================================================================

pathways <- pathways[
  lengths(pathways) >= min_size &
    lengths(pathways) <= max_size
]

if (length(pathways) == 0) {
  stop("No pathways remaining after filtering.")
}

# =============================================================================
# 4. Run GSVA
# =============================================================================

# gsva
if (method == "gsva") {

  param <- GSVA::gsvaParam(
    exprData = expr,
    geneSets = pathways,
    kcdf = kcdf
  )

  # ssgsea
} else if (method == "ssgsea") {

  param <- GSVA::ssgseaParam(
    exprData = expr,
    geneSets = pathways
  )

  # z-score
} else if (method == "zscore") {

  param <- GSVA::zscoreParam(
    exprData = expr,
    geneSets = pathways
  )

  # plage
} else if (method == "plage") {

  param <- GSVA::plageParam(
    exprData = expr,
    geneSets = pathways
  )

}

message(
  "[rna_gsva] ",
  length(pathways),
  " pathways retained after size filtering."
)

score_matrix <- GSVA::gsva(param)

# =============================================================================
# 5. Return object
# =============================================================================

# Parameters
params <- list(
  timestamp = Sys.time(),
  source = source,
  method = method,
  organism = organism,
  geneset_collection = geneset_collection,
  geneset_subcollection = geneset_subcollection,
  pathway_scope = pathway_scope,
  min_size = min_size,
  max_size = max_size
)

# Final object
obj <- list(
  params = params,
  pathway_scores = score_matrix,
  pathways = pathways,
  pathway_metadata = pathway_metadata,
  n_pathways = nrow(score_matrix),
  n_samples = ncol(score_matrix)
)

class(obj) <- "gsva_result"

# =============================================================================
# 6. Attach to project
# =============================================================================

if (save) {

  proj <- .attach_to_project(
    proj,
    obj,
    slot = "analyses",
    subtype = "gsva",
    log = list(
      source = source,
      method = method,
      n_pathways = nrow(score_matrix)
    )
  )
}

return(invisible(proj))

}

# ===============================================================================
# 7. Console summary
# =============================================================================
#' @method print gsva_result
#' @export

print.gsva_result <- function(x, ...) {

  cat("\nGSVA Analysis\n")
  cat("----------------------\n")
  cat("Method: ", x$params$method, "\n", sep = "")
  cat("Pathways: ", x$n_pathways, "\n", sep = "")
  cat("Samples: ", x$n_samples, "\n", sep = "")
  cat("\n")

}
