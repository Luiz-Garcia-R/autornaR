#' Gene Set Variation Analysis (GSVA)
#'
#' @description
#' Computes pathway activity scores for individual samples using GSVA-based
#' methods. Gene sets can be obtained directly from a previous
#' \code{rna.gsea()} analysis, MSigDB collections, or user-supplied pathways.
#'
#' Unlike GSEA, which evaluates pathway enrichment between groups,
#' GSVA estimates pathway activity for each sample individually.
#'
#' @param project rna_project object.
#' @param source Character. Source of gene sets:
#'   \code{gsea}, \code{msigdb} or \code{custom}.
#' @param gsea Optional gsea_result object when \code{source = 'gsea'}.
#' @param pathways Optional named list of gene sets when
#'   \code{source = 'custom'}.
#' @param pathway_scope Character. Determines search scope between:
#'   \code{top}, \code{significant}, or \code{all}
#' @param geneset_collection MSigDB collection when
#'   \code{source = "msigdb"}.
#' @param geneset_subcollection MSigDB subcollection.
#' @param method GSVA method:
#'   \code{gsva}, \code{ssgsea}, \code{zscore} or \code{plage}.
#' @param kcdf Kernel used by GSVA.
#' @param min_size Minimum pathway size.
#' @param max_size Maximum pathway size.
#' @param save Logical. Save result into project.
#'
#' @return
#' Object of class \code{gsva_result}.
#'
#' @export

rna.gsva <- function(
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

  # ---------------------------
  # Dependencies
  # ---------------------------
  bioc_pkgs <- c(
    "GSVA"
  )

  .check_dependencies(bioc_pkgs, bioc = TRUE)

  # ---------------------------
  # Load project
  # ---------------------------
  proj <- project

  expr <- as.matrix(.get_expr(proj))
  organism <- .get_organism(proj)
  gsea_obj <- .get_gsea(proj)

  # Harmonize gene identifiers
  expr <- .convert_expr_to_symbols(
    expr = expr,
    proj = proj
  )

  # ---------------------------
  # Obtain pathways
  # ---------------------------
  if (source == "gsea") {

    if (is.null(gsea_obj)) {
      stop(
        "No GSEA result found in project. Run rna.gsea() first."
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

  # msigdb
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

  # Custom
  if (source == "custom") {

    if (is.null(pathways)) {
      stop(
        "Provide 'pathways' when source = 'custom'."
      )
    }
  }

  # ---------------------------
  # Filter pathways
  # ---------------------------
  pathways <- pathways[
    lengths(pathways) >= min_size &
      lengths(pathways) <= max_size
  ]

  if (length(pathways) == 0) {
    stop("No pathways remaining after filtering.")
  }

  # ---------------------------
  # Run GSVA
  # ---------------------------
  if (method == "gsva") {

    param <- GSVA::gsvaParam(
      exprData = expr,
      geneSets = pathways,
      kcdf = kcdf
    )

  } else if (method == "ssgsea") {

    param <- GSVA::ssgseaParam(
      exprData = expr,
      geneSets = pathways
    )

  } else if (method == "zscore") {

    param <- GSVA::zscoreParam(
      exprData = expr,
      geneSets = pathways
    )

  } else if (method == "plage") {

    param <- GSVA::plageParam(
      exprData = expr,
      geneSets = pathways
    )

  }

  message(
    "[rna.gsva] ",
    length(pathways),
    " pathways retained after size filtering."
  )

  score_matrix <- GSVA::gsva(param)

  # ---------------------------
  # Output object
  # ---------------------------

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

  obj <- list(
    params = params,
    pathway_scores = score_matrix,
    pathways = pathways,
    pathway_metadata = pathway_metadata,
    n_pathways = nrow(score_matrix),
    n_samples = ncol(score_matrix)
  )

  class(obj) <- "gsva_result"

  # ---------------------------
  # Save into project
  # ---------------------------
  if (save) {

    proj <- .attach_to_project(
      proj,
      obj,
      slot = "analyses",
      subtype = "gsva",
      prefix = "gsva",
      log = list(
        source = source,
        method = method,
        n_pathways = nrow(score_matrix)
      )
    )
  }

  return(invisible(proj))

}


# ---------------------------
# Console summary
# ---------------------------
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
