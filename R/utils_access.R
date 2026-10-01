# ============================================================================ #
#                   ### Internal helpers for object access ###                 #
# ============================================================================ #

# =============================================================================
# 1. General project object accessor
# =============================================================================
#' @keywords internal

.get_project_object <- function(project,
                                slot,
                                subtype,
                                name = NULL) {

  if (is.null(project[[slot]])) {
    stop(sprintf("Project slot '%s' not found.", slot))
  }

  obj <- project[[slot]][[subtype]]

  if (is.null(obj)) {
    stop(sprintf(
      "No '%s' object found in project.",
      subtype
    ))
  }

  if (!is.null(name)) {

    if (!is.list(obj) || is.null(obj[[name]])) {
      stop(sprintf(
        "Field '%s' not found in '%s'.",
        name,
        subtype
      ))
    }

    return(obj[[name]])
  }

  obj
}

# =============================================================================
# 2. Normalized data
# =============================================================================

# Expr. data
.get_expr <- function(project) {
  .get_project_object(
    project,
    slot = "data",
    subtype = "normalized_data",
    name = "expr_matrix"
  )
}

# Metadata
.get_meta <- function(project) {
  .get_project_object(
    project,
    slot = "data",
    subtype = "normalized_data",
    name = "metadata"
  )
}

# Normalized datra
.get_norm_method <- function(project) {
  .get_project_object(
    project,
    slot = "data",
    subtype = "normalized_data",
    name = "method"
  )
}

# =============================================================================
# 3. Input data
# =============================================================================

# Imp data
.get_imp <- function(project) {
  .get_project_object(
    project,
    slot = "input",
    subtype = "imp_data"
  )
}

# Organism
.get_organism <- function(project) {
  .get_project_object(
    project,
    slot = "input",
    subtype = "imp_data",
    name = "organism"
  )
}

# Gene ID type
.get_gene_id_type <- function(project) {
  .get_project_object(
    project,
    slot = "input",
    subtype = "imp_data",
    name = "gene_id_type"
  )
}

# Annotation
.get_gene_annotation <- function(project) {
  .get_project_object(
    project,
    slot = "input",
    subtype = "imp_data",
    name = "gene_annotation"
  )
}

# =============================================================================
# 4. Analyses
# =============================================================================

# PCA / UMAP / T-SNE
.get_dimred <- function(project) {
  .get_project_object(
    project,
    slot = "analyses",
    subtype = "dimred"
  )
}

# DESeq2 / LIMMA
.get_comp <- function(project) {
  .get_project_object(
    project,
    slot = "analyses",
    subtype = "comparison"
  )
}

# GSEA
.get_gsea <- function(project) {
  .get_project_object(
    project,
    slot = "analyses",
    subtype = "gsea"
  )
}

# GSVA
.get_gsva <- function(project) {
  .get_project_object(
    project,
    slot = "analyses",
    subtype = "gsva"
  )
}

# Pathway scores
.get_gsva_scores <- function(project) {
  .get_project_object(
    project,
    slot = "analyses",
    subtype = "gsva",
    name = "pathway_scores"
  )
}

