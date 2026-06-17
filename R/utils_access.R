# ============================
# Internal helpers for object access
# ============================

# --- General helper ---
.get_from_container <- function(container, id = NULL, name = NULL) {

  if (is.null(container)) {
    stop("Container not found.")
  }

  # --- resolve id ---
  if (is.null(id) || identical(id, "last")) {

    if (is.null(container$last)) {
      stop("No 'last' entry available in container.")
    }

    obj <- container[[container$last]]

  } else if (is.numeric(id)) {

    if (length(id) != 1) {
      stop("'id' must be a single numeric index.")
    }

    ids <- setdiff(names(container), "last")

    if (id < 1 || id > length(ids)) {
      stop("Index out of bounds.")
    }

    obj <- container[[ids[id]]]

  } else if (is.character(id)) {

    if (!id %in% names(container)) {
      stop("Invalid ID.")
    }

    obj <- container[[id]]

  } else {
    stop("Invalid 'id' argument.")
  }

  # --- optional extraction ---
  if (!is.null(name)) {
    if (!name %in% names(obj)) {
      stop(paste0("Field '", name, "' not found in object."))
    }
    return(obj[[name]])
  }

  return(obj)
}

# ============================
# 1) Normalized data
# ============================
# --- expr_matrix ---
.get_expr <- function(project, id = NULL) {
  .get_from_container(project$data$normalized_data, id, "expr_matrix")
  }

# --- metadata ---
  .get_meta <- function(project, id = NULL) {
    .get_from_container(project$data$normalized_data, id, "metadata")
  }

# --- norm_method ---
  .get_norm_method <- function(project, id = NULL) {
    .get_from_container(project$data$normalized_data, id, "method")
  }

# ============================
# 2) Input data
# ============================
# --- imp_data ---
  .get_imp <- function(project, id = NULL) {
    .get_from_container(project$input$imp_data, id)
  }

# --- organism ---
  .get_organism <- function(project, id = NULL) {
    .get_from_container(project$input$imp_data, id, "organism")
  }

# --- gene_id_type ---
  .get_gene_id_type <- function(project, id = NULL) {
    .get_from_container(project$input$imp_data, id, "gene_id_type")
  }

# ============================
# 3) Analyses
# ============================
# --- dimred_data ---
  .get_dimred <- function(project, id = NULL) {
    .get_from_container(project$analyses$dimred, id)
  }

# --- comp_data (Whole container) ---
  .get_comp <- function(project) {
    project$analyses$comparison
  }

# --- comp_data (Unique object) ---
  .get_comp_obj <- function(project, id = NULL) {
    .get_from_container(project$analyses$comparison, id)
  }

# --- gene_annotation ---
  .get_gene_annotation <- function(project, id = NULL) {
    .get_from_container(project$input$imp_data, id, "gene_annotation")
  }

# --- gsea ---
  .get_gsea <- function(project, id = NULL) {
    .get_from_container(project$analyses$gsea, id)
  }

# --- gsva ---
  .get_gsva <- function(project, id = NULL) {
    .get_from_container(
      project$analyses$gsva, id)
  }

# --- gsva scores ---
  .get_gsva_scores <- function(project, id = NULL) {
    .get_from_container(
      project$analyses$gsva,
      id,
      "pathway_scores"
    )
  }
