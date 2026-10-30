
#' Construct a weighted gene co-expression network and detect modules
#'
#' @description
#' Constructs a weighted gene co-expression network using normalized
#' expression data stored in an \code{rna_project} object and identifies
#' gene co-expression modules using the WGCNA framework.
#'
#' The function filters genes with zero or invalid variance, optionally
#' retains the most variable genes, performs expression-data quality
#' control, constructs the network using a predefined soft-thresholding
#' power, and detects modules through hierarchical clustering of the
#' topological overlap matrix (TOM).
#'
#' Module assignments, module eigengenes, network information, and
#' analysis metadata are organized into a result object. Diagnostic
#' dendrogram plots can also be generated to visualize the identified
#' modules.
#'
#' @param project An \code{rna_project} object containing normalized
#'   gene expression data and sample metadata.
#'
#' @param power Numeric value specifying the soft-thresholding power
#'   used to construct the network. This value should be selected
#'   beforehand using \code{\link{rna_pick_power}} or another
#'   appropriate diagnostic approach. Must be finite and positive.
#'
#' @param genes_n Integer or \code{NULL}. Maximum number of genes to
#'   retain, selected by decreasing expression variance across samples.
#'   Genes with non-finite or zero variance are excluded before selection.
#'   If \code{NULL}, all eligible genes are retained. Using a subset of
#'   the most variable genes can reduce computational requirements.
#'   Defaults to \code{10000}.
#'
#' @param network_type Character string specifying the network type
#'   used by WGCNA. Supported values are \code{"signed"},
#'   \code{"unsigned"}, and \code{"signed hybrid"}. A signed network
#'   preserves the distinction between positive and negative
#'   correlations; an unsigned network uses the magnitude of
#'   correlations; and a signed hybrid network excludes negative
#'   correlations from the adjacency calculation.
#'   Defaults to \code{"signed"}.
#'
#' @param tom_type Character string specifying the type of topological
#'   overlap matrix (TOM) used for module detection. Supported values
#'   are \code{"signed"}, \code{"unsigned"}, \code{"unsigned 2"},
#'   \code{"signed Nowick"}, and \code{"signed Nowick 2"}.
#'   This parameter controls how topological overlap is calculated
#'   and is independent of the adjacency network type specified by
#'   \code{network_type}. Defaults to \code{"signed"}.
#'
#' @param cor_method Character string specifying the correlation
#'   method used during network construction. Currently supports
#'   \code{"pearson"}, which is also used in downstream calculations
#'   of module membership (kME) to ensure consistency with the
#'   original network. Defaults to \code{"pearson"}.
#'
#' @param max_block_size Integer specifying the maximum number of genes
#'   allowed in each block during network construction. When the number
#'   of genes exceeds this limit, WGCNA divides the data into multiple
#'   blocks and constructs networks separately. Larger values may
#'   increase memory requirements. Defaults to \code{20000}.
#'
#' @param min_module_size Integer specifying the minimum number of genes
#'   required to form a module during dynamic tree cutting.
#'   Smaller modules are generally not retained as independent modules.
#'   Defaults to \code{30}.
#'
#' @param deep_split Integer between 0 and 4 controlling the sensitivity
#'   of dynamic tree cutting. Higher values generally produce more
#'   modules by allowing finer subdivisions of the gene dendrogram.
#'   Higher values correspond to more aggressive module detection.
#'   Defaults to \code{2}.
#'
#' @param merge_cut_height Numeric value between 0 and 1 specifying
#'   the threshold for merging modules based on the dissimilarity
#'   between their module eigengenes. Lower values result in more
#'   conservative merging, whereas higher values allow more modules
#'   to be merged. Defaults to \code{0.25}.
#'
#' @param reassign_threshold Numeric value between 0 and 1 specifying
#'   the threshold used for gene reassignment based on module
#'   membership significance during module detection. A value of
#'   \code{0} disables reassignment based on this threshold.
#'   Defaults to \code{1e-6}.
#'
#' @param numeric_labels Logical. If \code{TRUE}, modules are
#'   represented by numeric labels in the network result. If
#'   \code{FALSE}, modules are represented by their assigned colors.
#'   The returned \code{genes_modules} table also includes module
#'   colors regardless of this setting. Defaults to \code{FALSE}.
#'
#' @param pam_respects_dendro Logical. Specifies whether the
#'   partitioning around medoids (PAM) stage of module detection
#'   must respect the dendrogram produced by hierarchical clustering.
#'   If \code{TRUE}, PAM assignments are constrained by the dendrogram;
#'   if \code{FALSE}, genes may be assigned to modules independently
#'   of their position in the dendrogram. Defaults to \code{FALSE}.
#'
#' @param threads Integer or \code{NULL}. Number of threads used
#'   during network construction. If \code{NULL}, the \code{nThreads}
#'   argument is omitted and WGCNA uses its default behavior.
#'   A positive integer specifies the number of threads passed to
#'   the network construction procedure. Defaults to \code{NULL}.
#'
#' @param plot Logical. If \code{TRUE}, generates diagnostic
#'   dendrogram plots showing the module color assignments for
#'   each network block. Defaults to \code{TRUE}.
#'
#' @param save Logical. If \code{TRUE}, attaches the network result
#'   to the project under the \code{analyses} slot with subtype
#'   \code{"wgcna_network"} and invisibly returns the updated project.
#'   If \code{FALSE}, does not attach the result and instead invisibly
#'   returns an object of class \code{wgcna_network_result}.
#'   Defaults to \code{TRUE}.
#'
#' @section Input data:
#' The expression matrix must be numeric, with genes in rows and
#' samples in columns. Gene identifiers and sample identifiers must
#' be unique.
#'
#' The metadata must identify samples through a \code{Sample} column
#' or through its row names. Sample identifiers are matched to the
#' expression matrix column names; metadata must contain an entry
#' for every expression sample.
#'
#' Genes with non-finite or zero variance are excluded before the
#' most variable genes are selected. WGCNA quality control is then
#' applied to identify and remove additional problematic samples
#' or genes.
#'
#' @section Network construction:
#' The network is constructed using the adjacency function and
#' topological overlap matrix defined by \code{network_type} and
#' \code{tom_type}, respectively. The soft-thresholding power
#' specified by \code{power} is applied to the adjacency calculation.
#'
#' When the number of selected genes exceeds \code{max_block_size},
#' WGCNA may divide the data into multiple blocks. Module detection
#' is performed independently within each block, and the resulting
#' module assignments are returned in a unified result object.
#'
#' Modules are identified through hierarchical clustering and
#' dynamic tree cutting. Closely related modules may subsequently
#' be merged according to \code{merge_cut_height}.
#'
#' @section Module eigengenes:
#' Module eigengenes are calculated by WGCNA as part of the
#' module detection procedure. They represent the first principal
#' component of the expression profiles of genes assigned to each
#' module and provide a summarized representation of module-level
#' expression patterns.
#'
#' The returned eigengene matrix is ordered using
#' \code{\link[WGCNA]{orderMEs}}. Additional analyses, such as
#' module membership (kME) estimation and associations between
#' module eigengenes and sample-level phenotypes, can be performed
#' in subsequent steps.
#'
#' @section Interpretation:
#' Identified modules represent groups of genes with similar
#' expression patterns across the samples included in the analysis.
#' Module detection depends on the selected soft-thresholding
#' power, network type, TOM definition, and module detection
#' parameters.
#'
#' Module colors are arbitrary identifiers assigned by WGCNA
#' and do not have intrinsic biological meaning. Biological
#' interpretation requires additional analyses, such as functional
#' enrichment, module membership evaluation, and associations
#' with relevant phenotypes.
#'
#' The resulting modules describe co-expression relationships
#' across the samples supplied in the project as a whole.
#' They do not directly represent differential co-expression
#' between experimental groups or establish causal relationships
#' between genes.
#'
#' @section Scope of the analysis:
#' This function constructs a single co-expression network using
#' all retained samples. It does not perform differential
#' co-expression analysis, consensus network construction,
#' module preservation analysis, or phenotype association testing.
#' These analyses require additional steps.
#'
#' @return
#' If \code{save = TRUE}, invisibly returns the updated
#' \code{rna_project} object containing the network result.
#'
#' If \code{save = FALSE}, invisibly returns an object of class
#' \code{wgcna_network_result}, containing:
#' \itemize{
#'   \item \code{network}: the complete network object returned
#'   by \code{\link[WGCNA]{blockwiseModules}}.
#'   \item \code{genes_modules}: a data frame containing gene
#'   identifiers, module labels, and module colors.
#'   \item \code{MEs}: the ordered module eigengene matrix.
#'   \item Network construction parameters, including the
#'   soft-thresholding power and network types.
#'   \item The number and identifiers of retained genes and samples.
#'   \item Identifiers of genes excluded due to invalid variance,
#'   gene selection, or WGCNA quality control, and samples removed
#'   during quality control.
#'   \item \code{timestamp}: the date and time of analysis.
#' }
#'
#' @seealso
#' \code{\link{rna_pick_power}},
#' \code{\link[WGCNA]{blockwiseModules}},
#' \code{\link[WGCNA]{goodSamplesGenes}},
#' \code{\link[WGCNA]{plotDendroAndColors}},
#' \code{\link[WGCNA]{orderMEs}}
#'
#' @references
#' Langfelder P, Horvath S. (2008).
#' WGCNA: an R package for weighted correlation network analysis.
#' \emph{BMC Bioinformatics}, 9, 559.
#' \doi{10.1186/1471-2105-9-559}
#'
#' @export

rna_wgcna <- function(
    project,
    power,
    genes_n = 10000,
    network_type = c("signed", "unsigned", "signed hybrid"),
    tom_type = c("signed", "unsigned", "unsigned 2", "signed Nowick", "signed Nowick 2"),
    cor_method = "pearson",
    max_block_size = 20000,
    min_module_size = 30,
    deep_split = 2,
    merge_cut_height = 0.25,
    reassign_threshold = 1e-6,
    numeric_labels = FALSE,
    pam_respects_dendro = FALSE,
    threads = NULL,
    plot = TRUE,
    save = TRUE
) {

# =============================================================================
# 0. Basic checks
# =============================================================================

.check_dependencies("WGCNA")

# project
if (!inherits(project, "rna_project")) {
  stop("`project` must be an object of class 'rna_project'.")
}

# power
if (!is.numeric(power) ||
    length(power) != 1L ||
    !is.finite(power) ||
    power <= 0) {
  stop("`power` must be a single finite positive numeric value.")
}

# genes_n
if (!is.null(genes_n)) {
  if (!is.numeric(genes_n) ||
      length(genes_n) != 1L ||
      !is.finite(genes_n) ||
      genes_n < 1 ||
      genes_n != floor(genes_n)) {
    stop("`genes_n` must be a positive integer or NULL.")
  }

  genes_n <- as.integer(genes_n)
}

# network_type
network_type <- match.arg(network_type)

# tom_type
tom_type <- match.arg(tom_type)

# cor_method
cor_method <- match.arg(cor_method)

# max_block_size
if (!is.numeric(max_block_size) ||
    length(max_block_size) != 1L ||
    !is.finite(max_block_size) ||
    max_block_size < 1 ||
    max_block_size != floor(max_block_size)) {
  stop("`max_block_size` must be a positive integer.")
}

max_block_size <- as.integer(max_block_size)

# min_module_size
if (!is.numeric(min_module_size) ||
    length(min_module_size) != 1L ||
    !is.finite(min_module_size) ||
    min_module_size < 1 ||
    min_module_size != floor(min_module_size)) {
  stop("`min_module_size` must be a positive integer.")
}

min_module_size <- as.integer(min_module_size)

if (min_module_size > max_block_size) {
  stop("`min_module_size` cannot exceed `max_block_size`.")
}

# deep_split
if (!is.numeric(deep_split) ||
    length(deep_split) != 1L ||
    !is.finite(deep_split) ||
    deep_split < 0 ||
    deep_split > 4 ||
    deep_split != floor(deep_split)) {
  stop("`deep_split` must be an integer between 0 and 4.")
}

deep_split <- as.integer(deep_split)

# merge_cut_height
if (!is.numeric(merge_cut_height) ||
    length(merge_cut_height) != 1L ||
    !is.finite(merge_cut_height) ||
    merge_cut_height < 0 ||
    merge_cut_height > 1) {
  stop("`merge_cut_height` must be a finite numeric value between 0 and 1.")
}

# reassign_threshold
if (!is.numeric(reassign_threshold) ||
    length(reassign_threshold) != 1L ||
    !is.finite(reassign_threshold) ||
    reassign_threshold < 0 ||
    reassign_threshold > 1) {
  stop("`reassign_threshold` must be a finite numeric value between 0 and 1.")
}

# numeric_labels
if (!is.logical(numeric_labels) ||
    length(numeric_labels) != 1L ||
    is.na(numeric_labels)) {
  stop("`numeric_labels` must be TRUE or FALSE.")
}

# pam_respects_dendro
if (!is.logical(pam_respects_dendro) ||
    length(pam_respects_dendro) != 1L ||
    is.na(pam_respects_dendro)) {
  stop("`pam_respects_dendro` must be TRUE or FALSE.")
}

# threads
if (!is.null(threads)) {
  if (!is.numeric(threads) ||
      length(threads) != 1L ||
      !is.finite(threads) ||
      threads < 1 ||
      threads != floor(threads)) {
    stop("`threads` must be a positive integer or NULL.")
  }

  threads <- as.integer(threads)
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
# 4. Build network and detect modules
# =============================================================================

cor <- WGCNA::cor

wgcna_args <- list(
  datExpr = datExpr,
  power = power,
  networkType = network_type,
  TOMType = tom_type,
  maxBlockSize = max_block_size,
  minModuleSize = min_module_size,
  deepSplit = deep_split,
  reassignThreshold = reassign_threshold,
  mergeCutHeight = merge_cut_height,
  numericLabels = numeric_labels,
  pamRespectsDendro = pam_respects_dendro,
  saveTOMs = FALSE,
  verbose = 3
)

if (!is.null(threads)) {
  wgcna_args$nThreads <- threads
}

net <- do.call(
  WGCNA::blockwiseModules,
  wgcna_args
)

# =============================================================================
# 5. Organize module assignments
# =============================================================================

if (numeric_labels) {
  module_colors <- WGCNA::labels2colors(net$colors)
} else {
  module_colors <- net$colors
}

genes_modules <- data.frame(
  Gene = colnames(datExpr),
  Module = net$colors,
  ModuleColor = module_colors,
  stringsAsFactors = FALSE
)

# =============================================================================
# 6. Prepare module eigengenes
# =============================================================================

MEs <- WGCNA::orderMEs(net$MEs)

# =============================================================================
# 7. Diagnostic plots
# =============================================================================

if (plot) {

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)

  n_blocks <- length(net$dendrograms)

  graphics::par(
    mfrow = c(n_blocks, 1)
  )

  for (b in seq_len(n_blocks)) {

    block_idx <- net$blockGenes[[b]]

    block_colors <- module_colors[block_idx]

    WGCNA::plotDendroAndColors(
      net$dendrograms[[b]],
      block_colors,
      "Module colors",
      dendroLabels = FALSE,
      hang = 0.03,
      main = paste("Co-expression network - Block", b)
    )
  }
}

# =============================================================================
# 8. Store result
# =============================================================================

obj <- list(
  network = net,
  datExpr = datExpr,
  genes_modules = genes_modules,
  MEs = MEs,
  power = power,
  network_type = network_type,
  tom_type = tom_type,
  cor_method = cor_method,
  n_genes = ncol(datExpr),
  n_samples = nrow(datExpr),
  genes_used = colnames(datExpr),
  samples_used = rownames(datExpr),
  genes_invalid_variance = prepared$genes_invalid_variance,
  genes_excluded_selection = prepared$genes_excluded_selection,
  genes_removed_qc = prepared$genes_removed_qc,
  samples_removed_qc = prepared$samples_removed_qc,
  timestamp = Sys.time()
)

class(obj) <- "wgcna_network_result"

# =============================================================================
# 9. Console summary
# =============================================================================

.print_header("WGCNA Network Results")

.print_block("Network Summary", function() {

  cat("Network type:          ", network_type, "\n")
  cat("TOM type:              ", tom_type, "\n")
  cat("Soft-threshold power:  ", power, "\n")
  cat("Samples used:          ", obj$n_samples, "\n")
  cat("Genes used:            ", obj$n_genes, "\n")
  cat("Blocks:                ", length(net$dendrograms), "\n")

})

# ===============================================
# Module detection
# ===============================================

module_table <- table(obj$genes_modules$ModuleColor)

# Exclude genes assigned to the unassigned module (grey).
module_sizes <- module_table[names(module_table) != "grey"]

n_modules <- length(module_sizes)

n_unassigned <- sum(
  obj$genes_modules$ModuleColor == "grey"
)

.print_block("Module Detection", function() {

  cat("Modules identified:    ", n_modules, "\n")
  cat("Unassigned genes:      ", n_unassigned, "\n")

})

# ===============================================
# Module sizes
# ===============================================

.print_block("Module Sizes", function() {

  if (n_modules > 0) {

    module_sizes <- sort(module_sizes, decreasing = TRUE)

    for (i in seq_along(module_sizes)) {

      cat(
        names(module_sizes)[i],
        ": ",
        module_sizes[i],
        " genes\n",
        sep = ""
      )

    }

  } else {

    cat("No co-expression modules were identified.\n")

  }

})

# =============================================================================
# 10. Attach to project
# =============================================================================

if (save) {

  proj <- .attach_to_project(
    proj,
    obj,
    slot = "analyses",
    subtype = "wgcna_network",
    log = list(
      power = power,
      network_type = network_type,
      tom_type = tom_type,
      n_genes = obj$n_genes,
      n_samples = obj$n_samples
    )
  )

  return(invisible(proj))
}

return(invisible(obj))
}
