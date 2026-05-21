#' Build DEG sets and structured gene membership from differential expression results
#'
#' @description
#' Constructs structured gene sets from differential expression results stored in
#' the active \code{rna_project}. This function formalizes DEGs as explicit sets,
#' enabling downstream analyses such as overlaps, enrichment, and custom visualization.
#'
#' Optionally, a summary plot of up- and down-regulated genes can be generated.
#'
#' @param project \code{rna_project} object created by \code{rna.project()}.
#' @param comparison Character. Name of a comparison stored in \code{rna_project}.
#'   If \code{NULL}, the most recent comparison is used.
#' @param padj_cutoff Numeric. Adjusted p-value threshold used to define significance
#'   (default: \code{0.05}).
#' @param log2fc_cutoff Numeric. Log2 fold-change threshold for directional classification
#'   (default: \code{0.5}).
#' @param direction Character. One of \code{"both"}, \code{"up"}, or \code{"down"}.
#'   Determines which genes are flagged as DEGs in the main output (default: \code{"both"}).
#' @param plot Logical. If \code{TRUE}, displays a barplot summarizing DEG counts per group
#'   (default: \code{FALSE}).
#' @param save Logical. Whether to store results in the active \code{rna_project}
#'   (default: \code{TRUE}).
#' @param verbose Logical. If \code{TRUE}, prints a summary to the console (default: \code{TRUE}).
#'
#' @return
#' An object of class \code{"rna_sets"} containing:
#' \describe{
#'   \item{deg}{Logical vector indicating DEG membership after applying \code{direction}.}
#'   \item{deg_all}{Logical vector indicating all DEGs (independent of direction).}
#'   \item{up}{Logical vector for genes upregulated in the test group.}
#'   \item{down}{Logical vector for genes upregulated in the reference group.}
#'   \item{one_hot}{Data frame encoding gene membership across groups (one-hot format).}
#'   \item{params}{List of parameters and thresholds used to define sets.}
#' }
#'
#' @details
#' This function transforms differential expression results into explicit gene sets
#' based on statistical significance and effect size thresholds.
#'
#' \strong{DEG definition}
#' Genes are considered differentially expressed if:
#' \itemize{
#'   \item \code{padj < padj_cutoff}
#'   \item absolute log2 fold-change exceeds \code{log2fc_cutoff}
#' }
#'
#' \strong{Directional classification}
#' \itemize{
#'   \item Upregulated genes: \code{log2FC > log2fc_cutoff}
#'   \item Downregulated genes: \code{log2FC < -log2fc_cutoff}
#' }
#'
#' \strong{Set representation}
#' The output includes multiple complementary representations:
#' \itemize{
#'   \item Binary DEG membership (\code{deg}, \code{deg_all})
#'   \item Directional sets (\code{up}, \code{down})
#'   \item One-hot encoded matrix (\code{one_hot}) for set operations
#' }
#'
#' This structure allows flexible downstream use, including intersections,
#' enrichment analysis, and visualization.
#'
#' @section Interpretation guide:
#' - The number of DEGs depends strongly on both \code{padj_cutoff} and \code{log2fc_cutoff}
#' - Increasing the fold-change threshold yields more conservative, high-effect genes
#' - Directional filtering is useful for asymmetric biological interpretations
#'
#' The distinction between "up" and "down" depends on the contrast definition
#' used in the original comparison.
#'
#' @section Pipeline context:
#' This function is typically used after \code{rna.compare()} to formalize DEG sets.
#' The resulting object is designed to be used with downstream tools such as:
#' \itemize{
#'   \item visualization (e.g., heatmaps, Venn diagrams)
#'   \item enrichment analysis
#'   \item custom set operations
#' }
#'
#' It represents a key step in separating statistical testing from biological
#' interpretation.
#'
#' @section Important considerations:
#' \strong{Threshold sensitivity}
#' DEG sets can vary substantially depending on cutoff choices. Always inspect
#' results across multiple thresholds when possible.
#'
#' \strong{Directionality depends on contrast}
#' The meaning of "up" and "down" is tied to the test vs reference groups defined
#' in the original comparison.
#'
#' \strong{Gene universe consistency}
#' Only genes present in both the expression matrix and comparison results
#' are retained.
#'
#' @section Side effects:
#' If \code{save = TRUE}, results are stored in the \code{rna_project} under the
#' \code{analyses$rna_sets} slot.
#'
#' @examples
#' \dontrun{
#' # Build DEG sets from last comparison
#' rna.sets(project = my_project)
#'
#' # More stringent thresholds
#' rna.sets(my_project,
#'          padj_cutoff = 0.01,
#'          log2fc_cutoff = 1
#' )
#'
#' # Only upregulated genes
#' rna.sets(my_project,
#'          direction = "up")
#'
#' # Plot DEG summary
#' rna.sets(my_project,
#'          plot = TRUE)
#'
#' }
#'
#' @export

rna.sets <- function(project,
                     comparison = NULL,
                     padj_cutoff = 0.05,
                     log2fc_cutoff = 0.5,
                     direction = c("both", "up", "down"),
                     plot = FALSE,
                     save = TRUE,
                     verbose = TRUE) {

  direction <- match.arg(direction)

  # ---------------------------
  # 0) Basic checks
  # ---------------------------
  .check_dependencies("ggplot2")

  # ---------------------------
  # 1) Get active project
  # ---------------------------
  proj <- project

  expr <- as.matrix(.get_expr(proj))
  metadata <- .get_meta(proj)
  comp_data <- .get_comp_obj(proj, comparison)

  # ---------------------------
  # 2) Validate input
  # ---------------------------
  if (is.null(proj$data$normalized_data)) {
    stop("No normalization results found in rna_project. Run rna.normalize() first.")
  }

  if (is.null(proj$analyses$comparison) ||
      length(proj$analyses$comparison) == 0) {
    stop("No comparisons found. Run rna.compare() first.")
  }

  # ---------------------------
  # 3) Accessing comp object
  # ---------------------------
  comparison_id <- if (is.null(comparison)) {
    proj$analyses$comparison$last
  } else {
    comparison
  }

  if (is.null(comp_data) || !inherits(comp_data, "rnaCompare")) {
    stop("Selected object is not a valid comparison.")
  }

  res <- as.data.frame(comp_data$res)
  res$Gene <- rownames(res)

  group_levels <- c(
    comp_data$groups$test,
    comp_data$groups$reference
  )

  if (length(group_levels) != 2) {
    stop("rna.sets() supports only 2-group comparisons.")
  }

  # Comparison label
  comparison_label <- paste(
    group_levels[1],
    "vs",
    group_levels[2]
  )

  # ---------------------------
  # 4) Align genes
  # ---------------------------
  common_genes <- intersect(rownames(expr), rownames(res))
  expr <- expr[common_genes, , drop = FALSE]
  res  <- res[common_genes, , drop = FALSE]

  # ---------------------------
  # 5) Detection filter
  # ---------------------------
  keep_genes <- rep(TRUE, length(common_genes))
  names(keep_genes) <- common_genes
  detection_matrix <- NULL

  # ---------------------------
  # 6) DEG classification
  # ---------------------------
  base_sig <- res$padj < padj_cutoff & !is.na(res$padj)

  up_sig   <- base_sig & res$log2FoldChange >  log2fc_cutoff
  down_sig <- base_sig & res$log2FoldChange < -log2fc_cutoff

  # --- Apply detection filter ---
  up_sig   <- up_sig   & keep_genes
  down_sig <- down_sig & keep_genes

  # --- Full DEG universe ---
  deg_all <- up_sig | down_sig

  # --- Optional directional filter for returned DEG column ---
  if (direction == "up") {
    deg_filtered <- up_sig
  } else if (direction == "down") {
    deg_filtered <- down_sig
  } else {
    deg_filtered <- deg_all
  }

  deg_membership <- data.frame(
    DEG = deg_filtered,
    row.names = common_genes
  )

  # ---------------------------
  # 7) One-hot by group detection
  # ---------------------------
  one_hot <- matrix(
    FALSE,
    nrow = length(common_genes),
    ncol = length(group_levels),
    dimnames = list(common_genes, group_levels)
  )

  one_hot[, group_levels[1]] <- up_sig   # test
  one_hot[, group_levels[2]] <- down_sig # reference

  one_hot_df <- as.data.frame(one_hot)

  # ---------------------------
  # 8) Optional DEG summary plot
  # ---------------------------
  if (plot) {

    n_up_test <- sum(up_sig & keep_genes)
    n_up_ref  <- sum(down_sig & keep_genes)

    plot_df <- data.frame(
      group = c(
        paste0("Up in ", group_levels[2]),
        paste0("Up in ", group_levels[1])
      ),
      value = c(n_up_test, n_up_ref)
    )

    vivid_colors <- scales::hue_pal()(length(unique(plot_df$group)))

    g <- ggplot2::ggplot(plot_df,
                         ggplot2::aes(x = .data$group,
                                      y = .data$value,
                                      fill = .data$group)) +
      ggplot2::geom_col(alpha = 0.8) +
      ggplot2::geom_text(
        ggplot2::aes(label = .data$value),
        vjust = -0.2,
        size = 4
      ) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::scale_fill_manual(values = vivid_colors) +
      ggplot2::labs(
        title = paste("DEG Summary -", comparison_label),
        x = "",
        y = "DEGs count"
      ) +
      ggplot2::theme(
        legend.position = "none",
        axis.text.x = ggplot2::element_text(
          angle = 45,
          hjust = 1,
          size = 12
        )
      )

    print(g)
  }

  # ---------------------------
  # 9) Output
  # ---------------------------
  obj <- list(
    params = list(
      timestamp = Sys.time(),
      comparison = comparison_id,
      comparison_label = comparison_label,
      total_genes = length(common_genes),
      total_deg = sum(deg_all),
      padj_cutoff = padj_cutoff,
      log2fc_cutoff = log2fc_cutoff,
      direction = direction
    ),
    membership = list(
      deg = deg_membership,
      deg_all = deg_all,
      up = up_sig,
      down = down_sig,
      one_hot = one_hot_df
    )
  )

  class(obj) <- "rna_sets"

  # ---------------------------
  # 10) Attach to project
  # ---------------------------
  if (save) {

    proj <- .attach_to_project(
      proj,
      obj,
      slot = "analyses",
      subtype = "sets",
      prefix = "sets",
      log = list(
        comparison = comparison_id,
        direction = direction,
        log2fc_cutoff = log2fc_cutoff,
        padj_cutoff = padj_cutoff,
        total_deg = sum(deg_all)
      )
    )
  }

  # ---------------------------
  # 11) Return
  # ---------------------------
  if (verbose) {
    .print_header("RNA DEG Sets")

    .print_block("Comparison", function() {
      cat(comparison_label, "\n")
    })

    .print_block("DEG Summary", function() {
      cat("Genes tested:      ", length(common_genes), "\n")
      cat("Total DEGs:        ", sum(deg_all), "\n")
      cat("Up in ", group_levels[2], ":           ", sum(up_sig), "\n", sep = "")
      cat("Up in ", group_levels[1], ":           ", sum(down_sig), "\n", sep = "")

      if (direction != "both") {
        cat("Directional filter:      ", direction, "\n")
        cat("Returned DEGs:           ", sum(deg_filtered), "\n")
      }
    })
  }

  return(invisible(proj))

}
