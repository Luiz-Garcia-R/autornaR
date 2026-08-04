#' Differential expression analysis using DESeq2 or limma-voom
#'
#' @description
#' Performs differential expression analysis between experimental groups
#' using RNA-seq count data stored in the active \code{rna_project}.
#' Supports both \pkg{DESeq2} and \pkg{limma} (voom) frameworks.
#'
#' The function handles design matrix construction, contrast specification,
#' and result extraction, returning standardized outputs across methods.
#'
#' @param project \code{rna_project} object created by \code{rna.project()}.
#' @param group_col Character. Column name in metadata defining experimental
#'   groups (default: \code{"Group"}).
#' @param batch_col Optional character. Column specifying batch effects.
#'   If present in metadata, it is included in the design formula.
#' @param method Character. Differential expression method:
#'   \code{"deseq2"} or \code{"limma"}.
#' @param contrast Specifies the comparison(s) to perform. Can be:
#'   \itemize{
#'     \item A character vector of length 2: \code{c("test", "reference")}
#'     \item A list of such vectors for multiple contrasts
#'   }
#' @param plot Logical. If \code{TRUE}, displays diagnostic plots:
#'   MA plot (both methods) and dispersion plot (DESeq2 only).
#' @param save Logical. Whether to store results in the active \code{rna_project}.
#' @param store_model Logical. Wheter to store the complete model results in the active \code{rna_project}.
#' @param verbose Logical. Whether to print progress messages.
#'
#' @details
#' The function requires raw count data and expects that normalization has been
#' performed with \code{rna.normalize(method = "none")} to preserve integer counts.
#'
#' \strong{Contrast definition:}
#' All contrasts are defined as \strong{test − reference}, which determines
#' the sign of the log2 fold-change:
#' \itemize{
#'   \item Positive log2FC -> higher expression in test group
#'   \item Negative log2FC -> higher expression in reference group
#' }
#'
#' \strong{Method-specific details:}
#' \itemize{
#'   \item \pkg{DESeq2}: Uses negative binomial GLM with dispersion estimation.
#'   \item \pkg{limma}: Uses voom transformation followed by linear modeling
#'   with empirical Bayes moderation.
#' }
#'
#' Multiple contrasts can be specified and will be computed independently.
#' Each result is stored with a unique identifier in:
#' \code{rna_project$analyses$comparison}.
#'
#' @return
#' Returns (invisibly) either:
#' \itemize{
#'   \item A single object of class \code{"rnaCompare"} if one contrast is specified
#'   \item A named list of such objects if multiple contrasts are provided
#' }
#'
#' Each \code{"rnaCompare"} object contains:
#' \itemize{
#'   \item \code{model}: Fitted model object (\code{DESeqDataSet} or \code{MArrayLM})
#'   \item \code{res}: Data frame with differential expression results
#'   \item \code{contrast}: Contrast specification used
#'   \item \code{groups}: List with \code{test} and \code{reference}
#'   \item \code{method}: Method used ("deseq2" or "limma")
#'   \item \code{plots}: List of plotting functions (MA and dispersion)
#' }
#'
#' @section Output summary:
#' A summary is printed to the console including:
#' \itemize{
#'   \item Total genes tested
#'   \item Number of significant genes (adjusted p-value < 0.05)
#'   \item Number of up- and down-regulated genes
#'   \item Top-ranked genes
#' }
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' rna.compare(my_project,
#'             contrast = c("treated", "control"))
#'
#' # Multiple contrasts
#' rna.compare(my_project,
#'             contrast = list(
#'                             c("treated", "control"),
#'                             c("treated2", "control")
#'   )
#' )
#'
#' # Using limma instead of DESeq2
#' rna.compare(my_project,
#'             contrast = c("treated", "control"),
#'             method = "limma")
#'
#' # Including batch effect
#' rna.compare(my_project,
#'             contrast = c("treated", "control"),
#'             batch_col = "Batch")
#' }
#'
#' @importFrom dplyr rename
#' @importFrom stats model.matrix aov
#'
#' @export

rna.compare <- function(project,
                        group_col = "Group",
                        batch_col = "Batch",
                        method = c("deseq2", "limma"),
                        contrast = NULL,
                        plot = TRUE,
                        save = TRUE,
                        store_model = FALSE,
                        verbose = TRUE
) {

  method <- match.arg(method)

  # ===========================================================================
  # 0) Basic checks
  # ===========================================================================
  if (method == "deseq2") {
    .check_dependencies(c("DESeq2", "SummarizedExperiment"), bioc = TRUE)
  }

  if (method == "limma") {
    .check_dependencies(c("limma", "edgeR"), bioc = TRUE)
  }

  # ===========================================================================
  # 1) get active project
  # ===========================================================================
  proj <- project

  expr_mat <- as.matrix(.get_expr(proj))
  metadata <- .get_meta(proj)
  norm_method <- .get_norm_method(proj)

  # ===========================================================================
  # 2) Validate input
  # ===========================================================================
  # --- Factor handling ---
  if (!is.factor(metadata[[group_col]])) {
    if (verbose)
      message(sprintf(
        "[rna.compare] '%s' is not a factor. Converting using order of appearance.",
        group_col
      ))
    metadata[[group_col]] <-
      factor(metadata[[group_col]], levels = unique(metadata[[group_col]]))
  }

  groups <- levels(metadata[[group_col]])

  if (is.null(contrast)) {

    if (length(groups) != 2) {
      stop(
        "More than two groups detected.\n",
        "Please specify contrasts manually using contrast = c('group1','group2') ",
        "or a list of contrasts."
      )
    }

    stop(
      "Contrast direction was not specified.\n\n",
      "Detected groups: ",
      paste(groups, collapse = " and "),
      ".\n\n",
      "Please specify the comparison direction explicitly.\n",
      "Examples:\n",
      "  contrast = c('", groups[2], "','", groups[1], "')\n",
      "  contrast = c('", groups[1], "','", groups[2], "')\n\n",
      "The first group is treated as the test group and the second as the reference group."
    )

  } else {

    if (is.character(contrast))
      contrast_list <- list(contrast)

    if (is.list(contrast))
      contrast_list <- contrast

  }

  # ===========================================================================
  # 3) Reading counting matrix
  # ===========================================================================
  counts <- round(expr_mat)

  if (norm_method != "none") {
    stop(
      "Differential expression requires raw counts.\n",
      "Please run rna.normalize(method = 'none') before rna.compare()."
    )
  }

  if (any(grepl("\\.", rownames(counts)))) {
    stop(
      "Gene IDs contain version suffixes.\n",
      "Please run rna.normalize(clean_gene_versions = TRUE) upstream."
    )
  }

  metadata <- metadata[match(colnames(counts), metadata$Sample), ]
  rownames(metadata) <- metadata$Sample

  # --- Design formula ---
  design_formula <- if (!is.null(batch_col) && batch_col %in% colnames(metadata)) {
    as.formula(paste("~", batch_col, "+", group_col))
  } else {
    as.formula(paste("~", group_col))
  }

  # ===========================================================================
  # 4) DESeq2
  # ===========================================================================
  if (method == "deseq2") {
    if (verbose) message("[rna.compare] Creating DESeq2 dataset...")

    dds <- DESeq2::DESeqDataSetFromMatrix(
      countData = counts,
      colData = metadata,
      design = design_formula
    )

    if (verbose) message("[rna.compare] Running DESeq2...")

    dds_final <- tryCatch(
      {
        DESeq2::DESeq(dds)
      },
      error = function(e) {
        if (verbose)
          message("[rna.compare] Parametric fit failed. Using gene-wise dispersions.")

        dds <- DESeq2::estimateSizeFactors(dds)
        dds <- DESeq2::estimateDispersionsGeneEst(dds)
        DESeq2::dispersions(dds) <- SummarizedExperiment::mcols(dds)$dispGeneEst
        dds
      }
    )

    results_all <- list()

    for (ct in contrast_list) {

      test_group <- ct[1]
      ref_group  <- ct[2]

      if (verbose) {
        message(sprintf("[rna.compare] Running contrast: %s vs %s",
                        test_group, ref_group))
      }

      contrast_use <- c(group_col, test_group, ref_group)

      res <- DESeq2::results(dds_final, contrast = contrast_use)
      res_df <- as.data.frame(res[order(res$pvalue), ])

      comparison_id <- paste(test_group, "vs", ref_group, sep = "_")

      ma_plot <- function() { DESeq2::plotMA(res, ylim=c(-5,5)) }
      disp_plot <- function() { DESeq2::plotDispEsts(dds_final) }


      # Store results
      if (store_model) {

        comparison_obj <- list(
          timestamp = Sys.time(),
          model = dds_final,
          res = res_df,
          contrast = contrast_use,
          groups = list(reference = ref_group, test = test_group),
          method = method,
          plots = list(MA = ma_plot, Dispersion = disp_plot)
        )

      } else {

        comparison_obj <- list(
          timestamp = Sys.time(),
          res = res_df,
          contrast = contrast_use,
          groups = list(reference = ref_group, test = test_group),
          method = method
        )

      }

      class(comparison_obj) <- "rnaCompare"
      results_all[[comparison_id]] <- comparison_obj
    }

  # ===========================================================================
  # 5) LIMMA
  # ===========================================================================
  } else if (method == "limma") {

    design <- model.matrix(
      as.formula(paste("~ 0 +",
                       if (!is.null(batch_col) && batch_col %in% colnames(metadata))
                         paste(batch_col, "+", group_col)
                       else
                         group_col)),
      data = metadata
    )

    dge <- edgeR::DGEList(counts = counts)
    dge <- edgeR::calcNormFactors(dge)

    v <- limma::voom(dge, design, plot = FALSE)

    fit <- limma::lmFit(v, design)

    results_all <- list()

    for (ct in contrast_list) {

      test_group <- ct[1]
      ref_group  <- ct[2]

      if (verbose) {
        message(sprintf("[rna.compare] Running contrast: %s vs %s",
                        test_group, ref_group))
      }

      # expected columns names
      coef_test <- paste0(group_col, test_group)
      coef_ref  <- paste0(group_col, ref_group)

      if (!(coef_test %in% colnames(design)))
        stop("Test group not found in design matrix.")

      if (!(coef_ref %in% colnames(design)))
        stop("Reference group not found in design matrix.")

      # create explicit contrast: test - ref
      contrast_formula <- paste0(coef_test, " - ", coef_ref)

      contrast_matrix <- limma::makeContrasts(
        contrasts = contrast_formula,
        levels = design
      )

      fit2 <- limma::contrasts.fit(fit, contrast_matrix)
      fit2 <- limma::eBayes(fit2)

      res_df <- limma::topTable(
        fit2,
        number = Inf,
        sort.by = "P"
      )

      res_df <- res_df |>
        dplyr::rename(
          log2FoldChange = "logFC",
          stat = "t",
          pvalue = "P.Value",
          padj = "adj.P.Val",
          log_odds = "B"
        )

      comparison_id <- paste(test_group, "vs", ref_group, sep = "_")

      ma_plot <- function() {
        limma::plotMD(fit2)
      }

      # Store results
      if (store_model) {

        comparison_obj <- list(
          timestamp = Sys.time(),
          model = fit2,
          res = res_df,
          contrast = contrast_formula,
          groups = list(reference = ref_group, test = test_group),
          method = method,
          plots = list(MA = ma_plot, Dispersion = NULL)
        )

      } else {

        comparison_obj <- list(
          timestamp = Sys.time(),
          res = res_df,
          contrast = contrast_formula,
          groups = list(reference = ref_group, test = test_group),
          method = method
        )

      }

      class(comparison_obj) <- "rnaCompare"

      results_all[[comparison_id]] <- comparison_obj
    }
  }

  # ===========================================================================
  # 6) Plot functions
  # ===========================================================================
  if (method == "deseq2") {
    ma_plot <- function() { DESeq2::plotMA(res, ylim=c(-5,5)) }
    disp_plot <- function() { DESeq2::plotDispEsts(dds_final) }
  } else {
    ma_plot <- function() {
      limma::plotMD(fit2)
    }
    disp_plot <- NULL
  }

  # --- Primary comparison ---
  primary_name <- names(results_all)[1]
  primary_obj  <- results_all[[1]]

  res_df    <- primary_obj$res
  ref_group <- primary_obj$groups$reference
  test_group <- primary_obj$groups$test

  if (plot) {

    suppressMessages(suppressWarnings(ma_plot()))

    if (!is.null(disp_plot))
      suppressMessages(suppressWarnings(disp_plot()))

  }

  # ===========================================================================
  # 7) Output and Attach to project
  # ===========================================================================
  if (save) {

    for (comparison_name in names(results_all)) {

      comp_obj <- results_all[[comparison_name]]

      proj <- .attach_to_project(
        proj,
        comp_obj,
        slot = "analyses",
        subtype = "comparison",
        prefix = paste0("compare_", comparison_name),
        id = if (store_model) "last" else NULL,
        log = list(
          method = method,
          design = deparse(design_formula),
          contrast = comparison_name
        )
      )
    }
  }

  # ===========================================================================
  # 8) Return
  # ===========================================================================
  # --- Prepare output ---
  if (length(results_all) == 1) {
    obj <- results_all[[1]]
  } else {
    obj <- results_all
  }

  # --- Print ---
  alpha <- 0.05

  sig <- res_df[!is.na(res_df$padj) & res_df$padj < alpha, ]
  up   <- sig[sig$log2FoldChange > 0, ]
  down <- sig[sig$log2FoldChange < 0, ]

  .print_header("RNA Differential Expression")

  .print_block("Comparison", function() {
    cat("Method:            ", toupper(method), "\n")
    cat("Design formula:    ", deparse(design_formula), "\n")
    cat("Batch included:    ",
        if (!is.null(batch_col) && batch_col %in% colnames(metadata)) "Yes" else "No",
        "\n")
  })

  .print_block("Results summary (padj < 0.05)", function() {
    cat("Total genes tested:", nrow(res_df), "\n")
    cat("Significant genes: ", nrow(sig), "\n")
    cat("Upregulated:       ", nrow(up), "\n")
    cat("Downregulated:     ", nrow(down), "\n")
  })

  cat("Stored as:         ", primary_name, "\n")

  .print_block("Top 5 genes", function() {
    print(utils::head(sig, 5))
  })

  cat("Alpha threshold:   ", alpha, "\n")

  return(invisible(proj))
}
