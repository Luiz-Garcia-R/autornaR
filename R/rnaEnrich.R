#' Gene Ontology (GO) enrichment analysis for differential expression results
#'
#' @description
#' Performs Gene Ontology (GO) enrichment analysis on differentially expressed
#' genes derived from \code{rna.compare()}, using the \pkg{clusterProfiler}
#' framework. Enrichment is computed separately for up- and down-regulated genes,
#' allowing directional biological interpretation.
#'
#' The function automatically retrieves the selected comparison from the active
#' \code{rna_project}, applies filtering based on adjusted p-values and
#' log2 fold-change thresholds, maps gene identifiers, and computes enrichment
#' across one or more ontology domains (BP, MF, CC).
#'
#' @param project \code{rna_project} object created by \code{rna.project()}.
#' @param use_comparison Character or numeric. Identifier of the comparison
#'   stored in \code{rna_project$analyses$comparison}. If \code{NULL},
#'   the most recent comparison is used.
#' @param plot_style Character. Wheter use \code{barplot} or \code{dotplot}.
#' @param ont Character vector. Ontology domains to test:
#'   \itemize{
#'     \item \code{"BP"}: Biological Process
#'     \item \code{"MF"}: Molecular Function
#'     \item \code{"CC"}: Cellular Component
#'   }
#' @param enrich_p_cutoff Numeric. P-value cutoff used internally by
#'   \code{clusterProfiler::enrichGO()} to filter enriched GO terms
#'   (default: \code{0.1}).
#' @param padj_cutoff Numeric. Adjusted p-value threshold used to filter
#'   differentially expressed genes prior to enrichment (default: \code{0.05}).
#' @param log2fc_cutoff Numeric. Absolute log2 fold-change threshold used
#'   to define up- and down-regulated genes (default: \code{0.5}).
#' @param plot Logical. Whether to display dot plots of enriched GO terms
#'   (default: \code{TRUE}).
#' @param top_terms Integer. Number of top enriched terms to display in plots
#'   (default: \code{10}).
#' @param save Logical. Whether to store results in the active
#'   \code{rna_project} (default: \code{TRUE}).
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Retrieves differential expression results from \code{rna_project}
#'   \item Filters genes using \code{padj_cutoff} and \code{log2fc_cutoff}
#'   \item Splits genes into up- and down-regulated sets
#'   \item Performs GO enrichment separately for each set and ontology
#'   \item Optionally visualizes results using dot plots
#' }
#'
#' Gene identifiers are automatically interpreted using the project metadata
#' (\code{gene_id_type}) and passed to \code{clusterProfiler::enrichGO()}.
#'
#' Supported annotation databases:
#' \itemize{
#'   \item Human: \code{org.Hs.eg.db}
#'   \item Mouse: \code{org.Mm.eg.db}
#'   \item Zebrafish: \code{org.Dr.eg.db}
#' }
#'
#' \strong{Important:}
#' \itemize{
#'   \item Gene IDs must be cleaned (no Ensembl version suffixes)
#'   \item Enrichment depends strongly on filtering thresholds
#'   \item Up- and down-regulated enrichments may reveal distinct biology
#' }
#'
#' @return
#' An object of class \code{"enrich_result"} containing:
#' \describe{
#'   \item{comparison}{Comparison identifier used}
#'   \item{contrast}{Group contrast (reference vs test)}
#'   \item{organism}{Organism used for annotation}
#'   \item{gene_id_type}{Gene identifier type used in enrichment}
#'   \item{up}{Named list of enrichment results for up-regulated genes (per ontology)}
#'   \item{down}{Named list of enrichment results for down-regulated genes (per ontology)}
#'   \item{params}{List of parameters used in the analysis}
#' }
#'
#' @section Stored output:
#' When \code{save = TRUE}, results are stored in:
#' \code{rna_project$analyses$enrich}, allowing reproducibility and reuse.
#'
#' @section Visualization:
#' When \code{plot = TRUE}, dot plots are generated using
#' \code{enrichplot::dotplot()}, showing the top enriched GO terms per ontology
#' and direction (up/down).
#'
#' @examples
#' \dontrun{
#' # Run enrichment using last comparison
#' rna.enrich(project = my_project)
#'
#' # Run enrichment with stricter thresholds
#' rna.enrich(my_project,
#'            padj_cutoff = 0.01,
#'            log2fc_cutoff = 1)
#'
#' # Run enrichment for molecular function
#' rna.enrich(my_project,
#'            ont = "MF")
#' }
#'
#' @importFrom ggplot2 ggtitle theme_minimal
#'
#' @export

rna.enrich <- function(project,
                       use_comparison = NULL,
                       plot_style = c("dotplot", "barplot"),
                       ont = c("BP", "MF", "CC"),
                       enrich_p_cutoff = 0.1,
                       padj_cutoff = 0.05,
                       log2fc_cutoff = 0.5,
                       plot = TRUE,
                       top_terms = 10,
                       save = TRUE
) {

  plot_style <- match.arg(plot_style)

  # ---------------------------
  # 0) Basic checks
  # ---------------------------
  ont <- match.arg(ont, several.ok = TRUE)

  .check_dependencies("clusterProfiler")

  if (plot) {
    .check_dependencies("ggplot2")
  }

  if (plot && plot_style == "dotplot") {
    .check_dependencies("enrichplot")
  }

  # ---------------------------
  # 1) Get active project
  # ---------------------------
  proj <- project

  organism <- .get_organism(proj)
  gene_id_type <- .get_gene_id_type(proj)
  comps <- .get_comp(proj)

  # ---------------------------
  # 2) Validate input
  # ---------------------------
  if (is.null(use_comparison)) {
    use_comparison <- tail(setdiff(names(comps), "last"), 1)
    msg <- "[rna.enrich] Using last comparison: "
  } else {
    if (is.numeric(use_comparison)) {
      ids <- setdiff(names(comps), "last")
      use_comparison <- ids[use_comparison]
    }
    msg <- "[rna.enrich] Using stored comparison: "
  }

  message(msg, use_comparison)

  comp_obj <- .get_comp_obj(proj, use_comparison)

  # ---------------------------
  # 3) Validate organism
  # ---------------------------
  message("[rna.enrich] Using organism from project: ", organism)

  org_pkg <- switch(
    organism,
    human = "org.Hs.eg.db",
    mouse = "org.Mm.eg.db",
    zebrafish = "org.Dr.eg.db",
    stop("`organism` must be 'human', 'mouse' or 'zebrafish'.")
  )

  .check_dependencies(c("clusterProfiler", "enrichplot", "ggplot2", org_pkg))

  # --- Select organism database ---
  OrgDb <- switch(
    organism,
    human = org.Hs.eg.db::org.Hs.eg.db,
    mouse = org.Mm.eg.db::org.Mm.eg.db,
    zebrafish = org.Dr.eg.db::org.Dr.eg.db
  )

  # ---------------------------
  # 3) Validate input
  # ---------------------------
  res <- as.data.frame(comp_obj$res)

  contrast <- c(
    comp_obj$groups$reference,
    comp_obj$groups$test
  )

  # --- Gene ID sanity check ---
  gene_ids <- rownames(res)

  if (any(grepl("\\.", gene_ids))) {
    stop(
      "[rna_enrich] Gene IDs contain version suffixes (e.g. ENSG...\\.6).\n",
      "This function expects cleaned Ensembl IDs.\n",
      "Please run rna.normalize(clean_gene_versions = TRUE) upstream."
    )
  }

  # ---------------------------
  # 4) Filter up- and down-regulated genes
  # ---------------------------
  up_genes <- res[!is.na(res$padj) &
                    res$padj < padj_cutoff &
                    res$log2FoldChange > log2fc_cutoff, ]

  down_genes <- res[!is.na(res$padj) &
                      res$padj < padj_cutoff &
                      res$log2FoldChange < -log2fc_cutoff, ]

  .run_direction <- function(gene_table, direction_label) {

    if (nrow(gene_table) == 0) {
      return(NULL)
    }

    res_list <- lapply(ont, function(o) {
      .run_go_enrichment(
        gene_ids = rownames(gene_table),
        OrgDb = OrgDb,
        from_type = gene_id_type,
        ont = o,
        p_cutoff = enrich_p_cutoff
      )
    })

    names(res_list) <- ont

    # Dotplot
    if (plot) {
      for (o in names(res_list)) {
        ego <- res_list[[o]]

        if (is.null(ego)) next

        df <- as.data.frame(ego)

        if (nrow(df) == 0) next

          if (plot_style == "dotplot") {

            print(
              enrichplot::dotplot(ego, showCategory = top_terms) +
                ggplot2::ggtitle(
                  paste0("GO ", o, " (", direction_label, " genes)")
                ) +
                ggplot2::theme_minimal()
            )

            # Barplot
          } else if (plot_style == "barplot") {

            df <- as.data.frame(ego)

            df <- df[order(df$p.adjust), ]
            df <- head(df, top_terms)

            df$Description <- factor(df$Description, levels = rev(df$Description))

            p <- ggplot2::ggplot(
              df,
              ggplot2::aes(x = -log10(p.adjust),
                           y = .data$Description,
                           fill = p.adjust)
            ) +

              ggplot2::geom_col() +

              ggplot2::scale_fill_gradient(
                low = "#DD6765",
                high = "#327EBA"
              ) +
              ggplot2::labs(
                title = paste0("GO ", o, " (", direction_label, " genes)"),
                x = expression(-log[10]("adjusted p-value")),
                y = NULL
              ) +

              ggplot2::theme_minimal() +

              ggplot2::theme(
                axis.text.y = ggplot2::element_text(size = 10)
              )

            print(p)
          }
        }
      }

    res_list
  }

  if (nrow(up_genes) == 0 && nrow(down_genes) == 0) {
    message("[rna_enrich] No significant genes found with the given thresholds.")
    return(invisible(NULL))
  }


  # ---------------------------
  # 5) enrichment + plotting
  # ---------------------------
  ego_up   <- .run_direction(up_genes, "up-regulated")
  ego_down <- .run_direction(down_genes, "down-regulated")

  # ---------------------------
  # 6) Output
  # ---------------------------

  params <- list(
    timestamp = Sys.time(),
    padj_cutoff = padj_cutoff,
    log2fc_cutoff = log2fc_cutoff,
    ont = ont,
    top_terms = top_terms
  )

  obj <- list(
    comparison = use_comparison,
    contrast = contrast,
    organism = proj$input$imp_data$organism,
    gene_id_type = gene_id_type,
    up = ego_up,
    down = ego_down,
    params = params
  )

  class(obj) <- "enrich_result"

  # ---------------------------
  # 7) Attach to project
  # ---------------------------
  if (save) {

    proj <- .attach_to_project(
      proj,
      obj,
      slot = "analyses",
      subtype = "enrichment",
      prefix = "enrich",
      log = list(
        comparison = use_comparison,
        organism = proj$input$imp_data$organism,
        padj_cutoff = padj_cutoff,
        log2fc_cutoff = log2fc_cutoff,
        ont = ont
      )
    )
  }

  return(invisible(proj))
}

# ---------------------------
# 8) Print S3
# ---------------------------
#' Print method for \code{enrich_result} objects
#'
#' Displays a concise summary of Gene Ontology (GO) enrichment results,
#' including comparison information, organism, and the number of
#' significantly enriched terms for up- and down-regulated genes.
#'
#' This method is automatically called when an \code{enrich_result}
#' object is printed in the console.
#'
#' @param x An object of class \code{enrich_result}.
#' @param ... Additional arguments (ignored).
#'
#' @return
#' Invisibly returns \code{x}.
#'
#' @method print enrich_result
#' @export
#' @rdname enrich_result

print.enrich_result <- function(x, ...) {

  .print_header("GO Enrichment Result")

  .print_block("Overview", function() {
    cat("Class: enrich_result\n")
    cat("Comparison: ", x$comparison, "\n", sep = "")
    cat("Contrast: ", paste(x$contrast, collapse = " vs "), "\n", sep = "")
    cat("Organism: ", x$organism, "\n", sep = "")

    up_n <- if (!is.null(x$up)) {
      sum(sapply(x$up, function(e) {
        if (is.null(e)) 0 else nrow(as.data.frame(e))
      }))
    } else 0

    down_n <- if (!is.null(x$down)) {
      sum(sapply(x$down, function(e) {
        if (is.null(e)) 0 else nrow(as.data.frame(e))
      }))
    } else 0

    cat("Significant GO terms (up): ", up_n, "\n", sep = "")
    cat("Significant GO terms (down): ", down_n, "\n", sep = "")
  })

  invisible(x)
}

#' Summary method for \code{enrich_result} objects
#'
#' Returns a table of enriched Gene Ontology (GO) terms for a given
#' direction (up- or down-regulated genes).
#'
#' @param object An object of class \code{enrich_result}.
#' @param direction Character. Which gene set to summarize:
#' \code{"up"} or \code{"down"}.
#' @param ont Optional character specifying the ontology to return
#' (e.g. \code{"BP"}, \code{"MF"}, \code{"CC"}). If \code{NULL},
#' returns results for all ontologies.
#' @param ... Additional arguments (ignored).
#'
#' @return
#' A data frame with enriched GO terms for the selected direction
#' and ontology. If \code{ont = NULL}, returns a list of data frames
#' (one per ontology).
#'
#' @method summary enrich_result
#' @export
#' @rdname enrich_result

summary.enrich_result <- function(object,
                                  direction = c("up", "down"),
                                  ont = NULL,
                                  ...) {

  direction <- match.arg(direction)
  ego_list <- object[[direction]]

  if (is.null(ego_list)) {
    message("No enrichment results for ", direction)
    return(invisible(NULL))
  }

  if (is.null(ont)) {
    return(ego_list)
  }

  if (!ont %in% names(ego_list)) {
    stop("Ontology not found. Available: ",
         paste(names(ego_list), collapse = ", "))
  }

  as.data.frame(ego_list[[ont]])
}
