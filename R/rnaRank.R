#' Rank genes by mean or variance
#'
#' @description
#' Rank genes based on expression mean or variance within two groups of samples.
#' Optionally highlight specific genes and plot rankings per group.
#'
#' @param norm_data Object returned by [rna.normalize()] containing `expr_matrix` and `metadata`.
#' @param group_col Character; column in metadata defining sample groups (default `"Group"`).
#' @param metric Character; either `"mean"` or `"var"` to rank by mean expression or variance (default `"mean"`).
#' @param top_n Integer; number of top genes to retain per group (default `NULL`, keeps all genes).
#' @param species Character; `"human"`, `"mouse"` or `"zebrafish"` for gene annotation (default `"mouse"`).
#' @param highlight_gene Character vector of gene symbols or Ensembl IDs to highlight in the plot (default `NULL`).
#' @param plot Logical; whether to plot the ranking (default `TRUE`).
#'
#' @return
#' A list containing:
#' \describe{
#'   \item{ranked_genes}{Data frame with statistics and rank per gene per group.}
#'   \item{plot}{ggplot2 object of gene rankings.}
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_text facet_wrap labs theme element_text theme_minimal
#' @export

rna.rank <- function(norm_data,
                     group_col = "Group",
                     metric = c("mean", "var"),
                     top_n = NULL,
                     organism = "mouse",
                     highlight_gene = NULL,
                     plot = TRUE) {

  # --- Validate input ---
  if (!inherits(norm_data, "normalized_data")) {
    stop("Input must be a 'normalized_data' object returned by rna.normalize().")
  }

  metric <- match.arg(metric)
  expr <- norm_data$expr_matrix
  meta <- norm_data$metadata

  if (!group_col %in% colnames(meta)) {
    stop(paste("Column", group_col, "not found in metadata."))
  }

  groups <- unique(meta[[group_col]])
  if (length(groups) != 2) stop("This function supports exactly 2 groups for comparison.")

  # --- Compute statistics per group ---
  get_stats <- function(expr_sub, group_name) {
    data.frame(
      Gene = rownames(expr_sub),
      Mean = rowMeans(expr_sub, na.rm = TRUE),
      Var  = apply(expr_sub, 1, var, na.rm = TRUE),
      Group = group_name,
      stringsAsFactors = FALSE
    )
  }

  stats_list <- lapply(groups, function(g) {
    samples <- meta$Sample[meta[[group_col]] == g]
    expr_sub <- expr[, samples, drop = FALSE]
    get_stats(expr_sub, g)
  })
  stats_df <- do.call(rbind, stats_list)

  # --- Prepare biomaRt dataset & symbol attribute based on organism ---
  org_lower <- tolower(organism)
  ds <- dplyr::case_when(
    org_lower %in% c("human","h","hsapiens","homo sapiens")   ~ "hsapiens_gene_ensembl",
    org_lower %in% c("mouse","m","mmusculus","mus musculus")  ~ "mmusculus_gene_ensembl",
    org_lower %in% c("zebrafish","drerio","danio","dre","danio rerio") ~ "drerio_gene_ensembl",
    TRUE ~ "mmusculus_gene_ensembl" # fallback to mouse
  )

  if (!org_lower %in% c("human","h","hsapiens","homo sapiens",
                        "mouse","m","mmusculus","mus musculus",
                        "zebrafish","drerio","danio","dre","danio rerio")) {
    message("[rna.rank] organism not recognized → using mouse as fallback.")
  }

  # choose symbol attribute name per species (works with biomaRt common attributes)
  symbol_attr <- switch(
    ds,
    hsapiens_gene_ensembl = "hgnc_symbol",
    mmusculus_gene_ensembl = "mgi_symbol",
    drerio_gene_ensembl = "external_gene_name",
    "external_gene_name"
  )

  # --- Map Ensembl IDs to gene symbols (safe attempt) ---
  stats_df$ensembl_gene_id <- sub("\\..*$", "", stats_df$Gene)
  stats_df$gene_symbol <- NA_character_

  if (requireNamespace("biomaRt", quietly = TRUE)) {
    mart <- tryCatch({
      biomaRt::useEnsembl(biomart = "genes", dataset = ds)
    }, error = function(e) {
      # try the GRCh38 mirror if useEnsembl fails (network quirks)
      tryCatch(biomaRt::useEnsembl(biomart = "genes", dataset = ds, mirror = "useast"),
               error = function(e2) NULL)
    })

    if (!is.null(mart)) {
      # request ensembl id + the chosen symbol attribute
      attrs <- c("ensembl_gene_id", symbol_attr)
      gene_info <- tryCatch({
        biomaRt::getBM(
          attributes = attrs,
          filters = "ensembl_gene_id",
          values = unique(stats_df$ensembl_gene_id),
          mart = mart
        )
      }, error = function(e) NULL)

      if (!is.null(gene_info) && nrow(gene_info) > 0) {
        # ensure column name for symbol is standardized
        colnames(gene_info)[colnames(gene_info) == symbol_attr] <- "symbol_tmp"
        idx <- match(stats_df$ensembl_gene_id, gene_info$ensembl_gene_id)
        stats_df$gene_symbol <- gene_info$symbol_tmp[idx]
      } else {
        warning("biomaRt returned no annotations for the requested IDs; continuing with Ensembl IDs.")
      }
    } else {
      warning("Could not connect to Ensembl via biomaRt; continuing with Ensembl IDs.")
    }
  } else {
    message("Package 'biomaRt' not installed: gene symbols will not be added.")
  }

  stats_df$GeneLabel <- ifelse(is.na(stats_df$gene_symbol) | stats_df$gene_symbol == "",
                               stats_df$Gene,
                               stats_df$gene_symbol)

  # --- Filter top N genes per group ---
  metric_col <- if (metric == "mean") "Mean" else "Var"
  if (!is.null(top_n)) {
    stats_df <- do.call(rbind, lapply(split(stats_df, stats_df$Group), function(d) {
      d[order(-d[[metric_col]]), , drop = FALSE][seq_len(min(nrow(d), top_n)), , drop = FALSE]
    }))
  }

  # --- Assign rank per group ---
  stats_df <- do.call(rbind, lapply(split(stats_df, stats_df$Group), function(d) {
    d <- d[order(-d[[metric_col]]), , drop = FALSE]
    d$Rank <- seq_len(nrow(d))
    d
  }))

  # --- Plot rankings ---
  p <- ggplot2::ggplot(stats_df, ggplot2::aes(x = Rank, y = !!rlang::sym(metric_col), color = Group)) +
    ggplot2::geom_point(alpha = 0.3, size = 1) +
    ggplot2::facet_wrap(~Group, scales = "free_x", nrow = 1) +
    ggplot2::labs(
      x = "Gene rank",
      y = paste0("Expression (", metric, ")"),
      title = paste("Gene ranking by", metric, "per group")
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      axis.text.x = ggplot2::element_blank(),
      legend.position = "none"
    )

  # --- Highlight specific genes ---
  if (!is.null(highlight_gene)) {
    highlight_df <- subset(stats_df, GeneLabel %in% highlight_gene)
    if (nrow(highlight_df) > 0) {
      p <- p +
        ggplot2::geom_point(
          data = highlight_df,
          ggplot2::aes(x = Rank, y = !!rlang::sym(metric_col)),
          color = "#666666",
          size = 3
        ) +
        ggplot2::geom_text(
          data = highlight_df,
          ggplot2::aes(x = Rank, y = !!rlang::sym(metric_col), label = GeneLabel),
          vjust = -1,
          color = "black",
          size = 3
        )
    }
  }

  if (plot) print(p)

  out <- list(ranked_genes = stats_df, plot = p)
  message("[rna.rank] top rows:")
  print(utils::head(stats_df, 10))

  invisible(out)

}
