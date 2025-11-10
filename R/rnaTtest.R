#' Perform t-tests or Mann-Whitney tests for selected genes
#'
#' @description
#' Performs statistical comparison between two groups for a set of genes using t-test (default)
#' or Mann-Whitney test (utest = TRUE), and plots boxplots with significance annotation.
#'
#' @param normalized_data Object of class `normalized_data` (from `rna.normalize()`).
#' @param genes Character vector of gene symbols or Ensembl IDs to analyze.
#' @param group_col Character; metadata column defining sample groups (default `"Group"`).
#' @param species Character; `"human"` or `"mouse"` for gene annotation (default `"human"`).
#' @param utest Logical; whether to use Mann-Whitney test instead of t-test (default `FALSE`).
#' @param return_type Character; `"htest"` returns only test results, `"all"` returns results and plots (default `"all"`).
#'
#' @return
#' List of t-test (or Mann-Whitney) results and plots, depending on `return_type`.
#'
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_jitter annotate theme labs element_text
#' @importFrom dplyr left_join
#' @importFrom patchwork plot_layout
#' @export
rna.ttest <- function(normalized_data, genes,
                      group_col = "Group",
                      species = "human",
                      utest = FALSE,
                      return_type = c("htest", "all")) {
  
  # --- Check packages ---
  required_pkgs <- c("ggplot2", "dplyr", "biomaRt", "patchwork")
  missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs)) stop("Please install packages: ", paste(missing_pkgs, collapse = ", "))
  
  return_type <- match.arg(return_type)
  
  # --- Map Ensembl -> gene symbol ---
  gene_ids <- rownames(normalized_data$expr_matrix)
  ensembl_gene_id <- sub("\\..*$", "", gene_ids)
  gene_map <- data.frame(GeneID = gene_ids, Ensembl = ensembl_gene_id, stringsAsFactors = FALSE)
  
  gene_map$Symbol <- NA_character_
  if (requireNamespace("biomaRt", quietly = TRUE)) {
    mart <- tryCatch({
      biomaRt::useEnsembl(
        biomart = "genes",
        dataset = ifelse(tolower(species) %in% c("human","h","hsapiens"),
                         "hsapiens_gene_ensembl",
                         "mmusculus_gene_ensembl")
      )
    }, error = function(e) NULL)
    
    if (!is.null(mart)) {
      bm <- tryCatch({
        biomaRt::getBM(
          attributes = c("ensembl_gene_id","hgnc_symbol"),
          filters = "ensembl_gene_id",
          values = unique(ensembl_gene_id),
          mart = mart
        )
      }, error = function(e) NULL)
      
      if (!is.null(bm) && nrow(bm) > 0) {
        gene_map$Symbol <- bm$hgnc_symbol[match(gene_map$Ensembl, bm$ensembl_gene_id)]
      } else {
        warning("biomaRt returned empty; symbols unavailable.")
      }
    } else {
      warning("Could not connect to Ensembl; symbols unavailable.")
    }
  } else {
    message("Package 'biomaRt' not installed; symbols will not be added.")
  }
  
  # --- Initialize results ---
  results_list <- list()
  plots_list <- list()
  
  # --- Loop over genes ---
  for (gene in genes) {
    idx <- which(gene_map$GeneID == gene | gene_map$Symbol == gene)
    if (length(idx) == 0) {
      warning(paste0("Gene '", gene, "' not found; skipping."))
      next
    }
    
    gene_use <- gene_map$GeneID[idx[1]]
    gene_label <- ifelse(is.na(gene_map$Symbol[idx[1]]) | gene_map$Symbol[idx[1]] == "",
                         gene_use, gene_map$Symbol[idx[1]])
    
    # --- Extract expression ---
    expr_mat <- normalized_data$expr_matrix
    df_long <- data.frame(Sample = colnames(expr_mat),
                          Value = as.numeric(expr_mat[gene_use, ]))
    df_long <- dplyr::left_join(df_long, normalized_data$metadata, by = "Sample")
    
    # --- Check two groups ---
    groups <- unique(df_long[[group_col]])
    if (length(groups) != 2) {
      warning(paste0("Gene '", gene_label, "' does not have exactly 2 groups; skipping."))
      next
    }
    
    g1 <- df_long$Value[df_long[[group_col]] == groups[1]]
    g2 <- df_long$Value[df_long[[group_col]] == groups[2]]
    g1 <- g1[is.finite(g1)]
    g2 <- g2[is.finite(g2)]
    if (length(g1) < 2 || length(g2) < 2) {
      warning(paste0("Gene '", gene_label, "' has insufficient values; skipping."))
      next
    }
    
    # --- Statistical test ---
    if (utest) {
      test_name <- "Mann-Whitney"
      warning("Mann-Whitney on normalized data is not recommended.")
      res <- suppressWarnings(stats::wilcox.test(g1, g2, exact = FALSE))
    } else {
      test_name <- "t-test"
      res <- stats::t.test(g1, g2)
    }
    
    # --- P-value annotation ---
    pval <- res$p.value
    signif_label <- if (pval < 0.001) "***" else if (pval < 0.01) "**" else if (pval < 0.05) "*" else "ns"
    
    # --- Boxplot ---
    y_max <- max(df_long$Value, na.rm = TRUE)
    y_pos <- y_max * 1.02
    p <- ggplot2::ggplot(df_long, ggplot2::aes(x = .data[[group_col]], y = Value, fill = .data[[group_col]])) +
      ggplot2::geom_boxplot(alpha = 0.75, outlier.shape = NA) +
      ggplot2::geom_jitter(width = 0.1, alpha = 0.65, color = "black") +
      ggplot2::annotate("text", x = 1.5, y = y_pos, label = signif_label, size = 4) +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(legend.position = "none",
                     axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                     plot.title = ggplot2::element_text(hjust = 0.5)) +
      ggplot2::labs(title = gene_label,
                    subtitle = paste0(test_name, " p-value: ", signif(pval, 3)),
                    x = "", y = "Normalized expression")
    
    # --- Store results ---
    results_list[[gene_label]] <- res
    plots_list[[gene_label]] <- p
  }
  
  # --- Combine plots ---
  if (length(plots_list) > 0) {
    combined_plot <- Reduce(`+`, plots_list) + patchwork::plot_layout(ncol = 2)
    print(combined_plot)
  }
  
  # --- Return ---
  if (return_type == "all") {
    invisible(list(tests = results_list, plots = plots_list))
  } else {
    invisible(results_list)
  }
}
