#' ROC analysis for selected genes
#'
#' @description
#' Compute ROC curve and AUC for selected genes in a normalized RNA-seq dataset.
#' Uses logistic regression for ≤3 genes or PCA1 for more than 3 genes.
#'
#' @param normalized_data Object of class `normalized_data` (from `rna.normalize()`).
#' @param genes Character vector of gene symbols or Ensembl IDs to analyze, or `"all"` for all genes (default `NULL`).
#' @param group_col Character; metadata column defining sample groups (default `"Group"`).
#' @param species Character; `"human"`, `"mouse"` or "`"zebrafish"` for gene annotation (default `"mouse"`).
#' @param plot Logical; whether to plot the ROC curve (default `TRUE`).
#'
#' @return
#' A list containing:
#' \describe{
#'   \item{roc_object}{pROC ROC object.}
#'   \item{auc}{AUC value.}
#'   \item{genes}{Genes used for the analysis (symbols if available).}
#'   \item{method}{Method used: `"logistic"` or `"pca1"`.}
#' }
#'
#' @importFrom pROC roc auc
#' @export

rna.roc <- function(normalized_data,
                    genes = NULL,
                    group_col = "Group",
                    species = "mouse",
                    plot = TRUE) {

  # --- Check required packages ---
  required_pkgs <- c("pROC", "ggplot2", "dplyr")
  missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs)) stop("Please install packages: ", paste(missing_pkgs, collapse = ", "))

  # --- Prepare Ensembl -> gene symbol mapping ---
  gene_ids <- rownames(normalized_data$expr_matrix)
  ensembl_gene_id <- sub("\\..*$", "", gene_ids)
  gene_map <- data.frame(GeneID = gene_ids, Ensembl = ensembl_gene_id, stringsAsFactors = FALSE)

  # --- Try biomaRt annotation ---
  gene_map$gene_symbol <- NA_character_

  if (requireNamespace("biomaRt", quietly = TRUE)) {

    ds <- dplyr::case_when(
      tolower(species) %in% c("mouse","m","mmusculus") ~ "mmusculus_gene_ensembl",
      tolower(species) %in% c("human","h","hsapiens") ~ "hsapiens_gene_ensembl",
      tolower(species) %in% c("zebrafish","drerio","danio","dre") ~ "drerio_gene_ensembl",
      TRUE ~ {
        message("[rna.roc] species não reconhecido → fallback para mouse (default).")
        "mmusculus_gene_ensembl"
      }
    )

    mart <- tryCatch(
      biomaRt::useEnsembl(biomart = "genes", dataset = ds),
      error = function(e) NULL
    )

    if (!is.null(mart)) {

      # atributo gene symbol ENSEMBL universal → "external_gene_name"
      bm <- tryCatch({
        biomaRt::getBM(
          attributes = c("ensembl_gene_id","external_gene_name"),
          filters = "ensembl_gene_id",
          values = unique(ensembl_gene_id),
          mart = mart
        )
      }, error = function(e) NULL)

      if (!is.null(bm) && nrow(bm) > 0) {
        gene_map$gene_symbol <- bm$external_gene_name[match(gene_map$Ensembl, bm$ensembl_gene_id)]
      } else {
        warning("biomaRt returned empty; gene symbols unavailable.")
      }

    } else {
      warning("Could not connect to Ensembl; gene symbols unavailable.")
    }

  } else {
    message("Package 'biomaRt' not installed; gene symbols will not be added.")
  }

  # --- Determine which genes to use ---
  if (is.null(genes) || (length(genes) == 1 && tolower(genes) == "all")) {
    genes_use <- gene_map$GeneID
    gene_labels <- ifelse(is.na(gene_map$gene_symbol) | gene_map$gene_symbol == "",
                          genes_use, gene_map$gene_symbol)
    message("Using all genes (", length(genes_use), " genes).")
    show_gene_names <- FALSE
  } else {
    valid_idx <- which(gene_map$GeneID %in% genes | gene_map$gene_symbol %in% genes)
    if (length(valid_idx) == 0) stop("No matching genes found.")
    genes_use <- gene_map$GeneID[valid_idx]
    gene_labels <- ifelse(is.na(gene_map$gene_symbol[valid_idx]) | gene_map$gene_symbol[valid_idx] == "",
                          genes_use, gene_map$gene_symbol[valid_idx])
    show_gene_names <- TRUE
  }

  # --- Prepare data ---
  expr_mat <- normalized_data$expr_matrix[genes_use, , drop = FALSE]
  df <- as.data.frame(t(expr_mat))
  colnames(df) <- gene_labels
  df[[group_col]] <- as.factor(normalized_data$metadata[[group_col]])

  # --- Select method ---
  method <- if (length(genes_use) <= 3) "logistic" else "pca1"

  # --- Compute predictions ---
  if (method == "logistic") {
    formula_str <- paste(group_col, "~", paste(gene_labels, collapse = " + "))
    logit_model <- stats::glm(as.formula(formula_str), data = df, family = stats::binomial)
    pred <- stats::predict(logit_model, type = "response")
  } else {
    pca <- stats::prcomp(df[, gene_labels, drop = FALSE], scale. = TRUE)
    pred <- pca$x[, 1]  # PC1
  }

  # --- ROC and AUC ---
  classe <- df[[group_col]]
  roc_obj <- pROC::roc(classe, pred, levels = levels(classe), direction = "<")
  auc_val <- pROC::auc(roc_obj)

  # --- Plot ---
  if (plot) {
    plot(1 - roc_obj$specificities, roc_obj$sensitivities,
         type = "l", col = "#1f77b4", lwd = 3,
         main = if(method == "logistic") "ROC - Logistic" else "ROC - PCA1",
         xlab = "False Positive Rate",
         ylab = "True Positive Rate",
         xlim = c(0,1), ylim = c(0,1), bty = "n", las = 1)
    abline(a = 0, b = 1, col = "#d3d3d3", lty = 2, lwd = 2)
    text(0.8, 0.05, paste("AUC =", round(auc_val, 3)), col = "#1f77b4", cex = 1.2, font = 2)

    if (show_gene_names && length(genes_use) <= 5) {
      mtext(side = 3, line = 0.5, at = mean(par("usr")[1:2]),
            text = paste("Genes:", paste(gene_labels, collapse = " + ")), cex = 1)
    } else if (show_gene_names) {
      message("Number of genes used: ", length(genes_use))
    }

    grid(col = "lightgray", lty = "dotted")
  }

  # --- Return result ---
  invisible(list(
    roc_object = roc_obj,
    auc = auc_val,
    genes = gene_labels,
    method = method
  ))
}
