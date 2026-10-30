# ============================================================================ #
#                     ### Auxiliary general functions ###                      #
# ============================================================================ #

# =============================================================================
# 1. Global cache environment
# =============================================================================
#' @keywords internal

.autornar_cache <- new.env(parent = emptyenv())
.autornar_cache$gene_maps <- list()

# =============================================================================
# 2. Auxiliary print functions
# =============================================================================
#' @keywords internal

.print_header <- function(title) {
  cat("\n")
  cat(strrep("=", 50), "\n")
  cat(title, "\n")
  cat(strrep("=", 50), "\n")
}

#' @keywords internal
.print_block <- function(title, content, width = 40) {
  cat("\n", title, "\n", sep = "")
  cat(strrep("-", width), "\n")
  content()
  cat(strrep("-", width), "\n")
}

# =============================================================================
# 3. Attach to project
# =============================================================================
#' @keywords internal

.attach_to_project <- function(proj,
                               obj,
                               slot,
                               subtype,
                               log = list()) {

  slot <- as.character(slot)

  if (!slot %in% c("input", "data", "analyses", "qc")) {
    stop("Invalid slot.")
  }

  # Ensure base structure
  if (is.null(proj[[slot]])) {
    proj[[slot]] <- list()
  }

  # Replace previous object
  proj[[slot]][[subtype]] <- obj

  # ===============================================
  # 3.1. Logs
  # ===============================================

  if (is.null(proj$logs)) {
    proj$logs <- list()
  }

  proj$logs[[subtype]] <- c(
    list(timestamp = Sys.time()),
    log
  )

  proj
}

# =============================================================================
# 4. Set seed
# =============================================================================
#' @keywords internal

.set_seed <- function(seed) {
  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) {
      get(".Random.seed", envir = .GlobalEnv)
    } else {
      NULL
    }
    set.seed(seed)
    return(old_seed)
  }
  return(NULL)
}

# =============================================================================
# 5. Reset seed
# =============================================================================
#' @keywords internal

.reset_seed <- function(old_seed) {
  if (!is.null(old_seed)) {
    .Random.seed <<- old_seed
  }
}

# =============================================================================
# 6. Get Biomart dataset
# =============================================================================
#' @keywords internal

.get_biomart_dataset <- function(organism) {

  organism <- tolower(organism)

  switch(organism,
         "human" = "hsapiens_gene_ensembl",
         "mouse" = "mmusculus_gene_ensembl",
         "zebrafish" = "drerio_gene_ensembl",
         stop("Unsupported organism: ", organism)
  )
}

# =============================================================================
# 7. Order contrast levels
# =============================================================================
#' @keywords internal

.order_contrast_levels <- function(levels_vec) {

  stopifnot(length(levels_vec) == 2)

  baseline_patterns <- c(
    "^control$", "^ctrl$", "^baseline$", "^vehicle$",
    "^untreated$", "^mock$", "^hd$", "^wt$"
  )

  lv_lower <- tolower(levels_vec)

  is_baseline <- sapply(
    lv_lower,
    function(x) any(grepl(paste(baseline_patterns, collapse = "|"), x))
  )

  if (sum(is_baseline) == 1) {
    return(levels_vec[c(which(is_baseline), which(!is_baseline))])
  }

  # deterministic fallback
  sort(levels_vec)
}

# =============================================================================
# 8. Check dependencies
# =============================================================================
#' @keywords internal

.check_dependencies <- function(pkgs, bioc = FALSE) {

  missing <- pkgs[
    !vapply(pkgs, requireNamespace,
            logical(1), quietly = TRUE)
  ]

  if (length(missing) > 0) {

    installer <- if (bioc) {
      paste0(
        "BiocManager::install(c(",
        paste(sprintf('"%s"', missing), collapse = ", "),
        "))"
      )
    } else {
      paste0(
        "install.packages(c(",
        paste(sprintf('"%s"', missing), collapse = ", "),
        "))"
      )
    }

    stop(
      paste0(
        "Missing required packages: ",
        paste(missing, collapse = ", "),
        "\n\nInstall using:\n",
        installer
      ),
      call. = FALSE
    )
  }
}

# =============================================================================
# 9. Smart pathway naming
# =============================================================================
#' @keywords internal

.smart_pathway_name <- function(x, max_words = 6) {

  vapply(
    x,
    function(xx) {

      if (is.na(xx) || !nzchar(xx)) {
        return("Unknown_Pathway")
      }

      # 1. Remove prefixes
      xx <- gsub("^(REACTOME_|KEGG_|HALLMARK_|GOBP_|GO_)", "", xx)

      # 2. Normalize
      xx <- gsub("_", " ", xx)
      xx <- tolower(xx)

      # 3. Smart substitution (semantic compression)
      xx <- gsub("positive regulation of", "+reg.", xx)
      xx <- gsub("negative regulation of", "-reg.", xx)
      xx <- gsub("regulation of", "reg.", xx)

      xx <- gsub("immune response", "immune resp.", xx)
      xx <- gsub("inflammatory response", "inflammation", xx)

      xx <- gsub("b cell", "b-cell", xx)
      xx <- gsub("t cell", "t-cell", xx)

      xx <- gsub("mediated", "", xx)

      words <- unlist(strsplit(xx, "\\s+"))
      words <- words[nzchar(words)]

      if (length(words) > 1) {
        words <- words[c(TRUE, words[-1] != words[-length(words)])]
      }

      priority <- c(
        "immune","inflammation","interferon","cytokine",
        "apoptosis","cell","dna","repair",
        "metabolism","glycolysis","lipid",
        "signaling","mtor","nfkb","jak","stat"
      )

      idx_priority <- which(words %in% priority)
      if (length(idx_priority) > 0) {
        words <- c(words[idx_priority], words[-idx_priority])
      }

      words <- head(words, max_words)

      out <- paste(words, collapse = " ")
      out <- tools::toTitleCase(out)

      out <- gsub("\\s+", " ", out)
      out <- trimws(out)

      if (!nzchar(out)) {
        out <- "Unknown_Pathway"
      }

      out
    },
    FUN.VALUE = character(1)
  )
}

# =============================================================================
# 10. Go enrichment
# =============================================================================
#' @keywords internal

.run_go_enrichment <- function(
    gene_ids,
    OrgDb,
    from_type = "ENSEMBL",
    ont = "BP",
    p_cutoff = 0.1
) {

  if (length(gene_ids) == 0) return(NULL)

  keytype_map <- c(
    ENSEMBL = "ENSEMBL",
    ENTREZ  = "ENTREZID",
    SYMBOL  = "SYMBOL"
  )

  if (!from_type %in% names(keytype_map)) {
    stop("Unsupported gene_id_type: ", from_type)
  }

  if (from_type != "ENTREZ") {

    mapped <- clusterProfiler::bitr(
      gene_ids,
      fromType = keytype_map[[from_type]],
      toType = "ENTREZID",
      OrgDb = OrgDb
    )

    if (nrow(mapped) == 0) return(NULL)

    gene_ids <- mapped$ENTREZID
  }

  ego <- clusterProfiler::enrichGO(
    gene = gene_ids,
    OrgDb = OrgDb,
    keyType = "ENTREZID",
    ont = ont,
    pAdjustMethod = "BH",
    qvalueCutoff = p_cutoff,
    readable = TRUE
  )

  return(ego)
}

# =============================================================================
# 11. Heatmap gene selection
# =============================================================================
#' @keywords internal

.select_heatmap_genes <- function(expr_mat,
                                  res_df = NULL,
                                  genes = NULL,
                                  use_variance = FALSE,
                                  top_n = 100) {

  if (!is.null(res_df)) {

    if ("padj" %in% colnames(res_df) && sum(!is.na(res_df$padj)) > 0) {
      score_col <- "padj"
    } else if ("pvalue" %in% colnames(res_df)) {
      score_col <- "pvalue"
    } else {
      stop("Comparison results must contain 'padj' or 'pvalue'.")
    }

    genes_sel <- rownames(
      head(res_df[order(res_df[[score_col]]), , drop = FALSE], top_n)
    )

  } else if (!is.null(genes)) {

    genes_sel <- genes

  } else if (use_variance) {

    var_genes <- apply(expr_mat, 1, var, na.rm = TRUE)
    genes_sel <- names(sort(var_genes, decreasing = TRUE))[1:min(top_n, length(var_genes))]

  } else {

    stop("No gene selection method provided.")
  }

  # Check if genes exist in matrix
  genes_sel <- genes_sel[genes_sel %in% rownames(expr_mat)]

  if (length(genes_sel) == 0)
    stop("No selected genes found in expression matrix.")

  missing <- setdiff(genes_sel, rownames(expr_mat))

  if (length(missing) > 0)
    message("[rna.heatmap] ", length(missing), " genes not found in expression matrix.")

  genes_sel
}

# =============================================================================
# 12. Convert gene ids
# =============================================================================
#' @keywords internal

.convert_gene_ids <- function(genes,
                              from = "SYMBOL",
                              to = "ENSEMBL",
                              organism = "human") {

  pkgs <- c("AnnotationDbi")

  if (organism == "human") pkgs <- c(pkgs, "org.Hs.eg.db")
  if (organism == "mouse") pkgs <- c(pkgs, "org.Mm.eg.db")
  if (organism == "zebrafish") pkgs <- c(pkgs, "org.Dr.eg.db")

  .check_dependencies(pkgs)
  db <- switch(
    organism,
    human = org.Hs.eg.db::org.Hs.eg.db,
    mouse = org.Mm.eg.db::org.Mm.eg.db,
    zebrafish = org.Dr.eg.db::org.Dr.eg.db,
    stop("Unsupported organism")
  )
  keytype <- switch(
    from,
    SYMBOL = "SYMBOL",
    ENSEMBL = "ENSEMBL",
    stop("Unsupported 'from' type")
  )
  column <- switch(
    to,
    SYMBOL = "SYMBOL",
    ENSEMBL = "ENSEMBL",
    stop("Unsupported 'to' type")
  )
  mapped <- AnnotationDbi::mapIds(
    db,
    keys = genes,
    keytype = keytype,
    column = column,
    multiVals = "first"
  )

  mapped <- unname(mapped)
  mapped[is.na(mapped)] <- genes[is.na(mapped)]

  return(mapped)
}


# =============================================================================
# 13. Map gene symbols
# =============================================================================
#' @keywords internal

.map_gene_annotation <- function(genes, organism = "human") {

  pkgs <- c("AnnotationDbi")

  if (organism == "human") pkgs <- c(pkgs, "org.Hs.eg.db")
  if (organism == "mouse") pkgs <- c(pkgs, "org.Mm.eg.db")
  if (organism == "zebrafish") pkgs <- c(pkgs, "org.Dr.eg.db")

  .check_dependencies(pkgs)

  db <- switch(
    organism,
    human = org.Hs.eg.db::org.Hs.eg.db,
    mouse = org.Mm.eg.db::org.Mm.eg.db,
    zebrafish = org.Dr.eg.db::org.Dr.eg.db
  )

  genes_clean <- sub("\\..*$", "", genes)

  map <- AnnotationDbi::select(
    db,
    keys = genes_clean,
    keytype = "ENSEMBL",
    columns = c("SYMBOL", "ENTREZID")
  )

  map <- map[!duplicated(map$ENSEMBL), ]

  gene_map <- data.frame(
    gene_id = genes_clean,
    symbol = map$SYMBOL[match(genes_clean, map$ENSEMBL)],
    entrez = map$ENTREZID[match(genes_clean, map$ENSEMBL)],
    stringsAsFactors = FALSE
  )

  return(gene_map)
}

# =============================================================================
# 14. Convert gene symbols
# =============================================================================
#' @keywords internal

.convert_expr_to_symbols <- function(expr, proj) {

  gene_map <- .get_gene_annotation(proj)

  gene_symbols <- gene_map$symbol[
    match(
      rownames(expr),
      gene_map$gene_id
    )
  ]

  gene_symbols[
    is.na(gene_symbols) |
      gene_symbols == ""
  ] <- rownames(expr)[
    is.na(gene_symbols) |
      gene_symbols == ""
  ]

  valid <- !is.na(gene_symbols) &
    gene_symbols != ""

  expr <- expr[valid, , drop = FALSE]

  rownames(expr) <- gene_symbols[valid]

  expr <- rowsum(
    expr,
    group = rownames(expr)
  )

  as.matrix(expr)
}

# =============================================================================
# 15. Get gene annotation
# =============================================================================
#' @keywords internal

.align_gene_annotation <- function(gene_annotation, expr_mat) {

  if (is.null(gene_annotation)) {
    stop("Gene annotation not found.")
  }

  genes_expr <- rownames(expr_mat)

  if (is.null(genes_expr)) {
    stop("Expression matrix must have rownames.")
  }

  genes_expr_clean <- sub("\\..*$", "", genes_expr)
  gene_annotation$gene_id <- sub("\\..*$", "", gene_annotation$gene_id)

  idx <- match(genes_expr_clean, gene_annotation$gene_id)

  aligned <- gene_annotation[idx, , drop = FALSE]
  aligned$gene_id <- genes_expr_clean

  return(aligned)
}

# =============================================================================
# 16. Detect organism
# =============================================================================
#' @keywords internal
#' @importFrom stats na.omit

.detect_organism <- function(ids) {
  ids <- na.omit(ids)
  ids <- ids[1:min(100, length(ids))]
  ids <- sub("\\..*$", "", ids)

  if (mean(grepl("^ENSG", ids)) > 0.8) return("human")
  if (mean(grepl("^ENSMUSG", ids)) > 0.8) return("mouse")
  if (mean(grepl("^ENSDARG", ids)) > 0.8) return("zebrafish")

  return("unknown")
}

# =============================================================================
# 17. Normalize sample names
# =============================================================================
#' @keywords internal

.normalize_sample_names <- function(x) {

  x <- trimws(x)

  x <- gsub(
    "\\.(bam|counts|txt|csv)$",
    "",
    x,
    ignore.case = TRUE
  )

  x
}

# =============================================================================
# 18. Get biomart gene info
# =============================================================================
#' @keywords internal

.get_biomart_gene_info <- function(gene_ids, organism) {

  dataset <- switch(
    tolower(organism),
    "human"     = "hsapiens_gene_ensembl",
    "mouse"     = "mmusculus_gene_ensembl",
    "zebrafish" = "drerio_gene_ensembl",
    "mmusculus_gene_ensembl"
  )

  ensembl <- biomaRt::useEnsembl(
    biomart = "genes",
    dataset = dataset
  )

  gene_info <- biomaRt::getBM(
    attributes = c(
      "ensembl_gene_id",
      "external_gene_name",
      "gene_biotype",
      "description",
      "transcript_biotype",
      "start_position",
      "end_position"
    ),
    filters = "ensembl_gene_id",
    values = unique(gene_ids),
    mart = ensembl
  )

  gene_info$gene_length <- gene_info$end_position - gene_info$start_position + 1

  return(gene_info)
}

# =============================================================================
# 19. Clean gene ids
# =============================================================================
#' @keywords internal

.clean_gene_ids <- function(expr_mat) {

  gene_ids <- rownames(expr_mat)
  gene_ids_clean <- sub("\\..*$", "", gene_ids)

  tibble::tibble(
    original_id = gene_ids,
    ensembl_gene_id = gene_ids_clean
  )

}

# =============================================================================
# 20. Bootstrap AUC
# =============================================================================
#' @keywords internal

.bootstrap_auc <- function(df, group_col, safe_gene_labels, method, score_method, n_boot = 1000) {

  aucs <- numeric(n_boot)

  for (i in seq_len(n_boot)) {

    idx <- unlist(
      lapply(split(seq_len(nrow(df)), df[[group_col]]),
             function(i) sample(i, length(i), replace = TRUE))
    )

    df_boot <- df[idx, , drop = FALSE]

    # ===============================================
    # 20.1. prediction
    # ===============================================

    if (method == "single_gene") {

      pred <- df_boot[[safe_gene_labels[1]]]

    } else if (method == "logistic") {

      formula_str <- paste(group_col, "~", paste(safe_gene_labels, collapse = " + "))

      model <- stats::glm(
        as.formula(formula_str),
        data = df_boot,
        family = stats::binomial
      )

      pred <- stats::predict(model, type = "response")

    } else if (method == "signature") {

      if (score_method == "mean") {

        pred <- rowMeans(df_boot[, safe_gene_labels, drop = FALSE])

      } else if (score_method == "pca") {

        pca <- prcomp(df_boot[, safe_gene_labels, drop = FALSE], scale. = TRUE)
        pc1 <- pca$x[, 1]

        # Force direction
        if (cor(pc1, df_boot[[group_col]]) < 0) {
          pc1 <- -pc1
        }

        pred <- pc1

      }
    }

    # ===============================================
    # 20.2. ROC
    # ===============================================

    if (length(unique(pred)) < 2) {
      aucs[i] <- NA
      next
    }

    roc_obj <- pROC::roc(
      df_boot[[group_col]],
      pred,
      quiet = TRUE
    )

    aucs[i] <- as.numeric(pROC::auc(roc_obj))
  }

  aucs <- aucs[!is.na(aucs)]

  ci_boot <- quantile(aucs, probs = c(0.025, 0.975), na.rm = TRUE)

  list(
    auc_boot = aucs,
    ci_boot = ci_boot
  )
}


# =============================================================================
# 21. Prepare expression data for WGCNA
# =============================================================================
#' @keywords internal

.prepare_wgcna_data <- function(
    expr_mat,
    metadata,
    genes_n = 10000,
    verbose = TRUE
) {

# ===============================================
# 21.1. Validate arguments
# ===============================================

if (!is.null(genes_n)) {

  if (length(genes_n) != 1L ||
      !is.numeric(genes_n) ||
      !is.finite(genes_n) ||
      genes_n < 1 ||
      genes_n != floor(genes_n) ||
      genes_n > .Machine$integer.max) {

    stop("`genes_n` must be a positive integer or NULL.")
  }

  genes_n <- as.integer(genes_n)
}

if (!is.logical(verbose) ||
    length(verbose) != 1L ||
    is.na(verbose)) {

  stop("`verbose` must be TRUE or FALSE.")
}

# ===============================================
# 21.2. Validate expression matrix
# ===============================================

if (!is.matrix(expr_mat) || !is.numeric(expr_mat)) {
  stop("Expression data must be a numeric matrix.")
}

if (nrow(expr_mat) < 2L || ncol(expr_mat) < 3L) {
  stop(
    "Expression matrix must contain at least 2 genes ",
    "and 3 samples."
  )
}

if (is.null(rownames(expr_mat)) ||
    is.null(colnames(expr_mat)) ||
    anyNA(rownames(expr_mat)) ||
    anyNA(colnames(expr_mat)) ||
    any(!nzchar(trimws(rownames(expr_mat)))) ||
    any(!nzchar(trimws(colnames(expr_mat)))) ||
    anyDuplicated(rownames(expr_mat)) ||
    anyDuplicated(colnames(expr_mat))) {

  stop(
    "Expression matrix must have unique, non-empty gene ",
    "and sample identifiers."
  )
}

# Work on local objects.
expr <- expr_mat

# ===============================================
# 21.3. Validate and align metadata
# ===============================================

if (is.null(metadata) || !is.data.frame(metadata)) {
  stop("Metadata must be a data.frame.")
}

if ("Sample" %in% colnames(metadata)) {

  sample_ids <- as.character(metadata$Sample)

} else {

  sample_ids <- rownames(metadata)

  if (is.null(sample_ids)) {
    stop(
      "Metadata must contain a 'Sample' column ",
      "or sample IDs as row names."
    )
  }
}

if (anyNA(sample_ids) ||
    any(!nzchar(trimws(sample_ids))) ||
    anyDuplicated(sample_ids)) {

  stop(
    "Metadata sample identifiers must be non-missing, ",
    "non-empty and unique."
  )
}

rownames(metadata) <- sample_ids

missing_samples <- setdiff(colnames(expr), sample_ids)

if (length(missing_samples) > 0L) {
  stop(
    "The following expression samples are missing from metadata: ",
    paste(missing_samples, collapse = ", ")
  )
}

# Preserve expression-matrix sample order.
metadata <- metadata[colnames(expr), , drop = FALSE]

# ===============================================
# 21.4. Record initial dimensions
# ===============================================

initial_genes <- rownames(expr)
initial_samples <- colnames(expr)

n_genes_initial <- nrow(expr)
n_samples_initial <- ncol(expr)

if (verbose) {
  message(
    "[WGCNA] Initial data: ",
    n_genes_initial, " genes and ",
    n_samples_initial, " samples."
  )
}

# ===============================================
# 21.5. Exclude genes with invalid or zero variance
# ===============================================

gene_var <- apply(expr, 1L, stats::var)

eligible <- is.finite(gene_var) & gene_var > 0

genes_invalid_variance <- rownames(expr)[!eligible]

expr <- expr[eligible, , drop = FALSE]
gene_var <- gene_var[eligible]

if (nrow(expr) < 2L) {
  stop("Fewer than 2 genes with non-zero variance remain.")
}

# ===============================================
# 21.6. Select most variable genes
# ===============================================

ord <- order(gene_var, decreasing = TRUE)

if (!is.null(genes_n)) {
  ord <- head(ord, min(genes_n, length(ord)))
}

expr <- expr[ord, , drop = FALSE]

genes_selected <- rownames(expr)

genes_excluded_selection <- setdiff(
  setdiff(initial_genes, genes_invalid_variance),
  genes_selected
)

if (verbose) {
  message(
    "[WGCNA] Selected ",
    nrow(expr), " variable genes."
  )
}

# ===============================================
# 21.7. Prepare data for WGCNA
# ===============================================

# WGCNA expects samples in rows and genes in columns.
datExpr <- t(expr)

# ===============================================
# 21.8. Quality control
# ===============================================

gsg <- WGCNA::goodSamplesGenes(
  datExpr,
  verbose = 0
)

# Record exclusions before filtering.
samples_removed_qc <- rownames(datExpr)[!gsg$goodSamples]
genes_removed_qc <- colnames(datExpr)[!gsg$goodGenes]

if (length(samples_removed_qc) > 0L ||
    length(genes_removed_qc) > 0L) {

  if (verbose) {
    message(
      "[WGCNA] Removing ",
      length(samples_removed_qc), " problematic samples and ",
      length(genes_removed_qc), " problematic genes."
    )
  }
}

datExpr <- datExpr[
  gsg$goodSamples,
  gsg$goodGenes,
  drop = FALSE
]

# Keep metadata synchronized with retained samples.
metadata <- metadata[
  rownames(datExpr),
  ,
  drop = FALSE
]

if (nrow(datExpr) < 3L || ncol(datExpr) < 2L) {
  stop(
    "Insufficient samples or genes remain after quality control."
  )
}

# ===============================================
# 21.9. Return prepared data and provenance
# ===============================================

list(
  datExpr = datExpr,
  metadata = metadata,

  n_genes_initial = n_genes_initial,
  n_samples_initial = n_samples_initial,

  n_genes = ncol(datExpr),
  n_samples = nrow(datExpr),

  genes_used = colnames(datExpr),
  samples_used = rownames(datExpr),

  genes_invalid_variance = genes_invalid_variance,
  genes_excluded_selection = genes_excluded_selection,
  genes_removed_qc = genes_removed_qc,
  samples_removed_qc = samples_removed_qc
)
}

# =============================================================================
# 22. Prepare expression data for WGCNA
# =============================================================================
#' @keywords internal

.calculate_wgcna_metrics <- function(datExpr, MEs, wgcna_obj) {

  net <- wgcna_obj$network
  genes_modules <- wgcna_obj$genes_modules

  # ===============================================
  # 22.1. Intramodular connectivity
  # ===============================================

  connectivity <- vector("list", length(net$blockGenes))

  for (b in seq_along(net$blockGenes)) {

    block_idx <- net$blockGenes[[b]]

    block_expr <- datExpr[, block_idx, drop = FALSE]
    block_colors <- net$colors[block_idx]

    adj <- WGCNA::adjacency(
      block_expr,
      power = wgcna_obj$power,
      type = wgcna_obj$network_type
    )

    connectivity[[b]] <- WGCNA::intramodularConnectivity(
      adj,
      colors = block_colors
    )

    rownames(connectivity[[b]]) <- colnames(block_expr)
  }

  connectivity <- do.call(rbind, connectivity)

  # ===============================================
  # 22.2. Module membership (kME)
  # ===============================================

  kME_matrix <- stats::cor(
    datExpr,
    MEs,
    use = "pairwise.complete.obs",
    method = wgcna_obj$cor_method
  )

  gene_module <- as.character(genes_modules$Module)

  gene_idx <- match(
    genes_modules$Gene,
    rownames(kME_matrix)
  )

  kME <- vapply(seq_along(gene_module), function(i) {

    if (is.na(gene_module[i]) ||
        gene_module[i] %in% c("0", "grey")) {
      return(NA_real_)
    }

    me_name <- paste0("ME", gene_module[i])

    if (!me_name %in% colnames(kME_matrix)) {
      return(NA_real_)
    }

    kME_matrix[gene_idx[i], me_name]

  }, numeric(1))

  # ===============================================
  # 22.3. Organize gene metrics
  # ===============================================

  connectivity_idx <- match(
    genes_modules$Gene,
    rownames(connectivity)
  )

  gene_metrics <- data.frame(
    Gene = genes_modules$Gene,
    Module = genes_modules$Module,
    ModuleColor = genes_modules$ModuleColor,
    kME = kME,
    kWithin = connectivity[connectivity_idx, "kWithin"],
    kTotal = connectivity[connectivity_idx, "kTotal"],
    kOut = connectivity[connectivity_idx, "kOut"],
    stringsAsFactors = FALSE
  )

  # ===============================================
  # 22.4. Safe summary functions
  # ===============================================

  safe_mean <- function(x) {
    if (all(is.na(x))) return(NA_real_)
    mean(x, na.rm = TRUE)
  }

  safe_median <- function(x) {
    if (all(is.na(x))) return(NA_real_)
    stats::median(x, na.rm = TRUE)
  }

  # ===============================================
  # 22.5. Module summary
  # ===============================================

  module_data <- gene_metrics[
    !is.na(gene_metrics$ModuleColor) &
      gene_metrics$ModuleColor != "grey",
    ,
    drop = FALSE
  ]

  module_groups <- split(
    module_data,
    module_data$ModuleColor
  )

  module_summary <- do.call(
    rbind,
    lapply(module_groups, function(x) {

      data.frame(
        Module = x$Module[1],
        ModuleColor = x$ModuleColor[1],
        NGenes = nrow(x),
        Mean_kME = safe_mean(x$kME),
        Median_kME = safe_median(x$kME),
        Mean_kWithin = safe_mean(x$kWithin),
        Median_kWithin = safe_median(x$kWithin),
        stringsAsFactors = FALSE
      )
    })
  )

  rownames(module_summary) <- NULL

  module_summary <- module_summary[
    order(module_summary$NGenes, decreasing = TRUE),
  ]

  list(
    gene_metrics = gene_metrics,
    module_summary = module_summary
  )
}
