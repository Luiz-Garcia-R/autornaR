# --- Auxiliary general functions ---

# ============================
# Global cache environment
# ============================
#' @keywords internal

.autornar_cache <- new.env(parent = emptyenv())
.autornar_cache$gene_maps <- list()

# ============================
# Attach to project
# ============================
#' @keywords internal

.attach_to_project <- function(proj,
                               obj,
                               slot,
                               subtype,
                               prefix,
                               log = list(),
                               id = NULL) {

  slot <- as.character(slot)

  if (!slot %in% c("input", "data", "analyses", "qc")) {
    stop("Invalid slot")
  }

  # Ensure base structure
  if (is.null(proj[[slot]])) {
    proj[[slot]] <- list()
  }

  if (is.null(proj[[slot]][[subtype]])) {
    proj[[slot]][[subtype]] <- list()
  }

  if (!is.list(proj[[slot]][[subtype]])) {
    proj[[slot]][[subtype]] <- list(proj[[slot]][[subtype]])
  }

  if (is.null(id)) {
    id <- paste0(prefix, "_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  }

  # Ensure no structural collision
  if (!is.list(proj[[slot]][[subtype]])) {
    proj[[slot]][[subtype]] <- list()
  }

  proj[[slot]][[subtype]][[id]] <- obj
  proj[[slot]][[subtype]]$last <- id

  # Safe logs
  if (is.null(proj$logs)) proj$logs <- list()
  if (is.null(proj$logs[[subtype]])) proj$logs[[subtype]] <- list()

  proj$logs[[subtype]][[id]] <- c(list(timestamp = Sys.time()), log)
  proj$logs[[subtype]]$last <- id

  proj
}

# ============================
# Get last or selected
# ============================
#' @keywords internal

.get_last_or_selected <- function(x, id = NULL, what = "element") {

  if (length(x) == 0) {
    stop(sprintf("No %s found.", what))
  }

  ids <- setdiff(names(x), "last")

  if (is.null(id)) {

    if (!"last" %in% names(x)) {
      stop(sprintf("No %s found.", what))
    }

    return("last")
  }

  if (is.numeric(id)) {

    if (length(id) != 1) {
      stop(sprintf("%s index must be a single number.", what))
    }

    if (id < 1 || id > length(ids)) {
      stop(sprintf("%s index out of bounds.", what))
    }

    return(ids[id])
  }

  if (is.character(id)) {

    if (length(id) != 1) {
      stop(sprintf("%s id must be a single string.", what))
    }

    if (!id %in% names(x)) {
      stop(sprintf("%s '%s' not found.", what, id))
    }

    return(id)
  }

  stop("Invalid 'id' argument.")
}

# ============================
# Set seed
# ============================
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

# ============================
# Reset seed
# ============================
#' @keywords internal

.reset_seed <- function(old_seed) {
  if (!is.null(old_seed)) {
    .Random.seed <<- old_seed
  }
}

# ============================
# Get Biomart dataset
# ============================
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

# ============================
# Order contrast levels
# ============================
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

# ============================
# Check dependencies
# ============================
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

# ============================
# Smart pathway naming
# ============================
#' @keywords internal

.smart_pathway_name <- function(x, max_words = 6) {

  vapply(
    x,
    function(xx) {

      if (is.na(xx) || !nzchar(xx)) {
        return("Unknown_Pathway")
      }

      # 1. Remove prefixos
      xx <- gsub("^(REACTOME_|KEGG_|HALLMARK_|GOBP_|GO_)", "", xx)

      # 2. Normaliza
      xx <- gsub("_", " ", xx)
      xx <- tolower(xx)

      # 3. Substituições inteligentes (compressão semântica)
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

# ============================
# Auxiliary print functions
# ============================
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

# ============================
# Go enrichment
# ============================
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

# ============================
# Heatmap gene selection
# ============================
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

# ============================
# Convert gene ids
# ============================
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


# ============================
# Map gene symbols
# ============================
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

# ============================
# Get gene annotation
# ============================
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

# ============================
# Detect organism
# ============================
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

# ============================
# Normalize sample names
# ============================
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

# ============================
# Get biomart gene info
# ============================
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

# ============================
# Clean gene ids
# ============================
#' @keywords internal

.clean_gene_ids <- function(expr_mat) {

  gene_ids <- rownames(expr_mat)
  gene_ids_clean <- sub("\\..*$", "", gene_ids)

  tibble::tibble(
    original_id = gene_ids,
    ensembl_gene_id = gene_ids_clean
  )

}

# ============================
# Bootstrap AUC
# ============================
#' @keywords internal

.bootstrap_auc <- function(df, group_col, safe_gene_labels, method, score_method, n_boot = 1000) {

  aucs <- numeric(n_boot)

  for (i in seq_len(n_boot)) {

    idx <- unlist(
      lapply(split(seq_len(nrow(df)), df[[group_col]]),
             function(i) sample(i, length(i), replace = TRUE))
    )

    df_boot <- df[idx, , drop = FALSE]

    # --- prediction ---
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

    # --- ROC ---
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

