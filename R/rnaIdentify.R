#' Identify and classify RNA types from normalized RNA-seq data
#'
#' This function uses Ensembl BioMart to annotate gene biotypes and classify
#' expressed transcripts (including lncRNAs) in normalized RNA-seq data. It can
#' optionally generate a summary plot of RNA type proportions per group.
#'
#' @param normalized_data A `normalized_data` object returned by \code{rna_normalize()},
#'   containing an expression matrix (`expr_matrix`) and associated metadata.
#' @param group_col Character. Column name in metadata specifying the sample groups
#'   (default: `"Group"`).
#' @param species Character. Either `"human"`, `"mouse"` or `"zebrafish"` to define the Ensembl dataset
#'   (default: `"mouse"`).
#' @param top_n_types Integer. Maximum number of RNA biotypes displayed individually
#'   in the summary plot; remaining types are grouped as `"Other"` (default: `5`).
#' @param plot Logical. Whether to display a stacked barplot summarizing RNA types
#'   per group (default: `TRUE`).
#' @param lnc_only Logical. If `TRUE`, filters and classifies only lncRNAs
#'   (default: `FALSE`).
#'
#' @return A list of class `"rnaIdentify"` containing:
#' \item{rna_lists}{A named list of gene identifiers by RNA biotype.}
#' \item{plot_bar}{A `ggplot` object of the RNA biotype distribution (if `plot = TRUE`).}
#' \item{info}{A data frame with Ensembl annotations, gene types, and lncRNA subclassification.}
#' \item{summary_by_group}{A data frame summarizing RNA biotype proportions per group.}
#'
#' @details
#' The function connects to Ensembl via \pkg{biomaRt} to retrieve gene biotypes and
#' additional metadata. It supports human, mouse and zebra fish datasets.
#' When `lnc_only = TRUE`, only long non-coding RNAs (lncRNAs) are kept and subclassified
#' into categories (e.g., *Antisense*, *Intergenic*, *Processed transcript*).
#'
#' @examples
#' \dontrun{
#' result <- rna.identify(normalized_data, group_col = "Condition")
#' result_lnc <- rna.identify(normalized_data, lnc_only = TRUE)
#' }
#'
#' @export

rna.identify <- function(normalized_data,
                         group_col = "Group",
                         organism = "mouse",
                         top_n_types = 5,
                         plot = TRUE,
                         lnc_only = FALSE) {

  # --- Necessary packages ---
  req_pkgs <- c("biomaRt", "dplyr", "tibble", "tidyr", "ggplot2", "stringr", "rlang", "glue")
  missing_pkgs <- req_pkgs[!vapply(req_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs)) stop("Please install packages: ", paste(missing_pkgs, collapse = ", "))

  # --- Data validation ---
  if (!is.list(normalized_data) || is.null(normalized_data$expr_matrix))
    stop("Input must be a list returned by rna.normalize() (containing expr_matrix).")
  if (is.null(normalized_data$metadata))
    stop("The 'normalized_data' object must contain a 'metadata' slot.")

  expr_mat <- as.matrix(normalized_data$expr_matrix)
  metadata <- normalized_data$metadata

  if (!group_col %in% colnames(metadata))
    stop(paste0("`metadata` must contain column '", group_col, "'"))

  # --- Clean gene IDs ---
  gene_ids <- rownames(expr_mat)
  gene_ids_clean <- sub("\\..*$", "", gene_ids)          # remove version suffix
  id_map <- tibble::tibble(original_id = gene_ids, ensembl_gene_id = gene_ids_clean)

  # --- UX auto-detect organism if not specified / or blank ---
  if(is.null(organism) || organism == ""){
    if(all(grepl("^ENSDARG", gene_ids_clean)))      organism <- "zebrafish"
    else if(all(grepl("^ENSMUSG", gene_ids_clean))) organism <- "mouse"
    else if(all(grepl("^ENSG", gene_ids_clean)))    organism <- "human"
    else                                            organism <- "mouse"   # fallback
    message(glue::glue("[rna.identify] Organism auto-detected → {organism}"))
  }

  # --- Choose dataset (fallback to mouse) ---
  dataset <- switch(
    tolower(organism),
    "human"     = "hsapiens_gene_ensembl",
    "h"         = "hsapiens_gene_ensembl",
    "mouse"     = "mmusculus_gene_ensembl",
    "m"         = "mmusculus_gene_ensembl",
    "zebrafish" = "drerio_gene_ensembl",
    "danio"     = "drerio_gene_ensembl",
    "dre"       = "drerio_gene_ensembl",
    "mmusculus_gene_ensembl" # fallback
  )

  # --- Connect to BioMart with mirror fallbacks ---
  ensembl <- tryCatch({
    biomaRt::useEnsembl(biomart = "genes", dataset = dataset, mirror = "www")
  }, error = function(e1) {
    tryCatch(biomaRt::useEnsembl(biomart = "genes", dataset = dataset, mirror = "uswest"),
             error = function(e2) {
               tryCatch(biomaRt::useMart("ENSEMBL_MART_ENSEMBL", dataset = dataset),
                        error = function(e3) NULL)
             })
  })
  if (is.null(ensembl)) stop("[rna.identify] Could not connect to Ensembl BioMart (check internet / dataset).")

  # --- Check available attributes and pick symbol attr ---
  attrs <- biomaRt::listAttributes(ensembl)$name
  if (tolower(organism) %in% c("zebrafish", "danio", "dre")) {
    symbol_attr <- if ("zfin_id_symbol" %in% attrs) "zfin_id_symbol" else "external_gene_name"
  } else {
    symbol_attr <- if ("hgnc_symbol" %in% attrs) "hgnc_symbol" else "external_gene_name"
  }

  # Build attributes to request (only those present)
  attrs_to_use <- c()
  if ("ensembl_gene_id" %in% attrs) attrs_to_use <- c(attrs_to_use, "ensembl_gene_id")
  if ("gene_biotype" %in% attrs) attrs_to_use <- c(attrs_to_use, "gene_biotype")
  if (symbol_attr %in% attrs) attrs_to_use <- c(attrs_to_use, symbol_attr)
  if ("description" %in% attrs) attrs_to_use <- c(attrs_to_use, "description")
  if ("transcript_biotype" %in% attrs) attrs_to_use <- c(attrs_to_use, "transcript_biotype")

  # --- Decide primary filter: ENSEMBL-like IDs -> use ensembl_gene_id; else use external_gene_name ---
  looks_like_ensembl <- all(grepl("^ENS", gene_ids_clean, ignore.case = TRUE))
  primary_filter <- if (looks_like_ensembl && "ensembl_gene_id" %in% attrs) "ensembl_gene_id" else "external_gene_name"

  # --- Query BioMart (try primary filter; if no hits, try the other filter) ---
  gene_info <- tryCatch({
    biomaRt::getBM(attributes = attrs_to_use, filters = primary_filter,
                   values = unique(gene_ids_clean), mart = ensembl)
  }, error = function(e) NULL)

  # If no results and primary was ensembl -> try external_gene_name
  if ((is.null(gene_info) || nrow(gene_info) == 0) && primary_filter == "ensembl_gene_id") {
    gene_info <- tryCatch({
      biomaRt::getBM(attributes = attrs_to_use, filters = "external_gene_name",
                     values = unique(gene_ids_clean), mart = ensembl)
    }, error = function(e) NULL)
  }

  # If no results and primary was external -> try ensembl
  if ((is.null(gene_info) || nrow(gene_info) == 0) && primary_filter == "external_gene_name") {
    gene_info <- tryCatch({
      biomaRt::getBM(attributes = attrs_to_use, filters = "ensembl_gene_id",
                     values = unique(gene_ids_clean), mart = ensembl)
    }, error = function(e) NULL)
  }

  if (is.null(gene_info) || nrow(gene_info) == 0) {
    stop("[rna.identify] BioMart did not return any matches. Check that your rownames are Ensembl IDs (no suffix) or gene symbols, and that 'organism' is correct.")
  }

  # --- Normalize column types and symbol name ---
  if ("ensembl_gene_id" %in% colnames(gene_info))
    gene_info$ensembl_gene_id <- as.character(gene_info$ensembl_gene_id)
  if (symbol_attr %in% colnames(gene_info) && symbol_attr != "hgnc_symbol") {
    # normalize to gene_symbol for downstream clarity
    gene_info <- gene_info |> dplyr::rename(gene_symbol = !!rlang::sym(symbol_attr))
  } else if ("hgnc_symbol" %in% colnames(gene_info)) {
    gene_info <- gene_info |> dplyr::rename(gene_symbol = hgnc_symbol)
  } else {
    gene_info$gene_symbol <- NA_character_
  }
  # make sure gene_biotype/description/transcript_biotype exist
  if (!"gene_biotype" %in% colnames(gene_info)) gene_info$gene_biotype <- NA_character_
  if (!"description" %in% colnames(gene_info)) gene_info$description <- NA_character_
  if (!"transcript_biotype" %in% colnames(gene_info)) gene_info$transcript_biotype <- NA_character_

  # --- Merge id_map with gene_info safely (both char) ---
  id_map$ensembl_gene_id <- as.character(id_map$ensembl_gene_id)
  gene_info$ensembl_gene_id <- as.character(gene_info$ensembl_gene_id)

  info <- dplyr::left_join(id_map, gene_info, by = "ensembl_gene_id") |>
    dplyr::distinct(ensembl_gene_id, .keep_all = TRUE) |>
    dplyr::mutate(
      description = dplyr::na_if(description, ""),
      gene_symbol = dplyr::na_if(gene_symbol, ""),
      transcript_biotype = if (any(!is.na(gene_info$transcript_biotype))) transcript_biotype else gene_biotype
    )

  # --- Structural classification (lnc_type) based on description/text heuristics ---
  info <- info |>
    dplyr::mutate(
      lnc_type = dplyr::case_when(
        stringr::str_detect(description, stringr::regex("antisense", ignore_case = TRUE)) ~ "Antisense",
        stringr::str_detect(description, stringr::regex("intergenic|LINC", ignore_case = TRUE)) ~ "Intergenic",
        stringr::str_detect(description, stringr::regex("sense overlapping", ignore_case = TRUE)) ~ "Sense overlapping",
        stringr::str_detect(description, stringr::regex("intronic", ignore_case = TRUE)) ~ "Sense intronic",
        stringr::str_detect(description, stringr::regex("bidirectional", ignore_case = TRUE)) ~ "Bidirectional",
        stringr::str_detect(description, stringr::regex("processed transcript", ignore_case = TRUE)) ~ "Processed transcript",
        stringr::str_detect(description, stringr::regex("novel transcript", ignore_case = TRUE)) ~ "Novel lncRNA",
        TRUE ~ "Unspecified lncRNA"
      ),
      # patch UX
      lnc_type = ifelse(
        gene_biotype == "lincRNA" & lnc_type == "Unspecified lncRNA",
        "Intergenic",
        lnc_type
      )
    )

  # --- Filter only classic lncRNA if requested ---
  if (lnc_only) {
    lncRNA_biotypes <- c("lncRNA", "antisense", "lincRNA", "processed_transcript")
    info <- dplyr::filter(info, gene_biotype %in% lncRNA_biotypes)
  }

  # --- RNA lists by type ---
  rna_types <- unique(info$gene_biotype)
  rna_lists <- lapply(rna_types, function(t) info$original_id[info$gene_biotype == t])
  names(rna_lists) <- rna_types

  # --- Prepare plotting data ---
  expr_detected <- expr_mat > 0
  expr_long <- tibble::as_tibble(expr_detected, rownames = "original_id") |>
    tidyr::pivot_longer(cols = -original_id, names_to = "Sample", values_to = "expressed") |>
    dplyr::filter(expressed) |>
    dplyr::left_join(info[, c("original_id", "gene_biotype", "lnc_type")], by = "original_id") |>
    dplyr::left_join(metadata[, c("Sample", group_col)], by = "Sample")

  # --- Summarize and plot ---
  type_counts <- table(expr_long$gene_biotype)
  top_types <- names(sort(type_counts, decreasing = TRUE))[1:min(top_n_types, length(type_counts))]
  expr_long$gene_biotype_plot <- ifelse(expr_long$gene_biotype %in% top_types, expr_long$gene_biotype, "Others")

  expr_summary <- expr_long |>
    dplyr::count(!!rlang::sym(group_col), gene_biotype_plot) |>
    dplyr::group_by(!!rlang::sym(group_col)) |>
    dplyr::mutate(prop = n / sum(n)) |>
    dplyr::ungroup()

  p_bar <- NULL
  if (plot & !lnc_only) {
    p_bar <- ggplot2::ggplot(expr_summary, ggplot2::aes(
      x = !!rlang::sym(group_col), y = prop, fill = gene_biotype_plot
    )) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = "RNA type distribution per group",
        x = group_col,
        y = "Proportion of expressed genes",
        fill = "RNA type"
      )
    print(p_bar)
  }

  # --- lnc summary when requested ---
  if (lnc_only) {
    lnc_summary <- table(info$lnc_type)
    message("[rna_identify] lncRNA classification summary:")
    for (n in names(lnc_summary)) {
      message(" - ", n, ": ", lnc_summary[[n]])
    }
  }

  # --- Return results ---
  result <- list(
    rna_lists = rna_lists,
    plot_bar = p_bar,
    info = info,
    summary_by_group = expr_summary,
    genes_by_group = NULL
  )

  if (lnc_only) assign("lncRNAs_identified", result, envir = .GlobalEnv)

  class(result) <- "rnaIdentify"
  invisible(result)
}
