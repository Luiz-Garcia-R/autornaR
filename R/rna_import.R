#' Import RNA-seq count data from multiple formats
#'
#' @description
#' Imports and standardizes RNA-seq count matrices from common formats
#' (HISAT2, STAR, featureCounts, tximport, or pre-cleaned matrices).
#'
#' The function automatically detects gene identifiers (either as a column
#' or row names), sample columns, and gene ID types. It performs basic
#' validation and cleaning steps to ensure compatibility with downstream
#' analysis.
#'
#' When ENSEMBL identifiers are detected, version numbers are removed and
#' duplicated IDs are collapsed by summation. Sample names can be optionally
#' cleaned or manually renamed.
#'
#' If metadata is provided, the function automatically detects the sample
#' identifier column and validates its correspondence with the count matrix.
#' Summary information (e.g., group sizes) is reported when available.
#'
#' @param project \code{rna_project} object created by \code{rna.project()}.
#' @param raw_data Data frame containing raw counts and gene identifiers.
#' @param metadata Optional data frame with sample information. The sample
#' column is automatically detected based on column names or matching values
#' with the count matrix. If detection fails, the function will return an error
#' (in strict mode).
#' @param format Character. Input format: \code{"clean"}, \code{"hisat2"},
#'   \code{"star"}, \code{"featureCounts"}, or \code{"tximport"}.
#' @param gene_col Optional character specifying the column containing gene
#' identifiers. If \code{NULL} (default), the function will attempt to detect
#' the gene column automatically. If no suitable column is found, row names
#' will be used when possible.
#' @param group_col Optional character specifying the column containing sample
#' groups. If \code{NULL} (default).
#' @param rename_samples Optional character vector to rename sample columns.
#' Must have the same length and order as the detected sample columns.
#' @param gene_id_type Character. Gene identifier type:
#'   \code{"auto"}, \code{"ENSEMBL"}, \code{"SYMBOL"}, \code{"ENTREZ"}.
#'   Default is \code{"auto"}.
#' @param organism Character. Organism name:
#'   \code{"auto"}, \code{"mouse"}, \code{"human"}, or \code{"zebrafish"}.
#' @param clean_names Logical. If \code{TRUE}, automatically cleans sample names
#'   by removing file extensions and special characters.
#' @param strict Logical. If \code{TRUE}, stops execution on validation errors.
#' @param save Logical. Whether to store results in the active \code{rna_project}.
#' @param envir Environment where to assign the object.
#'
#' @return
#' An object of class \code{"imp_data"}, containing:
#' \itemize{
#'   \item \code{data}: Processed count matrix.
#'   \item \code{metadata}: Metadata table (if provided).
#'   \item \code{n_genes}: Number of genes.
#'   \item \code{n_samples}: Number of samples.
#'   \item \code{gene_id_type}: Detected or specified gene ID type.
#'   \item \code{organism}: Organism used.
#'   \item \code{warnings}: Validation warnings generated during import.
#'   \item \code{groups}: (optional) Named list with group sizes if a
#'     \code{Group} column is present in metadata.
#'   \item \code{gene_annotation}: Matrix with annotation of all possible genes
#'   of the provided data set.
#'   \item \code{warnings}: Validation warnings generated during import,
#'   including automatic detection steps.
#'}
#'
#' @details
#' This function serves as the entry point of the RNA-seq analysis workflow.
#' It is designed to handle heterogeneous input formats with minimal user
#' intervention by automatically detecting key structural components such as
#' gene identifiers and sample labels.
#'
#' The function does not perform normalization or transformation of counts.
#'
#' @examples
#' \dontrun{
#' # Basic import with metadata
#'   rna.import(project = my_project,
#'              raw_data = counts_df,
#'              metadata = meta_df)
#'
#' # Import with metadata and specifying an organism
#'   rna.import(my_project,
#'              raw_data = counts_df,
#'              metadata = meta_df,
#'              organism = "mouse")
#'
#' # Import with a specific format
#'   rna.import(my_project,
#'              raw_data = counts_df,
#'              metadata = meta_df,
#'              format = "star")
#'  }
#'
#' @importFrom stats na.omit
#'
#' @export

rna.import <- function(
    project,
    raw_data,
    metadata = NULL,
    format = c("clean", "hisat2", "star", "featureCounts", "tximport"),
    gene_col = NULL,
    group_col = NULL,
    gene_id_type = c("auto", "ENSEMBL", "SYMBOL", "ENTREZ"),
    organism = c("auto", "mouse", "human", "zebrafish"),
    rename_samples = NULL,
    clean_names = TRUE,
    strict = TRUE,
    save = TRUE,
    envir = parent.frame()
) {

  gene_id_type <- match.arg(gene_id_type)
  organism <- match.arg(organism)
  format <- match.arg(format)

  # ---------------------------
  # 1) Get active project
  # ---------------------------
  if (is.null(project)) {
    stop("Project must be provided. Use rna.project() first.")
  }
  proj <- project

  # ---------------------------
  # 2) Setup validation
  # ---------------------------
  errors <- character()
  warnings_list <- character()

  add_error <- function(msg) errors <<- c(errors, msg)
  add_warning <- function(msg) warnings_list <<- c(warnings_list, msg)

  if (length(errors) > 0 && strict) {
    stop(paste(errors, collapse = "\n"))
  }

  # ---------------------------
  # 3) Validate input
  # ---------------------------
  if (!is.data.frame(raw_data)) {
    add_error("'raw_data' must be a data.frame.")
  }

  # ---------------------------
  # 3.5) Handle rownames as gene IDs
  # ---------------------------
  if (is.null(colnames(raw_data)) || ncol(raw_data) == 0) {
    add_error("raw_data has no columns.")
  }

  # Check column
  possible_gene_cols <- c(
    "Gene_ID",
    "Geneid",
    "geneid",
    "GeneID",
    "GENEID",
    "Ensembl_ID",
    "ENSEMBL_ID",
    "ensembl_id",
    "Gene",
    "gene",
    "ID",
    "id",
    "Symbol",
    "symbol"
  )

  found <- intersect(possible_gene_cols, colnames(raw_data))

  if (length(found) == 0) {

    rn <- rownames(raw_data)

    # Check rownames
    if (!is.null(rn) && length(rn) > 0) {

      rn_sample <- rn[1:min(100, length(rn))]

      looks_like_gene <- mean(
        grepl("^ENS|^[A-Za-z0-9\\-\\.]+$", rn_sample)
      ) > 0.8

      if (looks_like_gene) {

        if (is.null(gene_col)) {

          stop(
            "No gene identifier column detected. Please specify 'gene_col'."
          )

        }

        raw_data[[gene_col]] <- rn

        raw_data <- raw_data[, c(
          gene_col,
          setdiff(colnames(raw_data), gene_col)
        )]

        add_warning(
          "Gene IDs detected in rownames and moved to column."
        )
      }
    }
  }

  # ---------------------------
  # 4) Detect gene column
  # ---------------------------
  if (!is.null(gene_col)) {
    possible_gene_cols <- c(gene_col, possible_gene_cols)
  }

  found <- intersect(
    possible_gene_cols,
    colnames(raw_data)
  )

  if (length(found) > 0) {
    gene_col <- found[1]
  }


  # ---------------------------
  # 5) Detect gene ID type
  # ---------------------------
  .detect_gene_id_type <- function(ids) {
    ids <- na.omit(ids)
    ids <- ids[1:min(length(ids), 100)]

    ids_clean <- sub("\\..*$", "", ids)

    if (mean(grepl("^ENS[A-Z]*G[0-9]+$", ids_clean)) > 0.8)
      return("ENSEMBL")

    if (mean(grepl("^[0-9]+$", ids)) > 0.9)
      return("ENTREZ")

    if (mean(!grepl("^ENS", ids) &
             grepl("^[A-Za-z][A-Za-z0-9\\-\\.]*$", ids)) > 0.8)
      return("SYMBOL")

    return(NA_character_)
  }

  # ---------------------------
  # 6) Identify count columns
  # ---------------------------
  if (format %in% c("hisat2", "star", "featureCounts")) {
    count.cols <- grep("\\.bam$|\\.counts$", colnames(raw_data), value = TRUE)

    if (length(count.cols) == 0) {
      meta_cols <- c("Start", "START", "End", "END", "Length", "LENGTH", "Chr", "CHR", "Strand", "STRAND")
      count.cols <- setdiff(
        colnames(raw_data)[sapply(raw_data, is.numeric)],
        c(gene_col, meta_cols)
      )
    }

    if (length(count.cols) == 0)
      add_error("No count columns found.")
  }

  if (format %in% c("clean", "tximport")) {
    count.cols <- setdiff(colnames(raw_data), gene_col)

    if (length(count.cols) == 0)
      add_error("No count columns found.")
  }

  # ---------------------------
  # 7) Build matrix
  # ---------------------------
  if (length(errors) == 0) {

    raw_data_counts <- raw_data[, c(gene_col, count.cols), drop = FALSE]

    if (gene_id_type == "auto") {
      gene_id_type <- .detect_gene_id_type(raw_data_counts[[gene_col]])

      if (is.na(gene_id_type) && strict)
        stop("Could not detect gene_id_type.")
    }

    if (gene_id_type == "ENSEMBL") {

      raw_data_counts[[gene_col]] <- sub("\\..*$", "", raw_data_counts[[gene_col]])

      if (any(duplicated(raw_data_counts[[gene_col]]))) {

        add_warning("Duplicate ENSEMBL IDs collapsed.")

        collapsed <- aggregate(
          raw_data_counts[, -1],
          by = list(Gene = raw_data_counts[[gene_col]]),
          FUN = sum
        )

        colnames(collapsed)[1] <- gene_col
        raw_data_counts <- collapsed
      }
    }
  }

  # ---------------------------
  # 8) Resolve organism
  # ---------------------------
  if (organism == "auto") {

    if (gene_id_type %in% c("SYMBOL", "ENTREZ")) {
      add_error("Cannot infer organism from SYMBOL/ENTREZ IDs. Please specify 'organism'.")

    } else {
      organism <- .detect_organism(raw_data_counts[[gene_col]])

      if (organism == "unknown") {
        add_error("Could not infer organism from ENSEMBL IDs.")

      } else {
        add_warning(paste0(
          "Organism inferred as: ", organism,
          " (based on ENSEMBL prefix)"
        ))
      }
    }

  } else {
    # Manual definition
    if (!organism %in% c("human", "mouse", "zebrafish")) {
      add_error("Invalid organism. Must be 'human', 'mouse', or 'zebrafish'.")
    }
  }

  # ---------------------------
  # 8.1) Map annotation & Extract Length
  # ---------------------------
  gene_ids <- raw_data_counts[[gene_col]]
  gene_ids_clean <- sub("\\..*$", "", gene_ids)

  gene_annotation <- .map_gene_annotation(
    gene_ids_clean,
    organism = organism
  )

  if (is.null(gene_annotation)) {

    gene_annotation <- data.frame(gene_id = gene_ids_clean, stringsAsFactors = FALSE)

  }

  possible_length_cols <- c("Length", "LENGTH", "length")
  length_col_found <- intersect(possible_length_cols, colnames(raw_data))

  if (length(length_col_found) > 0) {

    idx <- match(raw_data_counts[[gene_col]], raw_data[[gene_col]])
    gene_annotation$gene_length <- raw_data[[length_col_found[1]]][idx]

  }

  # ---------------------------
  # 9) Clean sample names
  # ---------------------------
  if (clean_names && format %in% c("hisat2", "star", "featureCounts")) {

    clean_sample_names <- function(x) {
      x <- basename(x)
      x <- gsub("\\.(bam|counts|txt|csv)$", "", x, ignore.case = TRUE)
      x <- gsub("[^A-Za-z0-9_\\-]+", "", x)
      x
    }

    colnames(raw_data_counts)[-1] <- clean_sample_names(
      colnames(raw_data_counts)[-1]
    )
  }

  # ---------------------------
  # 10) Sample names
  # ---------------------------
  if (length(errors) == 0) {

    sample.cols <- setdiff(colnames(raw_data_counts), gene_col)

    if (!is.null(rename_samples)) {

      if (length(rename_samples) == length(sample.cols)) {

        colnames(raw_data_counts)[match(sample.cols, colnames(raw_data_counts))] <- rename_samples
        sample.cols <- rename_samples

      } else {
        add_warning("rename_samples length mismatch.")
      }
    }

  } else {
    sample.cols <- character(0)
  }

  # ---------------------------
  # 11) Metadata
  # ---------------------------
  if (!is.null(metadata)) {

    if (!is.data.frame(metadata))
      add_error("'metadata' must be data.frame.")

    detected_sample_col <- NULL

    # ---------------------------
    # 11.1) Try by column name
    # ---------------------------
    possible_sample_cols <- c(
      "Sample", "sample", "SampleID", "sample_id",
      "Sample_Id", "sampleId", "ID", "id",
      "Run", "run", "filename", "file", "File",
      "Sample_Name", "sample_name"
    )

    found <- intersect(possible_sample_cols, colnames(metadata))

    if (length(found) > 0) {
      detected_sample_col <- found[1]
    }

    # ---------------------------
    # 11.2) Try by matching values
    # ---------------------------
    if (is.null(detected_sample_col)) {

      matches <- sapply(metadata, function(col) {
        if (!is.character(col) && !is.factor(col)) return(0)
        sum(col %in% sample.cols)
      })

      best_col <- names(which.max(matches))

      if (length(best_col) > 0 && matches[best_col] > 0) {
        detected_sample_col <- best_col
        add_warning(paste0(
          "Sample column inferred from values: '", best_col, "'"
        ))
      }
    }

    # ---------------------------
    # 11.3) Final validation
    # ---------------------------
    if (is.null(detected_sample_col)) {
      add_error("Could not detect sample column in metadata.")
    } else {

      # Standard sample name
      metadata$Sample <- as.character(metadata[[detected_sample_col]])
      sample.cols <- as.character(sample.cols)

      # Normalization of sample names
      metadata$Sample <- .normalize_sample_names(metadata$Sample)
      sample.cols <- .normalize_sample_names(sample.cols)

      colnames(raw_data_counts)[
        colnames(raw_data_counts) != gene_col
      ] <- sample.cols

      # Trim
      metadata$Sample <- trimws(metadata$Sample)
      sample.cols <- trimws(sample.cols)

      # Check matching before reordering
      matched <- sum(metadata$Sample %in% sample.cols)

      if (matched == 0) {
        add_error("No matching samples between metadata and count matrix.")
      }

      # Preserve original metadata order
      metadata_original <- metadata

      # Reorder safely
      idx <- match(sample.cols, metadata$Sample)

      if (any(is.na(idx))) {

        missing_samples <- sample.cols[is.na(idx)]

        add_warning(
          paste0(
            "Samples in count matrix not found in metadata: ",
            paste(missing_samples, collapse = ", ")
          )
        )
      }

      metadata <- metadata[idx[!is.na(idx)], , drop = FALSE]

      if (detected_sample_col != "Sample") {
        add_warning(paste0(
          "Using '", detected_sample_col, "' as Sample column."
        ))
      }

      if (matched < length(sample.cols)) {
        missing_in_meta <- setdiff(sample.cols, metadata$Sample)
        missing_in_counts <- setdiff(metadata$Sample, sample.cols)

        add_warning(paste0(
          "Mismatch detected: ",
          length(missing_in_meta), " samples missing in metadata, ",
          length(missing_in_counts), " samples missing in count matrix."
        ))
      }
    }

      # ---------------------------
      # Match validation
      # ---------------------------
      missing <- setdiff(metadata$Sample, sample.cols)

      if (length(missing) > 0)
        add_warning("Metadata samples not in matrix.")
    }

  # ---------------------------
  # 12) Group column validation
  # ---------------------------
  group_info <- NULL

  if (!is.null(metadata)) {

    detected_group_col <- NULL
    original_group_levels <- NULL

    # ---------------------------
    # 12.1) Manual override
    # ---------------------------
    if (!is.null(group_col)) {

      if (group_col %in% colnames(metadata)) {
        detected_group_col <- group_col
      } else {
        add_error(paste0("group_col '", group_col, "' not found in metadata."))
      }

    } else {

      # ---------------------------
      # 12.2) Try common names
      # ---------------------------
      possible_group_cols <- c(
        "Group", "group", "Condition", "condition",
        "Treatment", "treatment", "Class", "class",
        "Phenotype", "phenotype"
      )

      found <- intersect(possible_group_cols, colnames(metadata))

      if (length(found) > 0) {
        detected_group_col <- found[1]
        add_warning(paste0(
          "Group column inferred as: '", detected_group_col, "'"
        ))
      }

      # ---------------------------
      # 12.3) Heuristic fallback
      # ---------------------------
      if (is.null(detected_group_col)) {

        candidates <- sapply(metadata, function(col) {
          if (!is.character(col) && !is.factor(col)) return(Inf)
          n_unique <- length(unique(col))
          if (n_unique > 1 && n_unique <= nrow(metadata) * 0.5) {
            return(n_unique)
          }
          return(Inf)
        })

        best <- names(which.min(candidates))

        if (length(best) > 0 && is.finite(candidates[best])) {
          detected_group_col <- best
          add_warning(paste0(
            "Group column heuristically inferred as: '", best, "'"
          ))
        }
      }
    }

    # ---------------------------
    # 12.4) Build group info
    # ---------------------------
    if (!is.null(detected_group_col)) {

      if (is.null(original_group_levels)) {

        original_group_levels <- unique(
          as.character(metadata_original[[detected_group_col]])
        )

      }

      detected_levels <- original_group_levels

      metadata$Group <- factor(
        metadata[[detected_group_col]],
        levels = detected_levels
      )

      group_info <- list(
        column = detected_group_col,
        levels = detected_levels,
        sizes = as.list(table(metadata$Group))
      )

      add_warning(
        paste0(
          "Group order detected: ",
          paste(detected_levels, collapse = " -> ")
        )
      )

    } else {
      add_warning("No group column detected.")
    }
  }

  # ---------------------------
  # 13) Output object
  # ---------------------------
  obj <- list(
    timestamp = Sys.time(),
    data = raw_data_counts,
    metadata = metadata,
    warnings = warnings_list,
    n_genes = nrow(raw_data_counts),
    n_samples = length(sample.cols),
    detected_format = format,
    gene_id_type = gene_id_type,
    organism = organism,
    groups = group_info,
    gene_annotation = gene_annotation
  )

  class(obj) <- "imp_data"

  # ---------------------------
  # 13) Attach to project
  # ---------------------------
  if (save) {

  proj <- .attach_to_project(
    proj,
    obj,
    slot = "input",
    subtype = "imp_data",
    prefix = "import",
    log = list(
      n_genes = obj$n_genes,
      n_samples = obj$n_samples,
      format = format,
      gene_id_type = gene_id_type,
      organism = organism,
      n_groups = length(group_info),
      group_sizes = group_info,
      warnings = warnings_list
    )
  )
}

  # ---------------------------
  # 14) Return
  # ---------------------------
  counts_only <- raw_data_counts[, setdiff(colnames(raw_data_counts), gene_col), drop = FALSE]

  .print_header("RNA Import")

  .print_block("Overview", function() {
    cat("Format:            ", format, "\n")
    cat("Gene ID type:      ", gene_id_type, "\n")
    cat("Genes:             ", nrow(raw_data_counts), "\n")
    cat("Samples:           ", length(sample.cols), "\n")
    cat("Organism:          ", organism, "\n")
  })

  .print_block("Data integrity", function() {

    # convert safely
    counts_only <- data.frame(lapply(counts_only, as.numeric))

    cat("Missing values:    ", sum(is.na(counts_only)), "\n")
    cat("Min count:         ", suppressWarnings(min(counts_only, na.rm = TRUE)), "\n")
    cat("Max count:         ", suppressWarnings(max(counts_only, na.rm = TRUE)), "\n")

    dup_genes <- sum(duplicated(raw_data_counts[[gene_col]]))
    cat("Duplicated genes:  ", dup_genes, "\n")
    zero_prop <- mean(counts_only == 0, na.rm = TRUE)
    cat("Zero proportion:   ", round(zero_prop, 3), "\n")
  })

  .print_block("Samples (preview)", function() {
    preview <- head(sample.cols, 5)
    cat("First samples:     ", paste(preview, collapse = ", "), "\n")

    if (!is.null(rename_samples)) {
      cat("Custom names:       Yes\n")
    } else {
      cat("Custom names:       No\n")
    }

    cat("Cleaned names:     ", ifelse(clean_names, "Yes", "No"), "\n")
  })

  .print_block("Metadata", function() {

    if (is.null(metadata)) {
      cat("Metadata:          None\n")
    } else {
      cat("Metadata samples:  ", nrow(metadata), "\n")
      if ("Group" %in% colnames(metadata)) {

        group_counts <- table(metadata$Group)

        cat("Groups:            ", length(group_counts), "\n")
        cat("Group sizes:\n")

        for (g in names(group_counts)) {
          cat("  -", g, ":", group_counts[g], "\n")
        }

      } else {
        cat("Groups:            Not found (no 'Group' column)\n")
      }

      matched <- sum(metadata$Sample %in% sample.cols)
      cat("Matched samples:   ", matched, "/", length(sample.cols), "\n")
    }
  })

  # Warnings block (only if exists)
  if (length(warnings_list) > 0) {
    .print_block("Warnings", function() {
      for (w in warnings_list) {
        cat("- ", w, "\n", sep = "")
      }
    })
  }
  if (any(counts_only < 0, na.rm = TRUE)) {
    cat("Negative values detected!\n")
  }

  return(invisible(proj))

}

