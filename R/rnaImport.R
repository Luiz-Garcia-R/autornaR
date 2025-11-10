#' Import RNA-seq count data from multiple formats
#'
#' @description
#' Imports and standardizes RNA-seq count matrices from common formats
#' (HISAT2, STAR, featureCounts, tximport, or clean matrices). Automatically
#' detects sample columns, handles metadata, and validates structure.
#'
#' @param df Data frame containing raw counts and gene identifiers.
#' @param metadata Optional data frame with sample information (must include a column named `Sample`).
#' @param format Character. Count file format: `"clean"`, `"hisat2"`, `"star"`, `"featureCounts"`, or `"tximport"`.
#' @param sample.names Optional vector of custom sample names.
#' @param gene.col Character. Column containing gene identifiers (default `"Geneid"`).
#' @param strict Logical. If `TRUE`, stops execution on validation errors (default `TRUE`).
#' @param clean_names Logical. If `TRUE`, automatically cleans sample names by removing paths,
#'   prefixes (e.g., `hisat2_bams.`), and extensions (e.g., `.bam`, `.counts`) (default `TRUE`).
#' @param assign_result Logical. If `TRUE`, assigns the imported object to the environment (default `TRUE`).
#' @param assign_name Character. Name of the assigned object (default `"imp_data"`).
#' @param envir Environment where to assign the object.
#'
#' @return An object of class `"imp_data"` containing:
#'   - `data`: standardized count matrix
#'   - `metadata`: optional metadata
#'   - `warnings`: collected warnings and errors
#'   - `n_genes`: number of genes
#'   - `n_samples`: number of samples
#'   - `detected_format`: format detected/used
#'
#' @examples
#' counts <- data.frame(
#'   Geneid = paste0("Gene", 1:5),
#'   Sample1 = c(100, 50, 300, 0, 120),
#'   Sample2 = c(80, 60, 250, 5, 100),
#'   Sample3 = c(90, 55, 280, 2, 110)
#' )
#'
#' meta <- data.frame(
#'   Sample = c("Sample1", "Sample2", "Sample3"),
#'   Group = c("Control", "Control", "Treatment")
#' )
#'
#' imp <- rna.import(df = counts, format = "featureCounts", metadata = meta)
#' print(imp)
#'
#' @export

rna.import <- function(df,
                       metadata = NULL,
                       format = c("clean", "hisat2", "star", "featureCounts", "tximport"),
                       sample.names = NULL,
                       gene.col = "Geneid",
                       strict = TRUE,
                       clean_names = TRUE,
                       assign_result = TRUE,
                       assign_name = "imp_data",
                       envir = parent.frame()) {

  # --- Setup ---
  format <- match.arg(format)
  errors <- character()
  warnings_list <- character()
  add_error <- function(msg) errors <<- c(errors, msg)
  add_warning <- function(msg) warnings_list <<- c(warnings_list, msg)

  # --- 1. Validate input ---
  if (!is.data.frame(df)) add_error("'df' must be a data.frame with count columns.")

  # --- 1b. Auto-detect gene column name ---
  possible_gene_cols <- c(gene.col,
                          "Geneid","geneid","GeneID","GENEID",
                          "Gene","ID","id","Id",
                          "gene_name","gene","symbol","Symbol")

  found <- intersect(possible_gene_cols, colnames(df))

  if (length(found) == 0) {
    add_error(paste0("No gene column detected. Tried: ",
                     paste(unique(possible_gene_cols), collapse=", ")))
  } else {
    gene.col <- found[1]  # FIRST match wins
  }

  # --- 2. Identify count columns ---
  if (format %in% c("hisat2", "star", "featureCounts")) {
    count.cols <- grep("\\.bam$|\\.counts$", colnames(df), value = TRUE)
    if (length(count.cols) == 0) {
      meta_cols <- c("Start", "End", "Length", "Chr", "Strand")
      count.cols <- setdiff(colnames(df)[sapply(df, is.numeric)], c(gene.col, meta_cols))
    }
    if (length(count.cols) == 0)
      add_error("No count columns found (expected *.bam or *.counts).")
  }

  if (format %in% c("clean", "tximport")) {
    count.cols <- setdiff(colnames(df), gene.col)
    if (length(count.cols) == 0)
      add_error("Clean/tximport format but no count columns found.")
  }

  # --- 3. Build count matrix ---
  if (length(errors) == 0) {
    df_counts <- df[, c(gene.col, count.cols), drop = FALSE]

    # --- Intelligent cleaning of sample names ---
    if (clean_names && format %in% c("hisat2", "star", "featureCounts")) {
      clean_sample_names <- function(x) {
        x <- basename(x)  # remove path
        x <- gsub("^hisat2[_\\.-]*bams[_\\.-]*", "", x, ignore.case = TRUE)
        x <- gsub("^star[_\\.-]*bams[_\\.-]*", "", x, ignore.case = TRUE)
        x <- gsub("\\.(bam|sam|counts|txt|tsv|csv)$", "", x, ignore.case = TRUE)
        x <- gsub("\\.$|_$", "", x)  # trailing . or _
        x <- gsub("[^A-Za-z0-9_\\-]+", "", x)
        x
      }
      colnames(df_counts)[-1] <- clean_sample_names(colnames(df_counts)[-1])
    }
  }

  # --- 4. Adjust sample names ---
  sample.cols <- setdiff(colnames(df_counts), gene.col)
  if (!is.null(sample.names)) {
    if (length(sample.names) == length(sample.cols)) {
      colnames(df_counts)[match(sample.cols, colnames(df_counts))] <- sample.names
      sample.cols <- sample.names
    } else {
      add_warning("Length of `sample.names` does not match the number of sample columns.")
    }
  }

  # --- 5. Metadata validation ---
  if (!is.null(metadata)) {
    if (!is.data.frame(metadata)) add_error("'metadata' must be a data.frame.")
    if (!"Sample" %in% colnames(metadata)) add_error("Metadata must contain a column named 'Sample'.")
    if (!"Group" %in% colnames(metadata)) add_warning("Metadata does not contain a 'Group' column.")

    if ("Sample" %in% colnames(metadata)) {
      missing_samples <- setdiff(metadata$Sample, sample.cols)
      if (length(missing_samples) > 0)
        add_warning(paste("Samples present in metadata but missing in matrix:",
                          paste(missing_samples, collapse = ", ")))
    }
  }

  # --- 6. Handle errors and warnings ---
  if (length(errors) > 0 && strict) {
    stop(
      paste("Validation error(s):\n-",
            paste(errors, collapse = "\n- ")),
      call. = FALSE
    )
  } else if (length(errors) > 0) {
    warnings_list <- c(warnings_list, errors)
  }

  # --- 7. Build final object ---
  obj <- list(
    data = df_counts,
    metadata = metadata,
    warnings = warnings_list,
    n_genes = nrow(df_counts),
    n_samples = length(sample.cols),
    detected_format = format
  )
  class(obj) <- "imp_data"

  # --- 8. Assign to environment ---
  if (assign_result) assign(assign_name, obj, envir = envir)
  return(obj)
}


#' Print Method for imp_data Objects
#'
#' Custom print method for objects of class \code{imp_data}.
#'
#' @param x An object of class \code{imp_data}.
#' @param ... Additional arguments (ignored).
#'
#' @export
print.imp_data <- function(x, ...) {
  msg <- paste0(
    "==============================\n",
    "Object of class 'imp_data'\n",
    "==============================\n",
    "Samples: ", x$n_samples, "\n",
    "Genes:   ", x$n_genes, "\n"
  )

  if (!is.null(x$metadata)) msg <- paste0(msg, "Metadata: ", nrow(x$metadata), " rows\n")

  if (length(x$warnings) > 0) {
    msg <- paste0(
      msg,
      "Warnings (", length(x$warnings), "): ",
      paste(x$warnings, collapse = "; "), "\n"
    )
  }

  msg <- paste0(msg, "==============================")
  message(msg)
  invisible(x)
}
