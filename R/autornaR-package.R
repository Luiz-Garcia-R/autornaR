#' autornaR: Quick Reference and Main Functions for RNA-seq Analysis
#'
#' autornaR provides a comprehensive set of functions for RNA-seq data analysis,
#' from raw count import to normalization, quality control, exploratory analysis,
#' differential expression, and ROC/volcano analyses. This help topic serves as
#' a quick reference and guide for users.
#'
#' ## Main Workflow
#'
#' The recommended workflow guides users from raw RNA-seq count data to QC and
#' exploratory/differential analysis:
#'  - rna.import() – Import and standardize raw RNA-seq count matrices (HISAT2, STAR, featureCounts, tximport, or clean matrices).
#'  - rna.normalize() – Normalize imported data, filter lowly expressed genes, impute missing values, and remove outliers.
#'  - rna.qc() – Generate QC plots including library size, expression distribution, density plots, and correlation heatmaps.
#'
#' ## Main Functions Overview
#' | Function | Description |
#' |----------|-------------|
#' | `rna.import()`      | Import and standardize RNA-seq count data |
#' | `rna.normalize()`   | Normalize, filter, and impute missing values |
#' | `rna.qc()`          | QC plots: library size, expression, density, correlation heatmaps |
#' | `rna.rank()`        | Rank genes by mean or variance per group |
#' | `rna.roc()`         | ROC curves for discriminant genes |
#' | `rna.ttest()`       | T-test or Mann-Whitney for individual genes |
#' | `rna.volcano()`     | Volcano plot and DE gene summary |
#' | `rna.upset()`       | UpSet plot for presence/absence or DE genes across groups |
#'
#' ## Contact and Contributions
#' For suggestions, bug reports, or contributions, see the
#' [GitHub repository](https://github.com/Luiz-Garcia-R/autornaR).
#'
#' ## Example Workflow
#'
#' A small synthetic RNA-seq dataset with 100 genes and 6 samples (3 per group).
#' This dataset is for demonstration purposes only.
#'
#' @examples
#' \dontrun{
#'
#' # Synthetic gene names and raw counts
#' genes <- paste0("G", sprintf("%03d", 1:100))
#' raw_counts <- data.frame(
#'   Sample = paste0("S", 1:6),
#'   matrix(rpois(100*6, lambda = 50), nrow = 6, ncol = 100, dimnames = list(NULL, genes))
#' )
#'
#' metadata <- data.frame(
#'   Sample = raw_counts$Sample,
#'   Group  = rep(c("Control","Treatment"), each = 3)
#' )
#'
#' # Import
#' imp_data <- rna.import(raw_counts, metadata)
#'
#' # Normalize
#' normalized <- rna.normalize(imp_data)
#'
#' # QC
#' rna.qc(normalized, metadata)
#'
#' # Gene ranking and exploratory analysis
#' rna.rank(normalized, metric = "mean")
#' rna.roc(normalized, genes = c("G001","G002"))
#' rna.ttest(normalized, genes = c("G001","G002"))
#' rna.volcano(normalized)
#' rna.upset(normalized)
#'
#' # Print summary of normalization
#' print(normalized)
#'}
#'
#' @name autornaR
#'
"_PACKAGE"
