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
#'  - rna_import() – Import and standardize raw RNA-seq count matrices (HISAT2, STAR, featureCounts, tximport, or clean matrices).
#'  - rna_normalize() – Normalize imported data, filter lowly expressed genes, impute missing values, and remove outliers.
#'  - rna_qc() – Generate QC plots including library size, expression distribution, density plots, and correlation heatmaps.
#'
#' ## Main Functions Overview
#' | Function | Description |
#' |----------|-------------|
#' | `rna_project()`     | Create a new rna_project object |
#' | `rna_import()`      | Import and standardize RNA-seq count data |
#' | `rna_normalize()`   | Normalize, filter, and impute missing values |
#' | `rna_qc()`          | QC plots: library size, expression, density, correlation heatmaps |
#' | `rna_compare()`     | Performs differential expression analysis between two groups |
#' | `rna_degs()`        | Rank genes by mean or variance per group |
#' | `rna_roc()`         | ROC curves for discriminant genes |
#' | `rna_boxplot()`     | Limma t-test for individual genes |
#' | `rna_volcano()`     | Volcano plot and DE gene summary |
#' | `rna_sets()`        | Evaluatte presence/absence or DE genes across groups |
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
#' # Create a new rna project
#' my_project <- rna_project()
#'
#' # Import
#' my_project <- rna_import(project = my_project,
#'                          raw_counts,
#'                          metadata)
#'
#' # Normalize
#' my_project <- rna_normalize(my_project)
#'
#' # QC
#' my_project <- rna_qc(my_project)
#'
#' # Gene ranking and exploratory analysis
#' rna_compare(my_project,
#'             method = "limma",
#'             contrast = "Treatment", "Control")
#'
#' rna_roc(my_project, genes = c("G001","G002"))
#'
#' rna_boxplot(my_project, genes = "G001")
#'
#' rna_volcano(my_project)
#'
#' rna_sets(my_project)
#'
#' # Save rna_project object
#' rna_save(rna_project)
#'}
#'
#' @name autornaR
#'
"_PACKAGE"
