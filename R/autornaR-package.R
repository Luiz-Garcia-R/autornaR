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
#' | `rna.project()`     | Create a new rna_project object |
#' | `rna.import()`      | Import and standardize RNA-seq count data |
#' | `rna.normalize()`   | Normalize, filter, and impute missing values |
#' | `rna.qc()`          | QC plots: library size, expression, density, correlation heatmaps |
#' | `rna.compare()`     | Performs differential expression analysis between two groups |
#' | `rna.degs()`        | Rank genes by mean or variance per group |
#' | `rna.roc()`         | ROC curves for discriminant genes |
#' | `rna.boxplot()`     | Limma t-test for individual genes |
#' | `rna.volcano()`     | Volcano plot and DE gene summary |
#' | `rna.sets()`        | Evaluatte presence/absence or DE genes across groups |
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
#' my_project <- rna.project()
#'
#' # Import
#' my_project <- rna.import(project = my_project,
#'                          raw_counts,
#'                          metadata)
#'
#' # Normalize
#' my_project <- rna.normalize(my_project)
#'
#' # QC
#' my_project <- rna.qc(my_project)
#'
#' # Gene ranking and exploratory analysis
#' rna.compare(my_project,
#'             method = "limma",
#'             contrast = "Treatment", "Control")
#'
#' rna.roc(my_project, genes = c("G001","G002"))
#'
#' rna.boxplot(my_project, genes = "G001")
#'
#' rna.volcano(my_project)
#'
#' rna.sets(my_project)
#'
#' # Save rna_project object
#' rna.save(rna_project)
#'}
#'
#' @name autornaR
#'
"_PACKAGE"
