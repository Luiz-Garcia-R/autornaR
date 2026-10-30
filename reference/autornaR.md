# autornaR: Quick Reference and Main Functions for RNA-seq Analysis

autornaR provides a comprehensive set of functions for RNA-seq data
analysis, from raw count import to normalization, quality control,
exploratory analysis, differential expression, and ROC/volcano analyses.
This help topic serves as a quick reference and guide for users.

## Details

### Main Workflow

The recommended workflow guides users from raw RNA-seq count data to QC
and exploratory/differential analysis:

- rna_import() – Import and standardize raw RNA-seq count matrices
  (HISAT2, STAR, featureCounts, tximport, or clean matrices).

- rna_normalize() – Normalize imported data, filter lowly expressed
  genes, impute missing values, and remove outliers.

- rna_qc() – Generate QC plots including library size, expression
  distribution, density plots, and correlation heatmaps.

### Main Functions Overview

|  |  |
|----|----|
| Function | Description |
| [`rna_project()`](https://luiz-garcia-r.github.io/autornaR/reference/rna_project.md) | Create a new rna_project object |
| [`rna_import()`](https://luiz-garcia-r.github.io/autornaR/reference/rna_import.md) | Import and standardize RNA-seq count data |
| [`rna_normalize()`](https://luiz-garcia-r.github.io/autornaR/reference/rna_normalize.md) | Normalize, filter, and impute missing values |
| [`rna_qc()`](https://luiz-garcia-r.github.io/autornaR/reference/rna_qc.md) | QC plots: library size, expression, density, correlation heatmaps |
| [`rna_compare()`](https://luiz-garcia-r.github.io/autornaR/reference/rna_compare.md) | Performs differential expression analysis between two groups |
| [`rna_degs()`](https://luiz-garcia-r.github.io/autornaR/reference/rna_degs.md) | Rank genes by mean or variance per group |
| [`rna_roc()`](https://luiz-garcia-r.github.io/autornaR/reference/rna_roc.md) | ROC curves for discriminant genes |
| [`rna_boxplot()`](https://luiz-garcia-r.github.io/autornaR/reference/rna_boxplot.md) | Limma t-test for individual genes |
| [`rna_volcano()`](https://luiz-garcia-r.github.io/autornaR/reference/rna_volcano.md) | Volcano plot and DE gene summary |
| [`rna_sets()`](https://luiz-garcia-r.github.io/autornaR/reference/rna_sets.md) | Evaluatte presence/absence or DE genes across groups |
| [`rna_wgcna()`](https://luiz-garcia-r.github.io/autornaR/reference/rna_wgcna.md) | Performs a weighted correlation network analysis |

### Contact and Contributions

For suggestions, bug reports, or contributions, see the [GitHub
repository](https://github.com/Luiz-Garcia-R/autornaR).

### Example Workflow

A small synthetic RNA-seq dataset with 100 genes and 6 samples (3 per
group). This dataset is for demonstration purposes only.

## See also

Useful links:

- <https://luiz-garcia-r.github.io/autornaR/>

## Author

**Maintainer**: Luiz Garcia <luiz.cardoso@ufpr.br>
([ORCID](https://orcid.org/0000-0002-9616-0927))

## Examples

``` r
if (FALSE) { # \dontrun{

# Synthetic gene names and raw counts
genes <- paste0("G", sprintf("%03d", 1:100))
raw_counts <- data.frame(
  Sample = paste0("S", 1:6),
  matrix(rpois(100*6, lambda = 50), nrow = 6, ncol = 100, dimnames = list(NULL, genes))
)

metadata <- data.frame(
  Sample = raw_counts$Sample,
  Group  = rep(c("Control","Treatment"), each = 3)
)

# Create a new rna project
my_project <- rna_project()

# Import
my_project <- rna_import(project = my_project,
                         raw_counts,
                         metadata)

# Normalize
my_project <- rna_normalize(my_project)

# QC
my_project <- rna_qc(my_project)

# Gene ranking and exploratory analysis
rna_compare(my_project,
            method = "limma",
            contrast = "Treatment", "Control")

my_project <- rna_roc(my_project, genes = c("G001","G002"))

my_project <- rna_boxplot(my_project, genes = "G001")

my_project <- rna_volcano(my_project)

my_project <- rna_sets(my_project)

my_project <- rna_wgcna(my_project, power = 6)

# Save rna_project object
rna_save(rna_project)
} # }
```
