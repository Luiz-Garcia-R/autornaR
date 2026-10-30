# RNA-seq Quality Control (QC)

Performs comprehensive quality control (QC) on normalized RNA-seq data
stored in the active `rna_project`. This function provides a structured
overview of sample quality, expression distributions, and overall
dataset consistency, helping identify technical artifacts and potential
outliers before downstream analyses.

## Usage

``` r
rna_qc(
  project,
  group_col = "Group",
  display_corr_values = FALSE,
  highlight_outliers = TRUE,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  plot = TRUE,
  save = TRUE
)
```

## Arguments

- project:

  `rna_project` object created by
  [`rna_project()`](https://luiz-garcia-r.github.io/autornaR/reference/rna_project.md).

- group_col:

  Character. Column in metadata indicating biological groups (default:
  `"Group"`).

- display_corr_values:

  Logical. If `TRUE`, overlays correlation coefficients on the heatmap
  (default: `FALSE`).

- highlight_outliers:

  Logical. If `TRUE`, detects and highlights outlier samples based on
  library size (default: `TRUE`).

- cluster_rows:

  Logical. Whether to cluster rows in the sample correlation heatmap
  (default: `FALSE`).

- cluster_cols:

  Logical. Whether to cluster columns in the sample correlation heatmap
  (default: `FALSE`).

- plot:

  Logical. If `TRUE`, generates diagnostic plots (default: `TRUE`).

- save:

  Logical. Whether to store results in the active `rna_project`
  (default: `TRUE`).

## Value

An object of class `"rna_QC"` containing:

- data_ready:

  Expression matrix used for QC.

- corr_matrix_samples:

  Sample-wise correlation matrix.

- plot_density:

  Density plot of expression values (if `plot = TRUE`).

- plot_expr_boxplot:

  Boxplot of expression distributions per sample.

- plot_lib_boxplot:

  Barplot of library sizes per sample.

- summary:

  List with basic QC metrics such as library size, genes detected, and
  identified outliers.

- qc_metrics:

  Detailed QC metrics including library variation, correlation
  structure, and classification of dataset quality.

## Details

This function evaluates RNA-seq data quality through multiple
complementary perspectives:

**1. Library size assessment**

- Computes total counts per sample

- Identifies outliers using IQR-based thresholds

- Calculates coefficient of variation (CV) across samples

**2. Expression distribution**

- Boxplots reveal global expression differences between samples

- Density plots highlight distribution shifts across groups

**3. Sample correlation structure**

- Computes pairwise Pearson correlation between samples

- Identifies low-correlation samples

- Compares intra-group vs inter-group similarity

**4. Biological signal evaluation**

- Estimates whether group structure is preserved in the data

- Stronger intra-group than inter-group correlation suggests biological
  signal

**5. Automated QC classification**

- Library quality (based on CV)

- Correlation consistency

- Overall dataset assessment with diagnostic flags

These metrics are summarized into an interpretable QC report, helping
guide decisions such as removing samples, revisiting normalization, or
proceeding with downstream analyses.

## Interpretation guide

- **High-quality data**: low library CV, high sample correlation, clear
  group structure

- **Potential issues**:

  - High library size variation -\> technical variability

  - Low correlation samples -\> possible outliers or batch effects

  - Weak group separation -\> low biological signal or noisy data

QC should always be interpreted in context: biological heterogeneity,
experimental design, and normalization method can all influence these
metrics.

## Pipeline context

This function is typically used after
[`rna_normalize()`](https://luiz-garcia-r.github.io/autornaR/reference/rna_normalize.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic QC
rna_qc(project = my_project)

# QC with custom grouping
rna_qc(my_project,
       group_col = "Condition")

# Display correlation values and highlight structure
rna_qc(my_project,
       display_corr_values = TRUE)

# Disable plotting (metrics only)
rna_qc(my_project,
       plot = FALSE)

} # }
```
