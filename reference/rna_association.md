# Association analysis between genes and/or metadata variables

Computes pairwise associations between features derived from normalized
RNA-seq expression data and/or sample metadata, returning both a
structured result object and a heatmap-style visualization.

## Usage

``` r
rna_association(
  project,
  x_features,
  y_features = NULL,
  x_type = c("auto", "gene", "metadata", "gsva"),
  y_type = c("auto", "gene", "metadata", "gsva"),
  format = "continuous",
  group = NULL,
  method = c("spearman", "pearson", "kendall"),
  compute_mi = FALSE,
  mi_method = c("knn", "discrete"),
  cluster = TRUE,
  save = TRUE,
  verbose = TRUE
)
```

## Arguments

- project:

  `rna_project` object created by
  [`rna_project()`](https://luiz-garcia-r.github.io/autornaR/reference/rna_project.md).

- x_features:

  Character vector of features to place on the x-axis. Features may
  represent genes or metadata variables depending on `x_type`.

- y_features:

  Optional character vector of features to place on the y-axis. If
  `NULL`, uses `x_features`.

- x_type:

  Character. Type of x-axis features. Currently supported: `auto`,
  `gene`, `metadata` and `gsva`.

- y_type:

  Character. Type of y-axis features. Currently supported: `auto`,
  `gene` `metadata` and `gsva`.

- format:

  Character. Data format expected for association analysis. Currently
  supports `continuous`.

- group:

  Optional character. Restricts analysis to samples belonging to the
  specified group in metadata.

- method:

  Correlation method used for association analysis. One of `pearson`,
  `spearman`, or `kendall`.

- compute_mi:

  Logical. Whether to additionally compute mutual information for each
  association pair (default: `FALSE`).

- mi_method:

  Method used for mutual information estimation. One of `knn` or
  `discrete`.

- cluster:

  Logical. Whether to hierarchically cluster rows and columns when the
  association matrix is symmetric (default: `TRUE`).

- save:

  Logical. Whether to store results in the active `rna_project`
  (default: `TRUE`).

- verbose:

  Logical. Whether to display the association heatmap and print a
  summary object in the console (default: `TRUE`).

## Value

Invisibly returns either:

- the updated `rna_project` object when `save = TRUE`;

- or an object of class `rna_association` when `save = FALSE`.

The `rna_association` object contains:

- matrix:

  Association matrix.

- long:

  Long-format association table containing correlation coefficients and
  p-values.

- plot:

  Generated ggplot2 heatmap object.

- x_features:

  Features used on the x-axis.

- y_features:

  Features used on the y-axis.

- x_type:

  Type of x-axis features.

- y_type:

  Type of y-axis features.

- params:

  List containing analysis parameters and metadata.

## Details

The function currently supports continuous variables using
correlation-based association metrics and is designed to be extensible
to future multi-omics and mixed-data analyses.

The function was designed as a generalized association engine capable of
integrating transcriptomic and metadata-derived variables within a
unified framework.

Currently, associations are computed using correlation-based statistics
applied to continuous variables.

Gene identifiers may be provided either as gene symbols or Ensembl IDs,
provided that annotation information is available in the project object.

When `cluster = TRUE` and the matrix is symmetric, hierarchical
clustering is applied using:

- distance: `1 - correlation`

- linkage: complete linkage clustering

If `compute_mi = TRUE`, mutual information values are additionally
estimated for each feature pair using either:

- discretization-based estimation;

- K-nearest neighbors estimation.

The resulting heatmap is generated using ggplot2.

Significant associations (`p < 0.05`) are marked with asterisks.

## Examples

``` r
if (FALSE) { # \dontrun{

# Gene-gene association matrix
rna_association(
  project = my_project,
  x_features = c("SERTAD3", "MFSD4B", "RNF139")
)

# Metadata vs genes
rna_association(
  project = my_project,
  x_features = c("Age", "BMI"),
  y_features = c("IL7R", "GZMB"),
  x_type = "metadata",
  y_type = "gene"
)

# Spearman correlation with MI
rna_association(
  project = my_project,
  x_features = c("GeneA", "GeneB"),
  method = "spearman",
  compute_mi = TRUE
)

# Restrict analysis to one group
rna_association(
  project = my_project,
  x_features = c("GeneA", "GeneB"),
  group = "Control"
)
} # }
```
