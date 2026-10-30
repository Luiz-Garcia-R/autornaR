# Gene-level comparison and boxplot visualization using limma

Performs differential expression analysis for a single gene using a
linear model fitted with limma, optionally including batch effects. The
function computes log2 fold-change (log2FC), 95% confidence interval,
p-value, and Cohen's d, and generates a combined violin + boxplot
visualization with individual sample points.

Gene identifiers can be provided as Ensembl (with or without version),
ENTREZID, or gene symbol. Annotation is performed automatically using
organism-specific org.\*.eg.db databases.

## Usage

``` r
rna_boxplot(
  project,
  gene,
  group_col = "Group",
  batch_col = NULL,
  style = c("boxplot", "violin"),
  colors = NULL,
  show_signif = TRUE,
  seed = NULL,
  save = TRUE,
  verbose = TRUE
)
```

## Arguments

- project:

  `rna_project` object created by
  [`rna_project()`](https://luiz-garcia-r.github.io/autornaR/reference/rna_project.md).

- gene:

  Character. Single gene identifier (Ensembl ID, ENTREZID, or SYMBOL).

- group_col:

  Character. metadata column defining experimental groups.

- batch_col:

  Character. Optional metadata column defining batch effects to include
  in the design matrix.

- style:

  Character. Whether use `boxplot` or `violin` for data visualization.

- colors:

  Character. Optional box colors.

- show_signif:

  Logical. Whether to display significance annotation (asterisks based
  on p-value) on the plot.

- seed:

  Optional. Set a seed for reproducibility for point `jitter`.

- save:

  Logical. Whether to store results in `rna_project$analyses$boxplot`.

- verbose:

  Logical. Whether to print informative messages.

## Value

Invisibly returns a list containing:

- `gene`: Gene identifier used in the expression matrix.

- `gene_label`: Gene symbol (if available).

- `contrast`: Comparison vector (`test`, `reference`).

- `groups`: Named list with `test` and `reference`.

- `statistics`: List with:

  - `logFC`: log2 fold-change (test − reference)

  - `CI_low`, `CI_high`: 95% confidence interval

  - `p_value`: p-value from limma

  - `cohen_d`: Effect size

  - `scale`: Expression scale used

- `design`: Model formula used

- `method`: Statistical method ("limma")

This function is typically used after
[`rna_normalize()`](https://luiz-garcia-r.github.io/autornaR/reference/rna_normalize.md)
and before downstream analyses such as
[`rna_compare()`](https://luiz-garcia-r.github.io/autornaR/reference/rna_compare.md),
or clustering.

## Details

The contrast is always defined as **test − reference**, which determines
the sign of the log2 fold-change. However, for visualization purposes,
groups are displayed in the order **reference → test** (left to right)
to improve interpretability.

The plot combines:

- Violin plot (distribution)

- Boxplot (median and IQR)

- Jittered points (individual samples)

- Optional significance annotation

Results can be stored in the active `rna_project` for reproducibility.

## Examples
