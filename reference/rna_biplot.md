# PCA Biplot with Feature Selection

Generates a PCA biplot using previously computed dimensionality
reduction results, highlighting genes based on their contribution to
principal components and/or association with a grouping variable.

## Usage

``` r
rna_biplot(
  project,
  use_dimred = NULL,
  group_col = "Group",
  style = c("loading", "anova", "hybrid"),
  colors = NULL,
  n_genes = 10,
  genes = NULL,
  point_size = 2,
  save = TRUE,
  verbose = TRUE
)
```

## Arguments

- project:

  `rna_project` object created by
  [`rna_project()`](https://luiz-garcia-r.github.io/autornaR/reference/rna_project.md).

- use_dimred:

  Optional `dimred_object` to use in `rna_biplot`.

- group_col:

  Character. Column name in metadata used for grouping samples.

- style:

  Character. Method used for gene selection. One of `"loading"`,
  `"anova"`, or `"hybrid"`.

- colors:

  Character. Optional point colors.

- n_genes:

  Integer. Number of top genes to display when `biplot_genes` is not
  provided.

- genes:

  Character vector. Optional vector of gene names to include. Overrides
  automatic selection.

- point_size:

  Numeric. Size of sample points in the plot.

- save:

  Logical. Whether to store results in the active `rna_project`.

- verbose:

  Logical. Whether to print the plot.

## Value

A list of class `biplot_result` containing:

- `biplot`: list with selected genes, coordinates, scores, and PCs used.

- `plot`: ggplot object of the biplot.

## Details

The function supports three modes for gene selection:

- **loading**: selects genes with highest contribution to selected PCs

- **anova**: selects genes most associated with the grouping variable

- **hybrid**: combines loading magnitude and ANOVA significance

Gene selection can be automatic (top-ranked) or user-defined.

The biplot is constructed from a previously computed PCA stored in
`rna_project`. Gene vectors are scaled to match the PCA coordinate
space.

In `hybrid` mode, both loading magnitude and ANOVA-derived significance
are normalized to the range `[0,1]` before combination.

## Examples

``` r
if (FALSE) { # \dontrun{
# Standard PCA biplot
rna_biplot(project = my_project)

# PCA biplot with a specific gene set
rna_biplot(project = my_project,
           n_genes = 2,
           genes = c("IL6", "HSP70"))
} # }
```
