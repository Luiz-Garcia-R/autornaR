# Build DEG sets and structured gene membership from differential expression results

Constructs structured gene sets from differential expression results
stored in the active `rna_project`. This function formalizes DEGs as
explicit sets, enabling downstream analyses such as overlaps,
enrichment, and custom visualization.

Optionally, a summary plot of up- and down-regulated genes can be
generated.

## Usage

``` r
rna_sets(
  project,
  padj_cutoff = 0.05,
  log2fc_cutoff = 0.5,
  direction = c("both", "up", "down"),
  colors = NULL,
  plot = FALSE,
  save = TRUE,
  verbose = TRUE
)
```

## Arguments

- project:

  `rna_project` object created by
  [`rna_project()`](https://luiz-garcia-r.github.io/autornaR/reference/rna_project.md).

- padj_cutoff:

  Numeric. Adjusted p-value threshold used to define significance
  (default: `0.05`).

- log2fc_cutoff:

  Numeric. Log2 fold-change threshold for directional classification
  (default: `0.5`).

- direction:

  Character. One of `"both"`, `"up"`, or `"down"`. Determines which
  genes are flagged as DEGs in the main output (default: `"both"`).

- colors:

  Character. Optional bars colors.

- plot:

  Logical. If `TRUE`, displays a barplot summarizing DEG counts per
  group (default: `FALSE`).

- save:

  Logical. Whether to store results in the active `rna_project`
  (default: `TRUE`).

- verbose:

  Logical. If `TRUE`, prints a summary to the console (default: `TRUE`).

## Value

An object of class `"rna_sets"` containing:

- deg:

  Logical vector indicating DEG membership after applying `direction`.

- deg_all:

  Logical vector indicating all DEGs (independent of direction).

- up:

  Logical vector for genes upregulated in the test group.

- down:

  Logical vector for genes upregulated in the reference group.

- one_hot:

  Data frame encoding gene membership across groups (one-hot format).

- params:

  List of parameters and thresholds used to define sets.

## Details

This function transforms differential expression results into explicit
gene sets based on statistical significance and effect size thresholds.

**DEG definition** Genes are considered differentially expressed if:

- `padj < padj_cutoff`

- absolute log2 fold-change exceeds `log2fc_cutoff`

**Directional classification**

- Upregulated genes: `log2FC > log2fc_cutoff`

- Downregulated genes: `log2FC < -log2fc_cutoff`

**Set representation** The output includes multiple complementary
representations:

- Binary DEG membership (`deg`, `deg_all`)

- Directional sets (`up`, `down`)

- One-hot encoded matrix (`one_hot`) for set operations

This structure allows flexible downstream use, including intersections,
enrichment analysis, and visualization.

## Interpretation guide

- The number of DEGs depends strongly on both `padj_cutoff` and
  `log2fc_cutoff`

- Increasing the fold-change threshold yields more conservative,
  high-effect genes

- Directional filtering is useful for asymmetric biological
  interpretations

The distinction between "up" and "down" depends on the contrast
definition used in the original comparison.

## Pipeline context

This function is typically used after
[`rna_compare()`](https://luiz-garcia-r.github.io/autornaR/reference/rna_compare.md)
to formalize DEG sets. The resulting object is designed to be used with
downstream tools such as:

- visualization (e.g., heatmaps, Venn diagrams)

- enrichment analysis

- custom set operations

It represents a key step in separating statistical testing from
biological interpretation.

## Important considerations

**Threshold sensitivity** DEG sets can vary substantially depending on
cutoff choices. Always inspect results across multiple thresholds when
possible.

**Directionality depends on contrast** The meaning of "up" and "down" is
tied to the test vs reference groups defined in the original comparison.

**Gene universe consistency** Only genes present in both the expression
matrix and comparison results are retained.

## Side effects

If `save = TRUE`, results are stored in the `rna_project` under the
`analyses$rna_sets` slot.

## Examples

``` r
if (FALSE) { # \dontrun{
# Build DEG sets from last comparison
rna_sets(project = my_project)

# More stringent thresholds
rna_sets(my_project,
         padj_cutoff = 0.01,
         log2fc_cutoff = 1
)

# Only upregulated genes
rna_sets(my_project,
         direction = "up")

# Plot DEG summary
rna_sets(my_project,
         plot = TRUE)

} # }
```
