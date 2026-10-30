# Volcano plot for differential expression results

Generates a volcano plot from differential expression results stored in
the active `rna_project`. The function visualizes the relationship
between effect size (log2 fold-change) and statistical significance
(-log10 p-value), highlighting significantly up- and down-regulated
genes.

## Usage

``` r
rna_volcano(
  project,
  genes = NULL,
  top_genes = NULL,
  padj_threshold = 0.05,
  log2fc_threshold = 0.5,
  colors = NULL,
  save = TRUE,
  verbose = TRUE
)
```

## Arguments

- project:

  `rna_project` object created by
  [`rna_project()`](https://luiz-garcia-r.github.io/autornaR/reference/rna_project.md).

- genes:

  Optional character vector of genes to label on the plot. Gene
  identifiers can be provided either as gene symbols (e.g. `"TP53"`) or
  as gene IDs matching the rownames of the differential expression
  results (e.g. Ensembl IDs). Matching is performed automatically using
  the project annotation, and both representations are supported
  transparently. When provided, this argument takes precedence over
  `top_genes`.

- top_genes:

  Optional integer specifying the number of genes to label automatically
  based on statistical significance and effect size. Genes are ranked
  primarily by adjusted p-value (`padj`) and secondarily by absolute
  log2 fold-change, prioritizing genes that are both highly significant
  and strongly regulated. Ignored if `genes` is provided.

- padj_threshold:

  Numeric. Adjusted p-value cutoff for statistical significance (default
  `0.05`).

- log2fc_threshold:

  Numeric. Absolute log2 fold-change cutoff for biological relevance
  (default `0.5`).

- colors:

  Character. Optional points colors.

- save:

  Logical. Whether to store results in the active `rna_project` (default
  `TRUE`).

- verbose:

  Logical. If `TRUE`, prints the plot and summary to the console
  (default `TRUE`).

## Value

An object of class `"rna_volcano"` containing:

- plot:

  The ggplot object representing the volcano plot.

- thresholds:

  List with `padj` and `log2fc` thresholds used.

- DEGs:

  List with vectors of up- and down-regulated genes.

- comparison:

  Comparison ID used to generate the plot.

## Details

This function retrieves results produced by
[`rna_compare()`](https://luiz-garcia-r.github.io/autornaR/reference/rna_compare.md)
and does not perform differential expression analysis itself.

Genes are classified into three categories:

- Upregulated genes (significant and above log2FC threshold)

- Downregulated genes (significant and below negative log2FC threshold)

- Non-significant genes

Significance is defined using the adjusted p-value threshold
(`padj_threshold`), while biological relevance is controlled via the
log2 fold-change threshold (`log2fc_threshold`).

The function supports optional labeling of genes directly on the plot
using `ggrepel`, which is useful for highlighting key markers.

## **Visualization features**

- Color-coded points (Up, Down, Non-significant)

- Horizontal line indicating adjusted p-value threshold

- Vertical lines indicating log2 fold-change thresholds

- Optional gene labeling for interpretation

## **Interpretation guidelines**

- Genes in the upper left and upper right quadrants are the most
  significant and biologically relevant.

- Symmetry around zero log2FC may indicate balanced regulation across
  conditions.

- Large numbers of significant genes may reflect strong biological
  signal or potential confounding factors (e.g., batch effects).

## **Side effects**

If `save = TRUE`, the result is stored in the active `rna_project` under
`analyses$volcano`, along with metadata describing thresholds and
comparison used.

## See also

[`rna_compare`](https://luiz-garcia-r.github.io/autornaR/reference/rna_compare.md),
[`rna_sets`](https://luiz-garcia-r.github.io/autornaR/reference/rna_sets.md)

## Examples

``` r
if (FALSE) { # \dontrun{

# Generate basic volcano plot
rna_volcano(project = my_project)

# Specific gene labeling
rna_volcano(my_project,
            genes = c('gene1', 'gene2', 'gene3')
)

# Top gene labeling
rna_volcano(my_project,
            top_genes = 10
)

# Custom thresholds
rna_volcano(my_project,
            padj_threshold = 0.01,
            log2fc_threshold = 1
)

} # }
```
