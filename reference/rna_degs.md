# Extract and visualize top differentially expressed genes (DEGs)

Identifies the top up- and down-regulated genes from a differential
expression analysis performed with
[`rna_compare()`](https://luiz-garcia-r.github.io/autornaR/reference/rna_compare.md),
and optionally generates a bar plot of log2 fold-changes using gene
symbols.

## Usage

``` r
rna_degs(
  project,
  top_n = 10,
  padj_cutoff = 0.05,
  plot = TRUE,
  save = TRUE,
  verbose = TRUE
)
```

## Arguments

- project:

  `rna_project` object created by
  [`rna_project()`](https://luiz-garcia-r.github.io/autornaR/reference/rna_project.md).

- top_n:

  Integer. Number of genes to select from each tail of the log2
  fold-change distribution (default: 10).

- padj_cutoff:

  Numeric. Adjusted p-value threshold used to define significant genes
  (default: 0.05).

- plot:

  Logical. If `TRUE`, produces a bar plot of top genes.

- save:

  Logical. Whether to store results in the active `rna_project`.

- verbose:

  Logical. Whether to print informative messages.

## Value

An object of class `"deg_result"` containing:

- `comparison`: Comparison identifier

- `deg_table`: Full differential expression table

- `top_genes`: Data frame of selected top genes

- `plot`: `ggplot` object (or `NULL`)

- `summary`: Summary statistics including:

  - Total genes

  - Number of up- and down-regulated genes

  - Number of significant genes

  - Distribution metrics (mean, median, extremes)

- `params`: Parameters used (e.g., `top_n`, `timestamp`)

## Details

The function ranks genes based on log2 fold-change and extracts the top
`top_n` most upregulated and downregulated genes.

**comparison_label interpretation:** The log2 fold-change is defined as
**test − reference**, inherited from
[`rna_compare()`](https://luiz-garcia-r.github.io/autornaR/reference/rna_compare.md).
Therefore:

- Positive values indicate higher expression in the test group

- Negative values indicate higher expression in the reference group

The bar plot displays:

- Top upregulated genes (positive log2FC)

- Top downregulated genes (negative log2FC)

- Color-coded regulation direction

## Output summary

A summary is printed to the console including:

- Total genes analyzed

- Number of up- and down-regulated genes

- Directional bias (percentage of upregulated genes)

- Number of significant genes (if available)

- Table of top genes

## Examples

``` r
if (FALSE) { # \dontrun{
# Use last comparison
rna_degs(project = my_project)

# Select more genes
rna_degs(my_project,
         top_n = 20)

# Disable plotting
rna_degs(my_project,
         plot = FALSE)
} # }
```
