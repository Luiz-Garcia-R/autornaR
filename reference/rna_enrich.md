# Gene Ontology (GO) enrichment analysis for differential expression results

Performs Gene Ontology (GO) enrichment analysis on differentially
expressed genes derived from
[`rna_compare()`](https://luiz-garcia-r.github.io/autornaR/reference/rna_compare.md),
using the clusterProfiler framework. Enrichment is computed separately
for up- and down-regulated genes, allowing directional biological
interpretation.

The function automatically retrieves the selected comparison from the
active `rna_project`, applies filtering based on adjusted p-values and
log2 fold-change thresholds, maps gene identifiers, and computes
enrichment across one or more ontology domains (BP, MF, CC).

## Usage

``` r
rna_enrich(
  project,
  plot_style = c("dotplot", "barplot"),
  ont = c("BP", "MF", "CC"),
  enrich_p_cutoff = 0.1,
  padj_cutoff = 0.05,
  log2fc_cutoff = 0.5,
  plot = TRUE,
  top_terms = 10,
  save = TRUE
)
```

## Arguments

- project:

  `rna_project` object created by
  [`rna_project()`](https://luiz-garcia-r.github.io/autornaR/reference/rna_project.md).

- plot_style:

  Character. Wheter use `barplot` or `dotplot`.

- ont:

  Character vector. Ontology domains to test:

  - `"BP"`: Biological Process

  - `"MF"`: Molecular Function

  - `"CC"`: Cellular Component

- enrich_p_cutoff:

  Numeric. P-value cutoff used internally by
  [`clusterProfiler::enrichGO()`](https://rdrr.io/pkg/clusterProfiler/man/enrichGO.html)
  to filter enriched GO terms (default: `0.1`).

- padj_cutoff:

  Numeric. Adjusted p-value threshold used to filter differentially
  expressed genes prior to enrichment (default: `0.05`).

- log2fc_cutoff:

  Numeric. Absolute log2 fold-change threshold used to define up- and
  down-regulated genes (default: `0.5`).

- plot:

  Logical. Whether to display dot plots of enriched GO terms (default:
  `TRUE`).

- top_terms:

  Integer. Number of top enriched terms to display in plots (default:
  `10`).

- save:

  Logical. Whether to store results in the active `rna_project`
  (default: `TRUE`).

## Value

An object of class `"enrich_result"` containing:

- contrast:

  Group contrast (test vs reference)

- organism:

  Organism used for annotation

- gene_id_type:

  Gene identifier type used in enrichment

- up:

  Named list of enrichment results for up-regulated genes (per ontology)

- down:

  Named list of enrichment results for down-regulated genes (per
  ontology)

- params:

  List of parameters used in the analysis

## Details

The function performs the following steps:

1.  Retrieves differential expression results from `rna_project`

2.  Filters genes using `padj_cutoff` and `log2fc_cutoff`

3.  Splits genes into up- and down-regulated sets

4.  Performs GO enrichment separately for each set and ontology

5.  Optionally visualizes results using dot plots

Gene identifiers are automatically interpreted using the project
metadata (`gene_id_type`) and passed to
[`clusterProfiler::enrichGO()`](https://rdrr.io/pkg/clusterProfiler/man/enrichGO.html).

Supported annotation databases:

- Human: `org.Hs.eg.db`

- Mouse: `org.Mm.eg.db`

- Zebrafish: `org.Dr.eg.db`

**Important:**

- Gene IDs must be cleaned (no Ensembl version suffixes)

- Enrichment depends strongly on filtering thresholds

- Up- and down-regulated enrichments may reveal distinct biology

## Stored output

When `save = TRUE`, results are stored in:
`rna_project$analyses$enrich`, allowing reproducibility and reuse.

## Visualization

When `plot = TRUE`, dot plots are generated using
[`enrichplot::dotplot()`](https://rdrr.io/pkg/enrichplot/man/dotplot.html),
showing the top enriched GO terms per ontology and direction (up/down).

## Examples

``` r
if (FALSE) { # \dontrun{
# Run enrichment using last comparison
rna_enrich(project = my_project)

# Run enrichment with stricter thresholds
rna_enrich(my_project,
           padj_cutoff = 0.01,
           log2fc_cutoff = 1)

# Run enrichment for molecular function
rna_enrich(my_project,
           ont = "MF")
} # }
```
