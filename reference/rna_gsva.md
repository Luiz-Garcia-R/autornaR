# Gene Set Variation Analysis (GSVA)

Computes pathway activity scores for individual samples using GSVA-based
methods. Gene sets can be obtained directly from a previous
[`rna_gsea()`](https://luiz-garcia-r.github.io/autornaR/reference/rna_gsea.md)
analysis, MSigDB collections, or user-supplied pathways.

Unlike GSEA, which evaluates pathway enrichment between groups, GSVA
estimates pathway activity for each sample individually.

## Usage

``` r
rna_gsva(
  project,
  source = c("gsea", "msigdb", "custom"),
  pathways = NULL,
  pathway_scope = c("top", "significant", "all"),
  geneset_collection = "H",
  geneset_subcollection = NULL,
  method = c("gsva", "ssgsea", "zscore", "plage"),
  kcdf = "Gaussian",
  min_size = 10,
  max_size = 500,
  save = TRUE
)
```

## Arguments

- project:

  rna_project object.

- source:

  Character. Source of gene sets: `gsea`, `msigdb` or `custom`.

- pathways:

  Optional named list of gene sets when `source = 'custom'`.

- pathway_scope:

  Character. Determines search scope between: `top`, `significant`, or
  `all`

- geneset_collection:

  MSigDB collection when `source = "msigdb"`.

- geneset_subcollection:

  MSigDB subcollection.

- method:

  GSVA method: `gsva`, `ssgsea`, `zscore` or `plage`.

- kcdf:

  Kernel used by GSVA.

- min_size:

  Minimum pathway size.

- max_size:

  Maximum pathway size.

- save:

  Logical. Save result into project.

- gsea:

  Optional gsea_result object when `source = 'gsea'`.

## Value

Object of class `gsva_result`.
