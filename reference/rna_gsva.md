# Gene Set Variation Analysis (GSVA)

Estimates pathway activity scores for individual samples using gene
expression data and predefined gene sets. Gene sets can be obtained from
a previous
[`rna_gsea()`](https://luiz-garcia-r.github.io/autornaR/reference/rna_gsea.md)
analysis, MSigDB collections, or user-supplied pathways.

Unlike GSEA, which evaluates pathway enrichment between experimental
conditions using a ranked gene list, GSVA-based methods estimate
pathway-level variation across individual samples, allowing downstream
comparisons and longitudinal analyses.

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

  An `rna_project` object containing normalized gene expression data.

- source:

  Character specifying the source of gene sets: `"gsea"` to reuse
  pathways from a previous GSEA analysis, `"msigdb"` to retrieve gene
  sets from MSigDB, or `"custom"` to use user-supplied pathways.

- pathways:

  Optional named list of gene sets, required when `source = "custom"`.
  Each element should contain gene identifiers corresponding to the
  expression matrix.

- pathway_scope:

  Character specifying which pathways to retrieve from a previous GSEA
  result: `"top"`, `"significant"`, or `"all"`. Only applicable when
  `source = "gsea"`.

- geneset_collection:

  Character specifying the MSigDB collection when `source = "msigdb"`
  (e.g., `"H"`, `"C2"`, `"C5"`).

- geneset_subcollection:

  Optional character specifying an MSigDB subcollection (e.g.,
  `"CP:REACTOME"`, `"GO:BP"`).

- method:

  Character specifying the pathway scoring method: `"gsva"`, `"ssgsea"`,
  `"zscore"`, or `"plage"`.

- kcdf:

  Character specifying the kernel distribution used by GSVA to estimate
  expression distributions.

- min_size:

  Integer specifying the minimum number of genes required for a gene set
  to be included in the analysis.

- max_size:

  Integer specifying the maximum number of genes allowed in a gene set.

- save:

  Logical. If `TRUE`, stores the GSVA result in the `rna_project`
  object.

- gsea:

  Optional object of class `gsea_result`, required when
  `source = "gsea"`.

## Value

An object of class `gsva_result` containing the pathway activity scores
and associated analysis information.

## Gene set sources

Gene sets can be obtained from three different sources:

- **GSEA:** reuses pathways from a previous
  [`rna_gsea()`](https://luiz-garcia-r.github.io/autornaR/reference/rna_gsea.md)
  analysis. The `pathway_scope` argument determines whether the
  top-ranked, statistically significant, or all available pathways are
  considered.

      \item \strong{MSigDB:} retrieves gene sets directly from Molecular
        Signatures Database collections, allowing pathway activity estimation
        independently of a previous GSEA analysis.

      \item \strong{Custom:} accepts user-defined gene sets, allowing the
        evaluation of specific biological signatures or experimental gene
        lists.

## Scoring methods

Four methods are available for estimating pathway activity:

- **GSVA:** estimates variation in pathway activity across samples using
  a non-parametric approach based on expression distributions.

      \item \strong{ssGSEA:} calculates enrichment scores for each gene set
        within individual samples based on gene expression rankings.

      \item \strong{z-score:} combines standardized expression values of genes
        within each gene set to estimate pathway-level scores.

      \item \strong{PLAGE:} uses singular value decomposition (SVD) to derive
        pathway activity scores from the expression profiles of genes within
        each set.

These methods differ in their scoring assumptions and numerical
properties. Scores obtained using different methods should not be
assumed to be directly comparable.

## Interpretation

The resulting scores represent relative pathway-level expression
patterns across samples and can be used in downstream analyses,
including differential comparisons, clustering, correlation analysis,
and longitudinal modeling.

Higher scores generally indicate stronger representation of the
corresponding gene expression signature, but do not necessarily
demonstrate increased biochemical activity or functional activation of
the pathway. For PLAGE, score direction is arbitrary and should be
interpreted with respect to the underlying expression patterns.

Results depend on the selected gene sets, gene identifier compatibility,
expression data, and scoring method. Gene sets containing too few or too
many genes after matching to the expression matrix are excluded
according to `min_size` and `max_size`.

## Pipeline context

GSVA is typically performed after expression normalization and before
downstream statistical analyses. Unlike GSEA, it does not require a
previous differential expression comparison and can therefore be used to
investigate pathway-level variation across all samples in a project.

## See also

[`rna_gsea`](https://luiz-garcia-r.github.io/autornaR/reference/rna_gsea.md),
[`rna_compare`](https://luiz-garcia-r.github.io/autornaR/reference/rna_compare.md),
[`gsva`](https://rdrr.io/pkg/GSVA/man/gsva.html)
