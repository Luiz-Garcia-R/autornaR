# ROC analysis for gene signatures

Computes Receiver Operating Characteristic (ROC) curves and Area Under
the Curve (AUC) for gene-level or multi-gene signatures using normalized
RNA-seq data stored in the active `rna_project`.

This function provides a flexible framework to evaluate the
discriminative performance of individual genes or gene sets in binary
classification problems.

## Usage

``` r
rna_roc(
  project,
  genes = NULL,
  group_col = "Group",
  score_method = c("auto", "pca", "mean"),
  show_gene_names = TRUE,
  line_color = "#1f77b4",
  plot = TRUE,
  save = TRUE
)
```

## Arguments

- project:

  `rna_project` object created by
  [`rna_project()`](https://luiz-garcia-r.github.io/autornaR/reference/rna_project.md).

- genes:

  Character vector of gene symbols or Ensembl IDs, or `"all"` to use all
  genes.

- group_col:

  Character. Column name in metadata defining binary groups.

- score_method:

  Character. One of `"auto"`, `"pca"`, or `"mean"`.

- show_gene_names:

  Logical. Wheter to display gene names in subtitle.

- line_color:

  Character. Optional plot line color (default: `"#1f77b4"`).

- plot:

  Logical. Whether to display the ROC curve.

- save:

  Logical. Whether to store results in the active `rna_project`.

## Value

An object of class `"roc_result"` containing:

- roc_object:

  ROC object from pROC.

- auc:

  Numeric AUC value.

- auc_ci:

  Confidence interval for AUC.

- genes:

  Genes used in the analysis.

- method:

  Scoring strategy used (`"single_gene"`, `"logistic"`, or
  `"signature"`).

- score_method:

  Method used for multi-gene scoring (`"pca"` or `"mean"`).

- best_threshold:

  Optimal classification threshold (Youden index).

- best_sens:

  Sensitivity at optimal threshold.

- best_spec:

  Specificity at optimal threshold.

## Details

The function adapts the prediction strategy based on the number of genes
provided:

**Single gene** Expression values are directly used as predictor scores.

**Small gene sets (\<= 5 genes)** A logistic regression model is fitted
using gene expression as predictors.

**Larger gene sets (\> 5 genes)** A composite signature score is
computed:

- `pca`: first principal component (default in `"auto"` mode)

- `mean`: average expression across genes

The ROC curve is computed using pROC, and the optimal threshold is
determined using the Youden index.

Gene identifiers can be provided as Ensembl IDs or gene symbols. When
possible, gene symbols are retrieved automatically using biomaRt.

## Interpretation guide

- AUC ~ 0.5 indicates no discriminative ability

- AUC 0.7–0.8 indicates moderate performance

- AUC 0.8–0.9 indicates strong performance

- AUC \> 0.9 may indicate very strong signal or potential overfitting

For multi-gene signatures, interpretation depends on how the score is
constructed. PCA-based scores capture dominant variance structure, while
mean-based scores assume consistent directionality across genes.

## Important considerations

**Binary outcome required** The grouping variable must contain exactly
two levels.

**Unsupervised signatures** PCA-based scoring does not use class labels
and may not yield optimal classification.

**Large gene sets** Using `"all"` genes may produce unstable or
uninterpretable results, especially without prior feature selection.

**Data leakage** When evaluating predictive performance, avoid using
genes selected from the same dataset without proper validation.

## Pipeline context

This function is typically used after
[`rna_normalize()`](https://luiz-garcia-r.github.io/autornaR/reference/rna_normalize.md)
and optionally after feature selection steps such as
[`rna_compare()`](https://luiz-garcia-r.github.io/autornaR/reference/rna_compare.md).
It complements visualization tools like
[`rna_heatmap()`](https://luiz-garcia-r.github.io/autornaR/reference/rna_heatmap.md)
by providing quantitative classification metrics.

## Side effects

If `save = TRUE`, results are stored in the `rna_project` under the
`analyses$roc` slot.

## Examples

``` r
if (FALSE) { # \dontrun{
# Single gene ROC
rna_roc(project = my_project,
        genes = "CD8A"
)

# Multi-gene logistic model
rna_roc(my_project,
        genes = c("CD3D", "CD8A", "GZMB")
)

# Signature using PCA score
rna_roc(my_project,
        genes = c("Gene1", "Gene2", "Gene3", "Gene4", "Gene5", "Gene6"),
        score_method = "pca"
)
} # }
```
