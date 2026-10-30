# Diagnose soft-thresholding power for WGCNA

Evaluates candidate soft-thresholding powers for weighted gene
co-expression network analysis (WGCNA) using normalized expression data
stored in an `rna_project` object. The function filters genes with zero
or invalid variance, optionally retains the most variable genes,
performs basic expression-data quality control, and evaluates scale-free
topology fit across the candidate powers.

Diagnostic plots and a summary of the evaluated powers are provided to
support the selection of a suitable soft-thresholding power before
network construction.

## Usage

``` r
rna_pick_power(
  project,
  genes_n = 10000,
  powers = 1:20,
  network_type = c("signed", "unsigned", "signed hybrid"),
  fit_threshold = 0.9,
  plot = TRUE,
  save = TRUE
)
```

## Arguments

- project:

  An `rna_project` object containing normalized gene expression data and
  sample metadata.

- genes_n:

  Integer or `NULL`. Maximum number of genes to retain, selected by
  decreasing expression variance across samples. Genes with non-finite
  or zero variance are excluded before selection. If `NULL`, all
  eligible genes are retained. Using a subset of the most variable genes
  can reduce computation time during power evaluation. Defaults to
  `10000`.

- powers:

  Numeric vector of candidate soft-thresholding powers to evaluate.
  Values must be finite, positive, and unique. Defaults to `1:20`.

- network_type:

  Character string specifying the network type used by WGCNA. Supported
  values are `"signed"`, `"unsigned"`, and `"signed hybrid"`. A signed
  network preserves the distinction between positive and negative
  correlations; an unsigned network uses the magnitude of correlations;
  and a signed hybrid network excludes negative correlations from the
  adjacency calculation. Defaults to `"signed"`.

- fit_threshold:

  Numeric value between 0 and 1 defining the minimum signed scale-free
  topology fit index used to flag candidate powers. The signed fit index
  is calculated from the scale-free topology fit and the sign of its
  regression slope. Powers meeting or exceeding this threshold are
  marked in the `meets_threshold` column of the returned diagnostic
  table. This threshold does not automatically select a final power or
  guarantee that the resulting network will be suitable. Defaults to
  `0.90`.

- plot:

  Logical. If `TRUE`, displays diagnostic plots showing the signed
  scale-free topology fit across candidate powers and the corresponding
  mean network connectivity. The horizontal reference line in the fit
  plot corresponds to `fit_threshold`. Defaults to `TRUE`.

- save:

  Logical. If `TRUE`, attaches the diagnostic result to the project
  under the `analyses` slot with subtype `"wgcna_power"` and invisibly
  returns the updated project. If `FALSE`, does not attach the result
  and instead invisibly returns an object of class `power_result`.
  Defaults to `TRUE`.

## Value

If `save = TRUE`, invisibly returns the updated `rna_project` object
containing the diagnostic result.

If `save = FALSE`, invisibly returns an object of class `power_result`,
containing the fit diagnostics, candidate powers, network type, fit
threshold, numbers of genes and samples used, identifiers of retained
genes and samples, identifiers removed during quality control, and the
analysis timestamp.

## Input data

The expression matrix must be numeric, with genes in rows and samples in
columns. Gene identifiers and sample identifiers must be unique. The
metadata must identify samples through a `Sample` column or through its
row names. Sample identifiers are matched to the expression matrix
column names; metadata must contain an entry for every expression
sample.

Expression values must be finite. Genes with non-finite or zero variance
are excluded before the most variable genes are selected. WGCNA quality
control is then applied to identify and remove additional problematic
samples or genes.

## Interpretation

The scale-free topology fit and mean connectivity plots provide
complementary diagnostics. A higher scale-free topology fit may support
a candidate power, but increasing the power can also reduce network
connectivity. Therefore, power selection should consider both
diagnostics and the biological objectives of the analysis.

The `fit_threshold` argument defines a diagnostic criterion rather than
a universal cutoff. If no candidate power reaches the specified
threshold, the function reports this result without automatically
relaxing the criterion. Alternative thresholds may be explored when
justified by the dataset and the intended analysis.

## Scope of the analysis

This function evaluates candidate powers using the samples supplied in
the project as a whole. It does not test differential co-expression
between experimental groups and does not construct the final network or
identify co-expression modules. Group-specific or consensus network
analyses require additional steps.

## References

Langfelder P, Horvath S. (2008). WGCNA: an R package for weighted
correlation network analysis. *BMC Bioinformatics*, 9, 559.
[doi:10.1186/1471-2105-9-559](https://doi.org/10.1186/1471-2105-9-559)

## See also

[`pickSoftThreshold`](https://rdrr.io/pkg/WGCNA/man/pickSoftThreshold.html),
[`goodSamplesGenes`](https://rdrr.io/pkg/WGCNA/man/goodSamplesGenes.html)
