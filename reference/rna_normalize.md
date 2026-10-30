# Normalize RNA-seq count data with integrated QC and filtering

Performs normalization of raw RNA-seq count data stored in the active
`rna_project` object (created with
[`rna_project()`](https://luiz-garcia-r.github.io/autornaR/reference/rna_project.md)
and populated via
[`rna_import()`](https://luiz-garcia-r.github.io/autornaR/reference/rna_import.md)).
The function integrates gene filtering, outlier detection, and multiple
normalization strategies into a single, reproducible workflow.

The resulting normalized dataset is stored internally and used as the
standardized input for downstream analyses (e.g., differential
expression, clustering, and pathway analysis).

## Usage

``` r
rna_normalize(
  project,
  method = c("none", "log2", "cpm", "tpm", "rlog", "vst", "quantile", "upper-quartile"),
  filter_low = TRUE,
  min_expr = 10,
  filter_min_prop = 0.1,
  remove_outlier_samples = TRUE,
  remove_outlier_genes = TRUE,
  outlier_method = c("iqr", "zscore"),
  clean_gene_versions = TRUE,
  gene_id_sep = "\\.",
  save = TRUE
)
```

## Arguments

- project:

  `rna_project` object created by
  [`rna_project()`](https://luiz-garcia-r.github.io/autornaR/reference/rna_project.md).

- method:

  Character. Normalization method to apply:

  - `"none"`: No normalization (raw counts retained)

  - `"log2"`: Log2 transformation (`log2(count + 1)`)

  - `"cpm"`: Counts per million (library size normalization)

  - `"tpm"`: Transcripts per million (length-normalized expression)

  - `"rlog"`: Regularized log transformation (DESeq2)

  - `"vst"`: Variance stabilizing transformation (DESeq2)

  - `"quantile"`: Quantile normalization (limma)

  - `"upper-quartile"`: Upper-quartile normalization

- filter_low:

  Logical. If `TRUE`, removes lowly expressed genes prior to
  normalization (default: `TRUE`).

- min_expr:

  Numeric. Minimum count threshold used to define expression (default:
  `10`).

- filter_min_prop:

  Numeric. Minimum proportion of samples in which a gene must be
  expressed (above `min_expr`) to be retained (default: `0.1`).

- remove_outlier_samples:

  Logical. If `TRUE`, detects and removes samples with abnormal library
  sizes (default: `TRUE`).

- remove_outlier_genes:

  Logical. If `TRUE`, removes genes with extreme variance profiles
  (default: `TRUE`).

- outlier_method:

  Character. Method for outlier detection: `"iqr"` (robust, default) or
  `"zscore"`.

- clean_gene_versions:

  Logical. If `TRUE`, removes version suffixes from gene identifiers
  (e.g., `ENSG000001.5 -> ENSG000001`) when using Ensembl IDs (default:
  `TRUE`).

- gene_id_sep:

  Character. Regular expression defining the separator used for gene
  version removal (default: `"\."`).

- save:

  Logical. Whether to store the normalized data inside the active
  `rna_project` (default: `TRUE`).

## Value

Invisibly returns an object of class `"normalized_data"` containing:

- expr_matrix:

  Normalized expression matrix (genes × samples).

- metadata:

  Aligned sample metadata.

- method:

  Normalization method used.

- removed_genes:

  Number of genes removed during filtering/outlier steps.

- removed_samples:

  Number of samples removed as outliers.

- QC_metrics:

  Data frame with per-sample QC metrics (library size and detected
  genes).

- gene_id_version_cleaned:

  Logical indicating whether gene IDs were cleaned.

- gene_id_type:

  Original gene identifier type (e.g., ENSEMBL, SYMBOL).

## Details

The normalization workflow is structured into sequential steps:

**1. Gene filtering**  
Genes with insufficient expression are removed based on the proportion
of samples exceeding `min_expr`. This reduces noise and improves
statistical stability in downstream analyses.

**2. Outlier detection**  
Outlier samples are identified based on library size distributions,
while outlier genes are detected using variance-based criteria. Two
strategies are available:

- `"iqr"`: Robust, recommended for most datasets

- `"zscore"`: Sensitive to extreme deviations

**3. Normalization**  
The selected normalization method is applied. Methods differ in their
assumptions:

- Count-based scaling (`cpm`, `upper-quartile`)

- Length-aware normalization (`tpm`)

- Variance stabilization (`vst`, `rlog`)

- Distribution alignment (`quantile`)

**4. Quality control metrics**  
Basic QC metrics are computed, including library size and number of
detected genes per sample, enabling rapid inspection of data quality.

**5. Gene ID harmonization**  
When Ensembl IDs include version suffixes, they are optionally removed
to ensure compatibility with annotation databases and downstream tools.
If duplicates arise after cleaning, counts are aggregated by gene.

**TPM requirement**  
TPM normalization requires a `gene_length` column in
`rna_project$input$imp_data$annotation`.

## Filtering strategies and use cases

Filtering parameters can strongly influence downstream analyses. The
default configuration is optimized for general RNA-seq workflows,
particularly differential expression analysis.

However, more permissive filtering may be required when the goal is to
preserve genes expressed in rare or specific cell populations.

**Standard analysis (default)**

- Removes lowly expressed genes to reduce noise

- Recommended for differential expression and global analyses

- Default parameters:

  - `min_expr = 10`

  - `filter_min_prop = 0.1`

  - `remove_outlier_genes = TRUE`

**Deconvolution / gene set scoring**

- Preserves low-abundance and cell-type-specific genes

- Recommended for immune profiling, cell signatures, and scoring methods

- Suggested parameters:

  - `min_expr = 1`

  - `filter_min_prop = 0.05`

  - `remove_outlier_genes = FALSE`

**Important:** More stringent filtering may remove biologically relevant
markers, especially for rare cell populations.

## See also

[`rna_import`](https://luiz-garcia-r.github.io/autornaR/reference/rna_import.md),
[`rna_project`](https://luiz-garcia-r.github.io/autornaR/reference/rna_project.md),
[`rna_qc`](https://luiz-garcia-r.github.io/autornaR/reference/rna_qc.md)

## Examples
