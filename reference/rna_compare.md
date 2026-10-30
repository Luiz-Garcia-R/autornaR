# Differential expression analysis using DESeq2 or limma-voom

Performs differential expression analysis between experimental groups
using RNA-seq count data stored in the active `rna_project`. Supports
both DESeq2 and limma (voom) frameworks.

The function handles design matrix construction, contrast specification,
and result extraction, returning standardized outputs across methods.

## Usage

``` r
rna_compare(
  project,
  group_col = "Group",
  batch_col = "Batch",
  method = c("deseq2", "limma"),
  contrast = NULL,
  plot = TRUE,
  save = TRUE,
  store_model = FALSE,
  verbose = TRUE
)
```

## Arguments

- project:

  `rna_project` object created by
  [`rna_project()`](https://luiz-garcia-r.github.io/autornaR/reference/rna_project.md).

- group_col:

  Character. Column name in metadata defining experimental groups
  (default: `"Group"`).

- batch_col:

  Optional character. Column specifying batch effects. If present in
  metadata, it is included in the design formula.

- method:

  Character. Differential expression method: `"deseq2"` or `"limma"`.

- contrast:

  Specifies the comparison(s) to perform. Can be:

  - A character vector of length 2: `c("test", "reference")`

  - A list of such vectors for multiple contrasts

- plot:

  Logical. If `TRUE`, displays diagnostic plots: MA plot (both methods)
  and dispersion plot (DESeq2 only).

- save:

  Logical. Whether to store results in the active `rna_project`.

- store_model:

  Logical. Wheter to store the complete model results in the active
  `rna_project`.

- verbose:

  Logical. Whether to print progress messages.

## Value

Returns (invisibly) either:

- A single object of class `"rnaCompare"` if one contrast is specified

- A named list of such objects if multiple contrasts are provided

Each `"rnaCompare"` object contains:

- `model`: Fitted model object (`DESeqDataSet` or `MArrayLM`)

- `res`: Data frame with differential expression results

- `contrast`: Contrast specification used

- `groups`: List with `test` and `reference`

- `method`: Method used ("deseq2" or "limma")

- `plots`: List of plotting functions (MA and dispersion)

## Details

The function requires raw count data and expects that normalization has
been performed with `rna_normalize(method = "none")` to preserve integer
counts.

**Contrast definition:** All contrasts are defined as **test −
reference**, which determines the sign of the log2 fold-change:

- Positive log2FC -\> higher expression in test group

- Negative log2FC -\> higher expression in reference group

**Method-specific details:**

- DESeq2: Uses negative binomial GLM with dispersion estimation.

- limma: Uses voom transformation followed by linear modeling with
  empirical Bayes moderation.

Multiple contrasts can be specified and will be computed independently.
Each result is stored with a unique identifier in:
`rna_project$analyses$comparison`.

## Output summary

A summary is printed to the console including:

- Total genes tested

- Number of significant genes (adjusted p-value \< 0.05)

- Number of up- and down-regulated genes

- Top-ranked genes

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic usage
rna_compare(my_project,
            contrast = c("treated", "control"))

# Multiple contrasts
rna_compare(my_project,
            contrast = list(
                            c("treated", "control"),
                            c("treated2", "control")
  )
)

# Using limma instead of DESeq2
rna_compare(my_project,
            contrast = c("treated", "control"),
            method = "limma")

# Including batch effect
rna_compare(my_project,
            contrast = c("treated", "control"),
            batch_col = "Batch")
} # }
```
