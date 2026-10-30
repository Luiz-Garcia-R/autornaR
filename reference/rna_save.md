# Save a lightweight RNA project for meta-analysis

Extracts the essential components from a `rna_project` object and saves
a simplified, memory-efficient version as an `.rds` file. The function
removes large intermediate objects (e.g., expression matrices, model
objects, and plots) while preserving key results required for downstream
analyses and cross-study comparisons.

This function is designed as a bridge between individual RNA-seq
analyses and meta-analysis workflows (e.g., integration via
`metaOmicsR`).

## Usage

``` r
rna_save(project, file, compress = "xz", verbose = TRUE)
```

## Arguments

- project:

  A `rna_project` object generated through the standard workflow
  ([`rna_project()`](https://luiz-garcia-r.github.io/autornaR/reference/rna_project.md),
  [`rna_import()`](https://luiz-garcia-r.github.io/autornaR/reference/rna_import.md),
  [`rna_normalize()`](https://luiz-garcia-r.github.io/autornaR/reference/rna_normalize.md),
  [`rna_compare()`](https://luiz-garcia-r.github.io/autornaR/reference/rna_compare.md),
  etc.).

- file:

  Character. Output file path where the `.rds` object will be saved.

- compress:

  Character. Compression method passed to
  [`saveRDS()`](https://rdrr.io/r/base/readRDS.html) (default: `"xz"`).

- verbose:

  Logical. If `TRUE`, prints a message upon successful save (default:
  `TRUE`).

## Value

Invisibly returns a cleaned `rna_project` object containing:

- study_info:

  General study metadata (organism, platform, sample size, gene
  identifier type, and contrast definition).

- gene_info:

  Gene annotation table (gene_id, symbol, entrez).

- de_results:

  Differential expression results with standardized column names.

- enrichment:

  Pathway enrichment results (e.g., GSEA output).

- qc:

  Quality control summary extracted from project logs.

- version:

  Version of the originating `rna_project`.

## Details

The function performs the following steps:

**1. Validation**  
Ensures the input object is a valid `rna_project`.

**2. Extraction of latest results**  
For each component (input, normalization, comparison, enrichment, QC),
only the most recent analysis (identified by the `"last"` pointer) is
retained to avoid redundancy and reduce file size.

**3. Contrast parsing**  
The comparison contrast is parsed to extract:

- `condition_tested` (e.g., treated, disease)

- `reference_condition` (e.g., control, baseline)

This information is critical for consistent interpretation of effect
directions across studies in meta-analysis.

**4. Object simplification**  
Heavy elements such as expression matrices, model objects, and
intermediate structures are removed. Only summary-level data required
for integration are preserved.

**5. Serialization**  
The cleaned object is saved as a compressed `.rds` file, enabling
efficient storage and fast reloading in R environments.

## See also

[`rna_project`](https://luiz-garcia-r.github.io/autornaR/reference/rna_project.md),
[`rna_import`](https://luiz-garcia-r.github.io/autornaR/reference/rna_import.md),
[`rna_compare`](https://luiz-garcia-r.github.io/autornaR/reference/rna_compare.md),
[`rna_gsea`](https://luiz-garcia-r.github.io/autornaR/reference/rna_gsea.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Save a processed project for meta-analysis
rna_save(project = my_project,
         file = "my_project.rds")

# Save with alternative compression
rna_save(my_project,
         file = "my_project.rds",
         compress = "gzip")
} # }
```
