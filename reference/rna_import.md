# Import RNA-seq count data from multiple formats

Imports and standardizes RNA-seq count matrices from common formats
(HISAT2, STAR, featureCounts, tximport, or pre-cleaned matrices).

The function automatically detects gene identifiers (either as a column
or row names), sample columns, and gene ID types. It performs basic
validation and cleaning steps to ensure compatibility with downstream
analysis.

When ENSEMBL identifiers are detected, version numbers are removed and
duplicated IDs are collapsed by summation. Sample names can be
optionally cleaned or manually renamed.

If metadata is provided, the function automatically detects the sample
identifier column and validates its correspondence with the count
matrix. Summary information (e.g., group sizes) is reported when
available.

## Usage

``` r
rna_import(
  project,
  count_matrix,
  metadata = NULL,
  format = c("clean", "hisat2", "star", "featureCounts", "tximport"),
  gene_col = NULL,
  group_col = NULL,
  gene_id_type = c("auto", "ENSEMBL", "SYMBOL", "ENTREZ"),
  organism = c("auto", "mouse", "human", "zebrafish"),
  rename_samples = NULL,
  clean_names = TRUE,
  strict = TRUE,
  save = TRUE,
  envir = parent.frame()
)
```

## Arguments

- project:

  `rna_project` object created by
  [`rna_project()`](https://luiz-garcia-r.github.io/autornaR/reference/rna_project.md).

- count_matrix:

  Data frame containing raw counts and gene identifiers.

- metadata:

  Optional data frame with sample information. The sample column is
  automatically detected based on column names or matching values with
  the count matrix. If detection fails, the function will return an
  error (in strict mode).

- format:

  Character. Input format: `"clean"`, `"hisat2"`, `"star"`,
  `"featureCounts"`, or `"tximport"`.

- gene_col:

  Optional character specifying the column containing gene identifiers.
  If `NULL` (default), the function will attempt to detect the gene
  column automatically. If no suitable column is found, row names will
  be used when possible.

- group_col:

  Optional character specifying the column containing sample groups. If
  `NULL` (default).

- gene_id_type:

  Character. Gene identifier type: `"auto"`, `"ENSEMBL"`, `"SYMBOL"`,
  `"ENTREZ"`. Default is `"auto"`.

- organism:

  Character. Organism name: `"auto"`, `"mouse"`, `"human"`, or
  `"zebrafish"`.

- rename_samples:

  Optional character vector to rename sample columns. Must have the same
  length and order as the detected sample columns.

- clean_names:

  Logical. If `TRUE`, automatically cleans sample names by removing file
  extensions and special characters.

- strict:

  Logical. If `TRUE`, stops execution on validation errors.

- save:

  Logical. Whether to store results in the active `rna_project`.

- envir:

  Environment where to assign the object.

## Value

An object of class `"imp_data"`, containing:

- `data`: Processed count matrix.

- `metadata`: Metadata table (if provided).

- `n_genes`: Number of genes.

- `n_samples`: Number of samples.

- `gene_id_type`: Detected or specified gene ID type.

- `organism`: Organism used.

- `warnings`: Validation warnings generated during import.

- `groups`: (optional) Named list with group sizes if a `Group` column
  is present in metadata.

- `gene_annotation`: Matrix with annotation of all possible genes of the
  provided data set.

- `warnings`: Validation warnings generated during import, including
  automatic detection steps.

## Details

This function serves as the entry point of the RNA-seq analysis
workflow. It is designed to handle heterogeneous input formats with
minimal user intervention by automatically detecting key structural
components such as gene identifiers and sample labels.

The function does not perform normalization or transformation of counts.

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic import with metadata
  rna_import(project = my_project,
             count_matrix = counts_df,
             metadata = meta_df)

# Import with metadata and specifying an organism
  rna_import(my_project,
             count_matrix = counts_df,
             metadata = meta_df,
             organism = "mouse")

# Import with a specific format
  rna_import(my_project,
             count_matrix = counts_df,
             metadata = meta_df,
             format = "star")
 } # }
```
