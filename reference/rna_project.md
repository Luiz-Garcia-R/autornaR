# Initialize an RNA-seq analysis project container

Creates a structured `rna_project` object that serves as the central
container for all RNA-seq analyses within the autornaR workflow.

The project object is designed to store raw inputs, processed data,
analysis results, and execution logs in a consistent and reproducible
structure, enabling seamless transitions between analysis steps.

Once created, the project is automatically assigned to the global
environment as `rna_project`, becoming the active context for all
subsequent functions (e.g.,
[`rna_import()`](https://luiz-garcia-r.github.io/autornaR/reference/rna_import.md),
[`rna_normalize()`](https://luiz-garcia-r.github.io/autornaR/reference/rna_normalize.md),
[`rna_compare()`](https://luiz-garcia-r.github.io/autornaR/reference/rna_compare.md)).

## Usage

``` r
rna_project(project_name = NULL)
```

## Arguments

- project_name:

  Character. Optional project identifier used for labeling outputs and
  logs.

## Value

Invisibly returns an object of class `"rna_project"`, containing:

- project_info:

  Metadata about the project (name, creation time, preprocessing steps,
  list of downstream analyses).

- input:

  Imported data (e.g., count matrices, annotations).

- data:

  Processed data objects (e.g., normalized expression matrices).

- logs:

  Execution logs for reproducibility and tracking.

- version:

  Package version used to create the project.

## Details

The `rna_project` object acts as a lightweight state manager for the
analysis workflow. Instead of passing data manually between functions,
each step reads from and writes to this centralized structure.

**Design principles:**

- *Reproducibility*: All steps are stored and traceable via logs.

- *Modularity*: Each function updates a specific part of the project.

- *State awareness*: The current analysis state is always explicit.

- *Safety*: Prevents loss of intermediate results across steps.

**Typical workflow:**

1.  Initialize project with `rna_project()`

2.  Import data using
    [`rna_import()`](https://luiz-garcia-r.github.io/autornaR/reference/rna_import.md)

3.  Normalize expression with
    [`rna_normalize()`](https://luiz-garcia-r.github.io/autornaR/reference/rna_normalize.md)

4.  Perform analyses (e.g.,
    [`rna_compare()`](https://luiz-garcia-r.github.io/autornaR/reference/rna_compare.md),
    [`rna_gsea()`](https://luiz-garcia-r.github.io/autornaR/reference/rna_gsea.md))

Calling `rna_project()` again will overwrite the current active project.

## See also

[`rna_import`](https://luiz-garcia-r.github.io/autornaR/reference/rna_import.md),
[`rna_normalize`](https://luiz-garcia-r.github.io/autornaR/reference/rna_normalize.md),
[`rna_qc`](https://luiz-garcia-r.github.io/autornaR/reference/rna_qc.md),
[`rna_compare`](https://luiz-garcia-r.github.io/autornaR/reference/rna_compare.md)

## Examples

``` r
if (FALSE) { # \dontrun{

# Initialize a new project
my_project <- rna_project(project_name = "my_project")

# Continue workflow updating project
my_project <- rna_import(...)
my_project <- rna_normalize(...)
my_project <- rna_compare(...)

} # }
```
