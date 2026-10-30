# Dimensionality reduction for RNA-seq data (PCA, UMAP, t-SNE)

Performs dimensionality reduction on normalized RNA-seq expression data
stored in the active `rna_project`, generating low-dimensional
representations of samples using Principal Component Analysis (PCA),
Uniform Manifold Approximation and Projection (UMAP), and optionally
t-distributed Stochastic Neighbor Embedding (t-SNE).

The function automatically adapts parameters based on dataset size and
structure, and returns both numerical results and publication-ready
plots.

## Usage

``` r
rna_dimred(
  project,
  group_col = "Group",
  run_tsne = FALSE,
  add_ellipse = FALSE,
  ellipse_level = 0.95,
  colors = NULL,
  ncomp_pca = NULL,
  n_neighbors = NULL,
  tsne_perplexity = 30,
  use_pca_preprocessing = TRUE,
  pca_method = "auto",
  pca_var_threshold = 0.8,
  pca_max_dims = 50,
  pca_min_dims = 10,
  point_size = 2,
  seed = NULL,
  save = TRUE,
  verbose = TRUE
)
```

## Arguments

- project:

  `rna_project` object created by
  [`rna_project()`](https://luiz-garcia-r.github.io/autornaR/reference/rna_project.md).

- group_col:

  Character. Metadata column used to define sample groups for coloring
  plots (default: `"Group"`).

- run_tsne:

  Logical. Whether to compute t-SNE embeddings (default: `FALSE`).

- add_ellipse:

  Logical. If `TRUE`, adds confidence ellipses to PCA plot when group
  sizes allow (default: `FALSE`).

- ellipse_level:

  Numeric. Confidence level for PCA ellipses (default: `0.95`).

- colors:

  Character. Optional point colors.

- ncomp_pca:

  Integer. Maximum number of PCA components to compute. If `NULL`, uses
  the maximum possible (default).

- n_neighbors:

  Integer. Number of neighbors for UMAP. If `NULL`, automatically set
  based on sample size.

- tsne_perplexity:

  Numeric. Initial perplexity for t-SNE (automatically adjusted if too
  large; default: `30`).

- use_pca_preprocessing:

  Logical. Whether to use PCA-reduced space as input for UMAP and t-SNE
  (default: `TRUE`).

- pca_method:

  Character. Strategy for selecting number of PCs:

  - `"auto"`: selects PCs based on cumulative variance threshold

  - `"variance"`: selects PCs reaching `pca_var_threshold`

  - `"fixed"`: uses `ncomp_pca` directly

- pca_var_threshold:

  Numeric. Target cumulative variance explained when using adaptive PCA
  selection (default: `0.8`).

- pca_max_dims:

  Integer. Maximum number of PCs allowed (default: `50`).

- pca_min_dims:

  Integer. Minimum number of PCs enforced (default: `10`).

- point_size:

  Numeric. Size of points in plots (default: `2`).

- seed:

  Optional integer used to set the random seed for reproducibility. If
  `NULL`, the existing random number generator state is used.

- save:

  Logical. Whether to store results in the active `rna_project`
  (default: `TRUE`).

- verbose:

  Logical. If `TRUE`, prints progress messages and displays plots
  (default: `TRUE`).

## Value

An object of class `"dimred_result"` containing:

- res_dimred:

  List with dimensionality reduction results:

  - `summary`: dataset overview (samples, genes, grouping, seed and rng
    state)

  - `preprocessing`: PCA preprocessing parameters

  - `PCA`: PCA model, variance explained, loadings, coordinates

  - `UMAP`: UMAP model, coordinates, and parameters

  - `tSNE`: t-SNE model and coordinates (if computed)

- plots:

  List of `ggplot2` objects: `PCA`, `UMAP`, and optionally `tSNE`

## Details

The function expects normalized expression data available in
`rna_project$data$normalized_data`, including:

- `expr_matrix`: gene expression matrix (genes x samples)

- `metadata`: sample metadata

Internally, the expression matrix is transposed so that samples are
treated as observations.

**Adaptive behavior:**

- UMAP neighbors are set to `floor(n_samples / 2)` when unspecified

- t-SNE perplexity is capped at `floor((n_samples - 1)/3)`

- t-SNE is skipped when `n_samples < 6`

- PCA dimensionality is automatically selected when
  `use_pca_preprocessing = TRUE`

**Interpretation notes:**

- PCA captures global variance structure and is deterministic

- UMAP emphasizes local structure and clustering

- t-SNE is stochastic and sensitive to parameter choice

## Stored output

When `save = TRUE`, results are stored in:
`rna_project$analyses$dimred`, enabling reproducibility and reuse.

## Examples

``` r
if (FALSE) { # \dontrun{
# Run default dimensionality reduction
rna_dimred(project = my_project)

# Run with t-SNE and ellipses
rna_dimred(my_project,
           run_tsne = TRUE,
           add_ellipse = TRUE)
} # }
```
