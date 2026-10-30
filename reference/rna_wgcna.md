# Construct a weighted gene co-expression network and detect modules

Constructs a weighted gene co-expression network using normalized
expression data stored in an `rna_project` object and identifies gene
co-expression modules using the WGCNA framework.

The function filters genes with zero or invalid variance, optionally
retains the most variable genes, performs expression-data quality
control, constructs the network using a predefined soft-thresholding
power, and detects modules through hierarchical clustering of the
topological overlap matrix (TOM).

Module assignments, module eigengenes, network information, and analysis
metadata are organized into a result object. Diagnostic dendrogram plots
can also be generated to visualize the identified modules.

## Usage

``` r
rna_wgcna(
  project,
  power,
  genes_n = 10000,
  network_type = c("signed", "unsigned", "signed hybrid"),
  tom_type = c("signed", "unsigned", "unsigned 2", "signed Nowick", "signed Nowick 2"),
  cor_method = "pearson",
  max_block_size = 20000,
  min_module_size = 30,
  deep_split = 2,
  merge_cut_height = 0.25,
  reassign_threshold = 1e-06,
  numeric_labels = FALSE,
  pam_respects_dendro = FALSE,
  threads = NULL,
  plot = TRUE,
  save = TRUE
)
```

## Arguments

- project:

  An `rna_project` object containing normalized gene expression data and
  sample metadata.

- power:

  Numeric value specifying the soft-thresholding power used to construct
  the network. This value should be selected beforehand using
  [`rna_pick_power`](https://luiz-garcia-r.github.io/autornaR/reference/rna_pick_power.md)
  or another appropriate diagnostic approach. Must be finite and
  positive.

- genes_n:

  Integer or `NULL`. Maximum number of genes to retain, selected by
  decreasing expression variance across samples. Genes with non-finite
  or zero variance are excluded before selection. If `NULL`, all
  eligible genes are retained. Using a subset of the most variable genes
  can reduce computational requirements. Defaults to `10000`.

- network_type:

  Character string specifying the network type used by WGCNA. Supported
  values are `"signed"`, `"unsigned"`, and `"signed hybrid"`. A signed
  network preserves the distinction between positive and negative
  correlations; an unsigned network uses the magnitude of correlations;
  and a signed hybrid network excludes negative correlations from the
  adjacency calculation. Defaults to `"signed"`.

- tom_type:

  Character string specifying the type of topological overlap matrix
  (TOM) used for module detection. Supported values are `"signed"`,
  `"unsigned"`, `"unsigned 2"`, `"signed Nowick"`, and
  `"signed Nowick 2"`. This parameter controls how topological overlap
  is calculated and is independent of the adjacency network type
  specified by `network_type`. Defaults to `"signed"`.

- cor_method:

  Character string specifying the correlation method used during network
  construction. Currently supports `"pearson"`, which is also used in
  downstream calculations of module membership (kME) to ensure
  consistency with the original network. Defaults to `"pearson"`.

- max_block_size:

  Integer specifying the maximum number of genes allowed in each block
  during network construction. When the number of genes exceeds this
  limit, WGCNA divides the data into multiple blocks and constructs
  networks separately. Larger values may increase memory requirements.
  Defaults to `20000`.

- min_module_size:

  Integer specifying the minimum number of genes required to form a
  module during dynamic tree cutting. Smaller modules are generally not
  retained as independent modules. Defaults to `30`.

- deep_split:

  Integer between 0 and 4 controlling the sensitivity of dynamic tree
  cutting. Higher values generally produce more modules by allowing
  finer subdivisions of the gene dendrogram. Higher values correspond to
  more aggressive module detection. Defaults to `2`.

- merge_cut_height:

  Numeric value between 0 and 1 specifying the threshold for merging
  modules based on the dissimilarity between their module eigengenes.
  Lower values result in more conservative merging, whereas higher
  values allow more modules to be merged. Defaults to `0.25`.

- reassign_threshold:

  Numeric value between 0 and 1 specifying the threshold used for gene
  reassignment based on module membership significance during module
  detection. A value of `0` disables reassignment based on this
  threshold. Defaults to `1e-6`.

- numeric_labels:

  Logical. If `TRUE`, modules are represented by numeric labels in the
  network result. If `FALSE`, modules are represented by their assigned
  colors. The returned `genes_modules` table also includes module colors
  regardless of this setting. Defaults to `FALSE`.

- pam_respects_dendro:

  Logical. Specifies whether the partitioning around medoids (PAM) stage
  of module detection must respect the dendrogram produced by
  hierarchical clustering. If `TRUE`, PAM assignments are constrained by
  the dendrogram; if `FALSE`, genes may be assigned to modules
  independently of their position in the dendrogram. Defaults to
  `FALSE`.

- threads:

  Integer or `NULL`. Number of threads used during network construction.
  If `NULL`, the `nThreads` argument is omitted and WGCNA uses its
  default behavior. A positive integer specifies the number of threads
  passed to the network construction procedure. Defaults to `NULL`.

- plot:

  Logical. If `TRUE`, generates diagnostic dendrogram plots showing the
  module color assignments for each network block. Defaults to `TRUE`.

- save:

  Logical. If `TRUE`, attaches the network result to the project under
  the `analyses` slot with subtype `"wgcna_network"` and invisibly
  returns the updated project. If `FALSE`, does not attach the result
  and instead invisibly returns an object of class
  `wgcna_network_result`. Defaults to `TRUE`.

## Value

If `save = TRUE`, invisibly returns the updated `rna_project` object
containing the network result.

If `save = FALSE`, invisibly returns an object of class
`wgcna_network_result`, containing:

- `network`: the complete network object returned by
  [`blockwiseModules`](https://rdrr.io/pkg/WGCNA/man/blockwiseModules.html).

- `genes_modules`: a data frame containing gene identifiers, module
  labels, and module colors.

- `MEs`: the ordered module eigengene matrix.

- Network construction parameters, including the soft-thresholding power
  and network types.

- The number and identifiers of retained genes and samples.

- Identifiers of genes excluded due to invalid variance, gene selection,
  or WGCNA quality control, and samples removed during quality control.

- `timestamp`: the date and time of analysis.

## Input data

The expression matrix must be numeric, with genes in rows and samples in
columns. Gene identifiers and sample identifiers must be unique.

The metadata must identify samples through a `Sample` column or through
its row names. Sample identifiers are matched to the expression matrix
column names; metadata must contain an entry for every expression
sample.

Genes with non-finite or zero variance are excluded before the most
variable genes are selected. WGCNA quality control is then applied to
identify and remove additional problematic samples or genes.

## Network construction

The network is constructed using the adjacency function and topological
overlap matrix defined by `network_type` and `tom_type`, respectively.
The soft-thresholding power specified by `power` is applied to the
adjacency calculation.

When the number of selected genes exceeds `max_block_size`, WGCNA may
divide the data into multiple blocks. Module detection is performed
independently within each block, and the resulting module assignments
are returned in a unified result object.

Modules are identified through hierarchical clustering and dynamic tree
cutting. Closely related modules may subsequently be merged according to
`merge_cut_height`.

## Module eigengenes

Module eigengenes are calculated by WGCNA as part of the module
detection procedure. They represent the first principal component of the
expression profiles of genes assigned to each module and provide a
summarized representation of module-level expression patterns.

The returned eigengene matrix is ordered using
[`orderMEs`](https://rdrr.io/pkg/WGCNA/man/orderMEs.html). Additional
analyses, such as module membership (kME) estimation and associations
between module eigengenes and sample-level phenotypes, can be performed
in subsequent steps.

## Interpretation

Identified modules represent groups of genes with similar expression
patterns across the samples included in the analysis. Module detection
depends on the selected soft-thresholding power, network type, TOM
definition, and module detection parameters.

Module colors are arbitrary identifiers assigned by WGCNA and do not
have intrinsic biological meaning. Biological interpretation requires
additional analyses, such as functional enrichment, module membership
evaluation, and associations with relevant phenotypes.

The resulting modules describe co-expression relationships across the
samples supplied in the project as a whole. They do not directly
represent differential co-expression between experimental groups or
establish causal relationships between genes.

## Scope of the analysis

This function constructs a single co-expression network using all
retained samples. It does not perform differential co-expression
analysis, consensus network construction, module preservation analysis,
or phenotype association testing. These analyses require additional
steps.

## References

Langfelder P, Horvath S. (2008). WGCNA: an R package for weighted
correlation network analysis. *BMC Bioinformatics*, 9, 559.
[doi:10.1186/1471-2105-9-559](https://doi.org/10.1186/1471-2105-9-559)

## See also

[`rna_pick_power`](https://luiz-garcia-r.github.io/autornaR/reference/rna_pick_power.md),
[`blockwiseModules`](https://rdrr.io/pkg/WGCNA/man/blockwiseModules.html),
[`goodSamplesGenes`](https://rdrr.io/pkg/WGCNA/man/goodSamplesGenes.html),
[`plotDendroAndColors`](https://rdrr.io/pkg/WGCNA/man/plotDendroAndColors.html),
[`orderMEs`](https://rdrr.io/pkg/WGCNA/man/orderMEs.html)
