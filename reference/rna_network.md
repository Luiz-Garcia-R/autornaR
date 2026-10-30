# Pathway gene coexpression network from RNA-seq data

Builds a gene coexpression network for genes belonging to a selected
pathway using normalized RNA-seq expression data stored in the active
`rna_project`.

The function translates pathway-level gene sets (from GSEA results) into
a correlation-based network, where edges represent gene-gene
coexpression within a defined experimental group.

Beyond network construction, the function supports two complementary
analytical views:

- **Global structure (`node_filter = "all"`)**: preserves the full
  pathway topology, enabling interpretation of modular organization and
  overall connectivity.

- **Hub-centered structure (`node_filter = "top"`)**: restricts
  visualization to the most central genes, emphasizing regulatory
  influence and core network drivers.

These two modes represent different biological perspectives rather than
simple visual filtering.

## Usage

``` r
rna_network(
  project,
  pathway = NULL,
  group = NULL,
  threshold = 0.8,
  cor_p = 0.05,
  node_filter = c("all", "top"),
  top_nodes = 0.2,
  cor_method = c("pearson", "spearman", "auto"),
  community_method = c("louvain", "leiden", "none"),
  layout = c("fr", "kk", "stress", "lgl", "graphopt"),
  color_by = c("expression", "log2fc"),
  labels = TRUE,
  seed = NULL,
  plot = TRUE,
  save = TRUE
)
```

## Arguments

- project:

  `rna_project` object created by
  [`rna_project()`](https://luiz-garcia-r.github.io/autornaR/reference/rna_project.md).

- pathway:

  Character or numeric. Pathway name or index from GSEA results.

- group:

  Character. Experimental group used to compute correlations.

- threshold:

  Numeric. Minimum absolute correlation required to define an edge.

- cor_p:

  Numeric. Maximum p-value allowed for correlations.

- node_filter:

  Character. Defines the structural view of the network:

  - `"all"`: full pathway network (recommended for structural
    interpretation)

  - `"top"`: restricted network containing only the most central genes

  **Important:** When `node_filter = "top"`, the network is no longer a
  full pathway graph but a *centrality-enriched subnetwork*,
  highlighting influential genes rather than full modular architecture.

- top_nodes:

  Numeric. Fraction of most central nodes retained when
  `node_filter = "top"` (default: `0.2`).

- cor_method:

  Character. Correlation method used:

  - `"pearson"`: linear coexpression structure

  - `"spearman"`: rank-based robustness

  - `"auto"`: adaptive choice based on sample size

- community_method:

  Character. Community detection strategy:

  - `"louvain"`: favors larger, coarse-grained modules (global view)

  - `"leiden"`: more stable and fine-grained community structure

  - `"none"`: no explicit modular partitioning

  **Guidance:**

  - Use `louvain` for full-network interpretation (`all`)

  - Use `leiden` for hub-focused or reduced networks (`top`)

  - Use `none` when focusing on continuous gradients of connectivity

- layout:

  Character. Graph layout algorithm:

  - `"fr"`: fast global layout, good for overview plots

  - `"kk"`: balanced structure-preserving layout

  - `"stress"`: best preservation of relational distances (recommended
    for top nodes)

  - `"lgl"`: scalable layout for larger networks

  - `"graphopt"`: compact structure visualization

  **Guidance:**

  - `fr`: publication-style overview of full networks

  - `stress`: best for interpreting hub-centered structures

  - `kk`: balanced exploratory analysis

- color_by:

  Method used to color network nodes:

  - `expression`: colors nodes according to the z-scored mean expression
    of each gene across samples in the selected group.

  - `log2fc`: colors nodes according to the log2 fold change obtained
    from the current comparison stored in the project.

- labels:

  Logical. Wheter to display gene labels in the network plot.

- seed:

  Optional integer used for reproducibility.

- plot:

  Logical. Whether to display the network plot.

- save:

  Logical. Whether to store results in `rna_project`.

## Value

An object of class `"rna_network"` containing:

- params:

  List of parameters used to build the network

- pathway:

  Selected pathway name

- genes:

  Genes belonging to the pathway

- graph:

  Network object from `igraph`

- edges:

  Edge table with correlation statistics

- nodes:

  Node table containing gene annotations, expression summaries, degree,
  betweenness centrality, harmonic centrality, community assignments,
  hub classification, and the composite hub score, where applicable.

The function invisibly returns a list containing:

- `result`: the network object

- `project`: updated project object (when `save = TRUE`)

## Details

The function expects pre-computed GSEA results and normalized expression
data.

The resulting network is sensitive to filtering choices:

- Lower `threshold` increases connectivity and density

- `node_filter = "top"` shifts interpretation toward central regulators

- Community structure becomes less stable in highly filtered networks

Therefore, parameter selection should reflect the biological question:

**Recommended usage patterns:**

- *Pathway structure exploration*: `all + leiden + stress`

- *Regulatory hub discovery*: `top + none/leiden + stress`

- *Publication figure*: `all + louvain + fr`

## Interpretation note

Networks generated with `node_filter = "top"` represent a
*centrality-enriched projection of the pathway*, not the full
coexpression topology.

## Network centrality

Network centrality is characterized using three complementary metrics:

- **Degree**: number of direct connections, representing local
  connectivity.

- **Betweenness centrality**: measures how frequently a node lies on
  shortest paths between other nodes, highlighting potential bridging
  positions in the network.

- **Harmonic centrality**: quantifies node accessibility using the
  inverse of shortest-path distances to other nodes.

Harmonic centrality is used instead of conventional closeness centrality
because pathway coexpression networks may contain disconnected
components. Unlike conventional closeness centrality, harmonic
centrality accounts for unreachable nodes by assigning them zero
contribution, allowing comparisons across nodes in disconnected
networks.

Shortest-path distances are derived from absolute correlation weights,
such that stronger correlations correspond to shorter distances.

The `hub_score` is calculated as the sum of the standardized degree,
betweenness centrality, and harmonic centrality values. It represents a
composite exploratory measure of network centrality rather than a
validated biological measure of regulatory influence.

## Stored output

When `save = TRUE`, the network object is stored in:

`rna_project$analyses$network`

allowing reproducibility and downstream reuse.

## Examples

``` r
if (FALSE) { # \dontrun{

# Build network for top enriched pathway
rna_network(project = my_project)

# Use a specific pathway from GSEA
rna_network(my_project,
            pathway = 'Cellular Response To Starvation')

# Restrict labels to top 20% genes
rna_network(my_project,
            node_filter = 'top',
            top_nodes = 0.2)

# Use Spearman correlation
rna_network(my_project,
            cor_method = "spearman")

# Change layout algorithm
rna_network(my_project,
            layout = "stress")

} # }
```
