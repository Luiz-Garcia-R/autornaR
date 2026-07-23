#' Pathway gene coexpression network from RNA-seq data
#'
#' @description
#' Builds a gene coexpression network for genes belonging to a selected pathway
#' using normalized RNA-seq expression data stored in the active \code{rna_project}.
#'
#' The function translates pathway-level gene sets (from GSEA results) into a
#' correlation-based network, where edges represent gene-gene coexpression
#' within a defined experimental group.
#'
#' Beyond network construction, the function supports two complementary
#' analytical views:
#'
#' \itemize{
#'   \item \strong{Global structure (\code{node_filter = "all"})}:
#'   preserves the full pathway topology, enabling interpretation of
#'   modular organization and overall connectivity.
#'
#'   \item \strong{Hub-centered structure (\code{node_filter = "top"})}:
#'   restricts visualization to the most central genes, emphasizing
#'   regulatory influence and core network drivers.
#' }
#'
#' These two modes represent different biological perspectives rather than
#' simple visual filtering.
#'
#' @param project \code{rna_project} object created by \code{rna.project()}.
#' @param pathway Character or numeric. Pathway name or index from GSEA results.
#' @param group Character. Experimental group used to compute correlations.
#' @param threshold Numeric. Minimum absolute correlation required to define an edge.
#' @param cor_p Numeric. Maximum p-value allowed for correlations.
#' @param node_filter Character. Defines the structural view of the network:
#' \itemize{
#'   \item \code{"all"}: full pathway network (recommended for structural interpretation)
#'   \item \code{"top"}: restricted network containing only the most central genes
#' }
#'
#' \strong{Important:}
#' When \code{node_filter = "top"}, the network is no longer a full pathway graph
#' but a \emph{centrality-enriched subnetwork}, highlighting influential genes
#' rather than full modular architecture.
#'
#' @param top_nodes Numeric. Fraction of most central nodes retained when
#' \code{node_filter = "top"} (default: \code{0.2}).
#'
#' @param cor_method Character. Correlation method used:
#' \itemize{
#'   \item \code{"pearson"}: linear coexpression structure
#'   \item \code{"spearman"}: rank-based robustness
#'   \item \code{"auto"}: adaptive choice based on sample size
#' }
#'
#' @param community_method Character. Community detection strategy:
#' \itemize{
#'   \item \code{"louvain"}: favors larger, coarse-grained modules (global view)
#'   \item \code{"leiden"}: more stable and fine-grained community structure
#'   \item \code{"none"}: no explicit modular partitioning
#' }
#'
#' \strong{Guidance:}
#' \itemize{
#'   \item Use \code{louvain} for full-network interpretation (\code{all})
#'   \item Use \code{leiden} for hub-focused or reduced networks (\code{top})
#'   \item Use \code{none} when focusing on continuous gradients of connectivity
#' }
#'
#' @param layout Character. Graph layout algorithm:
#' \itemize{
#'   \item \code{"fr"}: fast global layout, good for overview plots
#'   \item \code{"kk"}: balanced structure-preserving layout
#'   \item \code{"stress"}: best preservation of relational distances (recommended for top nodes)
#'   \item \code{"lgl"}: scalable layout for larger networks
#'   \item \code{"graphopt"}: compact structure visualization
#' }
#'
#' \strong{Guidance:}
#' \itemize{
#'   \item \code{fr}: publication-style overview of full networks
#'   \item \code{stress}: best for interpreting hub-centered structures
#'   \item \code{kk}: balanced exploratory analysis
#' }
#'
#' @param seed Optional integer used for reproducibility.
#' @param plot Logical. Whether to display the network plot.
#' @param save Logical. Whether to store results in \code{rna_project}.
#'
#' @details
#' The function expects pre-computed GSEA results and normalized expression data.
#'
#' The resulting network is sensitive to filtering choices:
#'
#' \itemize{
#'   \item Lower \code{threshold} increases connectivity and density
#'   \item \code{node_filter = "top"} shifts interpretation toward central regulators
#'   \item Community structure becomes less stable in highly filtered networks
#' }
#'
#' Therefore, parameter selection should reflect the biological question:
#'
#' \strong{Recommended usage patterns:}
#' \itemize{
#'   \item \emph{Pathway structure exploration}: \code{all + leiden + stress}
#'   \item \emph{Regulatory hub discovery}: \code{top + none/leiden + stress}
#'   \item \emph{Publication figure}: \code{all + louvain + fr}
#' }
#'
#' @section Interpretation note:
#' Networks generated with \code{node_filter = "top"} represent a
#' \emph{centrality-enriched projection of the pathway}, not the full
#' coexpression topology.
#'
#' @return
#' An object of class \code{"rna_network"} containing:
#' \describe{
#'   \item{params}{List of parameters used to build the network}
#'   \item{pathway}{Selected pathway name}
#'   \item{genes}{Genes belonging to the pathway}
#'   \item{graph}{Network object from \code{igraph}}
#'   \item{edges}{Edge table with correlation statistics}
#'   \item{nodes}{Node table with gene annotations and network metrics}
#' }
#'
#' The function invisibly returns a list containing:
#' \itemize{
#'   \item \code{result}: the network object
#'   \item \code{project}: updated project object (when \code{save = TRUE})
#' }
#'
#' @section Stored output:
#' When \code{save = TRUE}, the network object is stored in:
#'
#' \code{rna_project$analyses$network}
#'
#' allowing reproducibility and downstream reuse.
#'
#' @examples
#' \dontrun{
#'
#' # Build network for top enriched pathway
#' rna.network(project = my_project)
#'
#' # Use a specific pathway from GSEA
#' rna.network(my_project,
#'             pathway = 'Cellular Response To Starvation')
#'
#' # Restrict labels to top 20% genes
#' rna.network(my_project,
#'             node_filter = 'top',
#'             top_nodes = 0.2)
#'
#' # Use Spearman correlation
#' rna.network(my_project,
#'             cor_method = "spearman")
#'
#' # Change layout algorithm
#' rna.network(my_project,
#'             layout = "stress")
#'
#' }
#'
#' @importFrom stats cor pt
#' @importFrom igraph graph_from_data_frame degree betweenness closeness
#' @importFrom igraph cluster_louvain cluster_leiden membership edge_density
#' @importFrom tidygraph as_tbl_graph
#' @importFrom ggraph ggraph geom_edge_link geom_node_point geom_node_text scale_edge_width
#' @importFrom ggplot2 aes scale_size scale_colour_gradient2 theme_void ggtitle
#' @importFrom ggrepel geom_text_repel
#'
#' @export

rna.network <- function(project,
                        pathway = NULL,
                        group = NULL,
                        threshold = 0.8,
                        cor_p = 0.05,
                        node_filter = c("all", "top"),
                        top_nodes = 0.2,
                        cor_method = c("pearson","spearman","auto"),
                        community_method = c("louvain", "leiden", "none"),
                        layout = c("fr","kk","stress","lgl","graphopt"),
                        seed = NULL,
                        plot = TRUE,
                        save = TRUE) {

  cor_method <- match.arg(cor_method)
  node_filter <- match.arg(node_filter)
  community_method <- match.arg(community_method)
  layout <- match.arg(layout)

  # ===========================================================================
  # 0) Basic checks
  # ===========================================================================
  pkgs <- c(
    "igraph",
    "tidygraph",
    "ggraph",
    "ggrepel",
    "ggplot2"
  )

  .check_dependencies(pkgs)

  # --- Set seed ---
  old_seed <- .set_seed(seed)
  on.exit(.reset_seed(old_seed), add = TRUE)

  # ===========================================================================
  # 1) Get active project
  # ===========================================================================
  proj <- project

  expr_mat <- as.matrix(.get_expr(proj))
  meta <- .get_meta(proj)

  organism <- .get_organism(proj)
  gene_id_type <- .get_gene_id_type(proj)
  gsea_obj <- .get_gsea(proj)
  comp_obj <- .get_comp_obj(proj)

  # ===========================================================================
  # 2) Validate input
  # ===========================================================================
  # Validate layout
  valid_layouts <- c("fr","kk","stress","lgl","graphopt")

  if (!layout %in% valid_layouts) {
    stop("Invalid layout. Choose one of: ",
         paste(valid_layouts, collapse = ", "))
  }

  if (gene_id_type != "ENSEMBL") {
    warning("Network currently assumes ENSEMBL gene IDs.")
  }

  if (is.null(gsea_obj)) {
    stop("No GSEA result found in project. Run rna.gsea() first.")
  }

  if (is.null(comp_obj)) {
    warning("No comparison results found.")
  }

  # ===========================================================================
  # 3) Select pathway
  # ===========================================================================
  pathways <- gsea_obj$pathways
  gsea_top <- gsea_obj$gsea_top

  if (is.null(pathway)) {

    pathway <- gsea_top$pathway[1]
    message("[rna.network] Using top pathway: ", pathway)

  }

  else if (is.numeric(pathway)) {

    pathway <- gsea_top$pathway[pathway]

  }

  if (!pathway %in% names(pathways)) {
    stop("Pathway not found in GSEA results.")
  }

  # ===========================================================================
  # 4) Extract genes
  # ===========================================================================
  genes_path <- pathways[[pathway]]

  # convert SYMBOL -> ENSEMBL if necessary
  if (gene_id_type == "ENSEMBL" && !grepl("^ENSG", genes_path[1])) {

    genes_path <- .convert_gene_ids(
      genes = genes_path,
      from = "SYMBOL",
      to = "ENSEMBL",
      organism = organism
    )
  }

  genes_path <- intersect(genes_path, rownames(expr_mat))

  if (length(genes_path) < 3) {

    stop(
      "Too few genes from pathway present in expression matrix.\n",
      "Possible gene ID mismatch (SYMBOL vs ENSEMBL)."
    )
  }

  # ===========================================================================
  # 5) Select group
  # ===========================================================================
  if (is.null(group)) {

    group <- unique(meta$Group)[1]
    message("[rna.network] Using group: ", group)

  }

  samples_group <- meta$Sample[meta$Group == group]

  expr_sub <- expr_mat[genes_path, samples_group, drop = FALSE]

  # ===========================================================================
  # 6) Correlation
  # ===========================================================================
  edges <- weight <- layout_weight <- NULL

  n <- ncol(expr_sub)

  if (cor_method == "auto") {
    cor_method <- if (n < 8) "spearman" else "pearson"
  }

  message("[rna.network] Correlation method: ", cor_method)

  cor_mat <- cor(t(expr_sub), method = cor_method)

  # Calculate correlation p-value
  t_stat <- cor_mat * sqrt((n - 2) / (1 - cor_mat^2))
  p_mat <- 2 * pt(-abs(t_stat), df = n - 2)

  sel <- which(
    abs(cor_mat) >= threshold &
      p_mat < cor_p &
      upper.tri(cor_mat),
    arr.ind = TRUE
  )

  cor_mat[cor_mat > 0.999999] <- 0.999999
  cor_mat[cor_mat < -0.999999] <- -0.999999

  edges_df <- data.frame(
    source = rownames(cor_mat)[sel[,1]],
    target = rownames(cor_mat)[sel[,2]],
    weight = cor_mat[sel],
    pvalue = p_mat[sel],
    stringsAsFactors = FALSE
  )

  edges_df$distance <- (1 - abs(edges_df$weight)) + 1e-6
  edges_df$sign <- ifelse(edges_df$weight > 0, "positive", "negative")

  # ===========================================================================
  # 7) get nodes
  # ===========================================================================
  nodes_df <- data.frame(
    name = rownames(expr_sub),
    stringsAsFactors = FALSE
  )

  symbols <- .convert_gene_ids(
    genes = nodes_df$name,
    from = "ENSEMBL",
    to = "SYMBOL",
    organism = organism
  )

  if (length(symbols) == length(nodes_df$name)) {
    nodes_df$symbol <- symbols
  } else {
    nodes_df$symbol <- nodes_df$name
  }

  # Nodes color
  nodes_df$expr_mean <- rowMeans(expr_sub)
  nodes_df$expr_z <- as.numeric(scale(nodes_df$expr_mean))


  # ===========================================================================
  # 8) Build network
  # ===========================================================================
  g <- igraph::graph_from_data_frame(
    d = if (nrow(edges_df) == 0) NULL else edges_df,
    vertices = nodes_df,
    directed = FALSE
  )

  if (nrow(edges_df) > 0) {
    igraph::E(g)$distance <- edges_df$distance
    igraph::E(g)$weight <- edges_df$weight
    igraph::E(g)$layout_weight <- abs(edges_df$weight)
  }

  if (igraph::ecount(g) == 0) {
    warning("Network contains no edges at this threshold.")
  }

  # ===========================================================================
  # 8.1 Network metrics
  # ===========================================================================
  nodes_df$degree <- igraph::degree(g)

  nodes_df$betweenness <- igraph::betweenness(
    g,
    weights = igraph::E(g)$distance,
    normalized = TRUE
  )

  nodes_df$closeness <- igraph::closeness(
    g,
    weights = igraph::E(g)$distance,
    normalized = TRUE
  )

  # Composed hub score
  nodes_df$hub_score <- scale(nodes_df$degree) +
    scale(nodes_df$betweenness) +
    scale(nodes_df$closeness)

  nodes_df$hub_score <- as.numeric(nodes_df$hub_score)

  nodes_df$expr_mean <- rowMeans(expr_sub)
  nodes_df$expr_z <- as.numeric(scale(nodes_df$expr_mean))

  # hub detection (top 10% degree)
  deg_cut <- mean(nodes_df$degree) + sd(nodes_df$degree)
  nodes_df$hub <- nodes_df$degree >= deg_cut

  nodes_df$label <- nodes_df$symbol

  igraph::vertex_attr(g, "degree") <- nodes_df$degree
  igraph::vertex_attr(g, "betweenness") <- nodes_df$betweenness
  igraph::vertex_attr(g, "closeness") <- nodes_df$closeness
  igraph::vertex_attr(g, "hub") <- nodes_df$hub
  igraph::vertex_attr(g, "label") <- nodes_df$label

  tg <- tidygraph::as_tbl_graph(g)

  # ===========================================================================
  # 8.2 Node filter
  # ===========================================================================
  if (node_filter == "top") {

    nodes_df$importance <- nodes_df$hub_score
    cutoff <- quantile(nodes_df$importance, 1 - top_nodes, na.rm = TRUE)

    keep_nodes <- nodes_df$name[nodes_df$importance >= cutoff]
    keep_nodes <- intersect(keep_nodes, igraph::V(g)$name)

    if (length(keep_nodes) < 2) {
      stop("Too few nodes after filtering.")
    }

    # subgraph
    g <- igraph::induced_subgraph(g, vids = keep_nodes)

    # Rebuild graph
    el <- igraph::as_data_frame(g, what = "edges")

    # Ensure valid weights
    el$distance <- (1 - abs(el$weight)) + 1e-6
    el$distance[is.na(el$distance)] <- 1e-6
    el$distance[el$distance <= 0] <- 1e-6

    igraph::E(g)$distance <- el$distance
    igraph::E(g)$weight <- el$weight

    # rebuild nodes
    nodes_df <- data.frame(
      name = igraph::V(g)$name,
      stringsAsFactors = FALSE
    )

    nodes_df$degree <- igraph::degree(g)

    nodes_df$betweenness <- igraph::betweenness(
      g,
      weights = igraph::E(g)$distance
    )

    nodes_df$closeness <- igraph::closeness(
      g,
      weights = igraph::E(g)$distance
    )

    nodes_df$hub_score <- scale(nodes_df$degree) +
      scale(nodes_df$betweenness) +
      scale(nodes_df$closeness)

    nodes_df$hub_score <- as.numeric(nodes_df$hub_score)
    nodes_df$symbol <- nodes_df$name

    tg <- tidygraph::as_tbl_graph(g)
  }

  tg <- tg |>
    tidygraph::activate(edges) |>
    dplyr::mutate(layout_weight = abs(weight))

  # ===========================================================================
  # 8.3 Community detection
  # ===========================================================================
  if (community_method == "louvain") {
    comm <- igraph::cluster_louvain(
      g,
      weights = abs(igraph::E(g)$weight)
    )
  }

  else if (community_method == "leiden") {
    comm <- igraph::cluster_leiden(
      g,
      weights = abs(igraph::E(g)$weight)
    )
  }

  else {
    comm <- NULL
  }

  if (!is.null(comm)) {
    nodes_df$community <- igraph::membership(comm)
    igraph::vertex_attr(g, "community") <- nodes_df$community

  } else {

    nodes_df$community <- NA
  }

  # ===========================================================================
  # 9) Plot
  # ===========================================================================
  if (layout %in% c("fr","kk","stress","lgl")) {
    p <- ggraph::ggraph(tg, layout = layout, weights = .data$layout_weight)

  } else {
    p <- ggraph::ggraph(tg, layout = layout)

  }

  p <- p +
    ggraph::geom_edge_link(
      ggplot2::aes(width = abs(weight)),
      colour = "grey70",
      alpha = 0.7
    ) +
    ggraph::scale_edge_width(range = c(0.2, 2)) +
    ggraph::geom_node_point(
      ggplot2::aes(colour = .data$expr_z, size = degree)
    ) +
    ggplot2::scale_size(range = c(3,8), name = "Degree") +
    ggplot2::scale_colour_gradient2(
      low = "#3B4CC0",
      mid = "white",
      high = "#d7191c",
      midpoint = 0,
      name = "Expr (z-score)"
    ) +
    ggraph::geom_node_text(
      ggplot2::aes(label = label),
      repel = TRUE,
      size = 3.5,
      max.overlaps = Inf
    ) +
    ggplot2::theme_void() +
    ggplot2::ggtitle(paste0(pathway, " - ", group))

  if (plot && interactive()) {
    print(p)
  }

  # ===========================================================================
  # 10) RNG handling
  # ===========================================================================

  rng_state <- if (exists(".Random.seed", envir = .GlobalEnv)) .Random.seed else NULL

  seed <- if (is.null(seed)) {
    "not set"
  } else {
    as.character(seed)
  }

  # ===========================================================================
  # 11) Output
  # ===========================================================================
  params = list(
    timestamp = Sys.time(),
    threshold = threshold,
    group = group,
    node_filter = node_filter,
    top_nodes = top_nodes,
    seed = seed,
    rng_state = rng_state
  )

  obj <- list(
    params = params,
    pathway = pathway,
    genes = genes_path,
    graph = g,
    edges = edges_df,
    nodes = nodes_df
  )

  class(obj) <- "rna_network"

  # ===========================================================================
  # 12) Attach to project
  # ===========================================================================
  if (save) {

    proj <- .attach_to_project(
      proj,
      obj,
      slot = "analyses",
      subtype = "network",
      prefix = "network",
      log = list(
        pathway = pathway,
        group = group,
        threshold = threshold,
        n_genes = length(genes_path),
        n_edges = nrow(edges_df),
        seed = seed
      )
    )
  }

  # ===========================================================================
  # 13) Return
  # ===========================================================================
  .print_header("RNA Pathway Network")

  .print_block("Summary", function() {

    cat("Pathway:             ", pathway, "\n")
    cat("Organism:            ", organism, "\n")
    cat("Group used:          ", group, "\n")
    cat("Correlation cut:     ", threshold, "\n")
    cat("Layout:              ", layout, "\n")
    cat("Seed:                ", seed, "\n")
  })

  .print_block("Network properties", function() {

    cat("Nodes:               ", igraph::vcount(g), "\n")
    cat("Edges:               ", igraph::ecount(g), "\n")
    cat("Communities detected:", length(unique(nodes_df$community)), "\n")
    cat("Hub genes:           ", sum(nodes_df$hub), "\n")
    top_hubs <- nodes_df[order(-nodes_df$degree), ]
    top_hubs <- head(top_hubs$symbol, 5)

    dens <- igraph::edge_density(g)
    cat("Edge density:        ", round(dens, 3), "\n")

    cat("Top hubs:            ", paste(top_hubs, collapse = ", "), "\n")

  })

  return(invisible(proj))

}
