#' Pathway gene coexpression network from RNA-seq data
#'
#' @description
#' Builds a gene coexpression network for genes belonging to a selected
#' pathway using normalized RNA-seq expression data stored in the active
#' \code{rna_project}. Edges represent pairwise gene correlations within
#' a specified experimental group, allowing exploration of pathway-level
#' regulatory structure and hub genes.
#'
#' The function integrates pathway information from previous GSEA results,
#' constructs a correlation network, computes topological metrics, optionally
#' detects communities, and generates a publication-ready network plot.
#'
#' @param project \code{rna_project} object created by \code{rna.project()}.
#' @param pathway Character or numeric. Pathway name or index from the
#' previously computed GSEA results. If \code{NULL}, the top enriched
#' pathway is used automatically.
#' @param group Character. Experimental group used to compute correlations.
#' If \code{NULL}, the first group in the metadata is used.
#' @param threshold Numeric. Minimum absolute correlation required for an
#' edge to be included in the network (default: \code{0.6}).
#' @param cor_p Numeric. Maximum p-value allowed for correlations to be
#' retained (default: \code{0.05}).
#' @param only_hubs Logical. If \code{TRUE}, only hub genes are labeled
#' in the network plot (default: \code{FALSE}).
#' @param cor_method Character. Correlation method used to build the network:
#' \itemize{
#'   \item \code{"pearson"}: Pearson correlation
#'   \item \code{"spearman"}: Spearman rank correlation
#'   \item \code{"auto"}: automatically selects Spearman when
#'   sample size is small (\code{n < 8}), otherwise Pearson
#' }
#' @param community_method Character. Community detection algorithm applied
#' to the network:
#' \itemize{
#'   \item \code{"louvain"}: Louvain modularity optimization
#'   \item \code{"leiden"}: Leiden algorithm
#'   \item \code{"none"}: skip community detection
#' }
#' @param layout Character. Graph layout algorithm used for visualization:
#' \itemize{
#'   \item \code{"fr"}: Fruchterman-Reingold force-directed layout
#'   \item \code{"kk"}: Kamada-Kawai layout
#'   \item \code{"stress"}: stress majorization layout
#'   \item \code{"lgl"}: large graph layout
#'   \item \code{"graphopt"}: graph optimization layout
#' }
#' @param seed Optional integer used to set the random seed for reproducibility.
#' If \code{NULL}, the current random number generator state is preserved.
#' @param plot Logical. Whether to display the network plot when running
#' interactively (default: \code{TRUE}).
#' @param save Logical. Whether to store the results in the active
#' \code{rna_project} (default: \code{TRUE}).
#'
#' @details
#' The function expects normalized RNA-seq data and GSEA results already
#' stored in the project object.
#'
#' Required project components:
#' \itemize{
#'   \item \code{rna_project$data$normalized_data$expr_matrix}: gene expression
#'   matrix (genes x samples)
#'   \item \code{rna_project$data$normalized_data$metadata}: sample metadata
#'   containing group annotations
#'   \item \code{rna_project$analyses$gsea}: previously computed GSEA results
#' }
#'
#' The selected pathway genes are extracted from the GSEA results and matched
#' to the expression matrix. Pairwise correlations are computed across samples
#' within the chosen experimental group, and edges are retained according to
#' correlation magnitude and statistical significance.
#'
#' \strong{Network construction:}
#' \itemize{
#'   \item Edges represent pairwise gene correlations
#'   \item Correlations are filtered by magnitude and p-value
#'   \item Edge weights correspond to correlation strength
#'   \item Distances used for centrality metrics are defined as
#'   \code{1 - |correlation|}
#' }
#'
#' \strong{Topological metrics:}
#' \itemize{
#'   \item degree centrality
#'   \item betweenness centrality
#'   \item closeness centrality
#'   \item composite hub score
#' }
#'
#' Hub genes are defined as nodes with degree greater than
#' \code{mean(degree) + sd(degree)}.
#'
#' When enabled, community detection partitions the network into modules
#' using modularity-based algorithms.
#'
#' \strong{Visualization:}
#'
#' The resulting network plot displays:
#' \itemize{
#'   \item nodes sized by degree
#'   \item nodes colored by log2 fold-change from differential expression analysis
#'   \item edge thickness proportional to correlation magnitude
#'   \item gene labels optionally restricted to hub genes
#' }
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
#'             pathway = "Cellular Response To Starvation")
#'
#' # Restrict labels to hub genes
#' rna.network(my_project,
#'             only_hubs = TRUE)
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
                        only_hubs = FALSE,
                        cor_method = c("pearson","spearman","auto"),
                        community_method = c("louvain", "leiden", "none"),
                        layout = c("fr","kk","stress","lgl","graphopt"),
                        seed = NULL,
                        plot = TRUE,
                        save = TRUE) {

  cor_method <- match.arg(cor_method)
  community_method <- match.arg(community_method)
  layout <- match.arg(layout)

  # ---------------------------
  # 0) Basic checks
  # ---------------------------
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

  # ---------------------------
  # 1) Get active project
  # ---------------------------
  proj <- project

  expr_mat <- as.matrix(.get_expr(proj))
  meta <- .get_meta(proj)

  organism <- .get_organism(proj)
  gene_id_type <- .get_gene_id_type(proj)
  gsea_obj <- .get_gsea(proj)
  comp_obj <- .get_comp_obj(proj)

  # ---------------------------
  # 2) Validate input
  # ---------------------------
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
    warning("No comparison results found. logFC will not be available.")
  }

  # ---------------------------
  # 3) Select pathway
  # ---------------------------
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

  # ---------------------------
  # 4) Extract genes
  # ---------------------------
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

  # ---------------------------
  # 5) Select group
  # ---------------------------
  if (is.null(group)) {

    group <- unique(meta$Group)[1]
    message("[rna.network] Using group: ", group)

  }

  samples_group <- meta$Sample[meta$Group == group]

  expr_sub <- expr_mat[genes_path, samples_group, drop = FALSE]

  # ---------------------------
  # 6) Correlation
  # ---------------------------
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

  # ---------------------------
  # 7) get nodes
  # ---------------------------
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

  nodes_df$logFC <- NA_real_

  if (!is.null(comp_obj) && !is.null(comp_obj$res)) {

    if ("log2FoldChange" %in% colnames(comp_obj$res)) {

      logfc <- comp_obj$res$log2FoldChange
      names(logfc) <- rownames(comp_obj$res)

      nodes_df$logFC <- logfc[match(nodes_df$name, names(logfc))]

    }
  }

  # ---------------------------
  # 8) Build network
  # ---------------------------
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

  # ---------------------------
  # 8.1 Network metrics
  # ---------------------------

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

  if (only_hubs) {

    nodes_df$label[!nodes_df$hub] <- ""
    igraph::vertex_attr(g, "label") <- nodes_df$label
    tg <- tidygraph::as_tbl_graph(g)

  }

  tg <- tg |>
    tidygraph::activate(edges) |>
    dplyr::mutate(layout_weight = abs(weight))

  # ---------------------------
  # 8.2 Community detection
  # ---------------------------
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

  # ---------------------------
  # 9) Plot
  # ---------------------------
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
      ggplot2::aes(colour = .data$logFC, size = degree)
    ) +
    ggplot2::scale_size(range = c(3,8), name = "Degree") +
    ggplot2::scale_colour_gradient2(
      low = "#3B4CC0",
      mid = "white",
      high = "#d7191c",
      midpoint = 0,
      name = "log2FC"
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

  # ---------------------------
  # 10) RNG handling
  # ---------------------------

  rng_state <- if (exists(".Random.seed", envir = .GlobalEnv)) .Random.seed else NULL

  seed <- if (is.null(seed)) {
    "not set"
  } else {
    as.character(seed)
  }

  # ---------------------------
  # 11) Output
  # ---------------------------
  params = list(
    timestamp = Sys.time(),
    threshold = threshold,
    group = group,
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

  # ---------------------------
  # 12) Attach to project
  # ---------------------------
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

  # ---------------------------
  # 13) Return
  # ---------------------------
  .print_header("RNA Pathway Network")

  .print_block("Summary", function() {

    cat("Pathway:             ", pathway, "\n")
    cat("Organism:            ", organism, "\n")
    cat("Group used:          ", group, "\n")
    cat("Genes in pathway:    ", length(genes_path), "\n")
    cat("Edges detected:      ", nrow(edges_df), "\n")
    cat("Communities detected:", length(unique(nodes_df$community)), "\n")
    cat("Correlation cut:     ", threshold, "\n")
    cat("Layout:              ", layout, "\n")
    cat("Seed:                ", seed, "\n")
  })

  .print_block("Network properties", function() {

    cat("Nodes:               ", igraph::vcount(g), "\n")
    cat("Edges:               ", igraph::ecount(g), "\n")
    cat("Hub genes:           ", sum(nodes_df$hub), "\n")
    top_hubs <- nodes_df[order(-nodes_df$degree), ]
    top_hubs <- head(top_hubs$symbol, 5)

    cat("Top hubs:            ", paste(top_hubs, collapse = ", "), "\n")

    dens <- igraph::edge_density(g)
    cat("Edge density:        ", round(dens, 3), "\n")

  })

  return(invisible(proj))

}
