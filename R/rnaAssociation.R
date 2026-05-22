#' Association analysis between genes and/or metadata variables
#'
#' Computes pairwise associations between features derived from normalized
#' RNA-seq expression data and/or sample metadata, returning both a structured
#' result object and a heatmap-style visualization.
#'
#' The function currently supports continuous variables using correlation-based
#' association metrics and is designed to be extensible to future multi-omics
#' and mixed-data analyses.
#'
#' @param project \code{rna_project} object created by
#'   \code{rna.project()}.
#' @param x_features Character vector of features to place on the x-axis.
#'   Features may represent genes or metadata variables depending on
#'   \code{x_type}.
#' @param y_features Optional character vector of features to place on the
#'   y-axis. If \code{NULL}, uses \code{x_features}.
#' @param x_type Character. Type of x-axis features. Currently supported:
#'   \code{"gene"} and \code{"metadata"}.
#' @param y_type Character. Type of y-axis features. Currently supported:
#'   \code{"gene"} and \code{"metadata"}.
#' @param format Character. Data format expected for association analysis.
#'   Currently supports \code{"continuous"}.
#' @param group Optional character. Restricts analysis to samples belonging
#'   to the specified group in metadata.
#' @param method Correlation method used for association analysis.
#'   One of \code{"pearson"}, \code{"spearman"}, or \code{"kendall"}.
#' @param compute_mi Logical. Whether to additionally compute mutual
#'   information for each association pair (default: \code{FALSE}).
#' @param mi_method Method used for mutual information estimation.
#'   One of \code{"knn"} or \code{"discrete"}.
#' @param cluster Logical. Whether to hierarchically cluster rows and columns
#'   when the association matrix is symmetric
#'   (default: \code{TRUE}).
#' @param save Logical. Whether to store results in the active
#'   \code{rna_project} (default: \code{TRUE}).
#' @param verbose Logical. Whether to display the association heatmap and
#'   print a summary object in the console (default: \code{TRUE}).
#'
#' @return
#' Invisibly returns either:
#'
#' \itemize{
#'   \item the updated \code{rna_project} object when
#'   \code{save = TRUE};
#'   \item or an object of class \code{rna_association} when
#'   \code{save = FALSE}.
#' }
#'
#' The \code{rna_association} object contains:
#'
#' \describe{
#'   \item{matrix}{Association matrix.}
#'   \item{long}{Long-format association table containing correlation
#'   coefficients and p-values.}
#'   \item{plot}{Generated \pkg{ggplot2} heatmap object.}
#'   \item{x_features}{Features used on the x-axis.}
#'   \item{y_features}{Features used on the y-axis.}
#'   \item{x_type}{Type of x-axis features.}
#'   \item{y_type}{Type of y-axis features.}
#'   \item{params}{List containing analysis parameters and metadata.}
#' }
#'
#' @details
#' The function was designed as a generalized association engine capable of
#' integrating transcriptomic and metadata-derived variables within a unified
#' framework.
#'
#' Currently, associations are computed using correlation-based statistics
#' applied to continuous variables.
#'
#' Gene identifiers may be provided either as gene symbols or Ensembl IDs,
#' provided that annotation information is available in the project object.
#'
#' When \code{cluster = TRUE} and the matrix is symmetric, hierarchical
#' clustering is applied using:
#'
#' \itemize{
#'   \item distance: \code{1 - correlation}
#'   \item linkage: complete linkage clustering
#' }
#'
#' If \code{compute_mi = TRUE}, mutual information values are additionally
#' estimated for each feature pair using either:
#'
#' \itemize{
#'   \item discretization-based estimation;
#'   \item K-nearest neighbors estimation.
#' }
#'
#' The resulting heatmap is generated using \pkg{ggplot2}.
#'
#' Significant associations (\code{p < 0.05}) are marked with asterisks.
#'
#' @examples
#' \dontrun{
#'
#' # Gene-gene association matrix
#' rna.association(
#'   project = my_project,
#'   x_features = c("SERTAD3", "MFSD4B", "RNF139")
#' )
#'
#' # Metadata vs genes
#' rna.association(
#'   project = my_project,
#'   x_features = c("Age", "BMI"),
#'   y_features = c("IL7R", "GZMB"),
#'   x_type = "metadata",
#'   y_type = "gene"
#' )
#'
#' # Spearman correlation with MI
#' rna.association(
#'   project = my_project,
#'   x_features = c("GeneA", "GeneB"),
#'   method = "spearman",
#'   compute_mi = TRUE
#' )
#'
#' # Restrict analysis to one group
#' rna.association(
#'   project = my_project,
#'   x_features = c("GeneA", "GeneB"),
#'   group = "Control"
#' )
#' }
#'
#' @importFrom stats cor.test complete.cases hclust as.dist
#'
#' @export

rna.association <- function(project,
                            x_features,
                            y_features = NULL,
                            x_type = "gene",
                            y_type = NULL,
                            format = "continuous",
                            group = NULL,
                            method = c("spearman", "pearson", "kendall"),
                            compute_mi = FALSE,
                            mi_method = c("knn", "discrete"),
                            cluster = TRUE,
                            save = TRUE,
                            verbose = TRUE) {

  method <- match.arg(method)
  mi_method <- match.arg(mi_method)

  # ---------------------------
  # Load data
  # ---------------------------
  expr_mat <- as.matrix(.get_expr(project))
  metadata <- .get_meta(project)

  # ---------------------------
  # Validate input
  # ---------------------------
  if (is.null(y_features)) {
    y_features <- x_features
    y_type <- x_type
  }

  if (is.null(y_type)) {

    y_type <- x_type

  }

  # Group filtering
  if (!is.null(group)) {

    samples_keep <- metadata$Sample[
      metadata$Group == group
    ]

    expr_mat <- expr_mat[, samples_keep, drop = FALSE]

    metadata <- metadata[
      metadata$Sample %in% samples_keep,
      ,
      drop = FALSE
    ]
  }

  # ---------------------------
  # Resolve features
  # ---------------------------
  x_data <- .resolve_features(
    project = project,
    expr_mat = expr_mat,
    metadata = metadata,
    features = x_features,
    type = x_type,
    format = format
  )

  y_data <- .resolve_features(
    project = project,
    expr_mat = expr_mat,
    metadata = metadata,
    features = y_features,
    type = y_type,
    format = format
  )

  # ---------------------------
  # Association matrix
  # ---------------------------
  corr_mat <- matrix(
    NA,
    nrow = nrow(y_data),
    ncol = nrow(x_data)
  )

  rownames(corr_mat) <- rownames(y_data)
  colnames(corr_mat) <- rownames(x_data)

  long <- data.frame()

  for (i in seq_len(nrow(x_data))) {

    for (j in seq_len(nrow(y_data))) {

      x <- as.numeric(x_data[i, ])
      y <- as.numeric(y_data[j, ])

      ok <- complete.cases(x, y)

      x <- x[ok]
      y <- y[ok]

      if (length(x) < 3) next

      test <- suppressWarnings(
        stats::cor.test(x, y, method = method)
      )

      r <- unname(test$estimate)

      corr_mat[j, i] <- r

      long <- rbind(long, data.frame(
        x_feature = rownames(x_data)[i],
        y_feature = rownames(y_data)[j],
        r = r,
        p = test$p.value
      ))
    }
  }

  # ---------------------------
  # Optional MI
  # ---------------------------
  if (compute_mi) {

    long$mi <- NA_real_

    for (k in seq_len(nrow(long))) {

      x <- as.numeric(
        x_data[as.character(long$x_feature[k]), ]
      )

      y <- as.numeric(
        y_data[as.character(long$y_feature[k]), ]
      )

      ok <- complete.cases(x, y)

      long$mi[k] <- .compute_mi(
        x[ok],
        y[ok],
        method = mi_method
      )
    }
  }

  # ---------------------------
  # Optional clustering
  # ---------------------------
  if (cluster &&
      identical(rownames(corr_mat),
                colnames(corr_mat))) {

    ord <- hclust(
      as.dist(1 - corr_mat)
    )$order

    corr_mat <- corr_mat[ord, ord]

    levels_ord <- rownames(corr_mat)

    long$x_feature <- factor(
      long$x_feature,
      levels = levels_ord
    )

    long$y_feature <- factor(
      long$y_feature,
      levels = levels_ord
    )
  }

  # ---------------------------
  # Plot
  # ---------------------------
  g <- ggplot2::ggplot(
    long,
    ggplot2::aes(
      x_feature,
      y_feature,
      fill = r
    )
  ) +

    ggplot2::geom_tile(
      color = "grey85",
      linewidth = 0.3
    ) +

    ggplot2::geom_text(
      ggplot2::aes(
        label = ifelse(p < 0.05, "*", "")
      ),
      size = 5
    ) +

    ggplot2::scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-1, 1)
    ) +

    ggplot2::theme_minimal(base_size = 12) +

    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1
      ),
      panel.grid = ggplot2::element_blank()
    ) +

    ggplot2::labs(
      title = paste("Association matrix: ", group),
      x = NULL,
      y = NULL,
      fill = "r"
    )

  if (verbose) {
    print(g)
  }

  # ---------------------------
  # Parameters
  # ---------------------------
  params <- list(
    timestamp = Sys.time(),
    method = method,
    compute_mi = compute_mi,
    mi_method = mi_method,
    cluster = cluster,
    x_type = x_type,
    y_type = y_type,
    format = format,
    group = group
  )

  # ---------------------------
  # Output object
  # ---------------------------
  obj <- list(
    matrix = corr_mat,
    long = long,
    plot = g,
    x_features = x_features,
    y_features = y_features,
    x_type = x_type,
    y_type = y_type,
    params = params
  )

  class(obj) <- "rna_association"

  # ---------------------------
  # Optional save
  # ---------------------------
  if (save) {

    project <- .attach_to_project(
      project,
      obj,
      slot = "analyses",
      subtype = "association",
      prefix = "assoc",
      log = list(
        method = method,
        x_type = x_type,
        y_type = y_type,
        n_x = length(x_features),
        n_y = length(y_features)
      )
    )

    print(obj)

    return(invisible(project))
  }

  return(invisible(obj))
}


# ---------------------------
# Print S3
# ---------------------------
#' @method print rna_association
#' @export

print.rna_association <- function(x, ...) {

  .print_header("RNA Association")

  .print_block("Association info", function() {

    cat("X type: ", x$x_type, "\n", sep = "")
    cat("Y type: ", x$y_type, "\n", sep = "")

    cat("X features: ",
        length(x$x_features),
        "\n",
        sep = "")

    cat("Y features: ",
        length(x$y_features),
        "\n",
        sep = "")

    cat("Method: ",
        x$params$method,
        "\n",
        sep = "")

    cat("Mutual information: ",
        x$params$compute_mi,
        "\n",
        sep = "")

    cat("Clustered: ",
        x$params$cluster,
        "\n",
        sep = "")

    cat("Group: ",
        x$params$group,
        "\n",
        sep = "")
  })

  .print_block("Top associations", function() {

    top <- x$long

    top <- top[
      order(abs(top$r), decreasing = TRUE),
    ]

    print(utils::head(top, 10))
  })

  invisible(x)
}
