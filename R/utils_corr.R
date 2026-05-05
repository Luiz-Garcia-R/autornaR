# ============================
# Auxiliary correlation functions
# ============================

#' @keywords internal
.prepare_corr_data <- function(expr_mat, metadata, type,
                               sample_x = NULL,
                               sample_y = NULL) {

  if (type == "group") {

    groups <- unique(metadata$Group)

    if (length(groups) != 2)
      stop("Group correlation requires exactly two groups.")

    group_means <- sapply(groups, function(g) {
      samples <- metadata$Sample[metadata$Group == g]
      rowMeans(expr_mat[, samples, drop = FALSE])
    })

    x <- group_means[,1]
    y <- group_means[,2]

    labels <- groups

    structure_info <- list(
      level = "gene-level",
      n = length(x)
    )
  }

  else if (type == "sample") {

    if (is.null(sample_x) || is.null(sample_y))
      stop("For sample correlation, provide sample_x and sample_y.")

    if (sample_x == sample_y)
      stop("sample_x and sample_y must be different.")

    if (!all(c(sample_x, sample_y) %in% colnames(expr_mat)))
      stop("Selected samples not found in expression matrix.")

    x <- expr_mat[, sample_x]
    y <- expr_mat[, sample_y]

    labels <- c(sample_x, sample_y)

    structure_info <- list(
      level = "gene-level",
      n = length(x)
    )
  }

  list(
    x = x,
    y = y,
    labels = labels,
    structure = structure_info
  )
}

#' @keywords internal
.run_corr_diagnostics <- function(x, y, method) {

  # ---------------------------
  # Basic metrics
  # ---------------------------
  n_x <- length(x)

  has_ties <- anyDuplicated(x) > 0 || anyDuplicated(y) > 0

  normal_x <- n_x <= 5000 &&
    tryCatch(stats::shapiro.test(x)$p.value > 0.05,
             error = function(e) FALSE)

  normal_y <- length(y) <= 5000 &&
    tryCatch(stats::shapiro.test(y)$p.value > 0.05,
             error = function(e) FALSE)

  # ---------------------------
  # Outlier detection (IQR rule)
  # ---------------------------
  detect_outliers <- function(v) {
    q <- stats::quantile(v, probs = c(0.25, 0.75), na.rm = TRUE)
    iqr <- q[2] - q[1]
    which(v < (q[1] - 1.5 * iqr) | v > (q[2] + 1.5 * iqr))
  }

  outlier_fraction <- max(
    length(detect_outliers(x)) / length(x),
    length(detect_outliers(y)) / length(y)
  )

  # ---------------------------
  # Method selection
  # ---------------------------
  if (method == "auto") {

    # Extreme outliers → Kendall
    if (outlier_fraction > 0.10) {
      method_used <- "kendall"

      # Moderate outliers or ties → Spearman
    } else if (outlier_fraction > 0.03 || has_ties) {
      method_used <- "spearman"

      # Large n → Pearson acceptable
    } else if (n_x > 5000) {
      method_used <- "pearson"

      # Small/moderate n
    } else if (normal_x && normal_y) {
      method_used <- "pearson"

    } else {
      method_used <- "spearman"
    }

  } else {
    method_used <- method
  }

  # ---------------------------
  # Return structured object
  # ---------------------------
  return(list(
    method_used = method_used,
    diagnostics = list(
      normal_x = normal_x,
      normal_y = normal_y,
      ties = has_ties,
      outlier_fraction = outlier_fraction
    )
  ))
}

  # ============================
  # Compute mi
  # ============================
#' @importFrom stats complete.cases

  .compute_mi <- function(x, y, method = c("knn", "discrete"), bins = 10, k = 5) {


  method <- match.arg(method)

  # Remove NA
  df <- data.frame(x = x, y = y)
  df <- df[complete.cases(df), ]

  x <- df$x
  y <- df$y

  if (length(x) < 10) {
    warning("Too few observations to estimate mutual information reliably.")
    return(NA_real_)
  }

  if (method == "discrete") {

    # ---------------------------
    # Discretization-based MI
    # ---------------------------

    # Equal-frequency bins
    x_disc <- cut(x,
                  breaks = quantile(x, probs = seq(0, 1, length.out = bins + 1), na.rm = TRUE),
                  include.lowest = TRUE,
                  labels = FALSE)

    y_disc <- cut(y,
                  breaks = quantile(y, probs = seq(0, 1, length.out = bins + 1), na.rm = TRUE),
                  include.lowest = TRUE,
                  labels = FALSE)

    # Joint distribution
    joint <- table(x_disc, y_disc)
    joint <- joint / sum(joint)

    px <- rowSums(joint)
    py <- colSums(joint)

    mi <- 0

    for (i in seq_len(nrow(joint))) {
      for (j in seq_len(ncol(joint))) {

        if (joint[i, j] > 0) {
          mi <- mi + joint[i, j] * log(joint[i, j] / (px[i] * py[j]))
        }
      }
    }

    return(as.numeric(mi))
  }

  if (method == "knn") {

    # ---------------------------
    # KNN-based MI
    # ---------------------------
    if (!requireNamespace("FNN", quietly = TRUE)) {
      stop("Package 'FNN' is required for method = 'knn'.")
    }

    # KNN distance
    xy <- cbind(x, y)

    nn <- FNN::get.knn(xy, k = k)

    # Compute distance
    eps <- nn$nn.dist[, k]

    # Marginal count
    nx <- sapply(seq_along(x), function(i) {
      sum(abs(x - x[i]) < eps[i]) - 1
    })

    ny <- sapply(seq_along(y), function(i) {
      sum(abs(y - y[i]) < eps[i]) - 1
    })

    n <- length(x)

    # Kraskov estimator (simplified)
    mi <- digamma(k) + digamma(n) -
      mean(digamma(nx + 1) + digamma(ny + 1))

    return(as.numeric(mi))
  }
}

  # ============================
  # Classify relationship
  # ============================
  .classify_relationship <- function(r, mi_signal, non_linear_signal) {

    if (is.null(mi_signal)) return(NULL)

    r_abs <- abs(r)

    if (mi_signal < 0.05) {

      return("no clear association")

    } else if (r_abs > 0.8 && non_linear_signal < 0.1) {

      return("strong linear relationship")

    } else if (r_abs > 0.8 && non_linear_signal >= 0.1) {

      return("strong linear with non-linear component")

    } else if (r_abs > 0.4 && mi_signal > 0.2) {

      return("monotonic (possibly non-linear)")

    } else if (r_abs < 0.3 && mi_signal > 0.2) {

      return("non-linear complex relationship")

    } else {

      return("weak or ambiguous relationship")
    }
  }

  # ============================
  # Compute entropy
  # ============================
  .compute_entropy <- function(x, method = "discrete", bins = 10) {

    if (method == "discrete") {

      x_disc <- cut(x, breaks = bins, labels = FALSE)
      p <- table(x_disc) / length(x_disc)

    } else {
      stop("Entropy currently implemented only for 'discrete' method.")
    }

    p <- p[p > 0]

    -sum(p * log(p))
  }

  # ============================
  # Normalize mi
  # ============================
  .normalize_mi <- function(mi, x, y, method = "discrete") {

    hx <- .compute_entropy(x, method = method)
    hy <- .compute_entropy(y, method = method)

    if (hx == 0 || hy == 0) return(NA)

    mi / sqrt(hx * hy)
  }
