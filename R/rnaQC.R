#' RNA-seq Quality Control (QC)
#'
#' @description
#' Performs comprehensive quality control (QC) on normalized RNA-seq data
#' stored in the active \code{rna_project}. This function provides a structured
#' overview of sample quality, expression distributions, and overall dataset
#' consistency, helping identify technical artifacts and potential outliers
#' before downstream analyses.
#'
#' @param project \code{rna_project} object created by \code{rna.project()}.
#' @param group_col Character. Column in metadata indicating biological groups
#'   (default: \code{"Group"}).
#' @param display_corr_values Logical. If \code{TRUE}, overlays correlation coefficients
#'   on the heatmap (default: \code{FALSE}).
#' @param highlight_outliers Logical. If \code{TRUE}, detects and highlights outlier samples
#'   based on library size (default: \code{TRUE}).
#' @param cluster_rows Logical. Whether to cluster rows in the sample correlation
#'   heatmap (default: \code{FALSE}).
#' @param cluster_cols Logical. Whether to cluster columns in the sample correlation
#'   heatmap (default: \code{FALSE}).
#' @param plot Logical. If \code{TRUE}, generates diagnostic plots (default: \code{TRUE}).
#' @param save Logical. Whether to store results in the active \code{rna_project}
#'   (default: \code{TRUE}).
#'
#' @return
#' An object of class \code{"rna_QC"} containing:
#' \describe{
#'   \item{data_ready}{Expression matrix used for QC.}
#'   \item{corr_matrix_samples}{Sample-wise correlation matrix.}
#'   \item{plot_density}{Density plot of expression values (if \code{plot = TRUE}).}
#'   \item{plot_expr_boxplot}{Boxplot of expression distributions per sample.}
#'   \item{plot_lib_boxplot}{Barplot of library sizes per sample.}
#'   \item{summary}{List with basic QC metrics such as library size, genes detected,
#'   and identified outliers.}
#'   \item{qc_metrics}{Detailed QC metrics including library variation,
#'   correlation structure, and classification of dataset quality.}
#' }
#'
#' @details
#' This function evaluates RNA-seq data quality through multiple complementary perspectives:
#'
#' \strong{1. Library size assessment}
#' - Computes total counts per sample
#' - Identifies outliers using IQR-based thresholds
#' - Calculates coefficient of variation (CV) across samples
#'
#' \strong{2. Expression distribution}
#' - Boxplots reveal global expression differences between samples
#' - Density plots highlight distribution shifts across groups
#'
#' \strong{3. Sample correlation structure}
#' - Computes pairwise Pearson correlation between samples
#' - Identifies low-correlation samples
#' - Compares intra-group vs inter-group similarity
#'
#' \strong{4. Biological signal evaluation}
#' - Estimates whether group structure is preserved in the data
#' - Stronger intra-group than inter-group correlation suggests biological signal
#'
#' \strong{5. Automated QC classification}
#' - Library quality (based on CV)
#' - Correlation consistency
#' - Overall dataset assessment with diagnostic flags
#'
#' These metrics are summarized into an interpretable QC report, helping guide
#' decisions such as removing samples, revisiting normalization, or proceeding
#' with downstream analyses.
#'
#' @section Interpretation guide:
#' - \strong{High-quality data}: low library CV, high sample correlation, clear group structure
#' - \strong{Potential issues}:
#'   - High library size variation -> technical variability
#'   - Low correlation samples -> possible outliers or batch effects
#'   - Weak group separation -> low biological signal or noisy data
#'
#' QC should always be interpreted in context: biological heterogeneity,
#' experimental design, and normalization method can all influence these metrics.
#'
#' @section Pipeline context:
#' This function is typically used after \code{rna.normalize()}.
#'
#' @examples
#' \dontrun{
#' # Basic QC
#' rna.qc(project = my_project)
#'
#' # QC with custom grouping
#' rna.qc(my_project,
#'        group_col = "Condition")
#'
#' # Display correlation values and highlight structure
#' rna.qc(my_project,
#'        display_corr_values = TRUE)
#'
#' # Disable plotting (metrics only)
#' rna.qc(my_project,
#'        plot = FALSE)
#'
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_bar geom_text geom_boxplot geom_density theme_minimal labs element_text theme
#' @importFrom pheatmap pheatmap
#' @importFrom scales squish
#'
#' @export

rna.qc <- function(project,
                   group_col = "Group",
                   display_corr_values = FALSE,
                   highlight_outliers = TRUE,
                   cluster_rows = FALSE,
                   cluster_cols = FALSE,
                   plot = TRUE,
                   save = TRUE) {

  # ---------------------------
  # 0) Basic checks
  # ---------------------------
  pkgs <- c("ggplot2", "pheatmap", "data.table")

  .check_dependencies(pkgs)

  # ---------------------------
  # 1) Get active project
  # ---------------------------
  proj <- project

  expr_mat <- .get_expr(proj)
  metadata <- .get_meta(proj)

  # ---------------------------
  # 2) Validate input
  # ---------------------------
  if (is.null(expr_mat)) {
    stop("No expression matrix available. Run 'rna.normalize()' first.")
  }

  if (is.null(metadata)) {
    stop("No metadata available. Load a metadata file using 'rna.import()'")
  }

  if (!"Sample" %in% colnames(metadata)) {
    stop("metadata must contain a 'Sample' column.")
  }

  if (!group_col %in% colnames(metadata)) {
    stop(
      "Column '", group_col, "' not found in metadata.\n",
      "Available columns: ", paste(colnames(metadata), collapse = ", ")
    )
  }

  # ---------------------------
  # 3) Align samples
  # ---------------------------
  common_samples <- intersect(colnames(expr_mat), metadata$Sample)
  if (length(common_samples) == 0) stop("No matching samples between metadata and expression matrix.")
  metadata_sub <- metadata[match(common_samples, metadata$Sample), , drop = FALSE]
  expr_mat <- expr_mat[, common_samples, drop = FALSE]

  # ---------------------------
  # 4) Library sizes and genes detected
  # ---------------------------
  lib_sizes <- colSums(expr_mat)
  genes_detected <- colSums(expr_mat > 0)

  # ---------------------------
  # 5) Detect library size outliers
  # ---------------------------
  if (highlight_outliers) {
    q <- quantile(lib_sizes, probs = c(0.25, 0.75))
    iqr <- diff(q)
    outliers_lib <- names(lib_sizes)[
      lib_sizes < (q[1] - 1.5 * iqr) |
        lib_sizes > (q[2] + 1.5 * iqr)
    ]
  } else {
    outliers_lib <- integer(0)
  }

  # ---------------------------
  # 6) Percent of genes detected per group
  # ---------------------------
  percent_genes_detected <- sapply(unique(metadata_sub[[group_col]]), function(grp) {
    samples_in_group <- metadata_sub$Sample[metadata_sub[[group_col]] == grp]
    mean(colSums(expr_mat[, samples_in_group, drop = FALSE] > 0) / nrow(expr_mat)) * 100
  })

  # ---------------------------
  # 7) Prepare long format table
  # ---------------------------
  dt <- data.table::as.data.table(expr_mat, keep.rownames = "Gene")

  df_long <- data.table::melt(
    dt,
    id.vars = "Gene",
    variable.name = "Sample",
    value.name = "Expression"
  )

  # Merge metadata
  df_long <- merge(
    df_long,
    metadata_sub,
    by = "Sample",
    all.x = TRUE,
    sort = FALSE
  )

  # Preserve sample order
  df_long$Sample <- factor(df_long$Sample, levels = colnames(expr_mat))

  # ---------------------------
  # 8) Library size barplot
  # ---------------------------
  if (plot) {
  df_box <- data.frame(Sample = colnames(expr_mat), Library_size = lib_sizes,
                       Genes_detected = genes_detected, Group = metadata_sub[[group_col]])
  df_box$Sample <- factor(df_box$Sample, levels = colnames(expr_mat))
  df_box$Outlier <- df_box$Sample %in% outliers_lib
  p_boxplot <- ggplot2::ggplot(df_box, ggplot2::aes(x = Sample, y = Library_size, fill = Group)) +
    ggplot2::geom_bar(stat = "identity", alpha = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = ifelse(Outlier, "*", "")), vjust = -0.5, color = "red", size = 5) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Library size per sample", x = "", y = "Library size") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   legend.position = "bottom")
  print(p_boxplot)
  }

  # ---------------------------
  # 9) Expression boxplot
  # ---------------------------
  if (plot) {
    p_expr_box <- ggplot2::ggplot(df_long, ggplot2::aes(x = Sample, y = Expression, fill = Sample)) +
    ggplot2::geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Expression distribution per sample", x = "", y = "Expression") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   legend.position = "none")
  print(p_expr_box)
  }

  # ---------------------------
  # 10) Density plot per group
  # ---------------------------
  if (plot) {
  p_density <- ggplot2::ggplot(df_long, ggplot2::aes(x = Expression,
                                                     color = .data[[group_col]],
                                                     fill = .data[[group_col]])) +
    ggplot2::geom_density(alpha = 0.3, na.rm = TRUE) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Expression density per group", x = "Expression", y = "Density",
                  color = group_col, fill = group_col) +
    ggplot2::theme(legend.position = "bottom")
  print(p_density)
  }

  # ---------------------------
  # 11) Correlation metrics
  # ---------------------------
  corr_mat <- cor(expr_mat, use = "pairwise.complete.obs", method = "pearson")

  mean_corr <- colMeans(corr_mat)
  global_mean_corr <- mean(mean_corr)

  low_corr_threshold <- quantile(mean_corr, 0.1)
  low_corr_samples <- names(mean_corr)[mean_corr < low_corr_threshold]

  # ---------------------------
  # 12) Library size metrics
  # ---------------------------
  lib_cv <- sd(lib_sizes) / mean(lib_sizes) * 100
  lib_summary <- summary(lib_sizes)

  # ---------------------------
  # 13) Intra vs inter-group correlation
  # ---------------------------
  groups <- metadata_sub[[group_col]]
  if (length(unique(groups)) > 1) {
    same_group <- outer(groups, groups, "==")
    intra_corr <- mean(corr_mat[same_group])
    inter_corr <- mean(corr_mat[!same_group])
  } else {
    intra_corr <- NA
    inter_corr <- NA
  }

  # ---------------------------
  # 14) QC classification
  # ---------------------------
  qc_flags <- list()

  # Library size quality
  if (lib_cv < 5) {
    lib_quality <- "Excellent"
  } else if (lib_cv < 10) {
    lib_quality <- "Acceptable"
    qc_flags <- c(qc_flags, "Moderate library size variation")
  } else {
    lib_quality <- "Warning"
    qc_flags <- c(qc_flags, "High library size variation")
  }

  # Correlation quality
  if (global_mean_corr > 0.95) {
    corr_quality <- "High consistency"
  } else if (global_mean_corr > 0.90) {
    corr_quality <- "Moderate consistency"
    qc_flags <- c(qc_flags, "Moderate sample correlation")
  } else {
    corr_quality <- "Low consistency"
    qc_flags <- c(qc_flags, "Low sample correlation")
  }

  # Biological signal
  if (!is.na(intra_corr) && !is.na(inter_corr)) {
    if ((intra_corr - inter_corr) > 0.02) {
      signal_quality <- "Clear group structure"
    } else {
      signal_quality <- "Weak group separation"
    }
  } else {
    signal_quality <- NA
  }

  # Final verdict
  if (lib_quality == "Excellent" &&
      corr_quality == "High consistency") {
    qc_verdict <- "Overall QC assessment: Excellent"
  } else if (lib_quality == "Warning" ||
             corr_quality == "Low consistency") {
    qc_verdict <- "Overall QC assessment: Potential issues detected"
  } else {
    qc_verdict <- "Overall QC assessment: Acceptable"
  }

  # ---------------------------
  # 15) Sample correlation heatmap
  # ---------------------------
  Var1 <- Var2 <- Var1_num <- Var2_num <- NULL

  if (plot) {

    corr_long <- data.table::as.data.table(as.table(corr_mat))
    names(corr_long) <- c("Var1", "Var2", "Correlation")

    corr_long[, Var1_num := match(as.character(Var1), colnames(corr_mat))]
    corr_long[, Var2_num := match(as.character(Var2), colnames(corr_mat))]

    corr_long <- corr_long[!is.na(Var1_num) & !is.na(Var2_num)]
    corr_long <- corr_long[Var1_num >= Var2_num]

    corr_long <- corr_long[order(Var1_num, Var2_num)]

    # Factor adjust
    corr_long[, Var1 := factor(Var1, levels = rev(colnames(corr_mat)))]
    corr_long[, Var2 := factor(Var2, levels = colnames(corr_mat))]

    p_corr <- ggplot2::ggplot(
      corr_long,
      ggplot2::aes(x = .data$Var2, y = .data$Var1, fill = .data$Correlation)
    ) +
      ggplot2::geom_tile(color = "white") +
      ggplot2::scale_fill_gradient2(
        low = "#3B4CC0",
        mid = "white",
        high = "#d7191c",
        midpoint = 0.9,
        limits = c(0.8, 1),
        oob = scales::squish
      ) +
      ggplot2::coord_fixed() +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = "Sample correlation heatmap",
        x = "",
        y = ""
      ) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        axis.text = ggplot2::element_text(size = 10),
        plot.title = ggplot2::element_text(size = 12, hjust = 0.5),
        legend.position = "bottom"
      )

    if (display_corr_values) {
      p_corr <- p_corr +
        ggplot2::geom_text(
          ggplot2::aes(label = round(.data$Correlation, 2)),
          size = 3
        )
    }

    print(p_corr)
  }

  # ---------------------------
  # 16) Output
  # ---------------------------
  obj <- list(
    summary = list(
      timestamp = Sys.time(),
      n_genes = nrow(expr_mat),
      n_samples = ncol(expr_mat),
      percent_genes_detected = percent_genes_detected,
      library_size = lib_sizes,
      genes_detected = genes_detected,
      outliers_lib = outliers_lib
    ),
    qc_metrics = list(
      library = list(
        min = unname(lib_summary["Min."]),
        median = unname(lib_summary["Median"]),
        max = unname(lib_summary["Max."]),
        cv_percent = lib_cv
      ),
      correlation = list(
        global_mean = global_mean_corr,
        low_correlation_samples = low_corr_samples,
        intra_group_mean = intra_corr,
        inter_group_mean = inter_corr,
        plot_corr_heatmap = if (plot) p_corr else NULL
      ),
      classification = list(
        library = lib_quality,
        correlation = corr_quality,
        biological_signal = signal_quality,
        flags = qc_flags,
        verdict = qc_verdict
      )
    ),
      data_ready = expr_mat,
      corr_matrix_samples = corr_mat,
    plot = list(
      plot_density = if (plot) p_density else NULL,
      plot_expr_boxplot = if (plot) p_expr_box else NULL,
      plot_lib_boxplot = if (plot) p_boxplot else NULL
    )
  )

  class(obj) <- "rna_QC"

  # ---------------------------
  # 17) Attach to project
  # ---------------------------
  if (save) {

    proj <- .attach_to_project(
      proj,
      obj,
      slot = "analyses",
      subtype = "QC",
      prefix = "QC",
      log = list(
        n = ncol(expr_mat),
        group_col = group_col,
        outliers_detected = length(outliers_lib),
        low_corr_samples = length(low_corr_samples),
        qc_status = qc_verdict
      )
    )
  }

  # ---------------------------
  # 18) Return
  # ---------------------------
  .print_header("RNA quality control")

  .print_block("QC Summary", function() {
    cat("Samples: ", obj$summary$n_samples, "\n", sep = "")
    cat("Genes: ", obj$summary$n_genes, "\n\n", sep = "")

    cat("Library size:\n")
    cat("  Min: ", round(obj$qc_metrics$library$min, 2), "\n", sep = "")
    cat("  Median: ", round(obj$qc_metrics$library$median, 2), "\n", sep = "")
    cat("  Max: ", round(obj$qc_metrics$library$max, 2), "\n", sep = "")
    cat("  CV (%): ", round(obj$qc_metrics$library$cv_percent, 2), "\n\n", sep = "")

    cat("Library size outliers: ", length(obj$summary$outliers_lib), "\n\n", sep = "")

    cat("Mean sample correlation: ",
        round(obj$qc_metrics$correlation$global_mean, 3), "\n", sep = "")

    if (length(obj$qc_metrics$correlation$low_correlation_samples) > 0) {
      cat("Low-correlation samples: ",
          paste(obj$qc_metrics$correlation$low_correlation_samples,
                collapse = ", "),
          "\n", sep = "")
    }

    if (!is.na(obj$qc_metrics$correlation$intra_group_mean)) {
      cat("\nIntra-group correlation: ",
          round(obj$qc_metrics$correlation$intra_group_mean, 3), "\n", sep = "")
      cat("Inter-group correlation: ",
          round(obj$qc_metrics$correlation$inter_group_mean, 3), "\n", sep = "")
    }
  })

    .print_block("QC Assessment", function() {
    cat("Library size quality: ",
        obj$qc_metrics$classification$library,
        "\n", sep = "")

    cat("Sample correlation: ",
        obj$qc_metrics$classification$correlation,
        "\n", sep = "")

    if (!is.na(obj$qc_metrics$classification$biological_signal)) {
      cat("Biological signal: ",
          obj$qc_metrics$classification$biological_signal,
          "\n", sep = "")
    }

    cat(obj$qc_metrics$classification$verdict,
        "\n", sep = "")
  })

    return(invisible(proj))
}
