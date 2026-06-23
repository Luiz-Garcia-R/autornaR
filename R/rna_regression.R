#' Regression analysis between genes and/or metadata variables
#'
#' Fits linear regression models using normalized RNA-seq expression data
#' and/or sample metadata variables as predictors and outcomes.
#'
#' The function supports gene expression features and metadata-derived
#' variables through automatic feature resolution and generates both a
#' structured regression object and an optional visualization when a
#' single predictor is supplied.
#'
#' @param project \code{rna_project} object created by
#'   \code{rna.project()}.
#' @param outcome Character. Outcome variable to be modeled.
#'   May represent either a gene or metadata variable depending on
#'   \code{outcome_type}.
#' @param predictors Character vector of predictor variables.
#'   Predictors may represent genes and/or metadata variables depending on
#'   \code{predictor_type}.
#' @param covariates Optional character vector of metadata variables to be
#'   included as adjustment covariates in the model.
#' @param outcome_type Character. Type of outcome variable.
#'   One of \code{"auto"}, \code{"gene"}, or \code{"metadata"}.
#' @param predictor_type Character. Type of predictor variables.
#'   One of \code{"auto"}, \code{"gene"}, or \code{"metadata"}.
#' @param group Optional character. Restricts analysis to samples belonging
#'   to the specified group in metadata.
#' @param save Logical. Whether to store results in the active
#'   \code{rna_project} (default: \code{TRUE}).
#' @param verbose Logical. Whether to display generated plots and print
#'   a summary object in the console (default: \code{TRUE}).
#'
#' @return
#' Invisibly returns either:
#'
#' \itemize{
#'   \item the updated \code{rna_project} object when
#'   \code{save = TRUE};
#'   \item or an object of class \code{rna_regression} when
#'   \code{save = FALSE}.
#' }
#'
#' The \code{rna_regression} object contains:
#'
#' \describe{
#'   \item{formula}{Model formula used for fitting.}
#'   \item{fit}{Fitted \code{lm} object.}
#'   \item{coefficients}{Regression coefficient table.}
#'   \item{plot}{Generated \pkg{ggplot2} visualization (when available).}
#'   \item{outcome}{Outcome variable.}
#'   \item{predictors}{Predictor variables.}
#'   \item{covariates}{Adjustment covariates.}
#'   \item{params}{List containing analysis parameters and metadata.}
#' }
#'
#'
#' @details
#' The function was designed as a generalized regression engine for
#' transcriptomic and metadata-derived variables.
#'
#' Features may be specified using either Ensembl identifiers or gene
#' symbols when annotation information is available in the project object.
#'
#' Currently, linear regression models are fitted using
#' \code{\link[stats]{lm}}.
#'
#' Unlike correlation analyses, regression models explicitly distinguish
#' between an outcome variable and one or more predictor variables.
#' Consequently, the model estimates how changes in predictor values are
#' associated with changes in the outcome while preserving directionality
#' (\code{outcome ~ predictors}).
#'
#' For single-predictor models, the regression coefficient
#' (\eqn{\beta}) represents the expected change in the outcome for a
#' one-unit increase in the predictor. The coefficient of determination
#' (\eqn{R^2}) quantifies the proportion of outcome variability explained
#' by the model. When only one predictor is included, \eqn{R^2}
#' corresponds to the squared Pearson correlation coefficient between the
#' predictor and outcome.
#'
#' Regression analysis provides a flexible framework for studying
#' gene-gene and gene-phenotype relationships because additional
#' covariates can be incorporated directly into the model. This allows
#' users to evaluate associations while adjusting for potential
#' confounding variables such as age, sex, batch effects, or other sample
#' characteristics.
#'
#' Missing values are removed using complete-case filtering prior to
#' model fitting.
#'
#' When only a single predictor is supplied, a regression plot is
#' automatically generated showing:
#'
#' \itemize{
#'   \item observed data points;
#'   \item fitted regression line;
#'   \item 95\% confidence interval;
#'   \item model summary statistics (\eqn{\beta}, p-value and R^2).
#' }
#'
#'
#' \itemize{
#'   \item R^2;
#'   \item adjusted R^2;
#'   \item residual normality p-value (Shapiro-Wilk test);
#'   \item residual standard error (\code{sigma}).
#' }
#'
#' The generated visualization uses gene symbols whenever available,
#' while model fitting internally preserves the original feature
#' identifiers.
#'
#' @examples
#' \dontrun{
#'
#' # Simple gene-gene regression
#' rna_regression(
#'   project = my_project,
#'   outcome = "SERTAD3",
#'   predictors = "MFSD4B"
#' )
#'
#' # Multiple predictors
#' rna_regression(
#'   project = my_project,
#'   outcome = "SERTAD3",
#'   predictors = c(
#'     "MFSD4B",
#'     "RNF139",
#'     "NFKBIA"
#'   )
#' )
#'
#' # Adjusted model
#' rna_regression(
#'   project = my_project,
#'   outcome = "SERTAD3",
#'   predictors = "MFSD4B",
#'   covariates = c(
#'     "Age",
#'     "Sex"
#'   )
#' )
#'
#' # Metadata outcome
#' rna_regression(
#'   project = my_project,
#'   outcome = "BMI",
#'   predictors = c(
#'     "IL7R",
#'     "GZMB"
#'   ),
#'   outcome_type = "metadata"
#' )
#'
#' # Restrict analysis to one group
#' rna_regression(
#'   project = my_project,
#'   outcome = "SERTAD3",
#'   predictors = "MFSD4B",
#'   group = "Control"
#' )
#' }
#'
#' @importFrom stats lm as.formula residuals shapiro.test
#'
#' @export

rna_regression <- function(
    project,
    outcome,
    predictors,
    covariates = NULL,
    outcome_type = c("auto", "gene", "metadata"),
    predictor_type = c("auto", "gene", "metadata"),
    group = NULL,
    save = TRUE,
    verbose = TRUE
) {

  outcome_type <- match.arg(outcome_type)
  predictor_type <- match.arg(predictor_type)

  # ---------------------------
  # Get active project
  # ---------------------------
  expr_mat <- as.matrix(.get_expr(project))
  metadata <- .get_meta(project)

  gene_map <- .get_gene_annotation(project)
  gene_map <- .align_gene_annotation(
    gene_map,
    expr_mat
  )

  # ---------------------------
  # Resolve gene ID
  # ---------------------------
  get_label <- function(feature) {

    idx <- which(
      gene_map$gene_id == feature |
        gene_map$symbol == feature
    )

    if (length(idx) == 0)
      return(feature)

    symbol <- gene_map$symbol[idx[1]]

    if (is.na(symbol) || symbol == "")
      return(feature)

    symbol
  }

  # ---------------------------
  # Group filtering
  # ---------------------------

  if (!is.null(group)) {

    samples_keep <- metadata$Sample[
      metadata$Group == group
    ]

    expr_mat <- expr_mat[
      ,
      samples_keep,
      drop = FALSE
    ]

    metadata <- metadata[
      metadata$Sample %in% samples_keep,
      ,
      drop = FALSE
    ]
  }

  # ---------------------------
  # Outcome
  # ---------------------------

  outcome_data <- .resolve_features(
    project = project,
    expr_mat = expr_mat,
    metadata = metadata,
    features = outcome,
    type = outcome_type
  )

  if (nrow(outcome_data) != 1) {

    stop(
      "'outcome' must contain exactly one variable."
    )
  }

  # ---------------------------
  # Predictors
  # ---------------------------
  predictor_data <- .resolve_features(
    project = project,
    expr_mat = expr_mat,
    metadata = metadata,
    features = predictors,
    type = predictor_type
  )

  # ---------------------------
  # Covariates
  # ---------------------------
  cov_data <- NULL

  if (!is.null(covariates)) {

    cov_data <- .resolve_features(
      project = project,
      expr_mat = expr_mat,
      metadata = metadata,
      features = covariates,
      type = "metadata"
    )
  }

  # ---------------------------
  # Build model dataframe
  # ---------------------------
  model_df <- data.frame(

    outcome = as.numeric(
      outcome_data[1, ]
    ),

    t(predictor_data),
    check.names = FALSE
  )

  colnames(model_df) <- c(
    outcome,
    predictors
  )

  if (!is.null(cov_data)) {

    cov_df <- as.data.frame(
      t(cov_data),
      check.names = FALSE
    )

    model_df <- cbind(
      model_df,
      cov_df
    )
  }

  model_df <- model_df[
    complete.cases(model_df),
    ,
    drop = FALSE
  ]

  # ---------------------------
  # Formula
  # ---------------------------

  rhs <- c(
    predictors,
    covariates
  )

  formula_obj <- stats::as.formula(
    paste(
      outcome,
      "~",
      paste(rhs, collapse = " + ")
    )
  )

  # ---------------------------
  # Fit model
  # ---------------------------

  fit <- stats::lm(
    formula_obj,
    data = model_df
  )

  fit_sum <- summary(fit)

  coef_table <- as.data.frame(
    fit_sum$coefficients
  )

  coef_table$variable <- rownames(
    coef_table
  )

  rownames(coef_table) <- NULL

  names(coef_table) <- c(
    "estimate",
    "std_error",
    "t-statistic",
    "p-value",
    "variable"
  )

  coef_table <- coef_table[
    ,
    c(
      "variable",
      "estimate",
      "std_error",
      "t-statistic",
      "p-value"
    )
  ]

  # ---------------------------
  # Diagnostics
  # ---------------------------
  residuals_fit <- residuals(fit)

  diagnostics_obj <- list(

    r_squared =
      fit_sum$r.squared,

    adj_r_squared =
      fit_sum$adj.r.squared,

    residual_normality_p =
      tryCatch(
        shapiro.test(residuals_fit)$p.value,
        error = function(e) NA
      ),

    sigma =
      fit_sum$sigma
  )

  # ---------------------------
  # Plot
  # ---------------------------
  g <- NULL

  if (length(predictors) == 1) {

    pred <- predictors[1]

    pred_label <- get_label(pred)
    outcome_label <- get_label(outcome)

    # Prepare title text
    title_text <- paste0(
      outcome_label,
      " ~ ",
      pred_label,
      if (!is.null(group))
        paste0(" (", group, ")")
      else
        ""
    )

    coef_main <- coef_table[
      coef_table$variable == pred,
      ,
      drop = FALSE
    ]

    if (nrow(coef_main) != 1) {

      stop(
        "Could not identify predictor coefficient: ",
        pred
      )
    }

    beta <- round(
      coef_main$estimate,
      3
    )

    pval <- signif(
      coef_main$p,
      3
    )

    r2 <- round(
      diagnostics_obj$r_squared,
      3
    )

    subtitle_text <- paste0(
      "\u03B2 = ", beta,
      " | R^2 = ", r2,
      " | p = ", pval
    )

    line_col <- "#2F5D8A"
    ci_fill  <- "#2F5D8A"
    point_col <- "grey40"

    g <- ggplot2::ggplot(
      model_df,
      ggplot2::aes(
        x = .data[[pred]],
        y = .data[[outcome]]
      )
    ) +

      ggplot2::geom_point(
        alpha = 0.40,
        size = 2.3,
        color = point_col
      ) +

      ggplot2::geom_smooth(
        method = "lm",
        se = TRUE,
        level = 0.95,
        color = line_col,
        fill = ci_fill,
        alpha = 0.25,
        linewidth = 1.1
      ) +

      ggplot2::theme_minimal(
        base_size = 12
      ) +

      ggplot2::labs(
        title = title_text,
        subtitle = subtitle_text,
        x = pred_label,
        y = outcome_label
      )
  }

  n_samples = nrow(model_df)

  if (verbose && !is.null(g)) {
    print(g)
  }

  # ---------------------------
  # Parameters
  # ---------------------------
  params <- list(
    timestamp = Sys.time(),
    group = group,
    outcome_type = outcome_type,
    predictor_type = predictor_type
  )

  # ---------------------------
  # Output
  # ---------------------------
  obj <- list(

    params = params,

    outcome = outcome,
    predictors = predictors,
    covariates = covariates,

    n_samples = n_samples,

    formula = formula_obj,

    fit = fit,

    coefficients = coef_table,

    model_metrics = list(
      r_squared = fit_sum$r.squared,
      adj_r_squared = fit_sum$adj.r.squared,
      sigma = fit_sum$sigma,
      residual_normality_p =
        diagnostics_obj$residual_normality_p,
      AIC = stats::AIC(fit),
      BIC = stats::BIC(fit)
    ),

    residuals = residuals(fit),
    fitted_values = stats::fitted(fit),

    model_frame = model_df,

    plot = g
  )

  class(obj) <- "rna_regression"

  # ---------------------------
  # Save
  # ---------------------------
  if (save) {

    project <- .attach_to_project(
      project,
      obj,
      slot = "analyses",
      subtype = "regression",
      prefix = "reg",
      log = list(
        outcome = outcome,
        predictors =
          length(predictors),
        covariates =
          length(covariates)
      )
    )

    if (verbose) {
      print(obj)
    }

    return(invisible(project))}

  return(invisible(obj))

  }

# ---------------------------
# Print S3
# ---------------------------
#' @method print rna_regression
#' @export

print.rna_regression <- function(x, ...) {

  .print_header("RNA Regression")

  .print_block("Model info", function() {

    cat(
      "Outcome: ",
      x$outcome,
      "\n",
      sep = ""
    )

    cat(
      "Predictors: ",
      length(x$predictors),
      "\n",
      sep = ""
    )

    cat(
      "Covariates: ",
      length(x$covariates),
      "\n",
      sep = ""
    )

    cat(
      "Group: ",
      x$params$group,
      "\n",
      sep = ""
    )

    cat(
      "n: ",
      x$n_samples,
      "\n",
      sep = ""
    )

  })

  .print_block("Formula", function() {

    print(x$formula)
  })

  .print_block("Diagnostics", function() {

    cat(
      "R^2: ",
      round(x$model_metrics$r_squared, 3),
      "\n",
      sep = ""
    )

    cat(
      "Adjusted R^2: ",
      round(x$model_metrics$adj_r_squared, 3),
      "\n",
      sep = ""
    )

    cat(
      "Sigma: ",
      round(
        x$model_metrics$sigma,
        3
      ),
      "\n",
      sep = ""
    )
  })

  .print_block("Top coefficients", function() {

    coef_tab <- x$coefficients

    coef_tab <- coef_tab[
      coef_tab$variable != "(Intercept)",
      ,
      drop = FALSE
    ]

    coef_tab <- coef_tab[
      order(abs(coef_tab$t),
            decreasing = TRUE),
      ,
      drop = FALSE
    ]

    print(
      utils::head(
        coef_tab,
        10
      )
    )
  })

  invisible(x)
}
