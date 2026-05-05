# ============================================
# Print Methods for autornaR objects
# ============================================

# --------------------------------------------
#' Print Method for rnaQC Objects
#'
#' Custom print method for objects of class \code{rnaQC}.
#'
#' @param x An object of class \code{rnaQC}.
#' @param ... Additional arguments (ignored).
#'
#' @method print rnaQC
#' @export
#' @rdname print.rnaQC

print.rnaQC <- function(x, ...) {
  message("==============================")
  message("Object of class 'rnaQC'")
  message("==============================")
  message("Number of genes: ", x$summary$n_genes)
  message("Number of samples: ", x$summary$n_samples)
  message("Outlier samples: ", if(length(x$summary$outliers_lib) > 0) paste(x$summary$outliers_lib, collapse=", ") else "None")
  message("==============================")
  message("You can explore this object with functions like:")
  message(" - rna.rank()")
  message(" - rna.ttest()")
  message(" - rna.roc()")
  message(" - rna.volcano()")
  invisible(x)
}

# --------------------------------------------
#' Print Method for rnaROC Objects
#'
#' Custom print method for objects returned by \code{rna.roc()}.
#'
#' @param x An object returned by \code{rna.roc()}.
#' @param ... Additional arguments (ignored).
#'
#' @method print rnaROC
#' @export
#' @rdname print.rnaROC

print.rnaROC <- function(x, ...) {
  message("==============================")
  message("Object returned by 'rna.roc()'")
  message("==============================")
  message("ROC method used: ", x$method)
  message("Genes analyzed: ", paste(x$genes, collapse=", "))
  message("AUC: ", round(as.numeric(x$auc), 3))
  message("You can visualize the ROC curve using plot(x$roc_object)")
  invisible(x)
}

# --------------------------------------------
#' Print Method for rnaVolcano Objects
#'
#' Custom print method for objects returned by \code{rna.volcano()}.
#'
#' @param x An object returned by \code{rna.volcano()}.
#' @param ... Additional arguments (ignored).
#'
#' @method print rnaVolcano
#' @export
#' @rdname print.rnaVolcano

print.rnaVolcano <- function(x, ...) {
  message("==============================")
  message("Object returned by 'rna.volcano()'")
  message("==============================")
  message("Upregulated genes: ", length(x$up))
  message("Downregulated genes: ", length(x$down))
  message("Full results available in x$full_results")
  message("Volcano plot is stored in x$plot (if generated)")
  invisible(x)
}
