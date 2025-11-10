# R/imports.R

#' @keywords internal
# --- Imported functions ---

# stats
#' @importFrom stats aggregate as.formula var cor median quantile rnorm sd setNames reorder shapiro.test t.test aov prcomp glm predict p.adjust wilcox.test

# graphics
#' @importFrom graphics abline grid mtext par text plot

# utils
#' @importFrom utils globalVariables head

# ggplot2
#' @importFrom ggplot2 ggplot aes geom_point geom_text geom_bar geom_boxplot geom_jitter geom_density geom_line scale_color_manual scale_fill_manual labs theme theme_minimal element_text element_rect facet_wrap

# tidyr
#' @importFrom tidyr pivot_longer

# dplyr
#' @importFrom dplyr %>% left_join filter

# pheatmap
#' @importFrom pheatmap pheatmap

# pROC
#' @importFrom pROC roc auc

# biomaRt
#' @importFrom biomaRt useEnsembl getBM

# patchwork
#' @importFrom patchwork plot_layout

# ComplexUpset
#' @importFrom ComplexUpset upset

# ggrepel
#' @importFrom ggrepel geom_text_repel

NULL
