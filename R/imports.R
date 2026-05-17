# R/imports.R

#' @keywords internal
# --- Imported functions ---

# ======================
# Base / Stats
# ======================

#' @importFrom stats aggregate as.formula var cor median quantile rnorm
#'   sd setNames reorder shapiro.test t.test aov prcomp glm predict
#'   p.adjust wilcox.test

# ======================
# Utils
# ======================

#' @importFrom utils globalVariables head packageVersion
#' @importFrom data.table := as.data.table

# ======================
# Graphics / Devices
# ======================

#' @importFrom graphics abline grid mtext par text plot
#' @importFrom grDevices colorRampPalette

# ======================
# ggplot2 ecosystem
# ======================

#' @importFrom ggplot2 ggplot aes geom_point geom_text geom_bar geom_boxplot
#'   geom_jitter geom_density geom_line scale_color_manual scale_fill_manual
#'   labs theme theme_minimal element_text element_rect facet_wrap

#' @importFrom ggrepel geom_text_repel

# ======================
# Graph / Network
# ======================

#' @importFrom igraph graph_from_data_frame
#' @importFrom ggraph ggraph geom_edge_link geom_node_point geom_node_text
#' @importFrom tidygraph as_tbl_graph activate

# ======================
# Tidyverse helpers
# ======================

#' @importFrom tidyr pivot_longer
#' @importFrom dplyr %>% left_join filter

# ======================
# Bioinformatics
# ======================

#' @importFrom pheatmap pheatmap
#' @importFrom pROC roc auc
#' @importFrom ComplexUpset upset

# ======================
# Plot composition
# ======================

#' @importFrom patchwork plot_layout

NULL
