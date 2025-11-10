# R/globals.R
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    # Common variables
    "Sample", "Group", "Expression", "Library_size", "Genes_detected", "Outlier",
    "Gene", "Mean", "Var", "Rank", "hgnc_symbol", "GeneLabel",
    "pred", "classe", "roc_object", "auc", "method",
    "Value", "signif_label", "y_max", "y_pos",
    "padj", "log2FoldChange", "group_color",
    "gene_map", "ensembl_gene_id", "Symbol", "GeneID", "desc",

    # UpSet / helpers
    "presence_df", "presence_list", "intersect",

    # ggplot / tidyverse helpers
    ".data",

    # PCA / dimensionality reduction
    "PC1", "PC2", "Comp1", "Comp2",

    # tSNE / UMAP
    "UMAP1", "UMAP2", "tSNE1", "tSNE2",

    # Misc and identifiers
    "i", "g1", "g2", "g1_samples", "g2_samples", "genes_use", "highlight_df",

    # From rna_identify()
    "transcript_biotype", "description", "gene_biotype", "original_id",
    "expressed", "gene_biotype_plot", "prop", "total",

    # Differential expression / volcano
    "Regulation", "pvalue",

    # Additional helpers (if any)
    "Comp1", "Comp2"
  ))
}


