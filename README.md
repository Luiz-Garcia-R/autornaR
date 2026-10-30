# autornaR

<!-- badges: start --> <!-- badges: end -->

<table>
<tr>
<td width="70%" valign="top">

**RNA-seq analysis, from data to biology.**

autornaR is an R package designed to streamline RNA-seq data analysis, from raw count import to normalization, quality control, statistical testing, and visualization.

It provides functions for filtering, imputation, outlier detection, gene ranking, ROC curves, dimensionality reduction, correlation, and differential expression analysis.

</td>

<td width="30%" align="center" valign="middle">

<img src="man/figures/logo_autornaR.png" width="160">

</td>
</tr>
</table>

## Installation

You can install the development version directly from GitHub:

```r
# Install remotes if not already installed
install.packages("remotes")

# Install autornaR from GitHub
remotes::install_github("Luiz-Garcia-R/autornaR")
```

## Getting started

New to autornaR? **Start with the [Getting started with autornaR vignette](https://luiz-garcia-r.github.io/autornaR/articles/overview.html)** for a guided introduction to the package's main workflow, from data import and normalization to quality control and downstream analyses.

The online vignette is the easiest way to get started, with step-by-step explanations and code examples.

You can also access the vignette directly from R:

```r
browseVignettes("autornaR")
```

Or open the introductory vignette specifically:

```r
vignette("overview", package = "autornaR")
```



## Input Data Format (Very Important)

To use autornaR, you need two input data frames:

  1. raw_data - raw RNA-seq counts (from HISAT2, STAR, featureCounts, tximport, or clean matrices)
    - Must contain a first column with gene IDs
    - Other columns must contain numeric counts for samples

Example of raw_data:

| GeneID     | Sample1 | Sample2 | Sample3 | Sample4 |
| ---------- | ------- | ------- | ------- | ------- |
| ENSG000001 | 12      | 8       | 15      | 20      |
| ENSG000002 | 0       | 3       | 1       | 2       |
| ENSG000003 | 45      | 50      | 47      | 52      |
| ENSG000004 | 5       | 2       | 3       | 1       |


  2. metadata - describes your samples and experimental groups
    Must contain at least two columns:
    - Sample: matching exactly the sample columns in raw_data
    - Group: experimental condition

Example of metadata:

| Sample  | Group     |
| ------- | --------- |
| Sample1 | Control   |
| Sample2 | Control   |
| Sample3 | Treatment |
| Sample4 | Treatment |

## Main Workflow

The recommended workflow guides users from raw count import to quality control and initial gene evaluation:
- rna_project() - Create a new rna_project object to store the results from downstream analyses
- rna_import() - Import and validate raw RNA-seq count data (raw_data + metadata)
- rna_normalize() - Normalize data, filter genes, impute missing values, remove outliers
- rna_qc() - Performs comprehensive quality control (QC) on normalized RNA-seq data
- rna_save() - Extracts and save the essential components from a object

# Exploratory and Differential Analysis
For deeper insights and pairwise comparisons:
  - rna_compare() - Performs differential expression analysis (Limma or DESeq2) between two groups
  - rna_association() - Compute gene or clinical features correlations among samples
  - rna_enrich() - performs Gene Ontology (GO) enrichment analysis
  - rna_degs() - Identify differentially expressed genes
  - rna_dimred() - Perform PCA and UMAP for dimensionality reduction
  - rna_heatmap() - Visualize top variable genes with heatmaps
  - rna_gsea() - Perform a GSEA rank to visualize top enriched pathway
  - rna_gsva() - Perform a GSVA to visualize individual enriched pathway
  - rna_network() - Build a gene-gene network using pathways identified by rna_gsea()
  - rna_roc() - Assess discriminatory power of specific genes
  - rna_boxplot() - Perform Limma t-tests for gene-level differences
  - rna_sets() - Evaluate differential expressed genes distribution
  - rna_volcano() - Generate volcano plots for differential expression analysis
  - rna_wgcna() - Performs a weighted correlation network analysis

## Example Workflow (Minimal)
```r
# Example raw counts
raw_counts <- data.frame(
  GeneID = c("ENSG000001","ENSG000002","ENSG000003"),
  Control1 = c(12,0,45),
  Control2 = c(8,3,50),
  Treatment1 = c(15,1,47),
  Treatment2 = c(20,2,52)
)

meta_df <- data.frame(
  Sample = c("Control1","Control2","Treatment1","Treatment2"),
  Group  = c("Control","Control","Treatment","Treatment")
)

# Create a new rna_project object
my_project <- rna_project()

# Import counts
my_project <- rna_import(my_project,
                        raw_data = raw_counts,
                        metadata = meta_df)

# Normalize
my_project <- rna_normalize(my_project)

# Quality control
my_project <- rna_qc(my_project)

# Dimension reduction
my_project <- rna_dimred(my_project)

# Differential expression
my_project <- rna_compare(my_project,
                          comparison = c("Treatment", "Control"))

# Evaluating top DEGs
my_project <- rna_degs(my_project)

# Performs a WGCNA
my_project <- rna_wgcna(my_project, power = 6)
```

**Note:** The package is primarily optimized for pairwise analyses (e.g., Control vs. Treatment).  
Workflows involving multiple groups are possible, but most functions involving measures like log2FC are tuned for two-group comparisons yet.

## Contact

For questions, suggestions, or contributions, open an issue or pull request on GitHub:
https://github.com/Luiz-Garcia-R/autornaR

Thank you for using autornaR!
