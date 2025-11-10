# autornaR

<!-- badges: start --> <!-- badges: end -->

**autornaR** is an R package designed to streamline RNA-seq data analysis, from raw count import to normalization, quality control, statistical testing, and visualization.  
It provides functions for filtering, imputation, outlier detection, gene ranking, ROC curves, dimensionality reduction, correlation, and differential expression analysis.

**Note:** The package is primarily optimized for pairwise analyses (e.g., Control vs. Treatment).  
Workflows involving multiple groups are possible, but most functions are tuned for two-group comparisons.

## Installation

You can install the development version directly from GitHub:

```r
# Install devtools if not already installed
install.packages("devtools")

# Install autornaR from GitHub
devtools::install_github("Luiz-Garcia-R/autornaR")
```

## Input Data Format (Very Important)

To use autornaR, you need two input data frames:

  1. raw_data – raw RNA-seq counts (from HISAT2, STAR, featureCounts, tximport, or clean matrices)
    - Must contain a first column with gene IDs
    - Other columns must contain numeric counts for samples

Example of raw_data:

| GeneID     | Sample1 | Sample2 | Sample3 | Sample4 |
| ---------- | ------- | ------- | ------- | ------- |
| ENSG000001 | 12      | 8       | 15      | 20      |
| ENSG000002 | 0       | 3       | 1       | 2       |
| ENSG000003 | 45      | 50      | 47      | 52      |
| ENSG000004 | 5       | 2       | 3       | 1       |


  2. metadata – describes your samples and experimental groups
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
rna.import() – Import and validate raw RNA-seq count data (raw_data + metadata)
rna.normalize() – Normalize data, filter genes, impute missing values, remove outliers
rna.qc() – Generate QC plots including boxplots, PCA, and density distributions

# Exploratory and Differential Analysis
For deeper insights and pairwise comparisons:
  - rna.compare() - Performs differential expression analysis between two groups
  - rna.corr() – Compute correlations among samples or experimental groups
  - rna.enrich() - performs Gene Ontology (GO) enrichment analysis
  - rna.identify() - Annotates gene biotypes and classify expressed transcripts
  - rna.rank() – Rank genes by mean or variance within each group
  - rna.degs() – Identify differentially expressed genes
  - rna.dimred() – Perform PCA and UMAP for dimensionality reduction
  - rna.heatmap() – Visualize top variable genes with heatmaps
  - rna.roc() – Assess discriminatory power of specific genes
  - rna.ttest() – Perform t-tests or Mann-Whitney tests for gene-level differences
  - rna.uset() – Visualize overlaps between gene sets with Upset Plots
  - rna.volcano() – Generate volcano plots for differential expression analysis

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

metadata <- data.frame(
  Sample = c("Control1","Control2","Treatment1","Treatment2"),
  Group  = c("Control","Control","Treatment","Treatment")
)

# Import counts
imp_data <- rna.import(raw_counts, metadata)

# Normalize
normalized_data <- rna.normalize(imp_data)

# QC plots
rna.qc(normalized_data)

# Correlation matrix
corr_mat <- rna.corr(normalized_data, metadata)

# Differential expression
rna.dems(normalized_data, metadata)

# Gene ranking
ranked_genes <- rna.rank(normalized_data, metric = "mean", top_n = 20)
```

## Contact

For questions, suggestions, or contributions, open an issue or pull request on GitHub:
https://github.com/Luiz-Garcia-R/autornaR

Thank you for using autornaR!
