# autornaR

<!-- badges: start --> <!-- badges: end -->

**autornaR** is an R package designed to streamline RNA-seq data analysis, from raw count import to normalization, quality control, statistical testing, and visualization.  
It provides functions for filtering, imputation, outlier detection, gene ranking, ROC curves, dimensionality reduction, correlation, and differential expression analysis.

**Note:** The package is primarily optimized for pairwise analyses (e.g., Control vs. Treatment).  
Workflows involving multiple groups are possible, but most functions are tuned for two-group comparisons yet.

## Installation

You can install the development version directly from GitHub:

```r
# Install remotes if not already installed
install.packages("remotes")

# Install autornaR from GitHub
remotes::install_github("Luiz-Garcia-R/autornaR")
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
rna.project() - Create a new rna_project object to store the results from downstream analyses
rna.import() - Import and validate raw RNA-seq count data (raw_data + metadata)
rna.normalize() - Normalize data, filter genes, impute missing values, remove outliers
rna.qc() - Performs comprehensive quality control (QC) on normalized RNA-seq data
- rna.save() - Extracts and save the essential components from a object

# Exploratory and Differential Analysis
For deeper insights and pairwise comparisons:
  - rna.compare() - Performs differential expression analysis between two groups
  - rna.corr() - Compute correlations among samples or experimental groups
  - rna.enrich() - performs Gene Ontology (GO) enrichment analysis
  - rna.degs() - Identify differentially expressed genes
  - rna.dimred() - Perform PCA and UMAP for dimensionality reduction
  - rna.heatmap() - Visualize top variable genes with heatmaps
  - rna.gsea() - Perform a GSEA rank to visualize top enriched pathways
  - rna.network() - Build a gene-gene network using pathways identified by rna.gsea()
  - rna.roc() - Assess discriminatory power of specific genes
  - rna.boxplot() - Perform Limma t-tests for gene-level differences
  - rna.sets() - Evaluate differential expressed genes distribution
  - rna.volcano() - Generate volcano plots for differential expression analysis

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
my_project <- rna.project()

# Import counts
my_project <- rna.import(my_project,
                        raw_data = raw_counts,
                        metadata = meta_df)

# Normalize
my_project <- rna.normalize(my_project)

# Quality control
my_project <- rna.qc(my_project)

# Dimension reduction
my_project <- rna.dimred(my_project)

# Differential expression
my_project <- rna.compare(my_project)

# Evaluating top DEGs
my_project <- rna.degs(my_project)

```

## Contact

For questions, suggestions, or contributions, open an issue or pull request on GitHub:
https://github.com/Luiz-Garcia-R/autornaR

Thank you for using autornaR!
