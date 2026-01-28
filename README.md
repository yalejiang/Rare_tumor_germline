# Rare Tumor Germline Mutation Analysis and Visualization Pipeline
A comprehensive R pipeline for **analysis and visualization of germline mutation landscape in rare tumors**, including rare cancer type distribution, germline mutation landscape heatmap, KEGG enrichment analysis of germline mutated genes, and somatic second hit mutation analysis. This pipeline generates four publication-quality core figures for rare tumor oncology bioinformatics research.

## Project Overview
This script is designed for the **germline/somatic mutation data analysis of rare tumors**, with the following key analyses (all results are publication-ready):
1. Rare cancer type sample distribution and germline mutation rate calculation
2. Germline mutation landscape heatmap (mutation type/rate + clinical annotation)
3. KEGG functional enrichment analysis of rare tumor germline mutated genes (dot + sankey combined plot)
4. Somatic second hit mutation distribution in rare tumors (by germline gene/rare cancer type)

All output figures are **TIFF format (300dpi)** with LZW compression, which is the standard format for academic publication (small file size + high resolution).

## Environment Setup
### R Version Requirement
R ≥ 4.2.0 (tested on R 4.3.1)

### Install Dependencies
Run the following code in R to install all required packages (the script will also auto-install missing packages):
```r
required_packages <- c(
  "tidyverse", "readxl", "ggprism", "viridis", "scico", "paletteer",
  "ggsci", "RColorBrewer", "patchwork", "scales", "purrr", "writexl",
  "maftools", "ComplexHeatmap", "circlize", "grid", "org.Hs.eg.db",
  "clusterProfiler", "conflicted", "ggsankey"
)

invisible(lapply(required_packages, function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
    library(pkg, character.only = TRUE)
  }
}))


Data Preparation
Input Files
Place all input files in the project root directory (same as the script oncotree_mutation_vis.R). The required input files are all related to rare tumor germline mutation data:
File Name	Data Content
Figure1.rds	Contains clin_info (rare tumor clinical data), cancer_level (rare cancer type factor), total_counts (total sample counts)
figure2.rds	Contains mat (germline mutation matrix), clin_info (clinical data), mut_data (gene germline mutation data)
figure4.rdata	Contains somic_data (somatic second hit mutation data), germ_data (rare tumor germline mutation data)
gene_info.xlsx	Rare tumor germline mutated genes, with a single column named Gene (HGNC standard symbol, no NA values)

