# Xenium Human Breast Cancer Spatial Transcriptomics Analysis

## Overview
This repository contains a comprehensive spatial transcriptomics analysis of a human breast cancer dataset generated using the 10x Genomics Xenium platform. The workflow identifies cell types, spatial patterns, and marker genes associated with breast cancer microenvironments.

## Data Source
[10x Genomics Xenium In-Situ Preview Dataset](https://www.10xgenomics.com/products/xenium-in-situ/preview-dataset-human-breast)

## Analysis Workflow
1. **Data Preprocessing & Quality Control**
   - Filtering cells based on feature count and transcript count thresholds
   - Removal of low-quality cells and outliers

2. **Normalization**
   - SCTransform method for variance-stabilizing normalization

3. **Dimensionality Reduction**
   - PCA computation of top 30 principal components
   - UMAP projection for visualization

4. **Clustering**
   - Graph-based clustering with resolution 0.2
   - 11 distinct cell clusters identified

5. **Cell Type Annotation**
   - Differential expression analysis
   - Marker-based annotation of clusters
   - Identification of tumor, immune and stromal populations

6. **Spatial Visualization**
   - Integration of gene expression with tissue coordinates
   - Mapping of cell types and key marker genes in spatial context

## Key Findings
- Distinct spatial organization of invasive carcinoma cells (FASN+) and ductal carcinoma cells (CEACAM6+)
- Immune cell infiltration throughout the tumor microenvironment
- Concentrated expression of FASN in invasive tumor regions
- Ductal structures marked by CEACAM6 expression

## Requirements
```r
library(Seurat)
library(tidyverse)
library(SeuratObject)
```

## References
1. Swinnen, J. V., et al. (2006). "Fatty acid synthase drives the growth of prostate cancer cells." *Cancer Research*
2. Blumenthal, R. D., et al. (2007). "Carcinoembryonic antigen (CEA) and CEACAM6 in cancer progression." *Cancer Biology & Therapy*
3. Hanahan, D., & Weinberg, R. A. (2011). "Hallmarks of cancer: The next generation." *Cell*
