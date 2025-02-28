### ### Spatial Data Analysis using Seurat Xenium Human Breast Cancer Dataset
library(Seurat)
library(tidyverse)

### 1. Download and unzip the Xenium human-breast data
# www.10xgenomics.com/products/xenium-in-situ/preview-dataset-human-breast

### 2. Load and preprocess the data
# Load the data
xenium.obj <- LoadXenium(
  "~/Desktop/Projects/V1.0/Spatial-scRNA-Seq/Xenium Human Breast Cancer Dataset/data/outs/", 
  fov = "fov", assay = "Xenium")


# remove cells with 0 counts
xenium.obj <- subset(xenium.obj, subset = nCount_Xenium > 0)

# Genes/cell (nFeature_Xenium) and transcript counts/cell (nCount_Xenium)
VlnPlot(xenium.obj, features = c("nFeature_Xenium", "nCount_Xenium"), 
        ncol = 2, pt.size = 0)

# 166363 cells X 313 genes
xenium.obj <- subset(xenium.obj, 
                     subset = nFeature_Xenium > 5 & nFeature_Xenium < 200 & 
                       nCount_Xenium > 10 & nCount_Xenium < 1000)

# If i wanted to 163779 cells X 313 genes

### 3. SCTransform normalization followed by standard Seurat workflow analysis
xenium.obj <- SCTransform(xenium.obj, assay = "Xenium")
xenium.obj <- RunPCA(xenium.obj, npcs = 30, features = rownames(xenium.obj))
xenium.obj <- RunUMAP(xenium.obj, dims = 1:30)
xenium.obj <- FindNeighbors(xenium.obj, reduction = "pca", dims = 1:30)
xenium.obj <- FindClusters(xenium.obj, resolution = 0.2)
DimPlot(xenium.obj, group.by = "seurat_clusters", label = T)

### 4. Cell type annotation
AllMarkers <- FindAllMarkers(xenium.obj)

rownames(xenium.obj)

# B cells: MS4A1, CD79A, macrophages: ITGAX, T cells: CD3E
FeaturePlot(object = xenium.obj, features=c("MS4A1", "CD79A", "ITGAX", "CD3E"))

# T cells: CD4, CD8A, NK cells: NKG7, Mast cells: KIT
FeaturePlot(object = xenium.obj, features=c("CD4", "CD8A", "NKG7", "KIT"))

# Endothelial: PECAM1, myoepithelial: KRT15, Fibroblasts:LUM, Prolif: MKI67
FeaturePlot(object = xenium.obj, features=c("PECAM1", "KRT15", "LUM", "MKI67"))

# ductal carcinoma: CEACAM6, invasive tumor:FASN
FeaturePlot(object = xenium.obj, features=c("CEACAM6", "FASN"))

xenium.obj <- RenameIdents(xenium.obj, 
                           `0` = "invasive", `2` = "ductal", `11` = "ductal", `6` = "myoepithelial", 
                           `3` = "Fibroblasts", `4` = "macrophages", `5` = "Endothelial", 
                           `1` = "T/NK cells", `10` = "T/NK cells", `7` = "B cells", `8` = "B cells", 
                           `9` = "Mast cells")

xenium.obj$cell_type <- Idents(xenium.obj)

DimPlot(xenium.obj, group.by = "cell_type", label = T)

### 5. visualization on spatial coordinates
ImageDimPlot(xenium.obj, group.by = "seurat_clusters", size = 0.9)
ImageDimPlot(xenium.obj, group.by = "cell_type", size = 0.9)

# visualize the expression level of genes at the per-cell level
ImageFeaturePlot(xenium.obj, features = "FASN", cols = c("white", "red"))
# Visualize the expression level of genes at the per-cell level with different color

ImageFeaturePlot(xenium.obj, features = "CEACAM6", cols = c("blue", "red"))

# Save the plots
ggsave("Xenium_Human_Breast_Cancer_Dataset_genes_per-cell.png", width = 5, height = 5)


# Plot tumor specific marker genes
ImageDimPlot(xenium.obj, fov = "fov", nmols = 20000, axes = T,
             molecules = c("CEACAM6", "FASN"))

# Save the plots generated at specific location 
ggsave("Xenium_Human_Breast_Cancer_Dataset.png", width = 5, height = 5)


