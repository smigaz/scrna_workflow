# BiocManager::install("rhdf5r")
library(dplyr)
library(hdf5r)
library(Seurat)
library(patchwork)
library(ggplot2)

# Load the dataset
counts <- Read10X_h5(filename = "./20k_PBMC_3p_HT_nextgem_Chromium_X_raw_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
pbmc_seurat_odj <- CreateSeuratObject(counts = counts, project = "PBMC", min.cells = 3, min.features = 200)
str(pbmc_seurat_odj)

# QC
pbmc_seurat_odj[["percent.mt"]] <- PercentageFeatureSet(pbmc_seurat_odj, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(pbmc_seurat_odj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(pbmc_seurat_odj, feature1 = "nCount_RNA", feature2 = "percent.mt") + geom_smooth(method = "lm")
plot2 <- FeatureScatter(pbmc_seurat_odj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = "lm")
plot1 + plot2

pbmc_seurat_odj <- subset(pbmc_seurat_odj, subset = nFeature_RNA > 200 & nFeature_RNA < 2000 & percent.mt < 5)

plot1_filtered <- FeatureScatter(pbmc_seurat_odj, feature1 = "nCount_RNA", feature2 = "percent.mt") + geom_smooth(method = "lm")
plot2_filtered <- FeatureScatter(pbmc_seurat_odj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = "lm")
plot1_filtered + plot2_filtered

VlnPlot(pbmc_seurat_odj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

pbmc_seurat_odj <- NormalizeData(pbmc_seurat_odj, normalization.method = "LogNormalize", scale.factor = 10000)
# pbmc_seurat_odj <- NormalizeData(pbmc_seurat_odj)

pbmc_seurat_odj <- FindVariableFeatures(pbmc_seurat_odj, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc_seurat_odj), 10)
top10

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc_seurat_odj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

all.genes <- rownames(pbmc_seurat_odj)
pbmc_seurat_odj <- ScaleData(pbmc_seurat_odj, features = all.genes)

# linear dimensional reduction
pbmc_seurat_odj <- RunPCA(pbmc_seurat_odj, features = VariableFeatures(object = pbmc_seurat_odj))
# Examine and visualize PCA results a few different ways
print(pbmc_seurat_odj[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(pbmc_seurat_odj, dims = 1:2, reduction = "pca")
DimPlot(pbmc_seurat_odj, reduction = "pca")

DimHeatmap(pbmc_seurat_odj, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc_seurat_odj, dims = 1:15, cells = 500, balanced = TRUE)

# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
pbmc <- JackStraw(pbmc_seurat_odj, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc_seurat_odj, dims = 1:20)

ElbowPlot(pbmc_seurat_odj)

pbmc_seurat_odj <- FindNeighbors(pbmc_seurat_odj, dims = 1:15)
pbmc_seurat_odj <- FindClusters(pbmc_seurat_odj, resolution = c(0.1,0.3,0.5,0.7,1))

DimPlot(pbmc_seurat_odj, group.by = "RNA_snn_res.0.1",label =TRUE)
Idents(pbmc_seurat_odj) <- "RNA_snn_res.0.1"
# Look at cluster IDs of the first 5 cells
head(Idents(pbmc_seurat_odj), 5)

# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
pbmc_seurat_odj <- RunUMAP(pbmc_seurat_odj, dims = 1:15)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pbmc_seurat_odj, reduction = "umap")
