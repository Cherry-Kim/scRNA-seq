#!/usr/bin/env Rscript
library(dplyr)
library(Seurat)
library(cowplot)

### STEP1. Setup the Seurat objects
Normal.data <- Read10X(data.dir = "/hykim/colon_N/outs/filtered_feature_bc_matrix")
Tumor.data <- Read10X(data.dir = "/hykim/colon_T/outs/filtered_feature_bc_matrix")
# Set up control object
Normal <- CreateSeuratObject(counts = Normal.data, project = "COLON_Normal", min.cells = 5, min.features = 200)
Normal$stim <- "Normal"
Normal[["percent.mt"]] <- PercentageFeatureSet(object = Normal, pattern = "^MT-")
Normal <- subset(Normal, subset = percent.mt < 20)
Normal <- NormalizeData(Normal, normalization.method = "LogNormalize", scale.factor = 10000)
Normal <- FindVariableFeatures(Normal, selection.method = "vst", nfeatures = 2000)
# Set up Tumor object
Tumor <- CreateSeuratObject(counts = Tumor.data, project = "COLON_TUMOR", min.cells = 5, min.features = 200)
Tumor$stim <- "Tumor"
Tumor[["percent.mt"]] <- PercentageFeatureSet(Tumor, pattern = "^MT-")
Tumor <- subset(Tumor, subset = percent.mt < 20)
Tumor <- NormalizeData(Tumor, normalization.method = "LogNormalize", scale.factor = 10000)
Tumor <- FindVariableFeatures(Tumor, selection.method = "vst", nfeatures = 2000)

### STEP2. Perform integration
colon.anchors <- FindIntegrationAnchors(object.list = list(Normal, Tumor), dims = 1:20)
colon.combined <- IntegrateData(anchorset = colon.anchors, dims = 1:20)

### STEP3. Perform an integrated analysis
DefaultAssay(colon.combined) <- "integrated"
#Run the standard workflow for visualization and clustering
all.genes <- rownames(colon.combined)
colon.combined <- ScaleData(colon.combined, features = all.genes)
colon.combined <- RunPCA(colon.combined, features = VariableFeatures(object = colon.combined) )
# t-SNE and Clustering
colon.combined <- RunUMAP(colon.combined, reduction = "pca", dims = 1:20)
colon.combined <- FindNeighbors(colon.combined, reduction = "pca", dims = 1:20)
colon.combined <- FindClusters(colon.combined, resolution = 0.5)

saveRDS(colon.combined, file = "colon.combined.rds")

