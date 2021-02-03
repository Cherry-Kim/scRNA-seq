#!/usr/bin/env Rscript
library(dplyr)
library(Seurat)
library(cowplot)

colon.combined <-  readRDS('colon.combined.rds')
# Visualization
p1 <- DimPlot(colon.combined, reduction = "umap", group.by = "stim")
p2 <- DimPlot(colon.combined, reduction = "umap", label = TRUE)
plot_grid(p1, p2)
DimPlot(colon.combined, reduction = "umap", split.by = "stim")

### STEP4. Identify conserved cell type markers
DefaultAssay(colon.combined) <- "RNA"
nk.markers <- FindConservedMarkers(colon.combined, ident.1 = 7, grouping.var = "stim", verbose = FALSE)
markers.to.plot=c("AKNA", "BCL11A", "BTG2", "EAF2", "GATA1", "GFI1B", "IKZF1","IKZF3","IRF5","IRF8","KLF1","KLF2","MYB","NFE2","POU2F2","RHOXF2","SPI1","SPIB","TAL1","TCF7")

colon.combined <- RenameIdents(seurat_obj, `0`="CD4+ T-cells", `1`="CD8+ T-cells", `3`="CD8+ T-cells",`7`="CD8+ T-cells",`14`="CD8+ T-cells", '2'="B-cells","5"="B-cells", "4"="Epithelial cells","10"="Epithelial cells","11"="Epithelial cells","12"="Epithelial cells","15"="Epithelial cells","21"="Epithelial cells","8"="Fibroblasts","17"="Fibroblasts","19"="Fibroblasts","9"="Monocytes","16"="Monocytes","13"="Endothelial cells","22"="Endothelial cells","18"="HSC","20"="Astrocytes","6"="Macrophages")
DimPlot(colon.combined, label = TRUE)
DimPlot(colon.combined, reduction = "umap", split.by = "stim", label = TRUE)

DotPlot(colon.combined, features = markers.to.plot, cols=c("blue", "red"), dot.scale = 8) + RotatedAxis()
library(ggplot2)
t.cells <- subset(colon.combined, idents = "Epithelial cells")
Idents(t.cells) <- "stim"
avg.t.cells <- log1p(AverageExpression(t.cells, verbose = FALSE)$RNA)
avg.t.cells$gene <- rownames(avg.t.cells)
genes.to.label = c("IRF8","MYB","TCF7")
p1 <- ggplot(avg.t.cells, aes(CTRL, tumor)) + geom_point() + ggtitle("Epithelial cells")
LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, col="red")

FeaturePlot(colon.combined, features = markers.to.plot )
FeaturePlot(colon.combined, features = markers.to.plot, split.by = "stim", max.cutoff = 3, cols = c("grey", "red"))

plots <- VlnPlot(colon.combined, features = markers.to.plot, split.by = "stim", pt.size = 0, combine = FALSE)
CombinePlots(plots = plots, ncol = 1)
