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
markers.to.plot=c("AKNA")

colon.combined <- RenameIdents(seurat_obj, `0`="CD4+ T-cells", `1`="CD8+ T-cells")
DimPlot(colon.combined, label = TRUE)
DimPlot(colon.combined, reduction = "umap", split.by = "stim", label = TRUE)

DotPlot(colon.combined, features = markers.to.plot, cols=c("blue", "red"), dot.scale = 8) + RotatedAxis()
library(ggplot2)
t.cells <- subset(colon.combined, idents = "Epithelial cells")
Idents(t.cells) <- "stim"
avg.t.cells <- log1p(AverageExpression(t.cells, verbose = FALSE)$RNA)
avg.t.cells$gene <- rownames(avg.t.cells)
genes.to.label = c("IRF8","MYB")
p1 <- ggplot(avg.t.cells, aes(CTRL, tumor)) + geom_point() + ggtitle("Epithelial cells")
LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, col="red")

FeaturePlot(colon.combined, features = markers.to.plot )
FeaturePlot(colon.combined, features = markers.to.plot, split.by = "stim", max.cutoff = 3, cols = c("grey", "red"))

plots <- VlnPlot(colon.combined, features = markers.to.plot, split.by = "stim", pt.size = 0, combine = FALSE)
CombinePlots(plots = plots, ncol = 1)
