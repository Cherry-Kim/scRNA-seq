BiocManager::install(c("SingleR","Seurat","celldex","scran","SingleCellExperiment") )
library(scran)
library(celldex)
library(Seurat)
library(SingleR)
library(SingleCellExperiment)

seurat_obj <-  readRDS('colon.combined.rds')

clusters <- seurat_obj$seurat_clusters
test <- as.SingleCellExperiment(seurat_obj)

hpca.se <- celldex::HumanPrimaryCellAtlasData()
bpe <- celldex::BlueprintEncodeData()

pred.hesc <- SingleR(test = test, ref = hpca.se, genes = "de", de.method = "wilcox", labels = hpca.se$label.main, clusters = clusters)
pred.bpe <- SingleR(test = test, ref = bpe, genes = "de", de.method = "wilcox", labels = bpe$label.main, clusters = clusters)

hpca.clusters <- clusters
levels(hpca.clusters) <- pred.hesc$labels
bpe.clusters <- clusters
levels(bpe.clusters) <- pred.bpe$labels
#pred.hesc$labels
#unique(pred.hesc$labels)
#pred.bpe$labels
#unique(pred.bpe$labels)

hpca.clusters <- gsub("0","Pre-B_cell_CD34-",hpca.clusters)
hpca.clusters <- gsub("1|3|7|14","NK_cell",hpca.clusters)

bpe.clusters <- gsub("0","CD4+ T-cells",bpe.clusters)
bpe.clusters <- gsub("1|3|7|14","CD8+ T-cells",bpe.clusters)

seurat_obj[["SingleR.labels"]] <- hpca.clusters

seurat_obj.bpe <- seurat_obj
seurat_obj.bpe[["SingleR.labels"]] <- bpe.clusters

Idents(seurat_obj) <- seurat_obj[["SingleR.labels"]]
Idents(seurat_obj.bpe) <- seurat_obj.bpe[["SingleR.labels"]]

# Visualization
DimPlot(seurat_obj, reduction = "umap")
DimPlot(seurat_obj.bpe, reduction = "umap")

