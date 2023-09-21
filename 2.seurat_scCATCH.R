library(Matrix)
matrix_dir = "Cancer/Colon/scRNA-seq/count/colon_T/outs/filtered_feature_bc_matrix/"
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path, header = FALSE, stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, header = FALSE, stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V1
#################################################

library(dplyr)
library(Seurat)

# import input_data(cellranger count output)
pbmc.data <- Read10X(data.dir = "/Cancer/Colon/scRNA-seq/count/colon_T/outs/filtered_feature_bc_matrix")
# setup the seurat object
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "colon_T")
#pbmc <- CreateSeuratObject(counts = pbmc.data, project = "colon", min.cells = 3, min.features = 200)

### STEP0-1. Preprocessing(QC)  : Low quality cell, mitochondria genome percent check
pbmc[["percent.mt"]] <- PercentageFeatureSet(object = pbmc, pattern = "^MT-")
#we have as much as 97.5% mitochondrial reads in some of the cells, so it is quite unlikely that there is much celltype signature left in those.
pbmc <- subset(pbmc, subset = percent.mt < 20)
#plot1 <- FeatureScatter(pbmc, feature1 = "percent.mt", feature2 = "nCount_RNA")

### STEP0-2. Normalization
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

### STEP1. Feature Selection 
# identification of highly variable features
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 2000)
# scaling data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(object = pbmc, features = all.genes)

### STEP2. Dimension Reduction
# perform linear dimensional reduction
pbmc <- RunPCA(object = pbmc, features = VariableFeatures(object = pbmc))
#DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
# determine the "dimensionality" of the datasetdetermine the "dimensionality" of the dataset
#pbmc <- JackStraw(object = pbmc, num.replicate = 100)
#pbmc <- ScoreJackStraw(object = pbmc, dims = 1:20)

### STEP3. Clustering & Visualization
#cluster the cells
pbmc <- FindNeighbors(object = pbmc, dims = 1:10)
pbmc <- FindClusters(object = pbmc, resolution = 0.5)
# Run non-linear dimensional reduction(tSNE)
pbmc <- RunTSNE(object=pbmc)
pbmc <- RunUMAP(object = pbmc, dims = 1:10)
plot1 <- DimPlot(object = pbmc, reduction = "umap")
plot2 <- DimPlot(object = pbmc, reduction = "tsne")
p2 <- plot1 + plot2
png("PCA.png")
p2
dev.off()
saveRDS(pbmc, file = "colon_T.rds")
# Finding differentially expressed features(biomarkers)
install.packages("magrittr")
library(magrittr)

pbmc.markers <- FindAllMarkers(object = pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(pbmc.markers, "pbmc.markers.txt", col.names=T, row.names=T, quote=F, sep='\t')
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(pbmc, features = top10$gene) + NoLegend()
#Assigning cell type identity to clusters
new.cluster.ids <- c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono", 
    "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
saveRDS(pbmc, file = "colon_T.final.rds")
pbmc <- readRDS('colon_T.final.rds')
##################################################

library(scCATCH)
#cluster marker genes identification
clu_markers <- findmarkergenes(pbmc,
                               species = 'Human',
                               cluster = 'All',
                               match_CellMatch = TRUE,
                               cancer = c('Colon Cancer','Colorectal Cancer'),
                               tissue = c('Colon','Colorectum'),
                               cell_min_pct = 0.25,
                               logfc = 0.25,
                               pvalue = 0.05)
#https://github.com/ZJUFanLab/scCATCH/wiki/human_tissues
#Cluster annotation
clu_ann <- scCATCH(clu_markers,
                    species = 'Human',
		    cancer = c('Colon Cancer','Colorectal Cancer'), 
                    tissue = c('Colon','Colorectum'))
write.table(clu_ann, "clu_ann.txt", sep='\t', row.names=F, col.names=T, quote=F)

#Assignment of Cell Types
pbmc <- RenameIdents( pbmc, '0'= 'NA', '1'='NA','2'='NA', '3'= 'NA','4'= 'NA','5'= 'NA', '6'= 'NA','7'='NA','8'= 'NA','9'= 'NA','10'= 'NA','11'='Cancer-Associated Fibroblast','12'= 'NA','13'= 'NA','14'= 'NA','15'= 'Endothelial Cell','16'= 'NA','17'= 'Cancer-Associated Fibroblast, Stem Cell','18'= 'NA','19'= 'NA')
DimPlot(pbmc, reduction = "tsne", label=TRUE)
