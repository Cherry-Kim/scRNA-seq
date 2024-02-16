Step1_input <- function(){ 
        library(sceasy)
        library(reticulate)
        library(dplyr)
        library(Seurat)
        library(monocle3)
        library(patchwork)
        library(ggplot2)
        library(harmony)
        library(ggdag)
        library(ggplot2)
        library(Matrix)

        # 1. counts.data
        counts.data <- read.delim('E9.5_data_log_seurat.txt', header=T, row.names=1)
        colnames(counts.data) <- sub("^X", "", colnames(counts.data))
        head(counts.data[1:5,1:3])
        #              23365 23366 23367 23368 23369
        #0610005C13Rik      0      0      0      0      0
      
        pheno.data <- read.csv('E9.5.singlet_type.XY.prop.csv', header=T, row.names=1)
        head(pheno.data[1:5,1:5])
        #      spot_class first_type second_type                               celltype1
        #23365    singlet         18          18 Lateral plate 


        annotations_to_filter <- c("Facial mesenchyme")
        rows <- pheno.data[pheno.data$celltype1 %in% annotations_to_filter, ]   
        row_names_to_match <- rownames(rows)
        a <- as.data.frame(counts.data)
        matching_rows <- colnames(a) %in% row_names_to_match
        subset_counts.data <- counts.data[,matching_rows, drop = FALSE]
        dim(subset_counts.data) 

        mat = Matrix(as.matrix(subset_counts.data),sparse=TRUE)
        #4 x 4 sparse Matrix of class "dgCMatrix"
        #      23382 23399 23583 23618
        #0610005C13Rik     .     .     .     .
        writeMM(obj = mat, file="matrix.mtx")
        write(x = rownames(subset_counts.data), file = "genes.tsv")
        write(x = colnames(subset_counts.data), file = "barcodes.tsv"
}

STEP2_seurat <- function(){
        pbmc.data <- Read10X(data.dir = "/BIO3/STEP3_Seurat/", gene.column=1)
        pbmc <- CreateSeuratObject(counts = pbmc.data)

        pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = '^MT-')
        # FeatureScatter can be used for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
        plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")+ theme(legend.position = "none")
        plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+ theme(legend.position = "none")
        plot1 + plot2
        ## seu.obj.filtered <- subset(seu.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & mitopercent < 5)

        #Step2. Normalizing the data
        pbmc <- NormalizeData(pbmc)

        # Identification of highly variable features (feature selection)
        pbmc <- FindVariableFeatures(pbmc)

        #Scaling the data
        all.genes <- rownames(pbmc)
        pbmc <- ScaleData(pbmc, features = all.genes)

        #Perform linear dimensional reduction
        pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
        ElbowPlot(pbmc)

        #Cluster the cells
        pbmc <- FindNeighbors(pbmc)
        pbmc <- FindClusters(pbmc)      #resolution = 0.5
        head(Idents(pbmc), 5)

        #Run non-linear dimensional reduction (UMAP/tSNE)
        pbmc <- RunUMAP(pbmc, dims = 1:15) #pbmc <- RunTSNE(pbmc, dims = 1:8)   
        DimPlot(pbmc, reduction = 'umap')

        #Finding differentially expressed features (cluster biomarkers)
        pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)
        pbmc.markers %>%
                group_by(cluster) %>%
                dplyr::filter(avg_log2FC > 1)

        write.csv(pbmc.markers, "markers.csv", quote=F)
}

