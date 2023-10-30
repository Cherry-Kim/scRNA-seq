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

        pbmc <- readRDS("mouse.rds")
        ## An object of class Seurat
        ## 23761 features across 121767 samples within 1 assay
        ## Active assay: RNA (23761 features, 0 variable features)

        ## Set identity classes to an existing column in meta data
        Idents(object = pbmc) <- pbmc@meta.data$annotation
        a <- subset(pbmc, idents = "Liver")     
        head(a@meta.data)

