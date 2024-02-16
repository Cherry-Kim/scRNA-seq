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
