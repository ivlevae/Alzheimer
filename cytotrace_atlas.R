library(CytoTRACE2)
library(Seurat)
library(SeuratDisk)
library(Matrix)

brain_data <- LoadH5Seurat('/home/bnvlab2/anndata2.h5seurat',    misc = FALSE, tools =FALSE, meta.data = FALSE)
metadata <- read.csv("/home/bnvlab2/adata_metadata.csv", row.names = 1)

brain_data <- AddMetaData(brain_data, metadata)

brain_data$annot_scANVI_leiden_30_3_256NB_5_res2_V2 <- as.factor(brain_data$annot_scANVI_leiden_30_3_256NB_5_res2_V2)
DimPlot(brain_data, reduction = "umap", group.by = 'annot_scANVI_leiden_30_3_256NB_5_res2_V2')

### filter low expressed genes to speed up calculations
min_cells <- 3

expr_matrix <- GetAssayData(brain_data, slot = "counts")

genes_expressed <- rowSums(expr_matrix > 0) >= min_cells
genes_to_keep <- names(genes_expressed[genes_expressed])
seurat_filtered <- subset(brain_data, features = genes_to_keep)


cytotrace2_result_fullmodel <- cytotrace2(seurat_filtered,
                                          species = "human",
                                          is_seurat = TRUE,
                                          slot_type = "counts",
                                          #full_model = TRUE,
                                          batch_size = 10000, #2500
                                          smooth_batch_size = 1000,
                                          parallelize_models = TRUE,
                                          parallelize_smoothing = TRUE,
                                          ncores = 4, #ncores = 4
                                          #max_pcs = 200,
                                          seed = 14)


gc()

annotation <- data.frame(phenotype = seurat_filtered@meta.data$annot_scANVI_leiden_30_3_256NB_5_res2_V2) %>% 
  set_rownames(., colnames(seurat_filtered))

View(annotation)

plots5 <- plotData(cytotrace2_result = cytotrace2_result_fullmodel,
                   annotation = annotation,
                   is_seurat = TRUE)



plots5[["CytoTRACE2_UMAP"]][[1]][["data"]]["UMAP_1"] <- seurat_filtered@reductions[["umap"]]@cell.embeddings[,1]
plots5[["CytoTRACE2_UMAP"]][[1]][["data"]]["UMAP_2"]<- seurat_filtered@reductions[["umap"]]@cell.embeddings[,2]

plots5[["CytoTRACE2_Potency_UMAP"]][[1]][["data"]]["UMAP_1"] <- seurat_filtered@reductions[["umap"]]@cell.embeddings[,1]
plots5[["CytoTRACE2_Potency_UMAP"]][[1]][["data"]]["UMAP_2"]<- seurat_filtered@reductions[["umap"]]@cell.embeddings[,2]


annotation
plots5
cytotrace2_result_fullmodel$CytoTRACE2_Score_factor <- as.factor(cytotrace2_result_fullmodel$CytoTRACE2_Score)
View(cytotrace2_result_fullmodel)

plots5$CytoTRACE2_Boxplot_byPheno  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


meta_cyto <- cytotrace2_result_fullmodel@meta.data

write.csv(meta_cyto, '/home/bnvlab2/Documents/Kate/Alzheimer/cytotrace2_meta.csv')

