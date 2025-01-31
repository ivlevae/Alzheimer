### SeuObjProcessedList was a dictionary 
### Active assay RNA 
wholebrain_list <- unname(unlist(wholebrain_list))
print(seurat_objects)
print(length(wholebrain_list))  
print(names(wholebrain_list))

gc()
ProcessInt <- function(data.integrated){
  data.integrated <- ScaleData(data.integrated, verbose = T, vars.to.regress = c('percent.mt',"percent.rb","S.Score","G2M.Score"))
  data.integrated <- RunPCA(data.integrated, npcs = 30, verbose = T)
  data.integrated <- FindNeighbors(data.integrated, dims = 1:30)
  data.integrated <- FindClusters(data.integrated, resolution = 1)
  data.integrated <- RunUMAP(data.integrated, reduction = "pca", dims = 1:30)
}

features <- SelectIntegrationFeatures(object.list = wholebrain_list, nfeatures = 2000)
gc()
data.anchors <- FindIntegrationAnchors(object.list = wholebrain_list, anchor.features = features)
data.anchors

data.anchorsrpca <- FindIntegrationAnchors(object.list = wholebrain_list, anchor.features = features, reduction = 'rpca') 
data.anchorsrpca
gc()
wholebrain <- IntegrateData(anchorset = data.anchorsrpca )

gc()

print(DefaultAssay(wholebrain))

wholebrain <- ProcessInt(wholebrain)

ProcessInt <- function(data.integrated){
  data.integrated <- RunPCA(data.integrated, npcs = 150, verbose = T)
  data.integrated <- FindNeighbors(data.integrated, dims = 1:150)
  data.integrated <- FindClusters(data.integrated, resolution = 1)
  data.integrated <- RunUMAP(data.integrated, reduction = "pca", dims = 1:150)
}


options(future.globals.maxSize = 100000 * 1024^2)
wholebrain150 <- ProcessInt(wholebrain)



DimPlot(wholebrain150, label = T, repel = T, label.box = T, raster = T)
DimPlot(wholebrain,raster = T)

DimPlot(wholebrain, raster = T, group.by = 'GSM')

wholebrain$donor_id

dim(wholebrain)

wholebrain30@active.ident <- wholebrain30$seurat_clusters
wholebrain@active.ident <- wholebrain$seurat_clusters


markerswholebrain <- FindAllMarkers(wholebrain30, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)

gc()
markerswholebrain150 <- FindAllMarkers(wholebrain, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)


sorted_data <- markerswholebrain %>%
  filter(p_val_adj < 0.05) %>%
  arrange(cluster, desc(avg_log2FC))

write.csv(sorted_data, '/home/bnvlab2/Documents/Kate/Alzheimer/Hippocampus/Integration_results/markerswholebrain.csv')

sorted_data150 <- markerswholebrain150 %>%
  filter(p_val_adj < 0.05) %>%
  arrange(cluster, desc(avg_log2FC))

write.csv(markerswholebrain150, '/home/bnvlab2/Documents/Kate/Alzheimer/Hippocampus/Integration_results/markerswholebrain150.csv')
DefaultAssay(wholebrain) <- 'RNA'
DefaultAssay(wholebrain150) <- 'RNA'

library(SeuratDisk)
SaveH5Seurat(wholebrain150, '/home/bnvlab2/Documents/Kate/Alzheimer/Hippocampus/Integration_results/wholebrain150.h5Seurat', overwrite = TRUE)
SaveH5Seurat(wholebrain, '/home/bnvlab2/Documents/Kate/Alzheimer/Hippocampus/Integration_results/wholebrain.h5Seurat', overwrite = TRUE)

Convert('/home/bnvlab2/Documents/Kate/Alzheimer/Hippocampus/Integration_results/wholebrain.h5Seurat', dest = 'h5ad', assay = 'RNA')


saveRDS(wholebrain150, '/home/bnvlab2/Documents/Kate/Alzheimer/Hippocampus/Integration_results/wholebrain150.rds')
saveRDS(wholebrain, '/home/bnvlab2/Documents/Kate/Alzheimer/Hippocampus/Integration_results/wholebrain.rds')

dim(wholebrain150)
dim(prenatal150)
DefaultAssay(prenatal150) <- 'RNA'

forebrain_biomarkers <- list(
  Progenitor_Cells = c("PAX6", "SOX2", "HES5", "VIM"),
  Neurogenesis = c("DCX", "NEUROD1", "ASCL1"),
  Excitatory_Neurons = c("TBR1", "SATB2", "CTIP2"),
  Inhibitory_Neurons = c("GAD1", "GAD2", "DLX1", "DLX2", "ARX"),
  Glial_Cells = c("OLIG2", "GFAP", "S100B"),
  Cerebral_Cortex = c("EMX1", "EMX2", "FEZF2", "RELN"),
  Hippocampus = c("PROX1", "ZBTB20", "CALB1"),
  Basal_Ganglia = c("FOXP1", "GSX2", "NKX2-1")
)

# Print the list
print(forebrain_biomarkers)

DefaultAssay(wholebrain) <- 'RNA'
p <- DotPlot(wholebrain,
             features = forebrain_biomarkers,
             cols = c("lightgrey", "blue"))

p + theme(axis.text.x = element_text(angle = 90) )
