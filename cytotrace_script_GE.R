library(CytoTRACE2)
library(Seurat)
library(SeuratDisk)


DefaultAssay(ge_data_filtered) <- 'RNA'
#RNA >> FindVariableFeatures >> ScaleData (scaled data for 2k x n cells)

ge_data_filtered@active.ident <- as.factor(ge_data_filtered$seurat_clusters)
DimPlot(ge_data_filtered, reduction = "umap")  

cytotrace2_result_fullmodel <- cytotrace2(ge_data_filtered,
                                          species = "human",
                                          is_seurat = TRUE,
                                          slot_type = "counts",
                                          full_model = TRUE,
                                          batch_size = 10000, #2500
                                          smooth_batch_size = 1000,
                                          parallelize_models = TRUE,
                                          parallelize_smoothing = TRUE,
                                          ncores = 16, #ncores = 4
                                          max_pcs = 200,
                                          seed = 14)

gc()
#write.csv(cytotrace2_result_fullmodel, '/home/baranov_lab/Documents/Kate/Cornea/Cells_Subset/cytotrace2_result_fullmodel.csv')
SaveH5Seurat(cytotrace2_result_fullmodel, 'cytotrace2_ge_data_filtered.h5Seurat', overwrite = TRUE)

cytotrace2_ge_data <- LoadH5Seurat("/home/bnvlab2/Documents/Kate/Alzheimer/cytotrace2_ge_data_filtered.h5Seurat")

table(cytotrace2_ge_data$id)

#cell_cluster
#seurat_clusters
annotation <- data.frame(phenotype = ge_data_filtered@meta.data$seurat_clusters) %>% 
              set_rownames(., colnames(ge_data_filtered))

View(annotation)
# plotting
plots5 <- plotData(cytotrace2_result = cytotrace2_result_fullmodel,
                   annotation = annotation,
                   is_seurat = TRUE)

plots51 <- plotData(cytotrace2_result = cytotrace2_result_fullmodel,
                   annotation = annotation,
                   is_seurat = TRUE)

plots5[["CytoTRACE2_UMAP"]][[1]][["data"]]["UMAP_1"] <- ge_data_filtered@reductions[["umap"]]@cell.embeddings[,1]
plots5[["CytoTRACE2_UMAP"]][[1]][["data"]]["UMAP_2"]<- ge_data_filtered@reductions[["umap"]]@cell.embeddings[,2]

plots5[["CytoTRACE2_Potency_UMAP"]][[1]][["data"]]["UMAP_1"] <- ge_data_filtered@reductions[["umap"]]@cell.embeddings[,1]
plots5[["CytoTRACE2_Potency_UMAP"]][[1]][["data"]]["UMAP_2"]<- ge_data_filtered@reductions[["umap"]]@cell.embeddings[,2]

annotation
plots5
cytotrace2_result_fullmodel$CytoTRACE2_Score_factor <- as.factor(cytotrace2_result_fullmodel$CytoTRACE2_Score)
View(cytotrace2_result_fullmodel)

cytotrace2_result_fullmodel$id <- ge_data_filtered$id
multipotent_subset <- subset(cytotrace2_result_fullmodel, 
                                             seurat_clusters %in% c('23','13','4','11', '21', '7' ) &
                                             CytoTRACE2_Potency == 'Multipotent')

dim(multipotent_subset)

DimPlot(multipotent_subset)

dim(cytotrace2_result_fullmodel)

p <- DotPlot(
  seurat_obj,
  features = scRNA_seq_progenitor_markers,
  assay = NULL
)
p + theme(axis.text.x = element_text(angle = 90) )


progenitor_markers <- c("NES", "SOX2", "PAX6", "NEUROG2", "EOMES", "DCX", "ASCL1", "FOXG1", "EMX2", "LHX2")
hippocampal_progenitor_markers <- c(
  # Neural Stem Cells (NSCs)
  "SOX2", "NES", "VIM", "FABP7", "SLC1A3",
  
  # Intermediate Neural Progenitors (INPs)
  "EOMES", "PAX6", "NEUROG2",
  
  # Differentiating Neurons
  "DCX", "NEUROD1", "PROX1", "ASCL1",
  
  # Hippocampal Niche-Specific Markers
  "GFAP", "SOX9",
  
  # Proliferation Markers
  "MKI67", "PCNA",
  
  # Dentate Gyrus-Specific Markers
  "TBR1", "ZBTB20", "CR"
)

# To print the list
print(hippocampal_progenitor_markers)

scRNA_seq_progenitor_markers <- c(
  # Neural Stem Cells (NSCs)
  "MKI67", "TOP2A",  
  "ASCL1", "DLX1", "DLX2", "DCX", "EOMES", "PAX6",
  "SOX2", "NES", "FABP7", 
  
  # Radial Glial Cells (RGCs)
  "HES1", "HES5", "GLI3", "SLC1A3", "PROM1",

  # Early Differentiating Neurons
 
  
  # Region-Specific Progenitors
  "FOXG1", "EMX2", "OTX2"

  
)



DefaultAssay(multipotent_subset) <- 'integrated'
multipotent_subset <- ScaleData(multipotent_subset, verbose = T, vars.to.regress = c('percent.mt',"percent.rb","S.Score","G2M.Score"))

SaveH5Seurat(multipotent_subset, '/home/bnvlab2/Documents/Kate/Alzheimer/multipotent_subset.h5Seurat', overwrite = T)


ependymal_cell_markers <- c(
  # General Markers
  "FOXJ1",  # Transcription factor critical for motile cilia formation
  "S100B",  # Calcium-binding protein expressed in ependymal and astrocytic cells
  "VIM",    # Vimentin, expressed in ependymal progenitors and immature ependymal cells
  
  # Cilia and CSF Flow-Associated Markers
  "DNAH5",  # Dynein axonemal heavy chain, involved in cilia motility
  "TPPP",   # Tubulin polymerization-promoting protein, cilia-associated
  "TUBB4A", # Tubulin beta 4A, a component of cilia
  
  # Junction and Polarity Markers
  "CDH1",   # E-cadherin, involved in apical junctions
  "ZO-1",   # Tight junction protein (also in vascular cells)
  
  # Other Associated Markers
  "AQP4",   # Aquaporin-4, associated with fluid regulation
  "SLC4A4"  # Bicarbonate transporter involved in CSF regulation
)

print(ependymal_cell_markers)

print(vascular_cell_markers)


DefaultAssay(cytotrace2_result_fullmodel) <- 'RNA'
p <- DotPlot(cytotrace2_result_fullmodel,
  features = ependymal_cell_markers,
  assay = NULL)

p + theme(axis.text.x = element_text(angle = 90) )

cluster_potency_table <- table(cytotrace2_result_fullmodel$seurat_clusters, cytotrace2_result_fullmodel$CytoTRACE2_Potency)
cluster_potency_df <- as.data.frame.matrix(cluster_potency_table)
View(cluster_potency_percentage)
cluster_potency_percentage <- cluster_potency_df / rowSums(cluster_potency_df) * 100
plots5$CytoTRACE2_UMAP
plotData()
plots5$CytoTRACE2_UMAP
DimPlot(cytotrace2_result_fullmodel, group.by = 'CytoTRACE2_Potency', split.by = 'CytoTRACE2_Potency', raster = T)

cytotrace2_result_fullmodel@meta.data[["CytoTRACE2_Score"]]




