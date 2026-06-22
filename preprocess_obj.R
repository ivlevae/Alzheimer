library(Seurat)
library(DoubletFinder)
library(SeuratDisk)

RDoublet <- function(tmp){
  sweep.res.list <- paramSweep_v3(tmp, PCs = 1:30, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pKopt <- as.numeric(as.character(bcmvn$pK[bcmvn$BCmetric == max(bcmvn$BCmetric)]))
  pKopt <- pKopt[order(pKopt, decreasing = TRUE) ]
  pKopt <- pKopt[1]
  homotypic.prop <- modelHomotypic(tmp$seurat_clusters)
  nExp_poi <- round(0.05*length(colnames(tmp)))  ## Assuming 5% doublet formation rate
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  tmp <- doubletFinder_v3(tmp, PCs = 1:30, pN = 0.25, pK = pKopt, nExp = nExp_poi, reuse.pANN = FALSE)
  tmp <- doubletFinder_v3(tmp, PCs = 1:30, pN = 0.25, pK = pKopt, nExp = nExp_poi.adj, reuse.pANN = paste("pANN_0.25",pKopt,nExp_poi, sep="_"))
  return (tmp)
}

s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

convertHumanGeneList <- function(x){
  require("biomaRt")
  human = useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", mirror = "www", host = "dec2021.archive.ensembl.org")
  mouse = useEnsembl("ensembl", dataset = "mmusculus_gene_ensembl", mirror = "www", host = "dec2021.archive.ensembl.org")
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
}

m.s.genes <- convertHumanGeneList(cc.genes.updated.2019$s.genes)
m.g2m.genes <- convertHumanGeneList(cc.genes.updated.2019$g2m.genes)

ProcessSeu <- function(Seurat){
  Seurat <- NormalizeData(Seurat)
  Seurat <- FindVariableFeatures(Seurat, selection.method = "vst", nfeatures = 3000)
  Seurat <- ScaleData(Seurat, verbose = T, vars.to.regress = c('percent.mt', "percent.rb","S.Score","G2M.Score"))
  
  Seurat <- RunPCA(Seurat, npcs = 100)
  Seurat <- FindNeighbors(Seurat, dims = 1:100)
  Seurat <- FindClusters(Seurat, resolution = 1)
  Seurat <- RunUMAP(Seurat, dims = 1:100)
  DimPlot(object = Seurat, reduction = "umap")
  return (Seurat)
}



file_path <- paste0("/home/bnvlab2/Downloads/","counts_exprMatrix.tsv.gz")
count_matrix <- fread(file_path)  

rownames(count_matrix1) <- count_matrix1$V1
count_matrix1$V1 <- NULL


brainObject <- CreateSeuratObject(counts = count_matrix1)




brainObject[["percent.rb"]] <- PercentageFeatureSet(brainObject, pattern = "^RPS|^RPL|^MRPS|^MRPL", assay = 'RNA') #change to uppercase for human
brainObject[["percent.mt"]] <- PercentageFeatureSet(brainObject, pattern = "^MT-") #change to uppercase for human
brainObject <- CellCycleScoring(brainObject, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE, nbin = 12) #remove 'm.' if operating with human data
VlnPlot(brainObject, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", 'percent.rb'), ncol = 4)
dim(brainObject)
brainObject <- subset(brainObject, subset = nCount_RNA > 300 & nFeature_RNA > 400)
dim(brainObject)
brainObject <- subset(brainObject, subset = nCount_RNA > 300 & nCount_RNA < 13000 & nFeature_RNA > 400 & nFeature_RNA < 3800 & percent.mt < 30 & percent.rb < 40)
dim(brainObject)



brainObject <- ProcessSeu(brainObject)
brainObject <- RDoublet(brainObject)


DimPlot(brainObject)
View(brainObject)


brainObject <- subset(brainObject, cells = colnames(brainObject )[which(brainObject [[]][12] == 'Singlet')])
brainObject <- subset(brainObject , cells = colnames(brainObject )[which(brainObject [[]][13] == 'Singlet')])

brainObject <- ProcessSeu(brainObject)

DimPlot(brainObject)

dim(brainObject)


filename <- 'GSM5618238_6'

brainObject$GSE <- 'GSE185553'
brainObject$GSM  <- 'GSM5618238_6' 

SaveH5Seurat(brainObject, paste0('/home/bnvlab2/Documents/Kate/Alzheimer/Hippocampus/Output/', filename, '.h5Seurat'), overwrite = TRUE)
saveRDS(brainObject, paste0('/home/bnvlab2/Documents/Kate/Alzheimer/Hippocampus/Output/', filename, '.rds'))

