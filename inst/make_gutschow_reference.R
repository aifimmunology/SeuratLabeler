library(Seurat)
library(Matrix)

labeling_dir <- "G:/Shared drives/Imm - Computational Biology/Cell_Labelling_Project/"
pbmc <- readRDS(file.path(labeling_dir,"Processed10X/pbmc5surf_final.rds"))
labels <- read.csv(file.path(labeling_dir,"Antibody Labels based on 10X Cite Seq/20200106_CITE_5K_labels.csv"))

pbmc@meta.data$cell_type <- labels$Curation[match(colnames(pbmc), labels$Barcodes)]

pbmc <- NormalizeData(pbmc)

pbmc <- FindVariableFeatures(object = pbmc, nfeatures = 5000)

variable_pbmc_5k_cite <- pbmc

variable_pbmc_5k_cite@assays$RNA@counts <- variable_pbmc_5k_cite@assays$RNA@counts[VariableFeatures(pbmc),]
variable_pbmc_5k_cite@assays$RNA@data <- variable_pbmc_5k_cite@assays$RNA@data[VariableFeatures(pbmc),]
variable_pbmc_5k_cite@assays$RNA@scale.data <- variable_pbmc_5k_cite@assays$RNA@scale.data[0,]

save(variable_pbmc_5k_cite, file = "data/variable_pbmc_5k_cite.RData")
