library(Seurat)
library(Matrix)

download.file("http://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_10k_v3/pbmc_10k_v3_filtered_feature_bc_matrix.h5",
              "inst/pbmc_10k_v3_filtered_feature_bc_matrix.h5",
              mode = "wb")

pbmc <- Read10X_h5("inst/pbmc_10k_v3_filtered_feature_bc_matrix.h5")
pbmc <- CreateSeuratObject(counts = pbmc,
                           project = "pbmc_10k_v3",
                           min.cells = 3,
                           min.features = 200)

pbmc <- NormalizeData(pbmc)

pbmc <- FindVariableFeatures(object = pbmc, nfeatures = 5000)

variable_pbmc_10k_v3 <- pbmc

variable_pbmc_10k_v3@assays$RNA@counts <- variable_pbmc_10k_v3@assays$RNA@counts[VariableFeatures(pbmc),]
variable_pbmc_10k_v3@assays$RNA@data <- variable_pbmc_10k_v3@assays$RNA@data[VariableFeatures(pbmc),]

save(variable_pbmc_10k_v3, file = "data/variable_pbmc_10k_v3.RData")

file.remove("inst/pbmc_10k_v3_filtered_feature_bc_matrix.h5")
