# Make RDS labeling reference file of Hao et al 2021 PBMC CITE-seq dataset

# Set Up
library(Seurat)
library(SeuratDisk)

# Get data
download.file(
  url = 'https://atlas.fredhutch.org/data/nygc/multimodal/pbmc_multimodal.h5seurat',
  destfile = 'pbmc_multimodal.h5seurat'
)

# Load data
hao2021_pbmc <- SeuratDisk::LoadH5Seurat('pbmc_multimodal.h5seurat')

# Add L2.5 label
make_l2.5 <- function(l2, l3){
  l2.5 <- l2
  l2.5[l3 == "Treg Naive"] <- "Treg Naive"
  l2.5[l3 == "Treg Memory"] <- "Treg Memory"
  l2.5[l3 %in% c("CD8 TEM_4", "CD8 TEM_5")] <- "CD8 TEMRA"
  return(l2.5)
}

hao2021_pbmc$celltype.l2.5 <- make_l2.5(
  l2 = hao2021_pbmc$celltype.l2,
  l3 = hao2021_pbmc$celltype.l3
)

# Save Rda
saveRDS(hao2021_pbmc, file = "hao2021_pbmc.rds")

# Remove original download
file.remove("pbmc_multimodal.h5seurat")
