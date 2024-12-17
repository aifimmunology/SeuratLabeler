library(SeuratLabeler)
library(SeuratDisk)
library(SeuratData)

in_h5 <- system.file("testdata/pbmc_1k_v3.h5", package = "SeuratLabeler")
in_ref_test <- 
h5_list <- H5weaver::h5dump(in_h5)

# Super small example data, <200 cells
set.seed(3)
bc <- sample(h5_list$matrix$barcodes, 80, replace = FALSE)
h5_list_sm <- H5weaver::subset_h5_list_by_observations(
    h5_list = h5_list, 
    match_values = bc, 
    match_target = 'barcodes'
)

H5weaver::write_h5_list(
    h5_list = h5_list_sm, 
    h5_file = "./testdata/pbmc_1k_v3_80cells.h5", 
    overwrite = TRUE
)

# Downsampled Hao 2021 reference
## Download and format data
download.file(
  url = 'https://atlas.fredhutch.org/data/nygc/multimodal/pbmc_multimodal.h5seurat',
  destfile = 'pbmc_multimodal.h5seurat'
)

hao2021_pbmc <- SeuratDisk::LoadH5Seurat('pbmc_multimodal.h5seurat')

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

# Downsample to 100 cells ber l2.5 label
Idents(hao2021_pbmc) <- "celltype.l2.5"
set.seed(3)
small_ref <- subset(x = hao2021_pbmc, downsample = 100)

# Save RDS
saveRDS(small_ref, file = "./testdata/hao2021_pbmc_downsampled.rds")

# Save h5 Seurat
SeuratDisk::SaveH5Seurat(small_ref, filename = './testdata/hao2021_pbmc_downsampled.h5seurat')

file.remove('pbmc_multimodal.h5seurat')

# Downsampled BMMC reference
## Download and format data
inst_data <- SeuratData::InstalledData()
if(!("bmcite" %in% inst_data$Dataset)) {
    SeuratData::InstallData("bmcite")
}
bm <- SeuratData::LoadData(ds = 'bmcite')
DefaultAssay(bm) <- 'RNA'

# Downsample to 100 cells ber l2 label, recompute spca
set.seed(3)
Idents(bm) <- "celltype.l2"
bm_small <- subset(x = bm, downsample = 100)

DefaultAssay(object = bm_small) <- 'RNA'
bm_small <- NormalizeData(object = bm_small)
bm_small <- FindVariableFeatures(object = bm_small)
bm_small <- ScaleData(object = bm_small)
bm_small <- RunPCA(object = bm_small)
DefaultAssay(object = bm_small) <- 'ADT'
VariableFeatures(object = bm_small) <- rownames(x = bm_small[["ADT"]])
bm_small <- NormalizeData(object = bm_small, normalization.method = 'CLR', margin = 2)
bm_small <- ScaleData(object = bm_small)
bm_small <- RunPCA(object = bm_small, reduction.name = 'apca')
bm_small <- FindMultiModalNeighbors(
          object = bm_small, reduction.list = list("pca", "apca"),
          dims.list = list(1:30, 1:18), modality.weight.name = "RNA.weight")
bm_small <- RunSPCA(bm_small, assay = 'RNA', graph = 'wsnn')
bm_small <- DietSeurat(object = bm_small, scale.data = FALSE, assays = c("RNA", "ADT"),
          dimreducs = "spca", graphs = c("wknn", "wsnn"))
DefaultAssay(object = bm_small) <- 'RNA'

# Save RDS
saveRDS(bm_small, file = './testdata/bmcite_downsampled.rds')

# Save h5 Seurat
SeuratDisk::SaveH5Seurat(bm_small, filename = './testdata/bmcite_downsampled.h5seurat')

# Copy of small BM data
inst_data <- SeuratData::InstalledData()
if(!("hcabm40k" %in% inst_data$Dataset)) {
    SeuratData::InstallData('hcabm40k')
}
SeuratData::LoadData(ds = 'hcabm40k')
DefaultAssay(hcabm40k) <- 'RNA'
hcabm40k_BM1 <- subset(hcabm40k, idents = 'MantonBM1')
saveRDS(hcabm40k_BM1, './testdata/bm_hcabm40k_downsampled.rds')