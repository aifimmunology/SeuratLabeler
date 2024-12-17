library(Seurat)
library(SeuratData)

InstallData("bmcite")
bm <- LoadData(ds = 'bmcite')

DefaultAssay(bm) <- 'RNA'
saveRDS(bm, file = 'bmcite.rds')
