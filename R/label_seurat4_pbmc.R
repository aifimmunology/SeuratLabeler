.validate_future_cores <- function(ncores, verbose = TRUE){
  if (ncores > 1){
    if(ncores > (future::availableCores() - 1)){
      if(verbose){
        message(
          sprintf("%s available cores. Reducing future cores from %s to %s",
                  future::availableCores(), ncores,  future::availableCores()-1)
          )
      }
      ncores <- future::availableCores() - 1
    }
  }
  return(ncores)
}

.check_transfer_values <- function(reference_so, ref_data){
  imeta <- which(unlist(ref_data) %in% names(reference_so@meta.data))
  if(length(imeta) > 0){
    meta = unlist(ref_data)[imeta]
  } else {
    meta <- NULL
  }

  idata <- which(unlist(ref_data) %in% names(reference_so@assays))
  if(length(idata) > 0){
    dat = unlist(ref_data)[idata]
  } else {
    dat = NULL
  }

  missing <- setdiff(unlist(ref_data), c(meta, dat))
  if(length(missing) > 0){
    stop(sprintf("Reference data to transfer not found in reference object: %s",
                 paste(missing, collapse = ", ")))
  }

  return(list(meta = meta, data = dat))
}

#' Label Seurat with PBMC Reference
#'
#' @param query_so The query seurat object
#' @param ref_so The reference seurat object
#' @param ref_reduction Reduction value passed to Seurat::FindTransferAnchors(). Defualt 'spca'.
#' @param ref_data List of reference data passed to Seurat::TransferData(). Default celltype.l1,
#' @param dims Dimensionsfor FindTransferAnchors, defualt 1:50
#' @param verbose Default FALSE
#' @param ncores Number of cores for parallel processing via Seurat's future::plan implementation.
#' @param sct_method Method argument passed to Seurat::SCTransform, default 'glmGamPoi'
label_seurat4_pbmc <- function(
    query_so,
    ref_so,
    ref_reduction = "spca",
    ref_data = list(
      celltype.l1 = "celltype.l1",
      celltype.l2 = "celltype.l2",
      celltype.l2.5 = "celltype.l2.5",
      celltype.l3 = "celltype.l3",
      predicted_ADT = "ADT"
    ),
    dims = 1:50,
    verbose = FALSE,
    ncores = 1,
    sct_method = "glmGamPoi"
  )
  {

  if(!(ref_reduction %in% names(ref_so@reductions))) {
    stop(sprintf("Reference dataset reduction '%s' not found in reference data object",
                 ref_reduction))
  }

  # validate and group the ref_data types
  ref_data_types <- .check_transfer_values(reference_so = ref_so, ref_data = ref_data)

  # Set Cores
  if (verbose){ message("Setting multiple cores for SCTransform") }
  ncores <- .validate_future_cores(ncores = ncores, verbose = verbose)
  if (ncores > 1){
    future::plan(strategy = 'multicore', workers = ncores)
  }

  if(verbose){message("Finding anchors")}
  anchors <- Seurat::FindTransferAnchors(
    reference = ref_so,
    query = query_so,
    reference.assay = 'SCT',
    normalization.method = 'SCT',
    reference.reduction = ref_reduction,
    recompute.residuals = TRUE,  # function default
    dims = dims,
    verbose = verbose
  )

  query_so <- Seurat::TransferData(
    anchorset = anchors,
    query = query_so,
    reference = ref_so,
    refdata = ref_data,
    weight.reduction = 'pcaproject',
    verbose = verbose,
    l2.norm = FALSE,
    dims = NULL,
    k.weight = 50,
    sd.weight = 1,
    eps = 0,
    n.trees = 50,
    slot = "data",
    prediction.assay = FALSE,
    store.weights = TRUE
  )

  return(query_so)
}
