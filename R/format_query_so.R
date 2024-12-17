#' Format input H5 as a Seurat object 
#'
#' Formats the .H5 file as a Seurat object with option to
#' bootstrap if there are fewer than 200 cells. 
#'
#' For bootstrapping, the entire dataset will be duplicated 
#' enough times to obtain 200 cells prior to labeling. Duplicates
#' will be removed downstream.
#'
#' @param in_h5 Filepath to input h5 file
#' @param allow_few_cells Logical value, default FALSE. If TRUE allows proceeding when
#' input dataset has fewer thatn 200 cells but bootstrapping up to >=200 cells. If FALSE,
#' errors out if too few cells
#' @param dup_pattern String value prefix for bootstrapped cellnames when allow_few_cells
#' behavior is triggered. Default '__dup__'
format_query_so <- function(in_h5, allow_few_cells = FALSE, dup_pattern = "__dup__"){
    # Read in counts
    query_data <- H5weaver::read_h5_dgCMatrix(
        h5_file = in_h5,
        feature_names = "name"
    )
    rownames(query_data) <- as.vector(rhdf5::h5read(in_h5, "/matrix/features/name"))
    
    # Handle low cell count option
    if(allow_few_cells & ncol(query_data) < 200) {
      stm("WARNING: < 200 Cells as Input. Attempting compensation")
      stm("WARNING: Labels may be inaccurate")

      ndup <- ceiling(200 / ncol(query_data))
        
      query_data_orig <- query_data
      for(i in 1:(ndup-1)) {
        dup_data <- query_data_orig
        colnames(dup_data) <- paste0(dup_pattern, i, "_", colnames(dup_data))
        query_data <- cbind(query_data, dup_data)
      }
      rm(query_data_orig)
    } else if (ncol(query_data) < 200) {
      stm("ERROR: < 200 Cells as Input. Too Few for Labeling")
      stop("ERROR: < 200 Cells as Input. Too Few for Labeling")
    } 
    
    query_so <- Seurat::CreateSeuratObject(counts = query_data)
    
    return(query_so)
}