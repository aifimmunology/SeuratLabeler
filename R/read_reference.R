#' Read in and validate input reference object file
#'
#' Currently supports .rds or .h5seurat input file and enforces read in
#' object as type 'Seurat'.
#'
#' @param ref_fp File path to reference object, an .rds or .h5seurat object
#' @return Seurat object contained in input file
#' @export 
read_reference <- function(ref_fp){
    valid_extensions <- c("rds", "h5seurat")
    
    extension <- tolower(gsub(".*[.]", "", ref_fp))
    if (!(extension %in% valid_extensions)) {
        stop(sprintf("Reference file extension [%s] not in valid_extensions [%s]",
                    extension,
                    paste(valid_extensions, collapse = ", ")))
    }
    
    if(extension == "rds"){
        ref_so <- readRDS(ref_fp)
        if(!inherits(ref_so, 'Seurat')){
            stop(sprintf("Input reference RDS file does not contain a Seurat object. Object in file '%s' is '%s'",
                 ref_fp,
                 class(ref_so)))
        }
    } else if (extension == "h5seurat") {
        ref_so <- SeuratDisk::LoadH5Seurat(ref_fp)
    } else {
        stop(sprintf("Need to define a read method for valid reference extension '%s'", extension))
    }
    
    return(ref_so)
}