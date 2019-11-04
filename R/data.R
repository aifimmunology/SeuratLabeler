#' Seurat reference 10X PBMC dataset
#'
#' This dataset was retrieved from the Seurat website from this URL: https://www.dropbox.com/s/3f3p5nxrn5b3y4y/pbmc_10k_v3.rds?dl=1
#'
#' It was loaded and 5,000 variable features were selected.
#'
#' The data was then filtered to retain only these variable features for portability.
#'
#' The script used to generate this file is available within this package in inst/make_reference_data.R
#'
#' @format A Seurat object with data from 5,000 genes and 11,537 cells
#' @export
#'
"variable_pbmc_10k_v3"
