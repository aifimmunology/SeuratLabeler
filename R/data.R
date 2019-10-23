#' Seurat reference 10X PBMC dataset
#'
#' This dataset was retrieved from 10x Genomics here: http://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_10k_v3/pbmc_10k_v3_filtered_feature_bc_matrix.h5
#'
#' It was loaded as a Seurat object, normalized with default settings, and 5,000 variable features were selected.
#'
#' The data was then filtered to retain only these variable features for portability.
#'
#' The script used to generate this file is available within this package in inst/make_reference_data.R
#'
#' @format A Seurat object with data from 5,000 genes and 11,537 cells
#' @export
#'
"variable_pbmc_10k_v3"
