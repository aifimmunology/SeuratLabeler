#! /usr/bin/Rscript

library(SeuratLabeler)

# Test input files
fp_ref_so <- testthat::test_path('testdata','hao2021_pbmc_downsampled.rds')
fp_ref_h5 <- testthat::test_path('testdata','hao2021_pbmc_downsampled.h5seurat')

# Make example bad inputs
my_td <- tempdir()
df <- data.frame(x=1:10)
fp_bad_data <- file.path(my_td, "bad_data.RDS")
fp_bad_data_csv <- file.path(my_td, "bad_data.csv")
saveRDS(df, fp_bad_data)
write.csv(df, fp_bad_data_csv)

testthat::test_that(
    'read_reference succeeds for rds containing seurat object', 
    {
        testthat::expect_true(
            {
                so_ref <- read_reference(fp_ref_so);
                inherits(so_ref, "Seurat")
            }
        )
    }
)

testthat::test_that(
    'read_reference errors for rds containing non-seurat object', 
    {
        testthat::expect_error(
            so_ref <- read_reference(fp_bad_data),
            regexp = "Input reference RDS file does not contain a Seurat object",
        )
    }
)

fp_ref_h5 <- testthat::test_path('testdata','hao2021_pbmc_downsampled.rds')
testthat::test_that(
    'read_reference succeeds for h5seurat reference file type', 
    {
        testthat::expect_true(
            {
                so_ref <- read_reference(fp_ref_h5);
                inherits(so_ref, "Seurat")
            }
        )
    }
)

testthat::test_that(
    'read_reference errors for non-supported file type', 
    {
        testthat::expect_error(
            so_ref <- read_reference(fp_bad_data),
            regexp = "Input reference RDS file does not contain a Seurat object",
        )
    }
)

file.remove(fp_bad_data)
file.remove(fp_bad_data_csv)