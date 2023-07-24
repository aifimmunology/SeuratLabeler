#! /usr/bin/Rscript

library(SeuratLabeler)

fp_big <- system.file("testdata/pbmc_1k_v3.h5", package = "SeuratLabeler")
fp_small <- testthat::test_path('testdata','pbmc_1k_v3_80cells.h5')
dup_pattern = "__duptest__"

testthat::test_that(
    'Seurat formatting succeeds for h5 with >200 cells', 
    {
        testthat::expect_warning(  # ignore the Seurat non-unique rownames warning
            testthat::expect_true(
                {
                    so_query_big <- format_query_so(fp_big, allow_few_cells = FALSE);
                    inherits(so_query_big, "Seurat")
                }
            ), 
            regexp = 'Non-unique features [(]rownames[)] present in the input matrix, making unique'
        )
    }
)

testthat::test_that(
    'Seurat formatting errors for h5 with <200 cells if allow_few_cells = FALSE', 
    {
        # ignore the Seurat non-unique rownames warning
        testthat::expect_error(
            so_query_sm <- format_query_so(fp_small, allow_few_cells = FALSE),
            regexp = 'ERROR: < 200 Cells as Input. Too Few for Labeling'
        )
    }
)

testthat::test_that(
    'Seurat formatting succeeds for h5 with <200 cells if allow_few_cells = TRUE', 
    {
        # Output Seurat object
        testthat::expect_warning(  # ignore the Seurat non-unique rownames warning
            testthat::expect_true(
                {
                    set.seed(3)
                    so_query_sm <- format_query_so(fp_small, allow_few_cells = TRUE);
                    inherits(so_query_sm, "Seurat")
                }
            ), 
            regexp = 'Non-unique features [(]rownames[)] present in the input matrix, making unique'
        )
        
        # object passed, generate in memory
        set.seed(3)
        so_query_sm <- suppressWarnings(format_query_so(fp_small, allow_few_cells = TRUE, dup_pattern = dup_pattern))
         
        # Number of cells greater than 200
        testthat::expect_gte(ncol(so_query_sm), 200)
        
        # Duplicate regex used
        testthat::expect_true(any(grepl(dup_pattern, colnames(so_query_sm))))
    }
)
