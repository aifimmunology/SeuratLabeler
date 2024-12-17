library(SeuratLabeler)

# Test input files
fp_ref_so <- testthat::test_path('testdata','hao2021_pbmc_downsampled.rds')
fp_big <- system.file("testdata/pbmc_1k_v3.h5", package = "SeuratLabeler")
fp_small <- testthat::test_path('testdata','pbmc_1k_v3_80cells.h5')

ref_pbmc <- read_reference(fp_ref_so)
so_query_big <- suppressWarnings(format_query_so(fp_big, allow_few_cells = FALSE))
so_query_sm <- suppressWarnings(format_query_so(fp_small, allow_few_cells = TRUE))

exp_warn_pat_list <- list(
    messages = character(), 
    warnings = c('Keys should be one or more alphanumeric characters',
                 'Feature names cannot have underscores')
)

testthat::test_that(
    'label_seurat4_pbmc succeeds without error', 
    {
        testthat::expect_true({
                quiet_so_list <- purrr::quietly(SeuratLabeler::label_seurat4_pbmc)(
                    query_so = so_query_big, 
                    ref_so = ref_pbmc
                );
                res <- filter_quietly(quiet_so_list, patterns_suppress = exp_warn_pat_list);
                inherits(res, "Seurat")
            }
        )
    }
)

testthat::test_that(
    'label_seurat4_pbmc transfers all labels without error', 
    {
        quiet_so_list <- purrr::quietly(SeuratLabeler::label_seurat4_pbmc)(query_so = so_query_big, ref_so = ref_pbmc)
        res <- filter_quietly(quiet_so_list, patterns_suppress = exp_warn_pat_list)
        expected_cols = c('predicted.celltype.l1.score','predicted.celltype.l1',
                          'predicted.celltype.l2.score','predicted.celltype.l2',
                          'predicted.celltype.l2.5.score','predicted.celltype.l2.5',
                          'predicted.celltype.l3.score','predicted.celltype.l3')
        testthat::expect_true(all(expected_cols %in% names(res@meta.data)))
        
    }
)

testthat::test_that(
    'label_seurat4_pbmc transfers predicted ADT without error', 
    {
        quiet_so_list <- purrr::quietly(SeuratLabeler::label_seurat4_pbmc)(query_so = so_query_big, ref_so = ref_pbmc)
        res <- filter_quietly(quiet_so_list, patterns_suppress = exp_warn_pat_list)
        expected_assays = c('RNA','predicted_ADT')
        testthat::expect_true(all(expected_assays %in% Assays(res)))
        testthat::expect_true(nrow(res[['predicted_ADT']]@data) == nrow(ref_pbmc[['ADT']]@data))
    }
)

testthat::test_that(
    'label_seurat4_pbmc correctly errors on bad ref_reduction value', 
    {
        testthat::expect_error(
            {
                quiet_so_list <- SeuratLabeler::label_seurat4_pbmc(
                    query_so = so_query_big, 
                    ref_so = ref_pbmc, 
                    ref_reduction = 'nonexistent_reduction')
            },
            regexp = "Reference dataset reduction '.*' not found in reference data object"
        )
    }
)

testthat::test_that(
    'label_seurat4_pbmc correctly errors on bad transfer data label name value', 
    {
        ref_data_list <- list(celltype.l4='celltype.l4')
        testthat::expect_error(
            {
                quiet_so_list <- SeuratLabeler::label_seurat4_pbmc(
                    query_so = so_query_big, 
                    ref_so = ref_pbmc, 
                    ref_data = ref_data_list
                )
            },
            regexp = "Reference data to transfer not found in reference object: "
        )
    }
)

testthat::test_that(
    'label_seurat4_pbmc correctly errors on bad transfer data assay name value', 
    {
        ref_data_list <- list(predicted_CITE='CITE')
        testthat::expect_error(
            {
                quiet_so_list <- SeuratLabeler::label_seurat4_pbmc(
                    query_so = so_query_big, 
                    ref_so = ref_pbmc, 
                    ref_data = ref_data_list
                )
            },
            regexp = "Reference data to transfer not found in reference object: "
        )
    }
)

testthat::test_that(
    'label_seurat4_pbmc runs multicore without error', 
    {
        testthat::expect_true({
            quiet_so_list <- purrr::quietly(SeuratLabeler::label_seurat4_pbmc)(
                query_so = so_query_big, 
                ref_so = ref_pbmc, 
                ncores = 10
            );
            res <- filter_quietly(quiet_so_list, patterns_suppress = exp_warn_pat_list);
            inherits(res, "Seurat")            
        })
    }
)