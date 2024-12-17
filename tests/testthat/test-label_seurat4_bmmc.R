# test-label_seurat4_bmmc

library(SeuratLabeler)
library(SeuratData)

# Test input files
fp_ref_so <- testthat::test_path('testdata','bmcite_downsampled.rds') 
fp_qry <- testthat::test_path('testdata','bm_hcabm40k_downsampled.rds') 

ref_bm <- read_reference(fp_ref_so)
so_query <- readRDS(fp_qry)

exp_warn_pat_list <- list(
    messages = character(), 
    warnings = c('Keys should be one or more alphanumeric characters',
                 'Feature names cannot have underscores')
)

testthat::test_that(
    'label_seurat4_bmmc succeeds without error', 
    {
        testthat::expect_true({
                quiet_so_list <- purrr::quietly(SeuratLabeler::label_seurat4_bmmc)(
                    query_so = so_query, 
                    ref_so = ref_bm
                );
                res <- filter_quietly(quiet_so_list, patterns_suppress = exp_warn_pat_list);
                inherits(res, "Seurat")
            }
        )
    }
)

testthat::test_that(
    'label_seurat4_bmmc transfers all labels without error', 
    {
        quiet_so_list <- purrr::quietly(SeuratLabeler::label_seurat4_bmmc)(
            query_so = so_query, 
            ref_so = ref_bm
        )
        res <- filter_quietly(quiet_so_list, patterns_suppress = exp_warn_pat_list)
        expected_cols = c('predicted.celltype.l1.score','predicted.celltype.l1',
                          'predicted.celltype.l2.score','predicted.celltype.l2')
        testthat::expect_true(all(expected_cols %in% names(res@meta.data)))
        
    }
)

testthat::test_that(
    'label_seurat4_bmmc transfers predicted ADT without error', 
    {
        quiet_so_list <- purrr::quietly(SeuratLabeler::label_seurat4_bmmc)(
            query_so = so_query, 
            ref_so = ref_bm
        )
        res <- filter_quietly(quiet_so_list, patterns_suppress = exp_warn_pat_list)
        expected_assays = c('RNA','predicted_ADT')
        testthat::expect_true(all(expected_assays %in% Assays(res)))
        testthat::expect_true(nrow(res[['predicted_ADT']]@data) == nrow(ref_bm[['ADT']]@data))
    }
)

testthat::test_that(
    'label_seurat4_bmmc correctly errors on bad ref_reduction value', 
    {
        testthat::expect_error(
            {
                quiet_so_list <- SeuratLabeler::label_seurat4_bmmc(
                    query_so = so_query, 
                    ref_so = ref_bm,
                    ref_reduction = 'nonexistent_reduction')
            },
            regexp = "Reference dataset reduction '.*' not found in reference data object"
        )
    }
)

testthat::test_that(
    'label_seurat4_bmmc correctly errors on bad transfer data label name value', 
    {
        ref_data_list <- list(celltype.l4='celltype.l4')
        testthat::expect_error(
            {
                quiet_so_list <- SeuratLabeler::label_seurat4_bmmc(
                    query_so = so_query, 
                    ref_so = ref_bm,
                    ref_data = ref_data_list
                )
            },
            regexp = "Reference data to transfer not found in reference object: "
        )
    }
)

testthat::test_that(
    'label_seurat4_bmmc correctly errors on bad transfer data assay name value', 
    {
        ref_data_list <- list(predicted_CITE='CITE')
        testthat::expect_error(
            {
                quiet_so_list <- SeuratLabeler::label_seurat4_bmmc(
                    query_so = so_query, 
                    ref_so = ref_bm,
                    ref_data = ref_data_list
                )
            },
            regexp = "Reference data to transfer not found in reference object: "
        )
    }
)

