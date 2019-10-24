context("data")
library(SeuratLabeler)

test_data <- Read10X_h5(system.file("testdata/pbmc_1k_v3.h5",
                                    package = "SeuratLabeler"))

test_that(
  "test H5 data can be imported",
  {
    expect_true("dgCMatrix" %in% class(test_data))
  }
)
