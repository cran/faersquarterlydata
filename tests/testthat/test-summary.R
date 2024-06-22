testthat::test_that("summary_faersdata", {
  
  faers_ascii_data <- retrieve_faersascii(ascii_dir = "resources/test_zip_ex_dir/ASCII")
  faers_ascii_data_unified <- unify_tabular_ascii(ascii_list = faers_ascii_data)
  testthat::expect_true(is.list(summary_faersdata(faers_ascii_data_unified)))
})
