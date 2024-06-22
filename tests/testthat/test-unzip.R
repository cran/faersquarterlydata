test_that("unzip_faerszip", {
  testthat::expect_no_error(unzip_faerszip(zip_folders_dir = "resources/xml_files",
                                           ex_dir = "resources/test_zip_ex_dir/"))
  
  testthat::expect_no_error(unzip_faerszip(zip_folders_dir = "resources/ascii_files/",
                                           ex_dir = "resources/test_zip_ex_dir/"))
})
