testthat::test_that("faersxml_to_r", {
  testthat::expect_error(object = faersxml_to_r(xml_address = "any_other_xml_adress.txt"))
  testthat::expect_true( is.list(faersxml_to_r(xml_address = "resources/1_ADR22Q1_example.xml")))
})


testthat::test_that("retrieve_faersxml", {
  testthat::expect_error(object = retrieve_faersxml(xml_address = "resources/1_ADR22Q1_example.xml", drug_wise = TRUE, reaction_wise = TRUE))
  testthat::expect_error(object = retrieve_faersxml(xml_address = "resources/1_ADR22Q1_example.xml", drug_wise = FALSE, reaction_wise = FALSE))
  testthat::expect_true( object = is.data.frame(retrieve_faersxml(xml_address = "resources/1_ADR22Q1_example.xml", drug_indication_pattern = "Parkinson")) )
  testthat::expect_true( object = is.data.frame(retrieve_faersxml(xml_address = "resources/1_ADR22Q1_example.xml")) )
  testthat::expect_true( object = is.data.frame(retrieve_faersxml(xml_address = "resources/1_ADR22Q1_example.xml", drug_wise = TRUE, reaction_wise = FALSE)) )
})


testthat::test_that("retrieve_faersxml_all", {
  testthat::expect_error(object = retrieve_faersxml_all(xml_address_vector = c("file", "other_file")))
  testthat::expect_true( object = is.data.frame(retrieve_faersxml_all(xml_address = "resources/1_ADR22Q1_example.xml" )))
  testthat::expect_true( object = is.data.frame(retrieve_faersxml_all(xml_address = "resources/1_ADR22Q1_example.xml" , cache_path = "resources/cache/cache_xml.Rdata")))
})

testthat::test_that("retrieve_unique_info", {
  testthat::expect_error(object = retrieve_unique_info(xml_address_vector = c("file", "other_file")))
  testthat::expect_true( object = is.list(retrieve_unique_info(xml_address_vector = "resources/1_ADR22Q1_example.xml" )))
})
