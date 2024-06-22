test_that("estimate_ror", {
  testthat::expect_equal(estimate_ror(n11 = 20, n10 = 10, n01 = 200, n00 = 200)$estimate, 2)
  testthat::expect_equal(estimate_ror(n11 = 20, n10 = 10, n01 = 200, n00 = 50)$estimate, 0.5)
  testthat::expect_true(is.vector(estimate_ror(n11 = 20, n10 = 10, n01 = 200, n00 = 200)$ic))
  testthat::expect_true(is.list(estimate_ror(n11 = 20, n10 = 10, n01 = 200, n00 = 50)))
})

test_that("estimate_ppr", {
  testthat::expect_no_error(estimate_prr(n11 = 20, n10 = 10, n01 = 200, n00 = 200))
})

test_that("estimate_chisq", {
  testthat::expect_no_error(estimate_chisq(n11 = 20, n10 = 10, n01 = 200, n00 = 200))
})

test_that("estimate_infoc", {
  testthat::expect_no_error(estimate_infoc(n11 = 20, n10 = 10, n01 = 200, n00 = 200))
})

test_that("estimate_ror_bygroup", {
  
  faers_ascii_data <- retrieve_faersascii(ascii_dir = "resources/test_zip_ex_dir/ASCII")
  faers_ascii_data_unified <- unify_tabular_ascii(ascii_list = faers_ascii_data)
  
  testthat::expect_error(estimate_ror_bygroup(tabular_faers_data = faers_ascii_data,
                                    group_of_interest_col = "sex", 
                                    group_of_interest_ref = "1", 
                                    event_of_interest_col = "pt"))
  
  testthat::expect_error(estimate_ror_bygroup(tabular_faers_data = faers_ascii_data_unified %>% 
                                      filter(sex %in% c("F", "M")),
                                    group_of_interest_col = "Z", 
                                    group_of_interest_ref = "M", 
                                    event_of_interest_col = "pt"))
  
  log_v <- sample(c(TRUE, FALSE), size = nrow(faers_ascii_data_unified) ,replace = TRUE, prob = c(0.9, 0.1))
  testthat::expect_true(estimate_ror_bygroup(tabular_faers_data = faers_ascii_data_unified %>% 
                                      mutate(ei = log_v)  %>% 
                                      filter(sex %in% c("F", "M")) ,
                                    group_of_interest_col = "sex", 
                                    group_of_interest_ref = "M", 
                                    event_of_interest_col = "ei") %>% is.data.frame())
  
  testthat::expect_error(estimate_ror_bygroup(tabular_faers_data = faers_ascii_data_unified %>% 
                                               mutate(ei = log_v)  ,
                                             group_of_interest_col = "sex", 
                                             group_of_interest_ref = "M", 
                                             event_of_interest_col = "ei") %>% is.data.frame())
  
  testthat::expect_true(estimate_ror_bygroup(tabular_faers_data = faers_ascii_data_unified %>% 
                                               filter(sex %in% c("F", "M")) ,
                                   group_of_interest_col = "sex", 
                                   group_of_interest_ref = "M", 
                                   event_of_interest_col = "pt") %>% is.data.frame())
  
  
})
  