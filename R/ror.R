

#' Estimate Reporting Odds Ratio
#'
#' @param n11 Number of events of interest within the group of interest
#' @param n10 Number of events of interest from all groups
#' @param n01 Number of all events within the group of interest
#' @param n00 Number of all events from all groups
#' @param n00 Number of all events from all groups
#' @param ic_range Confidence Interval range
#'
#' @return list with ROR estimate and a vector with the IC boundaries
#' @export
#' @import dplyr
#'
#' @examples
#' estimate_ror(n11 = 20, n10 = 10, n01 = 200, n00 = 200, ic_range = 0.90)
estimate_ror <- function(n11, n10, n01, n00, ic_range = 0.95){

  n11 <- as.numeric(n11)
  n10 <- as.numeric(n10)
  n01 <- as.numeric(n01)
  n00 <- as.numeric(n00)
  
  ror <- list()
  ic_range <- stats::qnorm(ic_range + (1-ic_range)/2) * sqrt(1/n11 + 1/n01 + 1/n10 + 1/n00)

  ror$estimate <- round((n11/n10) / (n01/n00), 2)
  ror$ic <- c(lower = round(exp( log(ror$estimate) - ic_range ), 2),
              upper = round(exp( log(ror$estimate) + ic_range), 2))

  return(ror)

}


#' Estimate Proportional Reporting Odds Ratio
#'
#' @param n11 Number of events of interest within the group of interest
#' @param n10 Number of events of interest from all groups
#' @param n01 Number of all events within the group of interest
#' @param n00 Number of all events from all groups
#' @param n00 Number of all events from all groups
#' @param ic_range Confidence Interval range
#'
#' @return Proportional Reporting Odds Ratio
#' @export
#' @import dplyr
#'
#' @examples
#' estimate_prr(n11 = 20, n10 = 10, n01 = 200, n00 = 200)
estimate_prr <- function(n11, n10, n01, n00, ic_range = 0.95){
  
  n11 <- as.numeric(n11)
  n10 <- as.numeric(n10)
  n01 <- as.numeric(n01)
  n00 <- as.numeric(n00)
  
  prr <- list()
  ic_range <- stats::qnorm(ic_range + (1-ic_range)/2) * sqrt(1/n11 + 1/(n11 + n01) + 1/n10 + 1/(n10 + n00))
  
  prr$estimate <- (n11/(n11 + n01)) / (n10/ (n10 + n00))
  prr$ic <- c(lower = round(exp( log(prr$estimate) - ic_range ), 2),
              upper = round(exp( log(prr$estimate) + ic_range), 2)) 
  
  return(prr)
  
}



#' Estimate Chi-Squared test with yates correction
#'
#' @param n11 Number of events of interest within the group of interest
#' @param n10 Number of events of interest from all groups
#' @param n01 Number of all events within the group of interest
#' @param n00 Number of all events from all groups
#' @param n00 Number of all events from all groups
#'
#' @return list with Chi-squared statistic and p-value
#' @export
#' @import dplyr
#'
#' @examples
#' estimate_chisq(n11 = 20, n10 = 10, n01 = 200, n00 = 200)
estimate_chisq <- function(n11, n10, n01, n00){
  
  n11 <- as.numeric(n11)
  n10 <- as.numeric(n10)
  n01 <- as.numeric(n01)
  n00 <- as.numeric(n00)
  
  mx <- matrix(data = c(n11, n10, n01, n00), ncol = 2, byrow = TRUE)
  
  ch <- chisq.test(mx, correct = TRUE) %>% suppressWarnings()
  
  chisq <- unname(ch$statistic)
  
  return(chisq)
  
}

#' Estimate Information Component
#'
#' @param n11 Number of events of interest within the group of interest
#' @param n10 Number of events of interest from all groups
#' @param n01 Number of all events within the group of interest
#' @param n00 Number of all events from all groups
#' @param n00 Number of all events from all groups
#'
#' @return List with Information Component estimate and its 0.95 IC
#' @export
#' @import dplyr
#'
#' @examples
#' estimate_infoc(n11 = 20, n10 = 10, n01 = 200, n00 = 200)
estimate_infoc <- function(n11, n10, n01, n00){
  
  n11 <- as.numeric(n11)
  n10 <- as.numeric(n10)
  n01 <- as.numeric(n01)
  n00 <- as.numeric(n00)
  
  infoc <- list()
  
  n11_exp <- ((n11 + n10) * (n11 + n01)) / (n11 + n10 + n01 + n00)
  
  infoc$estimate <- log2((n11 + 0.5)/ (n11_exp + 0.5))
  infoc$ic <- c(lower = infoc$estimate - 3.3 * (n11 + 0.5)^(-1/2) - 2 * (n11 +0.5)^(-3/2),
                upper = infoc$estimate + 2.4 * (n11 + 0.5)^(-1/2) - 0.5 * (n11 +0.5)^(-3/2)) 
  
  return(infoc)
  
}


#' Estimate Measures of Association
#'
#' @param tabular_faers_data FAERS tabular format. Output of function \link{retrieve_faersxml} or \link{retrieve_faersxml_all}
#' @param group_of_interest_col  a string, specifying the group of interest. 
#' Must me a column name of `tabular_faers_data`, and this columns should only contain two unique values.
#' @param group_of_interest_ref a string, specifying the group of interest reference. Must me a value from the group of interest column.
#' @param rename_vector optional. named vector to rename the group of interest, in order to show up in a 
#' @param event_of_interest_col a string, specifying the event of interest. Must me a column name of `tabular_faers_data`.
#' @param ... arguments passed to `estimate_ror` like `ic_range`.
#'
#' @return tibble with the event of interest counts, group of interest counts and the respective estimated measures of association 
#' (ROR and its IC, PRR and its IC, Information Component and Chi-squared statisti with Yates correction.
#' @export
#' @import dplyr
#'
#' @examples 
#' estimate_ror_bygroup(tabular_faers_data = dplyr::filter(als_faers_data_unified, 
#' sex %in% c("M", "F")),
#' group_of_interest_col = "sex", 
#' group_of_interest_ref = "M", 
#' event_of_interest_col = "pt") 
#' 
#' 
estimate_ror_bygroup <- function(tabular_faers_data, 
                          group_of_interest_col = NULL, 
                          group_of_interest_ref = NULL,
                          rename_vector = NULL,
                          event_of_interest_col = NULL,
                          ...){
  
  if(length(unique(tabular_faers_data[[group_of_interest_col]]))>2){stop("group of interest should only have two values")}
  if(!group_of_interest_ref %in% unique(tabular_faers_data[[group_of_interest_col]]) & !is.logical(tabular_faers_data[[group_of_interest_col]])){
    stop("group of interest reference group must be one value from the group of interest")
  }
  
  if(is.logical(tabular_faers_data[[event_of_interest_col]])){
    
    vars_of_interest <- c(event_of_interest_col, 
                          group_of_interest_ref, 
                          setdiff(as.character(unique(tabular_faers_data[[group_of_interest_col]])),group_of_interest_ref ) )
    
    tabular_faers_data %>% 
      group_by_at(vars(c(group_of_interest_col,  event_of_interest_col))) %>%
      count() %>%
      ungroup() %>%
      tidyr::pivot_wider(names_from =group_of_interest_col, values_from = n, values_fill = 0) %>%
      select(all_of(vars_of_interest)) %>%
      {if(is.vector(rename_vector)){ rename(., rename_vector) }else{.}} %>% 
      mutate(total_events = rowSums(.[,c(2,3)]),
             n01 = sum(.[[2]]) - .[[2]],
             n00 = sum(.[[3]]) - .[[3]]) %>% 
      {names(.)[5] <-  paste0("other_events_", names(.)[2]); .} %>% 
      {names(.)[6] <-  paste0("other_events_", names(.)[3]); .} %>% 
      mutate(ror = apply(., MARGIN = 1, FUN = function(x) { estimate_ror(n11 = x[2], n10 = x[3], n01 = x[5], n00 = x[6])})) %>% 
      mutate(prr = apply(., MARGIN = 1, FUN = function(x) { estimate_prr(n11 = x[2], n10 = x[3], n01 = x[5], n00 = x[6])})) %>% 
      mutate(infoc = apply(., MARGIN = 1, FUN = function(x) { estimate_infoc(n11 = x[2], n10 = x[3], n01 = x[5], n00 = x[6])})) %>% 
      mutate(chisq = apply(., MARGIN = 1, FUN = function(x) { estimate_chisq(n11 = x[2], n10 = x[3], n01 = x[5], n00 = x[6])})) %>% 
      tidyr::unnest_wider(c("ror", "prr", "infoc"), names_sep = "_") %>% 
      tidyr::unnest_wider(ends_with("ic") , names_sep = "_") %>% 
      filter(.[,1] == TRUE) %>% 
      select(-1) %>% 
      return()
    
  }else{
    
    vars_of_interest <- c(event_of_interest_col, 
                          as.character(group_of_interest_ref), 
                          setdiff(as.character(unique(tabular_faers_data[[group_of_interest_col]])),
                                  as.character(group_of_interest_ref ) ))
    
    
    tabular_faers_data %>% 
      group_by_at(vars(c(group_of_interest_col,  event_of_interest_col))) %>%
      count() %>%
      ungroup() %>%  
      tidyr::pivot_wider(names_from =group_of_interest_col, values_from = n, values_fill = 0) %>%
      select(all_of(vars_of_interest)) %>%
      {if(is.vector(rename_vector)){ rename(., rename_vector) }else{.}} %>%
      mutate(total_events = rowSums(.[,c(2,3)])) %>% 
      mutate(n01 = sum(.[[2]]) - .[[2]]) %>% 
      mutate(n00 = sum(.[[3]]) - .[[3]]) %>% 
      {names(.)[5] <-  paste0("other_events_", names(.)[2]); .} %>% 
      {names(.)[6] <-  paste0("other_events_", names(.)[3]); .} %>% 
      mutate(ror = apply(., MARGIN = 1, FUN = function(x) { estimate_ror(n11 = x[2], n10 = x[3], n01 = x[5], n00 = x[6])})) %>% 
      mutate(prr = apply(., MARGIN = 1, FUN = function(x) { estimate_prr(n11 = x[2], n10 = x[3], n01 = x[5], n00 = x[6])})) %>% 
      mutate(infoc = apply(., MARGIN = 1, FUN = function(x) { estimate_infoc(n11 = x[2], n10 = x[3], n01 = x[5], n00 = x[6])})) %>% 
      mutate(chisq = apply(., MARGIN = 1, FUN = function(x) { estimate_chisq(n11 = x[2], n10 = x[3], n01 = x[5], n00 = x[6])})) %>% 
      tidyr::unnest_wider(c("ror", "prr", "infoc"), names_sep = "_") %>% 
      tidyr::unnest_wider(ends_with("ic") , names_sep = "_") %>% 
      arrange(desc("total_events")) %>% 
      return()
    
    
    
  }
  
  
  
  
}

