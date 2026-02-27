
#' Get duplicated caseIDs
#' 
#' Retrieve the duplicated caseIDs to remove from the analysis.
#'
#' @param duplicates_dir directory path where the text files with the duplicates information are.
#'
#' @return an integer vector with all the caseids to be removed
#' @export
get_duplicate_caseids <- function(duplicates_dir = NULL){
  
  if (is.null(duplicates_dir)){return(000000)}
  
  if (!dir.exists(duplicates_dir)) {
    stop("directory does not exist")
  }
  
  deleted_files <- list.files(duplicates_dir, full.names = TRUE, pattern = ".txt")
  caseid_to_be_removed <- c()
  
  for (d_file in deleted_files){
    read_file <- data.table::fread(d_file, sep = "$", quote = "", fill = FALSE)
    caseid_to_be_removed <- c(caseid_to_be_removed, unique(read_file$V1)) %>% 
      as.numeric()
    
  }
  
  return(caseid_to_be_removed)
  
}


#' Convert a date string into a date format
#'
#' @param date_string A string vector with multiple formats (8, 6 or 4 digits)
#'
#' @return A converted Date
#' @export
#' @examples arrange_date("2020")
#' arrange_date("202006")
#' arrange_date("20200601")
arrange_date <- function(date_string) {
  if (any(stringr::str_detect(
    string = date_string,
    pattern = "[:alpha:]"
  ), na.rm = TRUE)) {
    stop("string contains letters")
  }

  date_string_modif <- dplyr::case_when(
    (nchar(date_string) == 4) ~ paste0(date_string, "0630"),
    (nchar(date_string) == 6) ~ paste0(date_string, "15"),
    (nchar(date_string) == 8) ~ as.character(date_string),
    TRUE ~ as.character(date_string)
  )

  date_string_modif <- as.Date(date_string_modif, format = "%Y%m%d")
  return(unname(date_string_modif))
}

#' Read FAERS ascii files
#' 
#' Read ASCII files from a directory, removing the duplicates.
#'
#' @param ascii_dir directory path where ascii files are
#' @param cache_path (optional) a string. Must have a ".Rdata" extension to save the read tabular formats in each loop.
#' @param drug_indication_pattern (optional) a string.filter ADRs with a specific drug indication pattern (**stringr** sintax) 
#' @param drug_pattern (optional) a string. filter ADRs with a specific drug name pattern (**stringr** sintax) 
#' @param primary_suspect (optional) a string. 
#' @param ... directory with duplicate information to be passed to \link{get_duplicate_caseids}
#'
#' @return A list with binded tibbles retrieved from files.
#' @export
retrieve_faersascii <- function(ascii_dir, cache_path = NULL, drug_indication_pattern = NULL, drug_pattern = NULL,
                                 primary_suspect = TRUE, ...) {
  
  #setting global vars
  caseid <- role_cod <- drugname <- indi_pt <- primaryid <- indi_drug_seq <- 
    ichicsr <- . <- NULL
  
  
  if (!dir.exists(ascii_dir)) {
    stop("directory does not exist")
  }

  ascii_files <- list.files(ascii_dir, full.names = TRUE, pattern = ".txt")
  ascii_drug_files <- ascii_files[which(stringr::str_detect(ascii_files, pattern = "drug|DRUG|indi|INDI"))]
  
  
  #getduplicated caseIDs
  duplicated_caseids <- get_duplicate_caseids(...)

  # gather drug files
  drug_info <- tibble::tibble()
  indi_info <- tibble::tibble()
  drug_indi_info <- tibble::tibble()

  # iterating over drug and indication files to filter the chosen ADRs
  for (ascii_drug_file in ascii_drug_files) {
    message(paste(
      "Retrieving drug/indication information: binding", ascii_drug_file,
      which(ascii_drug_file == ascii_drug_files), "out of", length(ascii_drug_files)
    ))


    if (stringr::str_detect(ascii_drug_file, "drug|DRUG")) {
      
      #read drug_info
      read_file <- data.table::fread(ascii_drug_file, sep = "$", quote = "", fill = FALSE) %>% 
        mutate_at(c("primaryid", "caseid"), ~as.numeric(.)) %>% 
        filter(!caseid %in% duplicated_caseids) %>% 
        {
          if ("nda_num" %in% names(.)) {
            .$nda_num <- as.character(.$nda_num)
            .
          } else {
            .
          }
        } %>%
        {
          if ("dose_amt" %in% names(.)) {
            .$dose_amt <- as.character(.$dose_amt)
            .
          } else {
            .
          }
        } %>%
        mutate_at(vars(ends_with("_dt")), list(~ arrange_date(.))) 

      #filter drug_file and bind it
      drug_info <- bind_rows(drug_info, read_file  %>% 
        {
          if (primary_suspect) {
            filter(., role_cod == "PS")
          } else {
            .
          }
        } %>%
        {
          if (!is.null(drug_pattern)) {
            filter(., stringr::str_detect(string = drugname, pattern = drug_pattern))
          } else {
            .
          }
        })
      
      #remove file to free space
      rm(list= c("read_file"))
      gc()
      
    }

    if (stringr::str_detect(ascii_drug_file, "indi|INDI")) {
      
      #read indication_info
      read_file <- data.table::fread(ascii_drug_file, sep = "$", quote = "", fill = FALSE) %>% 
        mutate_at(c("primaryid", "caseid"), ~as.numeric(.)) %>% 
        filter(!caseid %in% duplicated_caseids) %>%
        mutate_at(vars(ends_with("_dt")), list(~ arrange_date(.)))

      #filter indi_file and bind it
      indi_info <- bind_rows(indi_info, read_file  %>%
        {
          if (!is.null(drug_indication_pattern)) {
            filter(., stringr::str_detect(string = indi_pt, pattern = drug_indication_pattern))
          } else {
            .
          }
        })
      
      #remove file to free space
      rm(list= c("read_file"))
      gc()
      
    }
  }

  #set final structure for drug_indi_info 
  if (!is.null(drug_indication_pattern) & !is.null(drug_pattern)) {
    drug_indi_info <- drug_info %>%
      inner_join(indi_info %>% group_by(primaryid, caseid, indi_drug_seq) %>%
        reframe(indi_pt = paste(sort(unique(indi_pt)), collapse = "/"), .groups = "keep") %>%
        ungroup(), by = c("primaryid", "caseid", "drug_seq" = "indi_drug_seq"))
  }

  if (!is.null(drug_indication_pattern) & is.null(drug_pattern)) {
    if (primary_suspect == TRUE) {
      drug_indi_info <- drug_info %>%
        inner_join(indi_info %>% group_by(primaryid, caseid, indi_drug_seq) %>%
          reframe(indi_pt = paste(sort(unique(indi_pt)), collapse = "/"), .groups = "keep") %>%
          ungroup(), by = c("primaryid", "caseid", "drug_seq" = "indi_drug_seq"))
    } else {
      drug_indi_info <- indi_info
    }
  }


  if (is.null(drug_indication_pattern) & !is.null(drug_pattern)) {
    drug_indi_info <- drug_info %>% 
      filter(!caseid %in% duplicated_caseids)
  }

  if (is.null(drug_indication_pattern) & is.null(drug_pattern)) {
    drug_indi_info <- drug_info %>%
      full_join(indi_info %>% group_by(primaryid, caseid, indi_drug_seq) %>%
        reframe(indi_pt = paste(sort(unique(indi_pt)), collapse = "/"), .groups = "keep") %>%
        ungroup(), by = c("primaryid", "caseid", "drug_seq" = "indi_drug_seq"))
  }



  #filter drug_info & indi_info
  drug_info <- drug_info %>% filter(primaryid %in% drug_indi_info$primaryid) %>% 
    mutate_at("drugname", ~stringr::str_replace_all(., pattern = "AND", replacement = "\\\\")) %>% 
    mutate_at("drugname", ~stringr::str_squish(.))
  indi_info <- indi_info %>% filter(primaryid %in% drug_indi_info$primaryid)



  # cache objects to be filled
  cache <- list()
  cache_files <- c()

  # files to be filled
  demo_info <- tibble::tibble()
  outc_info <- tibble::tibble()
  reac_info <- tibble::tibble()
  rpsr_info <- tibble::tibble()
  ther_info <- tibble::tibble()

  #read the files and remove the duplicated cases, impute date, and filter by caseIDs with specified drug and drug indication patterns
  for (ascii_file in ascii_files) {
    message(paste("binding", ascii_file, which(ascii_file == ascii_files), "out of", length(ascii_files)))
    
    #garbage collection
    gc()


    if (stringr::str_detect(ascii_file, "demo|DEMO")) {
      demo_info <- demo_info %>% bind_rows(data.table::fread(ascii_file, sep = "$", quote = "", fill = FALSE)%>% 
                                             mutate_at(c("primaryid", "caseid"), ~as.numeric(.)) %>% 
                                             filter(!caseid %in% duplicated_caseids)  %>%
        mutate_at(vars(ends_with("_dt")), list(~ arrange_date(.))) %>%
        {
          if ("nda_num" %in% names(.)) {
            .$nda_num <- as.character(.$nda_num)
            .
          } else {
            .
          }
        } %>%
        mutate_at(c("age", "wt"), ~ suppressWarnings(as.numeric(.))) %>%
        filter(primaryid %in% drug_indi_info$primaryid))
    }
    if (stringr::str_detect(ascii_file, "outc|OUTC")) {
      outc_info <- outc_info %>% bind_rows(data.table::fread(ascii_file, sep = "$", quote = "", fill = FALSE) %>% 
                                             mutate_at(c("primaryid", "caseid"), ~as.numeric(.)) %>% 
                                             filter(!caseid %in% duplicated_caseids) %>%
        {
          if ("outc_cod" %in% names(.)) {
            dplyr::rename(., "outc_code" = "outc_cod")
          } else {
            .
          }
        } %>%
        mutate_at(vars(ends_with("_dt")), list(~ arrange_date(.))) %>%
        filter(primaryid %in% drug_indi_info$primaryid))
    }
    if (stringr::str_detect(ascii_file, "reac|REAC")) {
      reac_info <- reac_info %>% bind_rows(data.table::fread(ascii_file, sep = "$", quote = "", fill = FALSE) %>% 
                                             mutate_at(c("primaryid", "caseid"), ~as.numeric(.)) %>% 
                                             filter(!caseid %in% duplicated_caseids) %>%
        mutate_at(vars(ends_with("_dt")), list(~ arrange_date(.))) %>%
        filter(primaryid %in% drug_indi_info$primaryid))
    }
    if (stringr::str_detect(ascii_file, "rpsr|RPSR")) {
      rpsr_info <- rpsr_info %>% bind_rows(data.table::fread(ascii_file, sep = "$", quote = "", fill = FALSE)%>% 
                                             mutate_at(c("primaryid", "caseid"), ~as.numeric(.)) %>% 
                                             filter(!caseid %in% duplicated_caseids)  %>%
        mutate_at(vars(ends_with("_dt")), list(~ arrange_date(.))) %>%
        filter(primaryid %in% drug_indi_info$primaryid)) 
    }
    if (stringr::str_detect(ascii_file, "ther|THER")) {
      ther_info <- ther_info %>% bind_rows(data.table::fread(ascii_file, sep = "$", quote = "", fill = FALSE)%>% 
                                             mutate_at(c("primaryid", "caseid"), ~as.numeric(.)) %>% 
                                             filter(!caseid %in% duplicated_caseids)  %>%
        mutate_at(vars(ends_with("_dt")), list(~ arrange_date(.))) %>%
        {
          if ("dur" %in% names(.)) {
            .$dur <- as.character(.$dur)
            .
          } else {
            .
          }
        } %>%
        filter(primaryid %in% drug_indi_info$primaryid))
    }


    if (!is.null(cache_path)) {
      message(paste("saving", ascii_file, "into cache:", cache_path))

      cache_files <- c(cache_files, ascii_file)
      ascii_list <- list(
        demographics = demo_info,
        drug = drug_info,
        indication = indi_info,
        outcome = outc_info,
        reaction = reac_info,
        report_source = rpsr_info,
        therapy = ther_info
      )
      cache$files <- cache_files
      save(cache, file = cache_path)
    }
  }
  return(list(
    demographics = demo_info,
    drug = drug_info,
    indication = indi_info,
    outcome = outc_info,
    reaction = reac_info,
    report_source = rpsr_info,
    therapy = ther_info
  ))
}



#' Unify the list to a tabular format
#' 
#' Turn the list elements returned from \link{retrieve_faersascii} into a tabular format
#'
#' @param ascii_list list from \link{retrieve_faersascii}
#' @examples
#' unify_tabular_ascii(ascii_list = als_faers_data)
#' 
#'
#' @return A data frame representing FAERS data, with all components from the list joined.
#' @export
unify_tabular_ascii <- function(ascii_list) {
  
  
  #setting global vars
  primaryid <- caseid <- indi_drug_seq <- indi_pt <- dsg_drug_seq <- drugname <- 
    role_cod <- start_dt <- reporter_country <- sex <- event_dt <- age <- age_YR <- 
    indi_pt_all <- drugname_all <- start_dt_ps <- caseversion <- NULL
  
  drug_indi_info <- ascii_list$drug %>%
    mutate_at("drugname", ~stringr::str_squish(.)) %>% 
    left_join(ascii_list$indication %>% group_by(primaryid, caseid, indi_drug_seq) %>%
      reframe(indi_pt = paste(sort(unique(indi_pt)), collapse = "/"), .groups = "keep") %>%
      ungroup(), by = c("primaryid", "caseid", "drug_seq" = "indi_drug_seq")) %>% 
    left_join(ascii_list$therapy %>% 
                distinct(primaryid, caseid, dsg_drug_seq, .keep_all = TRUE),
              by = c("primaryid", "caseid", "drug_seq" = "dsg_drug_seq"))

  outcome_info <- ascii_list$outcome %>%
    tidyr::pivot_wider(
      id_cols = c("primaryid", "caseid"), names_from = "outc_code",
      values_fill = "0", values_from = "outc_code", values_fn = ~ ifelse(. == "0", "0", "1"),
      names_prefix = "outcome_"
    )

  report_source_info <- ascii_list$report_source %>%
    tidyr::pivot_wider(
      id_cols = c("primaryid", "caseid"), names_from = "rpsr_cod",
      values_fill = "0", values_from = "rpsr_cod", values_fn = ~ ifelse(. == "0", "0", "1"),
      names_prefix = "report_source_"
    )
  
  unified_faers <- ascii_list$reaction %>%
    left_join(distinct(ascii_list$demographics, primaryid, caseid, .keep_all = TRUE), 
              by = c("primaryid", "caseid"), relationship = "many-to-many") %>%
    mutate(patient_drug = lapply(as.integer(primaryid), FUN = function(x) {
      drug_indi_info %>%
        filter(primaryid == x) %>%
        mutate_at(c("prod_ai", "drugname", "indi_pt"), ~stringr::str_squish(.)) %>% 
        suppressMessages()
    })) %>% 
    left_join(outcome_info, by = c("primaryid", "caseid")) %>%
    left_join(report_source_info, by = c("primaryid", "caseid")) %>%
    left_join(drug_indi_info %>%
      group_by(primaryid, caseid) %>%
      reframe(
        indi_pt_all = paste(sort(unique(indi_pt)), collapse = " + "),
        drugname_all = paste(sort(unique(drugname)), collapse = " + "),
        .groups = "keep"
      ) %>%
      ungroup(), by = c("primaryid", "caseid")) %>%
    left_join(drug_indi_info %>%
      filter(role_cod == "PS") %>%
      group_by(primaryid, caseid) %>%
      reframe(
        indi_pt_ps = sort(unique(indi_pt)),
        drugname_ps = paste(sort(unique(drugname)), collapse = " + "),
        start_dt_ps = unique(start_dt),
        .groups = "keep"
      ) %>%
      ungroup() %>% suppressWarnings(), by = c("primaryid", "caseid"))  %>% 
    mutate(age_YR = case_when(age_cod == "YR" ~ as.numeric(age),
                              age_cod == "DEC" ~ as.numeric(age)*10,
                              age_cod == "MON" ~ as.numeric(age)/12,
                              age_cod == "WK" ~ as.numeric(age)/(365/7),
                              age_cod == "DY" ~ as.numeric(age)/365,
                              age_cod == "HR" ~ as.numeric(age)/(365*24),
                              TRUE ~ as.numeric(age))) %>% 
    relocate(age_YR, .after = "age")

  
  #de-duplication
  unified_faers <- unified_faers %>% 
    mutate_at("pt", ~stringr::str_squish(.))
  
  #latest case version
  unified_faers <- unified_faers %>%
    group_by(caseid) %>%
    filter(caseversion == max(caseversion)) %>%
    ungroup()
  # 
  # #remove suspected duplicated cases
  index_to_remove <- which(duplicated(unified_faers [,c( "reporter_country",
                                                          "sex", "event_dt",
                                                          "age", "pt",
                                                          "indi_pt_all",
                                                          "drugname_all",
                                                          "start_dt_ps")]) &
                              !is.na(unified_faers$reporter_country) &
                              !is.na(unified_faers$sex) &
                              !is.na(unified_faers$event_dt) &
                              !is.na(unified_faers$age) &
                              !is.na(unified_faers$pt) &
                              !is.na(unified_faers$indi_pt_all) &
                              !is.na(unified_faers$drugname_all) &
                              !is.na(unified_faers$start_dt_ps)
  )

  if(length(index_to_remove) >=1){unified_faers <- unified_faers[-index_to_remove,]}


  #trim
  message("unification and de-duplication applied.")
  
  return(unified_faers)
}



