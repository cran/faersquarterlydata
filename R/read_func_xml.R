#' Convert FAERS xml to an R list
#'
#' @param xml_address XML address file
#'
#' @return a list containing all the elements from `xml_address`
#' @export
faersxml_to_r <- function(xml_address) {
  if (!stringr::str_detect(xml_address, pattern = "\\.xml")) {
    stop("xml address should have an .xml extension")
  }

  data_2 <- xml2::as_list(xml2::read_xml(x = xml_address))

  return(data_2)
}


#' Convert FAERS xml to tabular format
#'
#' @param xml_address XML address to be read
#' @param drug_indication_pattern filter by ADR with a specific drug indication pattern (**stringr** sintax)
#' @param reaction_wise each row corresponds to a reaction (if TRUE, drug_wise cannot be TRUE)
#' @param drug_wise each row corresponds to a drug (if TRUE, reaction_wise cannot be TRUE)
#'
#' @return A tibble corresponding to the XML file
#' @export
#'
retrieve_faersxml <- function(xml_address, reaction_wise = TRUE, drug_wise = FALSE, drug_indication_pattern = NULL) {
  
  #setting global variables
  ichicsr <- . <- NULL
  
  message(paste("retrieving", xml_address))


  if (reaction_wise == TRUE & drug_wise == TRUE) {
    stop("reaction wise == TRUE  or drug wise == TRUE, not both")
  }
  if (reaction_wise == FALSE & drug_wise == FALSE) {
    stop("reaction wise == FALSE  or drug wise == FALSE, not both")
  }


  # read xml
  data_2 <- faersxml_to_r(xml_address)

  # turn it into a tibble
  xml_df <- tibble::as_tibble(data_2) %>%
    slice(-1)

  # tidyr::unnest the tibble
  xml_df3 <- xml_df %>% tidyr::unnest_wider(ichicsr, names_repair = "unique")

  # tidyr::unnest more columns
  xml_df4 <- tibble()
  for (i in 1:nrow(xml_df3)) {

    d <- xml_df3 %>%
      slice(i)

    drug_indications <- c()
    if (!is.null(drug_indication_pattern)) {
      drugs <- d$patient[[1]][which(names(d$patient[[1]]) == "drug")]

      for (x in drugs) {
        drug_indications <- c(drug_indications, x$drugindication[[1]])
      }

      if (!any(stringr::str_detect(drug_indications, drug_indication_pattern))) {
        next
      }
    }

    d <- d %>%
      select(-starts_with("message")) %>%
      tidyr::unnest_wider(c("primarysource", "sender", "receiver"), names_sep = "_") %>%
      tidyr::unnest_wider("patient",
        names_sep = "_", simplify = FALSE,
        names_repair = function(x) {
          reaction_i <- stringr::str_which(x, "patient_reaction")
          x[reaction_i] <- paste0(x[reaction_i], "_", 1:length(reaction_i))

          drug_i <- stringr::str_which(x, "patient_drug")
          x[drug_i] <- paste0(x[drug_i], "_", 1:length(drug_i))

          return(x)
        }
      )

    if (reaction_wise == TRUE) {
      # set patient_reaction in one col
      d <- d %>%
        mutate(patient_reaction = apply(.[, stringr::str_detect(names(.), "reaction")],
          MARGIN = 1,
          FUN = function(x) {
            reactionmeddraversionpt <- c()
            reactionmeddrapt <- c()
            reactionoutcome <- c()

            for (i in x) {
              l <- unlist(i)

              reactionmeddraversionpt <- c(reactionmeddraversionpt, l["reactionmeddraversionpt"])
              reactionmeddrapt <- c(reactionmeddrapt, l["reactionmeddrapt"])
              reactionoutcome <- c(reactionoutcome, l["reactionoutcome"])
            }

            return(list(
              "reactionmeddraversionpt" = reactionmeddraversionpt,
              "reactionmeddrapt" = reactionmeddrapt,
              "reactionoutcome" = reactionoutcome
            ))
          }
        )) %>%
        mutate(patient_drug = apply(.[, stringr::str_detect(names(.), "drug")],
          MARGIN = 1,
          FUN = function(x) {
            drugcharacterization <- c()
            medicinalproduct <- c()
            drugauthorizationnumb <- c()
            drugstructuredosagenumb <- c()
            drugstructuredosageunit <- c()
            drugseparatedosagenumb <- c()
            drugintervaldosageunitnumb <- c()
            drugintervaldosagedefinition <- c()
            drugdosagetext <- c()
            drugdosageform <- c()
            drugadministrationroute <- c()
            drugstartdateformat <- c()
            drugstartdate <- c()
            drugenddateformat <- c()
            drugenddate <- c()
            drugindication <- c()
            actiondrug <- c()
            drugrecurreadministration <- c()
            drugadditional <- c()
            activesubstance <- c()

            for (i in x) {
              l <- unlist(i)

              drugcharacterization <- c(drugcharacterization, l["drugcharacterization"])
              medicinalproduct <- c(medicinalproduct, l["medicinalproduct"])
              drugauthorizationnumb <- c(drugauthorizationnumb, l["drugauthorizationnumb"])
              drugstructuredosagenumb <- c(drugstructuredosagenumb, l["drugstructuredosagenumb"])
              drugstructuredosageunit <- c(drugstructuredosageunit, l["drugstructuredosageunit"])
              drugseparatedosagenumb <- c(drugseparatedosagenumb, l["drugseparatedosagenumb"])
              drugintervaldosageunitnumb <- c(drugintervaldosageunitnumb, l["drugintervaldosageunitnumb"])
              drugintervaldosagedefinition <- c(drugintervaldosagedefinition, l["drugintervaldosagedefinition"])
              drugdosagetext <- c(drugdosagetext, l["drugdosagetext"])
              drugdosageform <- c(drugdosageform, l["drugdosageform"])
              drugadministrationroute <- c(drugadministrationroute, l["drugadministrationroute"])
              drugstartdateformat <- c(drugstartdateformat, l["drugstartdateformat"])
              drugstartdate <- c(drugstartdate, l["drugstartdate"])
              drugenddateformat <- c(drugenddateformat, l["drugenddateformat"])
              drugenddate <- c(drugenddate, l["drugenddate"])
              drugindication <- c(drugindication, l["drugindication"])
              actiondrug <- c(actiondrug, l["actiondrug"])
              drugrecurreadministration <- c(drugrecurreadministration, l["drugrecurreadministration"])
              drugadditional <- c(drugadditional, l["drugadditional"])
              activesubstance <- c(activesubstance, paste(l[names(l)[stringr::str_which(names(l), "activesubstancename")]], collapse = " + "))
            }
            data <- tibble::tibble(
              "drugcharacterization" = drugcharacterization,
              "medicinalproduct" = medicinalproduct,
              "drugauthorizationnumb" = drugauthorizationnumb,
              "drugstructuredosagenumb" = drugstructuredosagenumb,
              "drugstructuredosageunit" = drugstructuredosageunit,
              "drugseparatedosagenumb" = drugseparatedosagenumb,
              "drugintervaldosageunitnumb" = drugintervaldosageunitnumb,
              "drugintervaldosagedefinition" = drugintervaldosagedefinition,
              "drugdosagetext" = drugdosagetext,
              "drugdosageform" = drugdosageform,
              "drugadministrationroute" = drugadministrationroute,
              "drugstartdateformat" = drugstartdateformat,
              "drugstartdate" = drugstartdate,
              "drugenddateformat" = drugenddateformat,
              "drugenddate" = drugenddate,
              "drugindication" = drugindication,
              "actiondrug" = actiondrug,
              "drugrecurreadministration" = drugrecurreadministration,
              "drugadditional" = drugadditional,
              "activesubstance" = activesubstance
            )

            return(data)
          }
        )) %>%
        tidyr::unnest_wider(c("patient_reaction")) %>%
        tidyr::unnest_longer(col = c(
          "reactionmeddraversionpt", "reactionmeddrapt", "reactionoutcome"
          # "drugcharacterization","medicinalproduct","drugauthorizationnumb","drugstructuredosagenumb",
          # "drugstructuredosageunit" , "drugdosagetext" ,"drugadministrationroute" ,"drugindication","actiondrug",
          # "drugrecurreadministration" , "drugadditional", "activesubstance"
        )) %>%
        select(-names(.)[stringr::str_detect(names(.), pattern = "_id|duplicate|(\\.\\.\\.)")]) %>%
        select(-names(.)[stringr::str_detect(names(.), pattern = "patient_drug_\\d")]) %>%
        select(-names(.)[stringr::str_detect(names(.), pattern = "patient_reaction_\\d")]) %>%
        tidyr::unnest(names(.)[-which(names(.) == "patient_drug")])
    }



    if (drug_wise == TRUE) {
      # tidyr::unnest longer for reaction into one reaction per row, and a tibble of all the drugs for each reaction
      d <- d %>%
        mutate(patient_drug = apply(.[, stringr::str_detect(names(.), "drug")],
          MARGIN = 1,
          FUN = function(x) {
            drugcharacterization <- c()
            medicinalproduct <- c()
            drugauthorizationnumb <- c()
            drugstructuredosagenumb <- c()
            drugstructuredosageunit <- c()
            drugseparatedosagenumb <- c()
            drugintervaldosageunitnumb <- c()
            drugintervaldosagedefinition <- c()
            drugdosagetext <- c()
            drugdosageform <- c()
            drugadministrationroute <- c()
            drugstartdateformat <- c()
            drugstartdate <- c()
            drugenddateformat <- c()
            drugenddate <- c()
            drugindication <- c()
            actiondrug <- c()
            drugrecurreadministration <- c()
            drugadditional <- c()
            activesubstance <- c()

            for (i in x) {
              l <- unlist(i)


              drugcharacterization <- c(drugcharacterization, l["drugcharacterization"])
              medicinalproduct <- c(medicinalproduct, l["medicinalproduct"])
              drugauthorizationnumb <- c(drugauthorizationnumb, l["drugauthorizationnumb"])
              drugstructuredosagenumb <- c(drugstructuredosagenumb, l["drugstructuredosagenumb"])
              drugstructuredosageunit <- c(drugstructuredosageunit, l["drugstructuredosageunit"])
              drugseparatedosagenumb <- c(drugseparatedosagenumb, l["drugseparatedosagenumb"])
              drugintervaldosageunitnumb <- c(drugintervaldosageunitnumb, l["drugintervaldosageunitnumb"])
              drugintervaldosagedefinition <- c(drugintervaldosagedefinition, l["drugintervaldosagedefinition"])
              drugdosagetext <- c(drugdosagetext, l["drugdosagetext"])
              drugdosageform <- c(drugdosageform, l["drugdosageform"])
              drugadministrationroute <- c(drugadministrationroute, l["drugadministrationroute"])
              drugstartdateformat <- c(drugstartdateformat, l["drugstartdateformat"])
              drugstartdate <- c(drugstartdate, l["drugstartdate"])
              drugenddateformat <- c(drugenddateformat, l["drugenddateformat"])
              drugenddate <- c(drugenddate, l["drugenddate"])
              drugindication <- c(drugindication, l["drugindication"])
              actiondrug <- c(actiondrug, l["actiondrug"])
              drugrecurreadministration <- c(drugrecurreadministration, l["drugrecurreadministration"])
              drugadditional <- c(drugadditional, l["drugadditional"])
              activesubstance <- c(activesubstance, paste(l[names(l)[stringr::str_which(names(l), "activesubstancename")]], collapse = " + "))
            }
            data <- tibble::tibble(
              "drugcharacterization" = drugcharacterization,
              "medicinalproduct" = medicinalproduct,
              "drugauthorizationnumb" = drugauthorizationnumb,
              "drugstructuredosagenumb" = drugstructuredosagenumb,
              "drugstructuredosageunit" = drugstructuredosageunit,
              "drugseparatedosagenumb" = drugseparatedosagenumb,
              "drugintervaldosageunitnumb" = drugintervaldosageunitnumb,
              "drugintervaldosagedefinition" = drugintervaldosagedefinition,
              "drugdosagetext" = drugdosagetext,
              "drugdosageform" = drugdosageform,
              "drugadministrationroute" = drugadministrationroute,
              "drugstartdateformat" = drugstartdateformat,
              "drugstartdate" = drugstartdate,
              "drugenddateformat" = drugenddateformat,
              "drugenddate" = drugenddate,
              "drugindication" = drugindication,
              "actiondrug" = actiondrug,
              "drugrecurreadministration" = drugrecurreadministration,
              "drugadditional" = drugadditional,
              "activesubstance" = activesubstance
            )

            return(data)
          }
        )) %>%
        mutate(patient_reaction = apply(.[, stringr::str_detect(names(.), "reaction")],
          MARGIN = 1,
          FUN = function(x) {
            reactionmeddraversionpt <- c()
            reactionmeddrapt <- c()
            reactionoutcome <- c()

            for (i in x) {
              l <- unlist(i)

              reactionmeddraversionpt <- c(reactionmeddraversionpt, l["reactionmeddraversionpt"])
              reactionmeddrapt <- c(reactionmeddrapt, l["reactionmeddrapt"])
              reactionoutcome <- c(reactionoutcome, l["reactionoutcome"])
            }

            return(tibble::tibble(
              "reactionmeddraversionpt" = reactionmeddraversionpt,
              "reactionmeddrapt" = reactionmeddrapt,
              "reactionoutcome" = reactionoutcome
            ))
          }
        )) %>%
        tidyr::unnest_wider(c("patient_drug")) %>%
        tidyr::unnest_longer(col = c(
          "drugcharacterization", "medicinalproduct", "drugauthorizationnumb", "drugstructuredosagenumb",
          "drugstructuredosageunit", "drugseparatedosagenumb", "drugintervaldosageunitnumb", "drugintervaldosagedefinition",
          "drugdosagetext", "drugdosageform", "drugadministrationroute",
          "drugstartdateformat", "drugstartdate", "drugenddateformat", "drugenddate",
          "drugindication", "actiondrug",
          "drugrecurreadministration", "drugadditional", "activesubstance"
        )) %>%
        select(-names(.)[stringr::str_detect(names(.), pattern = "_id|duplicate|\\.\\.\\.")]) %>%
        select(-names(.)[stringr::str_detect(names(.), pattern = "patient_drug_\\d")]) %>%
        select(-names(.)[stringr::str_detect(names(.), pattern = "patient_reaction_\\d")]) %>%
        tidyr::unnest(names(.)[-which(names(.) == "patient_reaction")])
    }


    # build
    cols_to_unnest <- intersect(
      names(d),
      c(
        "safetyreportversion", "safetyreportid", "primarysourcecountry",
        "occurcountry", "transmissiondateformat", "transmissiondate", "reporttype",
        "serious", "seriousnessdeath", "seriousnesslifethreatening", "seriousnesshospitalization",
        "seriousnessdisabling", "seriousnesscongenitalanomali", "seriousnessother",
        "receivedateformat", "receivedate", "receiptdateformat", "receiptdate",
        "fulfillexpeditecriteria", "companynumb", "primarysource_reportercountry",
        "primarysource_qualification", "sender_sendertype", "sender_senderorganization",
        "receiver_receivertype", "receiver_receiverorganization", "authoritynumb",
        "patient_patientsex", "patient_patientonsetage",
        "patient_patientonsetageunit", "patient_summary", "patient_patientweight",
        "patient_patientagegroup", "primarysource_literaturereference"
      )
    )

    d <- d %>% tidyr::unnest(cols = cols_to_unnest)
    xml_df4 <- bind_rows(xml_df4, d)
  }



  return(xml_df4)
}



#' Convert FAERS a number of xml files to tabular format
#'
#' @param xml_address_vector Vector with XML addresses to be read
#' @param cache_path a string. Must have a ".Rdata" extension to save the read tabular formats in each loop.
#' @param ... arguments to be passed to \link{retrieve_faersxml}
#'
#' @return A binded tibble with all the tibbles returned from `retrieve_faersxml`
#' @export
retrieve_faersxml_all <- function(xml_address_vector, ..., cache_path = NULL) {
  if (!all(file.exists(xml_address_vector))) {
    stop("not all files exist")
  }

  data <- tibble::tibble()

  cache <- list()
  cache_adress <- c()


  for (address in xml_address_vector) {
    message(paste(which(address == xml_address_vector), "out of", length(xml_address_vector)))

    d <- retrieve_faersxml(xml_address = address, ...)

    data <- data %>% bind_rows(d)

    if (!is.null(cache_path)) {
      message(paste("saving", address, "into cache:", cache_path))

      cache_adress <- c(cache_adress, address)
      cache$data <- data
      cache$xml_adress <- cache_adress
      save(cache, file = cache_path)
    }
  }
  return(data)
}


#' Retrieve unique drug and ADR information values from XML files
#'
#' @param xml_address_vector Vector with XML addresses to be read
#' @param ... arguments to be passed to \link{retrieve_faersxml}
#'
#' @return A list with all the unique information on FAERS variables
#' @export
retrieve_unique_info <- function(xml_address_vector, ...) {
  if (!all(file.exists(xml_address_vector))) {
    stop("not all files exist")
  }

  data <- list()

  adr_unique <- c()
  medicinal_product_unique <- c()
  drug_indication_unique <- c()

  for (address in xml_address_vector) {
    message(paste(which(address == xml_address_vector), "out of", length(xml_address_vector)))

    d <- retrieve_faersxml(xml_address = address, reaction_wise = TRUE)

    adr_unique <- c(adr_unique, d$reactionmeddrapt) %>% unique()
    medicinal_product_unique <- c(medicinal_product_unique, unlist(lapply(d$patient_drug, FUN = function(x) {
      x$medicinalproduct
    }))) %>% unique()
    drug_indication_unique <- c(drug_indication_unique, unlist(lapply(d$patient_drug, FUN = function(x) {
      x$drugindication
    }))) %>% unique()
  }

  data$adr_unique <- adr_unique
  data$medicinal_product_unique <- medicinal_product_unique
  data$drug_indication_unique <- drug_indication_unique

  return(data)
}
