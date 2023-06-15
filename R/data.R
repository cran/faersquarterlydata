
#' List ASCII data example
#'
#' A list containing data from FDA website. The list only contains safety reports which the ADR primary suspect drug was indicated for ALS.
#' List originated from \link{retrieve_faersascii}
#'
#' @format A data frame with 200 rows and 38 columns:
#' @source <https://www.fda.gov/drugs/questions-and-answers-fdas-adverse-event-reporting-system-faers/fda-adverse-event-reporting-system-faers-latest-quarterly-data-files>
"als_faers_data"


#' Tabular ASCII data example
#'
#' A subset of data from FAERS data. One row corresponds to one adverse drug reaction. All the ADR in this subset have a primary suspect drug indicated for ALS.
#' Data frame originated from \link{unify_tabular_ascii}
#'
#' @format A data frame with 1635 rows and 40 columns.
"als_faers_data_unified"


#' List of approved products by FDA
#'
#' @format A data frame.
#' @source <https://www.fda.gov/drugs/drug-approvals-and-databases/drugsfda-data-files>
"products_fda"

