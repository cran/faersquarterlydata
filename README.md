<!-- README.md is generated from README.Rmd. Please edit that file -->

# `faersquarterlydata`

The goal of `faersquarterlydata` is to provide an easy framework to read
and analyse FAERS XML/ASCII files. The package faersquarterlydata for R
programming language provides easy access and analysis to FDA Adverse
Event Report System (FAERS) database. This database contains information
on the reported Adverse Drug Events (ADRs) in the United States since
2004. The available data format in FDA website is in XML or ASCII
format, and therefore, the users need to be familiar with creation of
relational databases. This package allows the reading of these files and
transform them into tabular format, computing summary counts and
estimating some useful statistics like the Reporting Odds Ratio (ROR)
and Proportional Reporting Ratio (PRR), and therefore, enabling
reproducible research on this topic.

## Installation

You can install the development version of `faersquarterlydata` like so:


    install.packages(faersquarterlydata)  

## License

The license is GPL-3 (<https://cran.r-project.org/web/licenses/GPL-3>).

## FDA Adverse Event Reporting System

The Latest Quarterly Data Files from FAERS can be retrieved here:
<https://www.fda.gov/drugs/questions-and-answers-fdas-adverse-event-reporting-system-faers/fda-adverse-event-reporting-system-faers-latest-quarterly-data-files>

## Basic Usage

### First Step: Unzipping the .zip folders

FAERS database files are typically distributed in .zip files which
contain text files within. In order to facilitate the opening of these
files, we provided here this function:

    unzip_faerszip(zip_folders_dir= "directory_with_zip_files", ex_dir = "directory_with_text_files")

### Second Step: Reading and Filtering files

Each quarterly ASCII file will result in seven tables containing diverse
information. In order for the Demographic information and others to be
binded into one single table, and the same for the other types of text
files, the following function is available:

    als_faers_data <- retrieve_faersascii(ascii_dir = "directory_with_text_files/ascii", drug_indication_pattern =  "Amyothrophic lateral sclerosis|Motor neuron disease", primary_suspect = TRUE, duplicates_dir = "directory_with_text_files/deleted" )

### Third Step: Unify text files into a single table

In order to merge all these seven tables into one, and therefore, allow
more meaningful analysis, the package makes available the following
function:

    als_faers_data_unified <- unify_tabular_ascii(ascii_list = als_faers_data)

### Fourth Step: Database description

The filtered database can be described based on demographic information,
drug-related characteristics, ADR description, report source, outcome or
counts based on the date of the event. This description is computed,
partly, by `tableone` package . The following code was used to describe
the filtered database:

    summary_faers <- summary_faersdata(als_faers_data_unified)
