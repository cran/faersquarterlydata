#' Unzip FAERS zip folders
#'
#' @param zip_folders_dir directory containing FAERS zip folders
#' @param ex_dir directory to be exported the unzipped files
#'
#' @return None. Just unzips the folders to a specified location.
#' @export
unzip_faerszip <- function(zip_folders_dir, ex_dir) {
  for (zip_folder in zip_folders_dir) {
    zip_files <- list.files(zip_folder, full.names = TRUE, 
                            pattern = ".zip")

    for (zipped_file in zip_files) {
      message(paste("unzipping", zipped_file))
      unzip(zipped_file, exdir = ex_dir)
    }
  }
}
