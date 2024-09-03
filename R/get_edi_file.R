#' Download file from EDI data portal
#'
#' @param edi_https https of the EDI package
#' @param file name of the file in the EDI package (not full path)
#' @param lake_directory full path to repository directory
#'
#' @export
#'
get_edi_file <- function(edi_https, file, lake_directory){ #, curl_timeout = 60){

  if(!file.exists(file.path(lake_directory, "data_raw", file))){
    if(!dir.exists(dirname(file.path(lake_directory, "data_raw", file)))){
      dir.create(dirname(file.path(lake_directory, "data_raw", file)))
    }
    url_download <- httr::RETRY("GET",edi_https, httr::timeout(1500), pause_base = 5, pause_cap = 20, pause_min = 5, times = 3, quiet = FALSE)
    test_bin <- httr::content(url_download,'raw')
    writeBin(test_bin, file.path(lake_directory, "data_raw", file))
  }
}
