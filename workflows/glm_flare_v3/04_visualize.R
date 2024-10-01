#renv::restore()

library(tidyverse)
library(lubridate)

lake_directory <- here::here()
files.sources <- list.files(file.path(lake_directory, "R"), full.names = TRUE)
sapply(files.sources, source)

configure_run_file <- "configure_run.yml"

config <- FLAREr::set_configuration(configure_run_file,lake_directory, config_set_name = config_set_name)

config <- FLAREr::get_restart_file(config, lake_directory)

FLAREr::get_targets(lake_directory, config)

pdf_file <-   FLAREr::plotting_general_2(file_name = saved_file,
                                         target_file = file.path(config$file_path$qaqc_data_directory, paste0(config$location$site_id, "-targets-insitu.csv")),
                                         ncore = 2,
                                         obs_csv = FALSE)

#if(config$run_config$use_s3){
#  success <- aws.s3::put_object(file = pdf_file,
#                                object = file.path(config$location$site_id, basename(pdf_file)),
#                                bucket = "analysis",
#                                region = "s3",
#                                base_url = "flare-forecast.org",
#                                use_https = as.logical(Sys.getenv("USE_HTTPS")))
#  if(success){
#    unlink(pdf_file)
#  }
#}

png_file_name <- manager_plot(file_name = config$run_config$restart_file,
                              target_file = file.path(config$file_path$qaqc_data_directory, paste0(config$location$site_id, "-targets-insitu.csv")),
                              focal_depths = c(1, 5, 8))

if(config$run_config$use_s3 & !is.na(png_file_name)){
  success <- aws.s3::put_object(file = png_file_name,
                                object = file.path(config$location$site_id, basename(png_file_name)),
                                bucket = "analysis",
                                region = "s3",
                                base_url = "flare-forecast.org",
                                use_https = as.logical(Sys.getenv("USE_HTTPS")))
  if(success){
    unlink(png_file_name)
  }
}

if(config$run_config$use_s3){
  unlink(file.path(config$file_path$qaqc_data_directory, paste0(config$location$site_id, "-targets-insitu.csv")))
  unlink(config$run_config$restart_file)
}
