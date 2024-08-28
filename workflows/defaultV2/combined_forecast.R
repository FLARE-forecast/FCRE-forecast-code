library(tidyverse)
library(lubridate)

remotes::install_github("FLARE-forecast/FLAREr")

lake_directory <- here::here()
setwd(lake_directory)
forecast_site <- "fcre"
configure_run_file <- "configure_run.yml"
config_set_name <- "defaultV2"

Sys.setenv("AWS_DEFAULT_REGION" = "renc",
           "AWS_S3_ENDPOINT" = "osn.xsede.org",
           "USE_HTTPS" = TRUE)

FLAREr::ignore_sigpipe()

message("Beginning generate targets")
source(file.path(lake_directory, "workflows", config_set_name, "01_generate_targets.R"))

noaa_ready <- TRUE

while(noaa_ready){

  # if (as.Date(config$run_config$forecast_start_datetime) == Sys.Date()){
  #   stop('Stop here to keep action from failing')
  # }

  # ## add temporary check to catch if noaa data DNE
  # forecast_dir <- arrow::s3_bucket(bucket = file.path(config$s3$drivers$bucket,
  #                                                     paste0("stage2/reference_datetime=",as.Date(config$run_config$forecast_start_datetime)),
  #                                                     paste0("site_id=",config$location$site_id)),
  #                                  endpoint_override =  config$s3$drivers$endpoint, anonymous = TRUE)
  #
  # t <- try(arrow::open_dataset(forecast_dir) |> collect(), silent = TRUE)
  # if("try-error" %in% class(t)){
  #   stop('NOAA DATA NOT FOUND FOR THIS DAY')
  # }


  message("Generating inflow forecast")
  source(file.path(lake_directory, "workflows", config_set_name, "02_run_inflow_forecast.R"))

  # Run FLARE
  output <- FLAREr::run_flare(lake_directory = lake_directory,
                              configure_run_file = configure_run_file,
                              config_set_name = config_set_name)


  forecast_start_datetime <- lubridate::as_datetime(config$run_config$forecast_start_datetime) + lubridate::days(1)
  start_datetime <- lubridate::as_datetime(config$run_config$forecast_start_datetime) - lubridate::days(1)
  restart_file <- paste0(config$location$site_id,"-", (lubridate::as_date(forecast_start_datetime)- days(1)), "-",config$run_config$sim_name ,".nc")

  FLAREr::update_run_config2(lake_directory = lake_directory,
                             configure_run_file = configure_run_file,
                             restart_file = restart_file,
                             start_datetime = start_datetime,
                             end_datetime = NA,
                             forecast_start_datetime = forecast_start_datetime,
                             forecast_horizon = config$run_config$forecast_horizon,
                             sim_name = config$run_config$sim_name,
                             site_id = config$location$site_id,
                             configure_flare = config$run_config$configure_flare,
                             configure_obs = config$run_config$configure_obs,
                             use_s3 = config$run_config$use_s3,
                             bucket = config$s3$warm_start$bucket,
                             endpoint = config$s3$warm_start$endpoint,
                             use_https = TRUE)

  RCurl::url.exists("https://hc-ping.com/a996a401-97b3-4884-a778-02243e056d2a", timeout = 5)

  noaa_ready <- FLAREr::check_noaa_present_arrow(lake_directory,
                                                 configure_run_file,
                                                 config_set_name = config_set_name)
}
