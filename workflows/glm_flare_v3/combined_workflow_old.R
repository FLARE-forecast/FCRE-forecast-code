#readRenviron("~/.Renviron") # MUST come first
library(tidyverse)
library(lubridate)
lake_directory <- here::here()
setwd(lake_directory)
forecast_site <<- "fcre"
configure_run_file <<- "configure_run.yml"
update_run_config <<- TRUE
config_set_name <<- "defaultV2"

Sys.setenv("AWS_DEFAULT_REGION" = "renc",
           "AWS_S3_ENDPOINT" = "osn.xsede.org",
           "USE_HTTPS" = TRUE)

config <- FLAREr::set_configuration(configure_run_file,lake_directory, config_set_name = config_set_name)

message("Generating targets")
source(file.path("workflows", config_set_name, "01_generate_targets.R"))

noaa_ready <- TRUE
while(noaa_ready){
  
  config <- FLAREr::set_configuration(configure_run_file,lake_directory, config_set_name = config_set_name)

  setwd(lake_directory)

  message("Generating inflow forecast")
  source(file.path("workflows", config_set_name, "02_run_inflow_forecast.R"))

  setwd(lake_directory)

  message("Generating forecast")
  source(file.path("workflows", config_set_name, "03_run_flarer_forecast.R"))

  setwd(lake_directory)

  RCurl::url.exists("https://hc-ping.com/a996a401-97b3-4884-a778-02243e056d2a", timeout = 5)
  
  noaa_ready <- FLAREr::check_noaa_present_arrow(lake_directory,
                                         configure_run_file,
                                         config_set_name = config_set_name)


}
