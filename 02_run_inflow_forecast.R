renv::restore()

library(tidyverse)
library(lubridate)

lake_directory <- here::here()
s3_mode <- TRUE

files.sources <- list.files(file.path(lake_directory, "R"), full.names = TRUE)
sapply(files.sources, source)

configuration_file <- "configure_flare.yml"
run_config <- yaml::read_yaml(file.path(lake_directory,"configuration","FLAREr","configure_run.yml"))
forecast_site <- run_config$forecast_site

if(file.exists("~/.aws")){
  warning(paste("Detected existing AWS credentials file in ~/.aws,",
                "Consider renaming these so that automated upload will work"))
}

Sys.setenv("AWS_DEFAULT_REGION" = "s3",
           "AWS_S3_ENDPOINT" = "flare-forecast.org")

#Note: lake_directory need to be set prior to running this script
config <- yaml::read_yaml(file.path(lake_directory,"configuration","FLAREr",configuration_file))
config$file_path$qaqc_data_directory <- file.path(lake_directory, "data_processed")
config$file_path$data_directory <- file.path(lake_directory, "data_raw")

if(s3_mode){
  restart_exists <- aws.s3::object_exists(object = file.path(forecast_site, "configure_run.yml"), bucket = "restart")
  if(restart_exists){
    aws.s3::save_object(object = file.path(forecast_site, "configure_run.yml"), bucket = "restart", file = file.path(lake_directory,"configuration","FLAREr","configure_run.yml"))
  }
  config$file_path$noaa_directory <- file.path(lake_directory, "drivers", "noaa")
  config$file_path$inflow_directory <- file.path(lake_directory, "drivers", "inflow")
}else{
  config$file_path$noaa_directory <- file.path(dirname(lake_directory), "drivers", "noaa")
  config$file_path$inflow_directory <- file.path(dirname(lake_directory), "drivers", "inflow")
}
run_config <- yaml::read_yaml(file.path(lake_directory,"configuration","FLAREr","configure_run.yml"))
config$run_config <- run_config


if(config$run_config$forecast_horizon > 0){
  # Set up timings
  #Weather Drivers
  start_datetime <- lubridate::as_datetime(config$run_config$start_datetime)
  if(is.na(config$run_config$forecast_start_datetime)){
    end_datetime <- lubridate::as_datetime(config$run_config$end_datetime)
    forecast_start_datetime <- end_datetime
  }else{
    forecast_start_datetime <- lubridate::as_datetime(config$run_config$forecast_start_datetime)
    end_datetime <- forecast_start_datetime + lubridate::days(config$run_config$forecast_horizon)
  }
  forecast_hour <- lubridate::hour(forecast_start_datetime)
  if(forecast_hour < 10){forecast_hour <- paste0("0",forecast_hour)}
  forecast_path <- file.path(config$file_path$noaa_directory, "NOAAGEFS_1hr",config$location$site_id,lubridate::as_date(forecast_start_datetime),forecast_hour)


  #Weather Drivers

  noaa_forecast_path <- file.path(config$file_path$noaa_directory, config$met$forecast_met_model,config$location$site_id,lubridate::as_date(forecast_start_datetime),forecast_hour)

  if(s3_mode){
    aws.s3::save_object(object = file.path(forecast_site, "fcre-targets-inflow.csv"), bucket = "targets", file = file.path(config$file_path$qaqc_data_directory, "fcre-targets-inflow.csv"))
    aws.s3::save_object(object = file.path(forecast_site, "observed-met_fcre.nc"), bucket = "targets", file = file.path(config$file_path$qaqc_data_directory, "observed-met_fcre.nc"))

    if(config$run_config$forecast_horizon > 0){
      noaa_forecast_path <- file.path(lake_directory,"drivers/noaa", config$met$forecast_met_model,config$location$site_id,lubridate::as_date(forecast_start_datetime),forecast_hour)

      download_s3_objects(lake_directory,
                          bucket = "drivers",
                          prefix = file.path("noaa", config$met$forecast_met_model,config$location$site_id,lubridate::as_date(forecast_start_datetime),forecast_hour))
    }
  }else{
    local_noaa_forecast_path <- file.path(config$file_path$noaa_directory, config$met$forecast_met_model,config$location$site_id,lubridate::as_date(forecast_start_datetime),forecast_hour)
    noaa_forecast_path <- file.path(lake_directory, "drivers/noaa", config$met$forecast_met_model,config$location$site_id,lubridate::as_date(forecast_start_datetime),forecast_hour)
    files <- list.files(noaa_forecast_path, full.names = TRUE)
    for(i in 1:length(files)){
      dir.create(noaa_forecast_path)
      file.copy(from = files[i], to = noaa_forecast_path)
    }
  }

  message("Forecasting inflow and outflows")
  # Forecast Inflows

  config$future_inflow_flow_coeff <- c(0.0010803, 0.9478724, 0.3478991)
  config$future_inflow_flow_error <- 0.00965
  config$future_inflow_temp_coeff <- c(0.20291, 0.94214, 0.04278)
  config$future_inflow_temp_error <- 0.943

  forecast_files <- list.files(noaa_forecast_path, full.names = TRUE)
  if(length(forecast_files) == 0){
    stop(paste0("missing forecast files at: ", noaa_forecast_path))
  }
  temp_flow_forecast <- forecast_inflows_outflows(inflow_obs = file.path(config$file_path$qaqc_data_directory, "fcre-targets-inflow.csv"),
                                                  forecast_files = forecast_files,
                                                  obs_met_file = file.path(config$file_path$qaqc_data_directory,"observed-met_fcre.nc"),
                                                  output_dir = config$file_path$inflow_directory,
                                                  inflow_model = config$inflow$forecast_inflow_model,
                                                  inflow_process_uncertainty = FALSE,
                                                  forecast_location = config$file_path$forecast_output_directory,
                                                  config = config,
                                                  s3_mode = s3_mode,
                                                  bucket = "drivers")

  #NEED TO COPY TO BUCKET

  if(config$model_settings$model_name == "glm_aed"){

    run_dir <- str_replace_all(dirname(list.files(temp_flow_forecast[[1]], full.names = TRUE)[1]), "INFLOW-FLOWS-NOAAGEFS-AR1","INFLOW-FLOWS-NOAAGEFS-AR1-AED")

    if(!dir.exists(run_dir)){
      dir.create(run_dir, recursive = TRUE)
    }

    historical_chemistry_weir <- read_csv("/Users/quinn/Downloads/FCRE-forecast-code/data_processed/FCR_weir_inflow_2013_2019_20200624_allfractions_2poolsDOC.csv")

    for(i in 1:length(inflow_files)){
      inflow_files <- list.files(temp_flow_forecast[[1]], full.names = TRUE)
      inflow_forecast <- read_csv(inflow_files[i])

      historical_chemistry_weir_mean <- historical_chemistry_weir %>%
        mutate(doy = lubridate::yday(time)) %>%
        select(-c("FLOW", "TEMP", "SALT", "time")) %>%
        group_by(doy) %>%
        summarise(across(everything(), mean))

      file_name <- str_replace_all(list.files(temp_flow_forecast[[1]], full.names = TRUE)[i], "INFLOW-FLOWS-NOAAGEFS-AR1","INFLOW-FLOWS-NOAAGEFS-AR1-AED")

      inflow_forecast_aed <- inflow_forecast %>%
        mutate(doy = lubridate::yday(time)) %>%
        left_join(historical_chemistry_weir_mean, by = "doy") %>%
        select(-doy) %>%
        write_csv(file = file_name)

      #NEED TO COPY TO BUCKET
    }
  }else if(config$model_settings$model_name == "glm_oxy"){

    #NEED TO GET WORKING WITH S3

    run_dir <- str_replace_all(dirname(list.files(temp_flow_forecast, full.names = TRUE)[1]), "INFLOW-FLOWS-NOAAGEFS-AR1","INFLOW-FLOWS-NOAAGEFS-AR1-AED")

    if(!dir.exists(run_dir)){
      dir.create(run_dir, recursive = TRUE)
    }

    historical_chemistry_weir <- read_csv(file.path(config$file_path$qaqc_data_directory, "FCR_weir_inflow_2013_2019_20200624_allfractions_2poolsDOC.csv"))

    for(i in 1:length(inflow_files)){
      inflow_files <- list.files(temp_flow_forecast, full.names = TRUE)
      inflow_forecast <- read_csv(inflow_files[i])

      historical_chemistry_weir_mean <- historical_chemistry_weir %>%
        mutate(doy = lubridate::yday(time)) %>%
        select(-c("FLOW", "TEMP", "SALT", "time")) %>%
        select("OXY_oxy") %>%
        group_by(doy) %>%
        summarise(across(everything(), mean))

      file_name <- str_replace_all(list.files(temp_flow_forecast, full.names = TRUE)[i], "INFLOW-FLOWS-NOAAGEFS-AR1","INFLOW-FLOWS-NOAAGEFS-AR1-AEDOXY")

      inflow_forecast_aed <- inflow_forecast %>%
        mutate(doy = lubridate::yday(time)) %>%
        left_join(historical_chemistry_weir_mean, by = "doy") %>%
        select(-doy) %>%
        write_csv(file = file_name)

      #NEED TO COPY TO BUCKET
    }

    unlink(noaa_forecast_path, recursive = TRUE)
    if(s3_mode){
      unlink(dirname((temp_flow_forecast[[1]])[1]), recursive = TRUE)
    }
  }

  message(paste0("successfully generated inflow forecats for: ", file.path(config$met$forecast_met_model,config$location$site_id,lubridate::as_date(forecast_start_datetime),forecast_hour)))

}else{

  message(paste0("no forecast days in configuration file: ", file.path(config$met$forecast_met_model,config$location$site_id,lubridate::as_date(forecast_start_datetime),forecast_hour)))

}


