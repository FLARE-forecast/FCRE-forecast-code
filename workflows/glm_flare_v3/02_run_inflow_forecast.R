#renv::restore()

library(tidyverse)
library(lubridate)

lake_directory <- here::here()

files.sources <- list.files(file.path(lake_directory, "R"), full.names = TRUE)
sapply(files.sources, source)

configure_run_file <- "configure_run.yml"
config <- FLAREr::set_configuration(configure_run_file,lake_directory, config_set_name = config_set_name)

FLAREr::get_targets(lake_directory, config)


message("Forecasting inflow and outflows")
# Forecast Inflows

if(config$run_config$forecast_horizon > 0){

  inflow_model_coeff <- NULL
  inflow_model_coeff$future_inflow_flow_coeff <- c(0.0010803, 0.9478724, 0.3478991)
  inflow_model_coeff$future_inflow_flow_error <- 0.00965
  inflow_model_coeff$future_inflow_temp_coeff <- c(0.20291, 0.94214, 0.04278)
  inflow_model_coeff$future_inflow_temp_error <- 0.943

  temp_flow_forecast <- forecast_inflows_outflows_arrow(inflow_obs = file.path(config$file_path$qaqc_data_directory, "fcre-targets-inflow.csv"),
                                                        obs_met_file = file.path(config$file_path$qaqc_data_directory,"observed-met_fcre.csv"),
                                                        inflow_model = config$inflow$forecast_inflow_model,
                                                        inflow_process_uncertainty = FALSE,
                                                        inflow_model_coeff = inflow_model_coeff,
                                                        site_id = config$location$site_id,
                                                        use_s3_met = TRUE,
                                                        use_s3_inflow = config$run_config$use_s3,
                                                        met_bucket = config$s3$drivers$bucket,
                                                        met_endpoint = config$s3$drivers$endpoint,
                                                        inflow_bucket = config$s3$inflow_drivers$bucket,
                                                        inflow_endpoint = config$s3$inflow_drivers$endpoint,
                                                        inflow_local_directory = file.path(lake_directory, "drivers/inflow"),
                                                        forecast_start_datetime = config$run_config$forecast_start_datetime,
                                                        forecast_horizon = config$run_config$forecast_horizon)

message(paste0("successfully generated inflow forecats for: ", file.path(config$met$forecast_met_model,config$location$site_id,lubridate::as_date(config$run_config$forecast_start_datetime))))

}else{

  message("An inflow forecasts was not needed because the forecast horizon was 0 in run configuration file")

}
