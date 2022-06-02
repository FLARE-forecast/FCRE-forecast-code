library(tidyverse)
library(lubridate)
set.seed(100)

remotes::install_github("rqthomas/GLM3r")
Sys.setenv("AWS_DEFAULT_REGION" = "s3",
           "AWS_S3_ENDPOINT" = "flare-forecast.org")

lake_directory <- here::here()

starting_index <- 1

files.sources <- list.files(file.path(lake_directory, "R"), full.names = TRUE)
sapply(files.sources, source)

sim_names <- "ms_glmead_flare"
config_set_name <- "glm_aed_ms"

config_files <- paste0("configure_flare_glm_aed.yml")

num_forecasts <- 52 * 3 - 3
#num_forecasts <- 1#19 * 7 + 1
days_between_forecasts <- 7
forecast_horizon <- 16 #32
starting_date <- as_date("2018-07-20")
#second_date <- as_date("2020-12-01") - days(days_between_forecasts)
#starting_date <- as_date("2018-07-20")
second_date <- as_date("2019-01-01") - days(days_between_forecasts)

start_dates <- rep(NA, num_forecasts)
start_dates[1:2] <- c(starting_date, second_date)
for(i in 3:(3 + num_forecasts)){
  start_dates[i] <- as_date(start_dates[i-1]) + days(days_between_forecasts)
}

start_dates <- as_date(start_dates)
forecast_start_dates <- start_dates + days(days_between_forecasts)
forecast_start_dates <- forecast_start_dates[-1]

configure_run_file <- "configure_aed_run.yml"

j = 1
sites <- "fcre"

#function(i, sites, lake_directory, sim_names, config_files, )

message(paste0("Running site: ", sites[j]))

if(starting_index == 1){
  run_config <- yaml::read_yaml(file.path(lake_directory, "configuration", config_set_name, configure_run_file))
  run_config$configure_flare <- config_files[j]
  run_config$sim_name <- sim_names
  yaml::write_yaml(run_config, file = file.path(lake_directory, "configuration", config_set_name, configure_run_file))
  if(file.exists(file.path(lake_directory, "restart", sites[j], sim_names, configure_run_file))){
    unlink(file.path(lake_directory, "restart", sites[j], sim_names, configure_run_file))
  }
}

##'
# Set up configurations for the data processing
config_obs <- FLAREr::initialize_obs_processing(lake_directory, observation_yml = "observation_processing.yml", config_set_name = config_set_name)
config <- FLAREr::set_configuration(configure_run_file,lake_directory, config_set_name = config_set_name)

use_s3 <- FALSE

#' Clone or pull from data repositories

FLAREr::get_git_repo(lake_directory,
                     directory = config_obs$realtime_insitu_location,
                     git_repo = "https://github.com/FLARE-forecast/FCRE-data.git")

FLAREr::get_git_repo(lake_directory,
                     directory = config_obs$realtime_met_station_location,
                     git_repo = "https://github.com/FLARE-forecast/FCRE-data.git")

FLAREr::get_git_repo(lake_directory,
                     directory = config_obs$realtime_inflow_data_location,
                     git_repo = "https://github.com/FLARE-forecast/FCRE-data.git")

#get_git_repo(lake_directory,
#             directory = config_obs$manual_data_location,
#             git_repo = "https://github.com/FLARE-forecast/FCRE-data.git")

#' Download files from EDI

FLAREr::get_edi_file(edi_https = "https://pasta.lternet.edu/package/data/eml/edi/389/5/3d1866fecfb8e17dc902c76436239431",
                     file = config_obs$met_raw_obs_fname[2],
                     lake_directory)

FLAREr::get_edi_file(edi_https = "https://pasta.lternet.edu/package/data/eml/edi/271/5/c1b1f16b8e3edbbff15444824b65fe8f",
                     file = config_obs$insitu_obs_fname[2],
                     lake_directory)

FLAREr::get_edi_file(edi_https = "https://pasta.lternet.edu/package/data/eml/edi/198/8/336d0a27c4ae396a75f4c07c01652985",
                     file = config_obs$secchi_fname,
                     lake_directory)

FLAREr::get_edi_file(edi_https = "https://pasta.lternet.edu/package/data/eml/edi/200/11/d771f5e9956304424c3bc0a39298a5ce",
                     file = config_obs$ctd_fname,
                     lake_directory)

FLAREr::get_edi_file(edi_https = "https://pasta.lternet.edu/package/data/eml/edi/199/8/da174082a3d924e989d3151924f9ef98",
                     file = config_obs$nutrients_fname,
                     lake_directory)


FLAREr::get_edi_file(edi_https = "https://pasta.lternet.edu/package/data/eml/edi/202/7/f5fa5de4b49bae8373f6e7c1773b026e",
                     file = config_obs$inflow_raw_file1[2],
                     lake_directory)

FLAREr::get_edi_file(edi_https = "https://pasta.lternet.edu/package/data/eml/edi/542/1/791ec9ca0f1cb9361fa6a03fae8dfc95",
                     file = "silica_master_df.csv",
                     lake_directory)

FLAREr::get_edi_file(edi_https = "https://pasta.lternet.edu/package/data/eml/edi/551/5/38d72673295864956cccd6bbba99a1a3",
                     file = "Dissolved_CO2_CH4_Virginia_Reservoirs.csv",
                     lake_directory)

https_file <- "https://raw.githubusercontent.com/cayelan/FCR-GLM-AED-Forecasting/master/FCR_2013_2019GLMHistoricalRun_GLMv3beta/inputs/FCR_SSS_inflow_2013_2021_20211102_allfractions_2DOCpools.csv"
download.file(https_file,
              file.path(config$file_path$execute_directory, basename(https_file)))


#' Clean up observed meterology

cleaned_met_file <- met_qaqc(realtime_file = file.path(config_obs$file_path$data_directory, config_obs$met_raw_obs_fname[1]),
                             qaqc_file = file.path(config_obs$file_path$data_directory, config_obs$met_raw_obs_fname[2]),
                             cleaned_met_file = file.path(config_obs$file_path$targets_directory, config_obs$site_id,paste0("observed-met_",config_obs$site_id,".nc")),
                             input_file_tz = "EST",
                             nldas = NULL)

#' Clean up observed inflow

cleaned_inflow_file <- inflow_qaqc(realtime_file = file.path(config_obs$file_path$data_directory, config_obs$inflow_raw_file1[1]),
                                   qaqc_file = file.path(config_obs$file_path$data_directory, config_obs$inflow_raw_file1[2]),
                                   nutrients_file = file.path(config_obs$file_path$data_directory, config_obs$nutrients_fname),
                                   silica_file = file.path(config_obs$file_path$data_directory,  config_obs$silica_fname),
                                   co2_ch4 = file.path(config_obs$file_path$data_directory, config_obs$ch4_fname),
                                   cleaned_inflow_file = file.path(config_obs$file_path$targets_directory, config_obs$site_id, paste0(config_obs$site_id,"-targets-inflow.csv")),
                                   input_file_tz = 'EST')

#' Clean up observed insitu measurements

cleaned_insitu_file <- in_situ_qaqc(insitu_obs_fname = file.path(config_obs$file_path$data_directory,config_obs$insitu_obs_fname),
                                    data_location = config_obs$file_path$data_directory,
                                    maintenance_file = file.path(config_obs$file_path$data_directory,config_obs$maintenance_file),
                                    ctd_fname = file.path(config_obs$file_path$data_directory, config_obs$ctd_fname),
                                    nutrients_fname =  file.path(config_obs$file_path$data_directory, config_obs$nutrients_fname),
                                    secchi_fname = file.path(config_obs$file_path$data_directory, config_obs$secchi_fname),
                                    ch4_fname = file.path(config_obs$file_path$data_directory, config_obs$ch4_fname),
                                    cleaned_insitu_file = file.path(config_obs$file_path$targets_directory, config_obs$site_id, paste0(config_obs$site_id,"-targets-insitu.csv")),
                                    lake_name_code = config_obs$site_id,
                                    config = config_obs)

##` Download NOAA forecasts`

message("    Downloading NOAA data")

cycle <- "00"

for(i in 1:length(forecast_start_dates)){
  noaa_forecast_path <- file.path(config$met$forecast_met_model, config$location$site_id, forecast_start_dates[i], "00")
  if(length(list.files(file.path(lake_directory,"drivers", noaa_forecast_path))) == 0){
  FLAREr::get_driver_forecast(lake_directory, forecast_path = noaa_forecast_path)
  }
}

available_dates <- list.files(file.path(lake_directory,"drivers","noaa","NOAAGEFS_1hr","fcre"))

if(starting_index == 1){
  config$run_config$start_datetime <- as.character(paste0(start_dates[1], " 00:00:00"))
  config$run_config$forecast_start_datetime <- as.character(paste0(start_dates[2], " 00:00:00"))
  config$run_config$forecast_horizon <- 0
  config$run_config$restart_file <- NA
  run_config <- config$run_config
  yaml::write_yaml(run_config, file = file.path(config$file_path$configuration_directory, configure_run_file))
}

#for(i in 1:1){
for(i in starting_index:length(forecast_start_dates)){

    config <- FLAREr::set_configuration(configure_run_file, lake_directory, config_set_name = config_set_name)

    num_dates_skipped <- 1
    while(!lubridate::as_date(config$run_config$forecast_start_datetime) %in% lubridate::as_date(available_dates) & i <= length(forecast_start_dates)){
      print("here")
      FLAREr::update_run_config(config, lake_directory, configure_run_file, saved_file = NA, new_horizon = forecast_horizon, day_advance = num_dates_skipped * days_between_forecasts, new_start_datetime = FALSE)
      config <- FLAREr::set_configuration(configure_run_file, lake_directory, config_set_name = config_set_name)
      num_dates_skipped <- num_dates_skipped + 1
      i <- i + 1
    }

      config <- FLAREr::set_configuration(configure_run_file,lake_directory, config_set_name = config_set_name)
      config <- FLAREr::get_restart_file(config, lake_directory)

      message(paste0("     Running forecast that starts on: ", config$run_config$start_datetime))

      if(config$run_config$forecast_horizon > 0){
        noaa_forecast_path <- FLAREr::get_driver_forecast_path(config,
                                                               forecast_model = config$met$forecast_met_model)
        forecast_dir <- file.path(config$file_path$noaa_directory, noaa_forecast_path)
      }else{
        forecast_dir <- NULL
      }

      inflow_forecast_path <- FLAREr::get_driver_forecast_path(config,
                                                               forecast_model = config$inflow$forecast_inflow_model)

      if(!is.null(inflow_forecast_path)){
        FLAREr::get_driver_forecast(lake_directory, forecast_path = inflow_forecast_path)
        inflow_file_dir <- file.path(config$file_path$noaa_directory,inflow_forecast_path)
      }else{
        inflow_file_dir <- NULL
      }


      config$future_inflow_flow_coeff <- c(0.0010803, 0.9478724, 0.3478991)
      config$future_inflow_flow_error <- 0.00965
      config$future_inflow_temp_coeff <- c(0.20291, 0.94214, 0.04278)
      config$future_inflow_temp_error <- 0.943

      forecast_files <- list.files(file.path(lake_directory, "drivers", noaa_forecast_path), full.names = TRUE)
      temp_flow_forecast <- forecast_inflows_outflows(inflow_obs = file.path(config$file_path$qaqc_data_directory, "fcre-targets-inflow.csv"),
                                                        forecast_files = forecast_files,
                                                        obs_met_file = file.path(config$file_path$qaqc_data_directory,"observed-met_fcre.nc"),
                                                        output_dir = config$file_path$inflow_directory,
                                                        inflow_model = config$inflow$forecast_inflow_model,
                                                        inflow_process_uncertainty = FALSE,
                                                        forecast_location = config$file_path$forecast_output_directory,
                                                        config = config,
                                                        use_s3 = config$run_config$use_s3,
                                                        bucket = "drivers",
                                                        model_name = config$model_settings$model_name)

      #Need to remove the 00 ensemble member because it only goes 16-days in the future

      #pars_config <- NULL #readr::read_csv(file.path(config$file_path$configuration_directory, "FLAREr", config$model_settings$par_config_file), col_types = readr::cols())
      pars_config <- readr::read_csv(file.path(config$file_path$configuration_directory, config$model_settings$par_config_file), col_types = readr::cols())
      obs_config <- readr::read_csv(file.path(config$file_path$configuration_directory, config$model_settings$obs_config_file), col_types = readr::cols())
      states_config <- readr::read_csv(file.path(config$file_path$configuration_directory, config$model_settings$states_config_file), col_types = readr::cols())


      #Download and process observations (already done)

      met_out <- FLAREr::generate_glm_met_files(obs_met_file = file.path(config$file_path$qaqc_data_directory, paste0("observed-met_",config$location$site_id,".nc")),
                                                out_dir = config$file_path$execute_directory,
                                                forecast_dir = forecast_dir,
                                                config = config)

      inflow_outflow_files <- FLAREr::create_glm_inflow_outflow_files(inflow_file_dir = inflow_file_dir,
                                                                      inflow_obs = file.path(config$file_path$qaqc_data_directory, paste0(config$location$site_id, "-targets-inflow.csv")),
                                                                      working_directory = config$file_path$execute_directory,
                                                                      config = config,
                                                                      state_names = states_config$state_names)

      management <- NULL

      if(config$model_settings$model_name == "glm_aed"){

        inflow_outflow_files$inflow_file_name <- cbind(inflow_outflow_files$inflow_file_name, rep(file.path(config$file_path$execute_directory,basename(https_file)), length(inflow_outflow_files$inflow_file_name)))
      }

      #Create observation matrix
      obs <- FLAREr::create_obs_matrix(cleaned_observations_file_long = file.path(config$file_path$qaqc_data_directory,paste0(config$location$site_id, "-targets-insitu.csv")),
                                       obs_config = obs_config,
                                       config)

      #obs[ ,2:dim(obs)[2], ] <- NA

      states_config <- FLAREr::generate_states_to_obs_mapping(states_config, obs_config)

      model_sd <- FLAREr::initiate_model_error(config, states_config)

      init <- FLAREr::generate_initial_conditions(states_config,
                                                  obs_config,
                                                  pars_config,
                                                  obs,
                                                  config,
                                                  restart_file = config$run_config$restart_file,
                                                  historical_met_error = met_out$historical_met_error)
      #Run EnKF
      da_forecast_output <- FLAREr::run_da_forecast(states_init = init$states,
                                                    pars_init = init$pars,
                                                    aux_states_init = init$aux_states_init,
                                                    obs = obs,
                                                    obs_sd = obs_config$obs_sd,
                                                    model_sd = model_sd,
                                                    working_directory = config$file_path$execute_directory,
                                                    met_file_names = met_out$filenames,
                                                    inflow_file_names = inflow_outflow_files$inflow_file_name,
                                                    outflow_file_names = inflow_outflow_files$outflow_file_name,
                                                    config = config,
                                                    pars_config = pars_config,
                                                    states_config = states_config,
                                                    obs_config = obs_config,
                                                    management,
                                                    da_method = config$da_setup$da_method,
                                                    par_fit_method = config$da_setup$par_fit_method)

      # Save forecast

      saved_file <- FLAREr::write_forecast_netcdf(da_forecast_output = da_forecast_output,
                                                  forecast_output_directory = config$file_path$forecast_output_directory,
                                                  use_short_filename = TRUE)

      #Create EML Metadata
      eml_file_name <- FLAREr::create_flare_metadata(file_name = saved_file,
                                                     da_forecast_output = da_forecast_output)

      #rm(da_forecast_output)
      #gc()
      message("Generating plot")
      FLAREr::plotting_general_2(file_name = saved_file,
                                 target_file = file.path(config$file_path$qaqc_data_directory, paste0(config$location$site_id, "-targets-insitu.csv")),
                                 ncore = 2,
                                 obs_csv = FALSE)

      FLAREr::put_forecast(saved_file, eml_file_name, config)

      new_time <- as.character(lubridate::as_datetime(config$run_config$forecast_start_datetime) +
                                                                   lubridate::days(days_between_forecasts))

      FLAREr::update_run_config(config, lake_directory, configure_run_file, saved_file, new_horizon = forecast_horizon, day_advance = days_between_forecasts)
}
