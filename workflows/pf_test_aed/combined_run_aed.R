library(tidyverse)
library(lubridate)
set.seed(201)

sim_name <- "rol_july19_pf"
config_set_name <- "pf_test_aed"
configure_run_file <- "configure_aed_run.yml"

options(future.globals.maxSize = 891289600)

lake_directory <- here::here()

use_s3 <- FALSE

walk(list.files(file.path(lake_directory, "R"), full.names = TRUE), source)

unlink(file.path(lake_directory, "restart", "fcre", sim_name, configure_run_file))

yml <- yaml::read_yaml(file.path(lake_directory, "configuration", config_set_name, configure_run_file))
yml$sim_name <- sim_name
yaml::write_yaml(yml, file.path(lake_directory, "configuration", config_set_name, configure_run_file))

site <- "fcre"

# Set up configurations for the data processing
config_obs <- FLAREr::initialize_obs_processing(lake_directory, observation_yml = "observation_processing.yml", config_set_name = config_set_name)

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

FLAREr::get_edi_file(edi_https = "https://pasta.lternet.edu/package/data/eml/edi/389/7/02d36541de9088f2dd99d79dc3a7a853",
                     file = config_obs$met_raw_obs_fname[2],
                     lake_directory)

FLAREr::get_edi_file(edi_https = "https://pasta.lternet.edu/package/data/eml/edi/271/7/71e6b946b751aa1b966ab5653b01077f",
                     file = config_obs$insitu_obs_fname[2],
                     lake_directory)

FLAREr::get_edi_file(edi_https = "https://pasta.lternet.edu/package/data/eml/edi/198/11/81f396b3e910d3359907b7264e689052",
                     file = config_obs$secchi_fname,
                     lake_directory)

FLAREr::get_edi_file(edi_https = "https://pasta.lternet.edu/package/data/eml/edi/200/13/27ceda6bc7fdec2e7d79a6e4fe16ffdf",
                     file = config_obs$ctd_fname,
                     lake_directory)

FLAREr::get_edi_file(edi_https = "https://pasta.lternet.edu/package/data/eml/edi/199/11/509f39850b6f95628d10889d66885b76",
                     file = config_obs$nutrients_fname,
                     lake_directory)


FLAREr::get_edi_file(edi_https = "https://pasta.lternet.edu/package/data/eml/edi/202/9/c065ff822e73c747f378efe47f5af12b",
                     file = config_obs$inflow_raw_file1[2],
                     lake_directory)

FLAREr::get_edi_file(edi_https = "https://pasta.lternet.edu/package/data/eml/edi/542/1/791ec9ca0f1cb9361fa6a03fae8dfc95",
                     file = "silica_master_df.csv",
                     lake_directory)

FLAREr::get_edi_file(edi_https = "https://pasta.lternet.edu/package/data/eml/edi/551/7/38d72673295864956cccd6bbba99a1a3",
                     file = "Dissolved_CO2_CH4_Virginia_Reservoirs.csv",
                     lake_directory)

#' Clean up observed meteorology

cleaned_met_file <- met_data_bind(realtime_file = file.path(config_obs$file_path$data_directory, config_obs$met_raw_obs_fname[1]),
                                  qaqc_file = file.path(config_obs$file_path$data_directory, config_obs$met_raw_obs_fname[2]),
                                  cleaned_met_file = file.path(config_obs$file_path$targets_directory, config_obs$site_id,paste0("observed-met_",config_obs$site_id,".csv")),
                                  input_file_tz = "EST",
                                  nldas = NULL,
                                  site_id = config_obs$site_id)

#' Clean up observed inflow

cleaned_inflow_file <- inflow_data_combine(realtime_file = file.path(config_obs$file_path$data_directory, config_obs$inflow_raw_file1[1]),
                                           qaqc_file = file.path(config_obs$file_path$data_directory, config_obs$inflow_raw_file1[2]),
                                           nutrients_file = file.path(config_obs$file_path$data_directory, config_obs$nutrients_fname),
                                           silica_file = file.path(config_obs$file_path$data_directory,  config_obs$silica_fname),
                                           co2_ch4 = file.path(config_obs$file_path$data_directory, config_obs$ch4_fname),
                                           cleaned_inflow_file = file.path(config_obs$file_path$targets_directory, config_obs$site_id, paste0(config_obs$site_id,"-targets-inflow.csv")),
                                           input_file_tz = 'EST',
                                           site_id = config_obs$site_id)






#read_csv(cleaned_inflow_file) |>
#  mutate(observation = ifelse(variable == "PHS_frp", observation * 100, observation)) |>
#  write_csv(cleaned_inflow_file)

#' Clean up observed insitu measurements

cleaned_insitu_file <- in_situ_qaqc_csv(insitu_obs_fname = file.path(config_obs$file_path$data_directory,config_obs$insitu_obs_fname),
                                        data_location = config_obs$file_path$data_directory,
                                        maintenance_file = file.path(config_obs$file_path$data_directory,config_obs$maintenance_file),
                                        ctd_fname = file.path(config_obs$file_path$data_directory, config_obs$ctd_fname),
                                        nutrients_fname =  file.path(config_obs$file_path$data_directory, config_obs$nutrients_fname),
                                        secchi_fname = file.path(config_obs$file_path$data_directory, config_obs$secchi_fname),
                                        ch4_fname = file.path(config_obs$file_path$data_directory, config_obs$ch4_fname),
                                        cleaned_insitu_file = file.path(config_obs$file_path$targets_directory, config_obs$site_id, paste0(config_obs$site_id,"-targets-insitu.csv")),
                                        lake_name_code = config_obs$site_id,
                                        config = config_obs)

d <- read_csv(cleaned_insitu_file, show_col_types = FALSE) |>
  filter((variable %in% c("NH4", "NO3NO2", "SRP", "TN", "TP") & depth == 1.5) |
           variable %in% c("temperature", "oxygen", "secchi", "chla", "fdom", "depth"),
          variable != "depth" & datetime == as_datetime("2021-01-01 00:00:00") | variable == "depth") |>
  write_csv(cleaned_insitu_file)

config <- FLAREr::set_configuration(configure_run_file, lake_directory, config_set_name = config_set_name)

file.copy(file.path(config$file_path$configuration_directory,"FCR_SSS_inflow_2013_2021_20220413_allfractions_2DOCpools.csv"),
          file.path(config$file_path$execute_directory,"FCR_SSS_inflow_2013_2021_20220413_allfractions_2DOCpools.csv"))

config <- FLAREr::set_configuration(configure_run_file = configure_run_file ,lake_directory = lake_directory, config_set_name = config_set_name)
config <- FLAREr::get_restart_file(config, lake_directory)

pars_config <- readr::read_csv(file.path(config$file_path$configuration_directory, config$model_settings$par_config_file), col_types = readr::cols())
obs_config <- readr::read_csv(file.path(config$file_path$configuration_directory, config$model_settings$obs_config_file), col_types = readr::cols())
states_config <- readr::read_csv(file.path(config$file_path$configuration_directory, config$model_settings$states_config_file), col_types = readr::cols())


# Inflows
hist_interp_inflow <- interpolate_targets(targets = paste0(config$location$site_id,"-targets-inflow.csv"),
                                          lake_directory = lake_directory,
                                          targets_dir = 'targets',
                                          site_id = config$location$site_id,
                                          #variables = c('FLOW', 'SALT', 'TEMP'),
                                          variables <- c("time", "FLOW", "TEMP", "SALT",
                                                         'OXY_oxy',
                                                         'CAR_dic',
                                                         'CAR_ch4',
                                                         'SIL_rsi',
                                                         'NIT_amm',
                                                         'NIT_nit',
                                                         'PHS_frp',
                                                         'OGM_doc',
                                                         'OGM_docr',
                                                         'OGM_poc',
                                                         'OGM_don',
                                                         'OGM_donr',
                                                         'OGM_pon',
                                                         'OGM_dop',
                                                         'OGM_dopr',
                                                         'OGM_pop',
                                                         'PHY_cyano',
                                                         'PHY_green',
                                                         'PHY_diatom'),
                                          depth = F,
                                          method = 'linear')

hist_interp_inflow <- hist_interp_inflow |>
  mutate(flow_number = 1,
       parameter = 1) |>
  rename(prediction = observation)

# Write the interpolated data as the historal file
arrow::write_dataset(hist_interp_inflow,
                     'drivers/inflow/historical/historical_interp/site_id=fcre/')

# generate a simple "forecast" that has ensemble members
forecast_date <- config$run_config$forecast_start_datetime
future_inflow <- hist_interp_inflow |>
  filter(datetime > forecast_date,
         datetime <= as_date(forecast_date) + config$run_config$forecast_horizon) |>
  mutate(parameter = 1,
         flow_number = 1) |>
  #reframe(prediction = rnorm(10, mean = observation, sd = 1),
  #        parameter = 1:10,
  #        .by = c(site_id, datetime, variable)) |>
  mutate(reference_datetime = as_date(forecast_date))

arrow::write_dataset(future_inflow,
                     'drivers/inflow/future/historical_interp/',
                     partitioning = c('reference_datetime', 'site_id'))

#==========================================#

# outflows
hist_interp_outflow <- interpolate_targets(targets = paste0(config$location$site_id,"-targets-inflow.csv"),
                                           lake_directory = lake_directory,
                                           targets_dir = 'targets',
                                           site_id = config$location$site_id,
                                           #variables = c('FLOW', 'SALT', 'TEMP'),
                                           variables <- c("time", "FLOW"),
                                           depth = F,
                                           method = 'linear')

hist_interp_outflow <- hist_interp_outflow |>
  mutate(flow_number = 1,
         parameter = 1) |>
  rename(prediction = observation)

# Write the interpolated data as the historal file
arrow::write_dataset(hist_interp_outflow,
                     'drivers/outflow/historical/historical_interp/site_id=fcre/')

# generate a simple "forecast" that has ensemble members
future_outflow <- hist_interp_outflow |>
  filter(datetime > forecast_date,
         datetime <= as_date(forecast_date) + config$run_config$forecast_horizon) |>
  mutate(parameter = 1,
         flow_number = 1) |>
  #reframe(prediction = rnorm(10, mean = observation, sd = 1),
  #        parameter = 1:10,
  #        .by = c(site_id, datetime, variable)) |>
  mutate(reference_datetime = as_date(forecast_date))

arrow::write_dataset(future_outflow,
                     'drivers/outflow/future/historical_interp/',
                     partitioning = c('reference_datetime', 'site_id'))



met_out <- FLAREr::generate_met_files_arrow(obs_met_file = file.path(config$file_path$qaqc_data_directory, paste0("observed-met_",config$location$site_id,".csv")),
                                            out_dir = config$file_path$execute_directory,
                                            start_datetime = config$run_config$start_datetime,
                                            end_datetime = config$run_config$end_datetime,
                                            forecast_start_datetime = config$run_config$forecast_start_datetime,
                                            forecast_horizon =  config$run_config$forecast_horizon,
                                            site_id = config$location$site_id,
                                            use_s3 = TRUE,
                                            bucket = config$s3$drivers$bucket,
                                            endpoint = config$s3$drivers$endpoint,
                                            local_directory = NULL,
                                            use_forecast = config$met$use_forecasted_met)

# if(config$model_settings$model_name == "glm_aed"){
#   variables <- c("time", "FLOW", "TEMP", "SALT",
#                  'OXY_oxy',
#                  'CAR_dic',
#                  'CAR_ch4',
#                  'SIL_rsi',
#                  'NIT_amm',
#                  'NIT_nit',
#                  'PHS_frp',
#                  'OGM_doc',
#                  'OGM_docr',
#                  'OGM_poc',
#                  'OGM_don',
#                  'OGM_donr',
#                  'OGM_pon',
#                  'OGM_dop',
#                  'OGM_dopr',
#                  'OGM_pop',
#                  'PHY_cyano',
#                  'PHY_green',
#                  'PHY_diatom')
# }else{
#   variables <- c("time", "FLOW", "TEMP", "SALT")
# }
#
# if(config$run_config$forecast_horizon > 0){
#   inflow_forecast_dir = file.path(config$inflow$forecast_inflow_model, config$location$site_id, "0", lubridate::as_date(config$run_config$forecast_start_datetime))
# }else{
#   inflow_forecast_dir <- NULL
# }

# inflow_outflow_files <- FLAREr::create_inflow_outflow_files_arrow(inflow_forecast_dir = NULL,
#                                                                   inflow_obs = file.path(config$file_path$qaqc_data_directory, paste0(config$location$site_id, "-targets-inflow.csv")),
#                                                                   variables = variables,
#                                                                   out_dir = config$file_path$execute_directory,
#                                                                   start_datetime = config$run_config$start_datetime,
#                                                                   end_datetime = config$run_config$end_datetime,
#                                                                   forecast_start_datetime = config$run_config$forecast_start_datetime,
#                                                                   forecast_horizon =  config$run_config$forecast_horizon,
#                                                                   site_id = config$location$site_id,
#                                                                   use_s3 = use_s3,
#                                                                   bucket = config$s3$inflow_drivers$bucket,
#                                                                   endpoint = config$s3$inflow_drivers$endpoint,
#                                                                   local_directory = file.path(lake_directory, "drivers", inflow_forecast_dir),
#                                                                   use_forecast = FALSE,
#                                                                   use_ler_vars = FALSE)


inflow_outflow_files <- FLAREr::create_inflow_outflow_files_arrow(config, config_set_name)


if(config$model_settings$model_name == "glm_aed"){
  inflow_outflow_files$inflow_file_names <- cbind(inflow_outflow_files$inflow_file_names, rep(file.path(config$file_path$execute_directory,"FCR_SSS_inflow_2013_2021_20220413_allfractions_2DOCpools.csv"), length(inflow_outflow_files$inflow_file_name)))
}

#Create observation matrix
obs <- FLAREr::create_obs_matrix(cleaned_observations_file_long = file.path(config$file_path$qaqc_data_directory,paste0(config$location$site_id, "-targets-insitu.csv")),
                                 obs_config = obs_config,
                                 config)

obs_non_vertical <- FLAREr::create_obs_non_vertical(cleaned_observations_file_long = file.path(config$file_path$qaqc_data_directory,paste0(config$location$site_id, "-targets-insitu.csv")),
                                                    obs_config,
                                                    start_datetime = config$run_config$start_datetime,
                                                    end_datetime = config$run_config$end_datetime,
                                                    forecast_start_datetime = config$run_config$forecast_start_datetime,
                                                    forecast_horizon =  config$run_config$forecast_horizon)


states_config <- FLAREr::generate_states_to_obs_mapping(states_config, obs_config)

model_sd <- FLAREr::initiate_model_error(config, states_config)

init <- FLAREr::generate_initial_conditions(states_config,
                                            obs_config,
                                            pars_config,
                                            obs,
                                            config,
                                            obs_non_vertical = obs_non_vertical)

states_init = init$states
pars_init = init$pars
aux_states_init = init$aux_states_init
obs = obs
obs_sd = obs_config$obs_sd
model_sd = model_sd
working_directory = config$file_path$execute_directory
met_file_names = met_out$filenames
inflow_file_names = inflow_outflow_files$inflow_file_names[,1]
outflow_file_names = inflow_outflow_files$outflow_file_names
config = config
pars_config = pars_config
states_config = states_config
obs_config = obs_config
management = NULL
da_method = config$da_setup$da_method
par_fit_method = config$da_setup$par_fit_method
obs_secchi = obs_non_vertical$obs_secchi
obs_depth = obs_non_vertical$obs_depth
debug = FALSE
log_wq = FALSE




da_forecast_output <- FLAREr::run_da_forecast(states_init = init$states,
                                              pars_init = init$pars,
                                              aux_states_init = init$aux_states_init,
                                              obs = obs,
                                              obs_sd = obs_config$obs_sd,
                                              model_sd = model_sd,
                                              working_directory = config$file_path$execute_directory,
                                              met_file_names = met_out$filenames,
                                              inflow_file_names = inflow_outflow_files$inflow_file_names[,1],
                                              outflow_file_names = inflow_outflow_files$outflow_file_names,
                                              config = config,
                                              pars_config = pars_config,
                                              states_config = states_config,
                                              obs_config = obs_config,
                                              da_method = config$da_setup$da_method,
                                              par_fit_method = config$da_setup$par_fit_method,
                                              obs_secchi = obs_non_vertical$obs_secchi,
                                              obs_depth = obs_non_vertical$obs_depth)

# Save forecast

saved_file <- FLAREr::write_forecast_netcdf(da_forecast_output = da_forecast_output,
                                            forecast_output_directory = config$file_path$forecast_output_directory,
                                            use_short_filename = TRUE)

forecast_df <- FLAREr::write_forecast_arrow(da_forecast_output = da_forecast_output,
                                            use_s3 = use_s3,
                                            bucket = config$s3$forecasts_parquet$bucket,
                                            endpoint = config$s3$forecasts_parquet$endpoint,
                                            local_directory = file.path(lake_directory, "forecasts/parquet"))



targets_df <- read_csv(file.path(config$file_path$qaqc_data_directory,paste0(config$location$site_id, "-targets-insitu.csv")))



FLAREr::plotting_general(forecast_df, targets_df, file_name = saved_file)


rm(da_forecast_output)
gc()

FLAREr::generate_forecast_score_arrow(targets_file = file.path(config$file_path$qaqc_data_directory,paste0(config$location$site_id, "-targets-insitu.csv")),
                                      forecast_df = forecast_df,
                                      use_s3 = use_s3,
                                      bucket = config$s3$scores$bucket,
                                      endpoint = config$s3$scores$endpoint,
                                      local_directory = file.path(lake_directory, "scores/parquet"),
                                      variable_type = c("state","parameter","diagnostic"))

rm(forecast_df)
gc()
