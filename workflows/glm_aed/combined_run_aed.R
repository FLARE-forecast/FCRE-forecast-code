library(tidyverse)
library(lubridate)
set.seed(100)

options(future.globals.maxSize= 891289600)

Sys.setenv("AWS_DEFAULT_REGION" = "renc",
           "AWS_S3_ENDPOINT" = "osn.xsede.org",
           "USE_HTTPS" = TRUE)

lake_directory <- here::here()

files.sources <- list.files(file.path(lake_directory, "R"), full.names = TRUE)
purrr::walk(files.sources, source)

config_set_name <- "glm_aed"
configure_run_file <- "configure_run.yml"

noaa_ready <- TRUE

#while(noaa_ready){

  config <- FLAREr::set_configuration(configure_run_file,lake_directory, config_set_name = config_set_name)
  config <- FLAREr::get_restart_file(config, lake_directory)

  #message(paste0("     Running forecast that starts on: ", config$run_config$start_datetime))

  #Need to remove the 00 ensemble member because it only goes 16-days in the future

  #pars_config <- NULL #readr::read_csv(file.path(config$file_path$configuration_directory, "FLAREr", config$model_settings$par_config_file), col_types = readr::cols())
  pars_config <- NULL #readr::read_csv(file.path(config$file_path$configuration_directory, config$model_settings$par_config_file), col_types = readr::cols())
  obs_config <- readr::read_csv(file.path(config$file_path$configuration_directory, config$model_settings$obs_config_file), col_types = readr::cols())
  states_config <- readr::read_csv(file.path(config$file_path$configuration_directory, config$model_settings$states_config_file), col_types = readr::cols())


  met_start_datetime <- lubridate::as_datetime(config$run_config$start_datetime)
  met_forecast_start_datetime <- lubridate::as_datetime(config$run_config$forecast_start_datetime)

  if(config$run_config$forecast_horizon > 16 & config$met$use_forecasted_met & !config$met$use_openmeteo){
    met_forecast_start_datetime <- met_forecast_start_datetime - lubridate::days(1)
    if(met_forecast_start_datetime < met_start_datetime){
      met_start_datetime <- met_forecast_start_datetime
      message("horizon is > 16 days so adjusting forecast_start_datetime in the met file generation to use yesterdays forecast. But adjusted forecast_start_datetime < start_datetime")
    }
  }




  met_out <- FLAREr::generate_met_files_arrow(out_dir = config$file_path$execute_directory,
                                              start_datetime = met_start_datetime,
                                              end_datetime = config$run_config$end_datetime,
                                              forecast_start_datetime = met_forecast_start_datetime,
                                              forecast_horizon =  config$run_config$forecast_horizon,
                                              site_id = config$location$site_id,
                                              use_s3 = TRUE,
                                              bucket = config$s3$drivers$bucket,
                                              endpoint = config$s3$drivers$endpoint,
                                              use_hive_met = TRUE)

  # met_out <- FLAREr::generate_met_files_openmet(out_dir = config$file_path$execute_directory,
  #                                               start_datetime = lubridate::as_datetime(config$run_config$start_datetime),
  #                                               end_datetime = config$run_config$end_datetime,
  #                                               forecast_start_datetime = lubridate::as_datetime(config$run_config$forecast_start_datetime),
  #                                               forecast_horizon =  config$run_config$forecast_horizon,
  #                                               latitude = config$location$latitude,
  #                                               longitude = config$location$longitude,
  #                                               site_id = config$location$site_id,
  #                                               openmeteo_api = config$met$openmeteo_api,
  #                                               model = config$met$openmeteo_model,
  #                                               use_archive = config$met$use_openmeteo_archive,
  #                                               bucket = config$s3$drivers$bucket,
  #                                               endpoint = config$s3$drivers$endpoint)

  variables <- c("datetime", "FLOW", "TEMP", "SALT",
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
                 'PHY_diatom')


  targets_vera <- readr::read_csv("https://renc.osn.xsede.org/bio230121-bucket01/vera4cast/targets/project_id=vera4cast/duration=P1D/daily-inflow-targets.csv.gz",
                                  show_col_types = FALSE)

  inflow_hist_dates <- tibble(datetime = seq(min(targets_vera$datetime), max(targets_vera$datetime), by = "1 day"))

  targets_vera |>
    filter(!variable %in% c("DN_mgL_sample", "DC_mgL_sample")) |>
    select(datetime, variable, observation) |>
    pivot_wider(names_from = variable, values_from = observation) |>
    right_join(inflow_hist_dates, by = "datetime") |>
    mutate(across(Flow_cms_mean:DIC_mgL_sample, imputeTS::na_interpolation)) |>
    tidyr::fill(Flow_cms_mean:DIC_mgL_sample, .direction = "up") |>
    tidyr::fill(Flow_cms_mean:DIC_mgL_sample, .direction = "down") |>
    dplyr::rename(TEMP = Temp_C_mean,
                  FLOW = Flow_cms_mean) |>
    dplyr::mutate(NIT_amm = NH4_ugL_sample*1000*0.001*(1/18.04),
                  NIT_nit = NO3NO2_ugL_sample*1000*0.001*(1/62.00), #as all NO2 is converted to NO3
                  PHS_frp = SRP_ugL_sample*1000*0.001*(1/94.9714),
                  OGM_doc = DOC_mgL_sample*1000*(1/12.01)* 0.10,  #assuming 10% of total DOC is in labile DOC pool (Wetzel page 753)
                  OGM_docr = 1.5*DOC_mgL_sample*1000*(1/12.01)* 0.90, #assuming 90% of total DOC is in recalcitrant DOC pool
                  TN_ugL = TN_ugL_sample*1000*0.001*(1/14),
                  TP_ugL = TP_ugL_sample*1000*0.001*(1/30.97),
                  OGM_poc = 0.1*(OGM_doc+OGM_docr), #assuming that 10% of DOC is POC (Wetzel page 755
                  OGM_don = (5/6)*(TN_ugL_sample-(NIT_amm+NIT_nit))*0.10, #DON is ~5x greater than PON (Wetzel page 220)
                  OGM_donr = (5/6)*(TN_ugL_sample-(NIT_amm+NIT_nit))*0.90, #to keep mass balance with DOC, DONr is 90% of total DON
                  OGM_pon = (1/6)*(TN_ugL_sample-(NIT_amm+NIT_nit)), #detemined by subtraction
                  OGM_dop = 0.3*(TP_ugL_sample-PHS_frp)*0.10, #Wetzel page 241, 70% of total organic P = particulate organic; 30% = dissolved organic P
                  OGM_dopr = 0.3*(TP_ugL_sample-PHS_frp)*0.90,#to keep mass balance with DOC & DON, DOPr is 90% of total DOP
                  OGM_pop = TP_ugL_sample-(OGM_dop+OGM_dopr+PHS_frp), # #In lieu of having the adsorbed P pool activated in the model, need to have higher complexed P
                  CAR_dic = DIC_mgL_sample*1000*(1/52.515),
                  OXY_oxy = rMR::Eq.Ox.conc(TEMP, elevation.m = 506, #creating OXY_oxy column using RMR package, assuming that oxygen is at 100% saturation in this very well-mixed stream
                                            bar.press = NULL, bar.units = NULL,
                                            out.DO.meas = "mg/L",
                                            salinity = 0, salinity.units = "pp.thou"),
                  OXY_oxy = OXY_oxy *1000*(1/32),
                  CAR_ch4 = CH4_umolL_sample,
                  PHY_cyano = 0,
                  PHY_green = 0,
                  PHY_diatom = 0,
                  SIL_rsi = DRSI_mgL_sample*1000*(1/60.08),
                  SALT = 0) |>
    dplyr::mutate_if(is.numeric, round, 4) |>
    dplyr::select(dplyr::any_of(variables)) |>
    tidyr::pivot_longer(-c("datetime"), names_to = "variable", values_to = "observation") |>
    dplyr::select(datetime, variable, observation) |>
    readr::write_csv(file.path(config$file_path$qaqc_data_directory, paste0(config$location$site_id, "-targets-inflow.csv")))

  inflow_forecast_dir <- "inflow"
  convert_vera4cast_inflow_old(reference_date = lubridate::as_date(config$run_config$forecast_start_datetime),
                           model_id = "inflow_gefsClimAED",
                           save_path = file.path(lake_directory, "drivers", inflow_forecast_dir))

  ## fix datetime object to match the format needed for inflow/outflow function below
  variables[1] <- 'time'

  inflow_outflow_files <- FLAREr::create_inflow_outflow_files_arrow(inflow_forecast_dir = inflow_forecast_dir,
                                                                    inflow_obs = file.path(config$file_path$qaqc_data_directory, paste0(config$location$site_id, "-targets-inflow.csv")),
                                                                    variables = variables,
                                                                    out_dir = config$file_path$execute_directory,
                                                                    start_datetime = config$run_config$start_datetime,
                                                                    end_datetime = config$run_config$end_datetime,
                                                                    forecast_start_datetime = config$run_config$forecast_start_datetime,
                                                                    forecast_horizon =  config$run_config$forecast_horizon,
                                                                    site_id = config$location$site_id,
                                                                    use_s3 = FALSE,
                                                                    bucket = config$s3$inflow_drivers$bucket,
                                                                    endpoint = config$s3$inflow_drivers$endpoint,
                                                                    local_directory = file.path(lake_directory, "drivers", inflow_forecast_dir),
                                                                    use_forecast = TRUE,
                                                                    use_ler_vars = FALSE)

  inflow_outflow_files$inflow_file_names <- inflow_outflow_files$inflow_file_names[which(stringr::str_detect(inflow_outflow_files$inflow_file_names, "ens32", negate = TRUE))]
  inflow_outflow_files$outflow_file_names <- inflow_outflow_files$outflow_file_names[which(stringr::str_detect(inflow_outflow_files$outflow_file_names, "ens32", negate = TRUE))]


  readr::read_csv("https://renc.osn.xsede.org/bio230121-bucket01/vera4cast/targets/project_id=vera4cast/duration=P1D/daily-insitu-targets.csv.gz", show_col_types = FALSE) |>
    dplyr::mutate(observation = ifelse(variable == "DO_mgL_mean", observation*1000*(1/32), observation),
                  observation = ifelse(variable == "fDOM_QSU_mean", -151.3407 + observation*29.62654, observation),
                  depth_m = ifelse(depth_m == 0.1, 0.0, depth_m)) |>
    dplyr::rename(depth = depth_m) |>
    dplyr::filter(site_id == "fcre",
                  datetime >= as_datetime(config$run_config$start_datetime)) |>
    readr::write_csv(file.path(config$file_path$qaqc_data_directory,paste0(config$location$site_id, "-targets-insitu.csv")))


  #Create observation matrix
  obs <- FLAREr::create_obs_matrix(cleaned_observations_file_long = file.path(config$file_path$qaqc_data_directory,paste0(config$location$site_id, "-targets-insitu.csv")),
                                   obs_config = obs_config,
                                   config)

  obs_secchi_depth <- get_obs_secchi_depth(obs_file = file.path(config$file_path$qaqc_data_directory,paste0(config$location$site_id, "-targets-insitu.csv")),
                                           start_datetime = config$run_config$start_datetime,
                                           end_datetime = config$run_config$end_datetime,
                                           forecast_start_datetime = config$run_config$forecast_start_datetime,
                                           forecast_horizon =  config$run_config$forecast_horizon,
                                           secchi_sd = 0.1)

  #obs_secchi_depth <- NULL
  #obs[ ,2:dim(obs)[2], ] <- NA

  states_config <- FLAREr::generate_states_to_obs_mapping(states_config, obs_config)

  model_sd <- FLAREr::initiate_model_error(config, states_config)

  init <- FLAREr::generate_initial_conditions(states_config,
                                              obs_config,
                                              pars_config,
                                              obs,
                                              config,
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
                                                management = NULL,
                                                da_method = config$da_setup$da_method,
                                                par_fit_method = config$da_setup$par_fit_method,
                                                obs_secchi = obs_secchi_depth$obs_secchi,
                                                obs_depth = obs_secchi_depth$obs_depth)

  # Save forecast

  saved_file <- FLAREr::write_forecast_netcdf(da_forecast_output = da_forecast_output,
                                              forecast_output_directory = config$file_path$forecast_output_directory,
                                              use_short_filename = TRUE)

  forecast_df <- FLAREr::write_forecast_arrow(da_forecast_output = da_forecast_output,
                                              use_s3 = config$run_config$use_s3,
                                              bucket = config$s3$forecasts_parquet$bucket,
                                              endpoint = config$s3$forecasts_parquet$endpoint,
                                              local_directory = file.path(lake_directory, "forecasts/parquet"))

  vera_variables <- c("Temp_C_mean","Chla_ugL_mean", "DO_mgL_mean", "fDOM_QSU_mean", "NH4_ugL_sample",
                      "NO3NO2_ugL_sample", "SRP_ugL_sample", "DIC_mgL_sample","Secchi_m_sample",
                      "Bloom_binary_mean","CH4_umolL_sample","IceCover_binary_max")

  # Calculate probablity of bloom
  bloom_binary <- forecast_df |>
    dplyr::filter(depth == 1.6 & variable == "Chla_ugL_mean") |>
    dplyr::mutate(over = ifelse(prediction > 20, 1, 0)) |>
    dplyr::summarize(prediction = sum(over) / n(), .by = c(datetime, reference_datetime, model_id, site_id, depth, variable)) |> #pubDate
    dplyr::mutate(family = "bernoulli",
                  parameter = "prob",
                  variable = "Bloom_binary_mean") |>
    dplyr::rename(depth_m = depth) |>
    dplyr::select(reference_datetime, datetime, model_id, site_id, depth_m, family, parameter, variable, prediction)

  # Calculate probablity of having ice
  ice_binary <- forecast_df |>
    dplyr::filter(variable == "ice_thickness") |>
    dplyr::mutate(over = ifelse(prediction > 0, 1, 0)) |>
    dplyr::summarize(prediction = sum(over) / n(), .by = c(datetime, reference_datetime, model_id, site_id, depth, variable)) |> #pubDate
    dplyr::mutate(family = "bernoulli",
                  parameter = "prob",
                  variable = "IceCover_binary_max",
                  depth = NA) |>
    dplyr::rename(depth_m = depth) |>
    dplyr::select(reference_datetime, datetime, model_id, site_id, depth_m, family, parameter, variable, prediction)

  # Calculate probablity of being mixed
  min_depth <- 1
  max_depth <- 8
  threshold <- 0.1

  temp_forecast <- forecast_df |>
    filter(depth %in% c(min_depth, max_depth),
           variable == "Temp_C_mean") |>
    pivot_wider(names_from = depth, names_prefix = 'wtr_', values_from = prediction)

  colnames(temp_forecast)[which(colnames(temp_forecast) == paste0('wtr_', min_depth))] <- 'min_depth'
  colnames(temp_forecast)[which(colnames(temp_forecast) == paste0('wtr_', max_depth))] <- 'max_depth'

  mix_binary <- temp_forecast |>
    mutate(min_depth = rLakeAnalyzer::water.density(min_depth),
           max_depth = rLakeAnalyzer::water.density(max_depth),
           mixed = ifelse((max_depth - min_depth) < threshold, 1, 0)) |>
    summarise(prediction = 100*(sum(mixed)/n()), .by = c(datetime, reference_datetime, model_id, site_id, variable)) |> #pubDate
    dplyr::mutate(family = "bernoulli",
                  parameter = "prob",
                  variable = "Mixed_binary_sample",
                  depth = NA) |>
    dplyr::rename(depth_m = depth) |>
    dplyr::select(reference_datetime, datetime, model_id, site_id, depth_m, family, parameter, variable, prediction)


  # Combine into a vera data frame
  vera4cast_df <- forecast_df |>
    dplyr::rename(depth_m = depth) |>
    dplyr::mutate(prediction = ifelse(variable == "DO_mgL_mean", prediction/1000*(32),prediction),
                  prediction = ifelse(variable == "fDOM_QSU_mean", (151.3407 + prediction)/29.62654,prediction),
                  prediction = ifelse(variable == "NIT_amm", prediction/1000/0.001/(1/18.04),prediction),
                  variable = ifelse(variable == "NIT_amm", "NH4_ugL_sample", variable),
                  prediction = ifelse(variable == "NIT_nit", prediction/1000/0.001/(1/62.00),prediction),
                  variable = ifelse(variable == "NIT_amm", "NO3NO2_ugL_sample", variable),
                  prediction = ifelse(variable == "PHS_frp", prediction/1000/0.001/(1/94.9714),prediction),
                  variable = ifelse(variable == "PHS_frp", "SRP_ugL_sample", variable),
                  prediction = ifelse(variable == "CAR_dic", prediction/1000/(1/52.515), prediction),
                  variable = ifelse(variable == "CAR_dic", "DIC_mgL_sample", variable),
                  variable = ifelse(variable == "CAR_ch4", "CH4_umolL_sample", variable),
                  variable = ifelse(variable == "secchi", "Secchi_m_sample", variable),
                  depth_m = ifelse(depth_m == 0.0, 0.1, depth_m)) |>
    dplyr::select(-forecast, -variable_type) |> #pubDate
    dplyr::mutate(parameter = as.character(parameter)) |>
    dplyr::bind_rows(bloom_binary) |>
    dplyr::bind_rows(ice_binary) |>
    dplyr::bind_rows(mix_binary) |>
    dplyr::filter(variable %in% vera_variables) |>
    mutate(project_id = "vera4cast",
           model_id = config$run_config$sim_name,
           family = "ensemble",
           site_id = "fcre",
           duration = "P1D",
           datetime = lubridate::as_datetime(datetime),
           reference_datetime = lubridate::as_datetime(reference_datetime)) |>
    filter(datetime >= reference_datetime)

  #vera4cast_df |>
  #  filter(depth_m == 1.6 | is.na(depth_m)) |>
  #  ggplot(aes(x = datetime, y = prediction, group = factor(parameter))) +
  #  geom_line() +
  #  facet_wrap(~variable, scale = "free")

  file_name <- paste0(config$run_config$sim_name,
                      "-",
                      lubridate::as_date(vera4cast_df$reference_datetime[1]),".csv.gz")

  readr::write_csv(vera4cast_df, file = file_name)



  #FLAREr::generate_forecast_score_arrow(targets_file = file.path(config$file_path$qaqc_data_directory,paste0(config$location$site_id, "-targets-insitu.csv")),
  #                                      forecast_df = forecast_df,
  #                                      use_s3 = config$run_config$use_s3,
  #                                      bucket = config$s3$scores$bucket,
  #                                      endpoint = config$s3$scores$endpoint,
  #                                      local_directory = file.path(lake_directory, "scores/parquet"))

  #message("Generating plot")
  #FLAREr::plotting_general_2(file_name = saved_file,
  #                           target_file = file.path(config$file_path$qaqc_data_directory, paste0(config$location$site_id, "-targets-insitu.csv")),
  #                           ncore = 2,
  #                           obs_csv = FALSE)

  FLAREr::put_forecast(saved_file, eml_file_name = NULL, config)

  forecast_start_datetime <- lubridate::as_datetime(config$run_config$forecast_start_datetime) + lubridate::days(1)
  start_datetime <- lubridate::as_datetime(config$run_config$forecast_start_datetime)
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

  var1 <- Sys.getenv("AWS_ACCESS_KEY_ID")
  var2 <- Sys.getenv("AWS_SECRET_ACCESS_KEY")
  Sys.unsetenv("AWS_ACCESS_KEY_ID")
  Sys.unsetenv("AWS_SECRET_ACCESS_KEY")
  vera4castHelpers::submit(file_name, first_submission = FALSE)

  Sys.setenv("AWS_ACCESS_KEY_ID" = var1,
             "AWS_SECRET_ACCESS_KEY" = var2)

 # noaa_ready <- FLAREr::check_noaa_present_arrow(lake_directory,
 #                                                configure_run_file,
 #                                                config_set_name = config_set_name)

#}




