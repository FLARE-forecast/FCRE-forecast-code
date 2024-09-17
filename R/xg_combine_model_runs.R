## function for creating inflow forecasts using xgboost model

xg_combine_model_runs <- function(site_id,
                             forecast_start_datetime,
                             use_s3_inflow = FALSE,
                             inflow_bucket = NULL,
                             inflow_endpoint = NULL,
                             inflow_local_directory = NULL,
                             forecast_horizon = NULL,
                             inflow_model = NULL){

  if((!is.null(forecast_start_datetime)) && (forecast_horizon > 0)){

    #print('inside of run_inflow_model')

    forecast_date <- lubridate::as_date(forecast_start_datetime) - lubridate::days(1)
    forecast_hour <- lubridate::hour(forecast_start_datetime)

    ## pull in observed met data
    observed_airtemp <- read_csv(file.path(lake_directory,'targets', config_obs$site_id,
                                         paste0('observed-met_',config_obs$site_id,'.csv')), show_col_types = FALSE) |>
      filter(variable == 'air_temperature') |>
      mutate(date = as.Date(datetime),
             temp = observation - 273.15) |>
      summarise(daily_temperature = median(temp, na.rm = TRUE), .by = c("date"))

    ## pull in future NOAA data

    met_s3_future <- arrow::s3_bucket(file.path("bio230121-bucket01/flare/drivers/met/gefs-v12/stage2",paste0("reference_datetime=",noaa_date),paste0("site_id=",site_id)),
                                      endpoint_override = endpoint,
                                      anonymous = TRUE)
    # old method
    # met_s3_future <- arrow::s3_bucket(file.path(config$s3$drivers$bucket,
    #                                             paste0("stage2/parquet/0/", forecast_date),
    #                                             paste0("", site_id)),
    #                                   endpoint_override = endpoint,
    #                                   anonymous = TRUE)

    df_future <- arrow::open_dataset(met_s3_future) |>
      dplyr::filter(variable %in% c("precipitation_flux","air_temperature")) |>
      collect() |>
      rename(ensemble = parameter) |>
      mutate(variable = ifelse(variable == "precipitation_flux", "precipitation", variable),
             variable = ifelse(variable == "air_temperature", "temperature_2m", variable),
             prediction = ifelse(variable == "temperature_2m", prediction - 273.15, prediction))

    min_datetime <- min(df_future$datetime)


    ## pull in past NOAA data
    met_s3_past <- arrow::s3_bucket(paste0("bio230121-bucket01/flare/drivers/met/gefs-v12/stage3/site_id=",site_id),
                                    endpoint_override = endpoint,
                                    anonymous = TRUE)

    years_prior <- forecast_start_datetime - lubridate::days(1825) # 5 years

    df_past <- arrow::open_dataset(met_s3_past) |>
      #select(datetime, parameter, variable, prediction) |>
      dplyr::filter(variable %in% c("precipitation_flux","air_temperature"),
             ((datetime <= min_datetime  & variable == "precipitation_flux") |
                datetime < min_datetime  & variable == "air_temperature"),
             datetime > years_prior) |>
      collect() |>
      rename(ensemble = parameter) |>
      mutate(variable = ifelse(variable == "precipitation_flux", "precipitation", variable),
             variable = ifelse(variable == "air_temperature", "temperature_2m", variable),
             prediction = ifelse(variable == "temperature_2m", prediction - 273.15, prediction)) |>
      select(-reference_datetime)


    # combine past and future noaa data
    df_combined <- bind_rows(df_future, df_past) |>
      arrange(variable, datetime, ensemble) #|>
      # mutate(date = as.Date(datetime)) |>
      # right_join(observed_airtemp, by = c('date'))

    forecast_precip <- df_combined |>
      dplyr::filter(variable == 'precipitation') |>
      summarise(precip_hourly = median(prediction, na.rm = TRUE), .by = c("datetime")) |> # get the median hourly precip across all EMs
      mutate(date = lubridate::as_date(datetime)) |>
      summarise(precip = sum(precip_hourly, na.rm = TRUE), .by = c("date")) |> # get the total precip for each day
      mutate(sevenday_precip = RcppRoll::roll_sum(precip, n = 7, fill = NA,align = "right")) |>
      mutate(doy = lubridate::yday(date))

    forecast_temp <- df_combined |>
      dplyr::filter(variable == 'temperature_2m') |>
      summarise(temp_hourly = median(prediction, na.rm = TRUE), .by = c("datetime")) |> # get the median hourly temp across all EMs
      mutate(date = lubridate::as_date(datetime)) |>
      summarise(temperature = median(temp_hourly, na.rm = TRUE), .by = c("date")) |> # get median temp across hours of the day
      mutate(threeday_temp = RcppRoll::roll_sum(temperature, n = 3, fill = NA,align = "right")) |>
      right_join(observed_airtemp, by = c('date'))


    forecast_met <- forecast_precip |>
      right_join(forecast_temp, by = c('date'))


    print('done setting up met data')

    ## RUN PREDICTIONS
    #sensorcode_df <- read_csv('configuration/default/sensorcode.csv', show_col_types = FALSE)
    inflow_targets <- read_csv(file.path(lake_directory,'targets', config_obs$site_id,
                                         paste0(config_obs$site_id,"-targets-inflow.csv")), show_col_types = FALSE)
    ## RUN FLOW PREDICTIONS
    print('Running Flow Inflow Forecast')

    flow_targets <- inflow_targets |>
      dplyr::filter(variable == 'FLOW') |>
      rename(date = datetime)

    flow_drivers <- forecast_met |>
      left_join(flow_targets, by = c('date')) |>
      drop_na(observation)

    flow_training_df <- flow_drivers |>
      dplyr::filter(date < reference_datetime)

    flow_rec <- recipe(observation ~ precip + sevenday_precip + doy + temperature, # daily_temperature, threeday_temp
                  data = flow_training_df)

    ### RUN THIS FUNCTION TO FIND ERROR
    flow_predictions <- xg_run_inflow_model(train_data = flow_training_df,
                                             model_recipe = flow_rec,
                                             met_combined = df_combined,
                                             targets_df = flow_targets,
                                             drivers_df = flow_drivers,
                                             var_name = 'FLOW')

    ## RUN TEMPERATURE PREDICTIONS
    print('Running Temperature Inflow Forecast')

    temp_targets <- inflow_targets |>
      dplyr::filter(variable == 'TEMP') |>
      rename(date = datetime)

    temp_drivers <- forecast_met |>
      left_join(temp_targets, by = c('date')) |>
      drop_na(observation)

    temp_training_df <- temp_drivers |>
      dplyr::filter(date < reference_datetime)

    temp_rec <- recipe(observation ~ doy + threeday_temp + temperature,
                       data = temp_training_df)

    temp_predictions <- xg_run_inflow_model(train_data = temp_training_df,
                                             model_recipe = temp_rec,
                                             met_combined = df_combined,
                                             targets_df = temp_targets,
                                             drivers_df = temp_drivers,
                                             var_name = 'TEMP')

    ## RUN SALINITY PREDICTIONS
    print('Running Salinity Inflow Forecast')

    salt_targets <- inflow_targets |>
      dplyr::filter(variable == 'SALT') |>
      rename(date = datetime)

    salt_drivers <- forecast_met |>
      left_join(salt_targets, by = c('date')) |>
      drop_na(observation)

    salt_training_df <- salt_drivers |>
      dplyr::filter(date < reference_datetime)

    salt_rec <- recipe(observation ~ doy + threeday_temp,
                       data = salt_training_df)

    salt_predictions <- xg_run_inflow_model(train_data = salt_training_df,
                                             model_recipe = salt_rec,
                                             met_combined = df_combined,
                                             targets_df = salt_targets,
                                             drivers_df = salt_drivers,
                                             var_name = 'SALT')


    ## ALL OTHER VARIABLES NEEDED FOR AED
    inflow_variables <- inflow_targets |>
      filter(!(variable %in% c('TEMP','SALT','FLOW'))) |>
      distinct(variable) |>
      pull(variable)

    var_prediction_build <- data.frame()

    for (i in inflow_variables){
      print(i)

      var_targets <- inflow_targets |>
        dplyr::filter(variable == i) |>
        rename(date = datetime)

      var_drivers <- forecast_met |>
        left_join(var_targets, by = c('date')) |>
        drop_na(observation)

      var_training_df <- var_drivers |>
        dplyr::filter(date < reference_datetime)

      var_rec <- recipe(observation ~ doy,
                         data = var_training_df)

      var_predictions <- xg_run_inflow_model(train_data = var_training_df,
                                              model_recipe = var_rec,
                                              met_combined = df_combined,
                                              targets_df = var_targets,
                                              drivers_df = var_drivers,
                                              var_name = i)

      var_prediction_build <- dplyr::bind_rows(var_prediction_build, var_predictions)
    }


    ## COMBINE ALL INFLOW PREDICTIONS

    inflow_combined <- bind_rows(flow_predictions, temp_predictions, salt_predictions, var_prediction_build) |>
      mutate(reference_datetime = as.Date(reference_datetime))

    outflow_df <- inflow_combined
    outflow_df$flow_type <- 'outflow'

    ## CREATE AND SAVE INFLOWS
    historic_inflow <- inflow_combined |>
      filter(datetime < forecast_date)

    arrow::write_dataset(historic_inflow, path = file.path(lake_directory,
                                                               "drivers/inflow/historic",paste0('model_id=',config$flows$forecast_inflow_model),"site_id=fcre"))

    future_inflow <- inflow_combined |>
      filter(datetime > forecast_date,
             datetime <= as_date(forecast_date) + config$run_config$forecast_horizon)

    arrow::write_dataset(future_inflow,
                         file.path(lake_directory, "drivers/inflow/future",paste0('model_id=',config$flows$forecast_inflow_model)),
                         partitioning = c('reference_datetime', 'site_id'))

    ## CREATE AND SAVE OUTFLOWS

    historic_outflow <- outflow_df |>
      filter(datetime < forecast_date)

    arrow::write_dataset(historic_outflow, path = file.path(lake_directory,
                                                            "drivers/outflow/historic",paste0('model_id=',config$flows$forecast_inflow_model),"site_id=fcre"))

    future_outflow <- outflow_df |>
      filter(datetime > forecast_date,
             datetime <= as_date(forecast_date) + config$run_config$forecast_horizon)

    arrow::write_dataset(future_outflow,
                         file.path(lake_directory, "drivers/outflow/future",paste0('model_id=',config$flows$forecast_inflow_model)),
                         partitioning = c('reference_datetime', 'site_id'))

    #flow_combined <- bind_rows(inflow_combined, outflow_df)

    #flow_combined <- inflow_combined

    # if (forecast_horizon > 0) {
    #   inflow_forecast_path <- file.path(inflow_model, site_id, forecast_hour, lubridate::as_date(forecast_start_datetime))
    # }else {
    #   inflow_forecast_path <- NULL
    # }

    # if(use_s3_inflow){
    #   #FLAREr:::arrow_env_vars()
    #   inflow_s3 <- arrow::s3_bucket(bucket = file.path(inflow_bucket, inflow_forecast_path), endpoint_override = inflow_endpoint)
    #   #on.exit(FLAREr:::unset_arrow_vars(vars))
    # }else{
    #   inflow_s3 <- arrow::SubTreeFileSystem$create(file.path(inflow_local_directory, inflow_forecast_path))
    # }
    # inflow_s3 <- arrow::SubTreeFileSystem$create(file.path(inflow_local_directory, inflow_forecast_path))
    # arrow::write_dataset(flow_combined, path = inflow_s3)
    #
    inflow_local_files <- list.files(file.path(inflow_local_directory), full.names = TRUE, recursive = TRUE)

    print("Inflow files processed...")


  }else{
    message("nothing to forecast")
    inflow_local_files <- NULL
  }# end if statement

  return(inflow_local_files)
}
