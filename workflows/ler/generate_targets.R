
# read in QAQC data
# this function is from the VERA targets generation folder
source('R/target_generation_ThermistorTemp_C_hourly.R')

# FCR
fcr_latest <- "https://raw.githubusercontent.com/FLARE-forecast/FCRE-data/fcre-catwalk-data-qaqc/fcre-waterquality_L1.csv"
fcr_edi <- "https://pasta.lternet.edu/package/data/eml/edi/271/8/fbb8c7a0230f4587f1c6e11417fe9dce"

fcr_thermistor_temp_hourly <- target_generation_ThermistorTemp_C_hourly(current_file = fcr_latest,
                                                                        historic_file = fcr_edi)


cuts <- tibble::tibble(cuts = as.integer(factor(config$model_settings$modeled_depths)),
                       depth = config$model_settings$modeled_depths)

cleaned_insitu_file <- file.path(lake_directory, "targets", config$location$site_id, config$da_setup$obs_filename)

fcr_thermistor_temp_hourly |>
  filter(site_id == forecast_site) |>
  dplyr::mutate(cuts = cut(depth,
                           breaks = config$model_settings$modeled_depths,
                           include.lowest = TRUE,
                           right = FALSE,
                           labels = FALSE),
                variable = ifelse(variable == 'ThermistorTemp_C', 'temperature', variable)) |>
  dplyr::filter(lubridate::hour(datetime) == 0) |>

  dplyr::group_by(cuts, variable, datetime, site_id) |>
  dplyr::summarize(observation = mean(observation, na.rm = TRUE), .groups = "drop") |>
  dplyr::left_join(cuts, by = "cuts") |>
  dplyr::select(site_id, datetime, variable, depth, observation) |>
  write_csv(cleaned_insitu_file)
