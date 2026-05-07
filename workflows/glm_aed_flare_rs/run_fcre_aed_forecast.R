run_fcre_aed_forecast <- function(config_set_name    = "glm_aed_flare_rs",
                                  configure_run_file = "configure_run.yml") {

  library(tidyverse)
  library(lubridate)
  set.seed(100)

  Sys.setenv('GLM_PATH' = 'GLM3r')
  options(future.globals.maxSize = 891289600)

  Sys.setenv("AWS_DEFAULT_REGION" = "amnh1",
             "AWS_S3_ENDPOINT"    = "osn.mghpcc.org",
             "USE_HTTPS"          = TRUE)

  marker <- file.path("configuration", config_set_name, configure_run_file)
  candidates <- list.files("/tmp/functions",
                           pattern = "configure_run\\.yml$",
                           recursive = TRUE, full.names = TRUE)
  candidates <- candidates[grepl(marker, candidates, fixed = TRUE)]
  lake_directory <- if (length(candidates) > 0) {
    sub(paste0("/", marker, "$"), "", candidates[1])
  } else {
    here::here()
  }
  setwd(lake_directory)

  source(file.path(lake_directory, "R/convert_vera4cast_inflow.R"), local = TRUE)
  source(file.path(lake_directory, "R/generate_forecast_score_arrow.R"), local = TRUE)

  config <- FLAREr::set_up_simulation(configure_run_file = configure_run_file,
                                      lake_directory     = lake_directory,
                                      config_set_name    = config_set_name)

  # WORKAROUND: FLAREr::check_noaa_present catches arrow IOError when the
  # NOAA partition is missing, but its handler does `message(error_message)`
  # with the condition object, which re-signals the error and aborts the
  # action instead of returning FALSE. Drop this wrapper once FLAREr's
  # check_noaa_present uses message(conditionMessage(e)).
  noaa_ready <- tryCatch(
    FLAREr::check_noaa_present(lake_directory,
                               configure_run_file,
                               config_set_name = config_set_name),
    error = function(e) {
      message("check_noaa_present errored (treating as not-ready): ",
              conditionMessage(e))
      FALSE
    }
  )

  reference_date <- lubridate::as_date(config$run_config$forecast_start_datetime)
  inflow_prefix <- "vera4cast/forecasts/archive-parquet/project_id=vera4cast/duration=P1D/variable=Temp_C_mean/model_id=inflow_gefsClimAED"
  s3 <- FLAREr::flare_arrow_s3_bucket(server_name  = "vera4cast_forecasts",
                                      faasr_prefix = inflow_prefix,
                                      config       = config)
  avail_dates <- gsub("reference_date=", "", s3$ls())
  inflow_ready <- reference_date %in% lubridate::as_date(avail_dates)

  message(paste0("noaa ready: ", noaa_ready))
  message(paste0("inflow ready: ", inflow_ready))

  loop_start <- Sys.time()
  loop_budget_seconds <- 5.5 * 60 * 60

  while (noaa_ready & inflow_ready) {

    if (as.numeric(Sys.time() - loop_start, units = "secs") > loop_budget_seconds) {
      message("Approaching action time budget; exiting loop cleanly with restart written.")
      break
    }

    source(file.path(lake_directory, "workflows", config_set_name, "generate_inflow_forecast.R"), local = TRUE)

    insitu_local <- file.path(tempdir(), "daily-insitu-targets.csv.gz")
    FLAREr::flare_get_file(local_file    = "daily-insitu-targets.csv.gz",
                           remote_file   = "daily-insitu-targets.csv.gz",
                           server_name   = "vera4cast_targets",
                           local_folder  = tempdir(),
                           remote_folder = "vera4cast/targets/project_id=vera4cast/duration=P1D",
                           config        = config)

    readr::read_csv(insitu_local, show_col_types = FALSE) |>
      dplyr::mutate(observation = ifelse(variable == "DO_mgL_mean", observation*1000*(1/32), observation),
                    observation = ifelse(variable == "fDOM_QSU_mean", -151.3407 + observation*29.62654, observation),
                    depth_m = ifelse(depth_m == 0.1, 0.0, depth_m)) |>
      dplyr::rename(depth = depth_m) |>
      dplyr::filter(site_id == "fcre",
                    datetime >= as_datetime(config$run_config$start_datetime)) |>
      dplyr::mutate(datetime = lubridate::as_datetime(datetime)) |>
      readr::write_csv(file.path(config$file_path$qaqc_data_directory,
                                 paste0(config$location$site_id, "-targets-insitu.csv")))

    source(file.path(lake_directory, "workflows", config_set_name, "getLST.R"), local = TRUE)
    start_iso <- format(lubridate::as_datetime(config$run_config$start_datetime),         "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
    end_iso   <- format(lubridate::as_datetime(config$run_config$forecast_start_datetime), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
    data <- tryCatch(
      get_lst(bbox, start_iso, end_iso),
      error = function(e) { message("get_lst failed: ", conditionMessage(e)); NULL }
    )
    vals <- if (!is.null(data)) {
      tryCatch(get_vals(points, data), error = function(e) NULL)
    } else NULL
    if (!is.null(vals)) {
      clean_data(vals) |>
        readr::write_csv(file.path(config$file_path$qaqc_data_directory,
                                   paste0(config$location$site_id, "-targets-rs.csv")))
    } else {
      message("No new RS data")
      data <- as.data.frame(cbind("datetime"    = paste0(substr(as.character(config$run_config$start_datetime), 0, 10), "T00:00:00Z"),
                                  "observation" = NA,
                                  "site_id"     = "fcre",
                                  "depth"       = 0,
                                  "variable"    = "temperature"))
      data |>
        readr::write_csv(file.path(config$file_path$qaqc_data_directory,
                                   paste0(config$location$site_id, "-targets-rs.csv")))
    }

    FLAREr::run_flare(lake_directory     = lake_directory,
                      configure_run_file = configure_run_file,
                      config_set_name    = config_set_name)

    ref_date <- as.character(lubridate::as_date(config$run_config$forecast_start_datetime))
    forecasts_s3 <- FLAREr::flare_arrow_s3_bucket(
      server_name  = "forecasts_parquet",
      faasr_prefix = paste0("flare/forecasts/parquet",
                            "/site_id=", config$location$site_id,
                            "/model_id=", config$run_config$sim_name,
                            "/reference_date=", ref_date),
      config       = config
    )
    forecast_df <- arrow::open_dataset(forecasts_s3) |>
      collect() |>
      mutate(datetime       = lubridate::as_datetime(datetime),
             site_id        = config$location$site_id,
             model_id       = config$run_config$sim_name,
             reference_date = ref_date)

    vera_variables <- c("Temp_C_mean","Chla_ugL_mean", "DO_mgL_mean", "fDOM_QSU_mean", "NH4_ugL_sample",
                        "NO3NO2_ugL_sample", "SRP_ugL_sample", "DIC_mgL_sample","Secchi_m_sample",
                        "Bloom_binary_mean","CH4_umolL_sample","IceCover_binary_max", "CO2flux_umolm2s_mean", "CH4flux_umolm2s_mean",
                        "Mixed_binary_mean")

    bloom_binary <- forecast_df |>
      dplyr::filter(depth == 1.6 & variable == "Chla_ugL_mean") |>
      dplyr::mutate(over = ifelse(prediction > 20, 1, 0)) |>
      dplyr::summarize(prediction = sum(over) / n(), .by = c(datetime, reference_datetime, model_id, site_id, depth, variable)) |>
      dplyr::mutate(family    = "bernoulli",
                    parameter = "prob",
                    variable  = "Bloom_binary_mean",
                    datetime  = lubridate::as_datetime(datetime)) |>
      dplyr::rename(depth_m = depth) |>
      dplyr::select(reference_datetime, datetime, model_id, site_id, depth_m, family, parameter, variable, prediction)

    ice_binary <- forecast_df |>
      dplyr::filter(variable == "ice_thickness") |>
      dplyr::mutate(over = ifelse(prediction > 0, 1, 0)) |>
      dplyr::summarize(prediction = sum(over) / n(), .by = c(datetime, reference_datetime, model_id, site_id, depth, variable)) |>
      dplyr::mutate(family    = "bernoulli",
                    parameter = "prob",
                    variable  = "IceCover_binary_max",
                    depth     = NA,
                    datetime  = lubridate::as_datetime(datetime)) |>
      dplyr::rename(depth_m = depth) |>
      dplyr::select(reference_datetime, datetime, model_id, site_id, depth_m, family, parameter, variable, prediction)

    min_depth <- 1
    max_depth <- 8
    threshold <- 0.1

    temp_forecast <- forecast_df |>
      filter(variable %in% c("temp_1.0m_mean","temp_8.0m_mean")) |>
      mutate(depth    = ifelse(variable == "temp_1.0m_mean", 1.0, 8.0),
             variable = "Temp_C_mean",
             datetime = lubridate::as_datetime(datetime - lubridate::days(1))) |>
      pivot_wider(names_from = depth, names_prefix = "wtr_", values_from = prediction)

    colnames(temp_forecast)[which(colnames(temp_forecast) == paste0("wtr_", min_depth))] <- "min_depth"
    colnames(temp_forecast)[which(colnames(temp_forecast) == paste0("wtr_", max_depth))] <- "max_depth"

    mix_binary <- temp_forecast |>
      mutate(min_depth = rLakeAnalyzer::water.density(min_depth),
             max_depth = rLakeAnalyzer::water.density(max_depth),
             mixed     = ifelse((max_depth - min_depth) < threshold, 1, 0)) |>
      summarise(prediction = (sum(mixed)/n()), .by = c(datetime, reference_datetime, model_id, site_id, variable)) |>
      dplyr::mutate(family    = "bernoulli",
                    parameter = "prob",
                    variable  = "Mixed_binary_mean",
                    depth     = NA,
                    datetime  = lubridate::as_datetime(datetime)) |>
      dplyr::rename(depth_m = depth) |>
      dplyr::select(reference_datetime, datetime, model_id, site_id, depth_m, family, parameter, variable, prediction)

    vera4cast_df <- forecast_df |>
      dplyr::rename(depth_m = depth) |>
      dplyr::mutate(variable   = ifelse(variable == "oxy_mean", "DO_mgL_mean", variable),
                    depth_m    = ifelse(variable == "DO_mgL_mean", 1.6, depth_m),
                    datetime   = ifelse(variable == "DO_mgL_mean", datetime - lubridate::days(1), datetime),
                    prediction = ifelse(variable == "DO_mgL_mean", prediction/1000*(32), prediction),
                    variable   = ifelse(variable == "Temp_C_mean", "Temp_C_mean_all_depth", variable),
                    variable   = ifelse(variable == "temp_1.6m_mean", "Temp_C_mean", variable),
                    depth_m    = ifelse(variable == "Temp_C_mean", 1.6, depth_m),
                    datetime   = ifelse(variable == "Temp_C_mean", datetime - lubridate::days(1), datetime),
                    prediction = ifelse(variable == "fDOM_QSU_mean", (151.3407 + prediction)/29.62654, prediction),
                    prediction = ifelse(variable == "NIT_amm", prediction/1000/0.001/(1/18.04), prediction),
                    variable   = ifelse(variable == "NIT_amm", "NH4_ugL_sample", variable),
                    prediction = ifelse(variable == "NIT_nit", prediction/1000/0.001/(1/62.00), prediction),
                    variable   = ifelse(variable == "NIT_amm", "NO3NO2_ugL_sample", variable),
                    prediction = ifelse(variable == "PHS_frp", prediction/1000/0.001/(1/94.9714), prediction),
                    variable   = ifelse(variable == "PHS_frp", "SRP_ugL_sample", variable),
                    prediction = ifelse(variable == "CAR_dic", prediction/1000/(1/52.515), prediction),
                    variable   = ifelse(variable == "CAR_dic", "DIC_mgL_sample", variable),
                    variable   = ifelse(variable == "CAR_ch4", "CH4_umolL_sample", variable),
                    variable   = ifelse(variable == "secchi", "Secchi_m_sample", variable),
                    prediction = ifelse(variable == "co2_flux_mean", prediction/0.001/86400, prediction),
                    variable   = ifelse(variable == "co2_flux_mean", "CO2flux_umolm2s_mean", variable),
                    prediction = ifelse(variable == "ch4_flux_mean", prediction/0.001/86400, prediction),
                    variable   = ifelse(variable == "ch4_flux_mean", "CH4flux_umolm2s_mean", variable),
                    depth_m    = ifelse(depth_m == 0.0, 0.1, depth_m),
                    datetime   = lubridate::as_datetime(datetime)) |>
      dplyr::select(-forecast, -variable_type) |>
      dplyr::mutate(parameter = as.character(parameter)) |>
      dplyr::bind_rows(bloom_binary) |>
      dplyr::bind_rows(ice_binary) |>
      dplyr::bind_rows(mix_binary) |>
      dplyr::filter(variable %in% vera_variables) |>
      mutate(project_id         = "vera4cast",
             model_id           = config$run_config$sim_name,
             family             = "ensemble",
             site_id            = "fcre",
             duration           = "P1D",
             datetime           = lubridate::as_datetime(datetime),
             reference_datetime = lubridate::as_datetime(reference_datetime)) |>
      filter(datetime >= reference_datetime) |>
      distinct(reference_datetime, datetime, variable, depth_m, parameter, model_id, .keep_all = TRUE)

    file_name <- paste0(config$run_config$sim_name, "-",
                        lubridate::as_date(vera4cast_df$reference_datetime[1]), ".csv.gz")
    readr::write_csv(vera4cast_df, file = file_name)

    message("Scoring forecasts")

    forecasts_s3 <- FLAREr::flare_arrow_s3_bucket(
      server_name  = "forecasts_parquet",
      faasr_prefix = paste0("flare/forecasts/parquet",
                            "/site_id=", config$location$site_id,
                            "/model_id=", config$run_config$sim_name,
                            "/reference_date=", as.character(lubridate::as_date(config$run_config$forecast_start_datetime))),
      config       = config
    )
    forecast_df <- arrow::open_dataset(forecasts_s3) |>
      dplyr::collect() |>
      dplyr::mutate(site_id        = config$location$site_id,
                    model_id       = config$run_config$sim_name,
                    reference_date = lubridate::as_date(config$run_config$forecast_start_datetime))

    if (config$output_settings$evaluate_past & config$run_config$use_s3) {
      past_days <- lubridate::as_date(lubridate::as_date(config$run_config$forecast_start_datetime) -
                                        lubridate::days(config$run_config$forecast_horizon))
      past_s3 <- FLAREr::flare_arrow_s3_bucket(
        server_name  = "forecasts_parquet",
        faasr_prefix = paste0("flare/forecasts/parquet",
                              "/site_id=", config$location$site_id,
                              "/model_id=", config$run_config$sim_name,
                              "/reference_date=", as.character(past_days)),
        config       = config
      )
      past_forecasts <- arrow::open_dataset(past_s3) |>
        dplyr::collect() |>
        dplyr::mutate(site_id        = config$location$site_id,
                      model_id       = config$run_config$sim_name,
                      reference_date = past_days)
    } else {
      past_forecasts <- NULL
    }

    combined_forecasts <- dplyr::bind_rows(forecast_df, past_forecasts)

    targets_df <- read_csv(file.path(config$file_path$qaqc_data_directory,
                                     paste0(config$location$site_id, "-targets-insitu.csv")),
                           show_col_types = FALSE)

    scoring <- generate_forecast_score_arrow(targets_df       = targets_df,
                                             forecast_df      = combined_forecasts,
                                             use_s3           = config$run_config$use_s3,
                                             bucket           = config$s3$scores$bucket,
                                             endpoint         = config$s3$scores$endpoint,
                                             local_directory  = "./scores/fcre",
                                             variable_types   = c("state", "parameter"),
                                             config           = config)

    forecast_start_datetime <- lubridate::as_datetime(config$run_config$forecast_start_datetime) + lubridate::days(1)
    start_datetime          <- lubridate::as_datetime(config$run_config$forecast_start_datetime) - lubridate::days(4)
    restart_file <- paste0(config$location$site_id, "-",
                           (lubridate::as_date(forecast_start_datetime) - days(1)),
                           "-", config$run_config$sim_name, ".nc")

    FLAREr::update_run_config(lake_directory          = lake_directory,
                              configure_run_file      = configure_run_file,
                              restart_file            = restart_file,
                              start_datetime          = start_datetime,
                              end_datetime            = NA,
                              forecast_start_datetime = forecast_start_datetime,
                              forecast_horizon        = config$run_config$forecast_horizon,
                              sim_name                = config$run_config$sim_name,
                              site_id                 = config$location$site_id,
                              configure_flare         = config$run_config$configure_flare,
                              configure_obs           = config$run_config$configure_obs,
                              use_s3                  = config$run_config$use_s3,
                              bucket                  = config$s3$restart$bucket,
                              endpoint                = config$s3$restart$endpoint,
                              config                  = config,
                              use_https               = TRUE)

    # POSSIBLE BUG FIX -- safe to remove if upstream behavior changes.
    # update_run_config just wrote a new forecast_start_datetime to
    # disk + S3, but the in-memory `config` still holds the previous
    # values. Sourced helpers on the next iteration
    # (generate_inflow_forecast.R -> convert_vera4cast_inflow) read
    # from this outer `config` and would write inflow partitions under
    # the previous reference_date, while run_flare re-reads config
    # from disk and looks under the new one -- producing an
    # arrow::open_dataset IOError. Refreshing here keeps the in-memory
    # view aligned with on-disk state. combined_run_aed.R doesn't need
    # this because each cron firing runs the loop exactly once.
    config <- FLAREr::set_up_simulation(configure_run_file = configure_run_file,
                                        lake_directory     = lake_directory,
                                        config_set_name    = config_set_name)

    var1 <- Sys.getenv("AWS_ACCESS_KEY_ID")
    var2 <- Sys.getenv("AWS_SECRET_ACCESS_KEY")
    Sys.unsetenv("AWS_ACCESS_KEY_ID")
    Sys.unsetenv("AWS_SECRET_ACCESS_KEY")
    tryCatch(
      vera4castHelpers::submit(file_name, first_submission = FALSE),
      error = function(e) message("vera4cast submit failed (continuing): ", conditionMessage(e))
    )
    Sys.setenv("AWS_ACCESS_KEY_ID"     = var1,
               "AWS_SECRET_ACCESS_KEY" = var2)

    # See note above on the initial call site.
    noaa_ready <- tryCatch(
      FLAREr::check_noaa_present(lake_directory,
                                 configure_run_file,
                                 config_set_name = config_set_name),
      error = function(e) {
        message("check_noaa_present errored (treating as not-ready): ",
                conditionMessage(e))
        FALSE
      }
    )

    reference_date <- lubridate::as_date(forecast_start_datetime)
    s3 <- FLAREr::flare_arrow_s3_bucket(server_name  = "vera4cast_forecasts",
                                        faasr_prefix = inflow_prefix,
                                        config       = config)
    avail_dates <- gsub("reference_date=", "", s3$ls())
    inflow_ready <- reference_date %in% avail_dates
  }

  invisible(NULL)
}
