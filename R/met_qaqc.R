met_qaqc <- function(realtime_file,
                     qaqc_file,
                     cleaned_met_file,
                     input_file_tz,
                     nldas = NULL){

  if(!is.na(qaqc_file)){
    d1 <- readr::read_csv(realtime_file,
                          col_names = c("TIMESTAMP","RECORD","BattV","PTemp_C","PAR_Den_Avg","PAR_Tot_Tot","BP_kPa_Avg","AirTC_Avg","RH","Rain_mm_Tot","WS_ms_Avg","WindDir","SR01Up_Avg","SR01Dn_Avg","IR01UpCo_Avg","IR01DnCo_Avg","NR01TK_Avg","Albedo_Avg"),
                          col_types = list(
                            TIMESTAMP = readr::col_datetime(format = ""),
                            RECORD = readr::col_integer(),
                            BattV = readr::col_double(),
                            PTemp_C = readr::col_double(),
                            PAR_Den_Avg = readr::col_double(),
                            PAR_Tot_Tot = readr::col_double(),
                            BP_kPa_Avg = readr::col_double(),
                            AirTC_Avg = readr::col_double(),
                            RH = readr::col_double(),
                            Rain_mm_Tot = readr::col_double(),
                            WS_ms_Avg = readr::col_double(),
                            WindDir = readr::col_double(),
                            SR01Up_Avg = readr::col_double(),
                            SR01Dn_Avg = readr::col_double(),
                            IR01UpCo_Avg = readr::col_double(),
                            IR01DnCo_Avg = readr::col_double(),
                            NR01TK_Avg = readr::col_double(),
                            Albedo_Avg = readr::col_double()), skip = 4)
    #d1 <- d1[-85572, ]

    TIMESTAMP_in <- lubridate::force_tz(d1$TIMESTAMP, tzone = input_file_tz)

    d1$TIMESTAMP <- lubridate::with_tz(TIMESTAMP_in,tz = "UTC")

    d2 <- readr::read_csv(qaqc_file,
                          show_col_types = FALSE, progress = FALSE)


    TIMESTAMP_in <- lubridate::force_tz(d2$DateTime, tzone = input_file_tz)


    d2$TIMESTAMP <- lubridate::with_tz(TIMESTAMP_in,tz = "UTC")

    #d3 <- read.csv( fname[3])
    #TIMESTAMP_in <- as.POSIXct(d3$time,
    #                           format= "%Y-%m-%d %H:%M",
    #                           tz = input_file_tz)


    #d3$TIMESTAMP <- with_tz(TIMESTAMP_in,tz = local_tzone)

    d1 <- data.frame(time = d1$TIMESTAMP, ShortWave = d1$SR01Up_Avg, LongWave = d1$IR01UpCo_Avg, AirTemp = d1$AirTC_Avg, RelHum = d1$RH, WindSpeed = d1$WS_ms_Avg, Rain = d1$Rain_mm_Tot, pressure = d1$BP_kPa_Avg)
    d2 <- data.frame(time = d2$TIMESTAMP, ShortWave = d2$ShortwaveRadiationUp_Average_W_m2, LongWave = d2$InfraredRadiationUp_Average_W_m2, AirTemp = d2$AirTemp_Average_C, RelHum = d2$RH_percent, WindSpeed = d2$WindSpeed_Average_m_s, Rain = d2$Rain_Total_mm, pressure = d2$BP_Average_kPa)

    d1 <- d1[which(d1$time > d2$time[nrow(d2)] | d1$time < d2$time[1]), ]

    #d3 <- d3[which(d3$TIMESTAMP < d2$TIMESTAMP[1])]

    d <- rbind(d2, d1)

  }else{

    d1 <- readr::read_csv(realtime_file,
                          col_names = c("TIMESTAMP","RECORD","BattV","PTemp_C","PAR_Den_Avg","PAR_Tot_Tot","BP_kPa_Avg","AirTC_Avg","RH","Rain_mm_Tot","WS_ms_Avg","WindDir","SR01Up_Avg","SR01Dn_Avg","IR01UpCo_Avg","IR01DnCo_Avg","NR01TK_Avg","Albedo_Avg"),
                          col_types = list(
                            TIMESTAMP = readr::col_datetime(format = ""),
                            RECORD = readr::col_integer(),
                            BattV = readr::col_double(),
                            PTemp_C = readr::col_double(),
                            PAR_Den_Avg = readr::col_double(),
                            PAR_Tot_Tot = readr::col_double(),
                            BP_kPa_Avg = readr::col_double(),
                            AirTC_Avg = readr::col_double(),
                            RH = readr::col_double(),
                            Rain_mm_Tot = readr::col_double(),
                            WS_ms_Avg = readr::col_double(),
                            WindDir = readr::col_double(),
                            SR01Up_Avg = readr::col_double(),
                            SR01Dn_Avg = readr::col_double(),
                            IR01UpCo_Avg = readr::col_double(),
                            IR01DnCo_Avg = readr::col_double(),
                            NR01TK_Avg = readr::col_double(),
                            Albedo_Avg = readr::col_double())) %>%
      dplyr::slice(-c(1,2,3,4))

    #d1 <- d1[-85572, ]

    TIMESTAMP_in <- as.POSIXct(d1$TIMESTAMP,
                               format= "%Y-%m-%d %H:%M",
                               tz = input_file_tz)

    d1$TIMESTAMP <- with_tz(TIMESTAMP_in,tz = "UTC")

    d <- data.frame(time = d1$TIMESTAMP, ShortWave = d1$SR01Up_Avg, LongWave = d1$IR01UpCo_Avg, AirTemp = d1$AirTC_Avg, RelHum = d1$RH, WindSpeed = d1$WS_ms_Avg, Rain = d1$Rain_mm_Tot, pressure = d1$BP_kPa_Avg)
  }


  wshgt <- 3
  roughlength <- 0.000114
  d$WindSpeed <- d$WindSpeed * log(10.00 / 0.000114) / log(wshgt / 0.000114)

  maxTempC = 41 # an upper bound of realistic temperature for the study site in deg C
  minTempC = -24 # an lower bound of realistic temperature for the study site in deg C

  equation_of_time <- function(doy) {
    stopifnot(doy <= 367)
    f      <- pi / 180 * (279.5 + 0.9856 * doy)
    et     <- (-104.7 * sin(f) + 596.2 * sin(2 * f) + 4.3 *
                 sin(4 * f) - 429.3 * cos(f) - 2 *
                 cos(2 * f) + 19.3 * cos(3 * f)) / 3600  # equation of time -> eccentricity and obliquity
    return(et)
  }


  cos_solar_zenith_angle <- function(doy, lat, lon, dt, hr) {
    et <- equation_of_time(doy)
    merid  <- floor(lon / 15) * 15
    merid[merid < 0] <- merid[merid < 0] + 15
    lc     <- (lon - merid) * -4/60  ## longitude correction
    tz     <- merid / 360 * 24  ## time zone
    midbin <- 0.5 * dt / 86400 * 24  ## shift calc to middle of bin
    t0   <- 12 + lc - et - tz - midbin  ## solar time
    h    <- pi/12 * (hr - t0)  ## solar hour
    dec  <- -23.45 * pi / 180 * cos(2 * pi * (doy + 10) / 365)  ## declination
    cosz <- sin(lat * pi / 180) * sin(dec) + cos(lat * pi / 180) * cos(dec) * cos(h)
    cosz[cosz < 0] <- 0
    return(cosz)
  }

  d <- d %>%
    dplyr::mutate(ShortWave = ifelse(ShortWave < 0, 0, ShortWave),
                  RelHum = ifelse(RelHum < 0, 0, RelHum),
                  RelHum = ifelse(RelHum > 100, 100, RelHum),
                  AirTemp = ifelse(AirTemp> maxTempC, NA, AirTemp),
                  AirTemp = ifelse(AirTemp < minTempC, NA, AirTemp),
                  LongWave = ifelse(LongWave < 0, NA, LongWave),
                  WindSpeed = ifelse(WindSpeed < 0, 0, WindSpeed)) %>%
    filter(is.na(time) == FALSE)

  d <- d %>%
    mutate(day = lubridate::day(time),
           year = lubridate::year(time),
           hour = lubridate::hour(time),
           month = lubridate::month(time),
           minute = lubridate::minute(time)) %>%
           #AirTemp = ifelse(minute <= 10, AirTemp, NA),
           #RelHum = ifelse(minute <= 10, RelHum, NA),
           #WindSpeed = ifelse(minute <= 10, WindSpeed, NA),
           #pressure = ifelse(minute <= 10, pressure, NA)) %>%
    group_by(day, year, hour, month) %>%
    summarize(ShortWave = mean(ShortWave, na.rm = TRUE),
              LongWave = mean(LongWave, na.rm = TRUE),
              AirTemp = mean(AirTemp, na.rm = TRUE),
              RelHum = mean(RelHum, na.rm = TRUE),
              WindSpeed = mean(WindSpeed, na.rm = TRUE),
              pressure = mean(pressure, na.rm = TRUE),
              Rain = sum(Rain), .groups = "drop") %>%
    mutate(day = as.numeric(day),
           hour = as.numeric(hour)) %>%
    mutate(day = ifelse(as.numeric(day) < 10, paste0("0",day),day),
           hour = ifelse(as.numeric(hour) < 10, paste0("0",hour),hour)) %>%
    mutate(time = lubridate::as_datetime(paste0(year,"-",month,"-",day," ",hour,":00:00"),tz = "UTC")) %>%
    dplyr::select(time,ShortWave,LongWave,AirTemp,RelHum,WindSpeed,Rain,pressure) %>%
    arrange(time)

  d <- d  %>%
    mutate(hr = lubridate::hour(time),
           doy = lubridate::yday(time) + hr/24.)

  dt <- median(diff(d$doy)) * 86400 # average number of seconds in time interval
  d$hr <- (d$doy - floor(d$doy)) * 24 # hour of day for each element of doy
  lat <- 37.27
  lon <- 360-79.9
  ## calculate potential radiation
  d$cosz <- cos_solar_zenith_angle(doy = d$doy, lat, lon, dt, hr = d$hr)
  d$rpot <- 1366 * d$cosz

  d |>
    select(time, ShortWave, rpot) |>
    pivot_longer(cols = -time, names_to = "variable", values_to = "value") |>
    mutate(hour = hour(time)) |>
    ggplot(aes(x = time, y = value, color = variable)) +
    geom_point() +
    facet_grid(~hour)

  d <- d %>%
    mutate(ShortWave = ifelse(ShortWave > rpot, NA, ShortWave))

  d <- d %>%
    rename(surface_downwelling_shortwave_flux_in_air = ShortWave,
           surface_downwelling_longwave_flux_in_air = LongWave,
           air_temperature = AirTemp,
           relative_humidity = RelHum,
           wind_speed = WindSpeed,
           precipitation_flux = Rain,
           air_pressure = pressure)

  d <- d %>%
    mutate(air_temperature = air_temperature + 273.15,
           relative_humidity = relative_humidity / 100)

  #Note that mm hr-1 is the same as kg m2 hr-1. Converting to kg m2 s-1
  d$precipitation_flux <- d$precipitation_flux / (60 * 60)

  d$air_pressure <- d$air_pressure * 1000

  d$specific_humidity <-  rh2qair(rh = d$relative_humidity,
                                                  T = d$air_temperature,
                                                  press = d$air_pressure)

  d <- d %>%
    select(time, air_temperature, air_pressure, relative_humidity, surface_downwelling_longwave_flux_in_air, surface_downwelling_shortwave_flux_in_air, precipitation_flux, specific_humidity, wind_speed)

  cf_var_names1 <- c("air_temperature", "air_pressure", "relative_humidity", "surface_downwelling_longwave_flux_in_air",
                     "surface_downwelling_shortwave_flux_in_air", "precipitation_flux","specific_humidity","wind_speed")

  cf_var_units1 <- c("K", "Pa", "1", "Wm-2", "Wm-2", "kgm-2s-1", "1", "ms-1")  #Negative numbers indicate negative exponents

  d$time <- lubridate::with_tz(d$time, tzone = "UTC")

  #d <- d %>%
  #  tidyr::drop_na()

  if(!is.null(nldas)){

    print("Gap filling with NLDAS")
    d_nldas <- readr::read_csv(nldas, col_type = readr::cols())

    d_nldas <- d_nldas %>%
      rename(air_temperature = AirTemp,
             surface_downwelling_shortwave_flux_in_air = ShortWave,
             surface_downwelling_longwave_flux_in_air = LongWave,
             relative_humidity = RelHum,
             wind_speed = WindSpeed,
             precipitation_flux = Rain) %>%
      mutate(air_pressure = NA,
             specific_humidity = NA,
             air_temperature = air_temperature + 273.15,
             relative_humidity = relative_humidity/ 100,
             precipitation_flux = precipitation_flux * 1000 / (60 * 60 * 24)) %>%
      select(all_of(names(d)))

    d_nldas$time <- lubridate::with_tz(d_nldas$time, tzone = "UTC")

    d_nldas <- d_nldas %>%
      filter(time >= min(d$time), time <= max(d$time))

    d_nldas_gaps <- d_nldas %>%
      filter(!(time %in% d$time))

    d_full <- rbind(d, d_nldas_gaps) %>%
      arrange(time)

  }else{

    doy_mean <- d %>%
      mutate(doy = yday(time),
             hour = hour(time)) %>%
      group_by(doy, hour) %>%
      summarize(filled_shortwave = mean(surface_downwelling_shortwave_flux_in_air, na.rm = TRUE),
                filled_longwave = mean(surface_downwelling_longwave_flux_in_air, na.rm = TRUE),
                filled_temp = mean(air_temperature, na.rm = TRUE),
                filled_relative_humidity = mean(relative_humidity, na.rm = TRUE),
                filled_precip = mean(precipitation_flux, na.rm = TRUE),
                filled_specific_humidity = mean(specific_humidity, na.rm = TRUE),
                filled_wind_speed = mean(wind_speed, na.rm = TRUE),
                filled_air_press = mean(air_pressure, na.rm = TRUE),
                .groups = "drop")


    t <- seq(min(d$time), max(d$time), by = "1 hour")
    cont_time <- tibble(time = t)

    d_full <- dplyr::left_join(cont_time, d, by = "time") %>%
      mutate(doy = yday(time),
             hour = hour(time)) %>%
      left_join(doy_mean, by = c("doy","hour"))

    d_full <- d_full %>%
      mutate(surface_downwelling_shortwave_flux_in_air = ifelse(is.na(surface_downwelling_shortwave_flux_in_air), filled_shortwave, surface_downwelling_shortwave_flux_in_air),
             surface_downwelling_longwave_flux_in_air= ifelse(is.na(surface_downwelling_longwave_flux_in_air), filled_longwave, surface_downwelling_longwave_flux_in_air),
             air_temperature = ifelse(is.na(air_temperature), filled_temp, air_temperature),
             relative_humidity = ifelse(is.na(relative_humidity), filled_relative_humidity, relative_humidity),
             precipitation_flux = ifelse(is.na(precipitation_flux), filled_precip, precipitation_flux),
             specific_humidity = ifelse(is.na(specific_humidity), filled_specific_humidity, specific_humidity),
             wind_speed = ifelse(is.na(wind_speed), filled_wind_speed, wind_speed),
             air_pressure = ifelse(is.na(air_pressure), filled_air_press, air_pressure))
  }

  lat <- 37.27
  lon <- 360-79.9
  start_time <- dplyr::first((d$time))
  end_time <- dplyr::last((d$time))
  cf_units <- cf_var_units1

  output_file <- cleaned_met_file

  if(!dir.exists(dirname(cleaned_met_file))){
    dir.create(dirname(cleaned_met_file), recursive = TRUE)
  }

  start_time <- min(d$time)
  end_time <- max(d$time)

  data <- d_full %>%
    dplyr::select(-time)

  diff_time <- as.numeric(difftime(d_full$time, d_full$time[1], units = "hours"))

  cf_var_names <- names(data)

  time_dim <- ncdf4::ncdim_def(name="time",
                               units = paste("hours since", format(start_time, "%Y-%m-%d %H:%M")),
                               diff_time, #GEFS forecast starts 6 hours from start time
                               create_dimvar = TRUE)
  lat_dim <- ncdf4::ncdim_def("latitude", "degree_north", lat, create_dimvar = TRUE)
  lon_dim <- ncdf4::ncdim_def("longitude", "degree_east", lon, create_dimvar = TRUE)

  dimensions_list <- list(time_dim, lat_dim, lon_dim)

  nc_var_list <- list()
  for (i in 1:length(cf_var_names)) { #Each ensemble member will have data on each variable stored in their respective file.
    nc_var_list[[i]] <- ncdf4::ncvar_def(cf_var_names[i], cf_units[i], dimensions_list, missval=NaN)
  }

  nc_flptr <- ncdf4::nc_create(output_file, nc_var_list, verbose = FALSE, )

  #For each variable associated with that ensemble
  for (j in 1:ncol(data)) {
    # "j" is the variable number.  "i" is the ensemble number. Remember that each row represents an ensemble
    ncdf4::ncvar_put(nc_flptr, nc_var_list[[j]], unlist(data[,j]))
  }

  ncdf4::nc_close(nc_flptr)  #Write to the disk/storage

  return(cleaned_met_file)

}


##' converts relative humidity to specific humidity
##' @title RH to SH
##' @param rh relative humidity (proportion, not percent)
##' @param T absolute temperature (Kelvin)
##' @param press air pressure (Pascals)
##' @noRd
##' @author Mike Dietze, Ankur Desai
##' @aliases rh2rv
rh2qair <- function(rh, T, press = 101325) {
  stopifnot(T[!is.na(T)] >= 0)
  Tc <- T - 273.15
  es <- 6.112 * exp((17.67 * Tc) / (Tc + 243.5))
  e <- rh * es
  p_mb <- press / 100
  qair <- (0.622 * e) / (p_mb - (0.378 * e))
  ## qair <- rh * 2.541e6 * exp(-5415.0 / T) * 18/29
  return(qair)
} # rh2qair

