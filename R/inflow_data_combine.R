###Information and extra lines that might be helpful###
#####
### This function aggregates historical and real time inflow data at FCR,
### converts flow to cms. This creates one output file named "inflow_postQAQC.csv"

#Inputs to function:
#  Historical flow data. Here, I need two historical files from github to fill all of the time from 2013 until diana was launched ('hist1' and hist2').
#  Diana data ('FCRweir.csv' from git)
#  Local timezone
#  Timezone of input files (just to be similar to other FLARE fucntions)

#Units:
#  flow = cms
#  wtr_temp = degC
#  Timestep = daily
#Missing data is assigned the prior days flow

##These lines may be useful for getting this into FLARE
#
#     data_location <- paste0("your working directory")
#     diana_location <- paste0(data_location, "/", "diana_data")

#     if(!file.exists(diana_location)){
#       setwd(data_location)
#       system("git clone -b diana-data --single-branch https://github.com/CareyLabVT/SCCData.git diana-data/")
#     }
#
#
#     setwd(diana_location)
#     system(paste0("git pull"))
#####
inflow_data_combine <- function(realtime_file,
                            qaqc_file,
                            nutrients_file,
                            silica_file,
                            co2_ch4,
                            cleaned_inflow_file,
                            input_file_tz,
                            site_id){

  ##Step 2: Read in historical flow data, clean, and aggregate to daily mean##

  flow <- readr::read_csv(qaqc_file, guess_max = 1000000, show_col_types = FALSE) %>%
    dplyr::rename("timestamp" = DateTime) %>%
    dplyr::mutate(timestamp = lubridate::as_datetime(timestamp, tz = input_file_tz),
                  timestamp = lubridate::with_tz(timestamp, tzone = "UTC")) %>%
    dplyr::select(timestamp, WVWA_Flow_cms, WVWA_Temp_C, VT_Flow_cms, VT_Temp_C) %>%
    dplyr::mutate(date = as_date(timestamp)) %>%
    dplyr::group_by(date) %>%
    dplyr:: summarize(WVWA_Flow_cms = mean(WVWA_Flow_cms, na.rm = TRUE),
                      WVWA_Temp_C = mean(WVWA_Temp_C, na.rm = TRUE),
                      VT_Flow_cms = mean(VT_Flow_cms, na.rm = TRUE),
                      VT_Temp_C = mean(VT_Temp_C, na.rm = TRUE), .groups = "drop") %>%
    dplyr::ungroup() %>%
    dplyr::select(date,WVWA_Flow_cms,WVWA_Temp_C,VT_Flow_cms,VT_Temp_C) %>%
    dplyr::mutate(VT_Flow_cms = ifelse(is.nan(VT_Flow_cms), NA, VT_Flow_cms),
                  VT_Temp_C = ifelse(is.nan(VT_Temp_C), NA, VT_Temp_C),
                  WVWA_Flow_cms = ifelse(is.nan(WVWA_Flow_cms), NA, WVWA_Flow_cms),
                  WVWA_Temp_C = ifelse(is.nan(WVWA_Temp_C), NA, WVWA_Temp_C)) %>%
    dplyr::rename("time" = date) %>%
    dplyr::filter(!is.na(time)) %>%
    dplyr::arrange(time)

  inflow_temp_flow <- tibble::tibble(time = seq(first(flow$time), last(flow$time), by = "1 day")) %>%
    dplyr::left_join(flow, by = "time") %>%
    dplyr::mutate(TEMP = ifelse(is.na(VT_Temp_C), WVWA_Temp_C, VT_Temp_C),
                  FLOW = ifelse(time > as_date("2019-06-07"), VT_Flow_cms, WVWA_Flow_cms),
                  SALT = 0) %>%
    dplyr::mutate(TEMP = imputeTS::na_interpolation(TEMP),
                  FLOW = imputeTS::na_interpolation(FLOW)) %>%
    dplyr::select(time, FLOW, TEMP, SALT)

  ##Step 3: Read in diana data, convert flow from PSI to CSM, calculations to
  #account for building new weir in June 2019 (FCR Specific), and
  #aggregate to daily mean.##

  if(!is.na(realtime_file)){
  inflow_realtime <- readr::read_csv(realtime_file, guess_max = 1000000, show_col_types = FALSE) %>%
    dplyr::rename("timestamp" = DateTime) %>%
    dplyr::mutate(timestamp = lubridate::as_datetime(timestamp, tz = input_file_tz),
                  timestamp = lubridate::with_tz(timestamp, tzone = "UTC")) %>%
    dplyr::select(timestamp, WVWA_Flow_cms, WVWA_Temp_C, VT_Flow_cms, VT_Temp_C) %>%
    dplyr::mutate(date = as_date(timestamp)) %>%
    dplyr::group_by(date) %>%
    dplyr:: summarize(WVWA_Flow_cms = mean(WVWA_Flow_cms, na.rm = TRUE),
                      WVWA_Temp_C = mean(WVWA_Temp_C, na.rm = TRUE),
                      VT_Flow_cms = mean(VT_Flow_cms, na.rm = TRUE),
                      VT_Temp_C = mean(VT_Temp_C, na.rm = TRUE), .groups = "drop") %>%
    dplyr::ungroup() %>%
    dplyr::select(date,WVWA_Flow_cms,WVWA_Temp_C,VT_Flow_cms,VT_Temp_C) %>%
    dplyr::mutate(VT_Flow_cms = ifelse(is.nan(VT_Flow_cms), NA, VT_Flow_cms),
                  VT_Temp_C = ifelse(is.nan(VT_Temp_C), NA, VT_Temp_C),
                  WVWA_Flow_cms = ifelse(is.nan(WVWA_Flow_cms), NA, WVWA_Flow_cms),
                  WVWA_Temp_C = ifelse(is.nan(WVWA_Temp_C), NA, WVWA_Temp_C)) %>%
    dplyr::rename("time" = date) %>%
    dplyr::filter(!is.na(time)) %>%
    dplyr::arrange(time) %>%
    dplyr::filter(time > max(inflow_temp_flow$time)) ##REMOVE LATER -- this is a placeholder until we do a PR for main repo

  inflow_realtime_temp_flow <- inflow_realtime %>%
    dplyr::mutate(TEMP = ifelse(is.na(VT_Temp_C), WVWA_Temp_C, VT_Temp_C),
                  FLOW = ifelse(time > as_date("2019-06-07"), VT_Flow_cms, WVWA_Flow_cms),
                  SALT = 0) %>%
    dplyr::mutate(TEMP = imputeTS::na_interpolation(TEMP),
                  FLOW = imputeTS::na_interpolation(FLOW)) %>%
    dplyr::select(time, FLOW, TEMP, SALT)

  inflow_combined <- dplyr::full_join(inflow_temp_flow, inflow_realtime_temp_flow, by = "time") %>%
    dplyr::mutate(FLOW = ifelse(is.na(FLOW.x), FLOW.y, FLOW.x),
                  TEMP = ifelse(is.na(TEMP.x), TEMP.y, TEMP.x),
                  SALT = ifelse(is.na(SALT.x), SALT.y, SALT.x)) %>%
    dplyr::select(time, FLOW, TEMP, SALT)

  # inflow_combined <- tibble::tibble(time = seq(first(flow$time), last(flow$time), by = "1 day")) %>%
  #   dplyr::left_join(flow, by = "time") %>%
  #   dplyr::mutate(TEMP = ifelse(is.na(VT_Temp_C), WVWA_Temp_C, VT_Temp_C),
  #                 FLOW = ifelse(time > as_date("2019-06-07"), VT_Flow_cms, WVWA_Flow_cms),
  #                 SALT = 0) %>%
  #   dplyr::mutate(TEMP = imputeTS::na_interpolation(TEMP),
  #                 FLOW = imputeTS::na_interpolation(FLOW)) %>%
  #   dplyr::select(time, FLOW, TEMP, SALT)
  }else{
    inflow_combined <- inflow_temp_flow |>
    dplyr::select(time, FLOW, TEMP, SALT)
  }


  #### BRING IN THE NUTRIENTS

  if(!is.na(nutrients_file)){

    FCRchem <- readr::read_csv(nutrients_file, show_col_types = FALSE) %>%
      dplyr::select(Reservoir:DIC_mgL) %>%
      dplyr::filter(Reservoir == "FCR") %>%
      dplyr::filter(Site == 100) %>% #inflow site code
      dplyr::rename(time = DateTime) %>%
      dplyr::mutate(time = lubridate::force_tz(time, tzone = input_file_tz),
                    time = lubridate::with_tz(time, tzone = "UTC"),
                    time = lubridate::as_date(time)) %>%
      #mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) %>%
      dplyr::filter(TP_ugL < 100) %>% #remove outliers
      dplyr::select(time:DIC_mgL)

    silica <- readr::read_csv(silica_file, show_col_types = FALSE) %>%
      dplyr::filter(Reservoir == "FCR") %>%
      dplyr::filter(Site == 100) %>% #100 = weir inflow site
      dplyr::select(DateTime, DRSI_mgL) %>%
      dplyr::rename(time = DateTime) %>%
      dplyr::mutate(time = lubridate::force_tz(time, tzone = input_file_tz),
                    time = lubridate::with_tz(time, tzone = "UTC"),
                    time = lubridate::as_date(time))

    all_data <- dplyr::left_join(inflow_combined, FCRchem, by = "time")

    ghg <- readr::read_csv(co2_ch4, show_col_types = FALSE) %>%
      dplyr::filter(Reservoir == "FCR") %>%
      dplyr::filter(Site == 100) %>% #weir inflow
      dplyr::select(DateTime, CH4_umolL) %>%
      dplyr::rename(time = DateTime,
                    CAR_ch4 = CH4_umolL) %>%
      dplyr::mutate(time = lubridate::force_tz(time, tzone = input_file_tz),
                    time = lubridate::with_tz(time, tzone = "UTC"),
                    time = lubridate::as_date(time)) %>%
      dplyr::group_by(time) %>%
      tidyr::drop_na() %>%
      dplyr::summarise(CAR_ch4 = mean(CAR_ch4)) %>%
      dplyr::filter(CAR_ch4 < 0.2) #remove outliers


    start_date<- lubridate::as_date("2015-07-07")
    end_date<- lubridate::as_date("2020-12-31")


    #start_date<-as.POSIXct(strptime("2015-07-07", "%Y-%m-%d", tz="EST"))
    #end_date<-as.POSIXct(strptime("2020-12-31", "%Y-%m-%d", tz="EST"))

    #creating new dataframe with list of all dates
    datelist<-seq(start_date,end_date, "1 day") #changed from May 15, 2013 because of NA in flow
    datelist<-tibble::tibble(time = datelist)

    ghg1 <- dplyr::left_join(datelist, ghg, by="time")
    #need to interpolate missing data, but first need to fill first & last values
    ghg1$CAR_ch4[1]<-ghg$CAR_ch4[which.min(abs(ghg$time - start_date))]
    ghg1$CAR_ch4[length(ghg1$CAR_ch4)]<-ghg$CAR_ch4[which.min(abs(ghg$time - end_date))]
    ghg1$CAR_ch4 <- zoo::na.fill(zoo::na.approx(ghg1$CAR_ch4), "extend")

    alldata<-dplyr::left_join(all_data, ghg1, by="time") %>%
      dplyr::group_by(time) %>%
      dplyr::summarise_all(mean, na.rm=TRUE)
    #merge chem with CH4 data, truncating to start and end date period of CH4 (not chem)

    lastrow <- length(alldata$time) #need for extend function below
    #now need to interpolate missing values in chem; setting 1st and last value in time series as medians
    #then linearly interpolating the middle missing values
    alldata$TN_ugL[1]<-median(na.exclude(alldata$TN_ugL))
    alldata$TN_ugL[lastrow]<-median(na.exclude(alldata$TN_ugL))
    alldata$TN_ugL<-zoo::na.fill(zoo::na.approx(alldata$TN_ugL),"extend")

    alldata$TP_ugL[1]<-median(na.exclude(alldata$TP_ugL))
    alldata$TP_ugL[lastrow]<-median(na.exclude(alldata$TP_ugL))
    alldata$TP_ugL<-zoo::na.fill(zoo::na.approx(alldata$TP_ugL),"extend")

    alldata$NH4_ugL[1]<-median(na.exclude(alldata$NH4_ugL))
    alldata$NH4_ugL[lastrow]<-median(na.exclude(alldata$NH4_ugL))
    alldata$NH4_ugL<-zoo::na.fill(zoo::na.approx(alldata$NH4_ugL),"extend")

    alldata$NO3NO2_ugL[1]<-median(na.exclude(alldata$NO3NO2_ugL))
    alldata$NO3NO2_ugL[lastrow]<-median(na.exclude(alldata$NO3NO2_ugL))
    alldata$NO3NO2_ugL<-zoo::na.fill(zoo::na.approx(alldata$NO3NO2_ugL),"extend")

    alldata$SRP_ugL[1]<-median(na.exclude(alldata$SRP_ugL))
    alldata$SRP_ugL[lastrow]<-median(na.exclude(alldata$SRP_ugL))
    alldata$SRP_ugL<-zoo::na.fill(zoo::na.approx(alldata$SRP_ugL),"extend")

    alldata$DOC_mgL[1]<-median(na.exclude(alldata$DOC_mgL))
    alldata$DOC_mgL[lastrow]<-median(na.exclude(alldata$DOC_mgL))
    alldata$DOC_mgL<-zoo::na.fill(zoo::na.approx(alldata$DOC_mgL),"extend")

    alldata$DIC_mgL[1]<-median(na.exclude(alldata$DIC_mgL))
    alldata$DIC_mgL[lastrow]<-median(na.exclude(alldata$DIC_mgL))
    alldata$DIC_mgL<-zoo::na.fill(zoo::na.approx(alldata$DIC_mgL),"extend")

    alldata <- alldata[(!duplicated(alldata$time)),] #remove duplicated dates

    #need to convert mass observed data into mmol/m3 units for two pools of organic carbon
    weir_inflow <- alldata %>%
      dplyr::select(-c(Depth_m, Rep)) %>%
      dplyr::mutate(NIT_amm = NH4_ugL*1000*0.001*(1/18.04),
                    NIT_nit = NO3NO2_ugL*1000*0.001*(1/62.00), #as all NO2 is converted to NO3
                    PHS_frp = SRP_ugL*1000*0.001*(1/94.9714),
                    OGM_doc = DOC_mgL*1000*(1/12.01)* 0.10,  #assuming 10% of total DOC is in labile DOC pool (Wetzel page 753)
                    OGM_docr = 1.5*DOC_mgL*1000*(1/12.01)* 0.90, #assuming 90% of total DOC is in recalcitrant DOC pool
                    TN_ugL = TN_ugL*1000*0.001*(1/14),
                    TP_ugL = TP_ugL*1000*0.001*(1/30.97),
                    OGM_poc = 0.1*(OGM_doc+OGM_docr), #assuming that 10% of DOC is POC (Wetzel page 755
                    OGM_don = (5/6)*(TN_ugL-(NIT_amm+NIT_nit))*0.10, #DON is ~5x greater than PON (Wetzel page 220)
                    OGM_donr = (5/6)*(TN_ugL-(NIT_amm+NIT_nit))*0.90, #to keep mass balance with DOC, DONr is 90% of total DON
                    OGM_pon = (1/6)*(TN_ugL-(NIT_amm+NIT_nit)), #detemined by subtraction
                    OGM_dop = 0.3*(TP_ugL-PHS_frp)*0.10, #Wetzel page 241, 70% of total organic P = particulate organic; 30% = dissolved organic P
                    OGM_dopr = 0.3*(TP_ugL-PHS_frp)*0.90,#to keep mass balance with DOC & DON, DOPr is 90% of total DOP
                    OGM_pop = TP_ugL-(OGM_dop+OGM_dopr+PHS_frp), # #In lieu of having the adsorbed P pool activated in the model, need to have higher complexed P
                    CAR_dic = DIC_mgL*1000*(1/52.515),
                    OXY_oxy = rMR::Eq.Ox.conc(TEMP, elevation.m = 506, #creating OXY_oxy column using RMR package, assuming that oxygen is at 100% saturation in this very well-mixed stream
                                              bar.press = NULL, bar.units = NULL,
                                              out.DO.meas = "mg/L",
                                              salinity = 0, salinity.units = "pp.thou"),
                    OXY_oxy = OXY_oxy *1000*(1/32),
                    PHY_cyano = 0,
                    PHY_green = 0,
                    PHY_diatom = 0,
                    SIL_rsi = median(silica$DRSI_mgL),
                    SIL_rsi = SIL_rsi*1000*(1/60.08),
                    SALT = 0,
                    CAR_ch4 = 0) %>%
      dplyr::mutate_if(is.numeric, round, 4) #round to 4 digits

    #Long-term median pH of FCR is 6.5, at which point CO2/HCO3 is about 50-50

    #given this disparity, using a 50-50 weighted molecular weight (44.01 g/mol and 61.02 g/mol, respectively)


    #reality check of mass balance
    # hist(weir_inflow$TP_ugL - (weir_inflow$PHS_frp + weir_inflow$OGM_dop + weir_inflow$OGM_dopr + weir_inflow$OGM_pop))
    # hist(weir_inflow$TN_ugL - (weir_inflow$NIT_amm + weir_inflow$NIT_nit + weir_inflow$OGM_don + weir_inflow$OGM_donr + weir_inflow$OGM_pon))
    VARS <- c("time", "FLOW", "TEMP", "SALT",
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

    #clean it up and get vars in order
    inflow_clean <- weir_inflow %>%
      dplyr::select(dplyr::all_of(VARS))
  }else{
    inflow_clean <- inflow_combined
  }

  #inflow_clean %>%
  #   pivot_longer(cols = -time, names_to = "Nutrient", values_to = "values") %>%
  #   ggplot(aes(x = time, y = values)) +
  #   geom_point() +
  #   geom_line() +
  #   facet_wrap(~Nutrient, scales = "free")

  dir.create(dirname(cleaned_inflow_file), recursive = TRUE, showWarnings = FALSE)

  inflow_clean |>
    tidyr::pivot_longer(-time, names_to = "variable", values_to = "observation") |>
    dplyr::rename(datetime = time) |>
    dplyr::mutate(site_id = site_id) |>
    dplyr::select(datetime, site_id, variable, observation) |>
    readr::write_csv(cleaned_inflow_file)

  return(cleaned_inflow_file)

}


