#renv::restore()

library(tidyverse)
library(lubridate)

message("Beginning generate targets")
message(config_set_name)

#' Set the lake directory to the repository directory

lake_directory <- here::here()


#' Source the R files in the repository

#files.sources <- list.files(file.path(lake_directory, "R"), full.names = TRUE)
#sapply(files.sources, source)

#' Generate the `config_obs` object and create directories if necessary

#config_obs <- FLAREr::initialize_obs_processing(lake_directory, observation_yml = "observation_processing.yml", config_set_name)
#config <- FLAREr::set_configuration(configure_run_file,lake_directory, config_set_name = config_set_name)


#' Clone or pull from data repositories

FLAREr:::get_git_repo(lake_directory,
                     directory = config_obs$realtime_insitu_location,
                     git_repo = "https://github.com/FLARE-forecast/FCRE-data.git")

FLAREr:::get_git_repo(lake_directory,
                     directory = config_obs$realtime_met_station_location,
                     git_repo = "https://github.com/FLARE-forecast/FCRE-data.git")

FLAREr:::get_git_repo(lake_directory,
                     directory = config_obs$realtime_inflow_data_location,
                     git_repo = "https://github.com/FLARE-forecast/FCRE-data.git")

#get_git_repo(lake_directory,
#             directory = config_obs$manual_data_location,
#             git_repo = "https://github.com/FLARE-forecast/FCRE-data.git")

#' Download files from EDI

get_edi_file(edi_https = "https://pasta.lternet.edu/package/data/eml/edi/389/7/02d36541de9088f2dd99d79dc3a7a853",
                     file = config_obs$met_raw_obs_fname[2],
                     lake_directory)

get_edi_file(edi_https = "https://pasta.lternet.edu/package/data/eml/edi/271/7/71e6b946b751aa1b966ab5653b01077f",
                     file = config_obs$insitu_obs_fname[2],
                     lake_directory)

get_edi_file(edi_https = "https://pasta.lternet.edu/package/data/eml/edi/198/11/81f396b3e910d3359907b7264e689052",
                     file = config_obs$secchi_fname,
                     lake_directory)

get_edi_file(edi_https = "https://pasta.lternet.edu/package/data/eml/edi/200/13/27ceda6bc7fdec2e7d79a6e4fe16ffdf",
                     file = config_obs$ctd_fname,
                     lake_directory)

get_edi_file(edi_https = "https://pasta.lternet.edu/package/data/eml/edi/199/11/509f39850b6f95628d10889d66885b76",
                     file = config_obs$nutrients_fname,
                     lake_directory)


get_edi_file(edi_https = "https://pasta.lternet.edu/package/data/eml/edi/202/10/c065ff822e73c747f378efe47f5af12b",
                     file = config_obs$inflow_raw_file1[2],
                     lake_directory)

get_edi_file(edi_https = "https://pasta.lternet.edu/package/data/eml/edi/542/1/791ec9ca0f1cb9361fa6a03fae8dfc95",
                     file = "silica_master_df.csv",
                     lake_directory)

get_edi_file(edi_https = "https://pasta.lternet.edu/package/data/eml/edi/551/7/38d72673295864956cccd6bbba99a1a3",
                     file = "Dissolved_CO2_CH4_Virginia_Reservoirs.csv",
                     lake_directory)

message("Clean up observed meterology")

cleaned_met_file <- met_data_bind(realtime_file = file.path(lake_directory,'data_raw',config_obs$met_raw_obs_fname[1]),
                                 qaqc_file = file.path(lake_directory,'data_raw', config_obs$met_raw_obs_fname[2]),
                                 cleaned_met_file = file.path(lake_directory,'targets', config_obs$site_id,paste0("observed-met_",config_obs$site_id,".csv")),
                                 input_file_tz = "EST",
                                 nldas = NULL,
                                 site_id = config_obs$site_id)

# cleaned_met_file <- met_qaqc_csv(realtime_file = file.path(config_obs$file_path$data_directory, config_obs$met_raw_obs_fname[1]),
#                                  qaqc_file = file.path(config_obs$file_path$data_directory, config_obs$met_raw_obs_fname[2]),
#                                  cleaned_met_file = file.path(config_obs$file_path$targets_directory, config_obs$site_id,paste0("observed-met_",config_obs$site_id,".csv")),
#                                  input_file_tz = "EST",
#                                  nldas = NULL,
#                                  site_id = config_obs$site_id)

message("Clean up observed inflow")

cleaned_inflow_file <- inflow_data_combine(realtime_file = file.path(lake_directory,'data_raw', config_obs$inflow_raw_file1[1]),
                                       qaqc_file = file.path(lake_directory,'data_raw', config_obs$inflow_raw_file1[2]),
                                       nutrients_file = file.path(lake_directory,'data_raw', config_obs$nutrients_fname),
                                       silica_file = file.path(lake_directory,'data_raw',  config_obs$silica_fname),
                                       co2_ch4 = file.path(lake_directory,'data_raw', config_obs$ch4_fname),
                                       cleaned_inflow_file = file.path(lake_directory,'targets', config_obs$site_id, paste0(config_obs$site_id,"-targets-inflow.csv")),
                                       input_file_tz = 'EST',
                                       site_id = config_obs$site_id)

message("Clean up observed insitu measurements")

cleaned_insitu_file <- in_situ_qaqc_csv(insitu_obs_fname = file.path(lake_directory,'data_raw',config_obs$insitu_obs_fname),
                                        data_location = file.path(lake_directory,'data_raw'),
                                        maintenance_file = file.path(lake_directory,'data_raw',config_obs$maintenance_file),
                                        ctd_fname = file.path(lake_directory,'data_raw', config_obs$ctd_fname),
                                        nutrients_fname =  file.path(lake_directory,'data_raw', config_obs$nutrients_fname),
                                        secchi_fname = file.path(lake_directory,'data_raw', config_obs$secchi_fname),
                                        ch4_fname = file.path(lake_directory,'data_raw', config_obs$ch4_fname),
                                        cleaned_insitu_file = file.path(lake_directory,'targets', config_obs$site_id, paste0(config_obs$site_id,"-targets-insitu.csv")),
                                        lake_name_code = config_obs$site_id,
                                        config = config_obs)
#' Move targets to s3 bucket

message("Successfully generated targets")

FLAREr:::put_targets(site_id = config_obs$site_id,
                    cleaned_insitu_file,
                    cleaned_met_file,
                    cleaned_inflow_file,
                    use_s3 = config$run_config$use_s3,
                    config = config)

if(config$run_config$use_s3){
  message("Successfully moved targets to s3 bucket")
}

