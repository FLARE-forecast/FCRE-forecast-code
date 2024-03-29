#remotes::install_github("rqthomas/cronR")
#remotes::install_deps()
library(cronR)

lake_directory <- here::here()
config_set_name <- "defaultV2"

home_dir <- file.path(lake_directory, "workflows", config_set_name)

cmd <- cronR::cron_rscript(rscript = file.path(home_dir, "combined_workflow.R"),
                           rscript_log = file.path(home_dir, "fcre.log"),
                           log_append = FALSE,
                           #cmd = "/usr/local/bin/r", # use litter, more robust on CLI
                           workdir = file.path(home_dir))
#trailing_arg = "curl -fsS -m 10 --retry 5 -o /dev/null https://hc-ping.com/cb249e47-f56b-45da-af7f-9c0c47db1a6c")
cronR::cron_add(command = cmd,  frequency = 'hourly', id = 'fcre_forecast',ask = FALSE)

