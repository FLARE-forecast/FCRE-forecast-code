extract_secchi <- function(fname,
                        input_file_tz,
                        focal_depths){

  d <- read_csv(fname,
                col_types = readr::cols()) %>%
    filter(Reservoir == "FCR" & Site == 50) %>%
    select(DateTime, Secchi_m) %>%
    mutate(DateTime = force_tz(DateTime, input_file_tz)) %>%
    mutate(DateTime = lubridate::with_tz(DateTime, tzone = "UTC")) %>%
    group_by(DateTime) %>%
    summarise(secchi = mean(Secchi_m, na.rm = TRUE), .groups = 'drop') %>%
    rename("time" = DateTime) %>%
    pivot_longer(cols = -c(time), names_to = "variable", values_to = "observed") %>%
    mutate(depth = NA) %>%
    filter(!is.na(observed)) %>%
    select(time , depth, observed, variable) %>%
    distinct()

  return(d)
}
