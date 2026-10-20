# Min max date and location for all SET data 
library(tidyverse)
library(VulnToolkit)

npsset_data <- read_csv("https://raw.githubusercontent.com/laura-feher/NPSETr/refs/heads/master/data/future_slr_projections_sweet_2022.csv") %>% 
  select(`PSMSL ID`, `NOAA ID`, `NOAA Name`) %>% 
  rename(psmsl_id = `PSMSL ID`,
         noaa_id = `NOAA ID`,
         noaa_name = `NOAA Name`) %>% 
  distinct_all() %>% 
  filter(noaa_id!=0)

gauge_queries <- list()

for (i in 1:length(npsset_data$noaa_id)) {
  
  noaa_id <- npsset_data$noaa_id[i]
  
  params <- VulnToolkit::noaa.parameters(noaa_id) %>% 
    group_by(station) %>% 
    summarise(startDate = min(startDate),
              endDate = max(endDate))
    # mutate(startDate = format(lubridate::ymd_hm(startDate), format = "%Y%m%d"),
    #        endDate = format(lubridate::ymd_hm(endDate), format = "%Y%m%d")
  
  gauge_queries[[i]] <- params

}


url <- paste0(
  "https://api.tidesandcurrents.noaa.gov/",
  "mdapi/prod/webapi/stations.json?type=historicwl"
)


noaa_gauges <- httr::GET(url) |>
  httr::content(as = "text", encoding = "UTF-8") |>
  jsonlite::fromJSON() |>
  purrr::pluck("stations") |>
  as_tibble() |>
  transmute(
    noaa_id = id,
    noaa_name = name,
    latitude = lat,
    longitude = lng,
    state = state
  )

gauge_table <- bind_rows(gauge_queries) %>% 
  rename(noaa_id = station) %>% 
  left_join(noaa_gauges)

write_csv(gauge_table, "noaa_us_gauge_info.csv")



