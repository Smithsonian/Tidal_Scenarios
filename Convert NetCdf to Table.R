# AR6 SLR projections
library(tidyverse)
library(ncdf4)
library(dplyr)
library(tidyr)
library(tibble)

psmsl_connector <- read_csv("data/npsset_data.csv")

# ------------------------------------------------------------
# Function: Convert an IPCC AR6 regional NetCDF to tidy format
#           while retaining tide-gauge projections only.
# ------------------------------------------------------------

tidy_ar6_nc <- function(nc_file) {
  
  nc <- nc_open(nc_file)
  on.exit(nc_close(nc))
  
  # ----------------------------------------------------------
  # Read dimensions / coordinates
  # ----------------------------------------------------------
  
  location_id <- ncvar_get(nc, "locations")
  lat         <- ncvar_get(nc, "lat")
  lon         <- ncvar_get(nc, "lon")
  year        <- ncvar_get(nc, "years")
  quantile    <- ncvar_get(nc, "quantiles")
  
  # ----------------------------------------------------------
  # Determine whether this is a "values" or "rates" file
  # ----------------------------------------------------------
  
  if ("sea_level_change" %in% names(nc$var)) {
    
    variable_name <- "sea_level_change"
    value_name    <- "sea_level_change_mm"
    
  } else if ("sea_level_change_rate" %in% names(nc$var)) {
    
    variable_name <- "sea_level_change_rate"
    value_name    <- "sea_level_change_rate_mm_yr"
    
  } else {
    
    stop(
      "Could not find sea_level_change or sea_level_change_rate."
    )
  }
  
  # ----------------------------------------------------------
  # Read projection array
  #
  # Documentation says R reads dimensions as:
  # locations x years x quantiles
  # ----------------------------------------------------------
  
  x <- ncvar_get(nc, variable_name)
  
  # ----------------------------------------------------------
  # Identify tide-gauge locations only
  #
  # location IDs < 1e9 = tide gauges
  # location IDs > 1e9 = 1x1 degree global grid
  # ----------------------------------------------------------
  
  gauge_index <- which(location_id < 1e9)
  
  # Restrict array to gauges before making the long table.
  # This also avoids unnecessarily expanding ~66,000 grid
  # locations into a very large dataframe.
  x <- x[gauge_index, , , drop = FALSE]
  
  location_id <- location_id[gauge_index]
  lat         <- lat[gauge_index]
  lon         <- lon[gauge_index]
  
  # ----------------------------------------------------------
  # Convert array to tidy form
  # ----------------------------------------------------------
  
  tidy <- expand_grid(
    location_index = seq_along(location_id),
    year_index     = seq_along(year),
    quantile_index = seq_along(quantile)
  ) %>%
    mutate(
      location_id = location_id[location_index],
      latitude    = lat[location_index],
      longitude   = lon[location_index],
      year         = year[year_index],
      quantile     = quantile[quantile_index],
      
      value = x[cbind(
        location_index,
        year_index,
        quantile_index
      )]
    ) %>%
    select(
      location_id,
      latitude,
      longitude,
      year,
      quantile,
      value
    ) %>% 
    rename(psmsl_id = location_id)
  
  names(tidy)[names(tidy) == "value"] <- value_name
  
  tidy
}

medium_ncs <- list.files("data/ar6-regional-confidence/regional/confidence_output_files/medium_confidence/",
                      pattern = "^total.*values\\.nc$",
                      full.names = T,
                      recursive = T
                      )

nc_list <- list()
for (i in 1:length(medium_ncs)) {
  
  scenario_name <- str_split(medium_ncs[i], "\\/")[[1]][7]
  
  tab_version <- tidy_ar6_nc(medium_ncs[i]) %>% 
    left_join(psmsl_connector) %>% 
    filter(complete.cases(noaa_id)) %>% 
    mutate(scenario_name = scenario_name) %>% 
    mutate(confidence = "medium")
  
  nc_list[[i]] <- tab_version
  
}

low_ncs <- list.files("data/ar6-regional-confidence/regional/confidence_output_files/low_confidence/",
                         pattern = "^total.*values\\.nc$",
                         full.names = T,
                         recursive = T
)

list_index <- (length(medium_ncs)+1):(length(medium_ncs)+length(low_ncs))

for (i in 1:length(length(low_ncs))) {
  
  scenario_name <- str_split(low_ncs[i], "\\/")[[1]][7]
  
  tab_version <- tidy_ar6_nc(low_ncs[i]) %>% 
    left_join(psmsl_connector) %>% 
    filter(complete.cases(noaa_id)) %>% 
    mutate(scenario_name = scenario_name) %>% 
    mutate(confidence = "low")
  
  nc_list[[list_index[i]]] <- tab_version
  
}

all_slr_scenarios <- bind_rows(nc_list)

arrow::write_parquet(all_slr_scenarios, "data/ar6_us_compiled_senarios.parquet")

test_plot <- all_slr_scenarios %>% 
  filter(noaa_name == "Annapolis") %>% 
  arrange(year, quantile)
library(plotly)
plot_this <- ggplot(test_plot, aes(x = year, y = sea_level_change_mm)) +
  geom_line(aes(group = quantile, color = scenario_name), alpha = 0.3) +
  facet_wrap(scenario_name~confidence)

ggplotly(plot_this)
