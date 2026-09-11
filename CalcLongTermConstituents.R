# Query and calc long term nodal cycles for all NOAA psmsl gauges

library(tidyverse)
library(VulnToolkit)
library(zoo)

npsset_data <- read_csv("https://raw.githubusercontent.com/laura-feher/NPSETr/refs/heads/master/data/future_slr_projections_sweet_2022.csv") %>% 
  select(`PSMSL ID`, `NOAA ID`, `NOAA Name`) %>% 
  rename(psmsl_id = `PSMSL ID`,
         noaa_id = `NOAA ID`,
         noaa_name = `NOAA Name`) %>% 
  distinct_all() %>% 
  filter(noaa_id!=0)

pb <- txtProgressBar(min = 1, max = nrow(npsset_data), style = 3)
for (i in 1:nrow(npsset_data)) {
  
  noaa_id <- npsset_data$noaa_id[i]
  
  raw_csv <- paste0('output/tabs/raw_wl/',
                       noaa_id,
                       ".parquet"
                       )
  
  hl_csv <- paste0('output/tabs/classified_hl/',
                   noaa_id,
                   ".csv"
  )
  
  datum_csv <- paste0('output/tabs/datums/',
                   noaa_id,
                   ".csv"
  )
  
  # int_csv <- paste0('output/tabs/annual_summaries/',
  #                   noaa_id,
  #                   ".csv"
  # )
  
  final_csv <- paste0('output/tabs/nodal_parameters/',
                    noaa_id,
                    ".csv"
  )
  
  vis_fig <- paste0('output/figs/nodal_cycle/',
                      noaa_id,
                      ".jpg"
  )
  
  
  if (file.exists(final_csv)) {
    print(paste0(npsset_data$noaa_name[i], " already done."))
    
  } else {
    
    print(paste0("Analysing ", npsset_data$noaa_name[i], "..."))
    
    try({
      
      params <- VulnToolkit::noaa.parameters(noaa_id) %>% 
        filter(params == "Verified High/Low Water Level") %>% 
        mutate(startDate = format(lubridate::ymd_hm(startDate), format = "%Y%m%d"),
               endDate = format(lubridate::ymd_hm(endDate), format = "%Y%m%d")
               )
      
      if (nrow(params) != 0) {
        
        if (file.exists(raw_csv)) {
          temp_hl_data <- arrow::read_parquet(raw_csv)
        } else {
          
          print(paste0("... downloading "))
          
          temp_hl_data <- VulnToolkit::noaa(begindate = min(params$startDate), 
                                            enddate = min("20251231",
                                                          max(params$endDate)
                                            ),
                                            station = noaa_id,
                                            interval = "hourly",
                                            continuous = T,
                                            datum = "MSL",
                                            units = "meters"
          )
          
          names(temp_hl_data) <- c("dateTime", "waterLevel", "station")
          
          arrow::write_parquet(temp_hl_data, raw_csv)
        } 
        
        
        if (file.exists(hl_csv)&file.exists(datum_csv)) {
          
          hl <- read_csv(hl_csv)
          datums <- read_csv(datum_csv)
        
        } else {
          
          temp_hl_data_2 <- temp_hl_data %>% 
            mutate(year = lubridate::year(dateTime))
          
          years_of_record <- unique(temp_hl_data_2$year)
          
          hl_list <- list()
          datum_list <- list()
          
          for (j in 2:(length(years_of_record)-1)) {
            
            this_year <- years_of_record[j]
            print(paste0("...", this_year))
            
            if (nrow(
              temp_hl_data_2 %>% filter(year == this_year & !is.na(waterLevel))
              ) < (8736*0.75)) {
              next
            }
            
            custom_datums <- fitCustomTidalDatum(wlTable=temp_hl_data_2,
                                                 startDate = paste0(this_year, "-01-01 00:00"),
                                                 endDate =  paste0(this_year, "-12-31 23:59"),
                                                 graph = T,
                                                 gauge_data = "output/figs/year_prediction_plots/",
                                                 out_fig_name = paste0(noaa_id, "_", this_year)
                                                 )
            
            hl_list[[j]] <- custom_datums[[1]] %>% 
              mutate(station_name = npsset_data$noaa_name[i])
            
            datum_list[[j]] <- custom_datums[[2]] %>% 
              mutate(year = this_year,
                     station_id = noaa_id,
                     station_name = npsset_data$noaa_name[i]
              )
          
        }  # end of j years loop
          
          hl <- bind_rows(hl_list)
          write_csv(hl, hl_csv)
          
          datums <- bind_rows(datum_list)
          write_csv(datums, datum_csv)
          
        } # end of if hl and datum files exist
        floods <- datums %>% 
          filter(! Datum %in% c("MSL", "LOT", "HOT"))
        
        MTL <- datums %>% 
          select(-c(n_pred, n_obs, fallingTime, risingTime)) %>% 
          filter(Datum %in% c("MSL")) %>% 
          pivot_wider(names_from = "Datum",
                      values_from = c("observed", "predicted"))

        amplitudes <- floods %>% 
          left_join(MTL) %>% 
          mutate(wl = observed - observed_MSL) %>% 
          mutate(station_id = noaa_id)
        
        output_fig <- ggplot(amplitudes, aes(x = year, y = wl)) +
          geom_point(aes(color = Datum)) +
          ggtitle(paste0(noaa_id, ": ", npsset_data$noaa_name[i])) +
          theme_minimal() +
          theme(legend.title = element_blank())
        
        (output_fig)
        
        ggsave(vis_fig, output_fig, width = 4.25, height = 3)
        
        P18 <- 18.61
        omega18 <- 2 * pi / P18
        P44 <- 4.4
        omega44 <- 2 * pi / P44
        
        amplitudes <- amplitudes |>
          dplyr::mutate(
            sin18 = sin(omega18 * year),
            cos18 = cos(omega18 * year),
            sin44 = sin(omega44 * year),
            cos44 = cos(omega44 * year)
          )
        
        unique_floods <- unique(amplitudes$Datum)
        
        n_floods <- length(unique_floods)
        
        output_df <- data.frame(station_id = rep(noaa_id, n_floods), 
                                station_name = rep(npsset_data$noaa_name[i], n_floods),
                                start_date = rep(min(params$startDate), n_floods), 
                                end_date = rep(min("20251231",max(params$endDate), n_floods)),
                                tide = unique_floods, 
                                # form_factor = rep(temp_f$form.number, n_floods),
                                mean_amp = rep(NA, n_floods),
                                amp18 = rep(NA, n_floods),
                                phase18 = rep(NA, n_floods),
                                amp44 = rep(NA, n_floods),
                                phase44 = rep(NA, n_floods),
                                rse = rep(NA, n_floods),
                                r2 = rep(NA, n_floods),
                                mean_ampb = rep(NA, n_floods),
                                amp18b = rep(NA, n_floods),
                                phase18b = rep(NA, n_floods),
                                rseb = rep(NA, n_floods),
                                r2b = rep(NA, n_floods))
        
        for (j in 1:length(unique_floods)) {
          
          temp_flood <- unique_floods[j]
          
          temp_amps <- amplitudes %>% filter(Datum == temp_flood)
          
          fit <- lm(
            wl ~
              sin18 + cos18 +
              sin44 + cos44,
            data = temp_amps
          )
          
          fit2 <- lm(
            wl ~
              sin18 + cos18,
            data = temp_amps
          )
          
          b <- coef(fit)
          rse <- summary(fit)$sigma
          r2 <- summary(fit)$adj.r.squared
          # p_values <- broom::tidy(fit)
          
          amp18 <- as.numeric(sqrt(b["sin18"]^2 + b["cos18"]^2))
          phase18 <- as.numeric(atan2(-b["cos18"], b["sin18"]) * 18.61 / (2 * pi))
          
          amp44<- as.numeric(sqrt(b["sin44"]^2 + b["cos44"]^2))
          phase44 <- as.numeric(atan2(-b["cos44"], b["sin44"]) * 4.4 / (2 * pi))
          
          b2 <- coef(fit2)
          rseb <- summary(fit2)$sigma
          r2b <- summary(fit2)$adj.r.squared
          
          amp18b <- as.numeric(sqrt(b2["sin18"]^2 + b2["cos18"]^2))
          phase18b <- as.numeric(atan2(-b2["cos18"], b2["sin18"]) * 18.61 / (2 * pi))
          
          output_df$mean_amp[j] <- b["(Intercept)"]
          output_df$amp18[j] <- amp18
          output_df$phase18[j] <- phase18
          output_df$amp44[j] <- amp44
          output_df$phase44[j] <- phase44
          output_df$rse[j] <- rse
          output_df$r2[j] <- r2
          output_df$mean_ampb[j] <- b2["(Intercept)"]
          output_df$amp18b[j] <- amp18b
          output_df$phase18b[j] <- phase18b
          output_df$rseb[j] <- rseb
          output_df$r2b[j] <- r2b
          
        } # end of unique floods iteration
        
        write_csv(output_df, final_csv)
        
      } # end of if nrow params == 0
    }) # end try statement
    
    
  } # end if else already done
  
  setTxtProgressBar(pb, i)
  
} # end iteration

close(pb)


output_files <- list.files("output/tabs/nodal_parameters/",
                          full.names = T
                          )
output_list <- list()
for (i in 1:length(output_files)) {
  output_list[[i]] <- read_csv(output_files[i])
}

output_tab <- bind_rows(output_list)

write_csv(output_tab, "data/long_term_tidal_constituents.csv")


output_files <- list.files("output/tabs/datums/",
                           full.names = T
)
output_list <- list()
for (i in 1:length(output_files)) {
  output_list[[i]] <- read_csv(output_files[i])
}

output_tab <- bind_rows(output_list)

write_csv(output_tab, "data/annual_compiled_datums.csv")

