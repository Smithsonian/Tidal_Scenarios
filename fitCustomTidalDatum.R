library(slider)
#Our transformation function
scaleFUN <- function(x) sprintf("%.1f", x)

assign_lunar_day <- function(times, new_day_threshold = 23) {
  
  n <- length(times)
  bin_id <- integer(n)
  bin_id[1] <- 1
  bin_start <- times[1]
  count_i <- 1
  
  for (i in seq_len(n)[-1]) {
    elapsed <- as.numeric(difftime(times[i], bin_start, units = "hours"))
    if (elapsed >= new_day_threshold | count_i == 4) {
      bin_id[i] <- bin_id[i - 1] + 1
      bin_start <- times[i]
      count_i <- 1
    } else {
      bin_id[i] <- bin_id[i - 1]
      count_i <- count_i+1 
    }
  }
  bin_id
}

fitCustomTidalDatum <- function(wlTable, startDate='2015-01-01 00:00',
                                bufferStart = 15,
                                bufferEnd = 15,
                                endDate='2015-12-31 23:59', graph=F, out_fig_name="temp",
                                gauge_data="biomass/agb_biomass_inundation/NOAA_water_levels") {
  

  bufferStartDate <- ymd_hm(startDate) - days(bufferStart)
  bufferStartDate <- toString(format(bufferStartDate, "%Y-%m-%d %H:%M"))
  
  bufferEndDate <- ymd_hm(endDate) + days(bufferEnd)
  bufferEndDate <- toString(format(bufferEndDate, "%Y-%m-%d %H:%M"))
  
  wlTable <- wlTable %>%
    filter(dateTime >= bufferStartDate & dateTime <= bufferEndDate)
  
  # Fit harmonics
  fitHamonics <- TideHarmonics::ftide(x=wlTable$waterLevel,
                                      dto = wlTable$dateTime,
                                      nodal = F)
  
  spring_neap_halfwidth <- 6
  
  # Amplitudes are usually in fit$Amp (named vector) or a coefficient data frame
  # is_semidiurnal <- ifelse(fitHamonics$fval <= 3, 1, 0)
  
  min_timestep <- min(
    diff(sort(wlTable$dateTime[!is.na(wlTable$dateTime)]))
  )
  
  predictFromHarmonics <- predict(fitHamonics, from = min(wlTable$dateTime), to = max(wlTable$dateTime),
                                  by = as.numeric(min_timestep, units = "mins")/60,
                                  msl = F) + mean(wlTable$waterLevel, na.rm=T)
  
  
  # detrend and add msl, IDK why but sometimes predict from harmonics is relative to 0, sometimes to MSL
  # predictFromHarmonics <- predictFromHarmonics - mean(predictFromHarmonics)
  # predictFromHarmonics <- predictFromHarmonics + fitHamonics$msl
  predicted_wls <- data.frame(dateTime = seq(min(wlTable$dateTime), 
                                             max(wlTable$dateTime), 
                                             by =  as.numeric(min_timestep, units = "mins")/60 * 3600),
                              predicted = predictFromHarmonics)
  
  allWaterLevels <- wlTable %>%
    full_join(predicted_wls, by="dateTime") %>% 
    rename(observed = waterLevel) %>% 
    arrange(dateTime)
  
  classifiedTides <- allWaterLevels %>%
    mutate(HL = ifelse(predicted > lag(predicted) & predicted > lead(predicted), "H",
                       ifelse(predicted < lag(predicted) & predicted < lead(predicted), "L", NA))) %>%
    filter(! is.na(HL))
  
  # quarter_lunar_day <- 24.8412 / 4  # ~6.2103 hours
  
  classifiedTides <- classifiedTides %>%
    arrange(dateTime) %>%
    mutate(year = year(dateTime),
           month = month(dateTime),
           day = day(dateTime)) %>% 
    mutate(lunar_day_id = assign_lunar_day(dateTime)) %>%
    group_by(lunar_day_id, HL) %>%
    mutate(n_in_group = n(),
      classifiedTides = case_when(
        n_in_group == 1 & HL == "H" ~ "HH",
        n_in_group == 1 & HL == "L" ~ "LL",
        n_in_group == 2 & HL == "H" & predicted == max(predicted) ~ "HH",
        n_in_group == 2 & HL == "H" ~ "LH",
        n_in_group == 2 & HL == "L" & predicted == min(predicted) ~ "LL",
        n_in_group == 2 & HL == "L" ~ "HL",
        n_in_group > 2 ~ NA_character_,  # flag for manual review rather than guess
        TRUE ~ NA_character_
      )) %>%
    ungroup() %>%
    # select(-n_in_group)
    mutate(
      localMaxTR = slide_index_dbl(
        predicted, lunar_day_id, max, na.rm = TRUE,
        .before = spring_neap_halfwidth,
        .after  = spring_neap_halfwidth,
        .complete = TRUE
      ),
      localMinTR = slide_index_dbl(
        predicted, lunar_day_id, min, na.rm = TRUE,
        .before = spring_neap_halfwidth,
        .after  = spring_neap_halfwidth,
        .complete = TRUE
      ),flagMax = case_when(
        predicted == localMaxTR ~ "monthlyMax",
        predicted == localMinTR ~ "monthlyMin",
        .default = "none"
      )
    ) %>%
    filter((dateTime >= ymd_hm(startDate, tz="UTC")) & (dateTime <= ymd_hm(endDate, tz="UTC"))) %>% 
    group_by() %>% 
    mutate(highestHigh = max(predicted),
           lowestLow = min(predicted)
           ) %>% 
    ungroup() %>% 
    mutate(flagMax2 = case_when(
      predicted == highestHigh ~ "yearlyMax",
      predicted == lowestLow ~ "yearlyMin",
      .default = "none"
    )) %>% 
    mutate(classifiedTides = case_when(flagMax2 == "yearlyMax" ~ paste0(classifiedTides, "A"),
                                       flagMax2 == "yearlyMin" ~ paste0(classifiedTides, "A"),
                                       flagMax == "monthlyMax" ~ paste0(classifiedTides, "S"),
                                       flagMax == "monthlyMin" ~ paste0(classifiedTides, "S"),
                                       .default = classifiedTides),
           risingTime = difftime(dateTime, lag(dateTime), units = "hours"),
           fallingTime = difftime(dateTime, lead(dateTime), units = "hours")
           )

  # simpleHL <- classifiedTides %>%
  #   group_by(HL) %>%
  #   summarise(wl=mean(predicted),
  #             n=n()) %>%
  #   mutate(HL = recode(HL, "H"="MHW", "L"="MLW")) %>%
  #   rename(classifiedTides = HL)
  
  datum_summaries <- classifiedTides %>%
    group_by(classifiedTides) %>%
    summarise(n_obs = sum(!is.na(observed)),
              observed=mean(observed),
              n_pred = sum(!is.na(predicted)),
              predicted = mean(predicted),
              risingTime = mean(risingTime, na.rm = T),
              fallingTime = mean(fallingTime, na.rm = T)
              ) %>%
    # bind_rows(simpleHL) %>% 
    filter(complete.cases(.)) %>%
    rename(Datum = classifiedTides)
    #mutate(Datum = recode(Datum, "LH"="MLHW", "HH"="MHHW", "HHS"="MHHWS", "HL"="MHLW", "LL"="MLLW", "LLS"="MLLWS", "L" = "MLW", "H"="MHW")
    # )
  
  wlTableSubset <- allWaterLevels %>% 
    filter((dateTime >= ymd_hm(startDate, tz="UTC")) & (dateTime <= ymd_hm(endDate, tz="UTC")))
  
  otherDatums1 <- as_tibble(data.frame(Datum = c("MSL"),
                                       predicted=c(mean(wlTableSubset$predicted,na.rm=T)),
                                       observed = c(mean(wlTableSubset$observed,na.rm=T)),
                                       n_pred=c(sum(!is.na(wlTableSubset$predicted))),
                                       n_obs = c(sum(!is.na(wlTableSubset$observed))),
                                       stringsAsFactors = F))
  
  datum_summaries <- datum_summaries %>% 
    bind_rows(otherDatums1) %>% 
    arrange(-observed)
  
  otherDatums2 <- as_tibble(data.frame(Datum = c("HOT", "LOT"), 
                                       observed=c(max(wlTableSubset$observed, na.rm=T), min(wlTableSubset$observed, na.rm=T)),
                                       n_obs=c(1, 1),
                                       stringsAsFactors = F))
  
  datum_summaries <- datum_summaries %>% 
    bind_rows(otherDatums2) %>% 
    arrange(-observed)
  
  if (graph == T) {
    # print("  ... graphing")
    # Coded datums, observed and predicted WL by fractional month
    tidesPlot <- classifiedTides %>%
      mutate(day_in_month = days_in_month(dateTime),
             f_day_of_month = ((day + (hour(dateTime)/24 + (minute(dateTime)/(60*24))))) /
               days_in_month(month)
      ) %>% 
      rename(wl = predicted) %>% 
      mutate(classifiedTides = factor(classifiedTides, levels = c("HHA", "HHS", "HA", "HS", "HH", "H", "LH",
                                              "HL", "L", "LL", "LS", "LA",  "LLS", "LLA"
                                              )))
    
    # Add predicted to wl table
    # Change waterLevel name to measured
    # Gather excluding datetime
    
    shape_values <- c(
      "HA" = 7,
      "HHA" = 3,
      "HS" = 25,
      "HHS" = 16, # filled circle
      "HH" = 17, # filled triangle
      "H" = 15, # filled square
      "LH" = 18, # filled diamond
      "HL" = 5,   # open diamond
      "L" = 0,  # open square
      "LL" = 2,  # open triangle
      "LLS" = 1,  # open circle
      "LS" = 6,
      "LLA" = 8,
      "LA" = 9
    )
    
    wlPlotting <- allWaterLevels %>% 
      dplyr::select(dateTime, observed, predicted) %>% 
      filter((dateTime >= ymd_hm(startDate, tz="UTC")) & (dateTime <= ymd_hm(endDate, tz="UTC"))) %>%
      # left_join(wlTable, by = "dateTime") %>%
      # mutate(anomaly = predicted - waterLevel) %>% 
      #select(-waterLevel) %>% 
      tidyr::gather(key = measuredOrModeled, value = wl, -dateTime) %>%
      arrange(dateTime, measuredOrModeled) %>%
      mutate(year = year(dateTime), 
             month = month(dateTime), 
             day_in_month = days_in_month(dateTime),
             f_day_of_month = ((day(dateTime) + (hour(dateTime)/24 + (minute(dateTime)/(60*24))))) /
               days_in_month(month)
      )
    
    wl_plot <- ggplot(data = wlPlotting, aes(x=f_day_of_month, y=wl)) +
      facet_wrap(year~month) +
      geom_line(aes(col = measuredOrModeled), alpha=0.66) +
      geom_point(data = tidesPlot, aes(shape = classifiedTides)) +
      ylab("Water Level (m)") +
      xlab("Month (fraction)") +
      theme_minimal() +
      scale_x_continuous(labels=scaleFUN) +
      theme(legend.title = element_blank()) + 
      scale_shape_manual(values = shape_values)
    
    out_fig <- paste(gauge_data, "/", out_fig_name, ".pdf", sep = "")
    ggsave(out_fig, wl_plot, height = 8.5, width = 11)
  }
  
  return(list(classifiedTides, datum_summaries))
  
}

getPctInundationProfile <- function(wlTable) {
  targetElevations <- seq(round(min(wlTable$waterLevel, na.rm=T),2),
                          round(max(wlTable$waterLevel, na.rm=T),2),
                          0.01)
  
  fInundation <- c()
  for (i in 1:length(targetElevations)) {
    
    inundationEvents <- filter(wlTable, waterLevel >= targetElevations[i])
    fInundation <- c(fInundation, 
                     nrow(inundationEvents)/nrow(wlTable))
  }
  
  middle_df <- data.frame(elevation = targetElevations,
                          fractionInundation = fInundation)
  
  # Lowest inundation step should be 1
  lowest_df <- data.frame(elevation =  min(targetElevations) - 0.01,
                          fractionInundation = 1)
  
  # Highest should be 0
  highest_df <- data.frame(elevation =  max(targetElevations) + 0.01,
                           fractionInundation = 0)
  
  
  output_df <- bind_rows(lowest_df, middle_df, highest_df)
  return(output_df)
  
}
