library(tidyverse)
# install.packages("arrow", repos = c("https://apache.r-universe.dev", "https://cloud.r-project.org"))


# Read the raw lines
raw_lines <- readLines("data/Kopp_2014/eft237-sup-0005-table08.tsv")

# --- Build column names from known file structure ---
# Each RCP block has 33 columns:
#   col 1:  minimum Latin Hypercube sample
#   cols 2-32: percentiles 0.1 through 99.9
#   col 33: maximum Latin Hypercube sample
# The RCP header row in the file is truncated (only 68 fields vs 100),
# so we build column names directly from the known structure instead.

pct_labels <- c(
  "0",
  "0.1", "0.5", "1.0", "2.5", "5.0", "10.0", "15.0", "16.7", "20.0",
  "25.0", "30.0", "33.3", "35.0", "40.0", "45.0", "50.0", "55.0", "60.0",
  "65.0", "66.7", "70.0", "75.0", "80.0", "83.3", "85.0", "90.0", "95.0",
  "97.5", "99.0", "99.5", "99.9",
  "100"
)  # 33 labels

rcps <- c("RCP8.5", "RCP4.5", "RCP2.6")

col_names <- paste(rep(rcps, each = length(pct_labels)),
                   rep(pct_labels, times = length(rcps)),
                   sep = "_")
# 99 column names total (3 RCPs x 33 percentiles)

# --- Parse data rows ---
data_list <- list()
current_gauge_name   <- NA_character_
current_gauge_number <- NA_integer_

for (line in raw_lines[-(1:2)]) {
  line <- trimws(line)
  if (nchar(line) == 0) next
  
  # COASTLINE rows — skip
  if (startsWith(line, "COASTLINE")) next
  
  # Gauge header rows: contain "[digits]" and no leading tab
  if (grepl("\\[\\d+\\]", line) && !startsWith(line, "\t")) {
    current_gauge_name   <- trimws(sub("\\s*\\[.*", "", line))
    current_gauge_number <- as.integer(sub(".*\\[(\\d+)\\].*", "\\1", line))
    next
  }
  
  # Data rows: first field is a year
  parts <- strsplit(line, "\t")[[1]]
  year  <- suppressWarnings(as.integer(trimws(parts[1])))
  if (is.na(year)) next
  
  values <- suppressWarnings(as.numeric(trimws(parts[-1])))
  
  if (length(values) != length(col_names)) {
    warning(sprintf("Skipping row for gauge '%s' year %d: expected %d values, got %d",
                    current_gauge_name, year, length(col_names), length(values)))
    next
  }
  
  row_df <- tibble(
    gauge_name   = current_gauge_name,
    gauge_number = current_gauge_number,
    year         = year
  )
  
  val_df <- as_tibble(setNames(as.list(values), col_names))
  
  data_list[[length(data_list) + 1]] <- bind_cols(row_df, val_df)
}

wide_df <- bind_rows(data_list)


npsset_data <- read_csv("https://raw.githubusercontent.com/laura-feher/NPSETr/refs/heads/master/data/future_slr_projections_sweet_2022.csv") %>% 
  select(`PSMSL ID`, `NOAA ID`, `NOAA Name`) %>% 
  rename(psmsl_id = `PSMSL ID`,
         noaa_id = `NOAA ID`,
         noaa_name = `NOAA Name`) %>% 
  distinct_all()

write_csv(npsset_data, "data/npsset_data.csv")

# --- Pivot to long form ---
long_df <- wide_df %>%
  pivot_longer(
    cols      = all_of(col_names),
    names_to  = "rcp_pct",
    values_to = "slr_cm",
  ) %>%
  separate(rcp_pct, into = c("rcp", "percentile"), sep = "_") %>%
  select(gauge_name, gauge_number, year, rcp, percentile, slr_cm) %>% 
  mutate(percentile = as.numeric(percentile)) %>% 
  rename(psmsl_id = gauge_number,
         psmsl_name = gauge_name) %>% 
  left_join(npsset_data, by = "psmsl_id")

# --- Preview& save ---
# write_csv(long_df, "data/Kopp_2014/Kopp_2014_projections_long.csv")
# write_rds(long_df, "data/Kopp_2014/Kopp_2014_projections_long.Rds")
arrow::write_parquet(long_df, "data/Kopp_2014/Kopp_2014_projections_long.parquet")
