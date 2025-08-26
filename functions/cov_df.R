cov_df <- function(start_date = "2001-01-01",
                   end_date = "2024-12-31",
                   stoplight_path = "R:/Fisheries Management/forecasting/2025/code_dev/all_forecasts/data/stoplight-raw-data.csv") {
  
  library(tidyverse)
  library(lubridate)
  library(curl)
  library(waterData)
  library(purrr)
  library(readr)
  
  
  # 1. River Qmax data
  station_ids <- c("12134500", "12150800", "12149000")
  station_names <- c("sky_qmax", "snoh_qmax", "snq_qmax")
  names(station_names) <- station_ids
  
  get_station_qmax <- function(station_id) {
    importDVs(station_id, code = "00060", stat = "00003",
              sdate = start_date, edate = end_date) %>%
      dplyr::mutate(
        mon = lubridate::month(dates),
        yr = lubridate::year(dates),
        Brood_Year = ifelse(mon <= 3, yr - 1, yr)
      ) %>%
      dplyr::filter(mon %in% c(1, 2, 8, 9, 10, 11, 12)) %>%
      dplyr::group_by(Brood_Year) %>%
      dplyr::summarise(Qmax = max(val, na.rm = TRUE), .groups = "drop") %>%
      dplyr::rename(!!station_names[[station_id]] := Qmax)
  }
  
  qmax_list <- purrr::map(station_ids, get_station_qmax)
  qmax_df <- purrr::reduce(qmax_list, dplyr::full_join, by = "Brood_Year")
  
  # 2. PDO
  url_pdo <- "https://psl.noaa.gov/pdo/data/pdo.timeseries.ersstv5.csv"
  pdo_data <- readr::read_csv(url_pdo, skip = 1, col_names = c("Date", "PDO"),
                              col_types = readr::cols(Date = readr::col_date(format = "%Y-%m-%d")),
                              na = c("-9999", "-9999.0", "-9999.00", "-9999.000")) %>%
    dplyr::mutate(Year = lubridate::year(Date), Month = lubridate::month(Date))
  
  pdo_summer <- pdo_data %>%
    dplyr::filter(Month >= 5 & Month <= 9) %>%
    dplyr::group_by(Year) %>%
    dplyr::summarise(PDO_Summer = mean(PDO, na.rm = TRUE), count = sum(!is.na(PDO))) %>%
    dplyr::filter(count == 5) %>%
    dplyr::select(Year, PDO_Summer)
  
  pdo_winter <- pdo_data %>%
    dplyr::filter(Month == 12 | Month <= 3) %>%
    dplyr::mutate(WinterYear = dplyr::if_else(Month == 12, Year, Year - 1)) %>%
    dplyr::group_by(WinterYear) %>%
    dplyr::summarise(PDO_Winter = mean(PDO, na.rm = TRUE), count = sum(!is.na(PDO))) %>%
    dplyr::filter(count == 4) %>%
    dplyr::select(Year = WinterYear, PDO_Winter)
  
  pdo_seasonal <- dplyr::full_join(pdo_summer, pdo_winter, by = "Year") %>% dplyr::arrange(Year)
  
  # 3. NPGO
  temp_npgo <- tempfile()
  curl::curl_download("https://www.o3d.org/npgo/data/NPGO.txt", temp_npgo)
  npgo_raw <- read.table(temp_npgo, skip = 10, header = FALSE,
                         col.names = c("Year", "Month", "NPGO")) %>%
    dplyr::mutate(
      Date = as.Date(paste(Year, Month, "01", sep = "-")),
      NPGO = as.numeric(NPGO),
      Year = lubridate::year(Date),
      Month = lubridate::month(Date)
    ) %>%
    dplyr::filter(!is.na(NPGO))
  
  npgo_summer <- npgo_raw %>%
    dplyr::filter(Month >= 5 & Month <= 9) %>%
    dplyr::group_by(Year) %>%
    dplyr::summarise(NPGO_Summer = mean(NPGO, na.rm = TRUE), count = sum(!is.na(NPGO))) %>%
    dplyr::filter(count == 5) %>%
    dplyr::select(Year, NPGO_Summer)
  
  npgo_winter <- npgo_raw %>%
    dplyr::filter(Month == 12 | Month <= 3) %>%
    dplyr::mutate(WinterYear = dplyr::if_else(Month == 12, Year, Year - 1)) %>%
    dplyr::group_by(WinterYear) %>%
    dplyr::summarise(NPGO_Winter = mean(NPGO, na.rm = TRUE), count = sum(!is.na(NPGO))) %>%
    dplyr::filter(count == 4) %>%
    dplyr::select(Year = WinterYear, NPGO_Winter)
  
  npgo_seasonal <- dplyr::full_join(npgo_summer, npgo_winter, by = "Year") %>% dplyr::arrange(Year)
  
  # 4. Sea level
  temp_sea <- tempfile()
  curl::curl_download("https://psmsl.org/data/obtaining/rlr.monthly.data/385.rlrdata", temp_sea)
  sea_level_data <- read.csv(temp_sea, header = FALSE, sep = ";",
                             col.names = c("Decimal_Year", "Sea_Level", "Flag1", "Flag2"))
  
  sea_level_data$Decimal_Year <- gsub(";", "", sea_level_data$Decimal_Year)
  sea_level_data$Decimal_Year <- as.numeric(sea_level_data$Decimal_Year)
  
  sea_level_clean <- sea_level_data %>%
    dplyr::filter(Sea_Level != -9999900000) %>%
    dplyr::mutate(
      Year = floor(Decimal_Year),
      Month = floor((Decimal_Year - Year) * 12) + 1
    )
  
  sea_summer <- sea_level_clean %>%
    dplyr::filter(Month >= 5 & Month <= 9) %>%
    dplyr::group_by(Year) %>%
    dplyr::summarise(Sea_Level_Summer = mean(Sea_Level, na.rm = TRUE))
  
  sea_winter <- sea_level_clean %>%
    dplyr::filter(Month == 12 | Month <= 3) %>%
    dplyr::mutate(WinterYear = dplyr::if_else(Month == 12, Year, Year - 1)) %>%
    dplyr::group_by(WinterYear) %>%
    dplyr::summarise(Sea_Level_Winter = mean(Sea_Level, na.rm = TRUE), count = sum(!is.na(Sea_Level))) %>%
    dplyr::filter(count == 4) %>%
    dplyr::select(Year = WinterYear, Sea_Level_Winter)
  
  sea_level_seasonal <- dplyr::full_join(sea_summer, sea_winter, by = "Year") %>% dplyr::arrange(Year)
  
  # 5. NOAA stoplight indicators
  indicators_raw <- read.csv(stoplight_path, header = FALSE)
  indicators_raw <- indicators_raw[-22, ]
  indicators_raw$V1[1] <- "year"
  indicators <- as.data.frame(t(indicators_raw))
  colnames(indicators) <- c(
    "year", "pdo_winter", "pdo_summer", "oni", "sst_summer", 
    "sst_u20_winter", "sst_u20_summer", "t_deep_summer", "s_deep_summer", 
    "cop_rich_summer", "cop_north_summer", "cop_south_summer", "bio_trans", 
    "icht_bio_winter", "icht_index_winter", "chin_juv_june", "coho_juv_june", 
    "mean_ranks", "rank_ranks", "pc1", "pc2",  "upw_phys_trans", 
    "upw_hydr_trans", "upw_anom_spring", "upw_length", "cop_summer_index"
  )
  indicators <- indicators[-1, ]
  indicators$year <- floor(as.numeric(indicators$year))
  indicators <- suppressWarnings(
    dplyr::mutate(indicators, dplyr::across(-year, as.numeric))
  )
  
  # 6. Join ocean & indicator data
  ocean_covariates <- dplyr::full_join(pdo_seasonal, npgo_seasonal, by = "Year") %>%
    dplyr::full_join(sea_level_seasonal, by = "Year") %>%
    dplyr::full_join(indicators, by = c("Year" = "year"))
  
  # 7. Brood Year shift
  ocean_covariates <- ocean_covariates %>%
    dplyr::mutate(Brood_Year = Year - 1) %>%
    dplyr::select(-Year)
  
  # 8. Merge all
  all_covariates <- dplyr::full_join(ocean_covariates, qmax_df, by = "Brood_Year") %>%
    dplyr::arrange(Brood_Year) %>%
    dplyr::rename_with(~ paste0("Cov_", .), -Brood_Year)
  
  return(all_covariates)
}

add_covs_per_stock <- function(data_list, covariates_df, covariate_map) {
  updated_list <- purrr::imap(data_list, function(df, name) {
    covariates_to_add <- covariate_map[[name]]
    
    if (is.null(covariates_to_add)) {
      warning(glue::glue("No covariates defined for '{name}', leaving as is."))
      return(df)
    }
    
    join_data <- covariates_df %>%
      dplyr::select(Brood_Year, all_of(covariates_to_add))
    
    dplyr::left_join(df, join_data, by = "Brood_Year")
  })
  
  return(updated_list)
}


# rescale.covars <- TRUE
# 
# if(rescale.covars){
#   
#   cov.idx <- which(grepl("Cov",dimnames(raw.data)[[2]]))
#   print(cov.idx)
#   
#   for(i in cov.idx){
#     raw.data[,i] <-   raw.data[,i] / max(abs(raw.data[,i] ),na.rm=TRUE)
#   }
# }
# 
# raw.data[,cov.idx]