term_run_df <- function(
    forecast_yr,
    base_dir = "R:/Fisheries Management/forecasting/2025/code_dev/all_forecasts"
) {
	# --- Load libraries ---
	library(dplyr)
	library(tidyr)
	library(readr)
	library(tibble)
	library(stringr)
  library(DBI)
  library(odbc)
  


  # --- Helper: add metadata columns up front (1st row filled, remainder blank) ---
  add_metadata <- function(df, stock_name) {
    n_rows <- nrow(df)
    meta_row <- tibble(
      Stock_Name       = stock_name,
      Stock_Species    = "Chinook",
      Stock_Abundance  = "Terminal_Run",
      Forecasting_Year = forecast_yr
    )
    meta_blanks <- tibble(
      Stock_Name       = rep("", n_rows - 1),
      Stock_Species    = rep("", n_rows - 1),
      Stock_Abundance  = rep("", n_rows - 1),
      Forecasting_Year = rep(NA_real_, n_rows - 1)
    )
    meta_df <- bind_rows(meta_row, meta_blanks)
    bind_cols(meta_df, df)
  }
  
  # --- Load base escapement (SKY + SNQ) ---
  esc_base <- read_csv(
    file.path(base_dir, "data/escapement.csv"),
    show_col_types = FALSE
  ) %>%
    dplyr::select(year, sky_nor_esc, snq_nor_esc) %>%
    filter(between(year, 2006, forecast_yr - 1))
  
  # --- Load Tulalip escapement & merge ---
  
  tul_escapement <- read_csv(
    file.path(base_dir, "data/tul_escapement.csv"),
    show_col_types = FALSE
  ) %>%
    filter(between(year, 2006, forecast_yr - 1)) %>%
    dplyr::select(year, tul_esc)
  
  
  escapement <- esc_base %>%
    left_join(tul_escapement, by = "year")
  
  # --- Load age structures ---
  # sky_age_structure <- read_csv(file.path(base_dir, "data/sky_age_structure.csv"), show_col_types = FALSE)
  # snq_age_structure <- read_csv(file.path(base_dir, "data/snq_age_structure.csv"), show_col_types = FALSE)
  # con <- dbConnect(odbc(),
  #                  Driver = "SQL Server",
  #                  Server = "SVR-WEBSQLCLUSTER",
  #                  Database = "nrfss",
  #                  Trusted_Connection = "Yes")
  # 
  # source("functions/age_structure.R")
  age_structure <- sky_snq_age_str_df(con)
  sky_age_structure <- age_structure$sky_summary
  sky_age_structure[sky_age_structure == 0] <- 1
  snq_age_structure <- age_structure$snq_summary
  snq_age_structure[snq_age_structure == 0] <- 1
  
  # dbDisconnect(con)
  
  # Tulalip: ReturnYear, Total Fish, Age2.Fingerling...Age5.Fingerling
  tul_age_structure <- read_csv(file.path(base_dir, "data/tul_age_structure.csv"), show_col_types = FALSE) %>%
    rename( `TotalSamples` = `Total Fish`)
  
  
  # keep ReturnYear as-is
  
  # --- Load exploitation rate data (both Unmarked & Marked) ---
  er_raw <- read_csv(file.path(base_dir, "data/brood_year_er.csv"), show_col_types = FALSE)
  
  # Generic ER extender: pick column by name
  extend_er <- function(er_df, er_col, forecast_year) {
    er_vals <- er_df[[er_col]]
    last_two_avg <- mean(tail(er_vals, 2), na.rm = TRUE)
    years <- min(er_df$BroodYear):(forecast_year - 1)
    tibble(
      Brood_Year = years,
      ER = ifelse(
        years %in% er_df$BroodYear,
        er_vals[match(years, er_df$BroodYear)],
        last_two_avg
      )
    )
  }
  
  # Build ER series
  full_er_unmarked <- er_raw %>%
    dplyr::select(BroodYear, SKY.Indicator.Unmarked.ER) %>%
    extend_er("SKY.Indicator.Unmarked.ER", forecast_yr)
  
  full_er_marked <- er_raw %>%
    dplyr::select(BroodYear, SKY.Indicator.Marked.ER) %>%
    extend_er("SKY.Indicator.Marked.ER", forecast_yr)
  
  # --- Abundance calc (Fingerling + Yearling; used by SKY & SNQ) ---
  abundance_df <- function(esc_df, age_df, prefix, er_df, forecast_yr) {
    total_col <- paste0(prefix, "_nor_esc")
    
    age_cols <- grep("^Age", names(age_df), value = TRUE)
    finger_ages <- sort(readr::parse_number(grep("Fingerling", age_cols, value = TRUE)))
    yearl_ages  <- sort(readr::parse_number(grep("Yearling",   age_cols, value = TRUE)))
    
    df <- dplyr::left_join(
      esc_df %>% dplyr::select(year, dplyr::starts_with(prefix)),
      age_df %>% dplyr::rename(year = ReturnYear),
      by = "year"
    ) %>%
      tidyr::pivot_longer(
        cols = dplyr::all_of(age_cols),
        names_to = "Age_Rearing",
        values_to = "AgeValue"
      ) %>%
      dplyr::mutate(
        Age_Class = readr::parse_number(Age_Rearing),
        Rearing   = ifelse(grepl("Fingerling", Age_Rearing), "Fingerling", "Yearling"),
        Prop      = AgeValue / `TotalSamples`,
        Average_Escapement = NA_real_,
        Brood_Year = dplyr::if_else(
          Rearing == "Fingerling",
          year - Age_Class,
          year - (Age_Class - 1)
        )
      ) %>%
      dplyr::rename(Run_Year = year) %>%
      dplyr::left_join(er_df, by = c("Brood_Year" = "Brood_Year")) %>%
      dplyr::mutate(
        Average_Terminal_Run = dplyr::if_else(
          Run_Year == forecast_yr,
          NA_real_,
          (Prop * .data[[total_col]]) / (1 - ER)
        )
      ) %>%
      dplyr::select(Run_Year, Brood_Year, Age_Class, Rearing, Average_Escapement, Average_Terminal_Run)
    
    make_clean <- function(df, age_set, rearing_type) {
      if (length(age_set) == 0) {
        return(dplyr::tibble(
          Run_Year = numeric(0),
          Brood_Year = numeric(0),
          Age_Class = numeric(0),
          Average_Escapement = numeric(0),
          Average_Terminal_Run = numeric(0)
        ))
      }
      
      valid_years <- df %>%
        dplyr::filter(!is.na(Age_Class)) %>%
        dplyr::count(Brood_Year) %>%
        dplyr::filter(n == length(age_set)) %>%
        dplyr::pull(Brood_Year)
      
      if (length(valid_years) == 0) valid_years <- min(df$Brood_Year, na.rm = TRUE)
      
      df <- df %>%
        dplyr::filter(Brood_Year >= min(valid_years)) %>%
        dplyr::select(-Rearing) %>%
        dplyr::arrange(Brood_Year, Age_Class)
      
      pad <- tibble::tibble(
        Run_Year = rep(forecast_yr, length(age_set)),
        Age_Class = age_set
      ) %>%
        dplyr::mutate(
          Brood_Year = Run_Year - Age_Class + ifelse(rearing_type == "Fingerling", 0, 1),
          Average_Escapement = NA_real_,
          Average_Terminal_Run = NA_real_
        )
      
      dplyr::bind_rows(df, pad)
    }
    
    list(
      fingerlings = df %>% dplyr::filter(Rearing == "Fingerling") %>% make_clean(finger_ages, "Fingerling"),
      yearlings   = df %>% dplyr::filter(Rearing == "Yearling")   %>% make_clean(yearl_ages, "Yearling")
    )
  }
  
  # --- Abundance calc (Fingerlings ONLY; used by Tulalip) ---
  abundance_df_fingerling <- function(esc_df, age_df, prefix, er_df, forecast_yr) {
    total_col <- if (prefix == "tul") "tul_esc" else paste0(prefix, "_nor_esc")
    age_cols <- grep("^Age", names(age_df), value = TRUE)
    finger_ages <- sort(readr::parse_number(age_cols))  # all are Fingerling
    
    df <- dplyr::left_join(
      esc_df %>% dplyr::select(year, dplyr::starts_with(prefix)),
      age_df %>% dplyr::rename(year = ReturnYear),
      by = "year"
    ) %>%
      tidyr::pivot_longer(
        cols = dplyr::all_of(age_cols),
        names_to = "Age_Rearing",
        values_to = "AgeValue"
      ) %>%
      dplyr::mutate(
        Age_Class = readr::parse_number(Age_Rearing),
        Prop      = AgeValue / `TotalSamples`,
        Average_Escapement = NA_real_,
        Brood_Year = year - Age_Class,  # all Fingerlings
        Rearing = "Fingerling"
      ) %>%
      dplyr::rename(Run_Year = year) %>%
      dplyr::left_join(er_df, by = c("Brood_Year" = "Brood_Year")) %>%
      dplyr::mutate(
        Average_Terminal_Run = dplyr::if_else(
          Run_Year == forecast_yr,
          NA_real_,
          (Prop * .data[[total_col]]) / (1 - ER)
        )
      ) %>%
      dplyr::select(Run_Year, Brood_Year, Age_Class, Average_Escapement, Average_Terminal_Run)
    
    # trim to first brood year with full age set
    valid_years <- df %>%
      dplyr::filter(!is.na(Age_Class)) %>%
      dplyr::count(Brood_Year) %>%
      dplyr::filter(n == length(finger_ages)) %>%
      dplyr::pull(Brood_Year)
    
    if (length(valid_years) == 0) valid_years <- min(df$Brood_Year, na.rm = TRUE)
    
    df <- df %>%
      dplyr::filter(Brood_Year >= min(valid_years)) %>%
      dplyr::arrange(Brood_Year, Age_Class)
    
    # pad forecast year
    pad <- tibble(
      Run_Year = rep(forecast_yr, length(finger_ages)),
      Age_Class = finger_ages,
      Brood_Year = forecast_yr - finger_ages,  # all Fingerlings
      Average_Escapement = NA_real_,
      Average_Terminal_Run = NA_real_
    )
    
    dplyr::bind_rows(df, pad)
  }
  
  # --- Create abundance datasets ---
  sky_abundance <- abundance_df(escapement, sky_age_structure, "sky", full_er_unmarked, forecast_yr)
  snq_abundance <- abundance_df(escapement, snq_age_structure, "snq", full_er_unmarked, forecast_yr)
  tul_fingerlings <- abundance_df_fingerling(escapement, tul_age_structure, "tul", full_er_marked, forecast_yr)
  
  # --- Add metadata ---
  sky_fingerlings <- add_metadata(sky_abundance$fingerlings, "Sky_fing")
  sky_yearlings   <- add_metadata(sky_abundance$yearlings,   "Sky_year")
  snq_fingerlings <- add_metadata(snq_abundance$fingerlings, "Snq_fing")
  snq_yearlings   <- add_metadata(snq_abundance$yearlings,   "Snq_year")
  tul_fingerlings <- add_metadata(tul_fingerlings,           "Tul_fing")
  
  # --- Write CSVs (overwrite) ---
  output_dir <- file.path(base_dir, "output")
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  write_csv(sky_fingerlings, file.path(output_dir, "sky_fingerlings.csv"), na = "NA")
  write_csv(sky_yearlings,   file.path(output_dir, "sky_yearlings.csv"),   na = "NA")
  write_csv(snq_fingerlings, file.path(output_dir, "snq_fingerlings.csv"), na = "NA")
  write_csv(snq_yearlings,   file.path(output_dir, "snq_yearlings.csv"),   na = "NA")
  write_csv(tul_fingerlings, file.path(output_dir, "tul_fingerlings.csv"), na = "NA")
  
  message("CSV files written to: ", output_dir)
  
  return(list(
    sky_fingerlings = sky_fingerlings,
    sky_yearlings   = sky_yearlings,
    snq_fingerlings = snq_fingerlings,
    snq_yearlings   = snq_yearlings,
    tul_fingerlings = tul_fingerlings
  ))
}

# --- Example Usage ---
# res <- forecastr_df(
#   forecast_yr = 2025
# )