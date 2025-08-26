#==============================================================================
# DATA PREPARATION SCRIPT - SALMON FORECASTING
#==============================================================================
#
# PURPOSE: This script prepares the data for the salmon forecasting analysis by
# loading individual stock CSV files and creating the necessary data structures
# for both custom walk-forward models and standard forecastR models.
#
# DATA STRUCTURE:
# 1. raw_data_short: List of dataframes for custom walk-forward models
#    - Contains historical run data with proper data types
#    - One dataframe per stock (sky_fingerlings, snq_fingerlings, etc.)
#    - Used by custom models that implement walk-forward validation
#
# 2. raw_data_forecastr: List of dataframes for standard forecastR models
#    - Contains both run data and environmental covariates
#    - One dataframe per stock with all required variables
#    - Used by forecastR package models (Naive, Sibling Regression)
#
# INPUT REQUIREMENTS:
# - Individual stock CSV files in data/ directory
# - Each CSV must contain: Stock, Brood_Year, Age_Class, Average_Terminal_Run
# - Environmental covariates specific to each stock
# - Data should be clean with no missing values in key fields
#
# DATA PROCESSING:
# - Ensures proper data types (integer for Brood_Year, Age_Class)
# - Handles missing values appropriately
# - Creates consistent structure across all stocks
# - Validates data integrity before modeling
#
# OUTPUT:
# - raw_data_short: List of dataframes for custom models
# - raw_data_forecastr: List of dataframes for forecastR models
# - Both lists have same stock names as keys
#
# DEPENDENCIES:
# - dplyr for data manipulation
# - readr for CSV file reading
# - Individual stock CSV files in data/ directory
#==============================================================================

library(tidyverse)

cat("Preparing data files for salmon forecasting...\n")

#==============================================================================
# 1. CREATE RAW_DATA_SHORT LIST
#==============================================================================

cat("Creating raw_data_short list...\n")

# List all short files
short_files <- list.files("data/", pattern = "*_short.csv", full.names = TRUE)
cat("Found", length(short_files), "short data files\n")

# Read each short file into a list
raw_data_short <- list()
for (file in short_files) {
  stock_name <- gsub(".*/(.+)_short.csv", "\\1", file)
  cat("  Reading", stock_name, "\n")
  
  data <- read.csv(file) %>%
    as_tibble() %>%
    mutate(
      Brood_Year = as.integer(Brood_Year),
      Age_Class = as.integer(Age_Class),
      Average_Terminal_Run = as.numeric(Average_Terminal_Run)
    )
  
  raw_data_short[[stock_name]] <- data
}

cat("Created raw_data_short list with", length(raw_data_short), "stocks\n")
cat("  Stocks:", paste(names(raw_data_short), collapse = ", "), "\n\n")

#==============================================================================
# 2. CREATE RAW_DATA_FORECASTR LIST
#==============================================================================

cat("Creating raw_data_forecastr list...\n")

# List all forecastr files
forecastr_files <- list.files("data/", pattern = "*_forecastr.csv", full.names = TRUE)
cat("Found", length(forecastr_files), "forecastr data files\n")

# Read each forecastr file into a list
raw_data_forecastr <- list()
for (file in forecastr_files) {
  stock_name <- gsub(".*/(.+)_forecastr.csv", "\\1", file)
  cat("  Reading", stock_name, "\n")
  
  data <- read.csv(file) %>%
    as_tibble() %>%
    mutate(
      Brood_Year = as.integer(Brood_Year),
      Age_Class = as.integer(Age_Class)
    )
  
  raw_data_forecastr[[stock_name]] <- data
}

cat("Created raw_data_forecastr list with", length(raw_data_forecastr), "stocks\n")
cat("  Stocks:", paste(names(raw_data_forecastr), collapse = ", "), "\n\n")

#==============================================================================
# 3. DATA SUMMARY
#==============================================================================

cat("=== DATA SUMMARY ===\n")
cat("raw_data_short:\n")
for (stock in names(raw_data_short)) {
  cat("  ", stock, ": ", nrow(raw_data_short[[stock]]), " rows, ", 
      paste(range(raw_data_short[[stock]]$Brood_Year), collapse = "-"), " years\n")
}

cat("\nraw_data_forecastr:\n")
for (stock in names(raw_data_forecastr)) {
  cat("  ", stock, ": ", nrow(raw_data_forecastr[[stock]]), " rows, ", 
      paste(range(raw_data_forecastr[[stock]]$Brood_Year), collapse = "-"), " years\n")
}

cat("\nData preparation completed successfully!\n")
cat("Two lists created: raw_data_short and raw_data_forecastr\n")
cat("You can now run the main forecasting script.\n")
