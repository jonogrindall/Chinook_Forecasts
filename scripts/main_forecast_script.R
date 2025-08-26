#==============================================================================
# SALMON FORECASTING WORKFLOW - MAIN SCRIPT
#==============================================================================
# 
# PURPOSE: This script runs a complete salmon forecasting analysis using the 
# forecastR package plus 7 custom walk-forward models with comprehensive
# retrospective validation and model ranking.
# 
# WORKFLOW OVERVIEW:
# 1. Data loading and preprocessing from individual stock CSV files
# 2. Running 5 standard forecastR models (Naive, Sibling Regression variants)
# 3. Running 7 custom walk-forward models with environmental covariates:
#    - Linear Models with covariates (3 variants: basic, AICc, Box-Cox)
#    - Box-Cox + AICc combined model
#    - Sibling Regression with Box-Cox transformation
#    - Generalized Linear Model (GLM) with covariates
#    - Generalized Additive Model (GAM) with smooth terms
# 4. Model ranking and performance evaluation across all age classes
# 5. Generation of forecast predictions for current year (2025)
#
# CUSTOM MODELS INCLUDED:
# - LinearCov: Standard linear regression with environmental covariates
# - LinearCov_AICc: Linear regression with AICc covariate selection
# - LinearCov_BoxCox: Linear regression with Box-Cox transformation
# - LinearCov_BoxCox_AICc: Box-Cox transformation + AICc selection
# - SibReg_BoxCox: Sibling regression with Box-Cox transformation
# - GLM: Generalized Linear Models with covariates
# - GAM: Generalized Additive Models with smooth terms
#
# WALK-FORWARD VALIDATION:
# - All custom models implement proper retrospective validation
# - Minimum 10 years of retrospective data required
# - Out-of-sample predictions for each validation year
# - Performance metrics: MRE, MAE, MPE, MAPE, MASE, RMSE
#
# CONFIGURATION:
# - forecast_yr: 2025 (current forecast year, updates annually)
# - min_retro_yrs: 10 (minimum years for retrospective analysis)
# - retro.min.yrs: 10 (minimum years for forecastR models)
#
# DEPENDENCIES:
# - forecastR package (from GitHub: SalmonForecastR/ForecastR-Package)
# - tidyverse, MASS, mgcv, abind packages
# - Custom functions in functions/ directory
#
# INPUT FILES:
# - Individual stock CSV files in data/ directory
# - Each CSV contains: Stock, Brood_Year, Age_Class, Average_Terminal_Run, covariates
#
# OUTPUT:
# - fc object with all model results and rankings
# - table.ptf with forecast predictions for 2025
# - ranked_forecasts with model performance rankings
# - retro.pm with retrospective performance metrics
#==============================================================================

# Clear workspace and load libraries
rm(list = ls())
gc()

# Load required libraries
library(tidyverse)    # Data manipulation and visualization
library(MASS)         # For Box-Cox transformations
library(mgcv)         # For GAM models
library(abind)        # For array operations

# Load forecastR package (make sure it's available)
if (!require(forecastR)) {
  stop("forecastR package is not installed. Please install it first.")
}

#==============================================================================
# 1. DATA LOADING AND PREPROCESSING
#==============================================================================

cat("Loading and preprocessing data...\n")

# Source the data preparation script to create the lists
source("scripts/prepare_data.R")

cat("Data loaded successfully!\n")
cat("Stocks found:", names(raw_data_short), "\n")
cat("Short data stocks:", length(raw_data_short), "\n")
cat("ForecastR data stocks:", length(raw_data_forecastr), "\n\n")

#==============================================================================
# 2. FORECASTING PARAMETERS
#==============================================================================

# Set forecast year (will change each year)
forecast_yr <- 2025

# Define model settings for forecastR
model.settings <- list(
  Naive3 = list(model.type = "Naive", settings = list(avg.yrs = 3)),
  Naive5 = list(model.type = "Naive", settings = list(avg.yrs = 5)),
  SibRegSimple = list(model.type = "SibRegSimple", settings = NULL),
  SibRegLogPower = list(model.type = "SibRegLogPower", settings = NULL),
  SibRegKalman = list(model.type = "SibRegKalman", settings = NULL)
)

cat("Forecast year:", forecast_yr, "\n")
cat("Models to run:", paste(names(model.settings), collapse = ", "), "\n\n")

#==============================================================================
# 3. RUN STANDARD FORECASTR MODELS
#==============================================================================

cat("Running standard forecastR models...\n")

# Run multiFC to get baseline models (Naive and SibReg variants)
fc <- lapply(names(raw_data_forecastr), function(name) {
  cat("Running multiFC on:", name, "\n")
  
  # Convert tibble to data.frame
  df <- as.data.frame(raw_data_forecastr[[name]])
  
  # Run multiFC
  fc <- multiFC(
    data.file = df,
    settings.list = model.settings,
    retro.min.yrs = 10,
    out.type = "full",
    int.type = "Retrospective",
    do.retro = TRUE,
    tracing = TRUE
  )
})

# Name the list with original tibble names
names(fc) <- names(raw_data_forecastr)

cat("Standard forecastR models completed!\n")
cat("Stocks processed:", names(fc), "\n\n")

#==============================================================================
# 4. RUN CUSTOM WALK-FORWARD MODELS
#==============================================================================

cat("Running custom walk-forward models...\n")

# Create diagnostics environment for storing model diagnostics
diagnostics_env <- new.env()

# Source custom model functions
source("functions/run_fc_models.R")

# Source custom ranking function
source("functions/get_ranked_fc.R")

# Run custom walk-forward models
fc <- run_fc_models(
  forecast_yr = forecast_yr,
  raw_data_short = raw_data_short,
  raw_data_forecastr = raw_data_forecastr,
  fc = fc,
  diagnostics_env = diagnostics_env
)

cat("Custom walk-forward models completed!\n\n")

#==============================================================================
# 5. MODEL RANKING AND EVALUATION
#==============================================================================

cat("Generating model rankings...\n")

# Get ranked forecasts using forecastR's ranking function
if (exists("get_ranked_fc")) {
  ranked_forecasts <- get_ranked_fc(fc)
} else {
  cat("Warning: get_ranked_fc function not found. Skipping ranking generation.\n")
  ranked_forecasts <- NULL
}

cat("Model rankings generated!\n\n")

#==============================================================================
# 6. RESULTS SUMMARY
#==============================================================================

cat("=== FORECASTING RESULTS SUMMARY ===\n\n")

# Show available models for each stock
for (stock_name in names(fc)) {
  cat("Stock:", stock_name, "\n")
  
  # Show models in table.ptf (forecast predictions)
  if (!is.null(fc[[stock_name]]$table.ptf)) {
    models_in_ptf <- rownames(fc[[stock_name]]$table.ptf)
    cat("  Models with 2025 forecasts:", length(models_in_ptf), "\n")
    cat("  Model names:", paste(models_in_ptf, collapse = ", "), "\n")
  }
  
  # Show models in ranksum (performance rankings)
  if (!is.null(ranked_forecasts$ranked_forecasts[[stock_name]]$ranksum)) {
    models_in_ranksum <- ranked_forecasts$ranked_forecasts[[stock_name]]$ranksum$Model
    cat("  Models in rankings:", length(models_in_ranksum), "\n")
    cat("  Top 3 models:", paste(head(models_in_ranksum, 3), collapse = ", "), "\n")
  }
  
  cat("\n")
}

#==============================================================================
# 7. SAVE RESULTS
#==============================================================================

cat("Saving results...\n")

# Save main results
saveRDS(fc, "output/forecast_results.rds")
saveRDS(ranked_forecasts, "output/ranked_forecasts.rds")

# Save diagnostics if any were generated
if (length(ls(diagnostics_env)) > 0) {
  saveRDS(as.list(diagnostics_env), "output/model_diagnostics.rds")
}

cat("Results saved to output/ directory\n\n")

#==============================================================================
# 8. QUICK ACCESS TO KEY RESULTS
#==============================================================================

cat("=== QUICK ACCESS COMMANDS ===\n")
cat("# View all models for a specific stock:\n")
cat("stock <- 'tul_fingerlings'  # Change stock name as needed\n")
cat("print(rownames(fc[[stock]]$table.ptf))\n\n")

cat("# View model rankings for a stock:\n")
cat("print(ranked_forecasts$ranked_forecasts[[stock]]$ranksum)\n\n")

cat("# View 2025 forecasts for a stock:\n")
cat("print(fc[[stock]]$table.ptf)\n\n")

cat("# View best models by age class:\n")
cat("print(ranked_forecasts$ranked_forecasts[[stock]]$bestmodel)\n\n")

cat("=== SCRIPT COMPLETED SUCCESSFULLY ===\n")
