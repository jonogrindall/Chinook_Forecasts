#==============================================================================
# WALK-FORWARD LINEAR MODEL WITH COVARIATES
#==============================================================================
#
# This function implements a linear regression model with environmental covariates
# using proper walk-forward retrospective validation. The model predicts salmon
# returns based on historical run data and environmental covariates.
#
# MODEL DESCRIPTION:
# - Response variable: Average_Terminal_Run (salmon returns)
# - Predictors: Environmental covariates (e.g., ocean conditions, climate indices)
# - Model type: Linear regression (lm)
# - Validation: Walk-forward retrospective analysis
#
# WALK-FORWARD VALIDATION PROCESS:
# 1. For each validation year in the retrospective period:
#    - Training data: All years before validation year
#    - Validation data: Single year being predicted
#    - Fit linear model on training data
#    - Predict validation year using fitted model
#    - Calculate prediction error
# 2. Aggregate errors across all validation years
# 3. Calculate performance metrics (MRE, MAE, MPE, MAPE, MASE, RMSE)
# 4. Generate final forecast for current year using all available data
#
# COVARIATE HANDLING:
# - Covariates are selected based on age class using get_covars_by_age()
# - Only covariates available in training period are used
# - Missing covariate values are handled gracefully
#
# INPUT PARAMETERS:
# - stock_name: Name of the salmon stock (e.g., "tul_fingerlings")
# - raw_data_short: List containing historical run data by stock
# - fc: Forecast object from forecastR (will be updated)
# - forecast_yr: Year to forecast (e.g., 2025)
# - covars_by_age: List of covariates by age class
# - diagnostics_env: Environment to store model diagnostics
# - min_retro_yrs: Minimum years required for retrospective analysis (default: 14)
#
# OUTPUT:
# - Updated fc object with LinearCov model results added
# - Performance metrics added to fc[[stock]]$retro.pm$retro.pm.all.varyrs
# - Forecast predictions added to fc[[stock]]$table.ptf
#
# DEPENDENCIES:
# - walk_forward_template.R: Template function for walk-forward validation
# - forecastR package: For data structures and utilities
#==============================================================================

update_forecast_w_lm_cov <- function(
    stock_name,
    raw_data_short,
    fc,
    forecast_yr,
    covars_by_age,
    diagnostics_env,
    min_retro_yrs = 10
) {
  
  #============================================================================
  # MODEL FITTING FUNCTION
  #============================================================================
  # This function defines how to fit the linear model given training data
  # and available covariates. It's passed to the walk-forward template.
  
  fit_lm_model <- function(train_data, valid_covars) {
    # Fit standard linear regression model
    # Formula: Average_Terminal_Run ~ covariate1 + covariate2 + ...
    lm(Average_Terminal_Run ~ ., data = train_data)
  }
  
  #============================================================================
  # EXECUTE WALK-FORWARD VALIDATION
  #============================================================================
  # Use the walk-forward template to handle the retrospective validation
  # and forecast generation process
  
  source("functions/walk_forward_template.R")
  
  return(walk_forward_template(
    stock_name = stock_name,
    raw_data_short = raw_data_short,
    fc = fc,
    forecast_yr = forecast_yr,
    covars_by_age = covars_by_age,
    diagnostics_env = diagnostics_env,
    model_name = "LinearCov",
    model_fitting_function = fit_lm_model,
    min_retro_yrs = min_retro_yrs
  ))
}
