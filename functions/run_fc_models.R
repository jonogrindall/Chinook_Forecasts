#==============================================================================
# CUSTOM WALK-FORWARD MODEL RUNNER - MAIN ORCHESTRATOR FUNCTION
#==============================================================================
#
# PURPOSE: This function orchestrates the execution of 7 custom walk-forward 
# forecasting models and integrates them with the forecastR framework.
#
# CUSTOM MODELS INCLUDED:
# 1. LinearCov_walkforward - Linear models with environmental covariates
# 2. LinearCov_AICc_walkforward - Linear models with AICc covariate selection
# 3. LinearCov_BoxCox_walkforward - Linear models with Box-Cox transformation
# 4. LinearCov_BoxCox_AICc_walkforward - Box-Cox transformation + AICc selection
# 5. SibReg_BoxCox_walkforward - Sibling regression with Box-Cox transformation
# 6. GLM_cov_walkforward - Generalized Linear Models with covariates
# 7. GAM_cov_walkforward - Generalized Additive Models with smooth terms
#
# WALK-FORWARD VALIDATION APPROACH:
# Each model implements proper retrospective validation using:
# - Training data: All years before validation year
# - Validation data: Single year being predicted
# - Rolling window: Process repeats for each year in retrospective period
# - Performance metrics: MRE, MAE, MPE, MAPE, MASE, RMSE calculated across all validation years
# - Final forecasts: Use all available data to predict current year (2025)
#
# INPUT PARAMETERS:
# - forecast_yr: Target forecast year (e.g., 2025)
# - raw_data_short: List of historical run data by stock (for custom models)
# - raw_data_forecastr: List of covariate data by stock (for forecastR models)
# - fc: Existing forecast object from forecastR (contains standard model results)
# - diagnostics_env: Environment to store model diagnostics and intermediate results
#
# OUTPUT:
# - Updated fc object with custom model results integrated
# - Each model adds performance metrics to fc[[stock]]$retro.pm$retro.pm.all.varyrs
# - Each model adds forecast predictions to fc[[stock]]$table.ptf
# - Custom models appear in ranking tables alongside forecastR models
#
# COVARIATE STRUCTURE:
# - sky_fingerlings/sky_yearlings: Cov_sky_qmax, Cov_pdo_summer, Cov_pc1
# - snq_fingerlings/snq_yearlings: Cov_snq_qmax, Cov_pdo_summer, Cov_pc1
# - tul_fingerlings: Cov_pc1 (Age 2-4), Cov_NPGO_Summer (Age 5)
#
# DEPENDENCIES:
# - All individual model functions in functions/ directory
# - forecastR package for integration
# - dplyr, tibble, purrr, abind for data manipulation
#==============================================================================

run_fc_models <- function(
  forecast_yr,
  raw_data_short,
  raw_data_forecastr,
  fc,
  diagnostics_env
) {
  
  #============================================================================
  # 1. INITIALIZATION AND VALIDATION
  #============================================================================
  
  cat("Initializing custom walk-forward model runner...\n")
  
  # Validate inputs
  if (is.null(fc) || length(fc) == 0) {
    stop("Error: fc object is NULL or empty. Run forecastR models first.")
  }
  
  if (is.null(raw_data_short) || length(raw_data_short) == 0) {
    stop("Error: raw_data_short is NULL or empty.")
  }
  
  if (is.null(raw_data_forecastr) || length(raw_data_forecastr) == 0) {
    stop("Error: raw_data_forecastr is NULL or empty.")
  }
  
  # Get list of stocks to process
  stocks <- names(fc)
  cat("Processing stocks:", paste(stocks, collapse = ", "), "\n")
  
  #============================================================================
  # 2. DEFINE MODEL CONFIGURATION
  #============================================================================
  
  # Define the 7 custom models to run
  methods <- c(
    "lm_cov_walkforward",           # Linear models with covariates
    "lm_cov_aicc_walkforward",      # Linear models with AICc selection
    "lm_cov_bxcx_walkforward",      # Linear models with Box-Cox transformation
    "lm_cov_bxcx_aicc_walkforward", # Linear models with Box-Cox + AICc
    "sibreg_bxcx_walkforward",      # Sibling regression with Box-Cox
    "glm_cov",                      # Generalized Linear Models
    "gam_cov"                       # Generalized Additive Models
  )
  
  cat("Running", length(methods), "custom walk-forward models:\n")
  cat(paste("  ", methods, collapse = "\n"), "\n\n")
  
  #============================================================================
  # 3. MODEL EXECUTION LOOP
  #============================================================================
  
  # Define covariates by stock and age (based on original analysis)
  covars_by_stock_age <- list(
    "sky_fingerlings" = list(
      "2" = c("Cov_pdo_summer"),
      "3" = c("Cov_sky_qmax", "Cov_pc1"),
      "4" = c("Cov_NPGO_Summer", "Cov_sky_qmax"),
      "5" = c("Cov_NPGO_Summer", "Cov_sky_qmax")
    ),
    "sky_yearlings" = list(
      "3" = c("Cov_pc1"),
      "4" = c("Cov_sky_qmax"),
      "5" = c("Cov_NPGO_Summer")
    ),
    "tul_fingerlings" = list(
      "2" = c("Cov_pc1"),
      "3" = c("Cov_pc1"),
      "4" = c("Cov_pc1"),
      "5" = c("Cov_NPGO_Summer")
    ),
    "snq_fingerlings" = list(
      "2" = c("Cov_NPGO_Summer"),
      "3" = c("Cov_snq_qmax", "Cov_pc1"),
      "4" = c("Cov_pc1"),
      "5" = c("Cov_pdo_summer")
    ), 
    "snq_yearlings" = list(
      "3" = c("Cov_snq_qmax"),
      "4" = c("Cov_NPGO_Summer"),
      "5" = c("Cov_NPGO_Summer","Cov_snq_qmax")
    )
  )

  # Process each stock
  for (stock in stocks) {
    cat("Processing stock:", stock, "\n")
    
    # Get covariate information for this stock
    covars <- covars_by_stock_age[[stock]]
    
    # Run each custom model
    for (method in methods) {
      cat("  Running method:", method, "\n")
      
      tryCatch({
        #========================================================================
        # 4. DYNAMIC FUNCTION SOURCING AND EXECUTION
        #========================================================================
        
        # Determine which source file to use for each method
        source_file <- switch(
          method,
          lm_cov_walkforward = "functions/lm_cov_walkforward.R",
          lm_cov_aicc_walkforward = "functions/lm_cov_aiccc_walkforward.R",
          lm_cov_bxcx_walkforward = "functions/lm_cov_bxcx_walkforward.R",
          lm_cov_bxcx_aicc_walkforward = "functions/lm_cov_bxcx_aicc_walkforward_fixed.R",
          sibreg_bxcx_walkforward = "functions/sibreg_bxcx_walkforward.R",
          glm_cov = "functions/glm_cov_walkforward.R",
          gam_cov = "functions/gam_cov_walkforward.R",
          stop("Unknown method: ", method)
        )
        
        # Source the appropriate function file
        source(source_file)
        
        # Execute the model function with appropriate parameters
        fc <- switch(
          method,
          lm_cov_walkforward = update_forecast_w_lm_cov(
            stock, raw_data_short, fc, forecast_yr, covars, diagnostics_env
          ),
          lm_cov_aicc_walkforward = update_fc_w_lm_cov_AICc(
            stock, raw_data_short, fc, forecast_yr, covars, diagnostics_env
          ),
          lm_cov_bxcx_walkforward = update_fc_w_lm_cov_bxcx(
            stock, raw_data_short, fc, forecast_yr, covars, diagnostics_env
          ),
          lm_cov_bxcx_aicc_walkforward = update_fc_w_lm_cov_bxcx_aicc(
            stock, raw_data_short, fc, forecast_yr, covars, diagnostics_env
          ),
          sibreg_bxcx_walkforward = update_fc_w_lm_cov_bxcx_sibreg(
            stock, raw_data_short, fc, forecast_yr, diagnostics_env
          ),
          glm_cov = update_fc_w_glm_cov(
            stock, raw_data_short, fc, forecast_yr, covars, diagnostics_env
          ),
          gam_cov = update_fc_w_gam_cov(
            stock, raw_data_short, fc, forecast_yr, covars, diagnostics_env
          )
        )
        
        cat("    ✓ Completed successfully\n")
        
      }, error = function(e) {
        # Error handling for individual model failures
        cat("    ✗ Error in", method, "for", stock, ":", e$message, "\n")
        cat("    Continuing with next model...\n")
      })
    }
    
    cat("  Completed all models for", stock, "\n\n")
  }
  
  #============================================================================
  # 5. VALIDATION AND SUMMARY
  #============================================================================
  
  cat("Custom walk-forward model execution completed!\n")
  
  # Validate that fc object is still intact
  if (is.null(fc)) {
    stop("Critical error: fc object became NULL during execution")
  }
  
  # Count models added to each stock
  for (stock in stocks) {
    if (!is.null(fc[[stock]]$table.ptf)) {
      model_count <- nrow(fc[[stock]]$table.ptf)
      cat("Stock", stock, "has", model_count, "models in table.ptf\n")
    }
  }
  
  cat("\nAll custom models processed successfully!\n")
  
  # Return the updated fc object
  return(fc)
}
