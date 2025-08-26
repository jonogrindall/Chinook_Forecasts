# WALK-FORWARD LINEAR COVARIATE MODEL WITH AICc SELECTION
# Replaces lm_cov_aiccc.R with proper retrospective validation

update_fc_w_lm_cov_AICc <- function(
    stock_name,
    raw_data_short,
    fc,
    forecast_yr,
    covars_by_age,
    diagnostics_env,
    min_retro_yrs = 10
) {
  library(MuMIn)
  options(na.action = "na.fail")
  
  # Model fitting function for AICc selection
  fit_aicc_model <- function(train_data, valid_covars) {
    # Fit full model
    full_model <- lm(Average_Terminal_Run ~ ., data = train_data)
    
    # Use MuMIn for model selection based on AICc
    tryCatch({
      # Perform automated model selection
      candidate_models <- dredge(full_model, rank = "AICc")
      best_model <- get.models(candidate_models, 1)[[1]]
      return(best_model)
    }, error = function(e) {
      # Fallback to full model if dredge fails
      cat("      AICc selection failed, using full model:", e$message, "\n")
      return(full_model)
    })
  }
  
  # Use the walk-forward template
  source("functions/walk_forward_template.R")
  
  return(walk_forward_template(
    stock_name = stock_name,
    raw_data_short = raw_data_short,
    fc = fc,
    forecast_yr = forecast_yr,
    covars_by_age = covars_by_age,
    diagnostics_env = diagnostics_env,
    model_name = "LinearCov_AICc",
    model_fitting_function = fit_aicc_model,
    min_retro_yrs = min_retro_yrs
  ))
}
