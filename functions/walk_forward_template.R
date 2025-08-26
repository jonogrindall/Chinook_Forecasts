#==============================================================================
# WALK-FORWARD VALIDATION TEMPLATE
#==============================================================================
#
# This is the core template function that implements proper retrospective validation
# for any statistical model. It provides a standardized framework for walk-forward
# analysis that can be used by all custom forecasting models.
#
# WALK-FORWARD VALIDATION METHODOLOGY:
# Walk-forward validation is a time series cross-validation technique that:
# 1. Uses a rolling window approach where the model is trained on historical data
# 2. Validates on a single future time point
# 3. Moves the window forward and repeats the process
# 4. Aggregates performance metrics across all validation periods
# 5. Generates final forecasts using all available data
#
# KEY FEATURES:
# - Proper temporal separation between training and validation data
# - Handles missing data and covariates gracefully
# - Calculates comprehensive error metrics (MRE, MAE, MPE, MAPE, MASE, RMSE)
# - Supports both per-age and total (intersectional) metrics
# - Generates forecasts for current year using latest available covariate data
# - Integrates seamlessly with forecastR package structure
#
# INPUT PARAMETERS:
# - stock_name: Name of the salmon stock
# - raw_data_short: Historical run data by stock
# - fc: Forecast object from forecastR (will be updated)
# - forecast_yr: Year to forecast (e.g., 2025)
# - covars_by_age: List of covariates by age class
# - diagnostics_env: Environment for storing model diagnostics
# - model_name: Name of the model (e.g., "LinearCov", "GLM_cov")
# - model_fitting_function: Function that fits the model (takes data, returns fitted model)
# - min_retro_yrs: Minimum years required for retrospective analysis
#
# OUTPUT:
# - Updated fc object with model results
# - Performance metrics added to retro.pm$retro.pm.all.varyrs
# - Forecast predictions added to table.ptf
#
# ERROR METRICS CALCULATED:
# - MRE: Mean Residual Error (bias)
# - MAE: Mean Absolute Error
# - MPE: Mean Percentage Error
# - MAPE: Mean Absolute Percentage Error
# - MASE: Mean Absolute Scaled Error
# - RMSE: Root Mean Square Error
#
# DEPENDENCIES:
# - dplyr, tibble, purrr, abind packages
# - forecastR package for data structures
#==============================================================================

walk_forward_template <- function(
    stock_name,
    raw_data_short,
    fc,
    forecast_yr,
    covars_by_age,
    diagnostics_env,
    model_name = "GenericModel",
    model_fitting_function,  # Function that takes (data, covars) and returns a fitted model
    min_retro_yrs = 10
) {
  
  # Load required libraries
  library(dplyr)
  library(tibble)
  library(purrr)
  library(abind)
  
  #============================================================================
  # 1. PREPARE STOCK DATA
  #============================================================================
  # Extract and clean data for the specific stock
  
  stock_df <- raw_data_short[[stock_name]] %>%
    filter(!is.na(Average_Terminal_Run))  # Remove rows with missing response
  all_ages <- sort(unique(stock_df$Age_Class))  # Get unique age classes
  metrics <- c("MRE", "MAE", "MPE", "MAPE", "MASE", "RMSE")  # Define error metrics
  
  #============================================================================
  # 2. DETERMINE RETROSPECTIVE YEARS BY AGE
  #============================================================================
  # Calculate which years have sufficient data for retrospective validation
  # Each age class needs at least min_retro_yrs of data
  
  retro_yrs_by_age <- lapply(all_ages, function(age_val) {
    age_years <- stock_df %>% filter(Age_Class == age_val) %>% pull(Brood_Year)
    # Only include years that have enough historical data for training
    age_years[age_years >= (min(age_years) + min_retro_yrs)]
  })
  names(retro_yrs_by_age) <- all_ages
  
  #============================================================================
  # 3. ERROR METRICS CALCULATION FUNCTION
  #============================================================================
  # Function to calculate comprehensive error metrics from actual vs predicted values
  
  eval_model <- function(actual, fitted) {
    if (length(actual) == 0 || length(fitted) == 0) {
      return(tibble(MRE = 0, MAE = 0, MPE = 0, MAPE = 0, MASE = 0, RMSE = 0))
    }
    
    residuals <- fitted - actual
    pct_error <- residuals / actual * 100
    denom <- mean(abs(actual - mean(actual)), na.rm = TRUE)
    
    tibble(
      MRE  = mean(residuals, na.rm = TRUE),                    # Mean Residual Error (bias)
      MAE  = mean(abs(residuals), na.rm = TRUE),               # Mean Absolute Error
      MPE  = mean(pct_error, na.rm = TRUE),                    # Mean Percentage Error
      MAPE = mean(abs(pct_error), na.rm = TRUE),               # Mean Absolute Percentage Error
      MASE = if (denom > 0) mean(abs(residuals), na.rm = TRUE) / denom else 0,  # Mean Absolute Scaled Error
      RMSE = sqrt(mean(residuals^2, na.rm = TRUE))             # Root Mean Square Error
    ) %>% round(2)
  }
  
  #============================================================================
  # 4. INITIALIZE STORAGE STRUCTURES
  #============================================================================
  # Set up storage for diagnostics and retrospective predictions
  
  if (!stock_name %in% names(diagnostics_env)) diagnostics_env[[stock_name]] <- list()
  diagnostics_env[[stock_name]][[model_name]] <- list()
  
  retro_preds <- list()  # Store retrospective predictions by year and age
  
  #============================================================================
  # 5. WALK-FORWARD RETROSPECTIVE VALIDATION LOOP
  #============================================================================
  # Core loop: for each age class, perform walk-forward validation
  
  cat("Processing stock:", stock_name, "with", model_name, "\n")
  
  for (age_val in all_ages) {
    retro_yrs_age <- retro_yrs_by_age[[as.character(age_val)]]
    covars <- covars_by_age[[as.character(age_val)]]
    
    # Skip if no covariates available for this age
    if (is.null(covars) || length(covars) == 0) {
      cat("  No covariates for Age", age_val, "\n")
      next
    }
    
    # Get data for this age class and validate covariates
    age_df <- stock_df %>% filter(Age_Class == age_val)
    valid_covars <- covars[covars %in% names(age_df)]
    
    if (length(valid_covars) == 0) {
      cat("  No valid covariates for Age", age_val, "\n")
      next
    }
    
    cat("  Processing Age", age_val, "with covariates:", paste(valid_covars, collapse = ", "), "\n")
    
    # Walk-forward validation for each validation year
    for (validation_yr in retro_yrs_age) {
      # Training data: all years before validation year
      train_data <- age_df %>% 
        filter(Brood_Year < validation_yr) %>%
        dplyr::select(Average_Terminal_Run, all_of(valid_covars)) %>%
        filter(complete.cases(.))  # Remove rows with missing values
      
      # Validation data: just the validation year
      validation_data <- age_df %>% 
        filter(Brood_Year == validation_yr) %>%
        dplyr::select(Average_Terminal_Run, all_of(valid_covars)) %>%
        filter(complete.cases(.))
      
      # Skip if insufficient training data or no validation data
      if (nrow(train_data) < 3 || nrow(validation_data) == 0) next
      
      # Fit model and make prediction
      tryCatch({
        fit <- model_fitting_function(train_data, valid_covars)
        
        # Predict on validation year
        pred_val <- predict(fit, newdata = validation_data)
        actual_val <- validation_data$Average_Terminal_Run
        
        # Store prediction results
        yr_key <- as.character(validation_yr)
        age_key <- paste0("Age ", age_val)
        
        if (!yr_key %in% names(retro_preds)) retro_preds[[yr_key]] <- list()
        
        retro_preds[[yr_key]][[age_key]] <- list(
          actual = actual_val,
          fitted = pred_val,
          run_years = validation_yr
        )
        
        cat("    Year", validation_yr, ": actual =", round(actual_val, 1), ", predicted =", round(pred_val, 1), "\n")
        
      }, error = function(e) {
        cat("    Error fitting model for Age", age_val, "Year", validation_yr, ":", e$message, "\n")
      })
    }
  }
  
  #============================================================================
  # 6. CALCULATE ERROR METRICS BY AGE CLASS
  #============================================================================
  # Aggregate predictions and calculate performance metrics for each age
  
  error_metrics_by_age <- list()
  fixed_ages <- c(paste0("Age ", all_ages), "Total")
  
  for (age_val in all_ages) {
    # Extract all actual and predicted values for this age
    all_actual <- unlist(lapply(retro_preds, function(x) {
      if (!is.null(x[[paste0("Age ", age_val)]])) x[[paste0("Age ", age_val)]]$actual else NULL
    }))
    all_fitted <- unlist(lapply(retro_preds, function(x) {
      if (!is.null(x[[paste0("Age ", age_val)]])) x[[paste0("Age ", age_val)]]$fitted else NULL
    }))
    
    if (length(all_actual) > 0 && length(all_fitted) > 0) {
      cat("Age", age_val, ": calculating metrics for", length(all_actual), "predictions\n")
      error_metrics_by_age[[paste0("Age ", age_val)]] <- eval_model(all_actual, all_fitted)
    }
  }
  
  #============================================================================
  # 7. CALCULATE TOTAL (INTERSECTIONAL) METRICS
  #============================================================================
  # Calculate metrics for total returns across all age classes
  # Uses intersectional approach: only years where all ages have predictions
  
  all_years_per_age <- lapply(all_ages, function(age_val) {
    unique(unlist(lapply(retro_preds, function(x) {
      if (!is.null(x[[paste0("Age ", age_val)]])) x[[paste0("Age ", age_val)]]$run_years else NULL
    })))
  })
  years_used <- Reduce(intersect, all_years_per_age)  # Years with predictions for all ages
  
  if (length(years_used) > 0) {
    total_actual <- numeric()
    total_fitted <- numeric()
    
    for (yr in years_used) {
      yr_key <- as.character(yr)
      if (yr_key %in% names(retro_preds)) {
        pred_list <- retro_preds[[yr_key]]
        
        # Only include years where all age classes have predictions
        if (all(sapply(all_ages, function(a) paste0("Age ", a) %in% names(pred_list)))) {
          sum_actual <- sum(sapply(all_ages, function(a) pred_list[[paste0("Age ", a)]]$actual))
          sum_fitted <- sum(sapply(all_ages, function(a) pred_list[[paste0("Age ", a)]]$fitted))
          
          total_actual <- c(total_actual, sum_actual)
          total_fitted <- c(total_fitted, sum_fitted)
        }
      }
    }
    
    if (length(total_actual) > 0) {
      error_metrics_by_age[["Total"]] <- eval_model(total_actual, total_fitted)
    }
  }
  
  #============================================================================
  # 8. GENERATE FORECAST FOR CURRENT YEAR
  #============================================================================
  # Use all available data to generate forecasts for the current year
  # Since current year won't have covariate data, use latest available year's covariates
  
  preds_by_age <- list()
  
  for (age_val in all_ages) {
    covars <- covars_by_age[[as.character(age_val)]]
    
    if (is.null(covars) || length(covars) == 0) {
      preds_by_age[[paste0("Age ", age_val)]] <- NA
      next
    }
    
    age_df <- stock_df %>% filter(Age_Class == age_val)
    valid_covars <- covars[covars %in% names(age_df)]
    
    if (length(valid_covars) == 0) {
      preds_by_age[[paste0("Age ", age_val)]] <- NA
      next
    }
    
    # Training data: all available years
    train_data <- age_df %>% 
      dplyr::select(Average_Terminal_Run, all_of(valid_covars)) %>%
      filter(complete.cases(.))
    
    # Forecast data: use most recent available covariate data
    # Since forecast_yr won't have covariate data, use the latest available year
    latest_year <- max(age_df$Brood_Year[complete.cases(age_df[valid_covars])])
    forecast_data <- age_df %>% 
      filter(Brood_Year == latest_year) %>%
      dplyr::select(all_of(valid_covars)) %>%
      filter(complete.cases(.))
    
    if (nrow(train_data) < 3 || nrow(forecast_data) == 0) {
      preds_by_age[[paste0("Age ", age_val)]] <- NA
      next
    }
    
    # Fit model on all training data and predict forecast year
    tryCatch({
      fit <- model_fitting_function(train_data, valid_covars)
      pred_val <- predict(fit, newdata = forecast_data)
      preds_by_age[[paste0("Age ", age_val)]] <- round(pred_val)
      cat("  Forecast for Age", age_val, "in", forecast_yr, "(using", latest_year, "covariates):", round(pred_val), "\n")
    }, error = function(e) {
      cat("  Error forecasting Age", age_val, "for", forecast_yr, ":", e$message, "\n")
      preds_by_age[[paste0("Age ", age_val)]] <- NA
    })
  }
  
  # Calculate total forecast
  valid_preds <- preds_by_age[!is.na(preds_by_age)]
  total_pred <- if(length(valid_preds) > 0) sum(unlist(valid_preds), na.rm = TRUE) else NA
  
  #============================================================================
  # 9. ADD FORECAST TO table.ptf
  #============================================================================
  # Add forecast predictions to the forecastR table.ptf structure
  
  if(length(valid_preds) > 0) {
    row <- c(Model = model_name, preds_by_age, Total = total_pred)
    df_tbl <- fc[[stock_name]]$`table.ptf`
    df_new <- dplyr::bind_rows(tibble::rownames_to_column(df_tbl, "Model"), tibble::as_tibble_row(row))
    fc[[stock_name]]$`table.ptf` <- tibble::column_to_rownames(df_new, "Model")
    
    cat("Added", model_name, "forecast to table.ptf - Total:", total_pred, "\n")
  }
  
  #============================================================================
  # 10. BUILD RETRO ARRAY AND APPEND TO fc
  #============================================================================
  # Add performance metrics to the forecastR retro.pm structure
  
  if (length(error_metrics_by_age) > 0) {
    metric_values <- map(error_metrics_by_age, ~ as.numeric(.x[1, metrics]))
    
    # Fill missing ages with zeros to maintain consistent structure
    for (age_name in fixed_ages) {
      if (!age_name %in% names(metric_values)) {
        metric_values[[age_name]] <- rep(0, length(metrics))
      }
    }
    
    # Create matrix and array for forecastR compatibility
    mat <- matrix(
      unlist(metric_values[fixed_ages]),
      nrow = length(metrics),
      dimnames = list(metrics, fixed_ages)
    )
    
    arr_new <- array(mat, dim = c(1, nrow(mat), ncol(mat)),
                     dimnames = list(model_name, metrics, fixed_ages))
    
    # Append to existing forecastR array
    fc_arr <- fc[[stock_name]]$`retro.pm`$`retro.pm.all.varyrs`
    fc[[stock_name]]$`retro.pm`$`retro.pm.all.varyrs` <- abind::abind(fc_arr, arr_new, along = 1)
    
    cat("Successfully added", model_name, "metrics for", stock_name, "\n\n")
  }
  
  # Return updated forecast object
  fc
}
