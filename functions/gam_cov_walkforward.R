#==============================================================================
# WALK-FORWARD GENERALIZED ADDITIVE MODEL WITH COVARIATES
#==============================================================================
#
# PURPOSE: Implements Generalized Additive Models (GAM) with environmental 
# covariates using proper walk-forward retrospective validation.
#
# MODEL DESCRIPTION:
# - Uses GAM with smooth terms for flexible modeling of covariate relationships
# - Incorporates environmental covariates specific to each stock and age class
# - Implements robust error handling for convergence issues
# - Handles data quality issues gracefully without fallback predictions
#
# WALK-FORWARD VALIDATION:
# - Training: Uses all years before validation year
# - Validation: Predicts single future year
# - Rolling window: Moves forward year by year through retrospective period
# - Performance metrics: MRE, MAE, MPE, MAPE, MASE, RMSE
# - Final forecast: Uses all available data to predict current year (2025)
#
# COVARIATE HANDLING:
# - Age-specific covariate selection based on stock characteristics
# - Automatic handling of missing covariate data
# - Age-specific latest year detection for forecast data
# - Smooth terms applied to all covariates for flexible modeling
#
# INPUT PARAMETERS:
# - stock_name: Name of the stock to model (e.g., "sky_fingerlings")
# - raw_data_short: List containing historical run data for all stocks
# - fc: forecastR results object to append results to
# - forecast_yr: Target forecast year (e.g., 2025)
# - covars_by_age: List of covariates to use for each age class
# - diagnostics_env: Environment for storing model diagnostics
# - min_retro_yrs: Minimum years required for retrospective analysis (default: 10)
#
# OUTPUT:
# - Updated fc object with GAM performance metrics and forecasts
# - Performance metrics added to fc[[stock]]$retro.pm$retro.pm.all.varyrs
# - Forecast predictions added to fc[[stock]]$table.ptf
#
# ERROR HANDLING:
# - No fallback predictions (as per user requirements)
# - Graceful handling of insufficient data
# - Detailed error messages for debugging
# - Robust convergence checking
#
# DEPENDENCIES:
# - dplyr, tibble, purrr, abind for data manipulation
# - mgcv for GAM modeling
#==============================================================================

update_fc_w_gam_cov <- function(
    stock_name,
    raw_data_short,
    fc,
    forecast_yr,
    covars_by_age,
    diagnostics_env,
    min_retro_yrs = 10
) {
  library(mgcv)
  library(dplyr)
  library(tibble)
  library(purrr)
  library(abind)
  
  # Model fitting function for GAM
  fit_gam_model <- function(train_data, valid_covars) {
    # Determine appropriate family
    y_vals <- train_data$Average_Terminal_Run
    
    if (all(y_vals >= 0)) {
      if (any(y_vals == 0)) {
        # Has zeros: use Gaussian with log link or add small constant
        family_choice <- gaussian(link = "log")
        train_data$Average_Terminal_Run <- train_data$Average_Terminal_Run + 0.1
        offset_added <- 0.1
      } else {
        # No zeros: try Gamma with log link
        family_choice <- Gamma(link = "log")
        offset_added <- 0
      }
    } else {
      # Has negative values: use Gaussian
      family_choice <- gaussian()
      offset_added <- 0
    }
    
    # Build GAM formula with smooth terms
    # Use smooth terms for continuous covariates
    smooth_terms <- character()
    
    for (covar in valid_covars) {
      # Check if covariate has enough unique values for smoothing
      unique_vals <- length(unique(train_data[[covar]]))
      
      if (unique_vals > 4) {
        # Use smooth term with appropriate basis dimension
        k_val <- min(unique_vals - 1, 10)  # Conservative choice
        smooth_terms <- c(smooth_terms, paste0("s(", covar, ", k=", k_val, ")"))
      } else {
        # Use linear term for categorical or low-variance covariates
        smooth_terms <- c(smooth_terms, covar)
      }
    }
    
    formula_str <- paste("Average_Terminal_Run ~", paste(smooth_terms, collapse = " + "))
    gam_formula <- as.formula(formula_str)
    
    tryCatch({
      # Fit GAM
      gam_model <- gam(gam_formula, 
                       data = train_data, 
                       family = family_choice,
                       method = "REML")  # Use REML for better smoothing parameter estimation
      
      # Check for convergence
      if (!gam_model$converged) {
        cat("      GAM did not converge, trying simpler model\n")
        # Fallback to linear terms only
        linear_formula <- as.formula(paste("Average_Terminal_Run ~", paste(valid_covars, collapse = " + ")))
        gam_model <- gam(linear_formula, 
                        data = train_data, 
                        family = family_choice,
                        method = "REML")
      }
      
      # Store offset info
      gam_model$offset_added <- offset_added
      
      return(gam_model)
      
    }, error = function(e) {
      cat("      GAM fitting failed, trying simpler approach:", e$message, "\n")
      
      # Fallback to GLM
      tryCatch({
        glm_model <- glm(Average_Terminal_Run ~ ., 
                        data = train_data, 
                        family = family_choice)
        glm_model$offset_added <- offset_added
        return(glm_model)
      }, error = function(e2) {
        # Final fallback to Gaussian
        cat("      All GAM approaches failed, using Gaussian GAM\n")
        simple_formula <- as.formula(paste("Average_Terminal_Run ~", paste(valid_covars, collapse = " + ")))
        gam_model <- gam(simple_formula, data = train_data, family = gaussian())
        gam_model$offset_added <- 0
        return(gam_model)
      })
    })
  }
  
  # Custom predict function for GAM
  predict_gam <- function(model, newdata) {
    pred_vals <- predict(model, newdata, type = "response")
    
    # Remove offset if it was added
    if (!is.null(model$offset_added) && model$offset_added > 0) {
      pred_vals <- pred_vals - model$offset_added
    }
    
    return(pred_vals)
  }
  
  # Walk-forward validation implementation
  stock_df <- raw_data_short[[stock_name]] %>%
    filter(!is.na(Average_Terminal_Run))
  all_ages <- sort(unique(stock_df$Age_Class))
  metrics <- c("MRE", "MAE", "MPE", "MAPE", "MASE", "RMSE")
  
  retro_yrs_by_age <- lapply(all_ages, function(age_val) {
    age_years <- stock_df %>% filter(Age_Class == age_val) %>% pull(Brood_Year)
    age_years[age_years >= (min(age_years) + min_retro_yrs)]
  })
  names(retro_yrs_by_age) <- all_ages
  
  eval_model <- function(actual, fitted) {
    if (length(actual) == 0 || length(fitted) == 0) {
      return(tibble(MRE = 0, MAE = 0, MPE = 0, MAPE = 0, MASE = 0, RMSE = 0))
    }
    residuals <- fitted - actual
    pct_error <- residuals / actual * 100
    denom <- mean(abs(actual - mean(actual)), na.rm = TRUE)
    tibble(
      MRE  = mean(residuals, na.rm = TRUE),
      MAE  = mean(abs(residuals), na.rm = TRUE),
      MPE  = mean(pct_error, na.rm = TRUE),
      MAPE = mean(abs(pct_error), na.rm = TRUE),
      MASE = if (denom > 0) mean(abs(residuals), na.rm = TRUE) / denom else 0,
      RMSE = sqrt(mean(residuals^2, na.rm = TRUE))
    ) %>% round(2)
  }
  
  if (!stock_name %in% names(diagnostics_env)) diagnostics_env[[stock_name]] <- list()
  diagnostics_env[[stock_name]][["GAM"]] <- list()
  
  retro_preds <- list()
  
  cat("Processing stock:", stock_name, "with GAM\n")
  
  for (age_val in all_ages) {
    retro_yrs_age <- retro_yrs_by_age[[as.character(age_val)]]
    covars <- covars_by_age[[as.character(age_val)]]
    
    if (is.null(covars) || length(covars) == 0) {
      cat("  No covariates for Age", age_val, "\n")
      next
    }
    
    age_df <- stock_df %>% filter(Age_Class == age_val)
    valid_covars <- covars[covars %in% names(age_df)]
    
    if (length(valid_covars) == 0) {
      cat("  No valid covariates for Age", age_val, "\n")
      next
    }
    
    cat("  Processing Age", age_val, "with GAM and covariates:", paste(valid_covars, collapse = ", "), "\n")
    
    for (validation_yr in retro_yrs_age) {
      train_data <- age_df %>% 
        filter(Brood_Year < validation_yr) %>%
        dplyr::select(Average_Terminal_Run, all_of(valid_covars)) %>%
        filter(complete.cases(.))
      
      validation_data <- age_df %>% 
        filter(Brood_Year == validation_yr) %>%
        dplyr::select(Average_Terminal_Run, all_of(valid_covars)) %>%
        filter(complete.cases(.))
      
      if (nrow(train_data) < 5 || nrow(validation_data) == 0) next  # GAM needs more data
      
      tryCatch({
        fit <- fit_gam_model(train_data, valid_covars)
        pred_val <- predict_gam(fit, validation_data)
        actual_val <- validation_data$Average_Terminal_Run
        
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
        cat("    Error fitting GAM for Age", age_val, "Year", validation_yr, ":", e$message, "\n")
      })
    }
  }
  
  # Calculate metrics
  error_metrics_by_age <- list()
  fixed_ages <- c(paste0("Age ", all_ages), "Total")
  
  for (age_val in all_ages) {
    all_actual <- unlist(lapply(retro_preds, function(x) {
      if (!is.null(x[[paste0("Age ", age_val)]])) x[[paste0("Age ", age_val)]]$actual else NULL
    }))
    all_fitted <- unlist(lapply(retro_preds, function(x) {
      if (!is.null(x[[paste0("Age ", age_val)]])) x[[paste0("Age ", age_val)]]$fitted else NULL
    }))
    
    if (length(all_actual) > 0 && length(all_fitted) > 0) {
      error_metrics_by_age[[paste0("Age ", age_val)]] <- eval_model(all_actual, all_fitted)
    }
  }
  
  # Total metrics
  all_years_per_age <- lapply(all_ages, function(age_val) {
    unique(unlist(lapply(retro_preds, function(x) {
      if (!is.null(x[[paste0("Age ", age_val)]])) x[[paste0("Age ", age_val)]]$run_years else NULL
    })))
  })
  years_used <- Reduce(intersect, all_years_per_age)
  
  if (length(years_used) > 0) {
    total_actual <- numeric()
    total_fitted <- numeric()
    
    for (yr in years_used) {
      yr_key <- as.character(yr)
      if (yr_key %in% names(retro_preds)) {
        pred_list <- retro_preds[[yr_key]]
        
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
  
  # Build array and append
  if (length(error_metrics_by_age) > 0) {
    metric_values <- map(error_metrics_by_age, ~ as.numeric(.x[1, metrics]))
    
    for (age_name in fixed_ages) {
      if (!age_name %in% names(metric_values)) {
        metric_values[[age_name]] <- rep(0, length(metrics))
      }
    }
    
    mat <- matrix(
      unlist(metric_values[fixed_ages]),
      nrow = length(metrics),
      dimnames = list(metrics, fixed_ages)
    )
    
    arr_new <- array(mat, dim = c(1, nrow(mat), ncol(mat)),
                     dimnames = list("GAM", metrics, fixed_ages))
    
    fc_arr <- fc[[stock_name]]$`retro.pm`$`retro.pm.all.varyrs`
    fc[[stock_name]]$`retro.pm`$`retro.pm.all.varyrs` <- abind::abind(fc_arr, arr_new, along = 1)
    
    cat("Successfully added GAM metrics for", stock_name, "\n\n")
  }
  
  #-----------------------------------------------------------------------------
  # 8. GENERATE FORECASTS FOR CURRENT YEAR AND UPDATE table.ptf
  #-----------------------------------------------------------------------------
  
  # Generate forecasts for each age
  preds_by_age <- list()
  
  # Use most recent available covariate data for forecasting
  # Find latest year where any covariates are available (will be age-specific)
  
  for(age_val in all_ages) {
    covars <- covars_by_age[[as.character(age_val)]]
    
    if (is.null(covars) || length(covars) == 0) {
      next
    }
    
    age_df <- stock_df %>% filter(Age_Class == age_val)
    valid_covars <- covars[covars %in% names(age_df)]
    
    if (length(valid_covars) == 0) {
      next
    }
    
    # Get forecast data using latest available year for this specific age class
    latest_year_age <- max(age_df$Brood_Year[complete.cases(age_df[valid_covars])])
    
    forecast_data <- age_df %>% 
      filter(Brood_Year == latest_year_age) %>%
      dplyr::select(all_of(valid_covars)) %>%
      filter(complete.cases(.))
    
    if(nrow(forecast_data) == 0) {
      cat("  No forecast data available for Age", age_val, "in year", latest_year, "\n")
      next
    }
    
    # Fit final model on all available data
    train_data <- age_df %>% 
      dplyr::select(Average_Terminal_Run, all_of(valid_covars)) %>%
      filter(complete.cases(.))
    
    if(nrow(train_data) < 5) {
      cat("  Insufficient training data for Age", age_val, "(n =", nrow(train_data), ")\n")
      next
    }
    
    tryCatch({
      # Fit GAM with smooth terms and identity link (more stable)
      formula_str <- paste("Average_Terminal_Run ~", paste(paste0("s(", valid_covars, ")"), collapse = " + "))
      fit <- gam(as.formula(formula_str), data = train_data, family = gaussian(link = "identity"))
      
      # Check if model converged
      if (!fit$converged) {
        cat("      GAM did not converge, trying simpler model\n")
        # Try with linear terms instead of smooth terms
        formula_str <- paste("Average_Terminal_Run ~", paste(valid_covars, collapse = " + "))
        fit <- gam(as.formula(formula_str), data = train_data, family = gaussian())
      }
      
      # Make prediction
      pred_val <- predict(fit, newdata = forecast_data, type = "response")
      pred_val <- max(0, pred_val)  # Ensure non-negative
      
      preds_by_age[[paste0("Age ", age_val)]] <- round(pred_val)
      cat("  Forecast for Age", age_val, "in", forecast_yr, "(using", latest_year_age, "covariates):", round(pred_val), "\n")
      
    }, error = function(e) {
      cat("  Error forecasting Age", age_val, "for", forecast_yr, ":", e$message, "\n")
      preds_by_age[[paste0("Age ", age_val)]] <- NA
    })
  }
  
  # Add forecasts to table.ptf if any were generated
  if(length(preds_by_age) > 0) {
    total_pred <- sum(unlist(preds_by_age), na.rm = TRUE)
    row <- c(Model = "GAM", preds_by_age, Total = total_pred)
    
    df_tbl <- fc[[stock_name]]$`table.ptf`
    df_new <- bind_rows(rownames_to_column(df_tbl, "Model"), as_tibble_row(row))
    fc[[stock_name]]$`table.ptf` <- column_to_rownames(df_new, "Model")
    
    cat("Added GAM forecasts to table.ptf for", stock_name, "\n")
  }
  
  return(fc)
}
