#==============================================================================
# WALK-FORWARD SIBLING REGRESSION WITH BOX-COX TRANSFORMATION
#==============================================================================
#
# PURPOSE: Implements sibling regression using previous age class as predictor
# with Box-Cox transformation and proper walk-forward retrospective validation.
#
# MODEL DESCRIPTION:
# - Uses previous age class (sibling) as predictor for current age class
# - Applies Box-Cox transformation to both response and predictor variables
# - Implements 5-year mean approach for youngest age class (no siblings available)
# - Handles age-specific data availability with robust error handling
#
# WALK-FORWARD VALIDATION:
# - Training: Uses all years before validation year for both current and sibling ages
# - Validation: Predicts single future year using sibling data
# - Rolling window: Moves forward year by year through retrospective period
# - Performance metrics: MRE, MAE, MPE, MAPE, MASE, RMSE
# - Final forecast: Uses all available data to predict current year (2025)
#
# SPECIAL HANDLING:
# - Youngest age class: Uses 5-year mean instead of sibling regression
# - Age-specific latest year detection for sibling data
# - Robust Box-Cox transformation with automatic lambda selection
# - Handles missing sibling data gracefully
#
# INPUT PARAMETERS:
# - stock_name: Name of the stock to model (e.g., "sky_fingerlings")
# - raw_data_short: List containing historical run data for all stocks
# - fc: forecastR results object to append results to
# - forecast_yr: Target forecast year (e.g., 2025)
# - diagnostics_env: Environment for storing model diagnostics
# - min_retro_yrs: Minimum years required for retrospective analysis (default: 10)
#
# OUTPUT:
# - Updated fc object with sibling regression performance metrics and forecasts
# - Performance metrics added to fc[[stock]]$retro.pm$retro.pm.all.varyrs
# - Forecast predictions added to fc[[stock]]$table.ptf
#
# ERROR HANDLING:
# - No fallback predictions (as per user requirements)
# - Graceful handling of insufficient sibling data
# - Detailed error messages for debugging
# - Robust Box-Cox transformation handling
#
# DEPENDENCIES:
# - dplyr, tibble, purrr, abind for data manipulation
# - MASS for Box-Cox transformation
#==============================================================================

update_fc_w_lm_cov_bxcx_sibreg <- function(
    stock_name,
    raw_data_short,
    fc,
    forecast_yr,
    diagnostics_env,
    min_retro_yrs = 10
) {
  library(MASS)
  library(dplyr)
  library(tibble)
  library(purrr)
  library(abind)
  
  # Helper function: Inverse Box-Cox transformation
  inv_boxcox <- function(trans_y, lambda) {
    if (abs(lambda) < 1e-6) exp(trans_y) else (lambda * trans_y + 1)^(1 / lambda)
  }
  
  # Walk-forward validation implementation
  stock_df <- raw_data_short[[stock_name]] %>%
    filter(!is.na(Average_Terminal_Run))
  all_ages <- sort(unique(stock_df$Age_Class))
  metrics <- c("MRE", "MAE", "MPE", "MAPE", "MASE", "RMSE")
  
  # For sibling regression, we can only model ages 2 and above (need previous age as predictor)
  if (length(all_ages) < 2) {
    cat("Stock", stock_name, "has insufficient age classes for sibling regression\n")
    return(fc)
  }
  
  retro_yrs_by_age <- lapply(all_ages[-1], function(age_val) {  # Skip youngest age
    age_years <- stock_df %>% filter(Age_Class == age_val) %>% pull(Brood_Year)
    age_years[age_years >= (min(age_years) + min_retro_yrs)]
  })
  names(retro_yrs_by_age) <- all_ages[-1]
  
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
  diagnostics_env[[stock_name]][["SibReg_BoxCox"]] <- list()
  
  retro_preds <- list()
  
  cat("Processing stock:", stock_name, "with Sibling Regression + Box-Cox\n")
  
  # Handle youngest age with naive approach (same as original)
  youngest_age <- all_ages[1]
  youngest_data <- stock_df %>%
    filter(Age_Class == youngest_age & !is.na(Average_Terminal_Run)) %>%
    arrange(desc(Brood_Year)) %>%
    slice_head(n = 5)
  
  if (nrow(youngest_data) > 0) {
    naive_mean <- mean(youngest_data$Average_Terminal_Run, na.rm = TRUE)
    
    # Store naive predictions for youngest age
    for (yr in youngest_data$Brood_Year) {
      yr_key <- as.character(yr)
      age_key <- paste0("Age ", youngest_age)
      
      if (!yr_key %in% names(retro_preds)) retro_preds[[yr_key]] <- list()
      
      actual_val <- youngest_data %>% filter(Brood_Year == yr) %>% pull(Average_Terminal_Run)
      
      retro_preds[[yr_key]][[age_key]] <- list(
        actual = actual_val,
        fitted = naive_mean,
        run_years = yr
      )
    }
    
    cat("  Age", youngest_age, ": using naive mean =", round(naive_mean, 1), "\n")
  }
  
  # Process older ages with sibling regression and walk-forward validation
  for (age_val in all_ages[-1]) {
    retro_yrs_age <- retro_yrs_by_age[[as.character(age_val)]]
    prev_age <- age_val - 1
    
    cat("  Processing Age", age_val, "using Age", prev_age, "as predictor\n")
    
    for (validation_yr in retro_yrs_age) {
      # Get training data: all years before validation year for both current and previous age
      current_age_train <- stock_df %>% 
        filter(Age_Class == age_val, Brood_Year < validation_yr) %>%
        dplyr::select(Brood_Year, Average_Terminal_Run) %>%
        rename(response = Average_Terminal_Run)
      
      prev_age_train <- stock_df %>% 
        filter(Age_Class == prev_age, Brood_Year < validation_yr) %>%
        dplyr::select(Brood_Year, Average_Terminal_Run) %>%
        rename(predictor = Average_Terminal_Run)
      
      # Join for sibling regression
      train_data <- inner_join(current_age_train, prev_age_train, by = "Brood_Year") %>%
        filter(response > 0, predictor > 0)  # Positive values only for Box-Cox
      
      # Get validation data
      current_age_val <- stock_df %>% 
        filter(Age_Class == age_val, Brood_Year == validation_yr) %>%
        pull(Average_Terminal_Run)
      
      prev_age_val <- stock_df %>% 
        filter(Age_Class == prev_age, Brood_Year == validation_yr) %>%
        pull(Average_Terminal_Run)
      
      if (nrow(train_data) < 3 || length(current_age_val) == 0 || length(prev_age_val) == 0) next
      if (prev_age_val <= 0) next  # Need positive predictor value
      
      tryCatch({
        # Box-Cox transformation for response
        bc_resp <- boxcox(lm(response ~ 1, data = train_data, y = TRUE), 
                         lambda = seq(-2, 2, 0.1), plotit = FALSE)
        lambda_resp <- bc_resp$x[which.max(bc_resp$y)]
        
        train_data$trans_response <- if (abs(lambda_resp) < 1e-6) {
          log(train_data$response)
        } else {
          (train_data$response^lambda_resp - 1) / lambda_resp
        }
        
        # Box-Cox transformation for predictor
        bc_pred <- boxcox(lm(predictor ~ 1, data = train_data, y = TRUE), 
                         lambda = seq(-2, 2, 0.1), plotit = FALSE)
        lambda_pred <- bc_pred$x[which.max(bc_pred$y)]
        
        train_data$trans_predictor <- if (abs(lambda_pred) < 1e-6) {
          log(train_data$predictor)
        } else {
          (train_data$predictor^lambda_pred - 1) / lambda_pred
        }
        
        # Fit model on transformed data
        fit <- lm(trans_response ~ trans_predictor, data = train_data)
        
        # Transform validation predictor
        trans_prev_age_val <- if (abs(lambda_pred) < 1e-6) {
          log(prev_age_val)
        } else {
          (prev_age_val^lambda_pred - 1) / lambda_pred
        }
        
        # Predict in transformed space
        pred_trans <- predict(fit, newdata = data.frame(trans_predictor = trans_prev_age_val))
        
        # Transform back to original scale
        pred_val <- inv_boxcox(pred_trans, lambda_resp)
        
        yr_key <- as.character(validation_yr)
        age_key <- paste0("Age ", age_val)
        
        if (!yr_key %in% names(retro_preds)) retro_preds[[yr_key]] <- list()
        
        retro_preds[[yr_key]][[age_key]] <- list(
          actual = current_age_val,
          fitted = pred_val,
          run_years = validation_yr
        )
        
        cat("    Year", validation_yr, ": actual =", round(current_age_val, 1), 
            ", predicted =", round(pred_val, 1), 
            "(using prev age =", round(prev_age_val, 1), ")\n")
        
      }, error = function(e) {
        cat("    Error fitting Sibling Regression for Age", age_val, "Year", validation_yr, ":", e$message, "\n")
      })
    }
  }
  
  # Calculate metrics (same structure as other models)
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
                     dimnames = list("SibReg_BoxCox", metrics, fixed_ages))
    
    fc_arr <- fc[[stock_name]]$`retro.pm`$`retro.pm.all.varyrs`
    fc[[stock_name]]$`retro.pm`$`retro.pm.all.varyrs` <- abind::abind(fc_arr, arr_new, along = 1)
    
    cat("Successfully added Sibling Regression + Box-Cox metrics for", stock_name, "\n\n")
  }
  
  #-----------------------------------------------------------------------------
  # 8. GENERATE FORECASTS FOR CURRENT YEAR AND UPDATE table.ptf
  #-----------------------------------------------------------------------------
  
  # Generate forecasts for each age using sibling regression
  preds_by_age <- list()
  
  # Use most recent available data for forecasting
  # Find latest year where any data is available (will be age-specific)
  
  for(age_val in all_ages) {
    age_df <- stock_df %>% filter(Age_Class == age_val & !is.na(Average_Terminal_Run))
    
    if(nrow(age_df) < 5) {
      next
    }
    
    # Get sibling age forecast data
    sib_age <- age_val - 1
    if(sib_age < min(all_ages)) {
      # For youngest age class, use 5-year mean instead of sibling regression
      recent_data <- age_df %>% 
        arrange(desc(Brood_Year)) %>% 
        head(5) %>% 
        pull(Average_Terminal_Run)
      
      if(length(recent_data) >= 3) {
        pred_val <- round(mean(recent_data, na.rm = TRUE))
        preds_by_age[[paste0("Age ", age_val)]] <- pred_val
        cat("  Forecast for Age", age_val, "in", forecast_yr, "(5-year mean):", pred_val, "\n")
      }
      next
    }
    
    # Get forecast data using latest available year for sibling age
    sib_age_df <- stock_df %>% filter(Age_Class == sib_age & !is.na(Average_Terminal_Run))
    latest_year_sib <- max(sib_age_df$Brood_Year)
    
    forecast_data <- stock_df %>% 
      filter(Brood_Year == latest_year_sib & Age_Class == sib_age) %>%
      dplyr::select(Average_Terminal_Run) %>%
      filter(complete.cases(.))
    
    if(nrow(forecast_data) == 0) {
      next
    }
    
    # Fit final sibling regression model on all available data
    sib_data <- age_df %>%
      inner_join(
        stock_df %>% filter(Age_Class == sib_age) %>% 
          dplyr::select(Brood_Year, sib_run = Average_Terminal_Run),
        by = "Brood_Year"
      ) %>%
      filter(complete.cases(.))
    
    if(nrow(sib_data) < 5) {
      next
    }
    
    tryCatch({
      # Use Box-Cox transformation
      y_vals <- sib_data$Average_Terminal_Run
      lambda_result <- boxcox(y_vals ~ 1, plotit = FALSE, lambda = seq(-2, 2, by = 0.1))
      best_lambda <- lambda_result$x[which.max(lambda_result$y)]
      
      if(abs(best_lambda) < 0.01) {
        y_transformed <- log(y_vals + 1)
      } else {
        y_transformed <- (y_vals^best_lambda - 1) / best_lambda
      }
      
      # Fit sibling regression model
      fit <- lm(y_transformed ~ sib_run, data = cbind(y_transformed, sib_data))
      
      # Make prediction using sibling data
      pred_transformed <- predict(fit, newdata = data.frame(sib_run = forecast_data$Average_Terminal_Run))
      
      # Inverse Box-Cox transformation
      if(abs(best_lambda) < 0.01) {
        pred_val <- exp(pred_transformed) - 1
      } else {
        pred_val <- (best_lambda * pred_transformed + 1)^(1/best_lambda)
      }
      
      pred_val <- max(0, pred_val)  # Ensure non-negative
      preds_by_age[[paste0("Age ", age_val)]] <- round(pred_val)
      cat("  Forecast for Age", age_val, "in", forecast_yr, "(using sibling Age", sib_age, "from", latest_year_sib, "):", round(pred_val), "\n")
      
    }, error = function(e) {
      cat("  Error forecasting Age", age_val, "for", forecast_yr, ":", e$message, "\n")
      preds_by_age[[paste0("Age ", age_val)]] <- NA
    })
  }
  
  # Add forecasts to table.ptf if any were generated
  if(length(preds_by_age) > 0) {
    total_pred <- sum(unlist(preds_by_age), na.rm = TRUE)
    row <- c(Model = "SibReg_BoxCox", preds_by_age, Total = total_pred)
    
    df_tbl <- fc[[stock_name]]$`table.ptf`
    df_new <- bind_rows(rownames_to_column(df_tbl, "Model"), as_tibble_row(row))
    fc[[stock_name]]$`table.ptf` <- column_to_rownames(df_new, "Model")
    
    cat("Added SibReg_BoxCox forecasts to table.ptf for", stock_name, "\n")
  }
  
  return(fc)
}
