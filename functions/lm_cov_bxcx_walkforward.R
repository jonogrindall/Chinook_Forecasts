# WALK-FORWARD LINEAR COVARIATE MODEL WITH BOX-COX TRANSFORMATION
# Replaces lm_cov_bxcx.R with proper retrospective validation

update_fc_w_lm_cov_bxcx <- function(
    stock_name,
    raw_data_short,
    fc,
    forecast_yr,
    covars_by_age,
    diagnostics_env,
    min_retro_yrs = 10
) {
  library(MASS)
  library(forecast)
  
  # Model fitting function for Box-Cox transformation
  fit_boxcox_model <- function(train_data, valid_covars) {
    # Ensure positive values for Box-Cox transformation
    y_vals <- train_data$Average_Terminal_Run
    if (any(y_vals <= 0)) {
      # Add small constant to make all values positive
      y_vals <- y_vals + abs(min(y_vals)) + 1
    }
    
    # Find optimal lambda for Box-Cox transformation
    tryCatch({
      # Create formula
      formula_str <- paste("y_vals ~", paste(valid_covars, collapse = " + "))
      formula_obj <- as.formula(formula_str)
      
      # Prepare data with transformed response
      bc_data <- train_data
      bc_data$y_vals <- y_vals
      
      # Find optimal lambda
      lambda_fit <- boxcox(lm(formula_obj, data = bc_data), plotit = FALSE)
      lambda_optimal <- lambda_fit$x[which.max(lambda_fit$y)]
      
      # Apply Box-Cox transformation
      if (abs(lambda_optimal) < 1e-6) {
        # Lambda ≈ 0: use log transformation
        bc_data$y_transformed <- log(y_vals)
      } else {
        # Lambda ≠ 0: use Box-Cox transformation
        bc_data$y_transformed <- (y_vals^lambda_optimal - 1) / lambda_optimal
      }
      
      # Fit model with transformed response
      transformed_model <- lm(y_transformed ~ ., 
                             data = bc_data %>% dplyr::select(y_transformed, all_of(valid_covars)))
      
      # Store lambda and original offset for back-transformation
      transformed_model$lambda <- lambda_optimal
      transformed_model$y_offset <- abs(min(train_data$Average_Terminal_Run)) + 1
      
      # Create custom predict method for back-transformation
      transformed_model$predict_original <- function(model, newdata) {
        # Predict on transformed scale
        pred_transformed <- predict(model, newdata)
        
        # Back-transform to original scale
        if (abs(model$lambda) < 1e-6) {
          # Log transformation
          pred_original <- exp(pred_transformed)
        } else {
          # Box-Cox transformation
          pred_original <- (model$lambda * pred_transformed + 1)^(1/model$lambda)
        }
        
        # Remove the offset that was added
        pred_original <- pred_original - model$y_offset
        return(pred_original)
      }
      
      return(transformed_model)
      
    }, error = function(e) {
      # Fallback to regular linear model if Box-Cox fails
      cat("      Box-Cox transformation failed, using linear model:", e$message, "\n")
      return(lm(Average_Terminal_Run ~ ., data = train_data))
    })
  }
  
  # Custom predict function that handles Box-Cox back-transformation
  predict_boxcox <- function(model, newdata) {
    if (!is.null(model$predict_original)) {
      return(model$predict_original(model, newdata))
    } else {
      return(predict(model, newdata))
    }
  }
  
  # Modified walk-forward template for Box-Cox
  library(dplyr)
  library(tibble)
  library(purrr)
  library(abind)
  
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
  diagnostics_env[[stock_name]][["LinearCov_BoxCox"]] <- list()
  
  retro_preds <- list()
  
  cat("Processing stock:", stock_name, "with LinearCov_BoxCox\n")
  
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
    
    cat("  Processing Age", age_val, "with Box-Cox and covariates:", paste(valid_covars, collapse = ", "), "\n")
    
    for (validation_yr in retro_yrs_age) {
      train_data <- age_df %>% 
        filter(Brood_Year < validation_yr) %>%
        dplyr::select(Average_Terminal_Run, all_of(valid_covars)) %>%
        filter(complete.cases(.))
      
      validation_data <- age_df %>% 
        filter(Brood_Year == validation_yr) %>%
        dplyr::select(Average_Terminal_Run, all_of(valid_covars)) %>%
        filter(complete.cases(.))
      
      if (nrow(train_data) < 3 || nrow(validation_data) == 0) next
      
      tryCatch({
        fit <- fit_boxcox_model(train_data, valid_covars)
        pred_val <- predict_boxcox(fit, validation_data)
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
        cat("    Error fitting Box-Cox model for Age", age_val, "Year", validation_yr, ":", e$message, "\n")
      })
    }
  }
  
  # Calculate metrics (same as template)
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
                     dimnames = list("LinearCov_BoxCox", metrics, fixed_ages))
    
    fc_arr <- fc[[stock_name]]$`retro.pm`$`retro.pm.all.varyrs`
    fc[[stock_name]]$`retro.pm`$`retro.pm.all.varyrs` <- abind::abind(fc_arr, arr_new, along = 1)
    
    cat("Successfully added LinearCov_BoxCox metrics for", stock_name, "\n\n")
  }
  
  #-----------------------------------------------------------------------------
  # 8. GENERATE FORECASTS FOR CURRENT YEAR AND UPDATE table.ptf
  #-----------------------------------------------------------------------------
  
  # Generate forecasts for each age
  preds_by_age <- list()
  
  # Use most recent available covariate data for forecasting
  # Find latest year where any covariates are available (will be age-specific)
  
  for(age_val in all_ages) {
    # Get covariates for this age
    covars <- covars_by_age[[as.character(age_val)]]
    if (is.null(covars) || length(covars) == 0) {
      next
    }
    
    age_df <- stock_df %>% filter(Age_Class == age_val & !is.na(Average_Terminal_Run))
    valid_covars <- covars[covars %in% names(age_df)]
    
    if (length(valid_covars) == 0) {
      next
    }
    
    if(nrow(age_df) < 5) {
      cat("  Skipping Age", age_val, "- insufficient data\n")
      next
    }
    
    # Get forecast data using latest available year for this specific age class
    latest_year_age <- max(age_df$Brood_Year[complete.cases(age_df[valid_covars])])
    
    forecast_data <- age_df %>% 
      filter(Brood_Year == latest_year_age) %>%
      dplyr::select(all_of(valid_covars)) %>%
      filter(complete.cases(.))
    
    if(nrow(forecast_data) == 0) {
      cat("  Skipping Age", age_val, "- no forecast data available\n")
      next
    }
    
    # Fit final model on all available data
    train_data <- age_df %>% 
      dplyr::select(Average_Terminal_Run, all_of(valid_covars)) %>%
      filter(complete.cases(.))
    
    tryCatch({
      # Use Box-Cox transformation
      y_vals <- train_data$Average_Terminal_Run
      
      # Add small offset if needed for Box-Cox
      if (any(y_vals <= 0)) {
        offset <- abs(min(y_vals)) + 1
        y_vals <- y_vals + offset
      } else {
        offset <- 0
      }
      
      # Fit basic model for Box-Cox
      basic_formula <- as.formula(paste("y_vals ~", paste(valid_covars, collapse = " + ")))
      temp_data <- train_data
      temp_data$y_vals <- y_vals
      
      # Box-Cox transformation
      bc_result <- boxcox(lm(basic_formula, data = temp_data), plotit = FALSE)
      lambda_optimal <- bc_result$x[which.max(bc_result$y)]
      
      # Transform response
      if (abs(lambda_optimal) < 1e-6) {
        y_transformed <- log(y_vals)
      } else {
        y_transformed <- (y_vals^lambda_optimal - 1) / lambda_optimal
      }
      
      train_data_transformed <- train_data
      train_data_transformed$Average_Terminal_Run <- y_transformed
      
      # Fit model
      formula_str <- paste("Average_Terminal_Run ~", paste(valid_covars, collapse = " + "))
      fit <- lm(as.formula(formula_str), data = train_data_transformed)
      
      # Make prediction
      pred_transformed <- predict(fit, newdata = forecast_data)
      
      # Inverse Box-Cox transformation
      if (abs(lambda_optimal) < 1e-6) {
        pred_val <- exp(pred_transformed)
      } else {
        pred_val <- (lambda_optimal * pred_transformed + 1)^(1/lambda_optimal)
      }
      
      # Remove offset
      pred_val <- pred_val - offset
      pred_val <- max(0, pred_val)  # Ensure non-negative
      preds_by_age[[paste0("Age ", age_val)]] <- round(pred_val)
      cat("  Forecast for Age", age_val, "in", forecast_yr, "(using", latest_year, "covariates):", round(pred_val), "\n")
      
    }, error = function(e) {
      cat("  Error forecasting Age", age_val, "for", forecast_yr, ":", e$message, "\n")
      preds_by_age[[paste0("Age ", age_val)]] <- NA
    })
  }
  
  # Add forecasts to table.ptf if any were generated
  if(length(preds_by_age) > 0) {
    total_pred <- sum(unlist(preds_by_age), na.rm = TRUE)
    row <- c(Model = "LinearCov_BoxCox", preds_by_age, Total = total_pred)
    
    df_tbl <- fc[[stock_name]]$`table.ptf`
    df_new <- bind_rows(rownames_to_column(df_tbl, "Model"), as_tibble_row(row))
    fc[[stock_name]]$`table.ptf` <- column_to_rownames(df_new, "Model")
    
    cat("Added LinearCov_BoxCox forecasts to table.ptf for", stock_name, "\n")
  }
  
  return(fc)
}
