setwd("R:/Fisheries Management/forecasting/2025/code_dev/all_forecasts")

# --- Load Libraries ---
library(forecastR)
library(tibble)
library(dplyr)
library(purrr)
library(DBI)
library(odbc)

forecast_yr = 2025

source("functions/cov_df.R")
source("functions/run_fc_models.R")
source("functions/get_ranked_fc.R") 
source("functions/diagnostic_plots.R")

con <- dbConnect(odbc(),
                 Driver = "SQL Server",
                 Server = "SVR-WEBSQLCLUSTER",
                 Database = "nrfss",
                 Trusted_Connection = "Yes")

source("functions/age_structure.R")
source("functions/terminal_run_df.R")

raw_data <- term_run_df(
  forecast_yr = forecast_yr
)

dbDisconnect(con)


cov_data <- cov_df()
  #names(cov_data)

# --- Create Short df With Covariates ----
cov_map <- list(
  sky_fingerlings = c("Cov_pdo_summer", "Cov_NPGO_Summer", "Cov_pc1", "Cov_pc2", "Cov_sky_qmax"),
  sky_yearlings = c("Cov_pdo_summer", "Cov_NPGO_Summer", "Cov_pc1", "Cov_pc2", "Cov_sky_qmax"),
  snq_fingerlings = c("Cov_pdo_summer", "Cov_NPGO_Summer", "Cov_pc1", "Cov_pc2", "Cov_snq_qmax"),
  snq_yearlings = c("Cov_pdo_summer", "Cov_NPGO_Summer", "Cov_pc1", "Cov_pc2","Cov_snq_qmax"),
  tul_fingerlings = c("Cov_pdo_summer", "Cov_NPGO_Summer", "Cov_pc1", "Cov_pc2")
)

raw_data_short <- lapply(raw_data, function(df) {
  df %>% dplyr::select(Run_Year, Brood_Year, Age_Class, Average_Terminal_Run)
})

raw_data_short <- raw_data_short %>%
  add_covs_per_stock(cov_data, cov_map)

# --- Create ForecastR df With Covariates ----
raw_data_forecastr <- raw_data
    # names(cov_data)
cov_map_forecastr <- list(
  sky_fingerlings = c("Cov_NPGO_Summer", "Cov_pc1", "Cov_sky_qmax"),
  sky_yearlings = c("Cov_NPGO_Summer", "Cov_pc1", "Cov_sky_qmax"),
  snq_fingerlings = c("Cov_pdo_summer", "Cov_NPGO_Summer", "Cov_pc1", "Cov_snq_qmax"),
  snq_yearlings = c("Cov_NPGO_Summer", "Cov_snq_qmax"),
  tul_fingerlings = c("Cov_NPGO_Summer", "Cov_pc1")
)


# --- Join Covariates Data ---
raw_data_forecastr <- raw_data_forecastr %>%
  add_covs_per_stock(cov_data, cov_map_forecastr)

# --- Define Model Settings ---
model.settings <- list(
  Naive3 = list(model.type = "Naive", settings = list(avg.yrs = 3)),
  Naive5 = list(model.type = "Naive", settings = list(avg.yrs = 5)),
  SibRegSimple = list(model.type = "SibRegSimple", settings = NULL),
  SibRegLogPower = list(model.type = "SibRegLogPower", settings = NULL),
  SibRegKalman = list(model.type = "SibRegKalman", settings = NULL)#,
  # SibRegComplex = list(
  #   model.type = "SibRegComplex",
  #   settings = list(
  #     tol.AIC = 0.75,
  #     tol.r.sq = 0.02,
  #     incl.base.eq = FALSE
  #   )
  # )
)

fc <- lapply(names(raw_data_forecastr), function(name) {
  cat("Running multiFC on:", name, "\n")
  
  # Convert tibble to data.frame
  df <- as.data.frame(raw_data_forecastr[[name]])
  
  # --- Run Forecasts ---
  fc <- multiFC(
    data.file = df,
    settings.list = model.settings,
    retro.min.yrs = 14,
    out.type = "full",
    int.type = "Retrospective",
    do.retro = TRUE,
    tracing = TRUE
  )
})

# Name the list with original tibble names
names(fc) <- names(raw_data_forecastr)

# source("functions/explore_cov_cor.R")
# explore_cov_cor(raw_data_short)

# this was chose based on the above explore_cov_cor plots
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


stock_names <- names(fc)  
diagnostics_env <- new.env(parent = emptyenv())


fc <- run_fc_models(forecast_yr)

ranked_forecasts <- get_ranked_fc(fc)
# ranked_forecasts



# ## Diagnostic Plot Options 
# # generate all plots
# make_diagnostic_plots(save = TRUE)
# # generate specific plots example (stock, model, and age names found in ranked_forecasts)
# diagnostic_plot_selection <- list(
#   plot1  = list(stock = "sky_fingerlings", model = "LinearCov", age = "Age 4"),
#   plot2  = list(stock = "sky_fingerlings", model = "Naive3", age = "Age 2")
# )
# 
# make_diagnostic_plots(plot_targets = diagnostic_plot_selection, save = TRUE)
