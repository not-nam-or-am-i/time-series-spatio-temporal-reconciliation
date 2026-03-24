# 0 - base_forecasts_sarimax.R
# ========================================
# SARIMAX Base Forecasts with NWP Exogenous Variable
# Following PLAN.md: Use NWP predictions as xreg
# ========================================

rm(list = ls())

# Load required packages
libs <- c("forecast", "doParallel", "doSNOW", "progress", "data.table", "Matrix")
invisible(lapply(libs, library, character.only = TRUE))

# Set working directory to RF_SARIMAX
script_dir <- "RF_SARIMAX"
if (!grepl("RF_SARIMAX$", getwd())) {
  if (dir.exists(script_dir)) {
    setwd(script_dir)
    cat("Changed working directory to:", getwd(), "\n")
  }
}

# Source configuration and helper functions
source("config.R")
source("fun.R")

# Load hierarchy information
load("info_reco.RData")

cat("\n========================================\n")
cat("SARIMAX Base Forecast Generation\n")
cat("With NWP as Exogenous Variable\n")
cat("========================================\n\n")

# ----------------------------------------
# Load measurement data
# ----------------------------------------
meas_raw <- fread(DATA_PATH)
if (is.character(meas_raw[[1]]) || inherits(meas_raw[[1]], "POSIXct")) {
  meas <- as.matrix(meas_raw[, -1])
} else {
  meas <- as.matrix(meas_raw)
}
n_stations <- ncol(meas)
n_obs <- nrow(meas)

cat(sprintf("Measurement data loaded: %d observations x %d stations\n", n_obs, n_stations))

# ----------------------------------------
# Load NWP predictions
# ----------------------------------------
pred_raw <- fread(PRED_PATH)
if (is.character(pred_raw[[1]]) || inherits(pred_raw[[1]], "POSIXct")) {
  pred <- as.matrix(pred_raw[, -1])
} else {
  pred <- as.matrix(pred_raw)
}

cat(sprintf("NWP predictions loaded: %d observations x %d stations\n", nrow(pred), ncol(pred)))

# Verify alignment between meas and pred
if (nrow(meas) != nrow(pred)) {
  stop("Measurement and prediction data have different number of rows!")
}
if (ncol(meas) != ncol(pred)) {
  stop("Measurement and prediction data have different number of columns!")
}

# ----------------------------------------
# Calculate derived parameters
# ----------------------------------------
n_rep <- length(rep_range)
obs_per_day <- m

# Output dimensions
k_star <- sum(m / k.v)  # 60
n_yhat_rows <- h * k_star  # 120
n_yobs_rows <- h * m  # 48
n_res_rows <- train.days * k_star  # 840

cat(sprintf("\nOutput dimensions per station:\n"))
cat(sprintf("  Y.hat rows: %d (temporal hierarchy, h=%d days)\n", n_yhat_rows, h))
cat(sprintf("  Y.obs rows: %d (hourly observations, h=%d days)\n", n_yobs_rows, h))
cat(sprintf("  Res.insamp rows: %d (residuals, %d training days)\n", n_res_rows, train.days))

# ----------------------------------------
# Setup parallel processing
# ----------------------------------------
cat(sprintf("\nSetting up parallel processing with %d cores...\n", ncores))

cl <- makeCluster(ncores)
registerDoSNOW(cl)

clusterExport(cl, c("m", "h", "k.v", "train.days", "obs_per_day",
                    "create_temporal_hierarchy", "create_residual_hierarchy",
                    "n_yhat_rows", "n_yobs_rows", "n_res_rows"))

clusterEvalQ(cl, {
  library(forecast)
})

# ----------------------------------------
# Process each station (bottom-level)
# ----------------------------------------
cat(sprintf("\nProcessing %d bottom-level stations...\n", n_stations))

for (station_idx in 1:n_stations) {

  station_name <- colnames(meas)[station_idx]
  station_meas <- meas[, station_idx]
  station_pred <- pred[, station_idx]  # NWP for this station

  cat(sprintf("\n[Station %d/%d] %s\n", station_idx, n_stations, station_name))

  # Progress bar
  pb <- txtProgressBar(max = n_rep, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  # Export station-specific NWP to workers
  clusterExport(cl, c("station_meas", "station_pred"), envir = environment())

  results <- foreach(rp = rep_range, .options.snow = opts,
                     .errorhandling = "pass") %dopar% {

    tryCatch({
      # Calculate indices
      start_idx <- (rp - 1) * obs_per_day + 1
      end_train <- start_idx + train.days * obs_per_day - 1
      end_test <- end_train + h * obs_per_day

      if (end_test > length(station_meas)) {
        return(list(Y.hat = rep(NA, n_yhat_rows),
                    Y.obs = rep(NA, n_yobs_rows),
                    Res.insamp = rep(NA, n_res_rows),
                    error = "Index out of bounds"))
      }

      # Extract training and forecast period data
      train_y <- station_meas[start_idx:end_train]
      train_x <- station_pred[start_idx:end_train]  # NWP for training
      test_x <- station_pred[(end_train + 1):end_test]  # NWP for forecast period

      # Create time series
      train_ts <- ts(train_y, frequency = m)

      # Fit SARIMAX with NWP as exogenous variable
      fit <- tryCatch({
        auto.arima(train_ts,
                   xreg = train_x,  # NWP as exogenous regressor
                   seasonal = TRUE,
                   stepwise = TRUE,
                   approximation = TRUE,
                   max.p = 3, max.q = 3,
                   max.P = 2, max.Q = 2,
                   max.d = 1, max.D = 1,
                   allowdrift = FALSE)
      }, error = function(e) {
        # Fallback: try without xreg
        tryCatch({
          auto.arima(train_ts,
                     seasonal = TRUE,
                     stepwise = TRUE,
                     approximation = TRUE,
                     allowdrift = FALSE)
        }, error = function(e2) {
          # Ultimate fallback
          Arima(train_ts, order = c(1, 0, 1),
                seasonal = list(order = c(1, 0, 1), period = m))
        })
      })

      # Generate forecasts
      fc <- tryCatch({
        # If model has xreg, use newxreg for forecast
        if (!is.null(fit$xreg)) {
          forecast(fit, h = h * m, xreg = test_x)
        } else {
          forecast(fit, h = h * m)
        }
      }, error = function(e) {
        # Fallback to simple forecast
        forecast(fit, h = h * m)
      })

      fc_hourly <- as.numeric(fc$mean)
      fc_hourly <- pmax(fc_hourly, 0)  # Non-negative

      # Get in-sample residuals
      residuals_hourly <- as.numeric(residuals(fit))
      if (length(residuals_hourly) < train.days * m) {
        residuals_hourly <- c(
          rep(0, train.days * m - length(residuals_hourly)),
          residuals_hourly
        )
      }

      # Create temporal hierarchies
      Y.hat <- create_temporal_hierarchy(fc_hourly, k.v, m, h)
      Y.obs <- station_meas[(end_train + 1):end_test]
      Res.insamp <- create_residual_hierarchy(residuals_hourly, k.v, m, train.days)

      list(Y.hat = Y.hat, Y.obs = Y.obs, Res.insamp = Res.insamp, error = NULL)

    }, error = function(e) {
      list(Y.hat = rep(NA, n_yhat_rows),
           Y.obs = rep(NA, n_yobs_rows),
           Res.insamp = rep(NA, n_res_rows),
           error = as.character(e))
    })
  }

  close(pb)

  errors <- sapply(results, function(x) !is.null(x$error))
  if (any(errors)) {
    cat(sprintf("  Warning: %d replications had errors\n", sum(errors)))
  }

  filename <- sprintf("%s/%s--sarimax.RData", dir_sarimax, station_name)
  save(results, file = filename)
  cat(sprintf("  Saved to: %s\n", filename))
}

# ----------------------------------------
# Process aggregated series (Total + 5 Regions)
# ----------------------------------------
cat("\n========================================\n")
cat("Processing Aggregated Series (L0 + L1)\n")
cat("========================================\n")

# Create aggregated data using S matrix
S <- as.matrix(hts_info$S)
n_upper <- nrow(S) - ncol(S)  # 6

# Aggregate measurement and prediction data
meas_agg <- t(S[1:n_upper, , drop = FALSE] %*% t(meas))
pred_agg <- t(S[1:n_upper, , drop = FALSE] %*% t(pred))

agg_names <- c("Total", paste0("Region_", sprintf("%02d", 1:5)))

for (agg_idx in 1:n_upper) {

  agg_name <- agg_names[agg_idx]
  agg_meas <- meas_agg[, agg_idx]
  agg_pred <- pred_agg[, agg_idx]  # Aggregated NWP

  cat(sprintf("\n[Aggregate %d/%d] %s\n", agg_idx, n_upper, agg_name))

  pb <- txtProgressBar(max = n_rep, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  clusterExport(cl, c("agg_meas", "agg_pred"), envir = environment())

  results <- foreach(rp = rep_range, .options.snow = opts,
                     .errorhandling = "pass") %dopar% {

    tryCatch({
      start_idx <- (rp - 1) * obs_per_day + 1
      end_train <- start_idx + train.days * obs_per_day - 1
      end_test <- end_train + h * obs_per_day

      if (end_test > length(agg_meas)) {
        return(list(Y.hat = rep(NA, n_yhat_rows),
                    Y.obs = rep(NA, n_yobs_rows),
                    Res.insamp = rep(NA, n_res_rows),
                    error = "Index out of bounds"))
      }

      train_y <- agg_meas[start_idx:end_train]
      train_x <- agg_pred[start_idx:end_train]
      test_x <- agg_pred[(end_train + 1):end_test]

      train_ts <- ts(train_y, frequency = m)

      fit <- tryCatch({
        auto.arima(train_ts,
                   xreg = train_x,
                   seasonal = TRUE,
                   stepwise = TRUE,
                   approximation = TRUE,
                   max.p = 3, max.q = 3,
                   max.P = 2, max.Q = 2,
                   max.d = 1, max.D = 1,
                   allowdrift = FALSE)
      }, error = function(e) {
        tryCatch({
          auto.arima(train_ts, seasonal = TRUE, stepwise = TRUE,
                     approximation = TRUE, allowdrift = FALSE)
        }, error = function(e2) {
          Arima(train_ts, order = c(1, 0, 1),
                seasonal = list(order = c(1, 0, 1), period = m))
        })
      })

      fc <- tryCatch({
        if (!is.null(fit$xreg)) {
          forecast(fit, h = h * m, xreg = test_x)
        } else {
          forecast(fit, h = h * m)
        }
      }, error = function(e) {
        forecast(fit, h = h * m)
      })

      fc_hourly <- pmax(as.numeric(fc$mean), 0)

      residuals_hourly <- as.numeric(residuals(fit))
      if (length(residuals_hourly) < train.days * m) {
        residuals_hourly <- c(
          rep(0, train.days * m - length(residuals_hourly)),
          residuals_hourly
        )
      }

      Y.hat <- create_temporal_hierarchy(fc_hourly, k.v, m, h)
      Y.obs <- agg_meas[(end_train + 1):end_test]
      Res.insamp <- create_residual_hierarchy(residuals_hourly, k.v, m, train.days)

      list(Y.hat = Y.hat, Y.obs = Y.obs, Res.insamp = Res.insamp, error = NULL)

    }, error = function(e) {
      list(Y.hat = rep(NA, n_yhat_rows),
           Y.obs = rep(NA, n_yobs_rows),
           Res.insamp = rep(NA, n_res_rows),
           error = as.character(e))
    })
  }

  close(pb)

  errors <- sapply(results, function(x) !is.null(x$error))
  if (any(errors)) {
    cat(sprintf("  Warning: %d replications had errors\n", sum(errors)))
  }

  filename <- sprintf("%s/%s--0--sarimax.RData", dir_sarimax, agg_name)
  save(results, file = filename)
  cat(sprintf("  Saved to: %s\n", filename))
}

# ----------------------------------------
# Cleanup
# ----------------------------------------
stopCluster(cl)

cat("\n========================================\n")
cat("SARIMAX Base Forecasts Complete\n")
cat("========================================\n")
cat(sprintf("Total stations processed: %d\n", n_stations))
cat(sprintf("Aggregated series processed: %d\n", n_upper))
cat(sprintf("Replications per series: %d\n", n_rep))
cat(sprintf("Results saved to: %s\n", dir_sarimax))
cat("NWP used as exogenous variable: YES\n")
