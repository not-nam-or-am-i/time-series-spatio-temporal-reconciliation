# 0 - base_forecasts_ets.R
# ========================================
# ETS Base Forecasts (Matching Original Paper)
# Uses thief:::th.forecast() to fit ETS at each temporal aggregation level
# L0+L1: Pure ETS at all k levels
# L2: Hybrid ETS+NWP (NWP replaces hourly forecasts)
# ========================================

rm(list = ls())

# Load required packages
libs <- c("thief", "forecast", "doParallel", "doSNOW", "progress", "data.table", "Matrix")
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
cat("ETS Base Forecast Generation\n")
cat("(Matching Original Paper Approach)\n")
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
# Load NWP predictions (used for hybrid at L2)
# ----------------------------------------
pred_raw <- fread(PRED_PATH)
if (is.character(pred_raw[[1]]) || inherits(pred_raw[[1]], "POSIXct")) {
  pred <- as.matrix(pred_raw[, -1])
} else {
  pred <- as.matrix(pred_raw)
}

cat(sprintf("NWP predictions loaded: %d observations x %d stations\n", nrow(pred), ncol(pred)))

# Verify alignment
if (nrow(meas) != nrow(pred)) {
  stop("Measurement and prediction data have different number of rows!")
}

# ----------------------------------------
# Calculate derived parameters
# ----------------------------------------
n_rep <- length(rep_range)
obs_per_day <- m

# Temporal aggregation setup
k.v <- c(1, 2, 3, 4, 6, 8, 12, 24)
k.s <- sum(k.v)  # 60
Mk.v <- m / k.v  # Periods per day at each k level

# Output dimensions
k_star <- k.s  # 60
n_yhat_rows <- h * k_star  # 120
n_yobs_rows <- h * m  # 48
n_res_rows <- train.days * k_star  # 840

cat(sprintf("\nOutput dimensions per station:\n"))
cat(sprintf("  Y.hat rows: %d (temporal hierarchy, h=%d days)\n", n_yhat_rows, h))
cat(sprintf("  Y.obs rows: %d (hourly observations, h=%d days)\n", n_yobs_rows, h))
cat(sprintf("  Res.insamp rows: %d (residuals, %d training days)\n", n_res_rows, train.days))

# ----------------------------------------
# Helper function: Convert thief residuals to flat day-major vector
# ----------------------------------------
convert_thief_residuals <- function(frc_residuals, k.v, m, train_days) {
  # frc_residuals is a list: [[1]]=hourly(k=1), ..., [[8]]=daily(k=24)
  # Output: flat vector in day-major format [day1_k24, day1_k12, ..., day1_k1, day2_k24, ...]
  res_vector <- c()

  for (day in 1:train_days) {
    # k_order descending: k=24, k=12, ..., k=1
    for (k_idx in length(k.v):1) {
      n_periods <- m / k.v[k_idx]
      res_at_level <- as.numeric(frc_residuals[[k_idx]])
      start <- (day - 1) * n_periods + 1
      end_pos <- day * n_periods
      res_vector <- c(res_vector, res_at_level[start:end_pos])
    }
  }

  return(res_vector)
}

# ----------------------------------------
# Helper function: Extract Y.hat from thief forecasts (matching original)
# ----------------------------------------
extract_yhat_from_thief <- function(frc_forecast, k.v, Mk.v, k.s, forecast_horizon) {
  # Constructs Y.hat in day-major format from thief forecast list
  # frc_forecast[[1]] = hourly (k=1), frc_forecast[[8]] = daily (k=24)
  Y.hat <- numeric(k.s * forecast_horizon)

  for (hh in 1:forecast_horizon) {
    tmp <- NULL
    for (k in 1:length(k.v)) {
      # Reverse index to go from k=24 down to k=1
      rev_idx <- length(k.v) - k + 1
      n_per <- Mk.v[rev_idx]
      start <- n_per * (hh - 1) + 1
      end_pos <- n_per * hh
      tmp <- c(tmp, as.numeric(frc_forecast[[rev_idx]])[start:end_pos])
    }
    Y.hat[(k.s * (hh - 1) + 1):(k.s * hh)] <- tmp
  }

  return(Y.hat)
}

# ----------------------------------------
# Setup parallel processing
# ----------------------------------------
cat(sprintf("\nSetting up parallel processing with %d cores...\n", ncores))

cl <- makeCluster(ncores)
registerDoSNOW(cl)

clusterExport(cl, c("m", "h", "k.v", "k.s", "Mk.v", "train.days", "obs_per_day",
                    "convert_thief_residuals", "extract_yhat_from_thief",
                    "n_yhat_rows", "n_yobs_rows", "n_res_rows"))

clusterEvalQ(cl, {
  library(thief)
  library(forecast)
})

# ----------------------------------------
# BLOCK 1: Process aggregated series (L0 + L1) — PURE ETS
# ----------------------------------------
cat("\n========================================\n")
cat("Processing Aggregated Series (L0 + L1) — Pure ETS\n")
cat("========================================\n")

S <- as.matrix(hts_info$S)
n_upper <- nrow(S) - ncol(S)

# Aggregate measurement data using S matrix
meas_agg <- t(S[1:n_upper, , drop = FALSE] %*% t(meas))

agg_names <- c("Total", paste0("Region_", sprintf("%02d", 1:5)))

for (agg_idx in 1:n_upper) {

  agg_name <- agg_names[agg_idx]
  agg_meas <- meas_agg[, agg_idx]

  cat(sprintf("\n[Aggregate %d/%d] %s (Pure ETS)\n", agg_idx, n_upper, agg_name))

  pb <- txtProgressBar(max = n_rep, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  clusterExport(cl, c("agg_meas"), envir = environment())

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

      # Extract training window as ts object
      y <- ts(agg_meas[start_idx:end_train], frequency = m)

      # Create temporal aggregates at all k levels
      yk <- thief::tsaggregates(y, m = m, align = "end")

      # Fit ETS at all aggregation levels using thief
      frc <- tryCatch({
        thief:::th.forecast(aggy = yk, h = m * h, usemodel = "ets",
                           forecastfunction = NULL)
      }, error = function(e) {
        # Fallback: try with simpler approach
        tryCatch({
          thief:::th.forecast(aggy = yk, h = m * h, usemodel = "ets",
                             forecastfunction = NULL)
        }, error = function(e2) {
          return(NULL)
        })
      })

      if (is.null(frc)) {
        return(list(Y.hat = rep(NA, n_yhat_rows),
                    Y.obs = rep(NA, n_yobs_rows),
                    Res.insamp = rep(NA, n_res_rows),
                    error = "thief forecast failed"))
      }

      # Recalculate MSE for bottom level
      frc$mse[1] <- mean(frc$residuals[[1]]^2)

      # Extract Y.hat in day-major format
      Y.hat <- extract_yhat_from_thief(frc$forecast, k.v, Mk.v, k.s, h)

      # Non-negative constraint
      Y.hat <- pmax(Y.hat, 0)

      # Extract observed values
      Y.obs <- agg_meas[(end_train + 1):end_test]

      # Convert residuals to flat day-major vector
      Res.insamp <- convert_thief_residuals(frc$residuals, k.v, m, train.days)

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

  filename <- sprintf("%s/%s--0--ets.RData", dir_ets, agg_name)
  save(results, file = filename)
  cat(sprintf("  Saved to: %s\n", filename))
}

# ----------------------------------------
# BLOCK 2: Process bottom-level stations (L2) — Hybrid ETS+NWP
# ----------------------------------------
cat("\n========================================\n")
cat("Processing Bottom-Level Stations (L2) — Hybrid ETS+NWP\n")
cat("========================================\n")

for (station_idx in 1:n_stations) {

  station_name <- colnames(meas)[station_idx]
  station_meas <- meas[, station_idx]
  station_pred <- pred[, station_idx]

  cat(sprintf("\n[Station %d/%d] %s (ETS+NWP)\n", station_idx, n_stations, station_name))

  pb <- txtProgressBar(max = n_rep, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  clusterExport(cl, c("station_meas", "station_pred"), envir = environment())

  results <- foreach(rp = rep_range, .options.snow = opts,
                     .errorhandling = "pass") %dopar% {

    tryCatch({
      start_idx <- (rp - 1) * obs_per_day + 1
      end_train <- start_idx + train.days * obs_per_day - 1
      end_test <- end_train + h * obs_per_day

      if (end_test > length(station_meas)) {
        return(list(Y.hat = rep(NA, n_yhat_rows),
                    Y.obs = rep(NA, n_yobs_rows),
                    Res.insamp = rep(NA, n_res_rows),
                    error = "Index out of bounds"))
      }

      # Extract training window as ts object
      y <- ts(station_meas[start_idx:end_train], frequency = m)

      # Create temporal aggregates at all k levels
      yk <- thief::tsaggregates(y, m = m, align = "end")

      # Fit ETS at all aggregation levels
      frc <- tryCatch({
        thief:::th.forecast(aggy = yk, h = m * h, usemodel = "ets",
                           forecastfunction = NULL)
      }, error = function(e) {
        return(NULL)
      })

      if (is.null(frc)) {
        return(list(Y.hat = rep(NA, n_yhat_rows),
                    Y.obs = rep(NA, n_yobs_rows),
                    Res.insamp = rep(NA, n_res_rows),
                    error = "thief forecast failed"))
      }

      # KEY: Replace bottom-level (hourly, k=1) with NWP
      # frc$forecast[[1]] is the hourly level
      nwp_forecast <- ts(station_pred[(end_train + 1):end_test], frequency = m)
      frc$forecast[[1]] <- nwp_forecast

      # Replace bottom-level residuals with NWP residuals
      nwp_train <- ts(station_pred[start_idx:end_train], frequency = m)
      frc$residuals[[1]] <- nwp_train - y
      frc$mse[1] <- mean(frc$residuals[[1]]^2)

      # Extract Y.hat in day-major format
      Y.hat <- extract_yhat_from_thief(frc$forecast, k.v, Mk.v, k.s, h)

      # Non-negative constraint
      Y.hat <- pmax(Y.hat, 0)

      # Extract observed values
      Y.obs <- station_meas[(end_train + 1):end_test]

      # Convert residuals to flat day-major vector
      Res.insamp <- convert_thief_residuals(frc$residuals, k.v, m, train.days)

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

  filename <- sprintf("%s/%s--ets.RData", dir_ets, station_name)
  save(results, file = filename)
  cat(sprintf("  Saved to: %s\n", filename))
}

# ----------------------------------------
# Cleanup
# ----------------------------------------
stopCluster(cl)

cat("\n========================================\n")
cat("ETS Base Forecasts Complete\n")
cat("========================================\n")
cat(sprintf("Total stations processed: %d\n", n_stations))
cat(sprintf("Aggregated series processed: %d\n", n_upper))
cat(sprintf("Replications per series: %d\n", n_rep))
cat(sprintf("Results saved to: %s\n", dir_ets))
cat("\nApproach:\n")
cat("  - L0+L1 (Aggregated): Pure ETS at all temporal aggregation levels\n")
cat("  - L2 (Bottom-level): Hybrid ETS+NWP (NWP replaces hourly forecasts)\n")
cat("  - ETS fitted independently at each k level via thief::th.forecast()\n")
