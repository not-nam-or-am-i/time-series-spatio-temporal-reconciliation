# 0-base_forecasts_sarimax.R
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

# for local machine
cl <- makeCluster(ncores)
# for SLURM cluster
# cl <- parallel::makeCluster(ncores, type = "PSOCK")
registerDoSNOW(cl)

clusterExport(cl, c("m", "h", "k.v", "train.days", "obs_per_day",
                    "aggregate_hourly_to_k",
                    "assemble_yhat_from_levels", "assemble_residuals_from_levels",
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
      start_idx <- (rp - 1) * obs_per_day + 1
      end_train <- start_idx + train.days * obs_per_day - 1
      end_test <- end_train + h * obs_per_day

      if (end_test > length(station_meas)) {
        return(list(Y.hat = rep(NA, n_yhat_rows),
                    Y.obs = rep(NA, n_yobs_rows),
                    Res.insamp = rep(NA, n_res_rows),
                    error = "Index out of bounds"))
      }

      train_y_hourly <- station_meas[start_idx:end_train]
      train_x_hourly <- station_pred[start_idx:end_train]
      test_x_hourly <- station_pred[(end_train + 1):end_test]

      k_order <- sort(k.v, decreasing = TRUE)
      fc_list <- list()
      res_list <- list()

      for (kk in k_order) {
        ppd <- m / kk
        n_fc <- h * ppd
        n_train_k <- train.days * ppd

        train_k <- aggregate_hourly_to_k(train_y_hourly, kk)
        train_xk <- aggregate_hourly_to_k(train_x_hourly, kk)
        test_xk <- aggregate_hourly_to_k(test_x_hourly, kk)

        train_ts_k <- ts(train_k, frequency = max(ppd, 1))

        fit_k <- tryCatch({
          auto.arima(train_ts_k,
                     xreg = train_xk,
                     seasonal = ppd >= 2,
                     stepwise = TRUE,
                     approximation = TRUE,
                     max.p = 3, max.q = 3,
                     max.P = 2, max.Q = 2,
                     max.d = 1, max.D = 1,
                     allowdrift = FALSE)
        }, error = function(e) {
          tryCatch({
            auto.arima(train_ts_k, seasonal = ppd >= 2,
                       stepwise = TRUE, approximation = TRUE,
                       allowdrift = FALSE)
          }, error = function(e2) {
            if (ppd >= 2) {
              Arima(train_ts_k, order = c(1, 0, 0),
                    seasonal = list(order = c(0, 0, 0), period = ppd))
            } else {
              Arima(train_ts_k, order = c(1, 0, 0))
            }
          })
        })

        fc_k <- tryCatch({
          if (!is.null(fit_k$xreg)) {
            pmax(as.numeric(forecast(fit_k, h = n_fc, xreg = test_xk)$mean), 0)
          } else {
            pmax(as.numeric(forecast(fit_k, h = n_fc)$mean), 0)
          }
        }, error = function(e) {
          pmax(as.numeric(forecast(fit_k, h = n_fc)$mean), 0)
        })

        res_k <- as.numeric(residuals(fit_k))
        res_k[is.na(res_k)] <- 0
        if (length(res_k) < n_train_k) {
          res_k <- c(rep(0, n_train_k - length(res_k)), res_k)
        } else if (length(res_k) > n_train_k) {
          res_k <- tail(res_k, n_train_k)
        }

        fc_list[[as.character(kk)]] <- fc_k
        res_list[[as.character(kk)]] <- res_k
      }

      Y.hat <- assemble_yhat_from_levels(fc_list, k.v, m, h)
      Y.obs <- station_meas[(end_train + 1):end_test]
      Res.insamp <- assemble_residuals_from_levels(res_list, k.v, m, train.days)

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

S <- as.matrix(hts_info$S)
n_upper <- nrow(S) - ncol(S)

meas_agg <- t(S[1:n_upper, , drop = FALSE] %*% t(meas))
pred_agg <- t(S[1:n_upper, , drop = FALSE] %*% t(pred))

agg_names <- c("Total", paste0("Region_", sprintf("%02d", 1:5)))

for (agg_idx in 1:n_upper) {

  agg_name <- agg_names[agg_idx]
  agg_meas <- meas_agg[, agg_idx]
  agg_pred <- pred_agg[, agg_idx]

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

      train_y_hourly <- agg_meas[start_idx:end_train]
      train_x_hourly <- agg_pred[start_idx:end_train]
      test_x_hourly <- agg_pred[(end_train + 1):end_test]

      k_order <- sort(k.v, decreasing = TRUE)
      fc_list <- list()
      res_list <- list()

      for (kk in k_order) {
        ppd <- m / kk
        n_fc <- h * ppd
        n_train_k <- train.days * ppd

        train_k <- aggregate_hourly_to_k(train_y_hourly, kk)
        train_xk <- aggregate_hourly_to_k(train_x_hourly, kk)
        test_xk <- aggregate_hourly_to_k(test_x_hourly, kk)

        train_ts_k <- ts(train_k, frequency = max(ppd, 1))

        fit_k <- tryCatch({
          auto.arima(train_ts_k,
                     xreg = train_xk,
                     seasonal = ppd >= 2,
                     stepwise = TRUE,
                     approximation = TRUE,
                     max.p = 3, max.q = 3,
                     max.P = 2, max.Q = 2,
                     max.d = 1, max.D = 1,
                     allowdrift = FALSE)
        }, error = function(e) {
          tryCatch({
            auto.arima(train_ts_k, seasonal = ppd >= 2,
                       stepwise = TRUE, approximation = TRUE,
                       allowdrift = FALSE)
          }, error = function(e2) {
            if (ppd >= 2) {
              Arima(train_ts_k, order = c(1, 0, 0),
                    seasonal = list(order = c(0, 0, 0), period = ppd))
            } else {
              Arima(train_ts_k, order = c(1, 0, 0))
            }
          })
        })

        fc_k <- tryCatch({
          if (!is.null(fit_k$xreg)) {
            pmax(as.numeric(forecast(fit_k, h = n_fc, xreg = test_xk)$mean), 0)
          } else {
            pmax(as.numeric(forecast(fit_k, h = n_fc)$mean), 0)
          }
        }, error = function(e) {
          pmax(as.numeric(forecast(fit_k, h = n_fc)$mean), 0)
        })

        res_k <- as.numeric(residuals(fit_k))
        res_k[is.na(res_k)] <- 0
        if (length(res_k) < n_train_k) {
          res_k <- c(rep(0, n_train_k - length(res_k)), res_k)
        } else if (length(res_k) > n_train_k) {
          res_k <- tail(res_k, n_train_k)
        }

        fc_list[[as.character(kk)]] <- fc_k
        res_list[[as.character(kk)]] <- res_k
      }

      Y.hat <- assemble_yhat_from_levels(fc_list, k.v, m, h)
      Y.obs <- agg_meas[(end_train + 1):end_test]
      Res.insamp <- assemble_residuals_from_levels(res_list, k.v, m, train.days)

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
cat("Temporal levels: independent SARIMAX at each k level\n")
