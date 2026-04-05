# fun.R
# ========================================
# Helper Functions for RF_SARIMAX
# ========================================

#' Load replication data for a specific forecast method
#'
#' Loads and reorganizes forecast data from Results_SARIMAX or Results_RF
#'
#' @param i Replication index (1-350)
#' @param method Character: "sarimax" or "rf"
#' @param bottom_only Logical: if TRUE, return only bottom-level series (318)
#'
#' @return List with Yhat, B1 (observations), E (residuals)
load_replication <- function(i = 1, method = "sarimax", bottom_only = FALSE) {

  if (i < 1 | i > 350) {
    stop("Replication index must be between 1 and 350", call. = FALSE)
  }

  # Determine directory based on method
  if (method == "sarimax") {
    dir_path <- "./Results_SARIMAX"
  } else if (method == "rf") {
    dir_path <- "./Results_RF"
  } else if (method == "rf_nwp") {
    dir_path <- "./Results_RF_NWP"
  } else if (method == "lgbm") {
    dir_path <- "./Results_LGBM"
  } else if (method == "ets") {
    dir_path <- "./Results_ETS"
  } else if (method == "ets_author") {
    dir_path <- "./Results_ETS_Author"
  } else if (method == "sarimax_nwp") {
    dir_path <- "./Results_SARIMAX_NWP"
  } else {
    stop("Method must be 'sarimax', 'rf', 'rf_nwp', 'lgbm', 'ets', 'ets_author', or 'sarimax_nwp'",
         call. = FALSE)
  }

  filesn <- list.files(dir_path, pattern = "\\.RData$", full.names = TRUE)

  if (length(filesn) == 0) {
    stop(sprintf("No results files found in %s", dir_path), call. = FALSE)
  }

  # Initialize matrices
  # Dimensions:
  # - Yhat: [k* x h] rows x n columns = 120 x 324
  # - B1 (obs): [m x h] rows x n columns = 48 x 324
  # - E (residuals): [N x k*] rows x n columns

  k <- c(24, 12, 8, 6, 4, 3, 2, 1)  # Temporal aggregation order
  m <- 24
  h <- 2

  n_bottom <- 318
  n_agg <- 6  # Total + 5 regions
  n_total <- n_bottom + n_agg

  if (bottom_only) {
    n_cols <- n_bottom
  } else {
    n_cols <- n_total
  }

  # Output dimensions
  Yhat <- matrix(NA, nrow = 120, ncol = n_cols)
  Yobs <- matrix(NA, nrow = 48, ncol = n_cols)
  Res <- matrix(NA, nrow = 840, ncol = n_cols)  # 14 days * 60 (k*) = 840

  col_idx <- 1

  for (file in filesn) {
    # Load the results
    load(file)  # loads 'results' list

    # Extract station info from filename
    station_name <- gsub("--.*\\.RData$", "", basename(file))

    # Get data for replication i
    rep_data <- results[[i]]

    if (!is.null(rep_data)) {
      Yhat[, col_idx] <- rep_data$Y.hat
      Yobs[, col_idx] <- rep_data$Y.obs
      Res[, col_idx] <- rep_data$Res.insamp
      colnames(Yhat)[col_idx] <- station_name
      col_idx <- col_idx + 1
    }
  }

  # Trim to actual data loaded
  Yhat <- Yhat[, 1:(col_idx - 1), drop = FALSE]
  Yobs <- Yobs[, 1:(col_idx - 1), drop = FALSE]
  Res <- Res[, 1:(col_idx - 1), drop = FALSE]

  return(list(
    Yhat = Yhat,
    B1 = Yobs,
    E = Res
  ))
}


#' Drop near-zero values for numerical stability
#'
#' @param x Numeric vector or matrix
#' @param tol Tolerance threshold
#'
#' @return Cleaned vector/matrix with near-zero values set to zero
drop_zeros <- function(x, tol = sqrt(.Machine$double.eps)) {
  if (is.vector(x) && !is.matrix(x)) {
    x[abs(x) < tol] <- 0
    return(x)
  }
  x1 <- as.matrix(Matrix::drop0(x, tol = tol))
  return(x1)
}


#' Extract specific temporal aggregation level
#'
#' @param data Matrix with rows organized by temporal hierarchy
#' @param aggregation Aggregation level (24, 12, 8, 6, 4, 3, 2, or 1)
#'
#' @return Subset of data for requested aggregation level
extract_agg <- function(data, aggregation = 1) {
  agg <- c(24, 12, 8, 6, 4, 3, 2, 1)
  h <- NROW(data) / sum(agg)
  return(data[rep(agg, h * max(agg) / agg) == aggregation, , drop = FALSE])
}


#' Create temporal hierarchy from hourly forecasts
#'
#' Aggregates hourly forecasts to all temporal levels
#'
#' @param fc_hourly Numeric vector of hourly forecasts (length = m * h = 48)
#' @param k.v Temporal aggregation levels (default: c(1,2,3,4,6,8,12,24))
#' @param m Seasonal period (default: 24)
#' @param h Forecast horizon in days (default: 2)
#'
#' @return Numeric vector organized by temporal hierarchy
create_temporal_hierarchy <- function(fc_hourly, k.v = c(1,2,3,4,6,8,12,24),
                                      m = 24, h = 2) {

  # Total hourly forecasts
  n_hourly <- length(fc_hourly)  # Should be m * h = 48

  if (n_hourly != m * h) {
    stop(sprintf("Expected %d hourly forecasts, got %d", m * h, n_hourly))
  }

  # Create list to store aggregated forecasts
  # Order: k=24 (daily), k=12, k=8, k=6, k=4, k=3, k=2, k=1 (hourly)
  k_order <- sort(k.v, decreasing = TRUE)  # 24, 12, 8, 6, 4, 3, 2, 1

  result <- c()

  for (day in 1:h) {
    # Extract hourly forecasts for this day
    start_idx <- (day - 1) * m + 1
    end_idx <- day * m
    day_hourly <- fc_hourly[start_idx:end_idx]

    # Aggregate to each level
    for (k in k_order) {
      n_periods <- m / k  # Number of periods at this aggregation

      if (k == 1) {
        # Hourly - no aggregation needed
        agg_fc <- day_hourly
      } else {
        # Aggregate by summing consecutive k hours
        agg_fc <- numeric(n_periods)
        for (p in 1:n_periods) {
          idx_start <- (p - 1) * k + 1
          idx_end <- p * k
          agg_fc[p] <- sum(day_hourly[idx_start:idx_end])
        }
      }

      result <- c(result, agg_fc)
    }
  }

  return(result)
}


#' Create residual hierarchy from in-sample residuals
#'
#' Organizes residuals by temporal aggregation levels
#'
#' @param residuals Numeric vector of hourly residuals
#' @param k.v Temporal aggregation levels
#' @param m Seasonal period
#' @param train_days Number of training days
#'
#' @return Numeric vector organized by temporal hierarchy
create_residual_hierarchy <- function(residuals, k.v = c(1,2,3,4,6,8,12,24),
                                      m = 24, train_days = 14) {

  k_order <- sort(k.v, decreasing = TRUE)
  n_hourly <- length(residuals)

  result <- c()

  # Process each training day
  for (day in 1:train_days) {
    start_idx <- (day - 1) * m + 1
    end_idx <- min(day * m, n_hourly)

    if (start_idx > n_hourly) break

    day_res <- residuals[start_idx:end_idx]

    # Pad if necessary
    if (length(day_res) < m) {
      day_res <- c(day_res, rep(0, m - length(day_res)))
    }

    # Aggregate to each level
    for (k in k_order) {
      n_periods <- m / k

      if (k == 1) {
        agg_res <- day_res
      } else {
        agg_res <- numeric(n_periods)
        for (p in 1:n_periods) {
          idx_start <- (p - 1) * k + 1
          idx_end <- p * k
          agg_res[p] <- sum(day_res[idx_start:idx_end])
        }
      }

      result <- c(result, agg_res)
    }
  }

  return(result)
}


# ========================================
# Metric Calculation Functions
# ========================================

#' Calculate Normalized RMSE
#'
#' @param actual Numeric vector of actual values
#' @param forecast Numeric vector of forecast values
#'
#' @return nRMSE as percentage
calc_nRMSE <- function(actual, forecast) {
  mse <- mean((actual - forecast)^2, na.rm = TRUE)
  rmse <- sqrt(mse)
  mean_actual <- mean(actual, na.rm = TRUE)

  if (mean_actual == 0) return(NA)

  return(rmse / mean_actual * 100)
}


#' Calculate Normalized Mean Bias Error
#'
#' @param actual Numeric vector of actual values
#' @param forecast Numeric vector of forecast values
#'
#' @return nMBE as percentage (positive = over-forecasting)
calc_nMBE <- function(actual, forecast) {
  mbe <- mean(forecast - actual, na.rm = TRUE)
  mean_actual <- mean(actual, na.rm = TRUE)

  if (mean_actual == 0) return(NA)

  return(mbe / mean_actual * 100)
}


#' Calculate Symmetric Mean Absolute Percentage Error
#'
#' @param actual Numeric vector of actual values
#' @param forecast Numeric vector of forecast values
#'
#' @return SMAPE as percentage (0-200 scale)
calc_SMAPE <- function(actual, forecast) {
  denominator <- (abs(actual) + abs(forecast))

  # Avoid division by zero
  valid <- denominator > 0

  if (sum(valid) == 0) return(NA)

  smape <- 100 * mean(abs(actual[valid] - forecast[valid]) / denominator[valid])

  return(smape)
}


#' Calculate Forecast Skill Score
#'
#' Vectorized function for use with dplyr mutate
#'
#' @param rmse_method RMSE of the forecast method (scalar or vector)
#' @param rmse_baseline RMSE of the baseline (scalar or vector)
#'
#' @return Skill score as percentage (positive = better than baseline)
calc_skill <- function(rmse_method, rmse_baseline) {
  # Vectorized calculation
  result <- (1 - rmse_method / rmse_baseline) * 100

  # Handle zero baseline (avoid division by zero)
  result[rmse_baseline == 0] <- NA

  return(result)
}


# ========================================
# Persistence (PERS) Baseline Functions
# ========================================

#' Create persistence forecast
#'
#' Uses the last h days of training as forecast for the next h days
#' (matching the original approach - actual patterns, not repeated)
#'
#' @param data Time series data (vector or single-column matrix)
#' @param train_end Index of last training observation
#' @param h Forecast horizon (days)
#' @param m Seasonal period (hours per day)
#'
#' @return Numeric vector of persistence forecasts (length = h * m)
create_pers_forecast <- function(data, train_end, h = 2, m = 24) {
  # Persistence: use actual last h days of training (not repeating one day)
  # For h=2: uses days (train_end - 47) to train_end as forecast
  last_hdays_start <- train_end - h * m + 1
  pers_fc <- data[last_hdays_start:train_end]

  return(pers_fc)
}


# ========================================
# Multi-Level Temporal Forecasting Helpers
# ========================================

#' Aggregate hourly data to k-hour blocks by summation
aggregate_hourly_to_k <- function(hourly_data, k) {
  if (k == 1) return(hourly_data)
  n <- length(hourly_data)
  n_periods <- n %/% k
  agg <- numeric(n_periods)
  for (p in 1:n_periods) {
    agg[p] <- sum(hourly_data[((p - 1) * k + 1):(p * k)])
  }
  agg
}

#' Assemble Y.hat in day-major format from per-level forecasts
#'
#' @param fc_list Named list: fc_list[["k"]] = forecast vector of length h*(m/k)
#' @param k.v Temporal aggregation levels
#' @param m Seasonal period
#' @param h Forecast horizon in days
assemble_yhat_from_levels <- function(fc_list, k.v, m, h) {
  k_order <- sort(k.v, decreasing = TRUE)
  k_star <- sum(m / k.v)
  Y.hat <- numeric(h * k_star)
  idx <- 1
  for (day in 1:h) {
    for (kk in k_order) {
      n_periods <- m / kk
      fc_k <- fc_list[[as.character(kk)]]
      start_fc <- (day - 1) * n_periods + 1
      end_fc <- day * n_periods
      Y.hat[idx:(idx + n_periods - 1)] <- fc_k[start_fc:end_fc]
      idx <- idx + n_periods
    }
  }
  Y.hat
}

#' Assemble Res.insamp in day-major format from per-level residuals
#'
#' @param res_list Named list: res_list[["k"]] = residual vector of length train_days*(m/k)
#' @param k.v Temporal aggregation levels
#' @param m Seasonal period
#' @param train_days Number of training days
assemble_residuals_from_levels <- function(res_list, k.v, m, train_days) {
  k_order <- sort(k.v, decreasing = TRUE)
  k_star <- sum(m / k.v)
  Res <- numeric(train_days * k_star)
  idx <- 1
  for (day in 1:train_days) {
    for (kk in k_order) {
      n_periods <- m / kk
      res_k <- res_list[[as.character(kk)]]
      start_r <- (day - 1) * n_periods + 1
      end_r <- day * n_periods
      Res[idx:(idx + n_periods - 1)] <- res_k[start_r:end_r]
      idx <- idx + n_periods
    }
  }
  Res
}

#' Create ML features at temporal resolution k
#'
#' @param data Aggregated data at level k
#' @param k Temporal aggregation level
#' @param m Seasonal period (hours per day)
#' @param nwp_data NWP data aggregated at level k (NULL to omit)
#' @param start_period Absolute period index (1-indexed) of first element
create_ml_features_at_k <- function(data, k, m, nwp_data = NULL, start_period = 1) {
  n <- length(data)
  ppd <- m / k

  features <- data.frame(
    lag_1 = c(NA, data[1:(n - 1)])
  )

  if (ppd <= n) {
    features$lag_1day <- c(rep(NA, ppd), data[1:(n - ppd)])
  } else {
    features$lag_1day <- rep(NA, n)
  }

  if (2 * ppd <= n) {
    features$lag_2days <- c(rep(NA, 2 * ppd), data[1:(n - 2 * ppd)])
  } else {
    features$lag_2days <- rep(NA, n)
  }

  if (7 * ppd <= n) {
    features$lag_1week <- c(rep(NA, 7 * ppd), data[1:(n - 7 * ppd)])
  } else {
    features$lag_1week <- rep(NA, n)
  }

  if (!is.null(nwp_data)) {
    features$nwp <- nwp_data
  }

  rolling_mean <- rep(NA, n)
  rolling_sd <- rep(NA, n)
  if (ppd >= 2) {
    for (i in ppd:n) {
      w <- data[(i - ppd + 1):i]
      rolling_mean[i] <- mean(w, na.rm = TRUE)
      rolling_sd[i] <- sd(w, na.rm = TRUE)
    }
  } else {
    rolling_mean <- data
    rolling_sd <- rep(0, n)
  }
  features$rolling_mean <- rolling_mean
  features$rolling_sd <- rolling_sd

  pod <- ((start_period - 1 + 0:(n - 1)) %% ppd)
  features$period_of_day <- pod

  day_idx <- ((start_period - 1 + 0:(n - 1)) %/% ppd) + 1
  features$day_of_week <- (day_idx - 1) %% 7
  features$month <- ((day_idx - 1) %/% 30) %% 12 + 1

  if (ppd > 1) {
    features$period_sin <- sin(2 * pi * pod / ppd)
    features$period_cos <- cos(2 * pi * pod / ppd)
  } else {
    features$period_sin <- rep(0, n)
    features$period_cos <- rep(1, n)
  }

  features$target <- data
  features
}

#' Create single-step ML prediction features at temporal resolution k
#'
#' @param current_data History vector at level k (training + already-predicted)
#' @param k Temporal aggregation level
#' @param m Seasonal period
#' @param nwp_value Single NWP value at this step (NULL to omit)
#' @param current_period Absolute period index (1-indexed) being forecast
create_ml_pred_features_at_k <- function(current_data, k, m,
                                         nwp_value = NULL, current_period = 1) {
  n <- length(current_data)
  ppd <- m / k

  features <- data.frame(
    lag_1 = current_data[n],
    lag_1day = ifelse(n >= ppd, current_data[n - ppd + 1], current_data[max(1, n)]),
    lag_2days = ifelse(n >= 2 * ppd, current_data[n - 2 * ppd + 1],
                       ifelse(n >= ppd, current_data[n - ppd + 1], current_data[max(1, n)])),
    lag_1week = ifelse(n >= 7 * ppd, current_data[n - 7 * ppd + 1],
                       ifelse(n >= ppd, current_data[n - ppd + 1], current_data[max(1, n)]))
  )

  if (!is.null(nwp_value)) {
    features$nwp <- nwp_value
  }

  if (ppd >= 2 && n >= ppd) {
    features$rolling_mean <- mean(tail(current_data, ppd))
    features$rolling_sd <- sd(tail(current_data, ppd))
  } else {
    features$rolling_mean <- mean(current_data)
    features$rolling_sd <- ifelse(n > 1, sd(current_data), 0)
  }

  pod <- (current_period - 1) %% ppd
  day_idx <- (current_period - 1) %/% ppd + 1
  features$period_of_day <- pod
  features$day_of_week <- (day_idx - 1) %% 7
  features$month <- ((day_idx - 1) %/% 30) %% 12 + 1

  if (ppd > 1) {
    features$period_sin <- sin(2 * pi * pod / ppd)
    features$period_cos <- cos(2 * pi * pod / ppd)
  } else {
    features$period_sin <- 0
    features$period_cos <- 1
  }

  features
}

cat("Helper functions loaded successfully\n")
