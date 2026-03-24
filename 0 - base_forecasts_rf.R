# 0 - base_forecasts_rf.R
# ========================================
# Random Forest Base Forecasts with Full Feature Engineering
# Following PLAN.md: NWP + Lags + Rolling Stats + Calendar Features
# ========================================

rm(list = ls())

# Load required packages
libs <- c("ranger", "doParallel", "doSNOW", "progress", "data.table", "Matrix")
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
cat("Random Forest Base Forecast Generation\n")
cat("With NWP + Full Feature Engineering\n")
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

# Verify alignment
if (nrow(meas) != nrow(pred)) {
  stop("Measurement and prediction data have different number of rows!")
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
cat(sprintf("  Y.hat rows: %d\n", n_yhat_rows))
cat(sprintf("  Y.obs rows: %d\n", n_yobs_rows))
cat(sprintf("  Res.insamp rows: %d\n", n_res_rows))

# ----------------------------------------
# Feature creation function (FULL - per PLAN.md)
# ----------------------------------------
create_rf_features <- function(data, nwp_data, start_hour = 1) {
  # Creates features for Random Forest training/prediction
  #
  # Features (per PLAN.md):
  # - Lag features: lag_1h, lag_24h, lag_48h, lag_168h (1 week)
  # - Rolling statistics: rolling_mean_24h, rolling_sd_24h
  # - Calendar features: hour_of_day, day_of_week, month
  # - NWP predictions

  n <- length(data)

  # Initialize feature dataframe
  features <- data.frame(
    # NWP prediction (key feature)
    nwp = nwp_data,

    # Lag features
    lag_1 = c(NA, data[1:(n - 1)]),
    lag_24 = c(rep(NA, 24), data[1:(n - 24)]),
    lag_48 = c(rep(NA, 48), data[1:(n - 48)]),
    lag_168 = c(rep(NA, 168), data[1:(n - 168)])  # 1 week lag
  )

  # Rolling statistics (24-hour window)
  rolling_mean <- rep(NA, n)
  rolling_sd <- rep(NA, n)
  for (i in 24:n) {
    window <- data[(i - 23):i]
    rolling_mean[i] <- mean(window, na.rm = TRUE)
    rolling_sd[i] <- sd(window, na.rm = TRUE)
  }
  features$rolling_mean_24 <- rolling_mean
  features$rolling_sd_24 <- rolling_sd

  # Calendar features (assuming hourly data starting from a known point)
  # Hour of day: cycles 0-23
  hour_of_day <- ((start_hour - 1 + 0:(n - 1)) %% 24)
  features$hour_of_day <- hour_of_day

  # For day_of_week and month, we use cyclic features
  # Since we don't have actual dates, use proxies based on position
  day_idx <- ((start_hour - 1 + 0:(n - 1)) %/% 24) + 1
  features$day_of_week <- (day_idx - 1) %% 7  # 0-6

  # Month approximation (assuming 365 days/year, ~30 days/month)
  features$month <- ((day_idx - 1) %/% 30) %% 12 + 1  # 1-12

  # Sin/cos encoding for cyclic features (better for RF)
  features$hour_sin <- sin(2 * pi * hour_of_day / 24)
  features$hour_cos <- cos(2 * pi * hour_of_day / 24)

  # Target variable
  features$target <- data

  return(features)
}

# Feature creation for prediction (single step)
create_prediction_features <- function(current_data, nwp_value, current_hour, current_day) {
  n <- length(current_data)

  features <- data.frame(
    nwp = nwp_value,
    lag_1 = current_data[n],
    lag_24 = ifelse(n >= 24, current_data[n - 23], current_data[n]),
    lag_48 = ifelse(n >= 48, current_data[n - 47], current_data[max(1, n - 23)]),
    lag_168 = ifelse(n >= 168, current_data[n - 167], current_data[max(1, n - 23)]),
    rolling_mean_24 = ifelse(n >= 24, mean(tail(current_data, 24)), mean(current_data)),
    rolling_sd_24 = ifelse(n >= 24, sd(tail(current_data, 24)), sd(current_data)),
    hour_of_day = current_hour %% 24,
    day_of_week = (current_day - 1) %% 7,
    month = ((current_day - 1) %/% 30) %% 12 + 1,
    hour_sin = sin(2 * pi * (current_hour %% 24) / 24),
    hour_cos = cos(2 * pi * (current_hour %% 24) / 24)
  )

  return(features)
}

# ----------------------------------------
# Setup parallel processing
# ----------------------------------------
cat(sprintf("\nSetting up parallel processing with %d cores...\n", ncores))

cl <- makeCluster(ncores)
registerDoSNOW(cl)

clusterExport(cl, c("m", "h", "k.v", "train.days", "obs_per_day",
                    "create_temporal_hierarchy", "create_residual_hierarchy",
                    "create_rf_features", "create_prediction_features",
                    "n_yhat_rows", "n_yobs_rows", "n_res_rows"))

clusterEvalQ(cl, {
  library(ranger)
})

# ----------------------------------------
# Process each station (bottom-level)
# ----------------------------------------
cat(sprintf("\nProcessing %d bottom-level stations...\n", n_stations))

for (station_idx in 1:n_stations) {

  station_name <- colnames(meas)[station_idx]
  station_meas <- meas[, station_idx]
  station_pred <- pred[, station_idx]

  cat(sprintf("\n[Station %d/%d] %s\n", station_idx, n_stations, station_name))

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

      # Need extra history for lag and rolling features (1 week = 168 hours)
      extended_start <- max(1, start_idx - 168)
      extended_meas <- station_meas[extended_start:end_train]
      extended_pred <- station_pred[extended_start:end_train]

      # Create features with extended data
      # start_hour for calendar features
      start_hour_offset <- extended_start
      features <- create_rf_features(extended_meas, extended_pred, start_hour_offset)

      # Remove rows with NA (due to lags)
      features_clean <- features[complete.cases(features), ]

      if (nrow(features_clean) < 100) {
        return(list(Y.hat = rep(NA, n_yhat_rows),
                    Y.obs = rep(NA, n_yobs_rows),
                    Res.insamp = rep(NA, n_res_rows),
                    error = "Insufficient training data"))
      }

      # Train Random Forest
      rf_model <- ranger(
        target ~ . - target,  # Exclude target from predictors
        data = features_clean,
        num.trees = 500,
        mtry = floor(sqrt(ncol(features_clean) - 1)),
        min.node.size = 5,
        importance = "none",  # Skip for speed
        seed = 42 + rp
      )

      # Recursive multi-step forecasting with NWP
      fc_hourly <- numeric(h * m)
      train_data <- station_meas[start_idx:end_train]
      current_data <- train_data
      forecast_nwp <- station_pred[(end_train + 1):end_test]

      # Track hour and day for calendar features
      current_hour <- end_train + 1  # Hour index
      current_day <- (end_train %/% 24) + 1  # Day index

      for (step in 1:(h * m)) {
        # Create prediction features
        new_features <- create_prediction_features(
          current_data,
          forecast_nwp[step],
          current_hour,
          current_day
        )

        # Predict
        pred_val <- predict(rf_model, new_features)$predictions
        pred_val <- max(pred_val, 0)  # Non-negative

        fc_hourly[step] <- pred_val

        # Update for next step
        current_data <- c(current_data, pred_val)
        current_hour <- current_hour + 1
        if ((current_hour - 1) %% 24 == 0) {
          current_day <- current_day + 1
        }
      }

      # In-sample residuals
      train_features <- features_clean[, names(features_clean) != "target"]
      train_pred <- predict(rf_model, train_features)$predictions
      residuals_hourly <- features_clean$target - train_pred

      # Pad/trim residuals
      n_expected_res <- train.days * m
      if (length(residuals_hourly) < n_expected_res) {
        residuals_hourly <- c(rep(0, n_expected_res - length(residuals_hourly)), residuals_hourly)
      } else if (length(residuals_hourly) > n_expected_res) {
        residuals_hourly <- tail(residuals_hourly, n_expected_res)
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

  filename <- sprintf("%s/%s--rf.RData", dir_rf, station_name)
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

# Aggregate both measurement and prediction data
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

      extended_start <- max(1, start_idx - 168)
      extended_meas <- agg_meas[extended_start:end_train]
      extended_pred <- agg_pred[extended_start:end_train]

      features <- create_rf_features(extended_meas, extended_pred, extended_start)
      features_clean <- features[complete.cases(features), ]

      if (nrow(features_clean) < 100) {
        return(list(Y.hat = rep(NA, n_yhat_rows),
                    Y.obs = rep(NA, n_yobs_rows),
                    Res.insamp = rep(NA, n_res_rows),
                    error = "Insufficient training data"))
      }

      rf_model <- ranger(
        target ~ . - target,
        data = features_clean,
        num.trees = 500,
        mtry = floor(sqrt(ncol(features_clean) - 1)),
        min.node.size = 5,
        importance = "none",
        seed = 42 + rp
      )

      fc_hourly <- numeric(h * m)
      train_data <- agg_meas[start_idx:end_train]
      current_data <- train_data
      forecast_nwp <- agg_pred[(end_train + 1):end_test]

      current_hour <- end_train + 1
      current_day <- (end_train %/% 24) + 1

      for (step in 1:(h * m)) {
        new_features <- create_prediction_features(
          current_data,
          forecast_nwp[step],
          current_hour,
          current_day
        )

        pred_val <- max(predict(rf_model, new_features)$predictions, 0)
        fc_hourly[step] <- pred_val
        current_data <- c(current_data, pred_val)
        current_hour <- current_hour + 1
        if ((current_hour - 1) %% 24 == 0) current_day <- current_day + 1
      }

      train_features <- features_clean[, names(features_clean) != "target"]
      train_pred <- predict(rf_model, train_features)$predictions
      residuals_hourly <- features_clean$target - train_pred

      n_expected_res <- train.days * m
      if (length(residuals_hourly) < n_expected_res) {
        residuals_hourly <- c(rep(0, n_expected_res - length(residuals_hourly)), residuals_hourly)
      } else if (length(residuals_hourly) > n_expected_res) {
        residuals_hourly <- tail(residuals_hourly, n_expected_res)
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

  filename <- sprintf("%s/%s--0--rf.RData", dir_rf, agg_name)
  save(results, file = filename)
  cat(sprintf("  Saved to: %s\n", filename))
}

# ----------------------------------------
# Cleanup
# ----------------------------------------
stopCluster(cl)

cat("\n========================================\n")
cat("Random Forest Base Forecasts Complete\n")
cat("========================================\n")
cat(sprintf("Total stations processed: %d\n", n_stations))
cat(sprintf("Aggregated series processed: %d\n", n_upper))
cat(sprintf("Replications per series: %d\n", n_rep))
cat(sprintf("Results saved to: %s\n", dir_rf))
cat("\nFeatures used:\n")
cat("  - NWP predictions: YES\n")
cat("  - Lag features: lag_1h, lag_24h, lag_48h, lag_168h\n")
cat("  - Rolling stats: rolling_mean_24h, rolling_sd_24h\n")
cat("  - Calendar features: hour_of_day, day_of_week, month (+ sin/cos)\n")
