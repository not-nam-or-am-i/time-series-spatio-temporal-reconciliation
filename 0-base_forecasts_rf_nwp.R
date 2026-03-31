# 0-base_forecasts_rf_nwp.R
# ========================================
# Random Forest Base Forecasts with L2 NWP Replacement
#
# Training features for ALL levels (L0, L1, L2):
#   - lag_1, lag_24, lag_48, lag_168
#   - rolling_mean_24, rolling_sd_24
#   - calendar features (hour/day/month + sin/cos)
# (NO NWP as an RF input feature)
#
# Forecast generation:
#   - L2: train RF on lag/rolling/calendar only, then
#        (a) generate RF recursive hourly forecasts for all 48 hours
#        (b) overwrite ONLY the k=1 (hourly) rows of Y.hat with NWP
#        (c) overwrite ONLY the k=1 segment of Res.insamp with NWP residuals
#   - All negatives are clipped to zero.
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
cat("L2: NWP hourly replacement (k=1 only)\n")
cat("L0/L1/L2 training: lag/rolling/calendar only (no NWP feature)\n")
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
# Load NWP predictions (used only for L2 replacement)
# ----------------------------------------
pred_raw <- fread(PRED_PATH)
if (is.character(pred_raw[[1]]) || inherits(pred_raw[[1]], "POSIXct")) {
  pred <- as.matrix(pred_raw[, -1])
} else {
  pred <- as.matrix(pred_raw)
}

cat(sprintf("NWP predictions loaded: %d observations x %d stations\n", nrow(pred), ncol(pred)))

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
# Feature creation function (NO NWP feature)
# ----------------------------------------
create_rf_features <- function(data, start_hour = 1) {
  # Creates features for Random Forest training.
  #
  # Features (per baseline RF plan, minus NWP):
  # - Lag features: lag_1, lag_24, lag_48, lag_168
  # - Rolling statistics: rolling_mean_24, rolling_sd_24
  # - Calendar features: hour_of_day, day_of_week, month (+ sin/cos)

  n <- length(data)

  features <- data.frame(
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

  # Calendar features
  hour_of_day <- ((start_hour - 1 + 0:(n - 1)) %% 24)
  features$hour_of_day <- hour_of_day

  # Proxies for day_of_week and month based on position
  day_idx <- ((start_hour - 1 + 0:(n - 1)) %/% 24) + 1
  features$day_of_week <- (day_idx - 1) %% 7
  features$month <- ((day_idx - 1) %/% 30) %% 12 + 1

  features$hour_sin <- sin(2 * pi * hour_of_day / 24)
  features$hour_cos <- cos(2 * pi * hour_of_day / 24)

  # Target variable
  features$target <- data

  features
}

# ----------------------------------------
# Prediction feature creation (single step, NO NWP feature)
# ----------------------------------------
create_prediction_features <- function(current_data, current_hour, current_day) {
  n <- length(current_data)

  features <- data.frame(
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

  features
}

# ----------------------------------------
# Setup parallel processing
# ----------------------------------------
cat(sprintf("\nSetting up parallel processing with %d cores...\n", ncores))

cl <- makeCluster(ncores)
# for SLURM cluster
# cl <- parallel::makeCluster(ncores, type = "PSOCK")
registerDoSNOW(cl)

clusterExport(
  cl,
  c(
    "m", "h", "k.v", "train.days", "obs_per_day",
    "create_temporal_hierarchy", "create_residual_hierarchy",
    "create_rf_features", "create_prediction_features",
    "n_yhat_rows", "n_yobs_rows", "n_res_rows",
    "k_star"
  )
)

clusterEvalQ(cl, {
  library(ranger)
})

# ----------------------------------------
# Process each station (bottom-level, L2)
# ----------------------------------------
cat(sprintf("\nProcessing %d bottom-level stations (RF training only + L2 NWP replacement)...\n", n_stations))

for (station_idx in 1:n_stations) {
  station_name <- colnames(meas)[station_idx]
  station_meas <- meas[, station_idx]
  station_pred <- pred[, station_idx]

  cat(sprintf("\n[Station %d/%d] %s\n", station_idx, n_stations, station_name))

  pb <- txtProgressBar(max = n_rep, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  clusterExport(cl, c("station_meas", "station_pred"), envir = environment())

  results <- foreach(rp = rep_range, .options.snow = opts, .errorhandling = "pass") %dopar% {
    tryCatch({
      start_idx <- (rp - 1) * obs_per_day + 1
      end_train <- start_idx + train.days * obs_per_day - 1
      end_test <- end_train + h * obs_per_day

      if (end_test > length(station_meas)) {
        return(list(
          Y.hat = rep(NA, n_yhat_rows),
          Y.obs = rep(NA, n_yobs_rows),
          Res.insamp = rep(NA, n_res_rows),
          error = "Index out of bounds"
        ))
      }

      # Extract training period data (RF features train on lag/rolling/calendar only)
      train_y <- station_meas[start_idx:end_train]
      train_x <- station_pred[start_idx:end_train] # only for residual replacement at k=1
      forecast_nwp <- station_pred[(end_train + 1):end_test] # only for Y.hat replacement at k=1

      # Need extra history for lag and rolling features (1 week = 168 hours)
      extended_start <- max(1, start_idx - 168)
      extended_meas <- station_meas[extended_start:end_train]

      # Feature creation with aligned calendar features
      start_hour_offset <- extended_start
      features <- create_rf_features(extended_meas, start_hour_offset)

      # Remove rows with NA (due to lags)
      features_clean <- features[complete.cases(features), ]

      if (nrow(features_clean) < 100) {
        return(list(
          Y.hat = rep(NA, n_yhat_rows),
          Y.obs = rep(NA, n_yobs_rows),
          Res.insamp = rep(NA, n_res_rows),
          error = "Insufficient training data"
        ))
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

      # Recursive multi-step forecasting (RF only)
      fc_hourly_rf <- numeric(h * m)
      current_data <- train_y

      current_hour <- end_train + 1
      current_day <- (end_train %/% 24) + 1

      for (step in 1:(h * m)) {
        new_features <- create_prediction_features(
          current_data,
          current_hour,
          current_day
        )

        pred_val <- predict(rf_model, new_features)$predictions
        pred_val <- max(pred_val, 0) # non-negative at generation time

        fc_hourly_rf[step] <- pred_val

        current_data <- c(current_data, pred_val)
        current_hour <- current_hour + 1
        if ((current_hour - 1) %% 24 == 0) current_day <- current_day + 1
      }

      # In-sample residuals (RF residuals at all k)
      train_features <- features_clean[, names(features_clean) != "target"]
      train_pred <- predict(rf_model, train_features)$predictions
      residuals_hourly_rf <- features_clean$target - train_pred

      # Pad/trim residuals to full train window length (train.days * m)
      n_expected_res <- train.days * m
      if (length(residuals_hourly_rf) < n_expected_res) {
        residuals_hourly_rf <- c(
          rep(0, n_expected_res - length(residuals_hourly_rf)),
          residuals_hourly_rf
        )
      } else if (length(residuals_hourly_rf) > n_expected_res) {
        residuals_hourly_rf <- tail(residuals_hourly_rf, n_expected_res)
      }

      # Build temporal hierarchy from RF hourly forecasts
      Y.hat <- create_temporal_hierarchy(fc_hourly_rf, k.v, m, h)
      Y.obs <- station_meas[(end_train + 1):end_test]

      # ----------------------------------------
      # L2 replacement: overwrite ONLY k=1 rows
      # ----------------------------------------
      nwp_hourly <- pmax(forecast_nwp, 0)

      hourly_offset <- k_star - m # start offset within each day block for k=1 rows
      for (day in 1:h) {
        y_day_offset <- (day - 1) * k_star
        row_start <- y_day_offset + hourly_offset + 1
        row_end <- y_day_offset + k_star
        Y.hat[row_start:row_end] <- nwp_hourly[(day - 1) * m + 1:day * m]
      }
      Y.hat <- pmax(Y.hat, 0)

      # Residual hierarchy from RF residuals (then overwrite k=1 segment)
      Res.insamp <- create_residual_hierarchy(residuals_hourly_rf, k.v, m, train.days)

      # NWP residuals at hourly level over the training window
      nwp_residuals_hourly <- train_x - train_y

      for (day in 1:train.days) {
        r_day_offset <- (day - 1) * k_star
        r_row_start <- r_day_offset + hourly_offset + 1
        r_row_end <- r_day_offset + k_star
        Res.insamp[r_row_start:r_row_end] <- nwp_residuals_hourly[(day - 1) * m + 1:day * m]
      }

      list(Y.hat = Y.hat, Y.obs = Y.obs, Res.insamp = Res.insamp, error = NULL)
    }, error = function(e) {
      list(
        Y.hat = rep(NA, n_yhat_rows),
        Y.obs = rep(NA, n_yobs_rows),
        Res.insamp = rep(NA, n_res_rows),
        error = as.character(e)
      )
    })
  }

  close(pb)

  errors <- sapply(results, function(x) !is.null(x$error))
  if (any(errors)) {
    cat(sprintf("  Warning: %d replications had errors\n", sum(errors)))
  }

  filename <- sprintf("%s/%s--rf_nwp.RData", dir_rf_nwp, station_name)
  save(results, file = filename)
  cat(sprintf("  Saved to: %s\n", filename))
}

# ----------------------------------------
# Process aggregated series (L0 + L1)
# RF forecasts trained without NWP feature (no L2 replacement here)
# ----------------------------------------
cat("\n========================================\n")
cat("Processing Aggregated Series (L0 + L1)\n")
cat("Random Forest training: lag/rolling/calendar only (no NWP feature)\n")
cat("========================================\n")

S <- as.matrix(hts_info$S)
n_upper <- nrow(S) - ncol(S)

meas_agg <- t(S[1:n_upper, , drop = FALSE] %*% t(meas))
pred_agg <- t(S[1:n_upper, , drop = FALSE] %*% t(pred)) # only needed to match baseline structure; not used for training features

agg_names <- c("Total", paste0("Region_", sprintf("%02d", 1:5)))

for (agg_idx in 1:n_upper) {
  agg_name <- agg_names[agg_idx]
  agg_meas <- meas_agg[, agg_idx]

  cat(sprintf("\n[Aggregate %d/%d] %s\n", agg_idx, n_upper, agg_name))

  pb <- txtProgressBar(max = n_rep, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  clusterExport(cl, c("agg_meas"), envir = environment())

  results <- foreach(rp = rep_range, .options.snow = opts, .errorhandling = "pass") %dopar% {
    tryCatch({
      start_idx <- (rp - 1) * obs_per_day + 1
      end_train <- start_idx + train.days * obs_per_day - 1
      end_test <- end_train + h * obs_per_day

      if (end_test > length(agg_meas)) {
        return(list(
          Y.hat = rep(NA, n_yhat_rows),
          Y.obs = rep(NA, n_yobs_rows),
          Res.insamp = rep(NA, n_res_rows),
          error = "Index out of bounds"
        ))
      }

      extended_start <- max(1, start_idx - 168)
      extended_meas <- agg_meas[extended_start:end_train]

      features <- create_rf_features(extended_meas, extended_start)
      features_clean <- features[complete.cases(features), ]

      if (nrow(features_clean) < 100) {
        return(list(
          Y.hat = rep(NA, n_yhat_rows),
          Y.obs = rep(NA, n_yobs_rows),
          Res.insamp = rep(NA, n_res_rows),
          error = "Insufficient training data"
        ))
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

      # Recursive multi-step forecasting (RF only)
      fc_hourly_rf <- numeric(h * m)
      train_y <- agg_meas[start_idx:end_train]
      current_data <- train_y

      current_hour <- end_train + 1
      current_day <- (end_train %/% 24) + 1

      for (step in 1:(h * m)) {
        new_features <- create_prediction_features(
          current_data,
          current_hour,
          current_day
        )

        pred_val <- predict(rf_model, new_features)$predictions
        pred_val <- max(pred_val, 0)
        fc_hourly_rf[step] <- pred_val

        current_data <- c(current_data, pred_val)
        current_hour <- current_hour + 1
        if ((current_hour - 1) %% 24 == 0) current_day <- current_day + 1
      }

      train_features <- features_clean[, names(features_clean) != "target"]
      train_pred <- predict(rf_model, train_features)$predictions
      residuals_hourly_rf <- features_clean$target - train_pred

      n_expected_res <- train.days * m
      if (length(residuals_hourly_rf) < n_expected_res) {
        residuals_hourly_rf <- c(
          rep(0, n_expected_res - length(residuals_hourly_rf)),
          residuals_hourly_rf
        )
      } else if (length(residuals_hourly_rf) > n_expected_res) {
        residuals_hourly_rf <- tail(residuals_hourly_rf, n_expected_res)
      }

      Y.hat <- create_temporal_hierarchy(fc_hourly_rf, k.v, m, h)
      Y.hat <- pmax(Y.hat, 0)
      Y.obs <- agg_meas[(end_train + 1):end_test]
      Res.insamp <- create_residual_hierarchy(residuals_hourly_rf, k.v, m, train.days)

      list(Y.hat = Y.hat, Y.obs = Y.obs, Res.insamp = Res.insamp, error = NULL)
    }, error = function(e) {
      list(
        Y.hat = rep(NA, n_yhat_rows),
        Y.obs = rep(NA, n_yobs_rows),
        Res.insamp = rep(NA, n_res_rows),
        error = as.character(e)
      )
    })
  }

  close(pb)

  errors <- sapply(results, function(x) !is.null(x$error))
  if (any(errors)) {
    cat(sprintf("  Warning: %d replications had errors\n", sum(errors)))
  }

  filename <- sprintf("%s/%s--0--rf_nwp.RData", dir_rf_nwp, agg_name)
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
cat(sprintf("Results saved to: %s\n", dir_rf_nwp))

cat("\nFeatures used:\n")
cat("  - NWP predictions as RF feature: NO (training excludes NWP)\n")
cat("  - Lag features: lag_1, lag_24, lag_48, lag_168\n")
cat("  - Rolling stats: rolling_mean_24h, rolling_sd_24h\n")
cat("  - Calendar features: hour_of_day, day_of_week, month (+ sin/cos)\n")

cat("\nReplacement policy:\n")
cat("  - L2 (k=1 hourly rows in Y.hat): replaced by NWP hourly (clipped at 0)\n")
cat("  - L2 (k=1 hourly segment in Res.insamp): replaced by NWP residuals\n")

