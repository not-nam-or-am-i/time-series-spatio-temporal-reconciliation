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
# Setup parallel processing
# ----------------------------------------
cat(sprintf("\nSetting up parallel processing with %d cores...\n", ncores))

# for local machine
cl <- makeCluster(ncores)
# for SLURM cluster
# cl <- parallel::makeCluster(ncores, type = "PSOCK")
registerDoSNOW(cl)

clusterExport(
  cl,
  c(
    "m", "h", "k.v", "train.days", "obs_per_day",
    "aggregate_hourly_to_k",
    "create_ml_features_at_k", "create_ml_pred_features_at_k",
    "assemble_yhat_from_levels", "assemble_residuals_from_levels",
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
        return(list(Y.hat = rep(NA, n_yhat_rows),
                    Y.obs = rep(NA, n_yobs_rows),
                    Res.insamp = rep(NA, n_res_rows),
                    error = "Index out of bounds"))
      }

      extended_start <- max(1, start_idx - 168)
      extended_meas_hourly <- station_meas[extended_start:end_train]
      train_y_hourly <- station_meas[start_idx:end_train]
      train_x_hourly <- station_pred[start_idx:end_train]
      test_nwp_hourly <- station_pred[(end_train + 1):end_test]

      k_order <- sort(k.v, decreasing = TRUE)
      fc_list <- list()
      res_list <- list()

      for (kk in k_order) {
        ppd <- m / kk
        n_fc <- h * ppd
        n_train_k <- train.days * ppd

        if (kk == 1) {
          # k=1: NWP replacement (matching ETS hybrid approach)
          fc_list[["1"]] <- pmax(test_nwp_hourly, 0)
          res_list[["1"]] <- train_x_hourly - train_y_hourly
        } else {
          # k>1: independent RF (no NWP feature)
          extended_k <- aggregate_hourly_to_k(extended_meas_hourly, kk)
          train_k <- aggregate_hourly_to_k(train_y_hourly, kk)

          ext_start_period <- (extended_start - 1) %/% kk + 1
          features <- create_ml_features_at_k(extended_k, kk, m,
                                              nwp_data = NULL,
                                              start_period = ext_start_period)
          features_clean <- features[complete.cases(features), ]

          if (nrow(features_clean) < 5) {
            fc_list[[as.character(kk)]] <- rep(0, n_fc)
            res_list[[as.character(kk)]] <- rep(0, n_train_k)
            next
          }

          feature_names <- setdiff(names(features_clean), "target")
          rf_model <- ranger(
            target ~ . - target,
            data = features_clean,
            num.trees = 500,
            mtry = min(floor(sqrt(length(feature_names))), length(feature_names)),
            min.node.size = min(5, max(1, nrow(features_clean) - 1)),
            importance = "none",
            seed = 42 + rp
          )

          fc_k <- numeric(n_fc)
          current_data <- train_k
          fc_start_period <- end_train %/% kk + 1

          for (step in 1:n_fc) {
            new_features <- create_ml_pred_features_at_k(
              current_data, kk, m,
              nwp_value = NULL,
              current_period = fc_start_period + step - 1
            )
            pred_val <- max(predict(rf_model, new_features)$predictions, 0)
            fc_k[step] <- pred_val
            current_data <- c(current_data, pred_val)
          }

          train_pred <- predict(rf_model, features_clean[, feature_names, drop = FALSE])$predictions
          residuals_k <- features_clean$target - train_pred

          if (length(residuals_k) < n_train_k) {
            residuals_k <- c(rep(0, n_train_k - length(residuals_k)), residuals_k)
          } else if (length(residuals_k) > n_train_k) {
            residuals_k <- tail(residuals_k, n_train_k)
          }

          fc_list[[as.character(kk)]] <- fc_k
          res_list[[as.character(kk)]] <- residuals_k
        }
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
        return(list(Y.hat = rep(NA, n_yhat_rows),
                    Y.obs = rep(NA, n_yobs_rows),
                    Res.insamp = rep(NA, n_res_rows),
                    error = "Index out of bounds"))
      }

      extended_start <- max(1, start_idx - 168)
      extended_meas_hourly <- agg_meas[extended_start:end_train]
      train_y_hourly <- agg_meas[start_idx:end_train]

      k_order <- sort(k.v, decreasing = TRUE)
      fc_list <- list()
      res_list <- list()

      for (kk in k_order) {
        ppd <- m / kk
        n_fc <- h * ppd
        n_train_k <- train.days * ppd

        extended_k <- aggregate_hourly_to_k(extended_meas_hourly, kk)
        train_k <- aggregate_hourly_to_k(train_y_hourly, kk)

        ext_start_period <- (extended_start - 1) %/% kk + 1
        features <- create_ml_features_at_k(extended_k, kk, m,
                                            nwp_data = NULL,
                                            start_period = ext_start_period)
        features_clean <- features[complete.cases(features), ]

        if (nrow(features_clean) < 5) {
          fc_list[[as.character(kk)]] <- rep(0, n_fc)
          res_list[[as.character(kk)]] <- rep(0, n_train_k)
          next
        }

        feature_names <- setdiff(names(features_clean), "target")
        rf_model <- ranger(
          target ~ . - target,
          data = features_clean,
          num.trees = 500,
          mtry = min(floor(sqrt(length(feature_names))), length(feature_names)),
          min.node.size = min(5, max(1, nrow(features_clean) - 1)),
          importance = "none",
          seed = 42 + rp
        )

        fc_k <- numeric(n_fc)
        current_data <- train_k
        fc_start_period <- end_train %/% kk + 1

        for (step in 1:n_fc) {
          new_features <- create_ml_pred_features_at_k(
            current_data, kk, m,
            nwp_value = NULL,
            current_period = fc_start_period + step - 1
          )
          pred_val <- max(predict(rf_model, new_features)$predictions, 0)
          fc_k[step] <- pred_val
          current_data <- c(current_data, pred_val)
        }

        train_pred <- predict(rf_model, features_clean[, feature_names, drop = FALSE])$predictions
        residuals_k <- features_clean$target - train_pred

        if (length(residuals_k) < n_train_k) {
          residuals_k <- c(rep(0, n_train_k - length(residuals_k)), residuals_k)
        } else if (length(residuals_k) > n_train_k) {
          residuals_k <- tail(residuals_k, n_train_k)
        }

        fc_list[[as.character(kk)]] <- fc_k
        res_list[[as.character(kk)]] <- residuals_k
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

cat("\nApproach:\n")
cat("  - Independent RF model at each temporal aggregation level k\n")
cat("  - Features: lag (1-period, 1-day, 2-day, 1-week), rolling stats, calendar (no NWP)\n")
cat("  - L2 k=1: NWP hourly replacement (forecasts + residuals)\n")
cat("  - L2 k>1: independent RF at that temporal level\n")
cat("  - L0+L1: independent RF at each k level (no NWP)\n")

