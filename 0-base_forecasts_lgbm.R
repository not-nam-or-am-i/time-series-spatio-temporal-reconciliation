# 0-base_forecasts_lgbm.R
# ========================================
# LightGBM Base Forecasts with Full Feature Engineering
# Same features as Random Forest: NWP + Lags + Rolling Stats + Calendar
# ========================================

rm(list = ls())

# Load required packages
libs <- c("lightgbm", "doParallel", "doSNOW", "progress", "data.table", "Matrix")
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
cat("LightGBM Base Forecast Generation\n")
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
# LightGBM hyperparameters
# ----------------------------------------
lgbm_params <- list(
  objective = "regression",
  metric = "rmse",
  learning_rate = 0.05,
  num_leaves = 31,
  min_data_in_leaf = 20,
  feature_fraction = 0.8,
  bagging_fraction = 0.8,
  bagging_freq = 5,
  verbose = -1
)
lgbm_nrounds <- 300

cat(sprintf("\nLightGBM parameters:\n"))
cat(sprintf("  num_leaves: %d\n", lgbm_params$num_leaves))
cat(sprintf("  learning_rate: %.2f\n", lgbm_params$learning_rate))
cat(sprintf("  nrounds: %d\n", lgbm_nrounds))
cat(sprintf("  min_data_in_leaf: %d\n", lgbm_params$min_data_in_leaf))

# ----------------------------------------
# Setup parallel processing
# ----------------------------------------
cat(sprintf("\nSetting up parallel processing with %d cores...\n", ncores))

# for local machine
# cl <- makeCluster(ncores)
# for SLURM cluster
cl <- parallel::makeCluster(ncores, type = "PSOCK")
registerDoSNOW(cl)

clusterExport(cl, c("m", "h", "k.v", "train.days", "obs_per_day",
                    "aggregate_hourly_to_k",
                    "create_ml_features_at_k", "create_ml_pred_features_at_k",
                    "assemble_yhat_from_levels", "assemble_residuals_from_levels",
                    "lgbm_params", "lgbm_nrounds",
                    "n_yhat_rows", "n_yobs_rows", "n_res_rows"))

clusterEvalQ(cl, {
  library(lightgbm)
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

      extended_start <- max(1, start_idx - 168)
      extended_meas_hourly <- station_meas[extended_start:end_train]
      extended_nwp_hourly <- station_pred[extended_start:end_train]
      train_y_hourly <- station_meas[start_idx:end_train]
      test_nwp_hourly <- station_pred[(end_train + 1):end_test]

      k_order <- sort(k.v, decreasing = TRUE)
      fc_list <- list()
      res_list <- list()

      for (kk in k_order) {
        ppd <- m / kk
        n_fc <- h * ppd
        n_train_k <- train.days * ppd

        extended_k <- aggregate_hourly_to_k(extended_meas_hourly, kk)
        extended_nwp_k <- aggregate_hourly_to_k(extended_nwp_hourly, kk)
        train_k <- aggregate_hourly_to_k(train_y_hourly, kk)
        test_nwp_k <- aggregate_hourly_to_k(test_nwp_hourly, kk)

        ext_start_period <- (extended_start - 1) %/% kk + 1
        features <- create_ml_features_at_k(extended_k, kk, m,
                                            nwp_data = extended_nwp_k,
                                            start_period = ext_start_period)
        features_clean <- features[complete.cases(features), ]

        if (nrow(features_clean) < 5) {
          fc_list[[as.character(kk)]] <- rep(0, n_fc)
          res_list[[as.character(kk)]] <- rep(0, n_train_k)
          next
        }

        feature_names <- setdiff(names(features_clean), "target")

        params_k <- lgbm_params
        params_k$min_data_in_leaf <- min(lgbm_params$min_data_in_leaf,
                                         max(1, nrow(features_clean) %/% 3))

        dtrain <- lgb.Dataset(
          data = as.matrix(features_clean[, feature_names]),
          label = features_clean$target
        )

        lgb_model <- lgb.train(
          params = params_k,
          data = dtrain,
          nrounds = lgbm_nrounds,
          verbose = -1
        )

        fc_k <- numeric(n_fc)
        current_data <- train_k
        fc_start_period <- end_train %/% kk + 1

        for (step in 1:n_fc) {
          new_features <- create_ml_pred_features_at_k(
            current_data, kk, m,
            nwp_value = test_nwp_k[step],
            current_period = fc_start_period + step - 1
          )
          pred_val <- max(predict(lgb_model, as.matrix(new_features[, feature_names])), 0)
          fc_k[step] <- pred_val
          current_data <- c(current_data, pred_val)
        }

        train_pred <- predict(lgb_model, as.matrix(features_clean[, feature_names]))
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

  filename <- sprintf("%s/%s--lgbm.RData", dir_lgbm, station_name)
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

      extended_start <- max(1, start_idx - 168)
      extended_meas_hourly <- agg_meas[extended_start:end_train]
      extended_nwp_hourly <- agg_pred[extended_start:end_train]
      train_y_hourly <- agg_meas[start_idx:end_train]
      test_nwp_hourly <- agg_pred[(end_train + 1):end_test]

      k_order <- sort(k.v, decreasing = TRUE)
      fc_list <- list()
      res_list <- list()

      for (kk in k_order) {
        ppd <- m / kk
        n_fc <- h * ppd
        n_train_k <- train.days * ppd

        extended_k <- aggregate_hourly_to_k(extended_meas_hourly, kk)
        extended_nwp_k <- aggregate_hourly_to_k(extended_nwp_hourly, kk)
        train_k <- aggregate_hourly_to_k(train_y_hourly, kk)
        test_nwp_k <- aggregate_hourly_to_k(test_nwp_hourly, kk)

        ext_start_period <- (extended_start - 1) %/% kk + 1
        features <- create_ml_features_at_k(extended_k, kk, m,
                                            nwp_data = extended_nwp_k,
                                            start_period = ext_start_period)
        features_clean <- features[complete.cases(features), ]

        if (nrow(features_clean) < 5) {
          fc_list[[as.character(kk)]] <- rep(0, n_fc)
          res_list[[as.character(kk)]] <- rep(0, n_train_k)
          next
        }

        feature_names <- setdiff(names(features_clean), "target")

        params_k <- lgbm_params
        params_k$min_data_in_leaf <- min(lgbm_params$min_data_in_leaf,
                                         max(1, nrow(features_clean) %/% 3))

        dtrain <- lgb.Dataset(
          data = as.matrix(features_clean[, feature_names]),
          label = features_clean$target
        )

        lgb_model <- lgb.train(
          params = params_k,
          data = dtrain,
          nrounds = lgbm_nrounds,
          verbose = -1
        )

        fc_k <- numeric(n_fc)
        current_data <- train_k
        fc_start_period <- end_train %/% kk + 1

        for (step in 1:n_fc) {
          new_features <- create_ml_pred_features_at_k(
            current_data, kk, m,
            nwp_value = test_nwp_k[step],
            current_period = fc_start_period + step - 1
          )
          pred_val <- max(predict(lgb_model, as.matrix(new_features[, feature_names])), 0)
          fc_k[step] <- pred_val
          current_data <- c(current_data, pred_val)
        }

        train_pred <- predict(lgb_model, as.matrix(features_clean[, feature_names]))
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

  filename <- sprintf("%s/%s--0--lgbm.RData", dir_lgbm, agg_name)
  save(results, file = filename)
  cat(sprintf("  Saved to: %s\n", filename))
}

# ----------------------------------------
# Cleanup
# ----------------------------------------
stopCluster(cl)

cat("\n========================================\n")
cat("LightGBM Base Forecasts Complete\n")
cat("========================================\n")
cat(sprintf("Total stations processed: %d\n", n_stations))
cat(sprintf("Aggregated series processed: %d\n", n_upper))
cat(sprintf("Replications per series: %d\n", n_rep))
cat(sprintf("Results saved to: %s\n", dir_lgbm))
cat("\nApproach:\n")
cat("  - Independent LightGBM model at each temporal aggregation level k\n")
cat("  - Features: NWP, lag (1-period, 1-day, 2-day, 1-week), rolling stats, calendar\n")
cat("  - Recursive multi-step forecasting at each k level\n")
cat(sprintf("  - LightGBM: %d rounds, lr=%.2f, %d leaves\n",
            lgbm_nrounds, lgbm_params$learning_rate, lgbm_params$num_leaves))
