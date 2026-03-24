# 3 - evaluate_forecasts.R
# ========================================
# Unified Forecast Evaluation
# Calculates nRMSE, nMBE, Skill, and SMAPE for all methods
# ========================================

rm(list = ls())

# Load required packages
libs <- c("data.table", "dplyr", "tidyr", "Matrix")
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
cat("Unified Forecast Evaluation\n")
cat("========================================\n\n")

# ----------------------------------------
# Load measurement data for actual values
# ----------------------------------------
# Load data and handle potential timestamp column
meas_raw <- fread(DATA_PATH)
# Check if first column is timestamp (character/POSIXct) or numeric data
if (is.character(meas_raw[[1]]) || inherits(meas_raw[[1]], "POSIXct")) {
  meas <- as.matrix(meas_raw[, -1])  # Remove timestamp column
} else {
  meas <- as.matrix(meas_raw)  # No timestamp, keep all columns
}
n_stations <- ncol(meas)
n_obs <- nrow(meas)

# Create aggregated actuals using S matrix
S <- as.matrix(hts_info$S)
n_upper <- nrow(S) - ncol(S)  # 6

meas_agg <- t(S[1:n_upper, , drop = FALSE] %*% t(meas))
meas_full <- cbind(meas_agg, meas)  # 324 columns

cat(sprintf("Loaded %d series (%d aggregated + %d bottom-level)\n",
            ncol(meas_full), n_upper, n_stations))

# ----------------------------------------
# Helper function to extract actual values for a replication
# ----------------------------------------
get_actual_values <- function(rp, k, all_meas) {
  obs_per_day <- m
  start_idx <- (rp - 1) * obs_per_day + 1
  test_start <- start_idx + train.days * obs_per_day
  test_end <- test_start + h * obs_per_day - 1

  if (test_end > nrow(all_meas)) {
    return(NULL)
  }

  # Get hourly actuals (use drop = FALSE to preserve matrix)
  hourly_actuals <- all_meas[test_start:test_end, , drop = FALSE]  # 48 x 324

  # Aggregate temporally if k > 1
  if (k > 1) {
    n_periods <- (h * m) / k  # Number of k-hour periods
    agg_actuals <- matrix(NA, nrow = n_periods, ncol = ncol(hourly_actuals))

    for (p in 1:n_periods) {
      idx_start <- (p - 1) * k + 1
      idx_end <- p * k
      agg_actuals[p, ] <- colSums(hourly_actuals[idx_start:idx_end, , drop = FALSE])
    }
    return(agg_actuals)
  }

  return(hourly_actuals)
}

# ----------------------------------------
# Helper function to extract forecasts at a specific k level
# ----------------------------------------
extract_forecasts_at_k <- function(reco_results, rp, k) {
  # reco_results is a list of replications
  # Each has $free (or $nn) as [120 x 324] matrix

  if (rp > length(reco_results) || is.null(reco_results[[rp]]$free)) {
    return(NULL)
  }

  fc_matrix <- reco_results[[rp]]$free

  # Extract rows corresponding to aggregation level k
  # The structure is: for each day, k levels are ordered 24, 12, 8, 6, 4, 3, 2, 1
  k_order <- c(24, 12, 8, 6, 4, 3, 2, 1)

  # Calculate row indices for k
  result <- c()

  for (day in 1:h) {
    # Starting row for this day
    day_start <- (day - 1) * sum(m / k_order) + 1

    # Find k in k_order
    k_idx <- which(k_order == k)
    if (length(k_idx) == 0) return(NULL)

    # Calculate offset to k
    if (k_idx > 1) {
      offset <- sum(m / k_order[1:(k_idx - 1)])
    } else {
      offset <- 0
    }

    # Number of periods at this k level
    n_periods <- m / k

    # Row indices for this day and k level
    row_start <- day_start + offset
    row_end <- row_start + n_periods - 1

    result <- c(result, row_start:row_end)
  }

  return(fc_matrix[result, , drop = FALSE])
}

# ----------------------------------------
# Calculate PERS (Persistence) baseline
# ----------------------------------------
cat("Calculating PERS baseline...\n")

pers_results <- list()

for (rp in rep_range) {
  obs_per_day <- m
  start_idx <- (rp - 1) * obs_per_day + 1
  train_end <- start_idx + train.days * obs_per_day - 1

  # Last h days of training = persistence forecast (matching original approach)
  # For rp=1: training is days 1-14, so persistence uses days 13-14 (last 2 days)
  # This uses ACTUAL measurements from both days, not repeating one day
  last_hdays_start <- train_end - h * m + 1

  # Check bounds
  if (last_hdays_start < 1 || train_end > nrow(meas_full)) {
    pers_results[[rp]] <- list(free = NULL, hourly = NULL)
    next
  }

  # Use drop = FALSE to preserve matrix dimensions
  pers_hourly <- meas_full[last_hdays_start:train_end, , drop = FALSE]  # h*m x 324 (48 x 324)

  # Store hourly persistence (no need to repeat - using actual 2-day pattern)
  pers_results[[rp]] <- list(free = NULL, hourly = pers_hourly)
}

# ----------------------------------------
# Methods to evaluate
# ----------------------------------------
methods <- list(
  list(name = "pers", type = "baseline", dir = NULL),
  list(name = "sarimax_base", type = "base", dir = dir_sarimax),
  list(name = "rf_base", type = "base", dir = dir_rf),
  list(name = "lgbm_base", type = "base", dir = dir_lgbm),
  list(name = "ets_base", type = "base", dir = dir_ets),
  list(name = "ctwlsv_sarimax", type = "reco", dir = dir_ctwlsv_sarimax),
  list(name = "ctwlsv_rf", type = "reco", dir = dir_ctwlsv_rf),
  list(name = "ctwlsv_lgbm", type = "reco", dir = dir_ctwlsv_lgbm),
  list(name = "ctwlsv_ets", type = "reco", dir = dir_ctwlsv_ets),
  list(name = "ctbu_sarimax", type = "reco", dir = dir_ctbu_sarimax),
  list(name = "ctbu_rf", type = "reco", dir = dir_ctbu_rf),
  list(name = "ctbu_lgbm", type = "reco", dir = dir_ctbu_lgbm),
  list(name = "ctbu_ets", type = "reco", dir = dir_ctbu_ets)
)

# ----------------------------------------
# POOLED RMSE CALCULATION (matching original approach)
# ----------------------------------------
# The original computes RMSE by pooling ALL replications together,
# not by averaging per-replication RMSEs. This is mathematically different.

# Helper function for pooled nRMSE
calc_pooled_nRMSE <- function(all_actuals, all_forecasts) {
  # all_actuals and all_forecasts are vectors of all values across replications
  mse <- mean((all_actuals - all_forecasts)^2, na.rm = TRUE)
  rmse <- sqrt(mse)
  mean_actual <- mean(all_actuals, na.rm = TRUE)
  if (mean_actual == 0) return(NA)
  return(rmse / mean_actual * 100)
}

# ----------------------------------------
# Initialize results tibble (for pooled RMSE approach)
# ----------------------------------------
pooled_results <- tibble()

# ----------------------------------------
# Evaluate each method using POOLED RMSE
# (matching the original approach where RMSE is computed
#  across ALL replications together, not averaged)
# ----------------------------------------
for (method_info in methods) {

  method_name <- method_info$name
  method_type <- method_info$type

  cat(sprintf("\nEvaluating: %s\n", method_name))

  # Load method results
  if (method_type == "baseline") {
    # Use pre-calculated PERS
    method_results <- pers_results

  } else if (method_type == "base") {
    # Load base forecasts
    dir_path <- method_info$dir
    files <- list.files(dir_path, pattern = "\\.RData$", full.names = TRUE)

    method_results <- list()
    for (rp in rep_range) {
      Yhat_full <- matrix(NA, nrow = 120, ncol = nrow(S))

      # IMPORTANT: Load files in the correct order to match S matrix structure
      # S matrix has: Total (col 1), Region_01 (col 2), ..., Region_05 (col 6), then bottom-level
      # Files must be loaded in this exact order!

      # Define the correct order for aggregated files
      agg_names <- c("Total", "Region_01", "Region_02", "Region_03", "Region_04", "Region_05")
      agg_files <- files[grepl("--0--", files)]
      bottom_files <- files[!grepl("--0--", files)]

      # Sort aggregated files in correct order
      agg_files_ordered <- character(length(agg_names))
      for (i in seq_along(agg_names)) {
        pattern <- paste0(agg_names[i], "--0--")
        matching <- agg_files[grepl(pattern, agg_files, fixed = TRUE)]
        if (length(matching) > 0) {
          agg_files_ordered[i] <- matching[1]
        }
      }
      agg_files_ordered <- agg_files_ordered[agg_files_ordered != ""]

      col_idx <- 1
      for (f in agg_files_ordered) {
        load(f)
        if (rp <= length(results) && !is.null(results[[rp]]$Y.hat)) {
          Yhat_full[, col_idx] <- results[[rp]]$Y.hat
        }
        col_idx <- col_idx + 1
      }
      for (f in bottom_files) {
        load(f)
        if (rp <= length(results) && !is.null(results[[rp]]$Y.hat)) {
          Yhat_full[, col_idx] <- results[[rp]]$Y.hat
        }
        col_idx <- col_idx + 1
      }

      method_results[[rp]] <- list(free = Yhat_full)
    }

  } else {
    # Load reconciled forecasts
    dir_path <- method_info$dir
    result_file <- list.files(dir_path, pattern = "\\.RData$", full.names = TRUE)[1]

    if (is.na(result_file) || !file.exists(result_file)) {
      cat(sprintf("  Warning: No results found for %s, skipping\n", method_name))
      next
    }

    load(result_file)  # loads 'results'
    method_results <- results
  }

  # For each k level, collect ALL data across replications and compute pooled metrics
  for (k in c(1, 24)) {  # Focus on hourly (k=1) and daily (k=24)

    freq_label <- ifelse(k == 1, "Hourly", "Daily")
    n_periods <- (h * m) / k  # 48 for k=1, 2 for k=24

    # Initialize storage for pooled data (per series)
    n_series <- ncol(meas_full)
    all_actuals <- vector("list", n_series)
    all_forecasts <- vector("list", n_series)
    for (s in 1:n_series) {
      all_actuals[[s]] <- numeric(0)
      all_forecasts[[s]] <- numeric(0)
    }

    # Collect data across all replications
    for (rp in rep_range) {

      # Get actual values
      actuals <- get_actual_values(rp, k, meas_full)
      if (is.null(actuals)) next

      # Ensure actuals is a matrix
      if (!is.matrix(actuals)) {
        actuals <- matrix(actuals, nrow = 1)
      }

      # Get forecasts
      if (method_type == "baseline") {
        # Check if rp exists in method_results
        if (rp > length(method_results) || is.null(method_results[[rp]]$hourly)) {
          next
        }
        pers_hourly <- method_results[[rp]]$hourly

        # Ensure pers_hourly is a matrix
        if (!is.matrix(pers_hourly)) {
          pers_hourly <- matrix(pers_hourly, ncol = n_series)
        }

        if (k > 1) {
          forecasts <- matrix(NA, nrow = n_periods, ncol = ncol(pers_hourly))
          for (p in 1:n_periods) {
            idx_start <- (p - 1) * k + 1
            idx_end <- p * k
            forecasts[p, ] <- colSums(pers_hourly[idx_start:idx_end, , drop = FALSE])
          }
        } else {
          forecasts <- pers_hourly
        }
      } else {
        forecasts <- extract_forecasts_at_k(method_results, rp, k)
        if (is.null(forecasts)) next
      }

      # Ensure forecasts is a matrix
      if (!is.matrix(forecasts)) {
        forecasts <- matrix(forecasts, nrow = 1)
      }

      # Validate dimensions match
      n_actual_cols <- ncol(actuals)
      n_forecast_cols <- ncol(forecasts)
      n_cols_to_use <- min(n_series, n_actual_cols, n_forecast_cols)

      # Append data for each series (only up to valid columns)
      for (s in 1:n_cols_to_use) {
        all_actuals[[s]] <- c(all_actuals[[s]], actuals[, s])
        all_forecasts[[s]] <- c(all_forecasts[[s]], forecasts[, s])
      }
    }

    # Compute pooled metrics for each series
    for (series_idx in 1:n_series) {

      actual_vec <- all_actuals[[series_idx]]
      forecast_vec <- all_forecasts[[series_idx]]

      # Skip if empty or has NA
      if (length(actual_vec) == 0 || any(is.na(actual_vec)) || any(is.na(forecast_vec))) next

      # Determine hierarchy level
      if (series_idx == 1) {
        level <- "L0"
      } else if (series_idx <= 6) {
        level <- "L1"
      } else {
        level <- "L2"
      }

      # Calculate POOLED metrics (all replications together)
      nrmse <- calc_pooled_nRMSE(actual_vec, forecast_vec)
      nmbe <- mean(forecast_vec - actual_vec, na.rm = TRUE) / mean(actual_vec, na.rm = TRUE) * 100
      smape <- calc_SMAPE(actual_vec, forecast_vec)

      # Add to results (one row per series per k level)
      result_row <- tibble(
        method = method_name,
        k = k,
        freq = freq_label,
        series = series_idx,
        level = level,
        nRMSE = nrmse,
        nMBE = nmbe,
        SMAPE = smape
      )

      pooled_results <- bind_rows(pooled_results, result_row)
    }
    cat(sprintf("  k=%d: processed %d series\n", k, n_series))
  }
}

cat("\n========================================\n")
cat("Calculating skill scores...\n")

# ----------------------------------------
# Calculate skill scores (vs PERS baseline)
# ----------------------------------------
pers_rmse <- pooled_results %>%
  filter(method == "pers") %>%
  select(k, freq, series, level, pers_nRMSE = nRMSE)

# Join and calculate skill
all_results <- pooled_results %>%
  left_join(pers_rmse, by = c("k", "freq", "series", "level")) %>%
  mutate(
    skill = calc_skill(nRMSE, pers_nRMSE)
  )

# ----------------------------------------
# Save full results
# ----------------------------------------
save(all_results, file = "output/evaluation_results.RData")
cat("Full results saved to output/evaluation_results.RData\n")

# ----------------------------------------
# Create summary statistics
# (Average nRMSE across series within each level)
# ----------------------------------------
cat("\nGenerating summary statistics...\n")

summary_stats <- all_results %>%
  group_by(method, level, freq) %>%
  summarise(
    mean_nRMSE = mean(nRMSE, na.rm = TRUE),
    median_nRMSE = median(nRMSE, na.rm = TRUE),
    sd_nRMSE = sd(nRMSE, na.rm = TRUE),
    min_nRMSE = min(nRMSE, na.rm = TRUE),
    max_nRMSE = max(nRMSE, na.rm = TRUE),
    mean_nMBE = mean(nMBE, na.rm = TRUE),
    median_nMBE = median(nMBE, na.rm = TRUE),
    sd_nMBE = sd(nMBE, na.rm = TRUE),
    mean_skill = mean(skill, na.rm = TRUE),
    median_skill = median(skill, na.rm = TRUE),
    mean_SMAPE = mean(SMAPE, na.rm = TRUE),
    median_SMAPE = median(SMAPE, na.rm = TRUE),
    n_series = n(),
    .groups = "drop"
  )

# Save summary
fwrite(summary_stats, "output/summary_all_methods.csv")
cat("Summary saved to output/summary_all_methods.csv\n")

# ----------------------------------------
# Print summary table
# ----------------------------------------
cat("\n========================================\n")
cat("SUMMARY: Mean nRMSE by Method and Level\n")
cat("========================================\n\n")

summary_wide <- summary_stats %>%
  select(method, level, freq, mean_nRMSE) %>%
  pivot_wider(names_from = c(level, freq), values_from = mean_nRMSE,
              names_sep = "_")

print(as.data.frame(summary_wide), row.names = FALSE)

cat("\n========================================\n")
cat("SUMMARY: Mean Skill Score by Method and Level\n")
cat("========================================\n\n")

skill_wide <- summary_stats %>%
  select(method, level, freq, mean_skill) %>%
  pivot_wider(names_from = c(level, freq), values_from = mean_skill,
              names_sep = "_")

print(as.data.frame(skill_wide), row.names = FALSE)

cat("\n========================================\n")
cat("Evaluation Complete\n")
cat("========================================\n")
cat(sprintf("Total observations evaluated: %d\n", nrow(all_results)))
cat(sprintf("Methods evaluated: %s\n", paste(unique(all_results$method), collapse = ", ")))
