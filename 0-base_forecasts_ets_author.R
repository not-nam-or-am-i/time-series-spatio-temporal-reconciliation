# 0-base_forecasts_ets_author.R
# ========================================
# ETS Base Forecasts — Faithful to Original Author Code
# Matches the logic in 0-base_forecasts_ets_author_code.R exactly:
#   - L0+L1: Pure ETS via thief:::th.forecast(), no clipping
#   - L2: Hybrid ETS+NWP (NWP replaces hourly), no clipping
#   - Residuals stored as list per replication (author format)
# Adapted only for: our file paths, parallel setup, and output naming
# ========================================

rm(list = ls(all = TRUE))

# Load required packages (same as author + doSNOW for progress)
libs <- c("thief", "data.table", "doParallel", "doSNOW", "hts")
invisible(lapply(libs, library, character.only = TRUE))

# Set working directory to RF_SARIMAX
script_dir <- "RF_SARIMAX"
if (!grepl("RF_SARIMAX$", getwd())) {
  if (dir.exists(script_dir)) {
    setwd(script_dir)
    cat("Changed working directory to:", getwd(), "\n")
  }
}

# Source configuration
source("config.R")

cat("\n========================================\n")
cat("ETS Base Forecast Generation (Author Version)\n")
cat("========================================\n\n")

# ----------------------------------------
# Load data (same as author)
# ----------------------------------------
meas_raw <- fread(DATA_PATH, header = TRUE, sep = ",")
if (is.character(meas_raw[[1]]) || inherits(meas_raw[[1]], "POSIXct")) {
  meas <- as.data.frame(meas_raw[, -1])
} else {
  meas <- as.data.frame(meas_raw)
}

# Build hts object (same as author)
level1 <- max(as.numeric(substr(colnames(meas), 1, 2)))
level2 <- rle(as.numeric(substr(colnames(meas), 1, 2)))$length
nodes <- list(level1, level2)
temp <- ts(meas)
xx <- hts(temp, nodes, characters = c(2, 4))

cat(sprintf("Stations: %d, Regions: %d\n", ncol(meas), level1))
cat(sprintf("Replications: %d to %d\n", min(rep_range), max(rep_range)))

# ----------------------------------------
# Constants (same as author)
# ----------------------------------------
k.v <- c(1, 2, 3, 4, 6, 8, 12, 24)
k.s <- sum(k.v)  # 60
Mk.v <- m / k.v
n_rep <- length(rep_range)

# Output directory (separate from our ETS to allow comparison)
dir.create(dir_ets_author, showWarnings = FALSE, recursive = TRUE)

#################################################################################
# BLOCK 1: L0 + L1 — Pure ETS (matches author lines 32-80)
#################################################################################
cat("\n========================================\n")
cat("Processing Aggregated Series (L0 + L1) — Pure ETS\n")
cat("========================================\n")

L01 <- aggts(xx, levels = c(0, 1))
agg_names <- c("Total", paste0("Region_", sprintf("%02d", 1:level1)))

for (st in 1:ncol(L01)) {

  agg_name <- agg_names[st]
  cat(sprintf("\n[Aggregate %d/%d] %s\n", st, ncol(L01), agg_name))

  Y <- ts(L01[, st], freq = m)

  # Setup parallel
  # for local machine
  cl <- makeCluster(ncores)
  # for SLURM cluster
  # cl <- parallel::makeCluster(ncores, type = "PSOCK")
  registerDoSNOW(cl)
  clusterEvalQ(cl, { library(thief) })
  clusterExport(cl, c("Y", "m", "k.v", "k.s", "Mk.v",
                       "train.days", "forecast.horizon"))

  pb <- txtProgressBar(max = n_rep, style = 3)
  progress_fn <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress_fn)

  results <- foreach(i = rep_range, .options.snow = opts,
                     .errorhandling = "pass") %dopar% {

    tryCatch({
      # Window training data (same as author line 52)
      y <- window(Y, start = c(i, 1), end = c(i + train.days - 1, 24))
      yk <- thief::tsaggregates(y, m = m, align = "end")

      # Fit ETS at all levels (author line 54)
      frc <- thief:::th.forecast(aggy = yk, h = m * forecast.horizon,
                                 usemodel = "ets", forecastfunction = NULL)
      frc$mse[1] <- mean(frc$residuals[[1]]^2)

      # Build Y.hat (author lines 57-64)
      Y.hat <- numeric(k.s * forecast.horizon)
      for (hh in 1:forecast.horizon) {
        tmp <- NULL
        for (k in 1:length(k.v)) {
          tmp <- c(tmp, frc$forecast[[length(k.v) - k + 1]][
            (Mk.v[length(k.v) - k + 1] * (hh - 1) + 1):
            (Mk.v[length(k.v) - k + 1] * hh)])
        }
        Y.hat[(k.s * (hh - 1) + 1):(k.s * hh)] <- tmp
      }

      # Y.obs (author line 65)
      Y.obs <- as.numeric(window(Y, start = c(i + train.days, 1),
                                 end = c(i + train.days + forecast.horizon - 1, 24)))

      # Res.insamp — convert from author's list format to flat day-major vector
      # (same values, just rearranged for compatibility with reconciliation scripts)
      Res.insamp <- c()
      for (day in 1:train.days) {
        for (k_idx in length(k.v):1) {
          n_periods <- m / k.v[k_idx]
          res_at_level <- as.numeric(frc$residuals[[k_idx]])
          start_r <- (day - 1) * n_periods + 1
          end_r <- day * n_periods
          Res.insamp <- c(Res.insamp, res_at_level[start_r:end_r])
        }
      }

      list(Y.hat = Y.hat, Y.obs = Y.obs, Res.insamp = Res.insamp, error = NULL)

    }, error = function(e) {
      list(Y.hat = rep(NA, k.s * forecast.horizon),
           Y.obs = rep(NA, m * forecast.horizon),
           Res.insamp = rep(NA, k.s * train.days),
           error = as.character(e))
    })
  }

  close(pb)
  stopCluster(cl)

  errors <- sapply(results, function(x) !is.null(x$error))
  if (any(errors)) {
    cat(sprintf("  Warning: %d replications had errors\n", sum(errors)))
  }

  filename <- sprintf("%s/%s--0--ets.RData", dir_ets_author, agg_name)
  save(results, file = filename)
  cat(sprintf("  Saved to: %s\n", filename))
}

#################################################################################
# BLOCK 2: L2 — Hybrid ETS+NWP (matches author lines 82-140)
#################################################################################
cat("\n========================================\n")
cat("Processing Bottom-Level Stations (L2) — Hybrid ETS+NWP\n")
cat("========================================\n")

pred_raw <- fread(PRED_PATH, header = TRUE, sep = ",")
if (is.character(pred_raw[[1]]) || inherits(pred_raw[[1]], "POSIXct")) {
  pred <- as.data.frame(pred_raw[, -1])
} else {
  pred <- as.data.frame(pred_raw)
}

for (st in 1:ncol(meas)) {

  station_name <- colnames(meas)[st]
  cat(sprintf("\n[Station %d/%d] %s\n", st, ncol(meas), station_name))

  Y <- ts(meas[, st], freq = m)
  NWP <- ts(pred[, st], freq = m)

  # Setup parallel
  cl <- makeCluster(ncores)
  registerDoSNOW(cl)
  clusterEvalQ(cl, { library(thief) })
  clusterExport(cl, c("Y", "NWP", "m", "k.v", "k.s", "Mk.v",
                       "train.days", "forecast.horizon"))

  pb <- txtProgressBar(max = n_rep, style = 3)
  progress_fn <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress_fn)

  results <- foreach(i = rep_range, .options.snow = opts,
                     .errorhandling = "pass") %dopar% {

    tryCatch({
      # Window training data (author line 103)
      y <- window(Y, start = c(i, 1), end = c(i + train.days - 1, 24))
      yk <- thief::tsaggregates(y, m = m, align = "end")

      # Fit ETS at all levels (author line 105)
      frc <- thief:::th.forecast(aggy = yk, h = m * forecast.horizon,
                                 usemodel = "ets", forecastfunction = NULL)

      # Replace hourly (k=1) with NWP (author lines 107-109)
      frc$forecast[[1]] <- window(NWP, start = c(i + train.days, 1),
                                  end = c(i + train.days + forecast.horizon - 1, 24))
      frc$residuals[[1]] <- window(NWP, start = c(i, 1),
                                   end = c(i + train.days - 1, 24)) - y
      frc$mse[1] <- mean(frc$residuals[[1]]^2)

      # Build Y.hat (author lines 111-116)
      Y.hat <- numeric(k.s * forecast.horizon)
      for (hh in 1:forecast.horizon) {
        tmp <- NULL
        for (k in 1:length(k.v)) {
          tmp <- c(tmp, frc$forecast[[length(k.v) - k + 1]][
            (Mk.v[length(k.v) - k + 1] * (hh - 1) + 1):
            (Mk.v[length(k.v) - k + 1] * hh)])
        }
        Y.hat[(k.s * (hh - 1) + 1):(k.s * hh)] <- tmp
      }

      # Y.obs (author line 121)
      Y.obs <- as.numeric(window(Y, start = c(i + train.days, 1),
                                 end = c(i + train.days + forecast.horizon - 1, 24)))

      # Res.insamp — convert from author's list format to flat day-major vector
      # (same values, just rearranged for compatibility with reconciliation scripts)
      Res.insamp <- c()
      for (day in 1:train.days) {
        for (k_idx in length(k.v):1) {
          n_periods <- m / k.v[k_idx]
          res_at_level <- as.numeric(frc$residuals[[k_idx]])
          start_r <- (day - 1) * n_periods + 1
          end_r <- day * n_periods
          Res.insamp <- c(Res.insamp, res_at_level[start_r:end_r])
        }
      }

      list(Y.hat = Y.hat, Y.obs = Y.obs, Res.insamp = Res.insamp, error = NULL)

    }, error = function(e) {
      list(Y.hat = rep(NA, k.s * forecast.horizon),
           Y.obs = rep(NA, m * forecast.horizon),
           Res.insamp = rep(NA, k.s * train.days),
           error = as.character(e))
    })
  }

  close(pb)
  stopCluster(cl)

  errors <- sapply(results, function(x) !is.null(x$error))
  if (any(errors)) {
    cat(sprintf("  Warning: %d replications had errors\n", sum(errors)))
  }

  filename <- sprintf("%s/%s--ets.RData", dir_ets_author, station_name)
  save(results, file = filename)
  cat(sprintf("  Saved to: %s\n", filename))
}

cat("\n========================================\n")
cat("ETS Base Forecasts Complete (Author Version)\n")
cat("========================================\n")
cat(sprintf("Aggregated series: %d\n", ncol(L01)))
cat(sprintf("Bottom-level stations: %d\n", ncol(meas)))
cat(sprintf("Replications per series: %d\n", n_rep))
cat(sprintf("Results saved to: %s\n", dir_ets_author))
cat("\nKey differences from our pipeline version:\n")
cat("  - No pmax(Y.hat, 0) clipping on base forecasts\n")
cat("  - Uses window() for ts indexing (same as author)\n")
cat("  - Uses hts::aggts() for L0+L1 aggregation (same as author)\n")
