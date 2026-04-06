# 1-reco_ctbu_twlsv.R
# ========================================
# Cross-Temporal Bottom-Up T-WLSV Reconciliation
# Two-step approach:
#   1. Temporal reconciliation (thfrec) on bottom-level series only
#   2. Spatial bottom-up aggregation using S matrix
# ========================================

rm(list = ls())

# Load required packages
libs <- c("FoReco", "doParallel", "doSNOW", "progress", "Matrix")
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
cat("CTBU T-WLSV Reconciliation\n")
cat("(Temporal reconciliation + Spatial bottom-up)\n")
cat("========================================\n\n")

# ----------------------------------------
# CRITICAL: Format transformation functions
# Our forecasts are in "day-major" order: [day1_all_k, day2_all_k]
# FoReco expects "k-major" order: [all_days_k24, all_days_k12, ..., all_days_k1]
# ----------------------------------------
create_foreco_indices <- function(m, n_days, k_order) {
  periods_per_k <- m / k_order
  k_star <- sum(periods_per_k)

  to_foreco <- numeric(k_star * n_days)

  foreco_idx <- 1
  for (k in k_order) {
    n_periods <- m / k
    k_start_in_day <- cumsum(c(0, periods_per_k[-length(periods_per_k)]))[which(k_order == k)]

    for (day in 1:n_days) {
      day_offset <- (day - 1) * k_star
      for (p in 1:n_periods) {
        our_row <- day_offset + k_start_in_day + p
        to_foreco[foreco_idx] <- our_row
        foreco_idx <- foreco_idx + 1
      }
    }
  }

  # Create inverse mapping
  from_foreco <- numeric(k_star * n_days)
  from_foreco[to_foreco] <- 1:(k_star * n_days)

  list(to_foreco = to_foreco, from_foreco = from_foreco)
}

# Pre-compute transformation indices for forecasts (h days) and residuals (train.days)
k_order <- c(24, 12, 8, 6, 4, 3, 2, 1)
idx_forecast <- create_foreco_indices(m, h, k_order)
idx_residuals <- create_foreco_indices(m, train.days, k_order)

cat("Format transformation indices created.\n")

# ----------------------------------------
# Helper function to load bottom-level forecasts only
# ----------------------------------------
load_bottom_forecasts <- function(rp, method = "sarimax") {
  # Determine directory
  if (method == "sarimax") {
    dir_path <- dir_sarimax
  } else if (method == "rf") {
    dir_path <- dir_rf
  } else if (method == "rf_nwp") {
    dir_path <- dir_rf_nwp
  } else if (method == "lgbm") {
    dir_path <- dir_lgbm
  } else if (method == "ets") {
    dir_path <- dir_ets
  } else if (method == "ets_author") {
    dir_path <- dir_ets_author
  } else if (method == "sarimax_nwp") {
    dir_path <- dir_sarimax_nwp
  } else {
    stop(sprintf("Unknown method: %s", method))
  }

  # Get bottom-level files only (exclude --0-- aggregated files)
  files <- list.files(dir_path, pattern = "\\.RData$", full.names = TRUE)
  bottom_files <- files[!grepl("--0--", files)]

  if (length(bottom_files) == 0) {
    stop(sprintf("No bottom-level result files found in %s", dir_path))
  }

  # Initialize matrices for 318 bottom-level series
  n_bottom <- length(bottom_files)
  n_yhat <- 120  # h * k_star
  n_res <- 840   # train.days * k_star

  Yhat <- matrix(NA, nrow = n_yhat, ncol = n_bottom)
  E <- matrix(NA, nrow = n_res, ncol = n_bottom)

  col_idx <- 1
  for (f in bottom_files) {
    load(f)  # loads 'results'
    if (rp <= length(results) && !is.null(results[[rp]]$Y.hat)) {
      Yhat[, col_idx] <- results[[rp]]$Y.hat
      E[, col_idx] <- results[[rp]]$Res.insamp
    }
    col_idx <- col_idx + 1
  }

  return(list(Yhat = Yhat, E = E, n_bottom = n_bottom))
}

# ----------------------------------------
# Get S matrix for spatial aggregation
# ----------------------------------------
S <- as.matrix(hts_info$S)  # 324 x 318 matrix
n_upper <- nrow(S) - ncol(S)  # 6 upper-level series
n_bottom <- ncol(S)  # 318 bottom-level series

cat(sprintf("S matrix dimensions: %d x %d\n", nrow(S), ncol(S)))
cat(sprintf("Upper-level series: %d\n", n_upper))
cat(sprintf("Bottom-level series: %d\n", n_bottom))

# ----------------------------------------
# Process all base methods
# ----------------------------------------
for (method in c("sarimax", "rf", "rf_nwp", "lgbm", "ets", "ets_author", "sarimax_nwp")) {

  cat(sprintf("\n========================================\n"))
  cat(sprintf("Processing CTBU T-WLSV for %s\n", toupper(method)))
  cat(sprintf("========================================\n\n"))

  if (method == "sarimax") {
    output_dir <- dir_ctbu_sarimax
  } else if (method == "rf") {
    output_dir <- dir_ctbu_rf
  } else if (method == "rf_nwp") {
    output_dir <- dir_ctbu_rf_nwp
  } else if (method == "lgbm") {
    output_dir <- dir_ctbu_lgbm
  } else if (method == "ets") {
    output_dir <- dir_ctbu_ets
  } else if (method == "ets_author") {
    output_dir <- dir_ctbu_ets_author
  } else if (method == "sarimax_nwp") {
    output_dir <- dir_ctbu_sarimax_nwp
  }

  # Setup parallel processing
  # for local machine
  # cl <- makeCluster(ncores)
  # for SLURM cluster
  cl <- parallel::makeCluster(ncores, type = "PSOCK")
  registerDoSNOW(cl)

  clusterEvalQ(cl, {
    library(FoReco)
    library(Matrix)
  })

  clusterExport(cl, c("hts_info", "thf_info", "m", "h", "k.v",
                      "drop_zeros",
                      "dir_sarimax", "dir_rf", "dir_rf_nwp", "dir_lgbm", "dir_ets",
                      "dir_ets_author", "dir_sarimax_nwp",
                      "load_bottom_forecasts",
                      "method", "S", "n_upper", "n_bottom",
                      "idx_forecast", "idx_residuals"))

  # Progress bar
  pb <- txtProgressBar(max = length(rep_range), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  # Process all replications
  results <- foreach(rp = rep_range, .options.snow = opts,
                     .errorhandling = "pass") %dopar% {

    tryCatch({
      # Load bottom-level forecasts only
      bottom_data <- load_bottom_forecasts(rp, method = method)

      # Check for NA values
      if (any(is.na(bottom_data$Yhat)) || any(is.na(bottom_data$E))) {
        return(list(free = NULL, nn = NULL, time = NA,
                    error = "NA values in input data"))
      }

      Start <- Sys.time()

      # Step 1: Temporal reconciliation on each bottom-level series
      # CRITICAL: Transform to FoReco k-major format before thfrec
      # Apply thfrec to each of the 318 series

      # Transform data to FoReco format
      Yhat_foreco <- bottom_data$Yhat[idx_forecast$to_foreco, ]
      E_foreco <- bottom_data$E[idx_residuals$to_foreco, ]

      thf_reco_foreco <- matrix(NA, nrow = nrow(Yhat_foreco), ncol = n_bottom)

      for (j in 1:n_bottom) {
        # thfrec expects data in FoReco k-major format
        thf_result <- FoReco::thfrec(
          basef = Yhat_foreco[, j],
          comb = "wlsv",
          m = m,
          res = E_foreco[, j],
          keep = "recf"
        )

        thf_reco_foreco[, j] <- as.numeric(thf_result)
      }

      # Transform back to our day-major format
      thf_reco <- thf_reco_foreco[idx_forecast$from_foreco, ]

      # Step 2: Spatial bottom-up aggregation
      # S is 324 x 318: maps 318 bottom-level to 324 total series
      # For each temporal level (row), multiply by S

      # thf_reco is [120 x 318], we need [120 x 324]
      # Each row of thf_reco is one temporal level
      # We need to aggregate each temporal level spatially
      reco_full <- matrix(NA, nrow = nrow(thf_reco), ncol = nrow(S))

      for (row in 1:nrow(thf_reco)) {
        # S is 324 x 318, thf_reco[row, ] is 318 x 1
        # Result is 324 x 1
        reco_full[row, ] <- as.numeric(S %*% thf_reco[row, ])
      }

      time_free <- difftime(Sys.time(), Start, units = "secs")

      # Clean up near-zero values
      reco_free <- drop_zeros(reco_full)

      # SNTZ heuristic for non-negativity (as in the paper):
      # Step 1 already done: temporal reco + spatial BU (coherent)
      # Step 2: set negatives to zero on bottom-level hourly
      # Step 3: re-aggregate via spatio-temporal bottom-up
      reco_nn <- NULL
      time_nn <- NA

      if (any(reco_free < 0, na.rm = TRUE)) {
        Start <- Sys.time()

        # Step 2: Set negatives to zero on bottom-level hourly
        # thf_reco is [120 x 318] in day-major format
        # Extract hourly (k=1) rows and clip
        k_order <- sort(k.v, decreasing = TRUE)
        periods_per_k <- m / k_order
        k_star <- sum(periods_per_k)
        bottom_nn <- pmax(thf_reco, 0)

        # Step 3: Re-aggregate from clipped bottom-level hourly
        reco_sntz <- matrix(NA,
                            nrow = nrow(reco_free),
                            ncol = nrow(S))

        for (day in 1:h) {
          day_offset <- (day - 1) * k_star

          # Hourly rows: last 24 of each day block
          hourly_offset <- sum(periods_per_k[
            1:(length(k_order) - 1)
          ])
          hourly_start <- day_offset + hourly_offset + 1
          hourly_end <- day_offset + k_star
          hourly_bottom <- bottom_nn[
            hourly_start:hourly_end, ]  # [24 x 318]

          # Temporal agg + spatial BU for each k level
          row_offset <- day_offset
          for (ki in seq_along(k_order)) {
            kk <- k_order[ki]
            n_per <- periods_per_k[ki]

            if (kk == 1) {
              for (p in 1:n_per) {
                reco_sntz[row_offset + p, ] <-
                  as.numeric(S %*% hourly_bottom[p, ])
              }
            } else {
              for (p in 1:n_per) {
                h_start <- (p - 1) * kk + 1
                h_end <- p * kk
                agg_b <- colSums(
                  hourly_bottom[h_start:h_end, ,
                                drop = FALSE])
                reco_sntz[row_offset + p, ] <-
                  as.numeric(S %*% agg_b)
              }
            }
            row_offset <- row_offset + n_per
          }
        }

        time_nn <- difftime(Sys.time(), Start, units = "secs")
        reco_nn <- drop_zeros(reco_sntz)
      }

      list(
        free = reco_free,
        nn = reco_nn,
        time = as.numeric(time_free),
        time_nn = as.numeric(time_nn),
        error = NULL
      )

    }, error = function(e) {
      list(free = NULL, nn = NULL, time = NA,
           time_nn = NA, error = as.character(e))
    })
  }

  close(pb)

  # Check for errors
  errors <- sapply(results, function(x) !is.null(x$error) && !is.na(x$error))
  if (any(errors)) {
    cat(sprintf("\nWarning: %d replications had errors\n", sum(errors)))
    error_msgs <- unique(sapply(results[errors], function(x) x$error))
    cat("Error types:\n")
    for (msg in error_msgs) cat(sprintf("  - %s\n", msg))
  }

  # Save results
  filename <- sprintf("%s/ctbu_twlsv_%s_results.RData", output_dir, method)
  save(results, file = filename)
  cat(sprintf("\nResults saved to: %s\n", filename))

  # Summary statistics
  times <- sapply(results, function(x) x$time)
  times <- times[!is.na(times)]

  cat(sprintf("\nTiming summary:\n"))
  cat(sprintf("  Mean: %.2f seconds\n", mean(times)))
  cat(sprintf("  Total: %.2f minutes\n", sum(times) / 60))

  n_nn <- sum(sapply(results, function(x) !is.null(x$nn)))
  cat(sprintf("\nReplications requiring clipping: %d / %d\n",
              n_nn, length(results)))

  stopCluster(cl)
}

cat("\n========================================\n")
cat("CTBU T-WLSV Reconciliation Complete\n")
cat("========================================\n")
