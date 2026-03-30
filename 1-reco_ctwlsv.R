# 1-reco_ctwlsv.R
# ========================================
# Cross-Temporal WLSV Reconciliation
# Joint optimal reconciliation for both spatial and temporal dimensions
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
cat("Cross-Temporal WLSV Reconciliation\n")
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
# Helper function to load all forecasts for a replication
# ----------------------------------------
load_all_forecasts <- function(rp, method = "sarimax") {
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

  # Get all result files
  files <- list.files(dir_path, pattern = "\\.RData$", full.names = TRUE)

  if (length(files) == 0) {
    stop(sprintf("No result files found in %s", dir_path))
  }

  # Separate aggregated (--0--) from bottom-level files
  agg_files <- files[grepl("--0--", files)]
  bottom_files <- files[!grepl("--0--", files)]

  # IMPORTANT: Sort aggregated files in correct order to match S matrix structure
  # S matrix has: Total (col 1), Region_01 (col 2), ..., Region_05 (col 6)
  agg_names <- c("Total", "Region_01", "Region_02", "Region_03", "Region_04", "Region_05")
  agg_files_ordered <- character(length(agg_names))
  for (i in seq_along(agg_names)) {
    pattern <- paste0(agg_names[i], "--0--")
    matching <- agg_files[grepl(pattern, agg_files, fixed = TRUE)]
    if (length(matching) > 0) {
      agg_files_ordered[i] <- matching[1]
    }
  }
  agg_files_ordered <- agg_files_ordered[agg_files_ordered != ""]

  # Initialize matrices
  # Order: 6 aggregated (Total + 5 regions) + 318 bottom-level = 324
  n_total <- length(agg_files_ordered) + length(bottom_files)
  n_yhat <- 120  # h * k_star
  n_res <- 840   # train.days * k_star

  Yhat <- matrix(NA, nrow = n_yhat, ncol = n_total)
  E <- matrix(NA, nrow = n_res, ncol = n_total)

  col_idx <- 1

  # Load aggregated series first (L0 + L1) in correct order
  for (f in agg_files_ordered) {
    load(f)  # loads 'results'
    if (rp <= length(results) && !is.null(results[[rp]]$Y.hat)) {
      Yhat[, col_idx] <- results[[rp]]$Y.hat
      E[, col_idx] <- results[[rp]]$Res.insamp
    }
    col_idx <- col_idx + 1
  }

  # Load bottom-level series (L2)
  for (f in bottom_files) {
    load(f)  # loads 'results'
    if (rp <= length(results) && !is.null(results[[rp]]$Y.hat)) {
      Yhat[, col_idx] <- results[[rp]]$Y.hat
      E[, col_idx] <- results[[rp]]$Res.insamp
    }
    col_idx <- col_idx + 1
  }

  return(list(Yhat = Yhat, E = E))
}

# ----------------------------------------
# Process all base methods
# ----------------------------------------
for (method in c("sarimax", "rf", "rf_nwp", "lgbm", "ets", "ets_author", "sarimax_nwp")) {

  cat(sprintf("\n========================================\n"))
  cat(sprintf("Processing CTWLSV for %s\n", toupper(method)))
  cat(sprintf("========================================\n\n"))

  if (method == "sarimax") {
    output_dir <- dir_ctwlsv_sarimax
  } else if (method == "rf") {
    output_dir <- dir_ctwlsv_rf
  } else if (method == "rf_nwp") {
    output_dir <- dir_ctwlsv_rf_nwp
  } else if (method == "lgbm") {
    output_dir <- dir_ctwlsv_lgbm
  } else if (method == "ets") {
    output_dir <- dir_ctwlsv_ets
  } else if (method == "ets_author") {
    output_dir <- dir_ctwlsv_ets_author
  } else if (method == "sarimax_nwp") {
    output_dir <- dir_ctwlsv_sarimax_nwp
  }

  # Setup parallel processing
  cl <- makeCluster(ncores)
  registerDoSNOW(cl)

  clusterEvalQ(cl, {
    library(FoReco)
    library(Matrix)
  })

  # Get S matrix for sntz bottom-up re-aggregation
  S <- as.matrix(hts_info$S)  # 324 x 318
  n_upper <- nrow(S) - ncol(S)  # 6
  n_bottom <- ncol(S)  # 318

  clusterExport(cl, c("hts_info", "thf_info", "m", "h", "drop_zeros",
                      "dir_sarimax", "dir_rf", "dir_lgbm", "dir_ets",
                      "dir_rf_nwp", "dir_ets_author", "dir_sarimax_nwp",
                      "load_all_forecasts", "method",
                      "S", "n_upper", "n_bottom", "k.v",
                      "idx_forecast", "idx_residuals"))

  # Progress bar
  pb <- txtProgressBar(max = length(rep_range), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  # Process all replications
  results <- foreach(rp = rep_range, .options.snow = opts,
                     .errorhandling = "pass") %dopar% {

    tryCatch({
      # Load all forecasts and residuals for this replication
      all_data <- load_all_forecasts(rp, method = method)

      # Check for NA values
      if (any(is.na(all_data$Yhat)) || any(is.na(all_data$E))) {
        return(list(free = NULL, nn = NULL, time = NA,
                    error = "NA values in input data"))
      }

      # CRITICAL: Transform from our "day-major" format to FoReco "k-major" format
      # Our format: [day1_all_k, day2_all_k] per row
      # FoReco format: [all_days_k24, all_days_k12, ..., all_days_k1] per row
      Yhat_foreco <- all_data$Yhat[idx_forecast$to_foreco, ]
      E_foreco <- all_data$E[idx_residuals$to_foreco, ]

      # Apply cross-temporal reconciliation using FoReco::octrec
      # basef: [n x k*h] matrix where n=series, k*h=temporal periods
      # res: [n x N] matrix where N=residual periods

      Start <- Sys.time()

      reco_raw <- FoReco::octrec(
        basef = t(Yhat_foreco),  # [324 x 120] in FoReco format
        comb = "wlsv",
        m = m,
        C = hts_info$C,
        res = t(E_foreco),       # [324 x 840] in FoReco format
        keep = "recf",
        type = "M"
      )

      time_free <- difftime(Sys.time(), Start, units = "secs")

      # Transform back: FoReco output [324 x 120] -> our format [120 x 324]
      reco_foreco <- t(reco_raw)  # [120 x 324] in FoReco k-major format
      reco_free <- reco_foreco[idx_forecast$from_foreco, ]  # [120 x 324] in our day-major format
      reco_free <- drop_zeros(reco_free)

      # SNTZ heuristic for non-negativity (3 steps, as in the paper):
      # Step 1 already done above: reconciliation (coherent, may have negatives)
      # Step 2: set negative to zero (non-negative, but incoherent)
      # Step 3: spatio-temporal bottom-up aggregation (coherent again)
      reco_nn <- NULL
      time_nn <- NA

      if (any(reco_free < 0, na.rm = TRUE)) {
        Start <- Sys.time()

        # Step 2: Set negatives to zero on bottom-level hourly forecasts
        # reco_free is [120 x 324], bottom-level = columns (n_upper+1):324
        # Extract bottom-level hourly (k=1) forecasts only
        # In day-major format, hourly values are the last 24 per day-block of 60
        bottom_reco <- reco_free[, (n_upper + 1):(n_upper + n_bottom)]
        bottom_nn <- pmax(bottom_reco, 0)

        # Step 3: Re-aggregate from bottom-level hourly via bottom-up
        # First, extract only the hourly (k=1) rows from bottom_nn
        # In day-major format: each day has k*=60 rows, hourly is last 24
        k_order <- sort(k.v, decreasing = TRUE)  # 24,12,8,6,4,3,2,1
        periods_per_k <- m / k_order  # 1,2,3,4,6,8,12,24

        # For each day, extract hourly rows and re-aggregate all temporal levels
        reco_sntz <- matrix(NA, nrow = nrow(reco_free), ncol = nrow(S))
        k_star <- sum(periods_per_k)  # 60

        for (day in 1:h) {
          day_offset <- (day - 1) * k_star

          # Hourly rows are the last 24 of each day block
          hourly_start <- day_offset + sum(periods_per_k[1:(length(k_order) - 1)]) + 1
          hourly_end <- day_offset + k_star
          hourly_bottom <- bottom_nn[hourly_start:hourly_end, ]  # [24 x 318]

          # Temporal aggregation: build all k levels from hourly
          row_offset <- day_offset
          for (ki in seq_along(k_order)) {
            kk <- k_order[ki]
            n_periods <- periods_per_k[ki]

            if (kk == 1) {
              # Hourly: no aggregation, just spatial bottom-up
              for (p in 1:n_periods) {
                reco_sntz[row_offset + p, ] <- as.numeric(S %*% hourly_bottom[p, ])
              }
            } else {
              # Aggregate kk consecutive hours, then spatial bottom-up
              for (p in 1:n_periods) {
                h_start <- (p - 1) * kk + 1
                h_end <- p * kk
                agg_bottom <- colSums(hourly_bottom[h_start:h_end, , drop = FALSE])
                reco_sntz[row_offset + p, ] <- as.numeric(S %*% agg_bottom)
              }
            }
            row_offset <- row_offset + n_periods
          }
        }

        time_nn <- difftime(Sys.time(), Start, units = "secs")
        reco_nn <- drop_zeros(reco_sntz)
      }

      list(
        free = reco_free,
        nn = reco_nn,
        time_free = as.numeric(time_free),
        time_nn = as.numeric(time_nn),
        error = NULL
      )

    }, error = function(e) {
      list(free = NULL, nn = NULL, time_free = NA, time_nn = NA,
           error = as.character(e))
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
  filename <- sprintf("%s/ctwlsv_%s_results.RData", output_dir, method)
  save(results, file = filename)
  cat(sprintf("\nResults saved to: %s\n", filename))

  # Summary statistics
  times_free <- sapply(results, function(x) x$time_free)
  times_free <- times_free[!is.na(times_free)]

  cat(sprintf("\nTiming summary (free solution):\n"))
  cat(sprintf("  Mean: %.2f seconds\n", mean(times_free)))
  cat(sprintf("  Total: %.2f minutes\n", sum(times_free) / 60))

  n_nn <- sum(sapply(results, function(x) !is.null(x$nn)))
  cat(sprintf("\nReplications requiring non-negativity: %d / %d\n",
              n_nn, length(results)))

  stopCluster(cl)
}

cat("\n========================================\n")
cat("CTWLSV Reconciliation Complete\n")
cat("========================================\n")
