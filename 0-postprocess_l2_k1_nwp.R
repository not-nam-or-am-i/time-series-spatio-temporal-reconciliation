# 0-postprocess_l2_k1_nwp.R
# ========================================
# Post-process base forecasts: replace L2 bottom-level k=1 (hourly) blocks in
# Y.hat and Res.insamp with NWP-based values (same hybrid as sarimax_nwp).
#
# Inputs: Results from 0-base_forecasts_sarimax.R, 0-base_forecasts_rf.R,
#         0-base_forecasts_lgbm.R (bottom-level + aggregates copied).
# Outputs: Results_SARIMAX_NWP, Results_RF_NWP, Results_LGBM_NWP with *_nwp filenames.
#
# Usage:
#   Rscript 0-postprocess_l2_k1_nwp.R
#   Rscript 0-postprocess_l2_k1_nwp.R sarimax
#   Rscript 0-postprocess_l2_k1_nwp.R sarimax rf lgbm
# ========================================

rm(list = ls())

libs <- c("data.table")
invisible(lapply(libs, library, character.only = TRUE))

script_dir <- "RF_SARIMAX"
if (!grepl("RF_SARIMAX$", getwd())) {
  if (dir.exists(script_dir)) {
    setwd(script_dir)
    cat("Changed working directory to:", getwd(), "\n")
  }
}

source("config.R")

# ----------------------------------------
# k=1 slice layout (day-major Y.hat / Res.insamp from assemble_*_from_levels)
# ----------------------------------------
patch_yhat_k1 <- function(yhat, test_nwp_hourly, k_star, m, h) {
  if (length(yhat) != h * k_star) {
    stop(sprintf("Y.hat length %d != h*k_star %d", length(yhat), h * k_star))
  }
  if (length(test_nwp_hourly) != h * m) {
    stop(sprintf("test NWP length %d != h*m %d", length(test_nwp_hourly), h * m))
  }
  out <- yhat
  fc1 <- pmax(test_nwp_hourly, 0)
  for (d in seq_len(h)) {
    i0 <- (d - 1) * k_star + (k_star - m) + 1
    i1 <- d * k_star
    out[i0:i1] <- fc1[((d - 1) * m + 1):(d * m)]
  }
  out
}

patch_res_k1 <- function(res, train_nwp_minus_y, k_star, m, train_days) {
  if (length(res) != train.days * k_star) {
    stop(sprintf("Res.insamp length %d != train.days*k_star %d",
                 length(res), train.days * k_star))
  }
  if (length(train_nwp_minus_y) != train.days * m) {
    stop(sprintf("train residual length %d != train.days*m %d",
                 length(train_nwp_minus_y), train.days * m))
  }
  out <- res
  for (td in seq_len(train.days)) {
    i0 <- (td - 1) * k_star + (k_star - m) + 1
    i1 <- td * k_star
    out[i0:i1] <- train_nwp_minus_y[((td - 1) * m + 1):(td * m)]
  }
  out
}

copy_aggregate_files <- function(in_dir, out_dir, base_suffix, nwp_suffix) {
  patt <- sprintf("^.*--0--%s\\.RData$", base_suffix)
  files <- list.files(in_dir, pattern = patt, full.names = TRUE)
  if (length(files) == 0) {
    warning(sprintf("No aggregate files matching --0--%s.RData in %s",
                    base_suffix, in_dir))
    return(invisible(0L))
  }
  n <- 0L
  for (f in files) {
    bn <- basename(f)
    new_bn <- sub(sprintf("--0--%s\\.RData$", base_suffix),
                  sprintf("--0--%s.RData", nwp_suffix),
                  bn)
    dest <- file.path(out_dir, new_bn)
    file.copy(f, dest, overwrite = TRUE)
    n <- n + 1L
    cat(sprintf("  Copied aggregate: %s -> %s\n", bn, new_bn))
  }
  invisible(n)
}

process_station_file <- function(path_in, out_dir, in_suffix, out_suffix,
                                   station_meas, station_pred) {
  load(path_in) # results
  if (!exists("results")) {
    stop(sprintf("No 'results' in %s", path_in))
  }
  k_star <- sum(m / k.v)
  obs_per_day <- m

  for (ii in seq_along(rep_range)) {
    rp <- rep_range[ii]
    if (ii > length(results)) break
    r <- results[[ii]]
    if (is.null(r) || !is.null(r$error)) next
    if (is.null(r$Y.hat) || all(is.na(r$Y.hat))) next

    start_idx <- (rp - 1) * obs_per_day + 1
    end_train <- start_idx + train.days * obs_per_day - 1
    end_test <- end_train + h * obs_per_day

    if (end_test > length(station_meas)) next

    train_y_hourly <- station_meas[start_idx:end_train]
    train_x_hourly <- station_pred[start_idx:end_train]
    test_x_hourly <- station_pred[(end_train + 1):end_test]

    r$Y.hat <- patch_yhat_k1(r$Y.hat, test_x_hourly, k_star, m, h)
    nwp_res <- train_x_hourly - train_y_hourly
    r$Res.insamp <- patch_res_k1(r$Res.insamp, nwp_res, k_star, m, train.days)
    results[[ii]] <- r
  }

  station_name <- sub(sprintf("--%s\\.RData$", in_suffix), "", basename(path_in))
  out_path <- file.path(out_dir, sprintf("%s--%s.RData", station_name, out_suffix))
  save(results, file = out_path)
  cat(sprintf("  Wrote: %s\n", out_path))
}

# ----------------------------------------
# Main
# ----------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  models_to_run <- c("sarimax", "rf", "lgbm")
} else {
  models_to_run <- tolower(args)
}

model_map <- list(
  sarimax = list(
    in_dir = dir_sarimax, out_dir = dir_sarimax_nwp,
    in_suffix = "sarimax", out_suffix = "sarimax_nwp"
  ),
  rf = list(
    in_dir = dir_rf, out_dir = dir_rf_nwp,
    in_suffix = "rf", out_suffix = "rf_nwp"
  ),
  lgbm = list(
    in_dir = dir_lgbm, out_dir = dir_lgbm_nwp,
    in_suffix = "lgbm", out_suffix = "lgbm_nwp"
  )
)

meas_raw <- fread(DATA_PATH)
if (is.character(meas_raw[[1]]) || inherits(meas_raw[[1]], "POSIXct")) {
  meas <- as.matrix(meas_raw[, -1])
} else {
  meas <- as.matrix(meas_raw)
}
pred_raw <- fread(PRED_PATH)
if (is.character(pred_raw[[1]]) || inherits(pred_raw[[1]], "POSIXct")) {
  pred <- as.matrix(pred_raw[, -1])
} else {
  pred <- as.matrix(pred_raw)
}

if (nrow(meas) != nrow(pred) || ncol(meas) != ncol(pred)) {
  stop("meas and pred dimensions must match.")
}

k_star <- sum(m / k.v)
cat(sprintf("\nk_star = %d (m=%d, h=%d, train.days=%d)\n", k_star, m, h, train.days))
cat(sprintf("Models: %s\n\n", paste(models_to_run, collapse = ", ")))

for (mod in models_to_run) {
  if (!mod %in% names(model_map)) {
    warning(sprintf("Unknown model '%s', skipping.", mod))
    next
  }
  cfg <- model_map[[mod]]
  if (!dir.exists(cfg$in_dir)) {
    warning(sprintf("Input dir missing: %s — skipping %s", cfg$in_dir, mod))
    next
  }
  if (!dir.exists(cfg$out_dir)) {
    dir.create(cfg$out_dir, recursive = TRUE)
    cat(sprintf("Created: %s\n", cfg$out_dir))
  }

  cat(sprintf("\n========================================\n"))
  cat(sprintf("Post-process: %s -> %s\n", mod, cfg$out_suffix))
  cat("========================================\n")

  all_r <- list.files(cfg$in_dir, pattern = "\\.RData$", full.names = TRUE)
  bottom_files <- all_r[!grepl("--0--", all_r, fixed = TRUE)]
  bottom_files <- bottom_files[grepl(sprintf("--%s\\.RData$", cfg$in_suffix), basename(bottom_files))]

  if (length(bottom_files) == 0) {
    warning(sprintf("No bottom-level files in %s for suffix %s", cfg$in_dir, cfg$in_suffix))
    next
  }

  for (f in bottom_files) {
    station_name <- sub(sprintf("--%s\\.RData$", cfg$in_suffix), "", basename(f))
    if (!station_name %in% colnames(meas)) {
      warning(sprintf("Station '%s' not in meas columns, skipping %s", station_name, f))
      next
    }
    sm <- meas[, station_name, drop = TRUE]
    sp <- pred[, station_name, drop = TRUE]
    process_station_file(f, cfg$out_dir, cfg$in_suffix, cfg$out_suffix, sm, sp)
  }

  copy_aggregate_files(cfg$in_dir, cfg$out_dir, cfg$in_suffix, cfg$out_suffix)
}

cat("\n========================================\n")
cat("Post-processing complete.\n")
cat("========================================\n")

# ----------------------------------------
# Optional validation: patch(base SARIMAX) vs reference from 0-base_forecasts_sarimax_nwp.R
# Requires: Results_SARIMAX/{station}--sarimax.RData and Results_SARIMAX_NWP/{station}--sarimax_nwp.RData
# ii indexes results[[ii]] for rep_range[ii].
# Example: validate_postprocess_vs_sarimax_nwp("MyStation", ii = 1L)
# ----------------------------------------
validate_postprocess_vs_sarimax_nwp <- function(station_name, ii = 1L) {
  e <- new.env()
  sys.source("config.R", envir = e)
  meas_raw <- data.table::fread(e$DATA_PATH)
  pred_raw <- data.table::fread(e$PRED_PATH)
  if (is.character(meas_raw[[1]]) || inherits(meas_raw[[1]], "POSIXct")) {
    meas <- as.matrix(meas_raw[, -1])
  } else {
    meas <- as.matrix(meas_raw)
  }
  if (is.character(pred_raw[[1]]) || inherits(pred_raw[[1]], "POSIXct")) {
    pred <- as.matrix(pred_raw[, -1])
  } else {
    pred <- as.matrix(pred_raw)
  }
  m <- e$m
  h <- e$h
  k.v <- e$k.v
  train.days <- e$train.days
  rep_range <- e$rep_range
  dir_sarimax <- e$dir_sarimax
  dir_sarimax_nwp <- e$dir_sarimax_nwp

  ref <- file.path(dir_sarimax_nwp, sprintf("%s--sarimax_nwp.RData", station_name))
  base_f <- file.path(dir_sarimax, sprintf("%s--sarimax.RData", station_name))
  if (!file.exists(ref) || !file.exists(base_f)) {
    message("Need reference sarimax_nwp and base sarimax files for this station.")
    return(invisible(NULL))
  }
  env_ref <- new.env(parent = emptyenv())
  env_base <- new.env(parent = emptyenv())
  load(ref, envir = env_ref)
  load(base_f, envir = env_base)
  r_ref <- env_ref$results[[ii]]
  r_base <- env_base$results[[ii]]

  k_star <- sum(m / k.v)
  obs_per_day <- m
  rp <- rep_range[ii]
  start_idx <- (rp - 1) * obs_per_day + 1
  end_train <- start_idx + train.days * obs_per_day - 1
  end_test <- end_train + h * obs_per_day
  sm <- meas[, station_name, drop = TRUE]
  sp <- pred[, station_name, drop = TRUE]
  train_y <- sm[start_idx:end_train]
  train_x <- sp[start_idx:end_train]
  test_x <- sp[(end_train + 1):end_test]
  yhat_patch <- patch_yhat_k1(r_base$Y.hat, test_x, k_star, m, h)
  res_patch <- patch_res_k1(r_base$Res.insamp, train_x - train_y, k_star, m, train.days)
  cat(sprintf("max |Y.hat diff|: %g\n", max(abs(yhat_patch - r_ref$Y.hat), na.rm = TRUE)))
  cat(sprintf("max |Res diff|: %g\n", max(abs(res_patch - r_ref$Res.insamp), na.rm = TRUE)))
  invisible(max(abs(yhat_patch - r_ref$Y.hat), na.rm = TRUE))
}
