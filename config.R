# config.R
# ========================================
# Shared Configuration for RF_SARIMAX
# SARIMAX (with NWP) and Random Forest (with features)
# Following PLAN.md implementation
# ========================================

# ----------------------------------------
# Data Paths
# ----------------------------------------
DATA_PATH <- "./data/meas.csv"
PRED_PATH <- "./data/pred.csv"  # NWP predictions

# ----------------------------------------
# Hierarchy Parameters
# ----------------------------------------
m <- 24                              # Seasonal period (hours per day)
h <- 2                               # Forecast horizon (days)
k.v <- c(1, 2, 3, 4, 6, 8, 12, 24)   # Temporal aggregation levels (divisors of 24)
train.days <- 14                     # Training window (days)
forecast.horizon <- 2                # Same as h (for clarity)

# ----------------------------------------
# Replication Range
# ----------------------------------------
# Full run: 350 replications (365 - 14 - 2 + 1 = 350)
rep_range <- 1:350

# Quick test: uncomment below for testing
# rep_range <- 1:3

# ----------------------------------------
# Parallel Processing
# ----------------------------------------
# for local machine
ncores <- parallel::detectCores() - 1 # for local machine
# for SLURM cluster
# ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = 8))
cat(sprintf("Using %d cores for parallel processing\n", ncores))

# ----------------------------------------
# Output Directories
# ----------------------------------------
dir_sarimax <- "./Results_SARIMAX"
dir_rf <- "./Results_RF"
dir_lgbm <- "./Results_LGBM"
dir_ets <- "./Results_ETS"
dir_ets_author <- "./Results_ETS_Author"
dir_sarimax_nwp <- "./Results_SARIMAX_NWP"
dir_ctwlsv_sarimax <- "./results_ctwlsv_sarimax"
dir_ctwlsv_rf <- "./results_ctwlsv_rf"
dir_ctwlsv_lgbm <- "./results_ctwlsv_lgbm"
dir_ctwlsv_ets <- "./results_ctwlsv_ets"
dir_ctbu_sarimax <- "./results_ctbu_sarimax"
dir_ctbu_rf <- "./results_ctbu_rf"
dir_ctbu_lgbm <- "./results_ctbu_lgbm"
dir_ctbu_ets <- "./results_ctbu_ets"
dir_ctwlsv_ets_author <- "./results_ctwlsv_ets_author"
dir_ctbu_ets_author <- "./results_ctbu_ets_author"
dir_ctwlsv_sarimax_nwp <- "./results_ctwlsv_sarimax_nwp"
dir_ctbu_sarimax_nwp <- "./results_ctbu_sarimax_nwp"
dir_output <- "./output"

# Create directories if they don't exist
dirs <- c(dir_sarimax, dir_rf, dir_lgbm, dir_ets,
          dir_ets_author, dir_sarimax_nwp,
          dir_ctwlsv_sarimax, dir_ctwlsv_rf, dir_ctwlsv_lgbm, dir_ctwlsv_ets,
          dir_ctbu_sarimax, dir_ctbu_rf, dir_ctbu_lgbm, dir_ctbu_ets,
          dir_ctwlsv_ets_author, dir_ctbu_ets_author,
          dir_ctwlsv_sarimax_nwp, dir_ctbu_sarimax_nwp,
          dir_output)
for (d in dirs) {
  if (!dir.exists(d)) {
    dir.create(d, recursive = TRUE)
    cat(sprintf("Created directory: %s\n", d))
  }
}

# ----------------------------------------
# Derived Parameters
# ----------------------------------------
# Total observations per day
obs_per_day <- m  # 24 hours

# Total forecast steps (hourly)
forecast_steps <- forecast.horizon * m  # 48 hours

# k* = sum of temporal aggregation factors
k.star <- sum(m / k.v)  # 1 + 2 + 3 + 4 + 6 + 8 + 12 + 24 = 60

# Number of forecast rows per horizon
# h=1: k* rows, h=2: additional forecasts
# Total rows in Y.hat depends on aggregation structure

cat("\n========================================\n")
cat("Configuration loaded successfully\n")
cat("========================================\n")
cat(sprintf("  Seasonal period (m): %d hours\n", m))
cat(sprintf("  Forecast horizon (h): %d days (%d hours)\n", h, forecast_steps))
cat(sprintf("  Training window: %d days\n", train.days))
cat(sprintf("  Replications: %d to %d\n", min(rep_range), max(rep_range)))
cat(sprintf("  Temporal levels: %s\n", paste(k.v, collapse = ", ")))
cat("========================================\n\n")
