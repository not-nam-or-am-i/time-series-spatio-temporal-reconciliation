# 0-info_reco.R
# ========================================
# Create Hierarchy Information for Reconciliation
# ========================================
# This script creates the cross-sectional (spatial) and temporal
# hierarchy matrices needed for forecast reconciliation.

rm(list = ls())

# Load required packages
library(hts)
library(FoReco)
library(Matrix)
library(data.table)

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
cat("Creating Hierarchy Information\n")
cat("========================================\n\n")

# ----------------------------------------
# Load measurement data to get station names
# ----------------------------------------
meas <- fread(DATA_PATH)

# Get station names (exclude timestamp column if present)
if (is.character(meas[[1]]) || inherits(meas[[1]], "POSIXct")) {
  station_names <- colnames(meas)[-1]
} else {
  station_names <- colnames(meas)
}

n_bottom <- length(station_names)
cat(sprintf("Number of bottom-level series: %d\n", n_bottom))

# ----------------------------------------
# Cross-sectional (Spatial) Hierarchy
# ----------------------------------------
# The hierarchy structure:
# - Level 0 (L0): Total (1 series) - sum of all stations
# - Level 1 (L1): 5 Regions (aggregates)
# - Level 2 (L2): 318 Bottom-level stations

# Station naming convention: "XX|YYY" where XX is region (01-05)
# Extract region codes from station names
regions <- substr(station_names, 1, 2)
unique_regions <- sort(unique(regions))

cat(sprintf("Number of regions: %d\n", length(unique_regions)))
cat("Region distribution:\n")

# Count stations per region
level2 <- sapply(unique_regions, function(r) sum(regions == r))
names(level2) <- unique_regions
print(level2)

# Create cross-sectional aggregation matrix C (na x nb)
# na = number of upper levels = 6 (Total + 5 regions)
# nb = number of bottom levels = 318

# Row 1: Total = sum of all 318 stations
# Rows 2-6: Each region sums its respective stations

C <- rbind(
  rep(1, n_bottom),  # Total row
  Matrix::bdiag(lapply(level2, function(x) t(rep(1, x))))  # Region rows
)

cat(sprintf("\nC matrix dimensions: %d x %d\n", nrow(C), ncol(C)))

# Create hts_info using FoReco
hts_info <- FoReco::hts_tools(C = C)

cat("hts_info created with components:\n")
print(names(hts_info))
cat(sprintf("  S matrix: %d x %d\n", nrow(hts_info$S), ncol(hts_info$S)))

# ----------------------------------------
# Temporal Hierarchy
# ----------------------------------------
# m = 24: Seasonal period (hours per day)
# h = 2: Forecast horizon (days)
# Aggregation levels k: {1, 2, 3, 4, 6, 8, 12, 24}

thf_info <- FoReco::thf_tools(m = m, h = h)

cat("\nthf_info created with components:\n")
print(names(thf_info))
cat(sprintf("  Seasonal period (m): %d\n", m))
cat(sprintf("  Forecast horizon (h): %d\n", h))
cat(sprintf("  Temporal aggregation levels: %s\n", paste(thf_info$kset, collapse = ", ")))

# ----------------------------------------
# Variable names for all series
# ----------------------------------------
# Order: Total, Region_01-05, then 318 bottom-level stations
varname <- c(
  "Total",
  paste0("Region_", unique_regions),
  station_names
)

cat(sprintf("\nTotal series in hierarchy: %d\n", length(varname)))
cat(sprintf("  L0 (Total): 1\n"))
cat(sprintf("  L1 (Regions): %d\n", length(unique_regions)))
cat(sprintf("  L2 (Stations): %d\n", n_bottom))

# ----------------------------------------
# Save hierarchy information
# ----------------------------------------
save(hts_info, thf_info, varname, m, h, k.v, train.days, level2,
     file = "info_reco.RData")

cat("\n========================================\n")
cat("Hierarchy info saved to info_reco.RData\n")
cat("========================================\n")
cat("\nContents:\n")
cat("  - hts_info: Cross-sectional hierarchy tools\n")
cat("  - thf_info: Temporal hierarchy tools\n")
cat("  - varname: Names of all 324 series\n")
cat("  - m, h: Temporal parameters\n")
cat("  - k.v: Aggregation levels\n")
cat("  - train.days: Training window\n")
cat("  - level2: Station counts per region\n")
