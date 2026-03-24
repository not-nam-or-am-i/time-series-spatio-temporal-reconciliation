# Bug Fixes Documentation

This document describes critical bugs discovered and fixed in the RF_SARIMAX forecast reconciliation pipeline.

## Overview

Two major bugs were identified that caused incorrect reconciliation results:

1. **File Loading Order Bug** - Aggregated forecast files were loaded in wrong order
2. **FoReco Format Mismatch Bug** - Data format incompatibility with FoReco package

---

## Bug #1: File Loading Order

### Problem

When loading base forecasts for evaluation and reconciliation, aggregated files (Total, Region_01-05) were loaded in **alphabetical order** instead of matching the S matrix structure.

**Incorrect order (alphabetical):**
```
Region_01 -> column 1
Region_02 -> column 2
Region_03 -> column 3
Region_04 -> column 4
Region_05 -> column 5
Total     -> column 6
```

**Correct order (matching S matrix):**
```
Total     -> column 1 (318 stations)
Region_01 -> column 2 (27 stations)
Region_02 -> column 3 (73 stations)
Region_03 -> column 4 (101 stations)
Region_04 -> column 5 (86 stations)
Region_05 -> column 6 (31 stations)
```

### Impact

- Base forecast nRMSE showed unrealistic values (89%+ at L0/L1)
- Forecasts were compared against wrong actuals
- SARIMAX and RF appeared to perform terribly when they actually performed well

### Fix Applied

Added explicit file ordering in all scripts that load aggregated files:

```r
# Define the correct order for aggregated files
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
```

### Files Modified

- `3 - evaluate_forecasts.R` (lines 231-245)
- `1 - reco_ctwlsv.R` (lines 55-66)

---

## Bug #2: FoReco Format Mismatch

### Problem

Our forecasts are stored in **"day-major"** order, but the FoReco package expects **"k-major"** order for multi-day forecasts.

**Our format (day-major):**
```
[day1_k24, day1_k12_1, day1_k12_2, ..., day1_k1_1...24, day2_k24, day2_k12_1, ...]
```
- Day 1 temporal hierarchy (rows 1-60)
- Day 2 temporal hierarchy (rows 61-120)

**FoReco format (k-major):**
```
[day1_k24, day2_k24, day1_k12_1, day1_k12_2, day2_k12_1, day2_k12_2, ...]
```
- All k=24 values first (2 daily totals)
- All k=12 values next (4 half-day values)
- ... continuing to k=1 (48 hourly values)

### Impact

- Reconciled forecasts were **incoherent**: daily totals ≠ sum of hourly values
- FoReco was applying temporal constraints to wrong indices
- Reconciliation appeared to make forecasts worse instead of better
- nRMSE values were in thousands/hundreds of thousands

### Diagnosis

```r
# Before fix - incoherent output:
# Day 1 k=24 (row 1): 804.31
# Day 1 k=1 sum (rows 37-60): 1881.75
# Match: FALSE

# After fix - coherent output:
# Day 1 k=24 (row 1): 1105.22
# Day 1 k=1 sum (rows 37-60): 1105.22
# Match: TRUE
```

### Fix Applied

Created format transformation functions to convert between our format and FoReco format:

```r
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

# Usage:
k_order <- c(24, 12, 8, 6, 4, 3, 2, 1)
idx_forecast <- create_foreco_indices(m, h, k_order)        # For h=2 day forecasts
idx_residuals <- create_foreco_indices(m, train.days, k_order)  # For 14-day residuals

# Transform to FoReco format before calling octrec/thfrec
Yhat_foreco <- Yhat[idx_forecast$to_foreco, ]
E_foreco <- E[idx_residuals$to_foreco, ]

# Call FoReco function...

# Transform back to our format after reconciliation
reco_our <- reco_foreco[idx_forecast$from_foreco, ]
```

### Files Modified

- `1 - reco_ctwlsv.R` (lines 33-71, 135-162, 169-196)
- `1 - reco_ctbu_twlsv.R` (lines 36-74, 140-160)

---

## Verification

After fixes, results show expected behavior:

### Base Forecasts (3 replications test)
| Level | PERS nRMSE | SARIMAX nRMSE | Skill |
|-------|------------|---------------|-------|
| L0 (Total) | 28.27% | 16.96% | +40.0% |
| L1 (Regions) | 28.57% | 21.96% | +23.1% |
| L2 (Stations) | 38.23% | 27.02% | +29.3% |

### Reconciled Forecasts
- Temporal coherence: Daily = sum(Hourly) ✓
- Spatial coherence: Total = sum(Regions) = sum(Stations) ✓

---

## Re-running the Pipeline

After the fixes, re-run scripts in this order:

```bash
# 1. Re-run CTWLSV reconciliation
Rscript "1 - reco_ctwlsv.R"

# 2. Re-run CTBU reconciliation
Rscript "1 - reco_ctbu_twlsv.R"

# 3. Re-run evaluation
Rscript "3 - evaluate_forecasts.R"

# 4. Generate comparison
Rscript "4 - compare_methods.R"
```

**Note:** Base forecasts do NOT need to be re-run. The forecasts were saved correctly; only the loading/processing had bugs.

---

## Date

Fixes applied: January 23, 2026
