# Temporal Hierarchy Data Format

This document explains the temporal hierarchy data structure used in this project and its compatibility with the FoReco package.

## Configuration

```r
m = 24           # Seasonal period (hours per day)
h = 2            # Forecast horizon (days)
train.days = 14  # Training window (days)
k_order = c(24, 12, 8, 6, 4, 3, 2, 1)  # Temporal aggregation levels
```

## Temporal Aggregation Levels

| k | Description | Periods per day | Period duration |
|---|-------------|-----------------|-----------------|
| 24 | Daily | 1 | 24 hours |
| 12 | Half-day | 2 | 12 hours |
| 8 | Tri-daily | 3 | 8 hours |
| 6 | Quarter-day | 4 | 6 hours |
| 4 | 6-hourly | 6 | 4 hours |
| 3 | 8-hourly | 8 | 3 hours |
| 2 | Bi-hourly | 12 | 2 hours |
| 1 | Hourly | 24 | 1 hour |

**k* (total periods per day):** 1 + 2 + 3 + 4 + 6 + 8 + 12 + 24 = **60**

## Our Data Format (Day-Major)

Forecasts are stored with each day's complete temporal hierarchy concatenated:

```
Vector structure (120 elements for h=2 days):

Day 1 (rows 1-60):
  Row 1:      k=24 (1 daily total)
  Rows 2-3:   k=12 (2 half-day values)
  Rows 4-6:   k=8  (3 tri-daily values)
  Rows 7-10:  k=6  (4 quarter-day values)
  Rows 11-16: k=4  (6 four-hourly values)
  Rows 17-24: k=3  (8 three-hourly values)
  Rows 25-36: k=2  (12 bi-hourly values)
  Rows 37-60: k=1  (24 hourly values)

Day 2 (rows 61-120):
  Row 61:      k=24 (1 daily total)
  Rows 62-63:  k=12 (2 half-day values)
  ... (same pattern as Day 1)
  Rows 97-120: k=1  (24 hourly values)
```

**Coherence Property:**
- `Y[1] = sum(Y[37:60])` (Day 1 daily = sum of Day 1 hourly)
- `Y[61] = sum(Y[97:120])` (Day 2 daily = sum of Day 2 hourly)

## FoReco Data Format (K-Major)

FoReco expects all periods for each aggregation level grouped together:

```
Vector structure (120 elements for h=2 days):

k=24 (rows 1-2):
  Row 1: Day 1 daily total
  Row 2: Day 2 daily total

k=12 (rows 3-6):
  Rows 3-4: Day 1 half-day values
  Rows 5-6: Day 2 half-day values

k=8 (rows 7-12):
  Rows 7-9:   Day 1 tri-daily values
  Rows 10-12: Day 2 tri-daily values

... (continuing pattern)

k=1 (rows 73-120):
  Rows 73-96:  Day 1 hourly values (24 values)
  Rows 97-120: Day 2 hourly values (24 values)
```

## Index Transformation

The `create_foreco_indices()` function creates mappings between formats:

```r
# Example for first 20 indices:
to_foreco = [1, 61, 2, 3, 62, 63, 4, 5, 6, 64, 65, 66, 7, 8, 9, 10, 67, 68, 69, 70, ...]

# Interpretation:
# FoReco position 1 <- Our position 1  (Day 1 k=24)
# FoReco position 2 <- Our position 61 (Day 2 k=24)
# FoReco position 3 <- Our position 2  (Day 1 k=12, period 1)
# FoReco position 4 <- Our position 3  (Day 1 k=12, period 2)
# FoReco position 5 <- Our position 62 (Day 2 k=12, period 1)
# FoReco position 6 <- Our position 63 (Day 2 k=12, period 2)
# ...
```

## Residuals Structure

Residuals follow the same pattern but for `train.days` instead of `h`:

```
n_residuals = k* × train.days = 60 × 14 = 840 elements
```

## Spatial Hierarchy

The S matrix defines spatial aggregation:

```
S dimensions: 324 × 318

Row structure:
  Row 1:     Total (aggregates all 318 bottom stations)
  Rows 2-6:  Regions (each aggregates a subset of stations)
  Rows 7-324: Bottom-level stations (identity mapping)

Column structure:
  Columns 1-318: Bottom-level stations

S[i,j] = 1 if bottom station j contributes to aggregate i
S[i,j] = 0 otherwise
```

## Full Forecast Matrix

Combined spatio-temporal forecast matrix:

```
Dimensions: 120 rows × 324 columns

Rows: Temporal hierarchy (k* × h = 60 × 2)
Columns: Spatial hierarchy (Total + 5 Regions + 318 Stations)

Column order MUST match S matrix:
  Column 1:     Total
  Columns 2-6:  Region_01 through Region_05
  Columns 7-324: Bottom-level stations (alphabetical order)
```

## Code Reference

### Loading with Correct Order

```r
# Aggregated files must be loaded in this order:
agg_names <- c("Total", "Region_01", "Region_02", "Region_03", "Region_04", "Region_05")

# Match files to names
for (i in seq_along(agg_names)) {
  pattern <- paste0(agg_names[i], "--0--")
  matching <- agg_files[grepl(pattern, agg_files, fixed = TRUE)]
  if (length(matching) > 0) {
    agg_files_ordered[i] <- matching[1]
  }
}
```

### Format Transformation

```r
# Create indices once at script start
idx_forecast <- create_foreco_indices(m, h, k_order)
idx_residuals <- create_foreco_indices(m, train.days, k_order)

# Transform to FoReco format
Yhat_foreco <- Yhat[idx_forecast$to_foreco, ]
E_foreco <- E[idx_residuals$to_foreco, ]

# Call FoReco
reco_raw <- FoReco::octrec(basef = t(Yhat_foreco), ...)

# Transform back to our format
reco_foreco <- t(reco_raw)
reco_our <- reco_foreco[idx_forecast$from_foreco, ]
```

## Verification

To verify coherence of reconciled forecasts:

```r
# Day 1 temporal coherence
stopifnot(abs(Y[1] - sum(Y[37:60])) < 1e-6)  # k=24 = sum(k=1)

# Day 2 temporal coherence
stopifnot(abs(Y[61] - sum(Y[97:120])) < 1e-6)

# Spatial coherence (for each temporal level)
stopifnot(abs(Y[1, 1] - sum(Y[1, 7:324])) < 1e-6)  # Total = sum(bottom)
```
