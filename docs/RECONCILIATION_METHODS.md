# Reconciliation Methods

This document describes the two forecast reconciliation methods implemented in the RF_SARIMAX pipeline: CTWLSV (Cross-Temporal WLSV) and CTBU (Cross-Temporal Bottom-Up).

## Why Reconciliation?

Base forecasts are generated independently — each station and temporal level is forecast separately. This means the resulting forecasts are generally **incoherent**: the sum of hourly forecasts does not equal the daily forecast, and the sum of station forecasts does not equal the regional or total forecast. Reconciliation adjusts the base forecasts to restore coherence across both the spatial and temporal dimensions simultaneously.

## Method 1: CTWLSV — Cross-Temporal WLSV Reconciliation

**Script:** `1-reco_ctwlsv.R`
**Package:** `FoReco::octrec()`

### Overview

CTWLSV performs **joint optimal reconciliation** across both spatial and temporal dimensions simultaneously. It finds the reconciled forecasts that are closest to the base forecasts (in a weighted least squares sense) while satisfying all coherence constraints.

### Algorithm

For each replication:

1. **Load base forecasts**: all 324 series (6 aggregated + 318 bottom-level), each with 120 temporal values (Y.hat) and 840 residual values (Res.insamp)
2. **Format transformation**: convert from day-major storage format to FoReco's k-major format using pre-computed index mappings
3. **Reconciliation**: call `FoReco::octrec()` with the full cross-temporal system
4. **Format back**: convert reconciled forecasts back to day-major format

### FoReco Parameters

```r
FoReco::octrec(
  basef = t(Yhat_foreco),  # [324 x 120] base forecasts
  comb = "wlsv",           # Weighted Least Squares (variance scaling)
  m = 24,                  # Seasonal period
  C = hts_info$C,          # Cross-sectional aggregation matrix (6 x 318)
  res = t(E_foreco),       # [324 x 840] in-sample residuals
  keep = "recf",           # Return reconciled forecasts only
  type = "M"               # Matrix approach
)
```

| Parameter | Value | Description |
|-----------|-------|-------------|
| `comb` | `"wlsv"` | Weighted Least Squares with variance scaling. Uses residual variances to weight each series inversely by its forecast uncertainty. |
| `m` | 24 | Seasonal period (hours per day) |
| `C` | 6 x 318 matrix | Cross-sectional aggregation constraints |
| `res` | 324 x 840 matrix | In-sample residuals used to estimate the covariance structure |
| `type` | `"M"` | Matrix computation approach |

### WLSV Weighting

WLSV (Weighted Least Squares with Variance scaling) uses the **diagonal of the residual covariance matrix** as weights. This means:
- Series with larger residual variance (less accurate base forecasts) receive less weight during reconciliation
- Series with smaller residual variance (more accurate base forecasts) are adjusted less

This is a practical compromise between:
- **OLS** (unweighted, treats all series equally)
- **MinT-shr** (full covariance estimation, computationally expensive for 324 x 8 = 2,592 cross-temporal series)

### Non-Negativity Handling

Solar power cannot be negative. If the free (unconstrained) reconciled forecasts contain negative values, a second pass is run:

```r
FoReco::octrec(
  ...,            # same parameters
  sol = "osqp",   # Quadratic programming solver
  nn = TRUE       # Non-negativity constraint
)
```

This solves a constrained optimization problem ensuring all reconciled forecasts are non-negative while still satisfying coherence.

### Format Transformation

Base forecasts are stored in **day-major** order:
```
[day1_k24, day1_k12, ..., day1_k1, day2_k24, day2_k12, ..., day2_k1]
```

FoReco expects **k-major** order:
```
[all_days_k24, all_days_k12, ..., all_days_k1]
```

Pre-computed index mappings (`idx_forecast$to_foreco` and `idx_forecast$from_foreco`) handle the conversion in both directions. See [BUG_FIXES.md](BUG_FIXES.md) for the history of this issue.

---

## Method 2: CTBU — Cross-Temporal Bottom-Up with T-WLSV

**Script:** `1-reco_ctbu_twlsv.R`
**Package:** `FoReco::thfrec()` + manual spatial aggregation

### Overview

CTBU is a **two-step heuristic** approach:
1. **Temporal reconciliation** on bottom-level series only (using `thfrec`)
2. **Spatial bottom-up aggregation** using the S matrix

This is computationally cheaper than CTWLSV because it avoids the full cross-temporal optimization. It only reconciles the temporal dimension and then constructs upper-level forecasts by simple summation.

### Algorithm

For each replication:

1. **Load bottom-level forecasts only**: 318 series, each with 120 temporal values and 840 residual values
2. **Format transformation**: convert to FoReco k-major format
3. **Temporal reconciliation**: apply `thfrec()` independently to each of the 318 bottom-level series
4. **Format back**: convert reconciled bottom-level forecasts to day-major format
5. **Spatial aggregation**: multiply by S matrix to produce all 324 series

### Step 1: Temporal Reconciliation (thfrec)

```r
FoReco::thfrec(
  basef = Yhat_foreco[, j],  # 120-element forecast vector for series j
  comb = "wlsv",             # Weighted Least Squares (variance scaling)
  m = 24,                    # Seasonal period
  res = E_foreco[, j],       # 840-element residual vector for series j
  keep = "recf"              # Return reconciled forecasts only
)
```

This ensures that for each bottom-level station, the temporal aggregations are coherent (e.g., daily = sum of hourly).

### Step 2: Spatial Bottom-Up Aggregation

After temporal reconciliation, upper-level forecasts are constructed by summing:

```r
reco_full[row, ] <- S %*% thf_reco[row, ]
```

For each temporal row (120 rows total), the 318 bottom-level values are multiplied by the S matrix (324 x 318) to produce all 324 series. This is a pure bottom-up approach — upper-level base forecasts are discarded entirely.

### Non-Negativity

Since bottom-up aggregation only sums non-negative values, negative reconciled values are rare. If they occur (from the temporal reconciliation step), simple clipping is applied:

```r
reco_nn <- pmax(reco_free, 0)
```

---

## Comparison of Methods

| Aspect | CTWLSV | CTBU |
|--------|--------|------|
| Approach | Joint optimal | Two-step heuristic |
| Spatial handling | Optimal reconciliation | Bottom-up aggregation |
| Temporal handling | Optimal reconciliation | WLSV reconciliation |
| Uses upper-level base forecasts | Yes | No (discarded) |
| Coherence | Fully coherent | Fully coherent |
| Computation | ~30-60 min | ~15-30 min |
| Non-negativity | OSQP solver | Simple clipping |
| FoReco function | `octrec()` | `thfrec()` + S matrix |

### When CTWLSV is Better

- When upper-level base forecasts contain useful information (e.g., Total or Regional forecasts capture aggregate patterns that bottom-level models miss)
- When cross-series correlations matter

### When CTBU is Better

- When bottom-level forecasts are high quality and upper-level forecasts are noisy
- When computational resources are limited
- As a simpler, more interpretable baseline

## Applied To

Both reconciliation methods are applied to all four base forecast methods, producing 8 reconciled method variants:

| Base Method | CTWLSV | CTBU |
|------------|--------|------|
| SARIMAX | CTWLSV_SARIMAX | CTBU_SARIMAX |
| Random Forest | CTWLSV_RF | CTBU_RF |
| LightGBM | CTWLSV_LGBM | CTBU_LGBM |
| ETS | CTWLSV_ETS | CTBU_ETS |

Together with the 4 base methods and persistence baseline, this gives **13 methods** in total for comparison.
