# Data Description

This document describes the input data, hierarchy structure, and evaluation setup used in the RF_SARIMAX pipeline.

## Input Data

### Measurement Data (`data/meas.csv`)

| Property | Value |
|----------|-------|
| File | `data/meas.csv` |
| Rows | 8,760 (365 days x 24 hours) |
| Columns | 318 (one per PV station) |
| Resolution | Hourly |
| Unit | Solar power generation (kW or normalized) |
| Coverage | 1 full year |

Each column represents one bottom-level photovoltaic (PV) station. Station names follow the convention `XX|YYY`, where `XX` is the two-digit region code (01-05) and `YYY` identifies the station within that region.

### NWP Predictions (`data/pred.csv`)

| Property | Value |
|----------|-------|
| File | `data/pred.csv` |
| Rows | 8,760 (aligned with `meas.csv`) |
| Columns | 318 (one per PV station) |
| Resolution | Hourly |
| Unit | Predicted solar power (same unit as measurements) |

Numerical Weather Prediction (NWP) forecasts, aligned row-by-row and column-by-column with the measurement data. These serve as exogenous inputs for SARIMAX (as `xreg`), as a feature for Random Forest and LightGBM, and as hourly forecast replacements for the ETS hybrid approach at L2.

## Spatial Hierarchy

The 318 stations are organized into a three-level spatial hierarchy:

```
Level 0 (L0): Total
                |
    +-----------+-----------+
    |     |     |     |     |
Level 1 (L1): Region_01  Region_02  Region_03  Region_04  Region_05
                |
Level 2 (L2): Individual stations (318 total)
```

| Level | Name | Series | Description |
|-------|------|--------|-------------|
| L0 | Total | 1 | Sum of all 318 stations |
| L1 | Regions | 5 | Sum of stations within each region |
| L2 | Stations | 318 | Individual bottom-level PV stations |
| **Total** | | **324** | All series in the hierarchy |

### Aggregation Matrix (S)

The S matrix (324 x 318) maps bottom-level stations to all levels of the hierarchy:

```
S = [ C ]    where C (6 x 318) is the upper-level aggregation matrix
    [ I ]    and I (318 x 318) is the identity matrix
```

- **Row 1 of C** (Total): all ones — sums all 318 stations
- **Rows 2-6 of C** (Regions): block-diagonal — each region sums its member stations

The S matrix and related hierarchy tools are created by `0 - info_reco.R` and stored in `info_reco.RData`.

## Temporal Hierarchy

Hourly data (m = 24) is aggregated into 8 temporal levels using all divisors of 24:

| k | Name | Periods per day | Forecasts (h=2 days) | Aggregation |
|---|------|----------------|---------------------|-------------|
| 24 | Daily | 1 | 2 | Sum of 24 consecutive hours |
| 12 | 12-hourly | 2 | 4 | Sum of 12 consecutive hours |
| 8 | 8-hourly | 3 | 6 | Sum of 8 consecutive hours |
| 6 | 6-hourly | 4 | 8 | Sum of 6 consecutive hours |
| 4 | 4-hourly | 6 | 12 | Sum of 4 consecutive hours |
| 3 | 3-hourly | 8 | 16 | Sum of 3 consecutive hours |
| 2 | Bi-hourly | 12 | 24 | Sum of 2 consecutive hours |
| 1 | Hourly | 24 | 48 | No aggregation (base level) |

**k\*** = total periods per day across all levels = 1 + 2 + 3 + 4 + 6 + 8 + 12 + 24 = **60**

### Coherence Constraints

A coherent set of forecasts satisfies:
- **Temporal**: the daily forecast equals the sum of its 24 hourly forecasts (and similarly for all intermediate levels)
- **Spatial**: the Total forecast equals the sum of all 318 station forecasts; each Region forecast equals the sum of its member stations

## Rolling-Origin Evaluation

The pipeline uses a rolling-origin (expanding origin) evaluation design:

| Parameter | Value |
|-----------|-------|
| Training window | 14 days (336 hours) |
| Forecast horizon | 2 days (48 hours) |
| Total days | 365 |
| Replications | 350 (365 - 14 - 2 + 1) |

For each replication `rp` (1 to 350):
1. **Training period**: hours `(rp-1)*24 + 1` through `(rp-1)*24 + 336`
2. **Test period**: the next 48 hours immediately after training
3. The window shifts forward by 24 hours (1 day) for each replication

## Output Format

Each base forecast script produces one `.RData` file per series, containing a `results` list of length 350 (one entry per replication). Each entry stores:

| Field | Dimensions | Description |
|-------|-----------|-------------|
| `Y.hat` | 120 values | Forecasts at all 8 temporal levels for 2 days (60 per day) |
| `Y.obs` | 48 values | Hourly actual observations for the 2-day test period |
| `Res.insamp` | 840 values | In-sample residuals at all temporal levels for 14 training days (60 per day) |

### Y.hat Day-Major Format

The 120-element `Y.hat` vector is stored in **day-major** order:

```
[day1_k24, day1_k12, day1_k8, day1_k6, day1_k4, day1_k3, day1_k2, day1_k1,
 day2_k24, day2_k12, day2_k8, day2_k6, day2_k4, day2_k3, day2_k2, day2_k1]
```

Each day contributes 60 values (k* = 60). Within each day, values are ordered from the most aggregated level (k=24, 1 value) down to the hourly level (k=1, 24 values).

See [TEMPORAL_HIERARCHY_FORMAT.md](TEMPORAL_HIERARCHY_FORMAT.md) for the full technical specification.

## Evaluation Metrics

All metrics are computed by pooling forecasts across all 350 replications, then calculating a single metric value per (method, level, frequency) combination.

| Metric | Formula | Interpretation |
|--------|---------|---------------|
| nRMSE | `RMSE / mean(actual) * 100` | Normalized RMSE (%). Lower is better. |
| nMBE | `mean(forecast - actual) / mean(actual) * 100` | Bias (%). Positive = over-forecasting. |
| Skill | `(1 - RMSE_method / RMSE_pers) * 100` | Improvement over persistence (%). Positive = better than PERS. |
| SMAPE | `100 * mean(|actual - forecast| / ((|actual| + |forecast|) / 2))` | Symmetric MAPE (%). Lower is better. |

### Persistence Baseline (PERS)

The persistence forecast uses the last `h` days of the training period as the forecast for the next `h` days. For h=2, the forecast for the test period is simply the actual values from the last 2 days of training, assuming tomorrow will look like today.
