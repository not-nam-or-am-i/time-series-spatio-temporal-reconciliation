# Base Forecast Models: Parameters and Features

This document describes the four base forecast models used in the RF_SARIMAX pipeline, including their parameters, features, and forecasting strategies.

## Common Settings

All models share the following configuration (defined in `config.R`):

| Parameter | Value | Description |
|-----------|-------|-------------|
| Training window | 14 days (336 hours) | Rolling window per replication |
| Forecast horizon | 2 days (48 hours) | Hourly forecasts |
| Seasonal period (m) | 24 | Hours per day |
| Temporal levels (k) | 1, 2, 3, 4, 6, 8, 12, 24 | Aggregation factors |
| Replications | 350 | Rolling-origin evaluations |
| Series | 324 total | 1 Total + 5 Regions + 318 Stations |

## 1. SARIMAX

**Script:** `0 - base_forecasts_sarimax.R`
**Package:** `forecast::auto.arima()`

### Model Selection

SARIMAX uses `auto.arima()` to automatically select the best ARIMA(p,d,q)(P,D,Q)[24] model for each series and replication independently. The model order varies across series depending on the data characteristics.

### Search Space

| Parameter | Range | Description |
|-----------|-------|-------------|
| max.p | 0-3 | Non-seasonal AR order |
| max.d | 0-1 | Non-seasonal differencing |
| max.q | 0-3 | Non-seasonal MA order |
| max.P | 0-2 | Seasonal AR order |
| max.D | 0-1 | Seasonal differencing |
| max.Q | 0-2 | Seasonal MA order |
| Seasonal period | 24 | Hourly seasonality |

### Fitting Options

| Option | Value | Description |
|--------|-------|-------------|
| stepwise | TRUE | Stepwise search (faster than exhaustive) |
| approximation | TRUE | Approximate likelihood for speed |
| allowdrift | FALSE | No drift term |
| seasonal | TRUE | Enable seasonal component |

### Exogenous Variable

- **NWP predictions** are used as the sole exogenous regressor (`xreg`)
- Training: `auto.arima(train_ts, xreg = train_nwp, ...)`
- Forecasting: `forecast(fit, h = 48, xreg = test_nwp)`

### Fallback Cascade

If model fitting fails, the following fallback order is used:

1. `auto.arima()` with NWP as `xreg`
2. `auto.arima()` without `xreg` (pure SARIMA)
3. Fixed ARIMA(1,0,1)(1,0,1)[24]

### Forecasting Strategy

Direct multi-step: `forecast(fit, h = 48)` produces all 48 hourly forecasts at once.

### Temporal Hierarchy Construction

Hourly forecasts are aggregated to all 8 temporal levels using `create_temporal_hierarchy()`, producing 120 values per replication (60 per forecast day).

---

## 2. Random Forest

**Script:** `0 - base_forecasts_rf.R`
**Package:** `ranger`

### Model Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| num.trees | 500 | Number of trees in the ensemble |
| mtry | floor(sqrt(12)) = 3 | Features sampled per split |
| min.node.size | 5 | Minimum observations per leaf |
| importance | none | Skipped for speed |
| seed | 42 + rp | Reproducible per replication |

### Features (12 total)

| # | Feature | Description |
|---|---------|-------------|
| 1 | `nwp` | NWP solar power prediction for the target hour |
| 2 | `lag_1` | Measured value 1 hour ago |
| 3 | `lag_24` | Measured value 24 hours ago (same hour yesterday) |
| 4 | `lag_48` | Measured value 48 hours ago (same hour 2 days ago) |
| 5 | `lag_168` | Measured value 168 hours ago (same hour last week) |
| 6 | `rolling_mean_24` | Rolling mean of last 24 hours |
| 7 | `rolling_sd_24` | Rolling standard deviation of last 24 hours |
| 8 | `hour_of_day` | Hour of day (0-23) |
| 9 | `day_of_week` | Day of week (0-6) |
| 10 | `month` | Month of year (1-12, approximated) |
| 11 | `hour_sin` | sin(2 * pi * hour / 24) — cyclic hour encoding |
| 12 | `hour_cos` | cos(2 * pi * hour / 24) — cyclic hour encoding |

### Forecasting Strategy

**Recursive multi-step forecasting:**

1. Train on 14-day window (336 observations, after removing NAs from lag features)
2. Predict hour t+1 using current features
3. Append prediction to history
4. Update lag and rolling features using appended history
5. Predict hour t+2, repeat until all 48 hours are forecast

This means forecast errors can accumulate, since each step uses its own predicted values as lag inputs.

### Temporal Hierarchy Construction

Same as SARIMAX: hourly forecasts aggregated to all 8 temporal levels.

---

## 3. LightGBM

**Script:** `0 - base_forecasts_lgbm.R`
**Package:** `lightgbm`

### Model Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| objective | regression | Squared loss |
| metric | rmse | Root mean squared error |
| learning_rate | 0.05 | Step size shrinkage |
| num_leaves | 31 | Maximum leaves per tree |
| min_data_in_leaf | 20 | Minimum observations per leaf |
| feature_fraction | 0.8 | Column subsampling ratio |
| bagging_fraction | 0.8 | Row subsampling ratio |
| bagging_freq | 5 | Bagging every 5 iterations |
| nrounds | 300 | Number of boosting rounds |
| verbose | -1 | Suppress output |

### Features (12 total)

Identical to Random Forest (see Section 2).

### Forecasting Strategy

Same recursive multi-step forecasting as Random Forest:
1. Train LightGBM model using `lgb.Dataset()` + `lgb.train()`
2. Predict recursively: `predict(lgb_model, as.matrix(features))`
3. Append predictions and update features for next step

### Temporal Hierarchy Construction

Same as RF and SARIMAX.

---

## 4. Exponential Smoothing (ETS)

**Script:** `0 - base_forecasts_ets.R`
**Packages:** `thief`, `forecast`

### Approach

ETS follows the original paper's methodology: independent ETS models are fitted at **each temporal aggregation level** using `thief:::th.forecast()`, rather than fitting a single hourly model and aggregating.

### Model Selection

`thief:::th.forecast(aggy, h, usemodel = "ets")` internally calls `forecast::ets()` at each aggregation level with automatic model selection (error type, trend type, seasonal type).

### Two Blocks

**Block 1: Aggregated series (L0 Total + L1 Regions) — Pure ETS**

1. Create temporal aggregates: `thief::tsaggregates(y, m = 24, align = "end")`
2. Fit independent ETS at each of 8 k levels: `thief:::th.forecast(yk, h = 48, usemodel = "ets")`
3. Collect forecasts and residuals from all levels

**Block 2: Bottom-level stations (L2, 318 series) — Hybrid ETS + NWP**

1. Same temporal aggregation and ETS fitting as Block 1
2. **Replace hourly (k=1) forecasts** with NWP predictions:
   ```
   frc$forecast[[1]] <- nwp_forecast
   ```
3. **Replace hourly (k=1) residuals** with NWP residuals:
   ```
   frc$residuals[[1]] <- nwp_train - y_train
   frc$mse[1] <- mean(frc$residuals[[1]]^2)
   ```
4. Higher aggregation levels (k=2 through k=24) retain their ETS forecasts

### Why Hybrid at L2?

NWP predictions are only available at the station level (L2), not for aggregated series (L0, L1). At the hourly resolution (k=1), NWP provides better forecasts than pure ETS for individual stations. At higher aggregation levels (k>=2), ETS captures the smoothed temporal patterns effectively.

### Temporal Hierarchy Construction

Forecasts from `thief:::th.forecast()` are already at all temporal levels. They are assembled into the 120-element Y.hat vector in day-major format using `extract_yhat_from_thief()`.

Residuals from each level are converted to the 840-element flat vector using `convert_thief_residuals()`.

### Fallback

If `thief:::th.forecast()` fails for a series, the script falls back to `forecast::ets()` on hourly data with manual temporal aggregation via `create_temporal_hierarchy()`.

---

## Summary Comparison

| Aspect | SARIMAX | Random Forest | LightGBM | ETS |
|--------|---------|---------------|----------|-----|
| Model type | Linear (ARIMA) | Ensemble (bagging) | Ensemble (boosting) | State space |
| NWP usage | Exogenous regressor | Input feature | Input feature | Replaces k=1 at L2 |
| Other features | None | 11 engineered | 11 engineered | None |
| Model selection | auto.arima() per series | Fixed hyperparameters | Fixed hyperparameters | ets() per level per series |
| Forecast strategy | Direct (h=48) | Recursive (step-by-step) | Recursive (step-by-step) | Direct (per level) |
| Temporal levels | Aggregated post-hoc | Aggregated post-hoc | Aggregated post-hoc | Fitted independently |
| Package | forecast | ranger | lightgbm | thief + forecast |
