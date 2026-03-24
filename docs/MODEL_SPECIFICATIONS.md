# Model Specifications: SARIMAX and Random Forest

This document describes the configuration and implementation details of the SARIMAX and Random Forest base forecast models used in the RF_SARIMAX pipeline.

## Shared Configuration

Both models share the following parameters (defined in `config.R`):

| Parameter | Value | Description |
|-----------|-------|-------------|
| `m` | 24 | Seasonal period (hours per day) |
| `h` | 2 | Forecast horizon (days) |
| `train.days` | 14 | Training window (days) |
| `forecast_steps` | 48 | Total hourly forecasts (h × m) |
| `k.v` | {1,2,3,4,6,8,12,24} | Temporal aggregation levels |
| `rep_range` | 1:350 | Number of train-test replications |

---

## SARIMAX Model

**Implementation File:** `0 - base_forecasts_sarimax.R`

### Model Type

SARIMAX (Seasonal ARIMA with eXogenous variables) using the `forecast` package's `auto.arima()` function.

### Model Specification

```
SARIMAX(p,d,q)(P,D,Q)[m]
```

Where:
- `(p,d,q)` = Non-seasonal ARIMA orders
- `(P,D,Q)` = Seasonal ARIMA orders
- `m = 24` = Seasonal period (hourly data, daily seasonality)

### Parameter Bounds

The `auto.arima()` function searches for optimal orders within these bounds:

| Parameter | Maximum Value | Description |
|-----------|---------------|-------------|
| `max.p` | 3 | Maximum non-seasonal AR order |
| `max.d` | 1 | Maximum non-seasonal differencing |
| `max.q` | 3 | Maximum non-seasonal MA order |
| `max.P` | 2 | Maximum seasonal AR order |
| `max.D` | 1 | Maximum seasonal differencing |
| `max.Q` | 2 | Maximum seasonal MA order |

### Model Selection Options

| Option | Value | Description |
|--------|-------|-------------|
| `seasonal` | TRUE | Include seasonal component |
| `stepwise` | TRUE | Use stepwise search (faster) |
| `approximation` | TRUE | Use approximation for speed |
| `allowdrift` | FALSE | Disable drift term |

### Exogenous Variable

**NWP (Numerical Weather Prediction)** is used as the exogenous regressor:

- **Training:** NWP predictions for the 14-day training window are passed via `xreg` parameter
- **Forecasting:** Future NWP predictions for the 48-hour forecast horizon are passed via `newxreg` parameter

```r
# Model fitting with NWP
fit <- auto.arima(train_ts, xreg = train_nwp, ...)

# Forecasting with future NWP
fc <- forecast(fit, h = 48, xreg = future_nwp)
```

### Fallback Strategy

If the SARIMAX model fails to fit (e.g., due to convergence issues), the following fallback cascade is applied:

1. **First fallback:** Try `auto.arima()` without exogenous variable
2. **Ultimate fallback:** Fit fixed ARIMA(1,0,1)(1,0,1)[24] model

### Output Format

For each replication and series:
- `Y.hat`: 120-element vector (temporal hierarchy for h=2 days)
- `Y.obs`: 48-element vector (hourly observations)
- `Res.insamp`: 840-element vector (residuals for 14 training days)

---

## Random Forest Model

**Implementation File:** `0 - base_forecasts_rf.R`

### Model Type

Random Forest regression using the `ranger` package (fast implementation of Random Forests).

### Hyperparameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `num.trees` | 500 | Number of trees in the forest |
| `mtry` | √(n_features) | Number of variables to sample at each split |
| `min.node.size` | 5 | Minimum node size |
| `importance` | "none" | Skip importance calculation (for speed) |
| `seed` | 42 + rp | Reproducibility seed (varies by replication) |

### Feature Engineering

Random Forest uses engineered features from historical data and NWP predictions:

#### 1. NWP Prediction (Exogenous)

| Feature | Description |
|---------|-------------|
| `nwp` | Numerical Weather Prediction value for target hour |

#### 2. Lag Features

| Feature | Description |
|---------|-------------|
| `lag_1` | Value 1 hour ago |
| `lag_24` | Value 24 hours ago (same hour yesterday) |
| `lag_48` | Value 48 hours ago (same hour 2 days ago) |
| `lag_168` | Value 168 hours ago (same hour 1 week ago) |

#### 3. Rolling Statistics

| Feature | Window | Description |
|---------|--------|-------------|
| `rolling_mean_24` | 24 hours | Mean of last 24 hourly values |
| `rolling_sd_24` | 24 hours | Standard deviation of last 24 hourly values |

#### 4. Calendar Features

| Feature | Description | Encoding |
|---------|-------------|----------|
| `hour_of_day` | Hour (0-23) | Raw integer |
| `day_of_week` | Day of week (0-6) | Raw integer |
| `month` | Month (1-12) | Raw integer |
| `hour_sin` | Sine of hour | sin(2π × hour/24) |
| `hour_cos` | Cosine of hour | cos(2π × hour/24) |

The sine/cosine encoding captures the cyclic nature of hours, allowing the model to understand that hour 23 is close to hour 0.

### Complete Feature List

```
Total features: 12

1. nwp           - NWP prediction (exogenous)
2. lag_1         - 1-hour lag
3. lag_24        - 24-hour lag (daily pattern)
4. lag_48        - 48-hour lag
5. lag_168       - 168-hour lag (weekly pattern)
6. rolling_mean_24 - 24-hour rolling mean
7. rolling_sd_24   - 24-hour rolling standard deviation
8. hour_of_day   - Hour of day (0-23)
9. day_of_week   - Day of week (0-6)
10. month        - Month (1-12)
11. hour_sin     - Cyclic encoding of hour (sine)
12. hour_cos     - Cyclic encoding of hour (cosine)
```

### Forecasting Strategy

**Recursive Multi-Step Forecasting:**

1. Train model on 14-day window with all features
2. For each forecast step (h=1 to h=48):
   - Create features using most recent data (including previous predictions)
   - Predict next hour's value
   - Append prediction to history
   - Update lag and rolling features for next step

```
Time: t  t+1  t+2  t+3  ...  t+48
Pred: ── p₁ → p₂ → p₃ → ... → p₄₈
              ↑
        Uses p₁ as lag_1
```

### Non-Negativity Constraint

All predictions are clipped to be non-negative:
```r
pred_val <- max(predict(rf_model, new_features)$predictions, 0)
```

### Output Format

Same as SARIMAX:
- `Y.hat`: 120-element vector (temporal hierarchy for h=2 days)
- `Y.obs`: 48-element vector (hourly observations)
- `Res.insamp`: 840-element vector (residuals for 14 training days)

---

## Comparison Summary

| Aspect | SARIMAX | Random Forest |
|--------|---------|---------------|
| **Model Type** | Parametric (ARIMA) | Non-parametric (Ensemble) |
| **Seasonality** | Explicit seasonal component | Captured via lag features |
| **NWP Usage** | Exogenous regressor (`xreg`) | Feature (among 12 total) |
| **Feature Engineering** | None (automatic) | Extensive (12 features) |
| **Model Selection** | Automatic (`auto.arima`) | Fixed hyperparameters |
| **Forecasting** | Direct (closed-form) | Recursive (multi-step) |
| **Package** | `forecast` | `ranger` |
| **Typical Runtime** | 2-4 hours (350 reps) | 1-2 hours (350 reps) |

---

## Temporal Hierarchy Construction

Both models produce hourly forecasts first, then aggregate to all temporal levels:

```
Hourly (k=1):  h₁, h₂, h₃, ..., h₄₈   (48 values)
                    ↓ sum
Daily (k=24):  D₁ = Σ(h₁..h₂₄), D₂ = Σ(h₂₅..h₄₈)   (2 values)
                    ↓ other k levels
...
```

The `create_temporal_hierarchy()` function aggregates hourly forecasts to all k levels following the day-major format documented in [TEMPORAL_HIERARCHY_FORMAT.md](TEMPORAL_HIERARCHY_FORMAT.md).

---

## Dependencies

```r
# SARIMAX
library(forecast)    # auto.arima, Arima, forecast

# Random Forest
library(ranger)      # Fast RF implementation

# Shared
library(doParallel)  # Parallel processing
library(doSNOW)      # Progress bars
library(data.table)  # Fast data I/O
library(Matrix)      # Sparse matrices
```

---

## References

- Hyndman, R.J., & Athanasopoulos, G. (2021). *Forecasting: Principles and Practice* (3rd ed). OTexts.
- Wright, M.N., & Ziegler, A. (2017). ranger: A Fast Implementation of Random Forests. *Journal of Statistical Software*, 77(1).
