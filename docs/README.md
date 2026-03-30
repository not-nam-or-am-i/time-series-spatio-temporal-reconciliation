# RF_SARIMAX: Solar Forecast Reconciliation Pipeline

This project contains a self-contained implementation of four base forecast methods with cross-temporal and spatial reconciliation.

## Current Implementation

- **SARIMAX:** Seasonal ARIMA with NWP as exogenous regressor via `auto.arima()`
- **Random Forest:** 500 trees with 12 engineered features (NWP, lags, rolling stats, calendar)
- **LightGBM:** Gradient boosting with same 12 features as RF, 300 rounds, lr=0.05
- **ETS:** Exponential Smoothing at each temporal aggregation level via `thief`, with NWP hybrid at L2

## Documentation

- [README.md](README.md) - This file (project overview)
- [DATA_DESCRIPTION.md](DATA_DESCRIPTION.md) - Input data, hierarchy structure, evaluation setup
- [BASE_FORECAST_MODELS.md](BASE_FORECAST_MODELS.md) - All base forecast model parameters and features
- [RECONCILIATION_METHODS.md](RECONCILIATION_METHODS.md) - CTWLSV and CTBU reconciliation methods
- [BUG_FIXES.md](BUG_FIXES.md) - Critical bug fixes applied (file ordering, FoReco format)
- [TEMPORAL_HIERARCHY_FORMAT.md](TEMPORAL_HIERARCHY_FORMAT.md) - Technical reference for data formats

## Folder Structure

```
RF_SARIMAX/
├── data/
│   ├── meas.csv              # Measurement data (318 stations, hourly)
│   └── pred.csv              # NWP predictions (318 stations, hourly)
├── docs/
│   ├── README.md             # This file
│   ├── DATA_DESCRIPTION.md      # Input data and hierarchy structure
│   ├── BASE_FORECAST_MODELS.md  # Model parameters and features
│   ├── RECONCILIATION_METHODS.md # CTWLSV and CTBU reconciliation
│   ├── BUG_FIXES.md          # Bug fix documentation
│   └── TEMPORAL_HIERARCHY_FORMAT.md  # Data format reference
├── config.R                  # Shared configuration parameters
├── fun.R                     # Helper functions
├── 0-info_reco.R           # Create hierarchy information
├── 0-base_forecasts_sarimax.R   # SARIMAX forecast generation
├── 0-base_forecasts_rf.R        # Random Forest forecast generation
├── 0-base_forecasts_lgbm.R      # LightGBM forecast generation
├── 0-base_forecasts_ets.R       # ETS forecast generation
├── 1-reco_ctwlsv.R         # Cross-temporal WLSV reconciliation
├── 1-reco_ctbu_twlsv.R     # Bottom-up temporal reconciliation
├── 3-evaluate_forecasts.R  # Calculate evaluation metrics
├── 4-compare_methods.R     # Generate comparison tables
├── Results_SARIMAX/          # SARIMAX forecast outputs
├── Results_RF/               # Random Forest forecast outputs
├── Results_LGBM/             # LightGBM forecast outputs
├── Results_ETS/              # ETS forecast outputs
├── results_ctwlsv_sarimax/   # CTWLSV reconciled SARIMAX
├── results_ctwlsv_rf/        # CTWLSV reconciled RF
├── results_ctwlsv_lgbm/      # CTWLSV reconciled LightGBM
├── results_ctwlsv_ets/       # CTWLSV reconciled ETS
├── results_ctbu_sarimax/     # CTBU reconciled SARIMAX
├── results_ctbu_rf/          # CTBU reconciled RF
├── results_ctbu_lgbm/        # CTBU reconciled LightGBM
├── results_ctbu_ets/         # CTBU reconciled ETS
└── output/                   # Final comparison results (CSV)
```

## Execution Order

Run scripts in this order:

```bash
# 1. Setup hierarchy information (run once)
Rscript "0-info_reco.R"

# 2. Generate base forecasts (can run in parallel on different machines)
Rscript "0-base_forecasts_sarimax.R"   # ~2-4 hours
Rscript "0-base_forecasts_rf.R"        # ~1-2 hours
Rscript "0-base_forecasts_lgbm.R"      # ~1-2 hours
Rscript "0-base_forecasts_ets.R"       # ~2-5 hours

# 3. Run reconciliation (after base forecasts complete)
Rscript "1-reco_ctwlsv.R"              # ~30-60 min
Rscript "1-reco_ctbu_twlsv.R"          # ~15-30 min

# 4. Evaluate and compare
Rscript "3-evaluate_forecasts.R"       # ~10 min
Rscript "4-compare_methods.R"          # ~1 min
```

## Quick Test

To test with a small subset before full run:

1. Edit `config.R`
2. Change `rep_range <- 1:350` to `rep_range <- 1:3`
3. Run all scripts
4. Verify output format

## Methods Evaluated

| Method | Description |
|--------|-------------|
| PERS | Persistence baseline (previous day's values) |
| SARIMAX_BASE | SARIMAX with NWP as exogenous variable |
| RF_BASE | Random Forest with NWP + 11 engineered features |
| LGBM_BASE | LightGBM with NWP + 11 engineered features |
| ETS_BASE | ETS per temporal level (NWP hybrid at L2) |
| CTWLSV_SARIMAX | Cross-temporal WLSV reconciled SARIMAX |
| CTWLSV_RF | Cross-temporal WLSV reconciled RF |
| CTWLSV_LGBM | Cross-temporal WLSV reconciled LightGBM |
| CTWLSV_ETS | Cross-temporal WLSV reconciled ETS |
| CTBU_SARIMAX | Bottom-up temporal reconciled SARIMAX |
| CTBU_RF | Bottom-up temporal reconciled RF |
| CTBU_LGBM | Bottom-up temporal reconciled LightGBM |
| CTBU_ETS | Bottom-up temporal reconciled ETS |

## Metrics

- **nRMSE**: Normalized RMSE (%)
- **nMBE**: Normalized Mean Bias Error (%)
- **Skill**: Improvement over PERS baseline (%)
- **SMAPE**: Symmetric Mean Absolute Percentage Error (%)

## Hierarchy Structure

- **L0**: Total (1 series) - Sum of all stations
- **L1**: Regions (5 series) - Regional aggregates
- **L2**: Stations (318 series) - Bottom-level

## Temporal Aggregation

| k | Description | Periods/day | Forecasts (h=2) |
|---|-------------|-------------|-----------------|
| 1 | Hourly | 24 | 48 |
| 2 | Bi-hourly | 12 | 24 |
| 3 | 3-hourly | 8 | 16 |
| 4 | 4-hourly | 6 | 12 |
| 6 | 6-hourly | 4 | 8 |
| 8 | 8-hourly | 3 | 6 |
| 12 | 12-hourly | 2 | 4 |
| 24 | Daily | 1 | 2 |

**k*** (total periods per day): 1+2+3+4+6+8+12+24 = **60**

## Dependencies

```r
# Install required packages
install.packages(c(
  "FoReco",      # Forecast reconciliation
  "forecast",    # ARIMA/SARIMAX/ETS models
  "thief",       # Temporal hierarchical forecasting (ETS)
  "ranger",      # Fast Random Forest
  "lightgbm",    # LightGBM gradient boosting
  "hts",         # Hierarchical time series
  "doParallel",  # Parallel processing
  "doSNOW",      # Progress bars
  "Matrix",      # Sparse matrices
  "data.table",  # Fast data manipulation
  "dplyr",       # Data wrangling
  "tidyr"        # Data reshaping
))
```

## Output Files

After running all scripts, the `output/` folder will contain:

- `evaluation_results.RData` - Full evaluation results
- `summary_all_methods.csv` - Summary statistics
- `comparison_nRMSE.csv` - nRMSE comparison table
- `comparison_nMBE.csv` - nMBE comparison table
- `comparison_skill.csv` - Skill score comparison table
- `comparison_SMAPE.csv` - SMAPE comparison table
- `comparison_report.txt` - Text summary report

## Known Issues Fixed (January 2026)

Two critical bugs were discovered and fixed:

1. **File Loading Order Bug** - Aggregated forecast files were loaded alphabetically instead of matching the S matrix structure, causing forecasts to be compared against wrong actuals.

2. **FoReco Format Mismatch** - Data format incompatibility between our "day-major" storage format and FoReco's expected "k-major" format caused reconciled forecasts to be incoherent.

See [BUG_FIXES.md](BUG_FIXES.md) for full details.

**Important:** If reconciliation results show extremely high nRMSE values or incoherent forecasts (daily != sum of hourly), verify that:
1. Aggregated files are loaded in correct order: Total, Region_01, ..., Region_05
2. Format transformation indices are applied before calling FoReco functions
