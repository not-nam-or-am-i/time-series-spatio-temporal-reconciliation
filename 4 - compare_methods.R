# 4 - compare_methods.R
# ========================================
# Method Comparison and Visualization
# Generates comparison tables and summary reports
# ========================================

rm(list = ls())

# Load required packages
libs <- c("data.table", "dplyr", "tidyr")
invisible(lapply(libs, library, character.only = TRUE))

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
cat("Method Comparison and Analysis\n")
cat("========================================\n\n")

# ----------------------------------------
# Load evaluation results
# ----------------------------------------
if (!file.exists("output/evaluation_results.RData")) {
  stop("Evaluation results not found. Run '3 - evaluate_forecasts.R' first.")
}

load("output/evaluation_results.RData")
summary_stats <- fread("output/summary_all_methods.csv")

cat(sprintf("Loaded %d evaluation records\n", nrow(all_results)))
cat(sprintf("Methods: %s\n", paste(unique(all_results$method), collapse = ", ")))

# ----------------------------------------
# Generate Comparison Tables
# ----------------------------------------

cat("\n----------------------------------------\n")
cat("Generating comparison tables...\n")
cat("----------------------------------------\n")

# Desired column order after pivot
col_order <- c("method", "L0_Hourly", "L0_Daily",
               "L1_Hourly", "L1_Daily",
               "L2_Hourly", "L2_Daily")

# Table 1: nRMSE by method and level
nrmse_table <- summary_stats %>%
  select(method, level, freq, mean_nRMSE) %>%
  mutate(mean_nRMSE = round(mean_nRMSE, 2)) %>%
  pivot_wider(names_from = c(level, freq), values_from = mean_nRMSE,
              names_sep = "_") %>%
  select(any_of(col_order))

fwrite(nrmse_table, "output/comparison_nRMSE.csv")
cat("  - comparison_nRMSE.csv\n")

# Table 2: nMBE by method and level
nmbe_table <- summary_stats %>%
  select(method, level, freq, mean_nMBE) %>%
  mutate(mean_nMBE = round(mean_nMBE, 2)) %>%
  pivot_wider(names_from = c(level, freq), values_from = mean_nMBE,
              names_sep = "_") %>%
  select(any_of(col_order))

fwrite(nmbe_table, "output/comparison_nMBE.csv")
cat("  - comparison_nMBE.csv\n")

# Table 3: Skill score by method and level
skill_table <- summary_stats %>%
  select(method, level, freq, mean_skill) %>%
  mutate(mean_skill = round(mean_skill, 2)) %>%
  pivot_wider(names_from = c(level, freq), values_from = mean_skill,
              names_sep = "_") %>%
  select(any_of(col_order))

fwrite(skill_table, "output/comparison_skill.csv")
cat("  - comparison_skill.csv\n")

# Table 4: SMAPE by method and level
smape_table <- summary_stats %>%
  select(method, level, freq, mean_SMAPE) %>%
  mutate(mean_SMAPE = round(mean_SMAPE, 2)) %>%
  pivot_wider(names_from = c(level, freq), values_from = mean_SMAPE,
              names_sep = "_") %>%
  select(any_of(col_order))

fwrite(smape_table, "output/comparison_SMAPE.csv")
cat("  - comparison_SMAPE.csv\n")

# ----------------------------------------
# Generate Summary Report
# ----------------------------------------

cat("\n========================================\n")
cat("FORECAST COMPARISON SUMMARY\n")
cat("========================================\n\n")

cat("Methods compared:\n")
cat("  - PERS: Persistence baseline (previous day's values)\n")
cat("  - SARIMAX_BASE: SARIMAX with NWP as exogenous variable\n")
cat("  - RF_BASE: Random Forest (NWP + lag + calendar features)\n")
cat("  - LGBM_BASE: LightGBM (NWP + lag + calendar features)\n")
cat("  - ETS_BASE: ETS (per-level via thief, NWP hybrid at L2)\n")
cat("  - ETS_AUTHOR_BASE: ETS (original author code, no clipping)\n")
cat("  - SARIMAX_NWP_BASE: SARIMAX with NWP hourly replacement at L2 + SNTZ\n")
cat("  - CTWLSV_SARIMAX: Cross-temporal WLSV reconciled SARIMAX\n")
cat("  - CTWLSV_RF: Cross-temporal WLSV reconciled RF\n")
cat("  - CTWLSV_LGBM: Cross-temporal WLSV reconciled LightGBM\n")
cat("  - CTWLSV_ETS: Cross-temporal WLSV reconciled ETS\n")
cat("  - CTWLSV_ETS_AUTHOR: Cross-temporal WLSV reconciled ETS (author)\n")
cat("  - CTWLSV_SARIMAX_NWP: Cross-temporal WLSV reconciled SARIMAX+NWP\n")
cat("  - CTBU_SARIMAX: Bottom-up temporal reconciled SARIMAX\n")
cat("  - CTBU_RF: Bottom-up temporal reconciled RF\n")
cat("  - CTBU_LGBM: Bottom-up temporal reconciled LightGBM\n")
cat("  - CTBU_ETS: Bottom-up temporal reconciled ETS\n")
cat("  - CTBU_ETS_AUTHOR: Bottom-up temporal reconciled ETS (author)\n")
cat("  - CTBU_SARIMAX_NWP: Bottom-up temporal reconciled SARIMAX+NWP\n\n")

# ----------------------------------------
# Print nRMSE Table
# ----------------------------------------
cat("----------------------------------------\n")
cat("Table 1: Mean nRMSE (%) by Method and Level\n")
cat("----------------------------------------\n\n")

# Format for printing
nrmse_print <- nrmse_table
colnames(nrmse_print) <- gsub("_", " ", colnames(nrmse_print))
print(as.data.frame(nrmse_print), row.names = FALSE)

# ----------------------------------------
# Print Skill Score Table
# ----------------------------------------
cat("\n----------------------------------------\n")
cat("Table 2: Mean Skill Score (%) vs PERS Baseline\n")
cat("(Positive = better than PERS)\n")
cat("----------------------------------------\n\n")

skill_print <- skill_table
colnames(skill_print) <- gsub("_", " ", colnames(skill_print))
print(as.data.frame(skill_print), row.names = FALSE)

# ----------------------------------------
# Print nMBE Table
# ----------------------------------------
cat("\n----------------------------------------\n")
cat("Table 3: Mean nMBE (%) by Method and Level\n")
cat("(Positive = over-forecasting, Negative = under-forecasting)\n")
cat("----------------------------------------\n\n")

nmbe_print <- nmbe_table
colnames(nmbe_print) <- gsub("_", " ", colnames(nmbe_print))
print(as.data.frame(nmbe_print), row.names = FALSE)

# ----------------------------------------
# Print SMAPE Table
# ----------------------------------------
cat("\n----------------------------------------\n")
cat("Table 4: Mean SMAPE (%) by Method and Level\n")
cat("----------------------------------------\n\n")

smape_print <- smape_table
colnames(smape_print) <- gsub("_", " ", colnames(smape_print))
print(as.data.frame(smape_print), row.names = FALSE)

# ----------------------------------------
# Best Method Analysis
# ----------------------------------------
cat("\n========================================\n")
cat("BEST METHOD ANALYSIS\n")
cat("========================================\n\n")

# Find best method for each level/freq combination (by nRMSE)
best_by_nrmse <- summary_stats %>%
  group_by(level, freq) %>%
  slice_min(mean_nRMSE, n = 1) %>%
  select(level, freq, method, mean_nRMSE) %>%
  ungroup()

cat("Best method by nRMSE:\n")
for (i in 1:nrow(best_by_nrmse)) {
  row <- best_by_nrmse[i, ]
  cat(sprintf("  %s %s: %s (%.2f%%)\n",
              row$level, row$freq, row$method, row$mean_nRMSE))
}

# Find best method for each level/freq combination (by Skill)
best_by_skill <- summary_stats %>%
  filter(method != "pers") %>%  # Exclude PERS (it's the baseline)
  group_by(level, freq) %>%
  slice_max(mean_skill, n = 1) %>%
  select(level, freq, method, mean_skill) %>%
  ungroup()

cat("\nBest method by Skill Score:\n")
for (i in 1:nrow(best_by_skill)) {
  row <- best_by_skill[i, ]
  cat(sprintf("  %s %s: %s (%.2f%%)\n",
              row$level, row$freq, row$method, row$mean_skill))
}

# ----------------------------------------
# Reconciliation Impact Analysis
# ----------------------------------------
cat("\n========================================\n")
cat("RECONCILIATION IMPACT ANALYSIS\n")
cat("========================================\n\n")

# Compare base vs reconciled for SARIMAX
sarimax_comparison <- summary_stats %>%
  filter(method %in% c("sarimax_base", "ctwlsv_sarimax", "ctbu_sarimax")) %>%
  select(method, level, freq, mean_nRMSE) %>%
  pivot_wider(names_from = method, values_from = mean_nRMSE)

cat("SARIMAX: Base vs Reconciled nRMSE:\n")
print(as.data.frame(sarimax_comparison), row.names = FALSE)

# Compare base vs reconciled for RF
rf_comparison <- summary_stats %>%
  filter(method %in% c("rf_base", "ctwlsv_rf", "ctbu_rf")) %>%
  select(method, level, freq, mean_nRMSE) %>%
  pivot_wider(names_from = method, values_from = mean_nRMSE)

cat("\nRandom Forest: Base vs Reconciled nRMSE:\n")
print(as.data.frame(rf_comparison), row.names = FALSE)

# Compare base vs reconciled for LightGBM
lgbm_comparison <- summary_stats %>%
  filter(method %in% c("lgbm_base", "ctwlsv_lgbm", "ctbu_lgbm")) %>%
  select(method, level, freq, mean_nRMSE) %>%
  pivot_wider(names_from = method, values_from = mean_nRMSE)

cat("\nLightGBM: Base vs Reconciled nRMSE:\n")
print(as.data.frame(lgbm_comparison), row.names = FALSE)

# Compare base vs reconciled for ETS
ets_comparison <- summary_stats %>%
  filter(method %in% c("ets_base", "ctwlsv_ets", "ctbu_ets")) %>%
  select(method, level, freq, mean_nRMSE) %>%
  pivot_wider(names_from = method, values_from = mean_nRMSE)

cat("\nETS: Base vs Reconciled nRMSE:\n")
print(as.data.frame(ets_comparison), row.names = FALSE)

# Compare base vs reconciled for ETS Author
ets_author_comparison <- summary_stats %>%
  filter(method %in% c("ets_author_base", "ctwlsv_ets_author",
                        "ctbu_ets_author")) %>%
  select(method, level, freq, mean_nRMSE) %>%
  pivot_wider(names_from = method, values_from = mean_nRMSE)

cat("\nETS (Author): Base vs Reconciled nRMSE:\n")
print(as.data.frame(ets_author_comparison), row.names = FALSE)

# Compare our ETS vs author ETS
ets_vs_author <- summary_stats %>%
  filter(method %in% c("ets_base", "ets_author_base")) %>%
  select(method, level, freq, mean_nRMSE) %>%
  pivot_wider(names_from = method, values_from = mean_nRMSE)

cat("\nETS Ours vs Author (base nRMSE):\n")
print(as.data.frame(ets_vs_author), row.names = FALSE)

# Compare base vs reconciled for SARIMAX+NWP
sarimax_nwp_comparison <- summary_stats %>%
  filter(method %in% c("sarimax_nwp_base", "ctwlsv_sarimax_nwp",
                        "ctbu_sarimax_nwp")) %>%
  select(method, level, freq, mean_nRMSE) %>%
  pivot_wider(names_from = method, values_from = mean_nRMSE)

cat("\nSARIMAX+NWP: Base vs Reconciled nRMSE:\n")
print(as.data.frame(sarimax_nwp_comparison), row.names = FALSE)

# Compare SARIMAX vs SARIMAX+NWP
sarimax_vs_nwp <- summary_stats %>%
  filter(method %in% c("sarimax_base", "sarimax_nwp_base")) %>%
  select(method, level, freq, mean_nRMSE) %>%
  pivot_wider(names_from = method, values_from = mean_nRMSE)

cat("\nSARIMAX vs SARIMAX+NWP (base nRMSE):\n")
print(as.data.frame(sarimax_vs_nwp), row.names = FALSE)

# ----------------------------------------
# Save Summary Report as Text
# ----------------------------------------
report_file <- "output/comparison_report.txt"

sink(report_file)
cat("========================================\n")
cat("FORECAST COMPARISON REPORT\n")
cat("SARIMAX, Random Forest, LightGBM, and ETS\n")
cat("========================================\n\n")

cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("METHODS COMPARED:\n")
cat("-----------------\n")
cat("1. PERS - Persistence baseline\n")
cat("2. SARIMAX_BASE - SARIMAX with NWP exogenous variable\n")
cat("3. RF_BASE - Random Forest (NWP + lag + calendar features)\n")
cat("4. LGBM_BASE - LightGBM (NWP + lag + calendar features)\n")
cat("5. ETS_BASE - ETS (per-level via thief, NWP hybrid at L2)\n")
cat("6. ETS_AUTHOR_BASE - ETS (original author code, no clipping)\n")
cat("7. SARIMAX_NWP_BASE - SARIMAX with NWP hourly replacement at L2 + SNTZ\n")
cat("8. CTWLSV_SARIMAX - Cross-temporal WLSV reconciled SARIMAX\n")
cat("9. CTWLSV_RF - Cross-temporal WLSV reconciled RF\n")
cat("10. CTWLSV_LGBM - Cross-temporal WLSV reconciled LightGBM\n")
cat("11. CTWLSV_ETS - Cross-temporal WLSV reconciled ETS\n")
cat("12. CTWLSV_ETS_AUTHOR - Cross-temporal WLSV reconciled ETS (author)\n")
cat("13. CTWLSV_SARIMAX_NWP - Cross-temporal WLSV reconciled SARIMAX+NWP\n")
cat("14. CTBU_SARIMAX - Bottom-up temporal reconciled SARIMAX\n")
cat("15. CTBU_RF - Bottom-up temporal reconciled RF\n")
cat("16. CTBU_LGBM - Bottom-up temporal reconciled LightGBM\n")
cat("17. CTBU_ETS - Bottom-up temporal reconciled ETS\n")
cat("18. CTBU_ETS_AUTHOR - Bottom-up temporal reconciled ETS (author)\n")
cat("19. CTBU_SARIMAX_NWP - Bottom-up temporal reconciled SARIMAX+NWP\n\n")

cat("METRICS:\n")
cat("--------\n")
cat("- nRMSE: Normalized RMSE (%)\n")
cat("- nMBE: Normalized Mean Bias Error (%)\n")
cat("- Skill: Improvement over PERS (%)\n")
cat("- SMAPE: Symmetric MAPE (%)\n\n")

cat("HIERARCHY LEVELS:\n")
cat("-----------------\n")
cat("- L0: Total (1 series)\n")
cat("- L1: Regional aggregates (5 series)\n")
cat("- L2: Bottom-level stations (318 series)\n\n")

cat("TEMPORAL RESOLUTION:\n")
cat("--------------------\n")
cat("- Hourly: k=1 (48 forecasts per replication)\n")
cat("- Daily: k=24 (2 forecasts per replication)\n\n")

cat("========================================\n")
cat("RESULTS TABLES\n")
cat("========================================\n\n")

cat("Table 1: Mean nRMSE (%)\n")
cat("-----------------------\n")
print(as.data.frame(nrmse_table), row.names = FALSE)

cat("\n\nTable 2: Mean Skill Score (%)\n")
cat("-----------------------------\n")
print(as.data.frame(skill_table), row.names = FALSE)

cat("\n\nTable 3: Mean nMBE (%)\n")
cat("----------------------\n")
print(as.data.frame(nmbe_table), row.names = FALSE)

cat("\n\nTable 4: Mean SMAPE (%)\n")
cat("-----------------------\n")
print(as.data.frame(smape_table), row.names = FALSE)

cat("\n\n========================================\n")
cat("KEY FINDINGS\n")
cat("========================================\n\n")

# Auto-generate findings
cat("Best performing methods:\n")
for (i in 1:nrow(best_by_nrmse)) {
  row <- best_by_nrmse[i, ]
  cat(sprintf("  - %s %s: %s (nRMSE: %.2f%%)\n",
              row$level, row$freq, row$method, row$mean_nRMSE))
}

sink()

cat(sprintf("\nReport saved to: %s\n", report_file))

# ----------------------------------------
# Final Summary
# ----------------------------------------
cat("\n========================================\n")
cat("COMPARISON COMPLETE\n")
cat("========================================\n\n")

cat("Output files generated:\n")
cat("  - output/comparison_nRMSE.csv\n")
cat("  - output/comparison_nMBE.csv\n")
cat("  - output/comparison_skill.csv\n")
cat("  - output/comparison_SMAPE.csv\n")
cat("  - output/comparison_report.txt\n")
cat("  - output/evaluation_results.RData (full results)\n")
cat("  - output/summary_all_methods.csv\n")
