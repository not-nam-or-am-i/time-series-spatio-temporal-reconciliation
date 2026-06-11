# export_thesis_tables.R
# Export comparison CSVs to LaTeX tables for the thesis (17 rows x 6 columns).
#
# Usage (from repo root):
#   Rscript export_thesis_tables.R

libs <- c("data.table")
invisible(lapply(libs, library, character.only = TRUE))

script_dir <- "RF_SARIMAX"
if (!grepl("RF_SARIMAX$", getwd())) {
  if (dir.exists(script_dir)) {
    setwd(script_dir)
  }
}

out_dir <- "Vorlage_aus_Masterarbeit/tables"
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

col_order <- c("L0_Hourly", "L0_Daily", "L1_Hourly", "L1_Daily", "L2_Hourly", "L2_Daily")

# Row order: PERS, then each base family (CTWLSV then CTBU)
row_spec <- list(
  list(label = "PERS", method = "pers"),
  list(label = "SARIMAX", method = "ctwlsv_sarimax"),
  list(label = "SARIMAX", method = "ctbu_sarimax"),
  list(label = "SARIMAX+NWP", method = "ctwlsv_sarimax_nwp"),
  list(label = "SARIMAX+NWP", method = "ctbu_sarimax_nwp"),
  list(label = "RF", method = "ctwlsv_rf"),
  list(label = "RF", method = "ctbu_rf"),
  list(label = "RF+NWP", method = "ctwlsv_rf_nwp"),
  list(label = "RF+NWP", method = "ctbu_rf_nwp"),
  list(label = "LightGBM", method = "ctwlsv_lgbm"),
  list(label = "LightGBM", method = "ctbu_lgbm"),
  list(label = "LightGBM+NWP", method = "ctwlsv_lgbm_nwp"),
  list(label = "LightGBM+NWP", method = "ctbu_lgbm_nwp"),
  list(label = "ETS (author)", method = "ctwlsv_ets_author"),
  list(label = "ETS (author)", method = "ctbu_ets_author"),
  list(label = "ETS", method = "ctwlsv_ets"),
  list(label = "ETS", method = "ctbu_ets")
)

reco_label <- function(method) {
  if (method == "pers") return("---")
  if (startsWith(method, "ctwlsv_")) return("CTWLSV")
  if (startsWith(method, "ctbu_")) return("CTBU")
  method
}

export_metric <- function(csv_file, value_col, tex_file, caption, label) {
  if (!file.exists(csv_file)) {
    warning(sprintf("Missing %s — skip", csv_file))
    return(invisible(NULL))
  }

  dt <- fread(csv_file)
  dt[, method := tolower(method)]

  rows <- lapply(row_spec, function(spec) {
    row <- dt[method == tolower(spec$method)]
    if (nrow(row) == 0) {
      stop(sprintf("Method '%s' not found in %s", spec$method, csv_file))
    }
    vals <- as.numeric(row[1, ..col_order])
    data.table(
      base = spec$label,
      reconciliation = reco_label(spec$method),
      vals = list(vals)
    )
  })
  tab <- rbindlist(rows)

  lines <- c(
    "\\begin{table}[htbp]",
    "\\centering",
    "\\small",
    sprintf("\\caption{%s}", caption),
    sprintf("\\label{%s}", label),
    "\\begin{tabular}{llrrrrrr}",
    "\\toprule",
    "Base method & Reconciliation & L0 H & L0 D & L1 H & L1 D & L2 H & L2 D \\\\",
    "\\midrule"
  )

  for (i in seq_len(nrow(tab))) {
    vals <- tab$vals[[i]]
    vals_fmt <- sprintf("%.2f", vals)
    lines <- c(lines, sprintf(
      "%s & %s & %s \\\\",
      tab$base[i], tab$reconciliation[i], paste(vals_fmt, collapse = " & ")
    ))
    if (i %in% c(1, 5, 9, 13)) {
      lines <- c(lines, "\\addlinespace")
    }
  }

  lines <- c(lines, "\\bottomrule", "\\end{tabular}", "\\end{table}", "")
  writeLines(lines, tex_file)
  cat(sprintf("  Wrote %s\n", tex_file))
}

cat("Exporting thesis tables...\n")

metrics <- list(
  list(
    csv = "output/comparison_nRMSE.csv",
    tex = file.path(out_dir, "table_nrmse.tex"),
    caption = "Mean nRMSE (\\%) by base method, reconciliation approach, spatial level, and temporal resolution.",
    label = "tab:nrmse"
  ),
  list(
    csv = "output/comparison_nMBE.csv",
    tex = file.path(out_dir, "table_nmbe.tex"),
    caption = "Mean nMBE (\\%) by base method, reconciliation approach, spatial level, and temporal resolution.",
    label = "tab:nmbe"
  ),
  list(
    csv = "output/comparison_skill.csv",
    tex = file.path(out_dir, "table_skill.tex"),
    caption = "Mean Skill score (\\%) relative to PERS by base method, reconciliation approach, spatial level, and temporal resolution.",
    label = "tab:skill"
  ),
  list(
    csv = "output/comparison_SMAPE.csv",
    tex = file.path(out_dir, "table_smape.tex"),
    caption = "Mean SMAPE (\\%) by base method, reconciliation approach, spatial level, and temporal resolution.",
    label = "tab:smape"
  )
)

for (m in metrics) {
  export_metric(m$csv, NULL, m$tex, m$caption, m$label)
}

cat("Done.\n")
