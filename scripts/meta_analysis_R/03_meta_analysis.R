# Random-effects meta-analysis of gene fusion detection tools across 10 benchmarks.
#
# Runs three parallel meta-analyses per tool:
#   1. Precision  — TP / (TP + FP), binomial variance
#   2. Sensitivity — TP / (TP + FN), binomial variance (FN recovered where needed)
#   3. F1 — bootstrap variance to correctly account for the correlation
#           between precision and sensitivity through shared TP.
#
# Outputs: forest plots + pooled-estimate tables in results/{figures,tables}/.
#
# Written in base R (no tidyverse) to keep dependency footprint minimal.

suppressPackageStartupMessages(library(metafor))

# --- Paths --------------------------------------------------------------------
script_path <- if (!is.null(sys.frames()[[1]]$ofile)) {
  normalizePath(sys.frames()[[1]]$ofile)
} else {
  normalizePath(commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1])
}
SCRIPT_DIR  <- dirname(script_path)
REPO_ROOT   <- normalizePath(file.path(SCRIPT_DIR, "..", ".."), mustWork = FALSE)
INPUT_CSV   <- file.path(REPO_ROOT, "data", "processed",
                         "benchmarks_metrics_enriched.csv")
FIG_DIR     <- file.path(REPO_ROOT, "results", "figures")
TAB_DIR     <- file.path(REPO_ROOT, "results", "tables")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TAB_DIR, showWarnings = FALSE, recursive = TRUE)

MIN_BENCHMARKS_PER_TOOL <- 3
BOOT_ITER <- 5000

# --- Load and prepare ---------------------------------------------------------
df <- read.csv(INPUT_CSV, stringsAsFactors = FALSE, check.names = FALSE)

# Use recovered FN when the stored FN is missing.
df$FN_use <- ifelse(is.na(df$FN), df$FN_recovered, df$FN)
df$n_pred <- df$TP + df$FP
df$n_true <- df$TP + df$FN_use

# --- Bootstrap variance for F1 ------------------------------------------------
# For each row, resample TP/FP/FN from binomials with observed proportions,
# recompute F1, and return the empirical variance. This automatically accounts
# for the positive correlation between precision and sensitivity through the
# shared TP.
bootstrap_f1_var <- function(TP, FP, FN, n_iter = BOOT_ITER, seed = 42) {
  if (is.na(TP) || is.na(FP) || is.na(FN)) return(NA_real_)
  # rbinom requires integer size, so round counts. FN recovered from
  # Sensitivity may be fractional by 0.001-0.5 due to rounding of the
  # reported sensitivity value.
  TP <- round(TP); FP <- round(FP); FN <- round(FN)
  if (TP + FP <= 0 || TP + FN <= 0) return(NA_real_)
  set.seed(seed)
  n_pred <- TP + FP
  n_true <- TP + FN
  p_prec <- TP / n_pred
  p_sens <- TP / n_true

  f1s <- numeric(n_iter)
  for (i in seq_len(n_iter)) {
    tp_p <- rbinom(1, n_pred, p_prec)
    tp_s <- rbinom(1, n_true, p_sens)
    tp   <- min(tp_p, tp_s)
    prec <- tp / n_pred
    sens <- tp / n_true
    denom <- prec + sens
    f1s[i] <- if (is.na(denom) || denom == 0) 0 else 2 * prec * sens / denom
  }
  var(f1s)
}

# --- Compute effect sizes -----------------------------------------------------
# Precision: raw proportion; variance = p(1-p)/n
mask_p <- !is.na(df$TP) & !is.na(df$FP) & (df$n_pred > 0)
df$prec_yi <- ifelse(mask_p, df$TP / df$n_pred, NA_real_)
df$prec_vi <- ifelse(mask_p, df$prec_yi * (1 - df$prec_yi) / df$n_pred, NA_real_)

# Sensitivity: raw proportion; variance = p(1-p)/n
mask_s <- !is.na(df$TP) & !is.na(df$FN_use) & (df$n_true > 0)
df$sens_yi <- ifelse(mask_s, df$TP / df$n_true, NA_real_)
df$sens_vi <- ifelse(mask_s, df$sens_yi * (1 - df$sens_yi) / df$n_true, NA_real_)

# F1 with bootstrap variance
mask_f <- mask_p & mask_s
df$f1_yi <- rep(NA_real_, nrow(df))
df$f1_vi <- rep(NA_real_, nrow(df))
cat("Computing F1 bootstrap variances for", sum(mask_f), "rows...\n")
for (i in which(mask_f)) {
  p <- df$prec_yi[i]; s <- df$sens_yi[i]
  df$f1_yi[i] <- if (p + s > 0) 2 * p * s / (p + s) else 0
  df$f1_vi[i] <- bootstrap_f1_var(df$TP[i], df$FP[i], df$FN_use[i])
}
cat("Done.\n")

# --- Pick tools with enough benchmark coverage --------------------------------
tool_bm <- unique(df[, c("Tool", "Benchmark")])
tool_counts <- table(tool_bm$Tool)
tools_ok <- names(tool_counts)[tool_counts >= MIN_BENCHMARKS_PER_TOOL]
cat("\nTools with >=", MIN_BENCHMARKS_PER_TOOL, "benchmarks:",
    length(tools_ok), "\n")
cat(paste(tools_ok, collapse = ", "), "\n")

# --- Meta-analysis helper -----------------------------------------------------
run_meta <- function(sub, yi_col, vi_col, tool_name, metric_name) {
  keep <- !is.na(sub[[yi_col]]) & !is.na(sub[[vi_col]]) & sub[[vi_col]] > 0
  d <- sub[keep, ]
  if (nrow(d) < 2) return(NULL)
  res <- tryCatch(
    rma(yi = d[[yi_col]], vi = d[[vi_col]], method = "DL"),
    error = function(e) NULL
  )
  if (is.null(res)) return(NULL)
  data.frame(
    tool   = tool_name,
    metric = metric_name,
    k      = res$k,
    pooled = as.numeric(res$b),
    ci_lb  = res$ci.lb,
    ci_ub  = res$ci.ub,
    I2     = res$I2,
    Q      = res$QE,
    Q_p    = res$QEp,
    tau2   = res$tau2,
    stringsAsFactors = FALSE
  )
}

# --- Per-tool meta-analyses ---------------------------------------------------
results <- list()
for (t in tools_ok) {
  sub <- df[df$Tool == t, ]
  results[[length(results) + 1]] <- run_meta(sub, "prec_yi", "prec_vi", t, "precision")
  results[[length(results) + 1]] <- run_meta(sub, "sens_yi", "sens_vi", t, "sensitivity")
  results[[length(results) + 1]] <- run_meta(sub, "f1_yi",   "f1_vi",   t, "f1")
}
pooled_results <- do.call(rbind, results)
write.csv(pooled_results, file.path(TAB_DIR, "pooled_estimates.csv"),
          row.names = FALSE)
cat("\nPooled estimates written to results/tables/pooled_estimates.csv\n\n")
print(pooled_results, digits = 3)

# --- Forest plots -------------------------------------------------------------
plot_forest <- function(sub, yi_col, vi_col, tool_name, metric_name) {
  keep <- !is.na(sub[[yi_col]]) & !is.na(sub[[vi_col]]) & sub[[vi_col]] > 0
  d <- sub[keep, ]
  if (nrow(d) < 2) return(invisible(NULL))
  res <- tryCatch(
    rma(yi = d[[yi_col]], vi = d[[vi_col]], method = "DL"),
    error = function(e) NULL
  )
  if (is.null(res)) return(invisible(NULL))
  fname <- sprintf("%s/forest_%s_%s.pdf", FIG_DIR,
                   gsub("[^A-Za-z0-9]", "_", tool_name), metric_name)
  pdf(fname, width = 8, height = 4 + 0.3 * nrow(d))
  forest(res,
         slab = paste0("B", d$Benchmark, ": ", d$Dataset),
         header = c("Benchmark / Dataset", sprintf("%s [95%% CI]", metric_name)),
         xlab = metric_name,
         mlab = sprintf("Pooled (DL): I2 = %.0f%%, Q p = %.3f",
                        res$I2, res$QEp))
  dev.off()
  invisible(NULL)
}

cat("\nGenerating forest plots...\n")
for (t in tools_ok) {
  sub <- df[df$Tool == t, ]
  plot_forest(sub, "prec_yi", "prec_vi", t, "precision")
  plot_forest(sub, "sens_yi", "sens_vi", t, "sensitivity")
  plot_forest(sub, "f1_yi",   "f1_vi",   t, "f1")
}
cat("  Forest plots saved to results/figures/\n")

# --- Subgroup analysis: real vs simulated (F1) -------------------------------
cat("\n=== SUBGROUP ANALYSIS: real vs simulated (F1) ===\n")
sub_res <- list()
for (t in tools_ok) {
  d <- df[df$Tool == t & !is.na(df$f1_yi) & !is.na(df$f1_vi) & df$f1_vi > 0, ]
  if (nrow(d) < 3) next
  if (length(unique(d$"Type of dataset")) < 2) next
  res <- tryCatch(
    rma(yi = d$f1_yi, vi = d$f1_vi, mods = ~ d$"Type of dataset",
        method = "DL"),
    error = function(e) NULL
  )
  if (is.null(res)) next
  sub_res[[length(sub_res) + 1]] <- data.frame(
    tool = t, k = res$k,
    QM = res$QM, QM_p = res$QMp,
    significant = res$QMp < 0.05,
    stringsAsFactors = FALSE
  )
}
sub_df <- do.call(rbind, sub_res)
if (!is.null(sub_df)) {
  print(sub_df, digits = 3)
  write.csv(sub_df, file.path(TAB_DIR, "subgroup_real_vs_sim.csv"),
            row.names = FALSE)
}

# --- Sensitivity analysis: with vs without Edgren -----------------------------
cat("\n=== SENSITIVITY ANALYSIS: with vs without Edgren datasets (F1) ===\n")
EDGREN_DATASETS <- c("Dataset2.3", "Dataset5.1", "Dataset6.1",
                     "Dataset6.2", "Dataset7.6", "Dataset7.7",
                     "Dataset8.3", "Dataset10.2")
is_edgren <- function(x) x %in% EDGREN_DATASETS

sens_res <- list()
for (t in tools_ok) {
  d_full  <- df[df$Tool == t & !is.na(df$f1_yi) & !is.na(df$f1_vi) & df$f1_vi > 0, ]
  d_noedg <- d_full[!is_edgren(d_full$Dataset), ]
  n_edg_removed <- nrow(d_full) - nrow(d_noedg)
  if (nrow(d_full) < 2 || nrow(d_noedg) < 2) next
  full <- tryCatch(rma(yi = d_full$f1_yi, vi = d_full$f1_vi, method = "DL"),
                   error = function(e) NULL)
  noe  <- tryCatch(rma(yi = d_noedg$f1_yi, vi = d_noedg$f1_vi, method = "DL"),
                   error = function(e) NULL)
  if (is.null(full) || is.null(noe)) next
  sens_res[[length(sens_res) + 1]] <- data.frame(
    tool             = t,
    edgren_rows      = n_edg_removed,
    f1_with          = as.numeric(full$b),
    I2_with          = full$I2,
    k_with           = full$k,
    f1_without       = as.numeric(noe$b),
    I2_without       = noe$I2,
    k_without        = noe$k,
    delta_f1         = as.numeric(noe$b) - as.numeric(full$b),
    delta_I2         = noe$I2 - full$I2,
    stringsAsFactors = FALSE
  )
}
sens_df <- do.call(rbind, sens_res)
if (!is.null(sens_df)) {
  print(sens_df, digits = 3)
  write.csv(sens_df, file.path(TAB_DIR, "sensitivity_edgren.csv"),
            row.names = FALSE)
}

# --- Edgren truth-set drift: Dataset7.6 (27 fusions) vs 7.7 (99 fusions) ------
# Same raw RNA-seq data, different reference truth sets. Any variation in F1
# between the two datasets for the same tool is attributable to the truth set
# alone (not the tool version, parameters, or preprocessing).
cat("\n=== EDGREN TRUTH-SET DRIFT: Dataset 7.6 (P27) vs 7.7 (P99) ===\n")

d76 <- df[df$Dataset == "Dataset7.6",
          c("Tool", "TP", "FP", "FN_use", "prec_yi", "sens_yi", "f1_yi")]
d77 <- df[df$Dataset == "Dataset7.7",
          c("Tool", "TP", "FP", "FN_use", "prec_yi", "sens_yi", "f1_yi")]
names(d76) <- c("Tool", "TP_27", "FP_27", "FN_27",
                "prec_27", "sens_27", "f1_27")
names(d77) <- c("Tool", "TP_99", "FP_99", "FN_99",
                "prec_99", "sens_99", "f1_99")

if (nrow(d76) > 0 && nrow(d77) > 0) {
  edgren_drift <- merge(d76, d77, by = "Tool")
  edgren_drift$delta_f1   <- edgren_drift$f1_99   - edgren_drift$f1_27
  edgren_drift$delta_prec <- edgren_drift$prec_99 - edgren_drift$prec_27
  edgren_drift$delta_sens <- edgren_drift$sens_99 - edgren_drift$sens_27
  print(edgren_drift[, c("Tool", "f1_27", "f1_99", "delta_f1",
                         "prec_27", "prec_99", "sens_27", "sens_99")],
        digits = 3)
  write.csv(edgren_drift, file.path(TAB_DIR, "edgren_truth_set_drift.csv"),
            row.names = FALSE)
  cat(sprintf("\nRange of |delta_F1| across %d tools: %.3f -- %.3f\n",
              nrow(edgren_drift),
              min(abs(edgren_drift$delta_f1), na.rm = TRUE),
              max(abs(edgren_drift$delta_f1), na.rm = TRUE)))
} else {
  cat("No overlap between Dataset7.6 and Dataset7.7 tools found.\n")
}

cat("\nDone. Results in results/figures/ and results/tables/.\n")
