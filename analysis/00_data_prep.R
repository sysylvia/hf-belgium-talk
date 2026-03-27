#!/usr/bin/env Rscript
# ============================================================================
# 00_data_prep.R — HF Belgium Talk Data Preparation
# ============================================================================
# Purpose: Read HF_main_anasam.dta, diagnose missingness, impute baseline
#          covariates, and export analysis-ready datasets.
#
# Key decision: Drop BLVAR_FEED_EBF from CATE covariates (39% structurally
#               missing — pregnancy enrollees had no breastfeeding baseline).
#               This recovers ~240 observations for the causal forest.
#
# Outputs:
#   data/hf_analysis_full.rds     — Full dataset with imputed covariates
#   data/hf_cate_ready.rds        — CF-ready subset (complete CATE covariates)
#   outputs/missingness_report.txt — Diagnostic summary
#
# Usage: Rscript analysis/00_data_prep.R
#        (or source from hf_analysis.qmd)
# ============================================================================

library(haven)
library(dplyr)
library(tidyr)
library(glue)

# Optional but recommended for imputation diagnostics
if (!requireNamespace("naniar",     quietly = TRUE)) install.packages("naniar")
if (!requireNamespace("missForest", quietly = TRUE)) install.packages("missForest")

library(naniar)
library(missForest)

# ── Paths ──────────────────────────────────────────────────────────────────
data_dir   <- here::here("data")
output_dir <- here::here("outputs")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

dta_path <- file.path(data_dir, "HF_main_anasam.dta")
stopifnot("Data file not found" = file.exists(dta_path))


# ============================================================================
# 1. LOAD AND INSPECT
# ============================================================================

cat("═══════════════════════════════════════════════════════════════\n")
cat("  HF Belgium Talk — Data Preparation\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

df <- read_dta(dta_path)
cat(glue("Loaded: {nrow(df)} observations × {ncol(df)} variables\n"))
cat(glue("Treatment: {sum(df$treat == 1)} treated, {sum(df$treat == 0)} control\n\n"))


# ============================================================================
# 2. DEFINE VARIABLE SETS
# ============================================================================

# Primary outcomes
primary_outcomes <- c(
  "FLVAR_CH_HEMO",   # Hemoglobin (g/dL)
  "FLVAR_FEED_EBF",  # Exclusive breastfeeding (0/1)
  "FLVAR_FEED_DDS"   # Dietary diversity score
)

# Domain indices (secondary outcomes)
domain_indices <- c(
  "IN_CH",           # Child Health Index
  "IN_BFI",          # Breastfeeding Initiation
  "IN_CF",           # Complementary Feeding
  "IN_PRIKNOW_CH",   # Caregiving Knowledge
  "IN_EPDS",         # Maternal Mental Health (EPDS)
  "IN_DASS",         # Caregiver Mental Health (DASS)
  "IN_SS"            # Perceived Social Support
)

# CATE covariates — REVISED: BLVAR_FEED_EBF DROPPED (structural missingness)
cate_covariates <- c(
  "BBVAR_AGE_MON",        # Child age at baseline (months)
  "BBVAR_BOY",            # Child sex (1 = boy)
  "HHVAR_MOM_AGE_BASE",   # Mother's age at baseline
  "HHVAR_MOM_HIEDU",      # Mother's education (high school+)
  "HHVAR_ASSET",          # Household asset index
  "HHVAR_SIB",            # Number of siblings
  "HHVAR_SIZE_BASE",      # Household size
  "HHVAR_MOM_MIGHIS",     # Mother migration history
  "HHVAR_MOM_HOMETOWN",   # Mother is from hometown
  "BLVAR_DASSPRI_DEP",    # Baseline depression (DASS)
  "BLVAR_DASSPRI_ANX",    # Baseline anxiety (DASS)
  "BLVAR_DASSPRI_STR",    # Baseline stress (DASS)
  "BLVAR_KNOWPRI_INX",    # Baseline knowledge index
  "BLVAR_SSPRI_ALL",      # Baseline social support
  "BLVAR_MFEED_DDS"       # Baseline dietary diversity
  # BLVAR_FEED_EBF EXCLUDED: 39% missing (structural — pregnancy enrollees)
)

# LASSO-selected controls (pre-generated in dataset)
lasso_cols <- df |>
  select(starts_with(c("_CHILD_", "_FAM_", "_CARE_", "_CID_"))) |>
  names()

# Core identifiers
id_vars <- c("treat", "TID", "CID")


# ============================================================================
# 3. MISSINGNESS DIAGNOSIS
# ============================================================================

cat("── Missingness Diagnosis ─────────────────────────────────────\n\n")

# 3a. Outcome missingness
cat("=== Outcome Variables ===\n")
all_outcomes <- c(primary_outcomes, domain_indices)
outcome_miss <- sapply(all_outcomes, function(v) {
  n_miss <- sum(is.na(df[[v]]))
  pct    <- round(100 * n_miss / nrow(df), 1)
  cat(sprintf("  %-20s: %4d missing (%5.1f%%)\n", v, n_miss, pct))
  c(n_miss = n_miss, pct = pct)
})

# 3b. CATE covariate missingness
cat("\n=== CATE Covariates (revised set, 15 vars) ===\n")
cate_miss <- sapply(cate_covariates, function(v) {
  n_miss <- sum(is.na(df[[v]]))
  pct    <- round(100 * n_miss / nrow(df), 1)
  if (n_miss > 0) {
    cat(sprintf("  %-25s: %4d missing (%5.1f%%)\n", v, n_miss, pct))
  }
  c(n_miss = n_miss, pct = pct)
})
cat(sprintf("  Variables with zero missing: %d / %d\n",
            sum(cate_miss["n_miss", ] == 0), length(cate_covariates)))

# 3c. Dropped variable report
cat("\n=== Structurally Excluded ===\n")
blvar_ebf_miss <- sum(is.na(df[["BLVAR_FEED_EBF"]]))
cat(sprintf("  BLVAR_FEED_EBF: %d missing (%.1f%%) — DROPPED from CATE covariates\n",
            blvar_ebf_miss, 100 * blvar_ebf_miss / nrow(df)))
cat("  Reason: Structurally absent for pregnancy enrollees (no BF at enrollment)\n")

# 3d. Complete-case comparison
cc_old <- complete.cases(df[, c("FLVAR_FEED_DDS", "treat", "TID", cate_covariates, "BLVAR_FEED_EBF")])
cc_new <- complete.cases(df[, c("FLVAR_FEED_DDS", "treat", "TID", cate_covariates)])

cat(glue("\n=== Complete-Case Comparison (DDS outcome) ===\n"))
cat(sprintf("  Old spec (16 vars incl. BLVAR_FEED_EBF): %d / %d (%.1f%%)\n",
            sum(cc_old), nrow(df), 100 * sum(cc_old) / nrow(df)))
cat(sprintf("  New spec (15 vars excl. BLVAR_FEED_EBF): %d / %d (%.1f%%)\n",
            sum(cc_new), nrow(df), 100 * sum(cc_new) / nrow(df)))
cat(sprintf("  Observations recovered: %d\n\n", sum(cc_new) - sum(cc_old)))


# ============================================================================
# 4. MISSINGNESS MECHANISM ASSESSMENT
# ============================================================================

cat("── Missingness Mechanism Assessment ─────────────────────────\n\n")

# Which CATE covariates have any missing?
vars_with_missing <- cate_covariates[sapply(cate_covariates, function(v) sum(is.na(df[[v]])) > 0)]

if (length(vars_with_missing) > 0) {
  cat(sprintf("Variables with missing data: %s\n", paste(vars_with_missing, collapse = ", ")))
  cat(sprintf("Total obs with any CATE covariate missing: %d\n\n",
              sum(!cc_new) - sum(is.na(df[["FLVAR_FEED_DDS"]]))))  # subtract outcome missingness

  # Test whether missingness correlates with treatment
  cat("=== Missingness × Treatment Balance ===\n")
  for (v in vars_with_missing) {
    miss_indicator <- is.na(df[[v]])
    treat_miss <- mean(miss_indicator[df$treat == 1])
    ctrl_miss  <- mean(miss_indicator[df$treat == 0])
    diff <- treat_miss - ctrl_miss
    cat(sprintf("  %-25s: Treat=%.3f, Control=%.3f, Diff=%.4f\n",
                v, treat_miss, ctrl_miss, diff))
  }
  cat("  (Small diffs suggest MAR/MCAR, not differential attrition)\n\n")
} else {
  cat("No CATE covariates have missing data — no imputation needed.\n\n")
}


# ============================================================================
# 5. IMPUTATION OF BASELINE COVARIATES
# ============================================================================

cat("── Covariate Imputation ──────────────────────────────────────\n\n")

if (length(vars_with_missing) > 0 && sum(sapply(vars_with_missing, function(v) sum(is.na(df[[v]])))) > 0) {

  # Extract covariate matrix for imputation
  impute_vars <- cate_covariates  # impute within the full covariate set
  impute_df   <- df[, impute_vars]

  # Convert all to numeric for missForest
  impute_df <- impute_df |> mutate(across(everything(), as.numeric))

  n_missing_total <- sum(is.na(impute_df))
  cat(sprintf("Total missing values to impute: %d (across %d variables)\n",
              n_missing_total, length(vars_with_missing)))

  if (n_missing_total <= 50) {
    # For very few missing values, use median/mode imputation (simpler, justified)
    cat("Strategy: Median imputation (few missing values, low risk of bias)\n\n")

    imputed_df <- impute_df
    for (v in vars_with_missing) {
      med_val <- median(impute_df[[v]], na.rm = TRUE)
      n_imp   <- sum(is.na(imputed_df[[v]]))
      imputed_df[[v]][is.na(imputed_df[[v]])] <- med_val
      cat(sprintf("  %-25s: %d values imputed with median = %.3f\n", v, n_imp, med_val))
    }

  } else {
    # For more missing values, use missForest
    cat("Strategy: missForest (nonparametric, preserves variable relationships)\n\n")

    set.seed(2026)
    mf_result <- missForest(as.data.frame(impute_df), maxiter = 10, ntree = 100,
                             verbose = FALSE)

    cat(sprintf("  OOB imputation error (NRMSE): %.4f\n", mf_result$OOBerror[1]))
    if (length(mf_result$OOBerror) > 1) {
      cat(sprintf("  OOB imputation error (PFC):   %.4f\n", mf_result$OOBerror[2]))
    }

    imputed_df <- as_tibble(mf_result$ximp)

    for (v in vars_with_missing) {
      n_imp <- sum(is.na(impute_df[[v]]))
      cat(sprintf("  %-25s: %d values imputed\n", v, n_imp))
    }
  }

  # Replace covariates in main dataframe
  df_imputed <- df
  for (v in impute_vars) {
    df_imputed[[v]] <- imputed_df[[v]]
  }

  # Add missing indicators for transparency
  for (v in vars_with_missing) {
    indicator_name <- paste0("mi_", v)
    df_imputed[[indicator_name]] <- as.integer(is.na(df[[v]]))
  }

  cat("\nMissing indicators created: ", paste0("mi_", vars_with_missing, collapse = ", "), "\n")

} else {
  cat("No imputation needed — all CATE covariates complete.\n")
  df_imputed <- df
}


# ============================================================================
# 6. CONSTRUCT ANALYSIS-READY DATASETS
# ============================================================================

cat("\n── Export Analysis-Ready Datasets ────────────────────────────\n\n")

# 6a. Full dataset with imputed covariates
full_path <- file.path(data_dir, "hf_analysis_full.rds")
saveRDS(df_imputed, full_path)
cat(sprintf("  hf_analysis_full.rds: %d obs × %d vars\n", nrow(df_imputed), ncol(df_imputed)))

# 6b. CATE-ready subset (complete on outcome + treatment + covariates)
cate_complete <- complete.cases(df_imputed[, c("FLVAR_FEED_DDS", "treat", "TID", cate_covariates)])
df_cate <- df_imputed[cate_complete, ]
cate_path <- file.path(data_dir, "hf_cate_ready.rds")
saveRDS(df_cate, cate_path)
cat(sprintf("  hf_cate_ready.rds:   %d obs (%.1f%% of full)\n",
            nrow(df_cate), 100 * nrow(df_cate) / nrow(df_imputed)))

# 6c. Compare old vs new CF sample sizes
cat(sprintf("\n  CF sample improvement: %d → %d (+%d observations, +%.0f%%)\n",
            sum(cc_old), nrow(df_cate),
            nrow(df_cate) - sum(cc_old),
            100 * (nrow(df_cate) - sum(cc_old)) / sum(cc_old)))


# ============================================================================
# 7. VALIDATION
# ============================================================================

cat("\n── Validation ────────────────────────────────────────────────\n\n")

# 7a. Treatment balance in CATE-ready sample
cat("=== Treatment Balance (CATE-ready sample) ===\n")
cat(sprintf("  Treated:  %d (%.1f%%)\n",
            sum(df_cate$treat == 1), 100 * mean(df_cate$treat == 1)))
cat(sprintf("  Control:  %d (%.1f%%)\n",
            sum(df_cate$treat == 0), 100 * mean(df_cate$treat == 0)))

# 7b. Key covariate means by treatment (quick balance check)
cat("\n=== Quick Balance Check (CATE-ready sample) ===\n")
balance_vars <- c("BBVAR_AGE_MON", "BBVAR_BOY", "HHVAR_MOM_AGE_BASE",
                  "HHVAR_ASSET", "BLVAR_DASSPRI_DEP", "BLVAR_KNOWPRI_INX")
for (v in balance_vars) {
  if (v %in% names(df_cate)) {
    t_mean <- mean(df_cate[[v]][df_cate$treat == 1], na.rm = TRUE)
    c_mean <- mean(df_cate[[v]][df_cate$treat == 0], na.rm = TRUE)
    t_sd   <- sd(df_cate[[v]], na.rm = TRUE)
    std_diff <- if (t_sd > 0) (t_mean - c_mean) / t_sd else NA
    cat(sprintf("  %-25s: T=%.3f  C=%.3f  StdDiff=%.3f\n", v, t_mean, c_mean, std_diff))
  }
}
cat("  (StdDiff < 0.1 = well-balanced)\n")


# ============================================================================
# 8. SAVE DIAGNOSTIC REPORT
# ============================================================================

report_path <- file.path(output_dir, "missingness_report.txt")

sink(report_path)
cat("HF Belgium Talk — Missingness Report\n")
cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

cat("Dataset: HF_main_anasam.dta\n")
cat(sprintf("Observations: %d\n", nrow(df)))
cat(sprintf("Variables: %d\n\n", ncol(df)))

cat("CATE Covariate Set (revised):\n")
cat(paste(" ", cate_covariates, collapse = "\n"), "\n\n")

cat("EXCLUDED: BLVAR_FEED_EBF (39.2% structural missingness)\n")
cat("Reason: Pregnancy enrollees had no breastfeeding baseline.\n\n")

cat("Imputation Summary:\n")
for (v in vars_with_missing) {
  cat(sprintf("  %s: %d values imputed\n", v, sum(is.na(df[[v]]))))
}

cat(sprintf("\nCF Sample: %d → %d (improvement: +%d)\n",
            sum(cc_old), nrow(df_cate), nrow(df_cate) - sum(cc_old)))

sink()
cat(sprintf("\n  Diagnostic report saved: %s\n", report_path))


# ============================================================================
# SUMMARY
# ============================================================================

cat("\n═══════════════════════════════════════════════════════════════\n")
cat("  Data preparation complete.\n")
cat(sprintf("  CF sample: %d → %d observations (+%.0f%%)\n",
            sum(cc_old), nrow(df_cate),
            100 * (nrow(df_cate) - sum(cc_old)) / sum(cc_old)))
cat("  Key change: BLVAR_FEED_EBF dropped from CATE covariates\n")
cat("═══════════════════════════════════════════════════════════════\n")
