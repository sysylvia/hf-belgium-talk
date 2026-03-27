#!/usr/bin/env Rscript
# ============================================================================
# 03_domain_indices.R — Anderson ICW Domain Summary Indices
# ============================================================================
# Constructs 3 inverse-covariance-weighted (ICW) summary indices from the
# 3-factor CFA structure identified in 02_summary_indices.R.
#
# Domains:
#   1. Mental Health    = IN_EPDS + IN_DASS
#   2. Caregiver Support = IN_SS + IN_PRIKNOW_CH
#   3. Child Nutrition   = IN_CH + IN_CF + IN_BFI
#
# Method: Anderson (2008, JASA) ICW index
#   - Standardize to control-group mean/SD
#   - Weight by inverse covariance matrix (control group)
#   - Index = normalized weighted sum
#
# Output: data/hf_with_indices.rds
# ============================================================================

library(haven)
library(dplyr)
library(glue)

data_dir <- here::here("data")
df <- readRDS(file.path(data_dir, "hf_analysis_full.rds"))
cat(glue("Loaded: {nrow(df)} obs\n\n"))

# ── Anderson ICW function ──────────────────────────────────────────────────

make_anderson_icw <- function(df, components, treat_var = "treat", index_name = "index") {
  ctrl <- df[[treat_var]] == 0

  # Step 1: Standardize to control-group mean and SD
  z_mat <- sapply(components, function(v) {
    x <- df[[v]]
    mu <- mean(x[ctrl], na.rm = TRUE)
    sg <- sd(x[ctrl], na.rm = TRUE)
    if (sg == 0 || is.na(sg)) return(rep(NA_real_, length(x)))
    (x - mu) / sg
  })

  # Step 2: Inverse covariance weights (pairwise complete, control group only)
  z_ctrl <- z_mat[ctrl, , drop = FALSE]
  Sigma <- cov(z_ctrl, use = "pairwise.complete.obs")

  # Handle potential singularity
  Sigma_inv <- tryCatch(
    solve(Sigma),
    error = function(e) {
      cat(glue("  Warning: Singular covariance for {index_name}, using pseudoinverse\n"))
      MASS::ginv(Sigma)
    }
  )

  # Step 3: Weights = normalized row sums of inverse covariance
  weights <- rowSums(Sigma_inv)
  weights_norm <- weights / sum(weights)

  # Step 4: Weighted sum
  # Handle NAs: compute for rows with all components available
  complete_rows <- complete.cases(z_mat)
  index_vals <- rep(NA_real_, nrow(df))
  index_vals[complete_rows] <- as.vector(z_mat[complete_rows, , drop = FALSE] %*% weights_norm)

  list(
    index = index_vals,
    weights = setNames(weights_norm, components),
    n_complete = sum(complete_rows),
    n_missing = sum(!complete_rows)
  )
}

# ── Construct 3 domain indices ─────────────────────────────────────────────

cat("═══════════════════════════════════════════════════════════════\n")
cat("  Anderson ICW Domain Index Construction\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

domains <- list(
  list(name = "idx_mental",  label = "Mental Health",      components = c("IN_EPDS", "IN_DASS")),
  list(name = "idx_support", label = "Caregiver Support",  components = c("IN_SS", "IN_PRIKNOW_CH")),
  list(name = "idx_child",   label = "Child Nutrition",    components = c("IN_CH", "IN_CF", "IN_BFI"))
)

for (d in domains) {
  cat(glue("--- {d$label} ({d$name}) ---\n"))
  cat(glue("  Components: {paste(d$components, collapse = ', ')}\n"))

  result <- make_anderson_icw(df, d$components, index_name = d$name)
  df[[d$name]] <- result$index

  cat(glue("  ICW weights: {paste(sprintf('%s=%.3f', names(result$weights), result$weights), collapse = ', ')}\n"))
  cat(glue("  Complete: {result$n_complete}, Missing: {result$n_missing}\n"))

  # Descriptives by treatment
  ctrl_mean <- mean(df[[d$name]][df$treat == 0], na.rm = TRUE)
  treat_mean <- mean(df[[d$name]][df$treat == 1], na.rm = TRUE)
  cat(glue("  Control mean: {round(ctrl_mean, 4)}, Treatment mean: {round(treat_mean, 4)}\n"))
  cat(glue("  Diff: {round(treat_mean - ctrl_mean, 4)}\n\n"))
}

# ── Save ───────────────────────────────────────────────────────────────────

out_path <- file.path(data_dir, "hf_with_indices.rds")
saveRDS(df, out_path)
cat(glue("Saved: {out_path} ({nrow(df)} obs × {ncol(df)} vars)\n"))
cat(glue("New columns: idx_mental, idx_support, idx_child\n"))

cat("\n═══════════════════════════════════════════════════════════════\n")
cat("  Domain index construction complete.\n")
cat("═══════════════════════════════════════════════════════════════\n")
