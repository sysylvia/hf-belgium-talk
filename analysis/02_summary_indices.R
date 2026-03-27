#!/usr/bin/env Rscript
# ============================================================================
# 02_summary_indices.R — EFA → CFA → Factor Scores
# ============================================================================
# Purpose: Construct summary indices from 7 domain outcomes using exploratory
#          factor analysis to determine structure, then confirmatory factor
#          analysis with FIML to handle missing data and extract scores.
#
# Outputs:
#   data/hf_with_factors.rds  — Full dataset with factor scores appended
#   outputs/tables/efa_results.csv
#   outputs/figures/BT_efa_scree.png
#   outputs/figures/BT_cfa_loadings.png
#
# Usage: Rscript analysis/02_summary_indices.R
# ============================================================================

library(haven)
library(dplyr)
library(psych)
library(lavaan)
library(ggplot2)
library(glue)

# ── Paths ──────────────────────────────────────────────────────────────────
data_dir   <- here::here("data")
fig_path   <- here::here("outputs", "figures")
table_path <- here::here("outputs", "tables")
dir.create(fig_path,   showWarnings = FALSE, recursive = TRUE)
dir.create(table_path, showWarnings = FALSE, recursive = TRUE)

# ── Load ───────────────────────────────────────────────────────────────────
df <- readRDS(file.path(data_dir, "hf_analysis_full.rds"))
cat(glue("Loaded: {nrow(df)} obs\n\n"))

domain_vars <- c("IN_CH", "IN_BFI", "IN_CF", "IN_PRIKNOW_CH", "IN_EPDS", "IN_DASS", "IN_SS")
domain_labels <- c("Child Health", "BF Initiation", "Complementary Feeding",
                    "Caregiving Knowledge", "Maternal MH (EPDS)",
                    "Caregiver MH (DASS)", "Social Support")

domain_df <- df[, domain_vars]
cat(glue("Complete cases across all 7 domains: {sum(complete.cases(domain_df))}\n"))
cat(glue("Pairwise available: {nrow(domain_df) - max(colSums(is.na(domain_df)))}\n\n"))


# ============================================================================
# 1. EXPLORATORY FACTOR ANALYSIS
# ============================================================================

cat("═══════════════════════════════════════════════════════════════\n")
cat("  Step 1: Exploratory Factor Analysis\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

# Use pairwise complete observations for EFA
domain_complete <- domain_df[complete.cases(domain_df), ]
cat(glue("EFA sample: {nrow(domain_complete)} complete cases\n\n"))

# 1a. Parallel analysis (Horn's method) — determine number of factors
cat("=== Parallel Analysis ===\n")
pa <- fa.parallel(domain_complete, fm = "pa", fa = "fa", n.iter = 100,
                  plot = FALSE, quant = 0.95)
cat(glue("Suggested factors (parallel analysis): {pa$nfact}\n"))
cat(glue("Eigenvalues: {paste(round(pa$fa.values, 3), collapse = ', ')}\n\n"))

# 1b. Multiple criteria (limit to 3 factors max to avoid singularity)
cat("=== Additional Criteria ===\n")
tryCatch({
  vss_result <- VSS(domain_complete, n = 3, fm = "pa", plot = FALSE, rotate = "oblimin")
  cat(glue("VSS complexity 1 optimal: {which.max(vss_result$cfit.1)}\n"))
  cat(glue("VSS complexity 2 optimal: {which.max(vss_result$cfit.2)}\n\n"))
}, error = function(e) {
  cat(glue("VSS encountered numerical issues (common with 7 vars): {e$message}\n"))
  cat("Proceeding with parallel analysis result.\n\n")
})

# 1c. Fit 2-factor and 3-factor EFA with oblique rotation
n_factors_to_try <- min(pa$nfact, 3)  # Use parallel analysis suggestion, cap at 3

for (nf in 2:n_factors_to_try) {
  cat(glue("=== {nf}-Factor EFA (oblimin rotation) ===\n"))
  efa_fit <- fa(domain_complete, nfactors = nf, fm = "pa", rotate = "oblimin")
  print(efa_fit$loadings, cutoff = 0.2)
  cat(glue("\nRMSEA: {round(efa_fit$RMSEA[1], 4)}, BIC: {round(efa_fit$BIC, 1)}\n"))
  var_explained <- sum(efa_fit$Vaccounted["Proportion Var", 1:nf])
  cat(glue("Cumulative variance: {round(var_explained * 100, 1)}%\n\n"))
  if (nf == 2) efa2 <- efa_fit
  if (nf == 3) efa3 <- efa_fit
}

# 1d. Save EFA loadings
efa_loadings <- as.data.frame(unclass(efa2$loadings))
efa_loadings$variable <- rownames(efa_loadings)
efa_loadings$label <- domain_labels
write.csv(efa_loadings, file.path(table_path, "efa_results.csv"), row.names = FALSE)

# 1e. Scree plot
png(file.path(fig_path, "BT_efa_scree.png"), width = 800, height = 500, res = 150)
plot(pa$fa.values, type = "b", pch = 19, col = "#1f4e79",
     xlab = "Factor Number", ylab = "Eigenvalue",
     main = "Parallel Analysis Scree Plot",
     ylim = c(0, max(pa$fa.values) * 1.2))
lines(pa$fa.sim, type = "b", pch = 17, col = "#ed7d31", lty = 2)
abline(h = 0, col = "gray", lty = 3)
legend("topright", c("Actual eigenvalues", "Simulated (95th percentile)"),
       col = c("#1f4e79", "#ed7d31"), pch = c(19, 17), lty = c(1, 2), cex = 0.8)
dev.off()
cat("Scree plot saved: BT_efa_scree.png\n\n")


# ============================================================================
# 2. CONFIRMATORY FACTOR ANALYSIS (PRIMARY: FIML)
# ============================================================================

cat("═══════════════════════════════════════════════════════════════\n")
cat("  Step 2: Confirmatory Factor Analysis (FIML)\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

# Try 2-factor model first, then 3-factor if parallel analysis suggests it
# Factor assignment based on EFA loadings
cfa_2f_model <- '
  caregiver =~ IN_EPDS + IN_DASS + IN_SS + IN_PRIKNOW_CH
  child     =~ IN_CH + IN_CF + IN_BFI
'

# Fit with FIML (handles missing data)
cfa_fit_2f <- cfa(cfa_2f_model, data = domain_df, missing = "fiml", estimator = "MLR")

# Also try 3-factor if suggested
if (n_factors_to_try >= 3) {
  cfa_3f_model <- '
    mental    =~ IN_EPDS + IN_DASS
    support   =~ IN_SS + IN_PRIKNOW_CH
    child     =~ IN_CH + IN_CF + IN_BFI
  '
  cfa_fit_3f <- tryCatch(
    cfa(cfa_3f_model, data = domain_df, missing = "fiml", estimator = "MLR"),
    error = function(e) { cat(glue("3-factor CFA failed: {e$message}\n")); NULL }
  )
}

# Compare and select best
cat("=== Model Comparison ===\n")
fit_2f_stats <- fitMeasures(cfa_fit_2f, c("cfi", "rmsea", "srmr", "aic", "bic"))
cat(sprintf("  2-Factor: CFI=%.3f, RMSEA=%.3f, AIC=%.1f\n", fit_2f_stats["cfi"], fit_2f_stats["rmsea"], fit_2f_stats["aic"]))

if (exists("cfa_fit_3f") && !is.null(cfa_fit_3f)) {
  fit_3f_stats <- fitMeasures(cfa_fit_3f, c("cfi", "rmsea", "srmr", "aic", "bic"))
  cat(sprintf("  3-Factor: CFI=%.3f, RMSEA=%.3f, AIC=%.1f\n", fit_3f_stats["cfi"], fit_3f_stats["rmsea"], fit_3f_stats["aic"]))
  # Select lower AIC
  if (fit_3f_stats["aic"] < fit_2f_stats["aic"]) {
    cat("  → 3-factor model preferred (lower AIC)\n\n")
    cfa_fit <- cfa_fit_3f
    n_factors_selected <- 3
  } else {
    cat("  → 2-factor model preferred (lower AIC)\n\n")
    cfa_fit <- cfa_fit_2f
    n_factors_selected <- 2
  }
} else {
  cfa_fit <- cfa_fit_2f
  n_factors_selected <- 2
  cat("  → Using 2-factor model\n\n")
}

cat("=== CFA Fit Indices ===\n")
fit_measures <- fitMeasures(cfa_fit, c("chisq", "df", "pvalue",
                                        "cfi", "tli", "rmsea", "rmsea.ci.lower",
                                        "rmsea.ci.upper", "srmr"))
for (m in names(fit_measures)) {
  cat(sprintf("  %-20s: %.4f\n", m, fit_measures[m]))
}

cat(glue("\n  CFI: {round(fit_measures['cfi'], 3)} (target: >0.90)\n"))
cat(glue("  RMSEA: {round(fit_measures['rmsea'], 3)} (target: <0.08)\n"))
cat(glue("  SRMR: {round(fit_measures['srmr'], 3)} (target: <0.08)\n\n"))

# Standardized loadings
cat("=== Standardized Loadings ===\n")
std_loadings <- standardizedSolution(cfa_fit)
std_loadings_display <- std_loadings[std_loadings$op == "=~", c("lhs", "rhs", "est.std", "pvalue")]
print(std_loadings_display, row.names = FALSE)

# Extract factor scores
cat("\n=== Factor Scores ===\n")
factor_scores <- lavPredict(cfa_fit, type = "lv")
cat(glue("Factor scores extracted for {nrow(factor_scores)} observations\n"))
cat(glue("Caregiver factor: mean={round(mean(factor_scores[,1]),3)}, sd={round(sd(factor_scores[,1]),3)}\n"))
cat(glue("Child factor: mean={round(mean(factor_scores[,2]),3)}, sd={round(sd(factor_scores[,2]),3)}\n"))
cat(glue("Factor correlation: {round(cor(factor_scores[,1], factor_scores[,2]),3)}\n\n"))


# ============================================================================
# 3. ALSO FIT 1-FACTOR MODEL (for comparison)
# ============================================================================

cat("=== 1-Factor Model (for comparison) ===\n")
cfa_1f_model <- 'overall =~ IN_EPDS + IN_DASS + IN_SS + IN_PRIKNOW_CH + IN_CH + IN_CF + IN_BFI'
cfa_1f_fit <- cfa(cfa_1f_model, data = domain_df, missing = "fiml", estimator = "MLR")
fit_1f <- fitMeasures(cfa_1f_fit, c("cfi", "rmsea", "srmr", "aic", "bic"))
fit_2f <- fitMeasures(cfa_fit, c("cfi", "rmsea", "srmr", "aic", "bic"))

cat(sprintf("  1-Factor: CFI=%.3f, RMSEA=%.3f, AIC=%.1f\n", fit_1f["cfi"], fit_1f["rmsea"], fit_1f["aic"]))
cat(sprintf("  2-Factor: CFI=%.3f, RMSEA=%.3f, AIC=%.1f\n", fit_2f["cfi"], fit_2f["rmsea"], fit_2f["aic"]))
cat(sprintf("  ΔAIC: %.1f (negative = 2-factor better)\n\n", fit_2f["aic"] - fit_1f["aic"]))

# Also extract 1-factor score
overall_scores <- lavPredict(cfa_1f_fit, type = "lv")


# ============================================================================
# 4. SAVE AUGMENTED DATASET
# ============================================================================

cat("═══════════════════════════════════════════════════════════════\n")
cat("  Saving augmented dataset\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

df$factor_caregiver <- factor_scores[, 1]
df$factor_child     <- factor_scores[, 2]
df$factor_overall   <- overall_scores[, 1]

saveRDS(df, file.path(data_dir, "hf_with_factors.rds"))
cat(glue("Saved: hf_with_factors.rds ({nrow(df)} obs × {ncol(df)} vars)\n"))
cat("  New columns: factor_caregiver, factor_child, factor_overall\n\n")

# CFA loadings figure
loadings_df <- std_loadings[std_loadings$op == "=~", ]
loadings_df$label <- domain_labels[match(loadings_df$rhs, domain_vars)]
loadings_df$label <- factor(loadings_df$label, levels = rev(domain_labels))

p <- ggplot(loadings_df, aes(x = est.std, y = label, fill = lhs)) +
  geom_col(position = "dodge", alpha = 0.85, width = 0.6) +
  geom_vline(xintercept = 0, linewidth = 0.5) +
  scale_fill_manual(values = c("caregiver" = "#ed7d31", "child" = "#1f4e79"),
                    labels = c("Caregiver Factor", "Child Factor")) +
  labs(title = "CFA Standardized Loadings",
       subtitle = glue("2-factor model · CFI={round(fit_measures['cfi'],3)} · RMSEA={round(fit_measures['rmsea'],3)}"),
       x = "Standardized Loading", y = NULL, fill = NULL) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", color = "#1f4e79"),
        plot.subtitle = element_text(color = "#767171"),
        legend.position = "bottom")

ggsave(file.path(fig_path, "BT_cfa_loadings.png"), p, width = 9, height = 5, dpi = 200)
cat("Saved: BT_cfa_loadings.png\n")

cat("\n═══════════════════════════════════════════════════════════════\n")
cat("  Summary index construction complete.\n")
cat("═══════════════════════════════════════════════════════════════\n")
