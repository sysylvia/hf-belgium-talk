#!/usr/bin/env Rscript
# ============================================================================
# 01_cate_all_outcomes.R — Comprehensive CATE Analysis
# ============================================================================
# Runs causal forests on all outcomes + factor scores, produces BLP screening
# table, and deep-dives into outcomes with BLP p < 0.20.
#
# Depends on: 00_data_prep.R (imputed covariates), 02_summary_indices.R (factor scores)
#
# Outputs:
#   outputs/tables/blp_screening.csv     — BLP screening table (all outcomes)
#   outputs/figures/BT_blp_screening.png — Visual summary
#   outputs/figures/BT_deepdive_*.png    — Deep-dive figures (if any qualify)
#   data/cate_results.rds                — Full results object
# ============================================================================

library(haven)
library(dplyr)
library(tidyr)
library(grf)
library(ggplot2)
library(patchwork)
library(glue)

# ── Paths ──────────────────────────────────────────────────────────────────
data_dir   <- here::here("data")
fig_path   <- here::here("outputs", "figures")
table_path <- here::here("outputs", "tables")
dir.create(fig_path, showWarnings = FALSE, recursive = TRUE)
dir.create(table_path, showWarnings = FALSE, recursive = TRUE)

# ── Load data with factor scores ───────────────────────────────────────────
df <- readRDS(file.path(data_dir, "hf_with_factors.rds"))
cat(glue("Loaded: {nrow(df)} obs × {ncol(df)} vars\n\n"))

# ── CATE covariates (15 vars, BLVAR_FEED_EBF excluded) ─────────────────────
cate_vars <- c(
  "BBVAR_AGE_MON", "BBVAR_BOY",
  "HHVAR_MOM_AGE_BASE", "HHVAR_MOM_HIEDU", "HHVAR_ASSET",
  "HHVAR_SIB", "HHVAR_SIZE_BASE", "HHVAR_MOM_MIGHIS", "HHVAR_MOM_HOMETOWN",
  "BLVAR_DASSPRI_DEP", "BLVAR_DASSPRI_ANX", "BLVAR_DASSPRI_STR",
  "BLVAR_KNOWPRI_INX", "BLVAR_SSPRI_ALL", "BLVAR_MFEED_DDS"
)

# ── Outcomes to analyze ────────────────────────────────────────────────────
outcomes <- list(
  # Primary outcomes
  list(var = "FLVAR_FEED_DDS",  label = "Dietary Diversity Score",    group = "primary"),
  list(var = "FLVAR_CH_HEMO",   label = "Hemoglobin Level",           group = "primary"),
  # Domain indices (secondary)
  list(var = "IN_CH",           label = "Child Health Index",          group = "domain"),
  list(var = "IN_BFI",          label = "BF Initiation Index",         group = "domain"),
  list(var = "IN_CF",           label = "Complementary Feeding Index", group = "domain"),
  list(var = "IN_PRIKNOW_CH",   label = "Caregiving Knowledge Index",  group = "domain"),
  list(var = "IN_EPDS",         label = "Maternal MH (EPDS) Index",    group = "domain"),
  list(var = "IN_DASS",         label = "Caregiver MH (DASS) Index",   group = "domain"),
  list(var = "IN_SS",           label = "Social Support Index",        group = "domain"),
  # Factor scores
  list(var = "factor_caregiver", label = "Factor: Mental Health",      group = "factor"),
  list(var = "factor_child",     label = "Factor: Support/Knowledge",  group = "factor"),
  list(var = "factor_overall",   label = "Factor: Overall (1-factor)", group = "factor")
)

# FLVAR_FEED_EBF excluded: N=210 too small for causal forest

cat("═══════════════════════════════════════════════════════════════\n")
cat("  CATE Analysis — All Outcomes\n")
cat(glue("  Outcomes: {length(outcomes)}\n"))
cat(glue("  Covariates: {length(cate_vars)} (BLVAR_FEED_EBF excluded)\n"))
cat(glue("  Trees: 4000 · Tuning: all · Seed: 2026\n"))
cat("═══════════════════════════════════════════════════════════════\n\n")


# ============================================================================
# STEP 1: Run causal forests for all outcomes
# ============================================================================

results <- list()

for (i in seq_along(outcomes)) {
  o <- outcomes[[i]]
  cat(glue("\n--- [{i}/{length(outcomes)}] {o$label} ({o$var}) ---\n"))

  # Prepare data: drop NAs on this outcome + covariates
  cf_df <- df |>
    select(all_of(c(o$var, "treat", "TID", cate_vars))) |>
    drop_na()

  n_obs <- nrow(cf_df)
  cat(glue("  N = {n_obs}\n"))

  if (n_obs < 200) {
    cat("  SKIPPED: too few observations\n")
    next
  }

  Y <- as.numeric(cf_df[[o$var]])
  W <- as.numeric(cf_df$treat)
  X <- cf_df |> select(all_of(cate_vars)) |> as.matrix()
  clusters <- as.integer(as.factor(cf_df$TID))

  # Fit causal forest with tuning
  set.seed(2026)
  cf <- causal_forest(
    X, Y, W,
    clusters     = clusters,
    num.trees    = 4000,
    tune.parameters = "all",
    seed         = 2026
  )

  # ATE
  ate <- average_treatment_effect(cf)

  # BLP calibration test
  blp <- test_calibration(cf)
  blp_mat <- as.matrix(blp)
  blp_coef <- blp_mat["differential.forest.prediction", "Estimate"]
  blp_se   <- blp_mat["differential.forest.prediction", "Std. Error"]
  blp_p    <- blp_mat["differential.forest.prediction", ncol(blp_mat)]

  # RATE
  cates <- predict(cf)$predictions
  rate  <- rank_average_treatment_effect(cf, cates)
  rate_p <- 2 * pnorm(-abs(rate$estimate / rate$std.err))

  # Variable importance
  vi <- variable_importance(cf)

  # Store results
  results[[o$var]] <- list(
    var        = o$var,
    label      = o$label,
    group      = o$group,
    n          = n_obs,
    ate_est    = ate[["estimate"]],
    ate_se     = ate[["std.err"]],
    blp_coef   = blp_coef,
    blp_se     = blp_se,
    blp_p      = blp_p,
    rate_est   = rate$estimate,
    rate_se    = rate$std.err,
    rate_p     = rate_p,
    cate_mean  = mean(cates),
    cate_sd    = sd(cates),
    vi         = setNames(as.numeric(vi), cate_vars),
    cates      = cates,
    cf         = cf,
    cf_df      = cf_df
  )

  cat(glue("  ATE = {round(ate[['estimate']], 4)} (SE={round(ate[['std.err']], 4)})\n"))
  cat(glue("  BLP p = {round(blp_p, 4)}\n"))
  cat(glue("  RATE = {round(rate$estimate, 4)} (p={round(rate_p, 4)})\n"))
}


# ============================================================================
# STEP 2: BLP Screening Table
# ============================================================================

cat("\n\n═══════════════════════════════════════════════════════════════\n")
cat("  BLP Screening Table\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

screening <- bind_rows(lapply(results, function(r) {
  tibble(
    Outcome   = r$label,
    Group     = r$group,
    N         = r$n,
    ATE       = round(r$ate_est, 4),
    ATE_SE    = round(r$ate_se, 4),
    BLP_coef  = round(r$blp_coef, 4),
    BLP_p     = round(r$blp_p, 4),
    RATE      = round(r$rate_est, 4),
    RATE_p    = round(r$rate_p, 4),
    CATE_SD   = round(r$cate_sd, 4)
  )
}))

screening <- screening |> arrange(BLP_p)
print(screening, n = 20)
write.csv(screening, file.path(table_path, "blp_screening.csv"), row.names = FALSE)
cat(glue("\nSaved: blp_screening.csv\n"))

# Count qualifying outcomes
n_qualify <- sum(screening$BLP_p < 0.20)
cat(glue("\nOutcomes with BLP p < 0.20: {n_qualify} / {nrow(screening)}\n"))
qualifying <- screening |> filter(BLP_p < 0.20)
if (n_qualify > 0) {
  cat("Qualifying for deep-dive:\n")
  for (i in 1:nrow(qualifying)) {
    cat(glue("  - {qualifying$Outcome[i]} (BLP p = {qualifying$BLP_p[i]})\n"))
  }
}

# BLP screening figure
add_stars <- function(p) case_when(p < 0.01 ~ "***", p < 0.05 ~ "**", p < 0.10 ~ "*", p < 0.20 ~ "†", TRUE ~ "")

screening_plot <- screening |>
  mutate(
    Outcome = factor(Outcome, levels = rev(Outcome)),
    sig = BLP_p < 0.20,
    blp_label = glue("{round(BLP_p, 3)}{add_stars(BLP_p)}")
  )

p_screen <- ggplot(screening_plot, aes(y = Outcome, x = -log10(BLP_p), fill = sig)) +
  geom_col(alpha = 0.85, width = 0.6) +
  geom_vline(xintercept = -log10(0.20), linetype = "dashed", color = "#c00000", linewidth = 0.7) +
  geom_vline(xintercept = -log10(0.10), linetype = "dotted", color = "#767171", linewidth = 0.5) +
  geom_text(aes(label = blp_label), hjust = -0.1, size = 3.5, color = "#333") +
  scale_fill_manual(values = c("FALSE" = "#cccccc", "TRUE" = "#2e75b6"), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.3))) +
  labs(
    title    = "BLP Heterogeneity Screening — All Outcomes",
    subtitle = "Red dashed = p < 0.20 threshold · † p<0.20, * p<0.10, ** p<0.05, *** p<0.01",
    x        = expression(-log[10](p)),
    y        = NULL,
    caption  = glue("Causal forest · 4,000 trees · tuned · 15 covariates · clustered by township")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", color = "#1f4e79", size = 13),
    plot.subtitle = element_text(color = "#767171", size = 10),
    plot.caption  = element_text(color = "#767171", size = 8.5),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(fig_path, "BT_blp_screening.png"), p_screen, width = 10, height = 7, dpi = 200)
cat("Saved: BT_blp_screening.png\n")


# ============================================================================
# STEP 3: Deep-Dive for BLP p < 0.20
# ============================================================================

if (n_qualify > 0) {
  cat("\n\n═══════════════════════════════════════════════════════════════\n")
  cat("  Deep-Dive Heterogeneity Analysis\n")
  cat("═══════════════════════════════════════════════════════════════\n")

  for (q_row in 1:nrow(qualifying)) {
    outcome_var <- names(results)[which(sapply(results, function(r) r$label) == qualifying$Outcome[q_row])]
    r <- results[[outcome_var]]
    cat(glue("\n\n=== Deep-Dive: {r$label} (BLP p = {round(r$blp_p, 4)}) ===\n\n"))

    cates <- r$cates
    ate_est <- r$ate_est

    # 1. CATE distribution
    p1 <- ggplot(tibble(cate = cates), aes(x = cate)) +
      geom_histogram(binwidth = diff(range(cates))/30, fill = "#2e75b6", color = "white", alpha = 0.85) +
      geom_vline(xintercept = ate_est, color = "#ed7d31", linewidth = 1.4) +
      annotate("text", x = ate_est, y = Inf, label = glue("ATE = {round(ate_est, 3)}"),
               vjust = 2, hjust = -0.1, color = "#ed7d31", fontface = "bold", size = 4) +
      labs(title = glue("CATE Distribution — {r$label}"),
           x = "Estimated CATE", y = "Count") +
      theme_minimal(base_size = 11) +
      theme(plot.title = element_text(face = "bold", color = "#1f4e79"))

    # 2. Variable importance
    vi_df <- tibble(variable = cate_vars, importance = r$vi) |>
      arrange(desc(importance)) |> slice_head(n = 8) |>
      mutate(variable = forcats::fct_reorder(variable, importance))

    p2 <- ggplot(vi_df, aes(x = importance, y = variable)) +
      geom_col(fill = "#1f4e79", alpha = 0.85) +
      labs(title = "Variable Importance", x = "Importance", y = NULL) +
      theme_minimal(base_size = 11) +
      theme(plot.title = element_text(face = "bold", color = "#1f4e79"))

    # 3. GATES (quintile ATEs)
    cf_df <- r$cf_df
    cf_df$cate <- cates
    cf_df$quintile <- cut(cates, breaks = quantile(cates, probs = 0:5/5),
                          include.lowest = TRUE, labels = paste0("Q", 1:5))

    gates <- cf_df |>
      group_by(quintile) |>
      summarise(
        ate = mean(.data[[r$var]][treat == 1], na.rm = TRUE) - mean(.data[[r$var]][treat == 0], na.rm = TRUE),
        n = n(),
        .groups = "drop"
      )

    p3 <- ggplot(gates, aes(x = quintile, y = ate, fill = quintile == "Q5")) +
      geom_col(alpha = 0.85, width = 0.6) +
      geom_hline(yintercept = ate_est, linetype = "dashed", color = "#ed7d31", linewidth = 0.8) +
      geom_text(aes(label = round(ate, 3)), vjust = -0.5, size = 3.5, fontface = "bold") +
      scale_fill_manual(values = c("FALSE" = "#cccccc", "TRUE" = "#2e75b6"), guide = "none") +
      labs(title = "GATES — Quintile ATEs",
           subtitle = "Q5 = highest predicted CATE",
           x = "CATE Quintile", y = "ATE") +
      theme_minimal(base_size = 11) +
      theme(plot.title = element_text(face = "bold", color = "#1f4e79"))

    # 4. TOC curve
    rate_obj <- rank_average_treatment_effect(r$cf, cates)
    toc <- as_tibble(rate_obj$TOC)
    names(toc) <- tolower(names(toc))
    toc <- toc |> rename_with(~case_when(
      . %in% c("fraction","q") ~ "fraction",
      . %in% c("estimate","toc") ~ "toc_est",
      . %in% c("std.err","se","stderr") ~ "toc_se",
      TRUE ~ .
    ))

    p4 <- ggplot(toc, aes(x = fraction, y = toc_est)) +
      geom_ribbon(aes(ymin = toc_est - 1.96*toc_se, ymax = toc_est + 1.96*toc_se),
                  fill = "#1f4e79", alpha = 0.15) +
      geom_line(color = "#1f4e79", linewidth = 1.2) +
      geom_hline(yintercept = ate_est, color = "#ed7d31", linewidth = 0.8, linetype = "dashed") +
      scale_x_continuous(labels = scales::percent) +
      labs(title = "TOC Curve",
           subtitle = glue("RATE = {round(r$rate_est, 3)} (p = {round(r$rate_p, 3)})"),
           x = "Fraction treated", y = "ATE for treated fraction") +
      theme_minimal(base_size = 11) +
      theme(plot.title = element_text(face = "bold", color = "#1f4e79"))

    # Combine
    combined <- (p1 | p2) / (p3 | p4) +
      plot_annotation(
        title = glue("Heterogeneity Deep-Dive: {r$label}"),
        subtitle = glue("N = {r$n} · BLP p = {round(r$blp_p, 4)} · ATE = {round(ate_est, 4)}"),
        theme = theme(
          plot.title = element_text(face = "bold", color = "#1f4e79", size = 14),
          plot.subtitle = element_text(color = "#767171", size = 11)
        )
      )

    fname <- glue("BT_deepdive_{gsub('[^a-zA-Z0-9]', '_', r$var)}.png")
    ggsave(file.path(fig_path, fname), combined, width = 14, height = 10, dpi = 200)
    cat(glue("  Saved: {fname}\n"))
  }
}

# ============================================================================
# STEP 4: Save results
# ============================================================================

# Strip large objects for storage (keep summaries, drop full CF objects)
results_slim <- lapply(results, function(r) {
  r$cf <- NULL
  r$cf_df <- NULL
  r
})
saveRDS(results_slim, file.path(data_dir, "cate_results.rds"))
cat(glue("\nSaved: cate_results.rds\n"))

cat("\n═══════════════════════════════════════════════════════════════\n")
cat("  CATE analysis complete.\n")
cat(glue("  Outcomes analyzed: {length(results)}\n"))
cat(glue("  Deep-dive threshold: BLP p < 0.20\n"))
cat(glue("  Qualifying outcomes: {n_qualify}\n"))
cat("═══════════════════════════════════════════════════════════════\n")
