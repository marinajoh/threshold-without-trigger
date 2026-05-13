# =============================================================================
# regression.R
#
# Bayesian negative binomial regressions for the main analysis. Estimates
# the baseline model, short-exposure-window models (0 to 3, 6, and 12
# months), a dynamic-exposure-window model (0-6, 6-12, 12-24, 24+ months),
# comparative regression models (NATO target versus non-NATO target), and
# appendix models for extended short windows and pre-recognition windows.
# All tables are printed to the console in LaTeX form for inclusion in
# the thesis and appendix.
# =============================================================================

library(dplyr)
library(tibble)
library(knitr)
library(kableExtra)
library(lubridate)
library(tidyr)
library(janitor)
library(splines)
library(rstanarm)

options(mc.cores = 1)

# ---- Read derived series ----------------------------------------------------

# Reads aggregated count series and event-level data produced by
# clean_data.R. The list also includes the NATO membership and non-NATO
# initiator vectors so they do not need to be redefined here.

series_path <- "data/derived_series.rds"

obj <- readRDS(series_path)

df_quarter <- obj$df_quarter
cutoff_date <- as.Date(obj$article5_cyber_date)
far_future <- as.Date("2100-01-01")

# ---- Helper: share of period overlap ----------------------------------------

# Computes the share of a period that falls within a target window. Used to
# build the post-2014 exposure indicator: for each period, the share of
# days falling on or after the cutoff date. Periods entirely after the
# cutoff return 1, periods entirely before return 0, and the period
# containing the cutoff returns a fractional value. This is duplicated in
# modeltest.R so that each script runs independently.

share_overlap <- function(start, end, win_start, win_end) {
  den <- as.numeric(end - start)
  num <- pmax(0, as.numeric(pmin(end, win_end) - pmax(start, win_start)))
  num / den
}

# ---- Main regression data ---------------------------------------------------

# Adds period start and end columns to the quarterly series and constructs
# the post-2014 exposure indicator along with several alternative exposure
# windows used for short-window and dynamic-window specifications. The
# natural cubic spline basis is then added as two columns (ns1, ns2) so
# the spline can be referred to in formulas without re-fitting.

dfq <- df_quarter %>%
  mutate(period_start = quarter_start,
         period_end = quarter_start %m+% months(3),
         intercept_1 = 1,
         post_article5_exact = share_overlap(period_start, period_end, cutoff_date, far_future),
         post_0_3m = share_overlap(period_start, period_end, cutoff_date, cutoff_date %m+% months(3)),
         post_0_6m = share_overlap(period_start, period_end, cutoff_date, cutoff_date %m+% months(6)),
         post_0_12m = share_overlap(period_start, period_end, cutoff_date, cutoff_date %m+% months(12)),
         post_6_12m = share_overlap(period_start, period_end, cutoff_date %m+% months(6), cutoff_date %m+% months(12)),
         post_12_24m = share_overlap(period_start, period_end, cutoff_date %m+% months(12), cutoff_date %m+% months(24)),
         post_24plus = share_overlap(period_start, period_end, cutoff_date %m+% months(24), far_future))

B <- ns(dfq$time_index, df = 2)

dfq <- dfq %>%
  mutate(ns1 = as.numeric(B[, 1]),
         ns2 = as.numeric(B[, 2])) %>%
  dplyr::select(quarter_start,
                attacks,
                time_index,
                intercept_1,
                post_article5_exact,
                post_0_3m,
                post_0_6m,
                post_0_12m,
                post_6_12m,
                post_12_24m,
                post_24plus,
                ns1,
                ns2)

stopifnot(!anyNA(dfq$attacks))
stopifnot(!anyNA(dfq$post_article5_exact))
stopifnot(!anyNA(dfq$post_0_3m))
stopifnot(!anyNA(dfq$post_0_6m))
stopifnot(!anyNA(dfq$post_0_12m))
stopifnot(!anyNA(dfq$post_6_12m))
stopifnot(!anyNA(dfq$post_12_24m))
stopifnot(!anyNA(dfq$post_24plus))
stopifnot(!anyNA(dfq$ns1))
stopifnot(!anyNA(dfq$ns2))

# Consistency check: the four dynamic windows should partition the full
# post-2014 exposure indicator. This verifies that the post-period dummy
# variables are constructed correctly without gaps or overlaps.
dynamic_sum <- with(dfq, post_0_6m + post_6_12m + post_12_24m + post_24plus)
stopifnot(max(abs(dynamic_sum - dfq$post_article5_exact)) < 1e-10)

# ---- Main regression models -------------------------------------------------

# Four specifications fit on the quarterly NATO-target series: a baseline
# model using the full post-2014 exposure indicator, three short-window
# variants restricting exposure to the first 3, 6, or 12 months, and a
# dynamic-window model splitting the post-period into four bins. All
# models use the same prior structure: weakly informative normal priors
# centred at zero for coefficients, with a tighter scale (0.5) on the
# post-2014 coefficient, and an exponential prior for the dispersion
# parameter.

prior_loc_baseline <- c(log(mean(dfq$attacks + 0.5)), 0, 0, 0)
prior_scale_baseline <- c(1, 0.5, 1, 1)

m_baseline <- stan_glm(formula = attacks ~ 0 + intercept_1 + post_article5_exact + ns1 + ns2,
                       data = dfq,
                       family = neg_binomial_2(link = "log"),
                       prior = normal(location = prior_loc_baseline, scale = prior_scale_baseline, autoscale = FALSE),
                       prior_aux = exponential(rate = 1, autoscale = FALSE),
                       chains = 4,
                       iter = 4000,
                       warmup = 2000,
                       seed = 1001,
                       refresh = 0)

prior_loc_short_3 <- c(log(mean(dfq$attacks + 0.5)), 0, 0, 0)
prior_scale_short_3 <- c(1, 0.5, 1, 1)

m_short_3 <- stan_glm(formula = attacks ~ 0 + intercept_1 + post_0_3m + ns1 + ns2,
                      data = dfq,
                      family = neg_binomial_2(link = "log"),
                      prior = normal(location = prior_loc_short_3, scale = prior_scale_short_3, autoscale = FALSE),
                      prior_aux = exponential(rate = 1, autoscale = FALSE),
                      chains = 4,
                      iter = 4000,
                      warmup = 2000,
                      seed = 1002,
                      refresh = 0)

prior_loc_short_6 <- c(log(mean(dfq$attacks + 0.5)), 0, 0, 0)
prior_scale_short_6 <- c(1, 0.5, 1, 1)

m_short_6 <- stan_glm(formula = attacks ~ 0 + intercept_1 + post_0_6m + ns1 + ns2,
                      data = dfq,
                      family = neg_binomial_2(link = "log"),
                      prior = normal(location = prior_loc_short_6, scale = prior_scale_short_6, autoscale = FALSE),
                      prior_aux = exponential(rate = 1, autoscale = FALSE),
                      chains = 4,
                      iter = 4000,
                      warmup = 2000,
                      seed = 1003,
                      refresh = 0)

prior_loc_short_12 <- c(log(mean(dfq$attacks + 0.5)), 0, 0, 0)
prior_scale_short_12 <- c(1, 0.5, 1, 1)

m_short_12 <- stan_glm(formula = attacks ~ 0 + intercept_1 + post_0_12m + ns1 + ns2,
                       data = dfq,
                       family = neg_binomial_2(link = "log"),
                       prior = normal(location = prior_loc_short_12, scale = prior_scale_short_12, autoscale = FALSE),
                       prior_aux = exponential(rate = 1, autoscale = FALSE),
                       chains = 4,
                       iter = 4000,
                       warmup = 2000,
                       seed = 1004,
                       refresh = 0)

# Dynamic model has four post-period exposure variables, so its prior
# vectors are seven elements long. The scale on the four dynamic
# coefficients is set to 0.5, matching the baseline tightness.
prior_loc_dynamic <- c(log(mean(dfq$attacks + 0.5)), 0, 0, 0, 0, 0, 0)
prior_scale_dynamic <- c(1, 0.5, 0.5, 0.5, 0.5, 1, 1)

m_dynamic <- stan_glm(formula = attacks ~ 0 + intercept_1 + post_0_6m + post_6_12m + post_12_24m + post_24plus + ns1 + ns2,
                      data = dfq,
                      family = neg_binomial_2(link = "log"),
                      prior = normal(location = prior_loc_dynamic, scale = prior_scale_dynamic, autoscale = FALSE),
                      prior_aux = exponential(rate = 1, autoscale = FALSE),
                      chains = 4,
                      iter = 4000,
                      warmup = 2000,
                      seed = 1005,
                      refresh = 0)

# ---- Baseline coefficient table ---------------------------------------------

# Extracts the posterior draws of the post-2014 coefficient and constructs
# the main summary table: posterior median, 95 percent credible interval,
# the posterior probability that the coefficient is positive, and the
# implied rate ratio.

draws_baseline_all <- as.matrix(m_baseline)
draws_baseline <- draws_baseline_all
beta_baseline <- draws_baseline[, "post_article5_exact"]
rr_baseline <- exp(beta_baseline)

baseline_table <- tibble(Term = "Post-exposure",
                         Estimate = median(beta_baseline),
                         `95% CrI` = paste0("$[", sprintf("%.3f", quantile(beta_baseline, 0.025)), ", ", sprintf("%.3f", quantile(beta_baseline, 0.975)), "]$"),
                         Pr_beta_gt_0 = mean(beta_baseline > 0),
                         RR = median(rr_baseline))

# kableExtra sometimes emits two consecutive \centering commands when
# combined with kable_styling(latex_options = "HOLD_position"). The gsub
# removes the duplicate so the resulting LaTeX is clean.
baseline_kable <- as.character(
  kable(baseline_table,
        format = "latex",
        booktabs = TRUE,
        digits = 3,
        align = c("l", "r", "c", "r", "r"),
        col.names = c("Term", "Estimate", "95\\% CrI", "$P(\\beta > 0)$", "Rate ratio"),
        caption = "Baseline Bayesian negative binomial model.",
        label = "reg_nb_quarter_stan",
        row.names = FALSE,
        escape = FALSE) %>%
    kable_styling(latex_options = "HOLD_position"))

baseline_kable <- gsub("\\\\centering\\s*\\\\centering", "\\\\centering", baseline_kable, perl = TRUE)
cat(baseline_kable, "\n")

# ---- Baseline dispersion parameter table ------------------------------------

# Dispersion is reported only for the baseline model. All models in this
# script use the same exponential(rate = 1) prior on the dispersion
# parameter and are conditioned on the same data, so the dispersion
# estimates are effectively the same across specifications. Reporting it
# once for baseline avoids repetition.

aux_name <- if ("reciprocal_dispersion" %in% colnames(draws_baseline_all)) {
  "reciprocal_dispersion"
} else if ("shape" %in% colnames(draws_baseline_all)) {
  "shape"
} else if ("aux" %in% colnames(draws_baseline_all)) {
  "aux"
} else {
  stop("Could not find the dispersion parameter in posterior draws.")
}

aux_draws <- draws_baseline_all[, aux_name]

baseline_aux_table <- tibble(Parameter = "Dispersion parameter",
                             Median = median(aux_draws),
                             `95% CrI` = paste0("$[", sprintf("%.3f", quantile(aux_draws, 0.025)), ", ", sprintf("%.3f", quantile(aux_draws, 0.975)), "]$"))

baseline_aux_kable <- as.character(
  kable(baseline_aux_table,
        format = "latex",
        booktabs = TRUE,
        digits = 3,
        align = c("l", "r", "c"),
        col.names = c("Parameter", "Median", "95\\% CrI"),
        caption = "Posterior summary of the dispersion parameter.",
        label = "reg_nb_quarter_aux",
        row.names = FALSE,
        escape = FALSE) %>%
    kable_styling(latex_options = "HOLD_position"))

baseline_aux_kable <- gsub("\\\\centering\\s*\\\\centering", "\\\\centering", baseline_aux_kable, perl = TRUE)
cat(baseline_aux_kable, "\n")

# ---- Baseline model information table ---------------------------------------

# Reported only for the baseline model. All models use the same MCMC
# specification (4 chains, 4000 iterations including 2000 warmup), so
# reporting it once is sufficient.

baseline_modelinfo_table <- tibble(N = nrow(dfq),
                                   Chains = 4,
                                   Iter = 4000,
                                   Warmup = 2000)

baseline_modelinfo_kable <- as.character(
  kable(baseline_modelinfo_table,
        format = "latex",
        booktabs = TRUE,
        digits = 0,
        align = c("r", "r", "r", "r"),
        col.names = c("N", "Chains", "Iter", "Warmup"),
        caption = "Model information for the baseline specification.",
        label = "reg_nb_quarter_modelinfo",
        row.names = FALSE,
        escape = FALSE) %>%
    kable_styling(latex_options = "HOLD_position"))

baseline_modelinfo_kable <- gsub("\\\\centering\\s*\\\\centering", "\\\\centering", baseline_modelinfo_kable, perl = TRUE)
cat(baseline_modelinfo_kable, "\n")

# ---- Short-exposure window coefficient table --------------------------------

# Combines the three short-window models into a single table. Each row
# reports the posterior summary for the relevant coefficient.

draws_short_3 <- as.matrix(m_short_3)
beta_short_3 <- draws_short_3[, "post_0_3m"]
rr_short_3 <- exp(beta_short_3)

short_3_row <- tibble(Term = "Post, 0 to 3 months",
                      Estimate = median(beta_short_3),
                      `95% CrI` = paste0("$[", sprintf("%.3f", quantile(beta_short_3, 0.025)), ", ", sprintf("%.3f", quantile(beta_short_3, 0.975)), "]$"),
                      Pr_beta_gt_0 = mean(beta_short_3 > 0),
                      RR = median(rr_short_3))

draws_short_6 <- as.matrix(m_short_6)
beta_short_6 <- draws_short_6[, "post_0_6m"]
rr_short_6 <- exp(beta_short_6)

short_6_row <- tibble(Term = "Post, 0 to 6 months",
                      Estimate = median(beta_short_6),
                      `95% CrI` = paste0("$[", sprintf("%.3f", quantile(beta_short_6, 0.025)), ", ", sprintf("%.3f", quantile(beta_short_6, 0.975)), "]$"),
                      Pr_beta_gt_0 = mean(beta_short_6 > 0),
                      RR = median(rr_short_6))

draws_short_12 <- as.matrix(m_short_12)
beta_short_12 <- draws_short_12[, "post_0_12m"]
rr_short_12 <- exp(beta_short_12)

short_12_row <- tibble(Term = "Post, 0 to 12 months",
                       Estimate = median(beta_short_12),
                       `95% CrI` = paste0("$[", sprintf("%.3f", quantile(beta_short_12, 0.025)), ", ", sprintf("%.3f", quantile(beta_short_12, 0.975)), "]$"),
                       Pr_beta_gt_0 = mean(beta_short_12 > 0),
                       RR = median(rr_short_12))

short_table <- bind_rows(short_3_row,
                         short_6_row,
                         short_12_row)

short_kable <- as.character(
  kable(short_table,
        format = "latex",
        booktabs = TRUE,
        digits = 3,
        align = c("l", "r", "c", "r", "r"),
        col.names = c("Term", "Estimate", "95\\% CrI", "$P(\\beta > 0)$", "Rate ratio"),
        caption = "Short-exposure window models.",
        label = "short_windows_quarter_stan",
        row.names = FALSE,
        escape = FALSE) %>%
    kable_styling(latex_options = "HOLD_position"))

short_kable <- gsub("\\\\centering\\s*\\\\centering", "\\\\centering", short_kable, perl = TRUE)
cat(short_kable, "\n")

# ---- Dynamic-window coefficient table ---------------------------------------

# Reports the four dynamic-window coefficients from m_dynamic in a single
# table. Together these partition the full post-2014 exposure into bins
# 0-6, 6-12, 12-24, and 24+ months after the recognition.

draws_dynamic <- as.matrix(m_dynamic)

beta_dynamic_0_6 <- draws_dynamic[, "post_0_6m"]
rr_dynamic_0_6 <- exp(beta_dynamic_0_6)

dynamic_0_6_row <- tibble(Term = "Post, 0 to 6 months",
                          Estimate = median(beta_dynamic_0_6),
                          `95% CrI` = paste0("$[", sprintf("%.3f", quantile(beta_dynamic_0_6, 0.025)), ", ", sprintf("%.3f", quantile(beta_dynamic_0_6, 0.975)), "]$"),
                          Pr_beta_gt_0 = mean(beta_dynamic_0_6 > 0),
                          RR = median(rr_dynamic_0_6))

beta_dynamic_6_12 <- draws_dynamic[, "post_6_12m"]
rr_dynamic_6_12 <- exp(beta_dynamic_6_12)

dynamic_6_12_row <- tibble(Term = "Post, 6 to 12 months",
                           Estimate = median(beta_dynamic_6_12),
                           `95% CrI` = paste0("$[", sprintf("%.3f", quantile(beta_dynamic_6_12, 0.025)), ", ", sprintf("%.3f", quantile(beta_dynamic_6_12, 0.975)), "]$"),
                           Pr_beta_gt_0 = mean(beta_dynamic_6_12 > 0),
                           RR = median(rr_dynamic_6_12))

beta_dynamic_12_24 <- draws_dynamic[, "post_12_24m"]
rr_dynamic_12_24 <- exp(beta_dynamic_12_24)

dynamic_12_24_row <- tibble(Term = "Post, 12 to 24 months",
                            Estimate = median(beta_dynamic_12_24),
                            `95% CrI` = paste0("$[", sprintf("%.3f", quantile(beta_dynamic_12_24, 0.025)), ", ", sprintf("%.3f", quantile(beta_dynamic_12_24, 0.975)), "]$"),
                            Pr_beta_gt_0 = mean(beta_dynamic_12_24 > 0),
                            RR = median(rr_dynamic_12_24))

beta_dynamic_24plus <- draws_dynamic[, "post_24plus"]
rr_dynamic_24plus <- exp(beta_dynamic_24plus)

dynamic_24plus_row <- tibble(Term = "Post, 24 months or later",
                             Estimate = median(beta_dynamic_24plus),
                             `95% CrI` = paste0("$[", sprintf("%.3f", quantile(beta_dynamic_24plus, 0.025)), ", ", sprintf("%.3f", quantile(beta_dynamic_24plus, 0.975)), "]$"),
                             Pr_beta_gt_0 = mean(beta_dynamic_24plus > 0),
                             RR = median(rr_dynamic_24plus))

dynamic_table <- bind_rows(dynamic_0_6_row,
                           dynamic_6_12_row,
                           dynamic_12_24_row,
                           dynamic_24plus_row)

dynamic_kable <- as.character(
  kable(dynamic_table,
        format = "latex",
        booktabs = TRUE,
        digits = 3,
        align = c("l", "r", "c", "r", "r"),
        col.names = c("Term", "Estimate", "95\\% CrI", "$P(\\beta > 0)$", "Rate ratio"),
        caption = "Dynamic exposure window model.",
        label = "dynamic_periods_quarter_stan",
        row.names = FALSE,
        escape = FALSE) %>%
    kable_styling(latex_options = "HOLD_position"))

dynamic_kable <- gsub("\\\\centering\\s*\\\\centering", "\\\\centering", dynamic_kable, perl = TRUE)
cat(dynamic_kable, "\n")

# ---- Main comparative regression --------------------------------------------

# Pools NATO-target and non-NATO-target quarterly series in a single
# panel. The model includes group-specific intercepts, group-specific
# spline coefficients, and a group-specific post-period coefficient. The
# interaction term `target_post_*` identifies the NATO-specific post-2014
# effect after controlling for parallel changes in the comparison series.

df_events <- obj$df_comp_events %>%
  mutate(quarter_start = floor_date(interactionstartdate, "quarter"))

df_quarter_group <- df_events %>%
  count(group, target_nato, quarter_start, name = "attacks") %>%
  arrange(group, quarter_start)

all_quarters <- df_quarter %>%
  distinct(quarter_start) %>%
  arrange(quarter_start)

all_groups <- tibble(group = c("NATO target", "Non-NATO target"),
                     target_nato = c(1L, 0L))

dfq_comp <- tidyr::crossing(all_groups, all_quarters) %>%
  left_join(df_quarter_group, by = c("group", "target_nato", "quarter_start")) %>%
  mutate(attacks = ifelse(is.na(attacks), 0L, attacks)) %>%
  arrange(group, quarter_start)

quarter_lookup <- all_quarters %>%
  mutate(time_index = row_number())

dfq_comp <- dfq_comp %>%
  left_join(quarter_lookup, by = "quarter_start") %>%
  mutate(period_start = quarter_start,
         period_end = quarter_start %m+% months(3),
         intercept_1 = 1,
         post_article5_exact = share_overlap(period_start, period_end, cutoff_date, far_future),
         post_0_3m = share_overlap(period_start, period_end, cutoff_date, cutoff_date %m+% months(3)),
         post_0_6m = share_overlap(period_start, period_end, cutoff_date, cutoff_date %m+% months(6)),
         post_0_12m = share_overlap(period_start, period_end, cutoff_date, cutoff_date %m+% months(12)))

B_comp <- ns(dfq_comp$time_index, df = 2)

dfq_comp <- dfq_comp %>%
  mutate(ns1 = as.numeric(B_comp[, 1]),
         ns2 = as.numeric(B_comp[, 2]),
         target_ns1 = target_nato * ns1,
         target_ns2 = target_nato * ns2,
         target_post_article5_exact = target_nato * post_article5_exact,
         target_post_0_3m = target_nato * post_0_3m,
         target_post_0_6m = target_nato * post_0_6m,
         target_post_0_12m = target_nato * post_0_12m) %>%
  dplyr::select(group,
                target_nato,
                quarter_start,
                attacks,
                time_index,
                intercept_1,
                post_article5_exact,
                post_0_3m,
                post_0_6m,
                post_0_12m,
                ns1,
                ns2,
                target_ns1,
                target_ns2,
                target_post_article5_exact,
                target_post_0_3m,
                target_post_0_6m,
                target_post_0_12m)

stopifnot(!anyNA(dfq_comp$attacks))
stopifnot(!anyNA(dfq_comp$target_nato))
stopifnot(!anyNA(dfq_comp$post_article5_exact))
stopifnot(!anyNA(dfq_comp$post_0_3m))
stopifnot(!anyNA(dfq_comp$post_0_6m))
stopifnot(!anyNA(dfq_comp$post_0_12m))
stopifnot(!anyNA(dfq_comp$ns1))
stopifnot(!anyNA(dfq_comp$ns2))
stopifnot(!anyNA(dfq_comp$target_ns1))
stopifnot(!anyNA(dfq_comp$target_ns2))
stopifnot(!anyNA(dfq_comp$target_post_article5_exact))
stopifnot(!anyNA(dfq_comp$target_post_0_3m))
stopifnot(!anyNA(dfq_comp$target_post_0_6m))
stopifnot(!anyNA(dfq_comp$target_post_0_12m))

# Comparative models. Prior vectors are eight elements long: intercept,
# group indicator, two spline coefficients, two group-by-spline
# interactions, post-period coefficient, and group-by-post interaction.
# The group-by-post interaction receives the tighter 0.5 scale because
# it is the NATO-specific effect of interest.

prior_loc_comp_baseline <- c(log(mean(dfq_comp$attacks + 0.5)), 0, 0, 0, 0, 0, 0, 0)
prior_scale_comp_baseline <- c(1, 0.5, 1, 1, 1, 1, 0.5, 0.5)

m_comp_baseline <- stan_glm(formula = attacks ~ 0 + intercept_1 + target_nato + ns1 + ns2 + target_ns1 + target_ns2 + post_article5_exact + target_post_article5_exact,
                            data = dfq_comp,
                            family = neg_binomial_2(link = "log"),
                            prior = normal(location = prior_loc_comp_baseline, scale = prior_scale_comp_baseline, autoscale = FALSE),
                            prior_aux = exponential(rate = 1, autoscale = FALSE),
                            chains = 4,
                            iter = 4000,
                            warmup = 2000,
                            seed = 2001,
                            refresh = 0)

prior_loc_comp_short_3 <- c(log(mean(dfq_comp$attacks + 0.5)), 0, 0, 0, 0, 0, 0, 0)
prior_scale_comp_short_3 <- c(1, 0.5, 1, 1, 1, 1, 0.5, 0.5)

m_comp_short_3 <- stan_glm(formula = attacks ~ 0 + intercept_1 + target_nato + ns1 + ns2 + target_ns1 + target_ns2 + post_0_3m + target_post_0_3m,
                           data = dfq_comp,
                           family = neg_binomial_2(link = "log"),
                           prior = normal(location = prior_loc_comp_short_3, scale = prior_scale_comp_short_3, autoscale = FALSE),
                           prior_aux = exponential(rate = 1, autoscale = FALSE),
                           chains = 4,
                           iter = 4000,
                           warmup = 2000,
                           seed = 2002,
                           refresh = 0)

prior_loc_comp_short_6 <- c(log(mean(dfq_comp$attacks + 0.5)), 0, 0, 0, 0, 0, 0, 0)
prior_scale_comp_short_6 <- c(1, 0.5, 1, 1, 1, 1, 0.5, 0.5)

m_comp_short_6 <- stan_glm(formula = attacks ~ 0 + intercept_1 + target_nato + ns1 + ns2 + target_ns1 + target_ns2 + post_0_6m + target_post_0_6m,
                           data = dfq_comp,
                           family = neg_binomial_2(link = "log"),
                           prior = normal(location = prior_loc_comp_short_6, scale = prior_scale_comp_short_6, autoscale = FALSE),
                           prior_aux = exponential(rate = 1, autoscale = FALSE),
                           chains = 4,
                           iter = 4000,
                           warmup = 2000,
                           seed = 2003,
                           refresh = 0)

prior_loc_comp_short_12 <- c(log(mean(dfq_comp$attacks + 0.5)), 0, 0, 0, 0, 0, 0, 0)
prior_scale_comp_short_12 <- c(1, 0.5, 1, 1, 1, 1, 0.5, 0.5)

m_comp_short_12 <- stan_glm(formula = attacks ~ 0 + intercept_1 + target_nato + ns1 + ns2 + target_ns1 + target_ns2 + post_0_12m + target_post_0_12m,
                            data = dfq_comp,
                            family = neg_binomial_2(link = "log"),
                            prior = normal(location = prior_loc_comp_short_12, scale = prior_scale_comp_short_12, autoscale = FALSE),
                            prior_aux = exponential(rate = 1, autoscale = FALSE),
                            chains = 4,
                            iter = 4000,
                            warmup = 2000,
                            seed = 2004,
                            refresh = 0)

# ---- Comparative coefficient table ------------------------------------------

# Reports the NATO-specific post-period coefficient (the group-by-post
# interaction) from each of the four comparative specifications.

draws_comp_baseline <- as.matrix(m_comp_baseline)
beta_comp_baseline <- draws_comp_baseline[, "target_post_article5_exact"]
rr_comp_baseline <- exp(beta_comp_baseline)

comp_baseline_row <- tibble(Term = "NATO-specific post-exposure",
                            Estimate = median(beta_comp_baseline),
                            `95% CrI` = paste0("$[", sprintf("%.3f", quantile(beta_comp_baseline, 0.025)), ", ", sprintf("%.3f", quantile(beta_comp_baseline, 0.975)), "]$"),
                            Pr_beta_gt_0 = mean(beta_comp_baseline > 0),
                            RR = median(rr_comp_baseline))

draws_comp_short_3 <- as.matrix(m_comp_short_3)
beta_comp_short_3 <- draws_comp_short_3[, "target_post_0_3m"]
rr_comp_short_3 <- exp(beta_comp_short_3)

comp_short_3_row <- tibble(Term = "NATO-specific post, 0 to 3 months",
                           Estimate = median(beta_comp_short_3),
                           `95% CrI` = paste0("$[", sprintf("%.3f", quantile(beta_comp_short_3, 0.025)), ", ", sprintf("%.3f", quantile(beta_comp_short_3, 0.975)), "]$"),
                           Pr_beta_gt_0 = mean(beta_comp_short_3 > 0),
                           RR = median(rr_comp_short_3))

draws_comp_short_6 <- as.matrix(m_comp_short_6)
beta_comp_short_6 <- draws_comp_short_6[, "target_post_0_6m"]
rr_comp_short_6 <- exp(beta_comp_short_6)

comp_short_6_row <- tibble(Term = "NATO-specific post, 0 to 6 months",
                           Estimate = median(beta_comp_short_6),
                           `95% CrI` = paste0("$[", sprintf("%.3f", quantile(beta_comp_short_6, 0.025)), ", ", sprintf("%.3f", quantile(beta_comp_short_6, 0.975)), "]$"),
                           Pr_beta_gt_0 = mean(beta_comp_short_6 > 0),
                           RR = median(rr_comp_short_6))

draws_comp_short_12 <- as.matrix(m_comp_short_12)
beta_comp_short_12 <- draws_comp_short_12[, "target_post_0_12m"]
rr_comp_short_12 <- exp(beta_comp_short_12)

comp_short_12_row <- tibble(Term = "NATO-specific post, 0 to 12 months",
                            Estimate = median(beta_comp_short_12),
                            `95% CrI` = paste0("$[", sprintf("%.3f", quantile(beta_comp_short_12, 0.025)), ", ", sprintf("%.3f", quantile(beta_comp_short_12, 0.975)), "]$"),
                            Pr_beta_gt_0 = mean(beta_comp_short_12 > 0),
                            RR = median(rr_comp_short_12))

comp_combined_table <- bind_rows(comp_baseline_row,
                                 comp_short_3_row,
                                 comp_short_6_row,
                                 comp_short_12_row)

comp_kable <- as.character(
  kable(comp_combined_table,
        format = "latex",
        booktabs = TRUE,
        digits = 3,
        align = c("l", "r", "c", "r", "r"),
        col.names = c("Term", "Estimate", "95\\% CrI", "$P(\\beta > 0)$", "Rate ratio"),
        caption = "Comparative regression models.",
        label = "comp_models_quarter_stan",
        row.names = FALSE,
        escape = FALSE) %>%
    kable_styling(latex_options = "HOLD_position"))

comp_kable <- gsub("\\\\centering\\s*\\\\centering", "\\\\centering", comp_kable, perl = TRUE)
cat(comp_kable, "\n")

# ---- Appendix regressions: Extended short-exposure windows ------------------

# Extends the short-window analysis to 18 and 24 months for the appendix.
# Uses the same prior structure as the main short-window models.

dfq_short_appendix <- df_quarter %>%
  mutate(period_start = quarter_start,
         period_end = quarter_start %m+% months(3),
         intercept_1 = 1,
         post_0_18m = share_overlap(period_start, period_end, cutoff_date, cutoff_date %m+% months(18)),
         post_0_24m = share_overlap(period_start, period_end, cutoff_date, cutoff_date %m+% months(24)))

B_short_appendix <- ns(dfq_short_appendix$time_index, df = 2)

dfq_short_appendix <- dfq_short_appendix %>%
  mutate(ns1 = as.numeric(B_short_appendix[, 1]),
         ns2 = as.numeric(B_short_appendix[, 2])) %>%
  dplyr::select(quarter_start,
                attacks,
                time_index,
                intercept_1,
                post_0_18m,
                post_0_24m,
                ns1,
                ns2)

stopifnot(!anyNA(dfq_short_appendix$attacks))
stopifnot(!anyNA(dfq_short_appendix$post_0_18m))
stopifnot(!anyNA(dfq_short_appendix$post_0_24m))
stopifnot(!anyNA(dfq_short_appendix$ns1))
stopifnot(!anyNA(dfq_short_appendix$ns2))

prior_loc_short_18 <- c(log(mean(dfq_short_appendix$attacks + 0.5)), 0, 0, 0)
prior_scale_short_18 <- c(1, 0.5, 1, 1)

m_short_18 <- stan_glm(formula = attacks ~ 0 + intercept_1 + post_0_18m + ns1 + ns2,
                       data = dfq_short_appendix,
                       family = neg_binomial_2(link = "log"),
                       prior = normal(location = prior_loc_short_18, scale = prior_scale_short_18, autoscale = FALSE),
                       prior_aux = exponential(rate = 1, autoscale = FALSE),
                       chains = 4,
                       iter = 4000,
                       warmup = 2000,
                       seed = 1006,
                       refresh = 0)

prior_loc_short_24 <- c(log(mean(dfq_short_appendix$attacks + 0.5)), 0, 0, 0)
prior_scale_short_24 <- c(1, 0.5, 1, 1)

m_short_24 <- stan_glm(formula = attacks ~ 0 + intercept_1 + post_0_24m + ns1 + ns2,
                       data = dfq_short_appendix,
                       family = neg_binomial_2(link = "log"),
                       prior = normal(location = prior_loc_short_24, scale = prior_scale_short_24, autoscale = FALSE),
                       prior_aux = exponential(rate = 1, autoscale = FALSE),
                       chains = 4,
                       iter = 4000,
                       warmup = 2000,
                       seed = 1007,
                       refresh = 0)

draws_short_18 <- as.matrix(m_short_18)
beta_short_18 <- draws_short_18[, "post_0_18m"]
rr_short_18 <- exp(beta_short_18)

short_18_row <- tibble(Term = "Post, 0 to 18 months",
                       Estimate = median(beta_short_18),
                       `95% CrI` = paste0("$[", sprintf("%.3f", quantile(beta_short_18, 0.025)), ", ", sprintf("%.3f", quantile(beta_short_18, 0.975)), "]$"),
                       Pr_beta_gt_0 = mean(beta_short_18 > 0),
                       RR = median(rr_short_18))

draws_short_24 <- as.matrix(m_short_24)
beta_short_24 <- draws_short_24[, "post_0_24m"]
rr_short_24 <- exp(beta_short_24)

short_24_row <- tibble(Term = "Post, 0 to 24 months",
                       Estimate = median(beta_short_24),
                       `95% CrI` = paste0("$[", sprintf("%.3f", quantile(beta_short_24, 0.025)), ", ", sprintf("%.3f", quantile(beta_short_24, 0.975)), "]$"),
                       Pr_beta_gt_0 = mean(beta_short_24 > 0),
                       RR = median(rr_short_24))

short_appendix_table <- bind_rows(short_18_row,
                                  short_24_row)

short_appendix_kable <- as.character(
  kable(short_appendix_table,
        format = "latex",
        booktabs = TRUE,
        digits = 3,
        align = c("l", "r", "c", "r", "r"),
        col.names = c("Term", "Estimate", "95\\% CrI", "$P(\\beta > 0)$", "Rate ratio"),
        caption = "Extended short-exposure window models.",
        label = "short_windows_extended_quarter_stan",
        row.names = FALSE,
        escape = FALSE) %>%
    kable_styling(latex_options = "HOLD_position"))

short_appendix_kable <- gsub("\\\\centering\\s*\\\\centering", "\\\\centering", short_appendix_kable, perl = TRUE)
cat(short_appendix_kable, "\n")

# ---- Appendix regressions: Pre-recognition exposure windows -----------------

# Mirror-image of the short-window analysis using pre-recognition windows
# (3, 6, 12, 18, 24 months before 5 September 2014). Tests whether a
# placebo cut-off shifted earlier in time produces apparent effects, which
# would suggest pre-existing trends rather than a true post-recognition
# shift.

dfq_pre <- df_quarter %>%
  mutate(period_start = quarter_start,
         period_end = quarter_start %m+% months(3),
         intercept_1 = 1,
         pre_0_24m = share_overlap(period_start, period_end, cutoff_date %m-% months(24), cutoff_date),
         pre_0_18m = share_overlap(period_start, period_end, cutoff_date %m-% months(18), cutoff_date),
         pre_0_12m = share_overlap(period_start, period_end, cutoff_date %m-% months(12), cutoff_date),
         pre_0_6m = share_overlap(period_start, period_end, cutoff_date %m-% months(6), cutoff_date),
         pre_0_3m = share_overlap(period_start, period_end, cutoff_date %m-% months(3), cutoff_date))

B_pre <- ns(dfq_pre$time_index, df = 2)

dfq_pre <- dfq_pre %>%
  mutate(ns1 = as.numeric(B_pre[, 1]),
         ns2 = as.numeric(B_pre[, 2])) %>%
  dplyr::select(quarter_start,
                attacks,
                time_index,
                intercept_1,
                pre_0_24m,
                pre_0_18m,
                pre_0_12m,
                pre_0_6m,
                pre_0_3m,
                ns1,
                ns2)

stopifnot(!anyNA(dfq_pre$attacks))
stopifnot(!anyNA(dfq_pre$pre_0_24m))
stopifnot(!anyNA(dfq_pre$pre_0_18m))
stopifnot(!anyNA(dfq_pre$pre_0_12m))
stopifnot(!anyNA(dfq_pre$pre_0_6m))
stopifnot(!anyNA(dfq_pre$pre_0_3m))
stopifnot(!anyNA(dfq_pre$ns1))
stopifnot(!anyNA(dfq_pre$ns2))

prior_loc_pre_24 <- c(log(mean(dfq_pre$attacks + 0.5)), 0, 0, 0)
prior_scale_pre_24 <- c(1, 0.5, 1, 1)

m_pre_24 <- stan_glm(formula = attacks ~ 0 + intercept_1 + pre_0_24m + ns1 + ns2,
                     data = dfq_pre,
                     family = neg_binomial_2(link = "log"),
                     prior = normal(location = prior_loc_pre_24, scale = prior_scale_pre_24, autoscale = FALSE),
                     prior_aux = exponential(rate = 1, autoscale = FALSE),
                     chains = 4,
                     iter = 4000,
                     warmup = 2000,
                     seed = 1101,
                     refresh = 0)

prior_loc_pre_18 <- c(log(mean(dfq_pre$attacks + 0.5)), 0, 0, 0)
prior_scale_pre_18 <- c(1, 0.5, 1, 1)

m_pre_18 <- stan_glm(formula = attacks ~ 0 + intercept_1 + pre_0_18m + ns1 + ns2,
                     data = dfq_pre,
                     family = neg_binomial_2(link = "log"),
                     prior = normal(location = prior_loc_pre_18, scale = prior_scale_pre_18, autoscale = FALSE),
                     prior_aux = exponential(rate = 1, autoscale = FALSE),
                     chains = 4,
                     iter = 4000,
                     warmup = 2000,
                     seed = 1102,
                     refresh = 0)

prior_loc_pre_12 <- c(log(mean(dfq_pre$attacks + 0.5)), 0, 0, 0)
prior_scale_pre_12 <- c(1, 0.5, 1, 1)

m_pre_12 <- stan_glm(formula = attacks ~ 0 + intercept_1 + pre_0_12m + ns1 + ns2,
                     data = dfq_pre,
                     family = neg_binomial_2(link = "log"),
                     prior = normal(location = prior_loc_pre_12, scale = prior_scale_pre_12, autoscale = FALSE),
                     prior_aux = exponential(rate = 1, autoscale = FALSE),
                     chains = 4,
                     iter = 4000,
                     warmup = 2000,
                     seed = 1103,
                     refresh = 0)

prior_loc_pre_6 <- c(log(mean(dfq_pre$attacks + 0.5)), 0, 0, 0)
prior_scale_pre_6 <- c(1, 0.5, 1, 1)

m_pre_6 <- stan_glm(formula = attacks ~ 0 + intercept_1 + pre_0_6m + ns1 + ns2,
                    data = dfq_pre,
                    family = neg_binomial_2(link = "log"),
                    prior = normal(location = prior_loc_pre_6, scale = prior_scale_pre_6, autoscale = FALSE),
                    prior_aux = exponential(rate = 1, autoscale = FALSE),
                    chains = 4,
                    iter = 4000,
                    warmup = 2000,
                    seed = 1104,
                    refresh = 0)

prior_loc_pre_3 <- c(log(mean(dfq_pre$attacks + 0.5)), 0, 0, 0)
prior_scale_pre_3 <- c(1, 0.5, 1, 1)

m_pre_3 <- stan_glm(formula = attacks ~ 0 + intercept_1 + pre_0_3m + ns1 + ns2,
                    data = dfq_pre,
                    family = neg_binomial_2(link = "log"),
                    prior = normal(location = prior_loc_pre_3, scale = prior_scale_pre_3, autoscale = FALSE),
                    prior_aux = exponential(rate = 1, autoscale = FALSE),
                    chains = 4,
                    iter = 4000,
                    warmup = 2000,
                    seed = 1105,
                    refresh = 0)

draws_pre_24 <- as.matrix(m_pre_24)
beta_pre_24 <- draws_pre_24[, "pre_0_24m"]
rr_pre_24 <- exp(beta_pre_24)

pre_24_row <- tibble(Term = "Pre, 24 to 0 months",
                     Estimate = median(beta_pre_24),
                     `95% CrI` = paste0("$[", sprintf("%.3f", quantile(beta_pre_24, 0.025)), ", ", sprintf("%.3f", quantile(beta_pre_24, 0.975)), "]$"),
                     Pr_beta_gt_0 = mean(beta_pre_24 > 0),
                     RR = median(rr_pre_24))

draws_pre_18 <- as.matrix(m_pre_18)
beta_pre_18 <- draws_pre_18[, "pre_0_18m"]
rr_pre_18 <- exp(beta_pre_18)

pre_18_row <- tibble(Term = "Pre, 18 to 0 months",
                     Estimate = median(beta_pre_18),
                     `95% CrI` = paste0("$[", sprintf("%.3f", quantile(beta_pre_18, 0.025)), ", ", sprintf("%.3f", quantile(beta_pre_18, 0.975)), "]$"),
                     Pr_beta_gt_0 = mean(beta_pre_18 > 0),
                     RR = median(rr_pre_18))

draws_pre_12 <- as.matrix(m_pre_12)
beta_pre_12 <- draws_pre_12[, "pre_0_12m"]
rr_pre_12 <- exp(beta_pre_12)

pre_12_row <- tibble(Term = "Pre, 12 to 0 months",
                     Estimate = median(beta_pre_12),
                     `95% CrI` = paste0("$[", sprintf("%.3f", quantile(beta_pre_12, 0.025)), ", ", sprintf("%.3f", quantile(beta_pre_12, 0.975)), "]$"),
                     Pr_beta_gt_0 = mean(beta_pre_12 > 0),
                     RR = median(rr_pre_12))

draws_pre_6 <- as.matrix(m_pre_6)
beta_pre_6 <- draws_pre_6[, "pre_0_6m"]
rr_pre_6 <- exp(beta_pre_6)

pre_6_row <- tibble(Term = "Pre, 6 to 0 months",
                    Estimate = median(beta_pre_6),
                    `95% CrI` = paste0("$[", sprintf("%.3f", quantile(beta_pre_6, 0.025)), ", ", sprintf("%.3f", quantile(beta_pre_6, 0.975)), "]$"),
                    Pr_beta_gt_0 = mean(beta_pre_6 > 0),
                    RR = median(rr_pre_6))

draws_pre_3 <- as.matrix(m_pre_3)
beta_pre_3 <- draws_pre_3[, "pre_0_3m"]
rr_pre_3 <- exp(beta_pre_3)

pre_3_row <- tibble(Term = "Pre, 3 to 0 months",
                    Estimate = median(beta_pre_3),
                    `95% CrI` = paste0("$[", sprintf("%.3f", quantile(beta_pre_3, 0.025)), ", ", sprintf("%.3f", quantile(beta_pre_3, 0.975)), "]$"),
                    Pr_beta_gt_0 = mean(beta_pre_3 > 0),
                    RR = median(rr_pre_3))

pre_window_table <- bind_rows(pre_24_row,
                              pre_18_row,
                              pre_12_row,
                              pre_6_row,
                              pre_3_row)

pre_window_kable <- as.character(
  kable(pre_window_table,
        format = "latex",
        booktabs = TRUE,
        digits = 3,
        align = c("l", "r", "c", "r", "r"),
        col.names = c("Term", "Estimate", "95\\% CrI", "$P(\\beta > 0)$", "Rate ratio"),
        caption = "Pre-recognition exposure window models.",
        label = "pre_windows_quarter_stan",
        row.names = FALSE,
        escape = FALSE) %>%
    kable_styling(latex_options = "HOLD_position"))

pre_window_kable <- gsub("\\\\centering\\s*\\\\centering", "\\\\centering", pre_window_kable, perl = TRUE)
cat(pre_window_kable, "\n")