# =============================================================================
# modeltest.R
#
# Robustness and diagnostic checks for the count regression model. Compares
# Poisson against negative binomial, the negative binomial against zero-
# inflated alternatives, spline degrees of freedom, linear against spline
# trends, post-cutoff trend breaks, residual autocorrelation, temporal
# granularity, and prior sensitivity for the quarterly Bayesian model. All
# tables are printed to the console in LaTeX form; figures are saved to
# Figures/.
# =============================================================================

library(dplyr)
library(tibble)
library(knitr)
library(kableExtra)
library(modelsummary)
library(MASS)
library(pscl)
library(splines)
library(lubridate)
library(ggplot2)
library(rethinking)

# ---- Palette and base theme -------------------------------------------------

# Palette and base theme are duplicated across descriptive.R, modeltest.R, and
# dag.R so that each script can be run independently of the others. Any change
# here should be mirrored in those files.

palette_full <- c(
  "#6C8FD4",
  "#A78FD6",
  "#7FC8B2",
  "#B8D986",
  "#E6A86C",
  "#D89595",
  "#C7B486")

color_sim <- palette_full[1]
color_sim_dark <- "#3D5FA8"
color_obs <- palette_full[5]

base_theme <- theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(color = "black", face = "bold", size = 12),
    axis.title = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    legend.title = element_text(color = "black"),
    legend.text = element_text(color = "black"),
    strip.text = element_text(color = "black", face = "bold"),
    strip.background = element_blank())

dir.create("Figures", showWarnings = FALSE, recursive = TRUE)

# ---- Read derived series ----------------------------------------------------

# Reads aggregated count series produced by clean_data.R.
in_path <- "data/derived_series.rds"
obj <- readRDS(in_path)

df_month <- obj$df_month
df_quarter <- obj$df_quarter
df_halfyear <- obj$df_halfyear
df_year <- obj$df_year

cutoff_date <- as.Date("2014-09-05")

# ---- Helper: share of period overlap ----------------------------------------

# Computes the share of a period that falls within a target window. Used to
# build the post-2014 exposure indicator: for each period, the share of days
# falling on or after the cutoff date. Periods entirely after the cutoff
# return 1, periods entirely before return 0, and the period containing the
# cutoff returns a fractional value. This is duplicated in regression.R so
# that each script runs independently.

share_overlap <- function(start, end, win_start, win_end) {
  den <- as.numeric(end - start)
  num <- pmax(0, as.numeric(pmin(end, win_end) - pmax(start, win_start)))
  num / den
}

# ---- Tabular and column mappings --------------------------------------------

# Coefficient label mapping for modelsummary tables. Splines of different
# degrees of freedom are rendered with a unified "Time spline k" naming so
# that tables across specifications remain comparable.

coef_map_main <- c(
  "(Intercept)" = "Intercept",
  "post_article5_exact" = "Post Article 5",
  "time_index" = "Time index",
  "ns(time_index, df = 2)1" = "Time spline 1",
  "ns(time_index, df = 2)2" = "Time spline 2",
  "ns(time_index, df = 3)1" = "Time spline 1",
  "ns(time_index, df = 3)2" = "Time spline 2",
  "ns(time_index, df = 3)3" = "Time spline 3",
  "ns(time_index, df = 4)1" = "Time spline 1",
  "ns(time_index, df = 4)2" = "Time spline 2",
  "ns(time_index, df = 4)3" = "Time spline 3",
  "ns(time_index, df = 4)4" = "Time spline 4")

gof_map_main <- tribble(
  ~raw,     ~clean,            ~fmt,
  "nobs",   "N",               0,
  "AIC",    "AIC",             2,
  "BIC",    "BIC",             2,
  "logLik", "Log likelihood",  2)

# ---- Monthly working data ---------------------------------------------------

# Adds period start and end columns and constructs the post-2014 exposure
# indicator as the share of each month falling on or after the cutoff date.
# This gives a fractional value for the period containing the recognition
# and 0 or 1 elsewhere.

dfm <- df_month %>%
  mutate(period_start = month,
         period_end = month %m+% months(1),
         post_article5_exact = share_overlap(period_start,
                                             period_end,
                                             cutoff_date,
                                             as.Date("2100-01-01")))

# ---- Overdispersion diagnostics ---------------------------------------------

# Checks whether monthly counts show overdispersion (variance greater than the
# mean). The Pearson dispersion statistic from a Poisson fit indicates whether
# a negative binomial model is required. Values substantially above 1 suggest
# overdispersion.

mean_attacks <- mean(dfm$attacks, na.rm = TRUE)
var_attacks <- var(dfm$attacks, na.rm = TRUE)
vm_ratio <- var_attacks / mean_attacks
zero_share <- mean(dfm$attacks == 0, na.rm = TRUE)

m_pois_diag <- glm(attacks ~ post_article5_exact + ns(time_index, df = 2),
                   family = poisson(link = "log"),
                   data = dfm)

pearson_chi2 <- sum(residuals(m_pois_diag, type = "pearson")^2, na.rm = TRUE)
dispersion_stat <- pearson_chi2 / m_pois_diag$df.residual

overdisp_table <- data.frame(mean = mean_attacks,
                             variance = var_attacks,
                             variance_to_mean = vm_ratio,
                             share_zero = zero_share,
                             pearson_dispersion = dispersion_stat)

kable(overdisp_table,
      format = "latex",
      booktabs = TRUE,
      digits = 3,
      col.names = c("Mean", "Variance", "Variance to mean", "Share zero", "Pearson dispersion"),
      caption = "Summary diagnostics for monthly attack counts.",
      label = "overdispersion_month",
      row.names = FALSE,
      escape = FALSE) %>%
  kable_styling(latex_options = "HOLD_position")

# ---- Poisson against negative binomial, monthly, spline df 2 ----------------

# Compares the two main count distributions. AIC favours the model with lower
# value. Confirms whether the overdispersion observed above is large enough
# to prefer the negative binomial.

m_pois_df2 <- glm(attacks ~ post_article5_exact + ns(time_index, df = 2),
                  family = poisson(link = "log"),
                  data = dfm)

m_nb_df2 <- glm.nb(attacks ~ post_article5_exact + ns(time_index, df = 2),
                   data = dfm)

modelsummary(list("Poisson, spline df 2" = m_pois_df2,
                  "Negative binomial, spline df 2" = m_nb_df2),
             coef_map = coef_map_main,
             gof_map = gof_map_main,
             statistic = "({std.error})",
             stars = TRUE,
             output = "kableExtra",
             title = "Monthly count models with spline df 2.",
             escape = FALSE) %>%
  kable_styling(latex_options = "HOLD_position")

aic_pois_vs_nb <- AIC(m_pois_df2, m_nb_df2) %>%
  rownames_to_column(var = "model") %>%
  mutate(model = recode(model,
                        m_pois_df2 = "Poisson, spline df 2",
                        m_nb_df2 = "Negative binomial, spline df 2"),
         delta_aic = AIC - min(AIC)) %>%
  arrange(AIC)

# ---- Negative binomial against zero-inflated, monthly, spline df 2 ----------

# Tests whether the share of zero-count months is larger than the negative
# binomial alone can accommodate. A substantially lower AIC for the zero-
# inflated model would suggest a separate zero-generating process.

m_zinb_df2 <- zeroinfl(attacks ~ post_article5_exact + ns(time_index, df = 2) | 1,
                       data = dfm,
                       dist = "negbin")

aic_nb_vs_zinb_df2 <- AIC(m_nb_df2, m_zinb_df2) %>%
  rownames_to_column(var = "model") %>%
  mutate(model = recode(model,
                        m_nb_df2 = "Negative binomial, spline df 2",
                        m_zinb_df2 = "Zero inflated negative binomial, spline df 2"),
         delta_aic = AIC - min(AIC)) %>%
  arrange(AIC)

# ---- Negative binomial against zero-inflated, monthly, spline df 3 ----------

# Same zero-inflation check but with a more flexible time trend. Verifies
# that the conclusion from the df 2 comparison is not driven by inadequate
# trend modelling.

m_nb_df3 <- glm.nb(attacks ~ post_article5_exact + ns(time_index, df = 3),
                   data = dfm)

m_zinb_df3 <- zeroinfl(attacks ~ post_article5_exact + ns(time_index, df = 3) | 1,
                       data = dfm,
                       dist = "negbin")

aic_nb_vs_zinb_df3 <- AIC(m_nb_df3, m_zinb_df3) %>%
  rownames_to_column(var = "model") %>%
  mutate(model = recode(model,
                        m_nb_df3 = "Negative binomial, spline df 3",
                        m_zinb_df3 = "Zero inflated negative binomial, spline df 3"),
         delta_aic = AIC - min(AIC)) %>%
  arrange(AIC)

# ---- Spline flexibility, monthly negative binomial --------------------------

# Compares spline degrees of freedom 2, 3, and 4. Df 2 is the minimum for a
# non-linear spline, df 4 approaches the upper limit of useful flexibility
# given the length of the monthly series, and df 3 covers the intermediate
# level. The aim is to identify the simplest specification that captures the
# underlying time trend.

m_nb_df4 <- glm.nb(attacks ~ post_article5_exact + ns(time_index, df = 4),
                   data = dfm)

aic_spline_df <- data.frame(spline_df = c(2, 3, 4),
                            df_model = c(attr(logLik(m_nb_df2), "df"),
                                         attr(logLik(m_nb_df3), "df"),
                                         attr(logLik(m_nb_df4), "df")),
                            aic = c(AIC(m_nb_df2), AIC(m_nb_df3), AIC(m_nb_df4))) %>%
  arrange(aic) %>%
  mutate(delta_aic = aic - min(aic))

# ---- Linear against spline trend, monthly negative binomial -----------------

# Checks whether a linear time trend is sufficient or whether the spline's
# extra flexibility is justified.

m_nb_linear <- glm.nb(attacks ~ post_article5_exact + time_index,
                      data = dfm)

modelsummary(list("Negative binomial, linear trend" = m_nb_linear,
                  "Negative binomial, spline df 2" = m_nb_df2),
             coef_map = coef_map_main,
             gof_map = gof_map_main,
             statistic = "({std.error})",
             stars = TRUE,
             output = "kableExtra",
             title = "Monthly negative binomial models with linear and spline time trends.",
             escape = FALSE) %>%
  kable_styling(latex_options = "HOLD_position")

aic_linear_vs_df2 <- AIC(m_nb_linear, m_nb_df2) %>%
  rownames_to_column(var = "model") %>%
  mutate(model = recode(model,
                        m_nb_linear = "Negative binomial, linear trend",
                        m_nb_df2 = "Negative binomial, spline df 2"),
         delta_aic = AIC - min(AIC)) %>%
  arrange(AIC)

# ---- Residual autocorrelation, monthly negative binomial --------------------

# Tests for residual autocorrelation. AR(0) is the null of no autocorrelation;
# AR(1) and AR(2) are the simplest alternatives. Higher orders showed no
# improvement and were considered unnecessarily complex for the length of the
# series. Ljung-Box tests are run at lags 12 and 24 to check for seasonal
# autocorrelation in the monthly residuals, corresponding to one-year and
# two-year cycles.

res_nb_df2 <- residuals(m_nb_df2, type = "pearson")

ar0_model <- arima(res_nb_df2, order = c(0, 0, 0), include.mean = FALSE)
ar1_model <- arima(res_nb_df2, order = c(1, 0, 0), include.mean = FALSE)
ar2_model <- arima(res_nb_df2, order = c(2, 0, 0), include.mean = FALSE)

ar_aic <- tibble(ar_order = c(0, 1, 2),
                 aic = c(AIC(ar0_model), AIC(ar1_model), AIC(ar2_model))) %>%
  arrange(aic) %>%
  mutate(delta_aic = aic - min(aic))

kable(ar_aic,
      format = "latex",
      booktabs = TRUE,
      digits = 3,
      col.names = c("AR order", "AIC", "Delta AIC"),
      caption = "AIC comparison for AR(p) models fitted to Pearson residuals from the monthly negative binomial model with spline df 2.",
      label = "aic_ar_residuals_month",
      row.names = FALSE,
      escape = FALSE) %>%
  kable_styling(latex_options = "HOLD_position")

lb12 <- Box.test(res_nb_df2, lag = 12, type = "Ljung-Box")
lb24 <- Box.test(res_nb_df2, lag = 24, type = "Ljung-Box")

lb_table <- data.frame(lag = c(12, 24),
                       statistic = c(as.numeric(lb12$statistic), as.numeric(lb24$statistic)),
                       p_value = c(lb12$p.value, lb24$p.value))

kable(lb_table,
      format = "latex",
      booktabs = TRUE,
      digits = 3,
      col.names = c("Lag", "Test statistic", "p value"),
      caption = "Ljung Box tests for residual autocorrelation in Pearson residuals from the monthly negative binomial model with spline df 2.",
      label = "ljung_box_month",
      row.names = FALSE,
      escape = FALSE) %>%
  kable_styling(latex_options = "HOLD_position")

# ---- Post-cutoff trend break, monthly ---------------------------------------

# Tests whether the shape of the time trend changes after the cutoff. The
# trend-break model interacts the spline basis with the post-2014 indicator,
# allowing a different trend shape before and after. A likelihood ratio test
# against the single-trend model assesses whether the additional flexibility
# is justified.

m_trend_only <- glm.nb(attacks ~ ns(time_index, df = 2),
                       data = dfm)

m_trend_break <- glm.nb(attacks ~ ns(time_index, df = 2) + post_article5_exact:ns(time_index, df = 2),
                        data = dfm)

aic_trend_break <- AIC(m_trend_only, m_trend_break) %>%
  rownames_to_column(var = "model") %>%
  mutate(model = recode(model,
                        m_trend_only = "Single trend, spline df 2",
                        m_trend_break = "Trend break after cutoff, interacted spline df 2"),
         delta_aic = AIC - min(AIC)) %>%
  arrange(AIC)

# ---- Combined monthly AIC table ---------------------------------------------

# Combines the AIC comparisons above into a single table for the appendix.

monthly_aic_combined <- bind_rows(
  aic_pois_vs_nb %>%
    transmute(test = "Poisson versus negative binomial, monthly, spline df 2",
              model = model,
              df = .data$df,
              AIC = AIC,
              delta_aic = delta_aic),
  aic_nb_vs_zinb_df2 %>%
    transmute(test = "Negative binomial versus zero inflated negative binomial, monthly, spline df 2",
              model = model,
              df = .data$df,
              AIC = AIC,
              delta_aic = delta_aic),
  aic_nb_vs_zinb_df3 %>%
    transmute(test = "Negative binomial versus zero inflated negative binomial, monthly, spline df 3",
              model = model,
              df = .data$df,
              AIC = AIC,
              delta_aic = delta_aic),
  aic_spline_df %>%
    transmute(test = "Spline flexibility, monthly negative binomial",
              model = paste("Spline df", spline_df),
              df = df_model,
              AIC = aic,
              delta_aic = delta_aic),
  aic_linear_vs_df2 %>%
    transmute(test = "Linear versus spline trend, monthly negative binomial",
              model = model,
              df = .data$df,
              AIC = AIC,
              delta_aic = delta_aic),
  aic_trend_break %>%
    transmute(test = "Single trend versus post cutoff trend break, monthly",
              model = model,
              df = .data$df,
              AIC = AIC,
              delta_aic = delta_aic))

kable(monthly_aic_combined,
      format = "latex",
      booktabs = TRUE,
      digits = 3,
      col.names = c("Test", "Model", "df", "AIC", "Delta AIC"),
      caption = "Combined AIC comparisons across monthly model checks.",
      label = "monthly_aic_combined",
      row.names = FALSE,
      escape = FALSE) %>%
  kable_styling(latex_options = "HOLD_position")

lr_stat <- 2 * (logLik(m_trend_break) - logLik(m_trend_only))
df_diff <- attr(logLik(m_trend_break), "df") - attr(logLik(m_trend_only), "df")
p_val <- pchisq(as.numeric(lr_stat), df = df_diff, lower.tail = FALSE)

trend_break_lr <- data.frame(lr_statistic = as.numeric(lr_stat),
                             df = df_diff,
                             p_value = p_val)

kable(trend_break_lr,
      format = "latex",
      booktabs = TRUE,
      digits = 3,
      col.names = c("LR statistic", "df", "p value"),
      caption = "Likelihood ratio test for a post cutoff break in spline shape in the monthly model.",
      label = "lr_trend_break_month",
      row.names = FALSE,
      escape = FALSE) %>%
  kable_styling(latex_options = "HOLD_position")

# ---- Temporal granularity comparison ----------------------------------------

# Refits the main model at four temporal resolutions (month, quarter, half-
# year, year) to check whether the substantive conclusion is sensitive to
# how the data are aggregated.

dfq <- df_quarter %>%
  mutate(period_start = quarter_start,
         period_end = quarter_start %m+% months(3),
         post_article5_exact = share_overlap(period_start,
                                             period_end,
                                             cutoff_date,
                                             as.Date("2100-01-01")))

dfh <- df_halfyear %>%
  mutate(period_start = half_start,
         period_end = half_start %m+% months(6),
         post_article5_exact = share_overlap(period_start,
                                             period_end,
                                             cutoff_date,
                                             as.Date("2100-01-01")))

dfy <- df_year %>%
  mutate(period_start = year_start,
         period_end = year_start %m+% years(1),
         post_article5_exact = share_overlap(period_start,
                                             period_end,
                                             cutoff_date,
                                             as.Date("2100-01-01")))

m_nb_month <- glm.nb(attacks ~ post_article5_exact + ns(time_index, df = 2),
                     data = dfm)

m_nb_quarter <- glm.nb(attacks ~ post_article5_exact + ns(time_index, df = 2),
                       data = dfq)

m_nb_halfyear <- glm.nb(attacks ~ post_article5_exact + ns(time_index, df = 2),
                        data = dfh)

m_nb_year <- glm.nb(attacks ~ post_article5_exact + ns(time_index, df = 2),
                    data = dfy)

m_pois_month <- glm(attacks ~ post_article5_exact + ns(time_index, df = 2),
                    family = poisson(link = "log"),
                    data = dfm)

m_pois_quarter <- glm(attacks ~ post_article5_exact + ns(time_index, df = 2),
                      family = poisson(link = "log"),
                      data = dfq)

m_pois_halfyear <- glm(attacks ~ post_article5_exact + ns(time_index, df = 2),
                       family = poisson(link = "log"),
                       data = dfh)

m_pois_year <- glm(attacks ~ post_article5_exact + ns(time_index, df = 2),
                   family = poisson(link = "log"),
                   data = dfy)

modelsummary(list("Month" = m_nb_month,
                  "Quarter" = m_nb_quarter,
                  "Half year" = m_nb_halfyear,
                  "Year" = m_nb_year),
             coef_map = coef_map_main,
             gof_map = gof_map_main,
             statistic = "({std.error})",
             stars = TRUE,
             output = "kableExtra",
             title = "Negative binomial models across temporal aggregation levels.",
             escape = FALSE) %>%
  kable_styling(latex_options = "HOLD_position")

gran_aic <- data.frame(granularity = c("Month", "Quarter", "Half year", "Year"),
                       n_periods = c(nrow(dfm), nrow(dfq), nrow(dfh), nrow(dfy)),
                       aic_poisson = c(AIC(m_pois_month), AIC(m_pois_quarter), AIC(m_pois_halfyear), AIC(m_pois_year)),
                       aic_nb = c(AIC(m_nb_month), AIC(m_nb_quarter), AIC(m_nb_halfyear), AIC(m_nb_year))) %>%
  mutate(delta_aic_poisson = aic_poisson - min(aic_poisson),
         delta_aic_nb = aic_nb - min(aic_nb)) %>%
  arrange(aic_nb)

kable(gran_aic,
      format = "latex",
      booktabs = TRUE,
      digits = 3,
      col.names = c("Granularity",
                    "Periods",
                    "AIC Poisson",
                    "AIC negative binomial",
                    "Delta AIC Poisson",
                    "Delta AIC negative binomial"),
      caption = "Model fit comparison across temporal aggregation levels with spline df 2.",
      label = "aic_granularity_compare",
      row.names = FALSE,
      escape = FALSE) %>%
  kable_styling(latex_options = "HOLD_position")

# ---- Prior sensitivity setup, quarterly model -------------------------------

# Three prior sets are used to test sensitivity. The sceptical set pulls the
# post-2014 coefficient tightly towards zero (sd_b = 0.25); the weak set is
# moderately diffuse (sd_b = 0.50); the wide set is least restrictive
# (sd_b = 1.00). Comparing posteriors across these sets shows whether the
# main results are sensitive to prior choice or driven by the data.

dfq <- df_quarter %>%
  mutate(period_start = quarter_start,
         period_end = quarter_start %m+% months(3),
         post_article5_exact = share_overlap(period_start,
                                             period_end,
                                             cutoff_date,
                                             as.Date("2100-01-01")))

basis_q <- ns(dfq$time_index, df = 2)

dfq <- dfq %>%
  mutate(ns1 = as.numeric(basis_q[, 1]),
         ns2 = as.numeric(basis_q[, 2])) %>%
  dplyr::select(attacks, time_index, post_article5_exact, ns1, ns2)

prior_grid <- data.frame(prior_set = c("sceptical", "weak", "wide"),
                         sd_b = c(0.25, 0.50, 1.00),
                         sd_g = c(0.50, 1.00, 2.00),
                         sd_a = c(0.75, 1.00, 1.50),
                         rate_phi = c(1.00, 1.00, 1.00))

prior_labels <- c(
  sceptical = "Sceptical",
  weak = "Weak",
  wide = "Wide")

# 4000 draws matches the default total number of post-warmup samples in
# rethinking (4 chains times 1000 post-warmup draws), and is sufficient to
# characterise the prior and posterior distributions with low Monte Carlo
# error.
n_draws <- 4000

X <- model.matrix(~ post_article5_exact + ns1 + ns2, data = dfq)
colnames(X) <- c("Intercept", "Post Article 5", "Time spline 1", "Time spline 2")

# The + 0.5 adjustment is a standard small-count fix that prevents log(0)
# when some periods have zero counts, while having a negligible effect on
# the location of the prior.
mu_a <- log(mean(dfq$attacks + 0.5))

# ---- Prior draws ------------------------------------------------------------

# Generates n_draws samples from each prior set. Seeds ensure reproducibility.
# Different seeds are used for prior draws, prior predictive simulations,
# model fits, and posterior draws so that each step is independently
# reproducible.

set.seed(2001)
sceptical_draws <- tibble(
  prior_set = "sceptical",
  a = rnorm(n_draws, mu_a, prior_grid$sd_a[1]),
  b_post = rnorm(n_draws, 0, prior_grid$sd_b[1]),
  g1 = rnorm(n_draws, 0, prior_grid$sd_g[1]),
  g2 = rnorm(n_draws, 0, prior_grid$sd_g[1]),
  phi = rexp(n_draws, rate = prior_grid$rate_phi[1]))

set.seed(2002)
weak_draws <- tibble(
  prior_set = "weak",
  a = rnorm(n_draws, mu_a, prior_grid$sd_a[2]),
  b_post = rnorm(n_draws, 0, prior_grid$sd_b[2]),
  g1 = rnorm(n_draws, 0, prior_grid$sd_g[2]),
  g2 = rnorm(n_draws, 0, prior_grid$sd_g[2]),
  phi = rexp(n_draws, rate = prior_grid$rate_phi[2]))

set.seed(2003)
wide_draws <- tibble(
  prior_set = "wide",
  a = rnorm(n_draws, mu_a, prior_grid$sd_a[3]),
  b_post = rnorm(n_draws, 0, prior_grid$sd_b[3]),
  g1 = rnorm(n_draws, 0, prior_grid$sd_g[3]),
  g2 = rnorm(n_draws, 0, prior_grid$sd_g[3]),
  phi = rexp(n_draws, rate = prior_grid$rate_phi[3]))

prior_draws <- bind_rows(
  sceptical_draws,
  weak_draws,
  wide_draws) %>%
  mutate(prior_label = recode(prior_set, !!!prior_labels))

prior_draws_long <- bind_rows(
  data.frame(prior_label = prior_draws$prior_label, param = "Intercept", value = prior_draws$a),
  data.frame(prior_label = prior_draws$prior_label, param = "Post Article 5", value = prior_draws$b_post),
  data.frame(prior_label = prior_draws$prior_label, param = "Time spline 1", value = prior_draws$g1),
  data.frame(prior_label = prior_draws$prior_label, param = "Time spline 2", value = prior_draws$g2))

# ---- Prior parameter densities figure ---------------------------------------

p_prior_params <- ggplot(prior_draws_long, aes(x = value)) +
  geom_density(color = color_sim, linewidth = 0.8) +
  facet_grid(param ~ prior_label, scales = "free") +
  labs(title = "Prior draws for model parameters",
       subtitle = "Quarterly specification with three alternative prior sets",
       x = NULL,
       y = "Density") +
  base_theme

p_prior_params
ggsave("Figures/prior_densities_parameters_quarter.pdf", p_prior_params, width = 8, height = 6)

# ---- Implied rate ratio figure ----------------------------------------------

# Shows the implied rate ratio (exp of the post coefficient) under each prior.
# Helps interpret the priors on the natural-scale magnitudes they imply.

prior_rr <- data.frame(prior_label = prior_draws$prior_label,
                       rr_post = exp(prior_draws$b_post))

p_prior_rr <- ggplot(prior_rr, aes(x = rr_post)) +
  geom_density(color = color_sim, linewidth = 0.8) +
  facet_wrap(~ prior_label, scales = "free_y") +
  labs(title = "Implied rate ratio under the prior on the post effect",
       subtitle = "Distribution of exp(b_post) by prior set",
       x = "Rate ratio implied by the prior on the post effect",
       y = "Density") +
  base_theme

p_prior_rr
ggsave("Figures/prior_densities_rate_ratio_quarter.pdf", p_prior_rr, width = 8, height = 4)

# ---- Prior predictive simulation --------------------------------------------

# Simulates counts from each prior set and compares to the observed quarterly
# distribution. Tells whether the prior puts mass on plausible count values
# or generates wildly implausible data.

sim_counts_sceptical <- data.frame(prior_set = character(0),
                                   prior_label = character(0),
                                   sim_y = integer(0))

set.seed(2011)
for (j in seq_len(nrow(sceptical_draws))) {
  beta <- c(sceptical_draws$a[j],
            sceptical_draws$b_post[j],
            sceptical_draws$g1[j],
            sceptical_draws$g2[j])
  mu <- as.numeric(exp(X %*% beta))
  y <- rnbinom(n = length(mu), mu = mu, size = sceptical_draws$phi[j])
  sim_counts_sceptical <- bind_rows(
    sim_counts_sceptical,
    data.frame(prior_set = "sceptical",
               prior_label = "Sceptical",
               sim_y = y))
}

sim_counts_weak <- data.frame(prior_set = character(0),
                              prior_label = character(0),
                              sim_y = integer(0))

set.seed(2012)
for (j in seq_len(nrow(weak_draws))) {
  beta <- c(weak_draws$a[j],
            weak_draws$b_post[j],
            weak_draws$g1[j],
            weak_draws$g2[j])
  mu <- as.numeric(exp(X %*% beta))
  y <- rnbinom(n = length(mu), mu = mu, size = weak_draws$phi[j])
  sim_counts_weak <- bind_rows(
    sim_counts_weak,
    data.frame(prior_set = "weak",
               prior_label = "Weak",
               sim_y = y))
}

sim_counts_wide <- data.frame(prior_set = character(0),
                              prior_label = character(0),
                              sim_y = integer(0))

set.seed(2013)
for (j in seq_len(nrow(wide_draws))) {
  beta <- c(wide_draws$a[j],
            wide_draws$b_post[j],
            wide_draws$g1[j],
            wide_draws$g2[j])
  mu <- as.numeric(exp(X %*% beta))
  y <- rnbinom(n = length(mu), mu = mu, size = wide_draws$phi[j])
  sim_counts_wide <- bind_rows(
    sim_counts_wide,
    data.frame(prior_set = "wide",
               prior_label = "Wide",
               sim_y = y))
}

sim_counts_df <- bind_rows(
  sim_counts_sceptical,
  sim_counts_weak,
  sim_counts_wide)

dist_df <- sim_counts_df %>%
  count(prior_label, sim_y, name = "n") %>%
  group_by(prior_label) %>%
  mutate(share = n / sum(n)) %>%
  ungroup()

obs_dist_df <- dfq %>%
  count(attacks, name = "n") %>%
  mutate(share = n / sum(n),
         sim_y = attacks) %>%
  dplyr::select(sim_y, share)

obs_dist_df <- data.frame(prior_label = rep(unique(dist_df$prior_label), each = nrow(obs_dist_df)),
                          sim_y = rep(obs_dist_df$sim_y, times = length(unique(dist_df$prior_label))),
                          share = rep(obs_dist_df$share, times = length(unique(dist_df$prior_label))))

p_prior_pred <- ggplot(dist_df, aes(x = sim_y, y = share)) +
  geom_col(aes(fill = "Prior predictive"),
           color = color_sim_dark,
           linewidth = 0.3) +
  geom_line(data = obs_dist_df,
            aes(x = sim_y, y = share, color = "Observed data", group = 1),
            linewidth = 0.8) +
  geom_point(data = obs_dist_df,
             aes(x = sim_y, y = share, color = "Observed data"),
             size = 1.2) +
  facet_wrap(~ prior_label) +
  scale_fill_manual(values = c("Prior predictive" = color_sim)) +
  scale_color_manual(values = c("Observed data" = color_obs)) +
  labs(title = "Prior predictive distribution for quarterly counts",
       subtitle = "Comparison of simulated counts and observed data",
       x = "Quarterly cyber incident counts",
       y = "Share",
       fill = NULL,
       color = NULL) +
  coord_cartesian(xlim = c(0, 50)) +
  base_theme

p_prior_pred
ggsave("Figures/prior_predictive_counts_quarter.pdf", p_prior_pred, width = 8, height = 4)

prior_pred_summary <- sim_counts_df %>%
  group_by(prior_label) %>%
  summarise(p_gt_20 = mean(sim_y > 20),
            p_gt_30 = mean(sim_y > 30),
            p_eq_0 = mean(sim_y == 0),
            .groups = "drop")

kable(prior_pred_summary,
      format = "latex",
      booktabs = TRUE,
      digits = 3,
      col.names = c("Prior set", "Pr(y > 20)", "Pr(y > 30)", "Pr(y = 0)"),
      caption = "Prior predictive checks for quarterly counts under alternative priors.",
      label = "prior_pred_quarter",
      row.names = FALSE,
      escape = FALSE) %>%
  kable_styling(latex_options = "HOLD_position")

# ---- Posterior sensitivity to prior choice ----------------------------------

# Fits the quarterly negative binomial model with each of the three prior
# sets and reports the posterior of the post-2014 coefficient. If the three
# posteriors are similar, the result is data-driven; if they diverge, the
# data are weak relative to the prior and the inference is prior-driven.

set.seed(2031)
fit_sceptical <- quap(alist(
  attacks ~ dgampois(lambda, phi),
  log(lambda) <- a + b_post * post_article5_exact + g1 * ns1 + g2 * ns2,
  a ~ dnorm(mu_a, prior_grid$sd_a[1]),
  b_post ~ dnorm(0, prior_grid$sd_b[1]),
  g1 ~ dnorm(0, prior_grid$sd_g[1]),
  g2 ~ dnorm(0, prior_grid$sd_g[1]),
  phi ~ dexp(prior_grid$rate_phi[1])),
  data = list(attacks = dfq$attacks,
              post_article5_exact = dfq$post_article5_exact,
              ns1 = dfq$ns1,
              ns2 = dfq$ns2,
              mu_a = mu_a))

set.seed(2032)
fit_weak <- quap(alist(
  attacks ~ dgampois(lambda, phi),
  log(lambda) <- a + b_post * post_article5_exact + g1 * ns1 + g2 * ns2,
  a ~ dnorm(mu_a, prior_grid$sd_a[2]),
  b_post ~ dnorm(0, prior_grid$sd_b[2]),
  g1 ~ dnorm(0, prior_grid$sd_g[2]),
  g2 ~ dnorm(0, prior_grid$sd_g[2]),
  phi ~ dexp(prior_grid$rate_phi[2])),
  data = list(attacks = dfq$attacks,
              post_article5_exact = dfq$post_article5_exact,
              ns1 = dfq$ns1,
              ns2 = dfq$ns2,
              mu_a = mu_a))

set.seed(2033)
fit_wide <- quap(alist(
  attacks ~ dgampois(lambda, phi),
  log(lambda) <- a + b_post * post_article5_exact + g1 * ns1 + g2 * ns2,
  a ~ dnorm(mu_a, prior_grid$sd_a[3]),
  b_post ~ dnorm(0, prior_grid$sd_b[3]),
  g1 ~ dnorm(0, prior_grid$sd_g[3]),
  g2 ~ dnorm(0, prior_grid$sd_g[3]),
  phi ~ dexp(prior_grid$rate_phi[3])),
  data = list(attacks = dfq$attacks,
              post_article5_exact = dfq$post_article5_exact,
              ns1 = dfq$ns1,
              ns2 = dfq$ns2,
              mu_a = mu_a))

set.seed(2041)
draws_sceptical_post <- extract.samples(fit_sceptical, n = 4000)$b_post

set.seed(2042)
draws_weak_post <- extract.samples(fit_weak, n = 4000)$b_post

set.seed(2043)
draws_wide_post <- extract.samples(fit_wide, n = 4000)$b_post

posterior_draws <- tibble(
  prior_label = c("Sceptical", "Weak", "Wide"),
  post_median = c(median(draws_sceptical_post),
                  median(draws_weak_post),
                  median(draws_wide_post)),
  post_l95 = c(quantile(draws_sceptical_post, 0.025),
               quantile(draws_weak_post, 0.025),
               quantile(draws_wide_post, 0.025)),
  post_u95 = c(quantile(draws_sceptical_post, 0.975),
               quantile(draws_weak_post, 0.975),
               quantile(draws_wide_post, 0.975)),
  p_gt0 = c(mean(draws_sceptical_post > 0),
            mean(draws_weak_post > 0),
            mean(draws_wide_post > 0)),
  rr_median = c(exp(median(draws_sceptical_post)),
                exp(median(draws_weak_post)),
                exp(median(draws_wide_post))),
  rr_l95 = c(exp(quantile(draws_sceptical_post, 0.025)),
             exp(quantile(draws_weak_post, 0.025)),
             exp(quantile(draws_wide_post, 0.025))),
  rr_u95 = c(exp(quantile(draws_sceptical_post, 0.975)),
             exp(quantile(draws_weak_post, 0.975)),
             exp(quantile(draws_wide_post, 0.975))))

kable(posterior_draws,
      format = "latex",
      booktabs = TRUE,
      digits = 3,
      col.names = c("Prior",
                    "$\\beta_{post}$",
                    "l95",
                    "u95",
                    "$Pr(\\beta > 0)$",
                    "Rate ratio",
                    "Rate ratio l95",
                    "Rate ratio u95"),
      caption = "Posterior sensitivity of the post coefficient under alternative priors.",
      label = "posterior_prior_sensitivity",
      row.names = FALSE,
      escape = FALSE) %>%
  kable_styling(latex_options = "HOLD_position")