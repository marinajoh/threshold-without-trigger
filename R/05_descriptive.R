# =============================================================================
# descriptive.R
#
# Produces descriptive figures and tables for the compositional analysis.
# Generates pre/post bar charts for method and severity, half-yearly
# composition time series for target type, coercive objective, damage type,
# and critical infrastructure sector, and comparative pre/post tables
# (NATO target versus non-NATO target) for four compositional dimensions.
# Output figures are saved to Figures/; tables are printed to the console
# in LaTeX form.
# =============================================================================

library(dplyr)
library(ggplot2)
library(lubridate)
library(tidyr)
library(scales)
library(knitr)

options(mc.cores = 1)

# ---- Read derived series ----------------------------------------------------

# Reads aggregated event-level data produced by clean_data.R. The descriptive
# analysis uses the event-level data, not the aggregated count series, so the
# focus is on df_filtered (NATO-target events) and df_events (full event set
# used for the comparative tables).

obj <- readRDS("data/derived_series.rds")

df <- obj$df_filtered %>%
  mutate(
    interactionstartdate = as.Date(interactionstartdate)) %>%
  filter(
    year(interactionstartdate) <= 2020)

cutoff_date <- as.Date(obj$article5_cyber_date)

dir.create("Figures", showWarnings = FALSE, recursive = TRUE)

# ---- Palette and base theme -------------------------------------------------

# Palette and base theme are duplicated across descriptive.R, modeltest.R,
# and dag.R so that each script can be run independently of the others.
# Any change here should be mirrored in those files.

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

palette_full <- c(
  "#6C8FD4",
  "#A78FD6",
  "#7FC8B2",
  "#B8D986",
  "#E6A86C",
  "#D89595",
  "#C7B486")

# ---- Colour mappings for compositional categories ---------------------------

# Each compositional dimension gets a named palette so figures use consistent
# colours across pre/post and time-series views.

method_colors <- setNames(
  palette_full[1:4],
  c("Vandalism", "Denial of service", "Network intrusion", "Network infiltration"))

target_colors <- setNames(
  palette_full[1:3],
  c("Private or non-state", "Government non-military", "Government military"))

objective_colors <- setNames(
  palette_full[1:4],
  c("Disruption", "Short-term espionage", "Long-term espionage", "Degradation"))

severity_colors <- setNames(
  palette_full[1:3],
  c("Low severity (1 to 2)", "Moderate severity (3 to 4)", "High severity (5)"))

damage_colors <- setNames(
  palette_full[1:4],
  c("Direct and immediate", "Direct and delayed", "Indirect and immediate", "Indirect and delayed"))

crit_inf_colors <- setNames(
  palette_full[1:7],
  c(
    "Government facilities",
    "Defence sector",
    "Energy",
    "Financial services",
    "Information infrastructure",
    "Election infrastructure and academia",
    "Other essential services"))

# Two-line labels for critical infrastructure used in the figure legend, so
# that long category names do not crowd the layout.
crit_inf_labels <- c(
  "Government facilities" = "Government\nfacilities",
  "Defence sector" = "Defence\nsector",
  "Energy" = "Energy",
  "Financial services" = "Financial\nservices",
  "Information infrastructure" = "Information\ninfrastructure",
  "Election infrastructure and academia" = "Election infrastructure\nand academia",
  "Other essential services" = "Other essential\nservices")

# ---- Compositional mappings from DCID codes to thesis categories ------------

# Each compositional dimension is mapped from DCID's numeric coding to the
# named categories used in the thesis. The grouping decisions follow
# Section 5 (Compositional Outcome Dimensions) in the thesis and the DCID
# codebook. Records with no valid code are dropped.

# Method of interaction. DCID codes 3 and 3.1 are grouped as network
# intrusion; 4 and its sub-codes (4.1 to 4.4) are grouped as network
# infiltration.
df_method_base <- df %>%
  mutate(
    method_main = case_when(
      method == 1 ~ "Vandalism",
      method == 2 ~ "Denial of service",
      method %in% c(3, 3.1) ~ "Network intrusion",
      method %in% c(4, 4.1, 4.2, 4.3, 4.4) ~ "Network infiltration",
      TRUE ~ NA_character_)) %>%
  filter(!is.na(method_main)) %>%
  mutate(
    method_main = factor(
      method_main,
      levels = c(
        "Vandalism",
        "Denial of service",
        "Network intrusion",
        "Network infiltration")))

# Type of target. Maps DCID targettype codes 1 to 3.
df_target_base <- df %>%
  mutate(
    target_main = case_when(
      targettype == 1 ~ "Private or non-state",
      targettype == 2 ~ "Government non-military",
      targettype == 3 ~ "Government military",
      TRUE ~ NA_character_)) %>%
  filter(!is.na(target_main)) %>%
  mutate(
    target_main = factor(
      target_main,
      levels = c(
        "Private or non-state",
        "Government non-military",
        "Government military")))

# Coercive objective. Maps DCID cyber_objective codes 1 to 4.
df_objective_base <- df %>%
  mutate(
    objective_main = case_when(
      cyber_objective == 1 ~ "Disruption",
      cyber_objective == 2 ~ "Short-term espionage",
      cyber_objective == 3 ~ "Long-term espionage",
      cyber_objective == 4 ~ "Degradation",
      TRUE ~ NA_character_)) %>%
  filter(!is.na(objective_main)) %>%
  mutate(
    objective_main = factor(
      objective_main,
      levels = c(
        "Disruption",
        "Short-term espionage",
        "Long-term espionage",
        "Degradation")))

# Severity. DCID severity scale 1 to 5, grouped into three bands.
df_severity_base <- df %>%
  mutate(
    severity_main = case_when(
      severity %in% 1:2 ~ "Low severity (1 to 2)",
      severity %in% 3:4 ~ "Moderate severity (3 to 4)",
      severity == 5 ~ "High severity (5)",
      TRUE ~ NA_character_)) %>%
  filter(!is.na(severity_main)) %>%
  mutate(
    severity_main = factor(
      severity_main,
      levels = c(
        "Low severity (1 to 2)",
        "Moderate severity (3 to 4)",
        "High severity (5)")))

# Type of damage. Maps DCID damage_type codes 1 to 4.
df_damage_base <- df %>%
  mutate(
    damage_main = case_when(
      damage_type == 1 ~ "Direct and immediate",
      damage_type == 2 ~ "Direct and delayed",
      damage_type == 3 ~ "Indirect and immediate",
      damage_type == 4 ~ "Indirect and delayed",
      TRUE ~ NA_character_)) %>%
  filter(!is.na(damage_main)) %>%
  mutate(
    damage_main = factor(
      damage_main,
      levels = c(
        "Direct and immediate",
        "Direct and delayed",
        "Indirect and immediate",
        "Indirect and delayed")))

# Critical infrastructure sector. DCID codes 17 underlying sectors based on
# the U.S. Department of Homeland Security framework. These are grouped into
# seven categories for analysis: six individual or paired sectors plus a
# residual category for the remaining sectors. The grouping is described in
# Section 5 (Compositional Outcome Dimensions) of the thesis.
df_critinf_base <- df %>%
  mutate(
    crit_inf_main = case_when(
      crit_inf == 11 ~ "Government facilities",
      crit_inf == 6 ~ "Defence sector",
      crit_inf == 8 ~ "Energy",
      crit_inf == 9 ~ "Financial services",
      crit_inf %in% c(3, 13) ~ "Information infrastructure",
      crit_inf == 17 ~ "Election infrastructure and academia",
      crit_inf %in% c(1, 2, 4, 5, 7, 10, 12, 14, 15, 16) ~ "Other essential services",
      TRUE ~ NA_character_)) %>%
  filter(!is.na(crit_inf_main)) %>%
  mutate(
    crit_inf_main = factor(
      crit_inf_main,
      levels = c(
        "Government facilities",
        "Defence sector",
        "Energy",
        "Financial services",
        "Information infrastructure",
        "Election infrastructure and academia",
        "Other essential services")))

# ---- Pre/post composition for method and severity ---------------------------

# Counts incidents per period (before/after the 5 September 2014 cutoff)
# and within each compositional category, then computes the share of
# incidents in each category by period.

df_method_prepost <- df_method_base %>%
  mutate(
    period = case_when(
      interactionstartdate < cutoff_date ~ "Before Article 5 recognition",
      interactionstartdate >= cutoff_date ~ "After Article 5 recognition"),
    period = factor(
      period,
      levels = c("Before Article 5 recognition", "After Article 5 recognition"))) %>%
  count(period, method_main, name = "n") %>%
  group_by(period) %>%
  mutate(share = n / sum(n)) %>%
  ungroup()

df_severity_prepost <- df_severity_base %>%
  mutate(
    period = case_when(
      interactionstartdate < cutoff_date ~ "Before Article 5 recognition",
      interactionstartdate >= cutoff_date ~ "After Article 5 recognition"),
    period = factor(
      period,
      levels = c("Before Article 5 recognition", "After Article 5 recognition"))) %>%
  count(period, severity_main, name = "n") %>%
  group_by(period) %>%
  mutate(share = n / sum(n)) %>%
  ungroup()

p_method_prepost <- ggplot(
  df_method_prepost,
  aes(x = period, y = share, fill = method_main)) +
  geom_col(position = "fill", width = 0.7, color = "white", linewidth = 0.3) +
  scale_fill_manual(values = method_colors) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title = "Method composition",
    x = NULL,
    y = "Share of incidents",
    fill = "Method") +
  base_theme

p_severity_prepost <- ggplot(
  df_severity_prepost,
  aes(x = period, y = share, fill = severity_main)) +
  geom_col(position = "fill", width = 0.7, color = "white", linewidth = 0.3) +
  scale_fill_manual(values = severity_colors) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title = "Severity composition",
    x = NULL,
    y = "Share of incidents",
    fill = "Severity") +
  base_theme

ggsave("Figures/method_composition_prepost.pdf", p_method_prepost, width = 9, height = 5)
ggsave("Figures/severity_composition_prepost.pdf", p_severity_prepost, width = 9, height = 5)

# ---- Half-yearly composition time series ------------------------------------

# For four dimensions (target, objective, damage, critical infrastructure),
# the share of each category is computed per half-year. The half-year index
# is calculated as the number of six-month periods between the cutoff date
# and the incident date; period_start aligns the periods to the cutoff so
# the visual reference at 5 September 2014 falls on a period boundary.
#
# Periods with fewer than two recorded incidents are filtered out as a
# relevance threshold for stability in sparse windows. Shares computed
# from a single event are not informative and can introduce volatile noise
# in periods where the data are thin.

df_target_share <- df_target_base %>%
  mutate(
    halfyear_index = floor(
      time_length(interval(cutoff_date, interactionstartdate), "months") / 6),
    period_start = cutoff_date %m+% months(6 * halfyear_index)) %>%
  count(period_start, target_main, name = "n") %>%
  complete(
    period_start = seq(min(period_start), max(period_start), by = "6 months"),
    target_main,
    fill = list(n = 0)) %>%
  group_by(period_start) %>%
  mutate(
    total_period = sum(n),
    share = ifelse(total_period >= 2, n / total_period, NA_real_)) %>%
  ungroup() %>%
  filter(!is.na(share))

df_objective_share <- df_objective_base %>%
  mutate(
    halfyear_index = floor(
      time_length(interval(cutoff_date, interactionstartdate), "months") / 6),
    period_start = cutoff_date %m+% months(6 * halfyear_index)) %>%
  count(period_start, objective_main, name = "n") %>%
  complete(
    period_start = seq(min(period_start), max(period_start), by = "6 months"),
    objective_main,
    fill = list(n = 0)) %>%
  group_by(period_start) %>%
  mutate(
    total_period = sum(n),
    share = ifelse(total_period >= 2, n / total_period, NA_real_)) %>%
  ungroup() %>%
  filter(!is.na(share))

df_damage_share <- df_damage_base %>%
  mutate(
    halfyear_index = floor(
      time_length(interval(cutoff_date, interactionstartdate), "months") / 6),
    period_start = cutoff_date %m+% months(6 * halfyear_index)) %>%
  count(period_start, damage_main, name = "n") %>%
  complete(
    period_start = seq(min(period_start), max(period_start), by = "6 months"),
    damage_main,
    fill = list(n = 0)) %>%
  group_by(period_start) %>%
  mutate(
    total_period = sum(n),
    share = ifelse(total_period >= 2, n / total_period, NA_real_)) %>%
  ungroup() %>%
  filter(!is.na(share))

df_critinf_share <- df_critinf_base %>%
  mutate(
    halfyear_index = floor(
      time_length(interval(cutoff_date, interactionstartdate), "months") / 6),
    period_start = cutoff_date %m+% months(6 * halfyear_index)) %>%
  count(period_start, crit_inf_main, name = "n") %>%
  complete(
    period_start = seq(min(period_start), max(period_start), by = "6 months"),
    crit_inf_main,
    fill = list(n = 0)) %>%
  group_by(period_start) %>%
  mutate(
    total_period = sum(n),
    share = ifelse(total_period >= 2, n / total_period, NA_real_)) %>%
  ungroup() %>%
  filter(!is.na(share))

# ---- Time series figures ----------------------------------------------------

p_target_share <- ggplot(
  df_target_share,
  aes(x = period_start, y = share, fill = target_main)) +
  geom_col(position = "fill", width = 160, color = "white", linewidth = 0.3) +
  geom_vline(xintercept = cutoff_date, colour = "black", linetype = "22", linewidth = 0.5) +
  scale_fill_manual(values = target_colors) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y", expand = expansion(mult = c(0.02, 0.04))) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(title = "Target type composition over time", x = NULL, y = "Share of incidents", fill = "Target type") +
  base_theme

p_objective_share <- ggplot(
  df_objective_share,
  aes(x = period_start, y = share, fill = objective_main)) +
  geom_col(position = "fill", width = 160, color = "white", linewidth = 0.3) +
  geom_vline(xintercept = cutoff_date, colour = "black", linetype = "22", linewidth = 0.5) +
  scale_fill_manual(values = objective_colors) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y", expand = expansion(mult = c(0.02, 0.04))) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(title = "Coercive objective composition over time", x = NULL, y = "Share of incidents", fill = "Objective") +
  base_theme

p_damage_share <- ggplot(
  df_damage_share,
  aes(x = period_start, y = share, fill = damage_main)) +
  geom_col(position = "fill", width = 160, color = "white", linewidth = 0.3) +
  geom_vline(xintercept = cutoff_date, colour = "black", linetype = "22", linewidth = 0.5) +
  scale_fill_manual(values = damage_colors) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y", expand = expansion(mult = c(0.02, 0.04))) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(title = "Damage type composition over time", x = NULL, y = "Share of incidents", fill = "Damage type") +
  base_theme

p_critinf_share <- ggplot(
  df_critinf_share,
  aes(x = period_start, y = share, fill = crit_inf_main)) +
  geom_col(position = "fill", width = 160, color = "white", linewidth = 0.3) +
  geom_vline(xintercept = cutoff_date, colour = "black", linetype = "22", linewidth = 0.5) +
  scale_fill_manual(
    values = crit_inf_colors,
    labels = crit_inf_labels) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y", expand = expansion(mult = c(0.02, 0.04))) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(title = "Critical infrastructure sector composition over time", x = NULL, y = "Share of incidents", fill = "Critical infrastructure sector") +
  base_theme +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

ggsave("Figures/target_composition_share.pdf", p_target_share, width = 9, height = 5)
ggsave("Figures/objective_composition_share.pdf", p_objective_share, width = 9, height = 5)
ggsave("Figures/damage_composition_share.pdf", p_damage_share, width = 9, height = 5)
ggsave("Figures/critinf_composition_share.pdf", p_critinf_share, width = 9, height = 5)

# ---- Comparative event-level data -------------------------------------------

# Builds the event-level frame used for the comparative tables: pools both
# NATO-target and non-NATO-target incidents initiated by non-NATO states,
# tags each event with its target group and pre/post period.

df_comp <- obj$df_events %>%
  mutate(
    interactionstartdate = as.Date(interactionstartdate)) %>%
  filter(
    initiator %in% obj$non_nato_cow,
    year(interactionstartdate) <= 2020) %>%
  mutate(
    target_nato = ifelse(state_a %in% obj$nato | state_b %in% obj$nato, 1L, 0L),
    group = ifelse(target_nato == 1L, "NATO target", "Non-NATO target"),
    period = ifelse(
      interactionstartdate < cutoff_date,
      "Before Article 5 recognition",
      "After Article 5 recognition"),
    period = factor(period, levels = c("Before Article 5 recognition", "After Article 5 recognition")),
    group = factor(group, levels = c("NATO target", "Non-NATO target")))

# ---- Comparative compositional mappings -------------------------------------

# The same DCID-to-category mappings used above, applied here to the
# comparative event-level frame. Repeated rather than abstracted into a
# helper function so the script remains self-contained and the mapping
# logic is visible at each use.

df_target_comp_base <- df_comp %>%
  mutate(
    target_main = case_when(
      targettype == 1 ~ "Private or non-state",
      targettype == 2 ~ "Government non-military",
      targettype == 3 ~ "Government military",
      TRUE ~ NA_character_)) %>%
  filter(!is.na(target_main)) %>%
  mutate(
    target_main = factor(
      target_main,
      levels = c("Private or non-state", "Government non-military", "Government military")))

df_objective_comp_base <- df_comp %>%
  mutate(
    objective_main = case_when(
      cyber_objective == 1 ~ "Disruption",
      cyber_objective == 2 ~ "Short-term espionage",
      cyber_objective == 3 ~ "Long-term espionage",
      cyber_objective == 4 ~ "Degradation",
      TRUE ~ NA_character_)) %>%
  filter(!is.na(objective_main)) %>%
  mutate(
    objective_main = factor(
      objective_main,
      levels = c("Disruption", "Short-term espionage", "Long-term espionage", "Degradation")))

df_damage_comp_base <- df_comp %>%
  mutate(
    damage_main = case_when(
      damage_type == 1 ~ "Direct and immediate",
      damage_type == 2 ~ "Direct and delayed",
      damage_type == 3 ~ "Indirect and immediate",
      damage_type == 4 ~ "Indirect and delayed",
      TRUE ~ NA_character_)) %>%
  filter(!is.na(damage_main)) %>%
  mutate(
    damage_main = factor(
      damage_main,
      levels = c("Direct and immediate", "Direct and delayed", "Indirect and immediate", "Indirect and delayed")))

df_critinf_comp_base <- df_comp %>%
  mutate(
    crit_inf_main = case_when(
      crit_inf == 11 ~ "Government facilities",
      crit_inf == 6 ~ "Defence sector",
      crit_inf == 8 ~ "Energy",
      crit_inf == 9 ~ "Financial services",
      crit_inf %in% c(3, 13) ~ "Information infrastructure",
      crit_inf == 17 ~ "Election infrastructure and academia",
      crit_inf %in% c(1, 2, 4, 5, 7, 10, 12, 14, 15, 16) ~ "Other essential services",
      TRUE ~ NA_character_)) %>%
  filter(!is.na(crit_inf_main)) %>%
  mutate(
    crit_inf_main = factor(
      crit_inf_main,
      levels = c(
        "Government facilities",
        "Defence sector",
        "Energy",
        "Financial services",
        "Information infrastructure",
        "Election infrastructure and academia",
        "Other essential services")))

# ---- Comparative tables -----------------------------------------------------

# Each comparative table reports incident counts and within-cell shares for
# the four combinations of group (NATO target, Non-NATO target) and period
# (before, after the recognition). Cell entries are formatted as "n (s%)".

target_comp_table <- df_target_comp_base %>%
  count(group, period, target_main, name = "n") %>%
  group_by(group, period) %>%
  mutate(share = 100 * n / sum(n)) %>%
  ungroup() %>%
  mutate(
    group_short = ifelse(group == "NATO target", "nato", "nonnato"),
    period_short = ifelse(period == "Before Article 5 recognition", "pre", "post"),
    cell = sprintf("%d (%.1f\\%%)", n, share)) %>%
  dplyr::select(group_short, period_short, target_main, cell) %>%
  pivot_wider(
    names_from = c(group_short, period_short),
    values_from = cell,
    names_glue = "{group_short}_{period_short}",
    values_fill = "0 (0.0)") %>%
  transmute(
    Category = target_main,
    `NATO pre` = nato_pre,
    `NATO post` = nato_post,
    `Non-NATO pre` = nonnato_pre,
    `Non-NATO post` = nonnato_post)

objective_comp_table <- df_objective_comp_base %>%
  count(group, period, objective_main, name = "n") %>%
  group_by(group, period) %>%
  mutate(share = 100 * n / sum(n)) %>%
  ungroup() %>%
  mutate(
    group_short = ifelse(group == "NATO target", "nato", "nonnato"),
    period_short = ifelse(period == "Before Article 5 recognition", "pre", "post"),
    cell = sprintf("%d (%.1f\\%%)", n, share)) %>%
  dplyr::select(group_short, period_short, objective_main, cell) %>%
  pivot_wider(
    names_from = c(group_short, period_short),
    values_from = cell,
    names_glue = "{group_short}_{period_short}",
    values_fill = "0 (0.0)") %>%
  transmute(
    Category = objective_main,
    `NATO pre` = nato_pre,
    `NATO post` = nato_post,
    `Non-NATO pre` = nonnato_pre,
    `Non-NATO post` = nonnato_post)

damage_comp_table <- df_damage_comp_base %>%
  count(group, period, damage_main, name = "n") %>%
  group_by(group, period) %>%
  mutate(share = 100 * n / sum(n)) %>%
  ungroup() %>%
  mutate(
    group_short = ifelse(group == "NATO target", "nato", "nonnato"),
    period_short = ifelse(period == "Before Article 5 recognition", "pre", "post"),
    cell = sprintf("%d (%.1f\\%%)", n, share)) %>%
  dplyr::select(group_short, period_short, damage_main, cell) %>%
  pivot_wider(
    names_from = c(group_short, period_short),
    values_from = cell,
    names_glue = "{group_short}_{period_short}",
    values_fill = "0 (0.0)") %>%
  transmute(
    Category = damage_main,
    `NATO pre` = nato_pre,
    `NATO post` = nato_post,
    `Non-NATO pre` = nonnato_pre,
    `Non-NATO post` = nonnato_post)

critinf_comp_table <- df_critinf_comp_base %>%
  count(group, period, crit_inf_main, name = "n") %>%
  group_by(group, period) %>%
  mutate(share = 100 * n / sum(n)) %>%
  ungroup() %>%
  mutate(
    group_short = ifelse(group == "NATO target", "nato", "nonnato"),
    period_short = ifelse(period == "Before Article 5 recognition", "pre", "post"),
    cell = sprintf("%d (%.1f\\%%)", n, share)) %>%
  dplyr::select(group_short, period_short, crit_inf_main, cell) %>%
  pivot_wider(
    names_from = c(group_short, period_short),
    values_from = cell,
    names_glue = "{group_short}_{period_short}",
    values_fill = "0 (0.0)") %>%
  transmute(
    Category = crit_inf_main,
    `NATO pre` = nato_pre,
    `NATO post` = nato_post,
    `Non-NATO pre` = nonnato_pre,
    `Non-NATO post` = nonnato_post)

# ---- Print comparative tables to LaTeX --------------------------------------

kable(
  target_comp_table,
  format = "latex",
  booktabs = TRUE,
  align = c("l", "r", "r", "r", "r"),
  caption = "Comparative target type composition before and after the Article 5 recognition.",
  label = "target_composition_comparative_prepost",
  row.names = FALSE,
  escape = FALSE)

kable(
  objective_comp_table,
  format = "latex",
  booktabs = TRUE,
  align = c("l", "r", "r", "r", "r"),
  caption = "Comparative coercive objective composition before and after the Article 5 recognition.",
  label = "objective_composition_comparative_prepost",
  row.names = FALSE,
  escape = FALSE)

kable(
  damage_comp_table,
  format = "latex",
  booktabs = TRUE,
  align = c("l", "r", "r", "r", "r"),
  caption = "Comparative damage type composition before and after the Article 5 recognition.",
  label = "damage_composition_comparative_prepost",
  row.names = FALSE,
  escape = FALSE)

kable(
  critinf_comp_table,
  format = "latex",
  booktabs = TRUE,
  align = c("l", "r", "r", "r", "r"),
  caption = "Comparative critical infrastructure sector composition before and after the Article 5 recognition.",
  label = "critinf_composition_comparative_prepost",
  row.names = FALSE,
  escape = FALSE)