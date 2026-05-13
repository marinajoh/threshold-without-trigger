# =============================================================================
# clean_data.R
#
# Reads the Dyadic Cyber Incident and Campaign Dataset (DCID v2.0), filters
# it to the thesis population, and produces aggregated incident-count series
# at four temporal resolutions (month, quarter, half-year, year). Output is
# a single RDS file consumed by descriptive.R, regression.R, and modeltest.R.
# =============================================================================

library(readxl)
library(dplyr)
library(janitor)
library(lubridate)
library(tidyr)
library(tibble)

# ---- Paths ------------------------------------------------------------------

data_path <- "data/DCID_2.0_Release_update_February_2023.xlsx"
out_path  <- "data/derived_series.rds"

# ---- Constants --------------------------------------------------------------

# Wales Summit Declaration, paragraph 72.
article5_cyber_date <- as.Date("2014-09-05")

# NATO members appearing in DCID v2.0 within the observation window.
# Used to identify the NATO-target series.
nato <- c(
  "US", "Canada", "UK", "France", "Germany", "Poland", "Turkey",
  "Estonia", "Lithuania")

# Non-NATO initiator states (Correlates of War country codes).
# Includes the four dominant cyber initiators discussed in the thesis
# (Russia, China, Iran, North Korea) and other non-NATO initiators recorded
# in DCID. Restricting initiators to non-NATO states means the analysis
# captures adversary cyber behaviour towards NATO, not NATO-initiated
# activity.
non_nato_cow <- c(
  365, 369, 372, 630, 652, 666, 710, 713, 731, 732, 740, 750, 770, 816)

# ---- Read and clean ---------------------------------------------------------

# Read raw DCID and standardise column names to snake_case.
df_raw <- read_excel(data_path) %>%
  clean_names(case = "snake")

# Cast start date to Date and restrict to the 2000 to 2020 observation window.
# DCID v2.0 includes some incidents with 2021 start dates; these fall outside
# the thesis window and are excluded.
df_events <- df_raw %>%
  mutate(
    interactionstartdate = as.Date(interactionstartdate)) %>%
  filter(
    interactionstartdate < as.Date("2021-01-01"))

# ---- Main NATO-target series ------------------------------------------------

# Filter to incidents initiated by non-NATO states and targeting at least one
# NATO member. This is the primary analytical population (N = 186 in the
# thesis). The post-2014 exposure indicator is constructed downstream in
# modeltest.R and regression.R using share_overlap() on the aggregated
# series, so no binary post-period column is added here.
df_filtered <- df_events %>%
  filter(
    initiator %in% non_nato_cow,
    state_a %in% nato | state_b %in% nato)

stopifnot(!anyNA(df_filtered$interactionstartdate))

# ---- Comparative event-level data -------------------------------------------

# Pools the NATO-target series with a comparison series of incidents between
# non-NATO states. Aggregated to quarters downstream in regression.R for
# the comparative regression models, and used pre/post in descriptive.R for
# the comparative composition tables.
df_comp_events <- df_events %>%
  filter(
    initiator %in% non_nato_cow) %>%
  mutate(
    target_nato = ifelse(state_a %in% nato | state_b %in% nato, 1L, 0L),
    group = ifelse(target_nato == 1L, "NATO target", "Non-NATO target"))

stopifnot(!anyNA(df_comp_events$interactionstartdate))

# ---- Aggregated count series ------------------------------------------------

# Each aggregated series counts the number of recorded incidents per period
# in the main NATO-target series. A time_index column gives a 1-based integer
# index for time-trend specifications. Periods with zero recorded incidents
# are filled with explicit zeros so the resulting series is dense in time.

# Monthly.
df_month <- df_filtered %>%
  mutate(month = floor_date(interactionstartdate, "month")) %>%
  count(month, name = "attacks") %>%
  arrange(month)

all_months <- tibble(
  month = seq(min(df_month$month, na.rm = TRUE),
              max(df_month$month, na.rm = TRUE),
              by = "1 month"))

df_month <- all_months %>%
  left_join(df_month, by = "month") %>%
  mutate(
    attacks = ifelse(is.na(attacks), 0L, attacks),
    time_index = row_number()) %>%
  dplyr::select(month, attacks, time_index)

stopifnot(!anyDuplicated(df_month$month))

# Quarterly.
df_quarter <- df_filtered %>%
  mutate(quarter_start = floor_date(interactionstartdate, "quarter")) %>%
  count(quarter_start, name = "attacks") %>%
  arrange(quarter_start)

all_quarters <- tibble(
  quarter_start = seq(min(df_quarter$quarter_start, na.rm = TRUE),
                      max(df_quarter$quarter_start, na.rm = TRUE),
                      by = "3 months"))

df_quarter <- all_quarters %>%
  left_join(df_quarter, by = "quarter_start") %>%
  mutate(
    attacks = ifelse(is.na(attacks), 0L, attacks),
    time_index = row_number()) %>%
  dplyr::select(quarter_start, attacks, time_index)

stopifnot(!anyDuplicated(df_quarter$quarter_start))

# Half-yearly.
df_halfyear <- df_filtered %>%
  mutate(half_start = floor_date(interactionstartdate, "halfyear")) %>%
  count(half_start, name = "attacks") %>%
  arrange(half_start)

all_halfyears <- tibble(
  half_start = seq(min(df_halfyear$half_start, na.rm = TRUE),
                   max(df_halfyear$half_start, na.rm = TRUE),
                   by = "6 months"))

df_halfyear <- all_halfyears %>%
  left_join(df_halfyear, by = "half_start") %>%
  mutate(
    attacks = ifelse(is.na(attacks), 0L, attacks),
    time_index = row_number()) %>%
  dplyr::select(half_start, attacks, time_index)

stopifnot(!anyDuplicated(df_halfyear$half_start))

# Yearly.
df_year <- df_filtered %>%
  mutate(year_start = floor_date(interactionstartdate, "year")) %>%
  count(year_start, name = "attacks") %>%
  arrange(year_start)

all_years <- tibble(
  year_start = seq(min(df_year$year_start, na.rm = TRUE),
                   max(df_year$year_start, na.rm = TRUE),
                   by = "1 year"))

df_year <- all_years %>%
  left_join(df_year, by = "year_start") %>%
  mutate(
    attacks = ifelse(is.na(attacks), 0L, attacks),
    time_index = row_number()) %>%
  dplyr::select(year_start, attacks, time_index)

stopifnot(!anyDuplicated(df_year$year_start))

# ---- Save -------------------------------------------------------------------

# Save all derived objects in one RDS file. Downstream scripts read the list
# and extract what they need.
saveRDS(
  list(
    df_raw = df_raw,
    df_events = df_events,
    df_filtered = df_filtered,
    df_comp_events = df_comp_events,
    df_month = df_month,
    df_quarter = df_quarter,
    df_halfyear = df_halfyear,
    df_year = df_year,
    nato = nato,
    non_nato_cow = non_nato_cow,
    article5_cyber_date = article5_cyber_date),
  out_path)

message("Saved: ", out_path)