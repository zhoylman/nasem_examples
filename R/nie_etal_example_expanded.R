# ============================================================
# Title: Simulated Drought Index + Drought Days Figure (100-year version; expanding climatology)
# Author: Zachary H. Hoylman
# Date: 1-29-2026
# Description:
#   Recreates the layout of a 4-panel figure (Nie et al., 2025 style):
#     - Left: normalized daily terrestrial water storage (TWS) "drought index" time series
#     - Right: stacked annual bars of drought days by category (D0–D4)
#
#   Key updates to Nie:
#     - Uses a 100 year climatology instead of a 30 year
#     - Uses an expanding reference climatology for standardization & drought classes.
#     - Minimum climatology length = 30 years.
#     - Years 1–30 are standardized using a fixed 30-year baseline (years 1–30).
#     - For year y >= 30, z-scores for year y are computed relative to YEARS 1..y
#       (so year 31 uses 31-year clim, year 60 uses 60-year clim, ..., year 100 uses 100-year clim).
#   Note:
#     - This is a z-score standardization (mean/sd). It mimics a standardized index but is not SPI/SPEI.
#     - Drought classes use USDM-like z thresholds:
#         D0: [-0.8, -0.5]
#         D1: [-1.3, -0.8]
#         D2: [-1.6, -1.3]
#         D3: [-2.0, -1.6]
#         D4: <= -2.0
#
# Reference:
#   Nie et al., 2025. One Earth. https://doi.org/10.1016/j.oneear.2025.101196
# ============================================================

# ---- Libraries ----
library(tidyverse)
library(cowplot)
library(scales)
library(shadowtext)

# ---- User settings ----
years_total = 100
days_per_year = 365
min_climatology_years = 30
seed = 42

# ---- Drought categories in PERCENTILE SPACE ----
# Percentile p in [0,1], smaller p = drier.
# D0: 20-30%, D1: 10-20%, D2: 5-10%, D3: 2-5%, D4: 0-2%
drought_bins_pct = tibble::tibble(
  cat   = factor(c("D0", "D1", "D2", "D3", "D4"), levels = c("D0","D1","D2","D3","D4")),
  lower = c(0.20, 0.10, 0.05, 0.02, -Inf),
  upper = c(0.30, 0.20, 0.10, 0.05, 0.02)
)

# Colors roughly matching USDM palette
drought_colors = c(
  "D0" = "#ffff05",
  "D1" = "#fdd37f",
  "D2" = "#ffaa01",
  "D3" = "#e60001",
  "D4" = "#730000"
)

# ---- Helpers ----
make_time_index = function(years_total, days_per_year) {
  n = years_total * days_per_year
  tibble::tibble(
    t = seq_len(n),
    year = ((t - 1) %/% days_per_year) + 1,
    day_of_year = ((t - 1) %% days_per_year) + 1,
    year_frac = year + (day_of_year - 1) / days_per_year
  )
}

# ---- Percentile-based expanding climatology ----
# Returns p in [0,1] where small p = unusually low (drier) values.
percentile_with_min_then_expanding_climatology = function(x, year, min_years = 30) {
  stopifnot(length(x) == length(year))
  p = rep(NA_real_, length(x))
  idx_by_year = split(seq_along(year), year)
  years_unique = sort(unique(year))
  
  # fixed baseline years 1..min_years
  fixed_years = seq_len(min_years)
  fixed_idx = unlist(idx_by_year[names(idx_by_year) %in% fixed_years])
  ref_fixed = x[fixed_idx]
  ecdf_fixed = stats::ecdf(ref_fixed)
  
  # apply fixed baseline to years <= min_years
  for (y in years_unique[years_unique <= min_years]) {
    target_idx = idx_by_year[[as.character(y)]]
    p[target_idx] = ecdf_fixed(x[target_idx])
  }
  
  # expanding reference for years > min_years
  for (y in years_unique[years_unique > min_years]) {
    base_years = seq_len(y)
    base_idx = unlist(idx_by_year[names(idx_by_year) %in% base_years])
    ref_vals = x[base_idx]
    ecdf_y = stats::ecdf(ref_vals)
    target_idx = idx_by_year[[as.character(y)]]
    p[target_idx] = ecdf_y(x[target_idx])
  }
  p
}

assign_drought_category_from_percentile = function(p, drought_bins_pct) {
  out = rep(NA_character_, length(p))
  for (i in seq_len(nrow(drought_bins_pct))) {
    in_bin = (p <= drought_bins_pct$upper[i]) & (p > drought_bins_pct$lower[i])
    out[in_bin] = as.character(drought_bins_pct$cat[i])
  }
  # D4 tail (<= 0.02)
  out[p <= 0.02] = "D4"
  factor(out, levels = levels(drought_bins_pct$cat))
}

summarize_drought_days = function(df_daily) {
  df_daily |>
    dplyr::filter(!is.na(drought_cat)) |>
    dplyr::count(year, drought_cat, name = "days") |>
    dplyr::group_by(year) |>
    dplyr::mutate(days = as.integer(days)) |>
    dplyr::ungroup()
}

summarize_drought_class_percent_lastN <- function(df_daily,
                                                  drought_bins_pct,
                                                  days_per_year,
                                                  years_total,
                                                  last_years = 10) {
  # which years to summarize
  start_year <- years_total - last_years + 1
  df_last <- df_daily |> dplyr::filter(year >= start_year)
  total_days <- last_years * days_per_year
  
  # count days by class in the last N years
  counts <- df_last |>
    dplyr::filter(!is.na(drought_cat)) |>
    dplyr::count(drought_cat, name = "days") |>
    dplyr::mutate(
      drought_cat = factor(drought_cat, levels = levels(drought_bins_pct$cat))
    ) |>
    tidyr::complete(
      drought_cat = factor(levels(drought_bins_pct$cat), levels = levels(drought_bins_pct$cat)),
      fill = list(days = 0L)
    ) |>
    dplyr::arrange(drought_cat)
  
  # ensure order from worst to best for cumulative sums
  # drought_bins_pct$cat is D0..D4; we want D4 worst first
  class_desc <- rev(levels(drought_bins_pct$cat)) # e.g., "D4","D3","D2","D1","D0"
  
  # create cumulative groups: D4, D3-D4, D2-D4, D1-D4, D0-D4
  cum_defs <- list(
    "D4"     = class_desc[1],
    "D3-D4"  = class_desc[1:2],
    "D2-D4"  = class_desc[1:3],
    "D1-D4"  = class_desc[1:4],
    "D0-D4"  = class_desc[1:5]
  )
  
  # compute days and pct for each cumulative bucket
  out <- purrr::imap_dfr(cum_defs, function(classes, label) {
    days_sum <- counts |> dplyr::filter(drought_cat %in% classes) |> dplyr::pull(days) |> sum(na.rm = TRUE)
    tibble::tibble(
      cum_label = label,
      days = as.integer(days_sum),
      pct_total = 100 * days_sum / total_days
    )
  })
  
  # order rows from worst to best (D4 first)
  out <- out |> dplyr::slice(match(names(cum_defs), out$cum_label))
  out
}

# ---- Scenario generators (daily TWS analogs) ----
scenario_A = function(df, seed = 1, days_per_year = 365) {
  set.seed(seed)
  seasonal = 1.0 * sin(2 * pi * df$day_of_year / days_per_year - pi / 2)
  noise = stats::rnorm(n = nrow(df), mean = 0, sd = sqrt(0.4))
  seasonal + noise
}

scenario_B = function(df, seed = 2, days_per_year = 365) {
  set.seed(seed)
  seasonal = 1.0 * sin(2 * pi * df$day_of_year / days_per_year - pi / 2)
  trend_per_year = -0.02
  trend = trend_per_year * (df$year - 1)
  noise = stats::rnorm(n = nrow(df), mean = 0, sd = sqrt(0.4))
  seasonal + trend + noise
}

scenario_C = function(df,
                      seed = 3,
                      days_per_year = 365,
                      decadal_period_years = 10,
                      decadal_amp = 0.5) {
  set.seed(seed)
  # --- base seasonal cycle (unchanged) ---
  seasonal =
    1.0 * sin(2 * pi * df$day_of_year / days_per_year - pi / 2)
  # --- decadal oscillation (low-frequency modulation) ---
  # Uses year_frac so oscillation is smooth through time
  decadal =
    decadal_amp * sin(2 * pi * df$year_frac / decadal_period_years)
  # --- stochastic noise ---
  noise =
    stats::rnorm(n = nrow(df), mean = 0, sd = sqrt(0.4))
  seasonal + decadal + noise
}

scenario_D = function(df, seed = 4, days_per_year = 365, years_total = 100) {
  set.seed(seed)
  seasonal = 1.0 * sin(2 * pi * df$day_of_year / days_per_year - pi / 2)
  sd0 = sqrt(0.4)
  years_to_double = 100
  sd_t = sd0 * 2 ^ ((df$year - 1) / years_to_double)
  noise = stats::rnorm(n = nrow(df), mean = 0, sd = sd_t)
  seasonal + noise
}

# ---- Build daily dataset and compute percentiles + classes ----
df_time = make_time_index(years_total = years_total, days_per_year = days_per_year)

scenarios = list(
  "A: Stationary Climate" = scenario_A(df = df_time, seed = seed, days_per_year = days_per_year),
  "B: Nonstationary Climate with Linear Trend" = scenario_B(df = df_time, seed = seed, days_per_year = days_per_year),
  "C: Nonstationary Climate with Decadal Oscillation" = scenario_C(df = df_time, seed = seed, days_per_year = days_per_year),
  "D: Nonstationary Climate with Increasing Variability Trend" = scenario_D(df = df_time, seed = seed, days_per_year = days_per_year, years_total = years_total)
)

df_all = purrr::imap_dfr(scenarios, function(tws_raw, scen) {
  
  # percentile relative to expanding climatology (min years)
  p = percentile_with_min_then_expanding_climatology(
    x = tws_raw,
    year = df_time$year,
    min_years = min_climatology_years
  )
  
  tibble::tibble(
    scenario = scen,
    year = df_time$year,
    day_of_year = df_time$day_of_year,
    year_frac = df_time$year_frac,
    tws_raw = tws_raw,
    tws_pct = p,                     # percentile in [0,1]
    tws_norm_plot = stats::qnorm(p), # for plotting as "standardized" index (NA where p==0 or 1)
    drought_cat = assign_drought_category_from_percentile(p = p, drought_bins_pct = drought_bins_pct)
  )
})

# ---- Annual drought-day summaries (for stacked bars) ----
df_days = df_all |>
  dplyr::group_by(scenario) |>
  dplyr::group_modify(~summarize_drought_days(.x)) |>
  dplyr::ungroup()

# ---- Plot builders ----
plot_timeseries = function(df_scen, label, days_per_year = 365) {
  
  # REQUIRE: df_scen has day_of_year (carry it through into df_all)
  df_plot = df_scen |>
    dplyr::mutate(
      seasonal_ref = sin(2 * pi * day_of_year / days_per_year - pi / 2)
    )
  
  # Fit trend ONLY on finite values (qnorm can create +/-Inf)
  df_fit = df_plot |>
    dplyr::filter(is.finite(tws_norm_plot), is.finite(year_frac))
  
  # Trend prediction across full record (set NA where can't compute)
  if (nrow(df_fit) < 2) {
    trend_hat = rep(0, nrow(df_plot))
  } else {
    fit = stats::lm(tws_norm_plot ~ year_frac, data = df_fit)
    trend_hat = stats::predict(fit, newdata = df_plot)
    trend_hat[!is.finite(trend_hat)] = NA_real_
  }
  
  # Seasonal wave that follows the fitted trend (B will tilt the seasonal line)
  df_plot = df_plot |>
    dplyr::mutate(seasonal_trended = seasonal_ref + trend_hat)
  
  ggplot2::ggplot(df_plot, ggplot2::aes(x = year_frac)) +
    
    # Daily drought index
    ggplot2::geom_line(
      ggplot2::aes(y = tws_norm_plot),
      linewidth = 0.25,
      alpha = 0.35,
      color = "#2F5D8A",
      na.rm = TRUE
    ) +
    
    # Seasonal reference signal riding on the trend
    ggplot2::geom_line(
      ggplot2::aes(y = seasonal_trended),
      linewidth = 0.25,
      color = "black",
      alpha = 0.6,
      na.rm = TRUE
    ) +
    
    # Linear trend line
    ggplot2::geom_smooth(
      ggplot2::aes(y = tws_norm_plot),
      method = "lm",
      linetype = "dashed",
      linewidth = 0.3,
      color = "grey40",
      se = FALSE,
      na.rm = TRUE
    ) +
    
    ggplot2::labs(x = "Years", y = "Drought Indicator", title = label) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank()
    )
}

plot_drought_bars = function(df_daily_scen, df_days_scen,
                             years_total, drought_bins_pct,
                             drought_colors, days_per_year,
                             last_years = 10) {
  
  # ---- stacked annual bars (ensure complete) ----
  df_complete =
    tidyr::expand_grid(
      year = seq_len(years_total),
      drought_cat = factor(levels(drought_bins_pct$cat), levels = levels(drought_bins_pct$cat))
    ) |>
    dplyr::left_join(df_days_scen, by = c("year", "drought_cat")) |>
    dplyr::mutate(days = tidyr::replace_na(days, 0L))
  
  # ---- percent-of-time summary over LAST `last_years` YEARS (cumulative) ----
  df_pct =
    summarize_drought_class_percent_lastN(
      df_daily = df_daily_scen |> dplyr::filter(!is.na(drought_cat)),
      drought_bins_pct = drought_bins_pct,
      days_per_year = days_per_year,
      years_total = years_total,
      last_years = last_years
    ) |>
    dplyr::mutate(
      # enforce desired order explicitly
      cum_label = factor(
        cum_label,
        levels = c("D0-D4", "D1-D4", "D2-D4", "D3-D4", "D4")
      )
    ) |>
    dplyr::arrange(cum_label) |>
    dplyr::mutate(
      label = sprintf("%s: %4.1f%%", cum_label, pct_total),
      
      color_key = dplyr::case_when(
        cum_label == "D4"    ~ "D4",
        cum_label == "D3-D4" ~ "D3",
        cum_label == "D2-D4" ~ "D2",
        cum_label == "D1-D4" ~ "D1",
        cum_label == "D0-D4" ~ "D0",
        TRUE ~ "D0"
      ),
      
      x = 2,
      y = 265 - (dplyr::row_number() - 1) * (days_per_year * 0.08)
    )
  
  ggplot2::ggplot(df_complete, ggplot2::aes(x = year, y = days, fill = drought_cat)) +
    ggplot2::geom_col(width = 0.9) +
    ggplot2::scale_fill_manual(values = drought_colors, drop = FALSE) +
    ggplot2::scale_y_continuous(limits = c(0, days_per_year), breaks = seq(0, days_per_year, 60)) +
    ggplot2::scale_x_continuous(limits = c(1, years_total), breaks = seq(0, years_total, 20)) +
    ggplot2::labs(x = "Year", y = "Days in Drought", fill = NULL) +
    
    # header text
    ggplot2::annotate(
      "text",
      x = 2,
      y = 300,
      label = paste0("Percent of Time in Drought\n(Last ", last_years, " Years)"),
      hjust = 0, vjust = 0,
      fontface = "bold", size = 3, color = "black"
    ) +
    
    # class percentages
    shadowtext::geom_shadowtext(
      data = df_pct,
      ggplot2::aes(x = x, y = y, label = label, color = color_key),
      inherit.aes = FALSE,
      hjust = 0,
      size = 3,
      fontface = "bold",
      bg.color = "black",   # outline color
      bg.r = 0.05           
    ) +
    ggplot2::scale_color_manual(values = drought_colors, guide = "none") +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "none",
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(t = 8, r = 40, b = 5, l = 5)
    )
}

# ---- Assemble 4x2 layout ----
scenario_order = c(
  "A: Stationary Climate",
  "B: Nonstationary Climate with Linear Trend",
  "C: Nonstationary Climate with Decadal Oscillation",
  "D: Nonstationary Climate with Increasing Variability Trend"
)

plots = purrr::map(scenario_order, function(scen) {
  df_s = df_all |> dplyr::filter(scenario == scen)
  df_d = df_days |> dplyr::filter(scenario == scen)
  
  p_left = plot_timeseries(df_scen = df_s, label = scen)
  p_right = plot_drought_bars(
    df_daily_scen = df_s,   # daily records (for last-10-year percent summary)
    df_days_scen  = df_d,   # annual aggregated days for stacked bars
    years_total   = years_total,
    drought_bins_pct = drought_bins_pct,
    drought_colors = drought_colors,
    days_per_year = days_per_year,
    last_years = 10
  )
  cowplot::plot_grid(p_left, p_right, ncol = 2, rel_widths = c(2.9, 1))
})

# ---- Shared legend ----
legend_plot = ggplot2::ggplot(
  tibble::tibble(`Drought Class` = factor(levels(drought_bins_pct$cat), levels = levels(drought_bins_pct$cat)),
                 x = 1, y = 1),
  ggplot2::aes(x = x, y = y, fill = `Drought Class`)
) +
  ggplot2::geom_col() +
  ggplot2::scale_fill_manual(values = drought_colors, drop = FALSE) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "right")

legend = cowplot::get_legend(legend_plot)

final = cowplot::plot_grid(
  cowplot::plot_grid(plotlist = plots, ncol = 1, align = "v")
)

# ---- Save ----
out_file = "~/nasem_examples/figs/nie_etal_example_expanded.png"
dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
ggplot2::ggsave(filename = out_file, plot = final, width = 12, height = 8, dpi = 300, bg = "white")

message("Saved: ", out_file)