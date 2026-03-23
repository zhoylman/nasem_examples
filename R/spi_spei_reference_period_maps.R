################################################################################
# Title: SPI & SPEI Reference Period Comparison Maps
# Description: Creates two figures comparing contemporary (30-year) vs full
#              period-of-record baselines:
#   Figure 9-1: 90-day SPI over the mid-south US
#   Figure 9-2: 90-day SPEI over Nebraska
#
# Inputs (GeoTIFFs):
#   ~/nasem_examples/data/
#     - spi_90d_rolling-30_2026-03-18.tif
#     - spi_90d_full_2026-03-18.tif
#     - spei_90d_rolling-30_2026-03-18.tif
#     - spei_90d_full_2026-03-18.tif
#
# Output (PNG):
#   ~/nasem_examples/figs/fig_9-1_spi_reference_periods.png
#   ~/nasem_examples/figs/fig_9-2_spei_reference_periods.png
#
# Author: Zachary H. Hoylman
# Date: 3-21-2026
################################################################################

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(terra)
  library(ggplot2)
  library(tidyterra)
  library(cowplot)
  library(tigris)
})

options(tigris_use_cache = TRUE)

# ------------------------------------------------------------------------------
# PATHS
# ------------------------------------------------------------------------------
DATA_DIR = path.expand("~/nasem_examples/data")
OUT_DIR  = path.expand("~/nasem_examples/figs")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# READ RASTERS
# ------------------------------------------------------------------------------
spi_rolling30  = terra::rast(file.path(DATA_DIR, "spi_90d_rolling-30_2026-03-18.tif"))
spi_full       = terra::rast(file.path(DATA_DIR, "spi_90d_full_2026-03-18.tif"))
spei_rolling30 = terra::rast(file.path(DATA_DIR, "spei_90d_rolling-30_2026-03-18.tif"))
spei_full      = terra::rast(file.path(DATA_DIR, "spei_90d_full_2026-03-18.tif"))

names(spi_rolling30)  = "value"
names(spi_full)       = "value"
names(spei_rolling30) = "value"
names(spei_full)      = "value"

# ------------------------------------------------------------------------------
# BOUNDARIES
# ------------------------------------------------------------------------------
message("Downloading county and state boundaries...")
counties_sf = tigris::counties(cb = TRUE, resolution = "20m", year = 2022) |>
  sf::st_transform(4326)

states_sf = tigris::states(cb = TRUE, resolution = "20m", year = 2022) |>
  dplyr::filter(!STUSPS %in% c("HI", "AK", "PR", "VI", "GU", "AS", "MP")) |>
  sf::st_transform(4326)

# ------------------------------------------------------------------------------
# USDM-STYLE CATEGORICAL CLASSIFICATION
# ------------------------------------------------------------------------------
usdm_labels = c(
  "D4\n< -2.05", "D3\n-2.05 \u2013 -1.64", "D2\n-1.64 \u2013 -1.28",
  "D1\n-1.28 \u2013 -0.84", "D0\n-0.84 \u2013 -0.52",
  "Normal\n-0.52 \u2013 0.52",
  "W0\n0.52 \u2013 0.84", "W1\n0.84 \u2013 1.28",
  "W2\n1.28 \u2013 1.64", "W3\n1.64 \u2013 2.05", "W4\n\u2265 2.05"
)

usdm_colors = c(
  "#730000",
  "#e60001",
  "#ffaa01",
  "#fdd37f",
  "#ffff05",
  "#f5f5f5",
  "#b2e5f6",
  "#62c8f5",
  "#1d6fe3",
  "#1232a8",
  "#0c1a5e"
)

classify_usdm = function(r) {
  rcl = matrix(c(
    -Inf,   -2.054, 1,
    -2.054, -1.644, 2,
    -1.644, -1.281, 3,
    -1.281, -0.842, 4,
    -0.842, -0.524, 5,
    -0.524,  0.524, 6,
     0.524,  0.842, 7,
     0.842,  1.281, 8,
     1.281,  1.644, 9,
     1.644,  2.054, 10,
     2.054,  Inf,   11
  ), ncol = 3, byrow = TRUE)
  classified = terra::classify(r, rcl, include.lowest = TRUE)
  levels(classified) = data.frame(id = 1:11, category = usdm_labels)
  classified
}

# ------------------------------------------------------------------------------
# PLOT FUNCTION
# ------------------------------------------------------------------------------
make_panel = function(raster, counties, states, xlim, ylim, subtitle) {

  r_cat = classify_usdm(raster)

  ggplot() +
    tidyterra::geom_spatraster(data = r_cat) +
    geom_sf(data = counties, fill = NA, color = "grey50", linewidth = 0.08) +
    geom_sf(data = states, fill = NA, color = "grey20", linewidth = 0.4) +
    coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
    scale_fill_manual(
      values = setNames(usdm_colors, usdm_labels),
      na.value = "grey90",
      drop = FALSE,
      guide = "none"
    ) +
    labs(title = subtitle) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(size = 18, hjust = 0.5, margin = margin(t = 0, b = 2)),
      axis.text  = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "grey30", linewidth = 0.5),
      plot.margin = margin(0, 4, 0, 4)
    )
}

# ------------------------------------------------------------------------------
# LEGEND BAR
# ------------------------------------------------------------------------------
make_legend_bar = function() {
  legend_df = data.frame(
    x = 1:11,
    fill = factor(usdm_labels, levels = usdm_labels)
  )

  ggplot(legend_df, aes(x = x, y = 1, fill = fill)) +
    geom_tile(width = 1, height = 2, color = "grey40", linewidth = 0.3) +
    scale_fill_manual(values = setNames(usdm_colors, usdm_labels), drop = FALSE) +
    scale_x_continuous(
      breaks = 1:11,
      labels = usdm_labels,
      expand = c(0.02, 0)
    ) +
    theme_void(base_size = 14) +
    theme(
      axis.text.x = element_text(color = "black", size = 10, angle = 30,
                                 hjust = 1, vjust = 1, lineheight = 0.85,
                                 face = "bold",
                                 margin = margin(t = 4)),
      legend.position = "none",
      plot.margin = margin(0, 20, 6, 20)
    )
}

# ------------------------------------------------------------------------------
# HELPER: assemble a two-panel figure with title, subtitle, and legend
# ------------------------------------------------------------------------------
make_figure = function(p_left, p_right, main_title, sub_title,
                       legend_rel = 0.2) {

  title_grob = ggdraw() +
    draw_label(main_title, fontface = "bold", size = 24, hjust = 0.5)

  subtitle_grob = ggdraw() +
    draw_label(sub_title, fontface = "italic", size = 14, hjust = 0.5,
               color = "grey30")

  maps = plot_grid(p_left, p_right, ncol = 2, align = "hv")
  legend = make_legend_bar()

  plot_grid(
    title_grob,
    subtitle_grob,
    maps,
    legend,
    ncol = 1,
    rel_heights = c(0.08, 0.07, 1, legend_rel)
  )
}

# ==============================================================================
# FIGURE 9-1: SPI — Mid-South US
# ==============================================================================
message("Creating Figure 9-1: SPI Reference Period Comparison...")

spi_xlim = c(-94.5, -83.5)
spi_ylim = c(33.0, 39.5)

p_spi_contemp = make_panel(
  raster = spi_rolling30, counties = counties_sf, states = states_sf,
  xlim = spi_xlim, ylim = spi_ylim,
  subtitle = "Contemporary Baseline (1997\u20132026)"
)

p_spi_full = make_panel(
  raster = spi_full, counties = counties_sf, states = states_sf,
  xlim = spi_xlim, ylim = spi_ylim,
  subtitle = "Full Record Baseline (1979\u20132026)"
)

fig_9_1 = make_figure(
  p_spi_contemp, p_spi_full,
  main_title = "90-day SPI on March 18, 2026",
  sub_title  = "Standardized Precipitation Index \u2014 Mid-South United States",
  legend_rel = 0.3
)

ggsave(
  filename = file.path(OUT_DIR, "fig_9-1_spi_reference_periods.png"),
  plot = fig_9_1, width = 10, height = 6, dpi = 300, bg = "white"

)
message("Saved: fig_9-1_spi_reference_periods.png")

# ==============================================================================
# FIGURE 9-2: SPEI — Nebraska
# ==============================================================================
message("Creating Figure 9-2: SPEI Reference Period Comparison...")

ne_bbox   = states_sf |> filter(STUSPS == "NE") |> st_bbox()
spei_xlim = c(ne_bbox["xmin"] - 0.5, ne_bbox["xmax"] + 0.5)
spei_ylim = c(ne_bbox["ymin"] - 0.5, ne_bbox["ymax"] + 0.5)

p_spei_contemp = make_panel(
  raster = spei_rolling30, counties = counties_sf, states = states_sf,
  xlim = spei_xlim, ylim = spei_ylim,
  subtitle = "Contemporary Baseline (1997\u20132026)"
)

p_spei_full = make_panel(
  raster = spei_full, counties = counties_sf, states = states_sf,
  xlim = spei_xlim, ylim = spei_ylim,
  subtitle = "Full Record Baseline (1979\u20132026)"
)

fig_9_2 = make_figure(
  p_spei_contemp, p_spei_full,
  main_title = "90-day SPEI on March 18, 2026",
  sub_title  = "Standardized Precipitation Evapotranspiration Index \u2014 Nebraska",
  legend_rel = 0.3
)

ggsave(
  filename = file.path(OUT_DIR, "fig_9-2_spei_reference_periods.png"),
  plot = fig_9_2, width = 10, height = 5.5, dpi = 300, bg = "white"
)
message("Saved: fig_9-2_spei_reference_periods.png")

message("Done!")
