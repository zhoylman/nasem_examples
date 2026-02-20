################################################################################
# Title: TerraClimate 30-Year Moving Percentile Trends — CONUS Maps (3x2)
# Description: Reads locally-downloaded GeoTIFF trend rasters (P, PET, DEF, VPD,
#              TMMN, TMMX), reprojects to EPSG:5070, overlays state boundaries,
#              and generates a 3-row x 2-column panel figure with per-panel
#              legends (horizontal), transparent NA, and units in legend titles.
#
# Inputs (GeoTIFFs):
#   ~/nasem_examples/data/
#     - TC_p05_p_trend_30yr_CONUS.tif
#     - TC_p95_pet_trend_30yr_CONUS.tif
#     - TC_p05_def_trend_30yr_CONUS.tif
#     - TC_p95_vpd_trend_30yr_CONUS.tif
#     - TC_p95_tmmn_trend_30yr_CONUS.tif
#     - TC_p95_tmmx_trend_30yr_CONUS.tif
#
# Output (PNG):
#   ~/nasem_examples/figs/drought_trends_3x2.png
#
# Author: Zachary H. Hoylman
# Date: 2-20-2026
################################################################################

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(terra)
  library(ggplot2)
  library(tidyterra)
  library(scales)
  library(cowplot)
})

# ------------------------------------------------------------------------------
# PATHS
# ------------------------------------------------------------------------------
DATA_DIR = path.expand("~/nasem_examples/data")
OUT_DIR  = path.expand("~/nasem_examples/figs")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

FILES = list(
  def  = file.path(DATA_DIR, "TC_p05_def_trend_30yr_CONUS.tif"),
  p    = file.path(DATA_DIR, "TC_p05_p_trend_30yr_CONUS.tif"),
  pet  = file.path(DATA_DIR, "TC_p95_pet_trend_30yr_CONUS.tif"),
  vpd  = file.path(DATA_DIR, "TC_p95_vpd_trend_30yr_CONUS.tif"),
  tmmn = file.path(DATA_DIR, "TC_p95_tmmn_trend_30yr_CONUS.tif"),
  tmmx = file.path(DATA_DIR, "TC_p95_tmmx_trend_30yr_CONUS.tif")
)

stopifnot(all(file.exists(unlist(FILES))))

# ------------------------------------------------------------------------------
# READ + PROJECT (EPSG:5070)
# ------------------------------------------------------------------------------
PROJ_5070 = "EPSG:5070"

read_and_project = function(path, band_name, method = "bilinear") {
  r = terra::rast(path)
  names(r) = band_name
  terra::project(r, PROJ_5070, method = method)
}

r_def_5070  = read_and_project(FILES$def,  "def")
r_p_5070    = read_and_project(FILES$p,    "p")
r_pet_5070  = read_and_project(FILES$pet,  "pet")
r_vpd_5070  = read_and_project(FILES$vpd,  "vpd")
r_tmmn_5070 = read_and_project(FILES$tmmn, "tmmn")
r_tmmx_5070 = read_and_project(FILES$tmmx, "tmmx")

# ------------------------------------------------------------------------------
# BASEMAP (states) — EPSG:5070
# ------------------------------------------------------------------------------
st_lines =
  sf::read_sf("https://eric.clst.org/assets/wiki/uploads/Stuff/gz_2010_us_040_00_20m.json") |>
  dplyr::filter(!NAME %in% c("Virgin Islands", "Hawaii", "Alaska", "Puerto Rico")) |>
  sf::st_make_valid() |>
  sf::st_transform(PROJ_5070)

# ------------------------------------------------------------------------------
# COLORS + PLOT HELPERS
# ------------------------------------------------------------------------------
pal_rywcb = c("#b30000", "#ffff00", "#ffffff", "#00ffff", "#0000b3") # P/DEF/etc.
pal_rev   = rev(pal_rywcb)                                           # PET (and VPD if desired)

plot_trend = function(r, title, pal, units_label, vmin, vmax) {
  ggplot() +
    tidyterra::geom_spatraster(data = r, na.rm = FALSE) +
    geom_sf(data = st_lines, fill = NA, color = "grey25", linewidth = 0.25) +
    coord_sf(crs = sf::st_crs(5070)) +
    scale_fill_gradientn(
      colors = pal,
      limits = c(vmin, vmax),
      oob = scales::squish,
      na.value = NA,  # transparent
      name = paste0("Trend (", units_label, "/yr)"),
      guide = guide_colorbar(
        title.position = "top",
        direction = "horizontal",
        barwidth = unit(10, "cm"),
        barheight = unit(0.4, "cm")
      )
    ) +
    labs(title = title) +
    theme_minimal(base_size = 11) +
    theme(
      axis.title = element_blank(),
      axis.text  = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      legend.position = "bottom",
      legend.justification = "center",
      legend.title = element_text(hjust = 0.5),
      legend.key.width  = unit(6, "cm"),
      legend.key.height = unit(0.4, "cm"),
      plot.margin = margin(8, 6, 6, 6),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
}

# ------------------------------------------------------------------------------
# PANELS (edit vmin/vmax as desired)
# ------------------------------------------------------------------------------
p1 = plot_trend(
  r = r_p_5070,
  title = "p05 Precipitation (P)",
  pal = pal_rywcb,
  units_label = "mm",
  vmin = -5, vmax = 5
)

p2 = plot_trend(
  r = r_pet_5070,
  title = "p95 Reference ET (ET\u2092)",   # ETₒ
  pal = pal_rev,
  units_label = "mm",
  vmin = -3, vmax = 3
)

p3 = plot_trend(
  r = r_def_5070,
  title = "p05 Water Deficit (P \u2212 ET\u2092)",  # P − ETₒ
  pal = pal_rywcb,
  units_label = "mm",
  vmin = -5, vmax = 5
)

p4 = plot_trend(
  r = r_vpd_5070,
  title = "p95 Vapor Pressure Deficit (VPD)",
  pal = pal_rev,
  units_label = "kPa",
  vmin = -0.005, vmax = 0.005
)

p5 = plot_trend(
  r = r_tmmn_5070,
  title = "p95 Minimum Temperature (Tmin)",
  pal = pal_rev,
  units_label = "\u00B0C",  # °C
  vmin = -0.05, vmax = 0.05
)

p6 = plot_trend(
  r = r_tmmx_5070,
  title = "p95 Maximum Temperature (Tmax)",
  pal = pal_rev,
  units_label = "\u00B0C",  # °C
  vmin = -0.05, vmax = 0.05
)

# ------------------------------------------------------------------------------
# 3 ROWS x 2 COLUMNS (cowplot) — legends kept per panel
# ------------------------------------------------------------------------------
final_plot = cowplot::plot_grid(
  p5, p6,
  p1, p4,
  p2, p3,
  ncol = 2,
  align = "hv"
)

# ------------------------------------------------------------------------------
# SAVE (white background)
# ------------------------------------------------------------------------------
ggsave(
  filename = file.path(OUT_DIR, "drought_trends_3x2.png"),
  plot = final_plot,
  width = 12,
  height = 14,
  bg = "white",
  dpi = 300
)
