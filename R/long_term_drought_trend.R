# ---- Libraries ----
library(sf)
library(terra)
library(ggplot2)
library(tidyterra)
library(patchwork)
library(tigris)
library(scales)

# ---- Local files ----
data_dir <- path.expand("~/nasem_examples/data")

f_def <- file.path(data_dir, "TC_p05_def_trend_30yr_CONUS.tif")
f_p   <- file.path(data_dir, "TC_p05_p_trend_30yr_CONUS.tif")
f_pet <- file.path(data_dir, "TC_p95_pet_trend_30yr_CONUS.tif")

stopifnot(file.exists(f_def), file.exists(f_p), file.exists(f_pet))

# ---- Read rasters ----
r_def <- terra::rast(f_def)
r_p   <- terra::rast(f_p)
r_pet <- terra::rast(f_pet)

names(r_def) <- "def"
names(r_p)   <- "p"
names(r_pet) <- "pet"

# ---- Project to EPSG:5070 (NAD83 / CONUS Albers) ----
# Use bilinear for smoother continuous fields; use "near" if you want exact values preserved.
r_def_5070 <- terra::project(r_def, "EPSG:5070", method = "bilinear")
r_p_5070   <- terra::project(r_p,   "EPSG:5070", method = "bilinear")
r_pet_5070 <- terra::project(r_pet, "EPSG:5070", method = "bilinear")

# ---- Palettes ----
pal_rywcb <- c("#b30000", "#ffff00", "#ffffff", "#00ffff", "#0000b3") # P & DEF
pal_rev   <- rev(pal_rywcb)                                           # PET

# ---- Map range ----
vmin <- -5
vmax <-  5

# ---- State lines (EPSG:5070) ----
st_lines <- read_sf('https://eric.clst.org/assets/wiki/uploads/Stuff/gz_2010_us_040_00_20m.json') |>
  filter(!NAME %in% c('Virgin Island', 'Hawaii', 'Alaska', 'Puerto Rico')) |>
  sf::st_as_sf() |>
  sf::st_transform(5070)

# ---- Plot helper ----
# ---- Plot helper ----
plot_trend <- function(r, title, pal, vmin = -5, vmax = 5) {
  ggplot() +
    tidyterra::geom_spatraster(
      data = r,
      na.rm = FALSE
    ) +
    geom_sf(data = st_lines, fill = NA, color = "grey25", linewidth = 0.25) +
    coord_sf(crs = sf::st_crs(5070)) +
    scale_fill_gradientn(
      colors = pal,
      limits = c(vmin, vmax),
      oob = scales::squish,
      na.value = NA,
      name = "30-Year Moving Percentile Trend (TerraClimate 1958–2024)",
      guide = guide_colorbar(
        title.position = "top",
        direction = "horizontal",
        barwidth = unit(10, "cm"),
        barheight = unit(0.4, "cm")
      )
    ) +
    labs(
      title = title
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.title = element_blank(),
      axis.text  = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5),
      legend.position = "bottom",
      legend.justification = "center",
      legend.title = element_text(hjust = 0.5),
      legend.key.width  = unit(6, "cm"),
      legend.key.height = unit(0.4, "cm"),
      plot.margin = margin(10, 5, 5, 5)
    )
}
p1 <- plot_trend(r_p_5070,   "5th Percentile Precipitaion (P) Drought Trend", pal_rywcb, vmin = -5, vmax = 5)
p2 <- plot_trend(r_pet_5070, "95th Percentile Reference ET (ETₒ) Drought Trend", pal_rev, vmin = -3, vmax = 3)
p3 <- plot_trend(r_def_5070, "5th Percentile Water Deficit (P - ETₒ) Drought Trend", pal_rywcb, vmin = -5, vmax = 5)


library(cowplot)

final_plot <- cowplot::plot_grid(
  p1 + theme(legend.position = "bottom"),
  p2 + theme(legend.position = "bottom"),
  p3 + theme(legend.position = "bottom"),
  ncol = 1,
  align = "v",
  rel_heights = c(1, 1, 1)
)

ggsave(
  filename = "~/nasem_examples/figs/drought_trends.png",
  plot = final_plot,
  width = 8,
  height = 14,
  bg = "white"
)

