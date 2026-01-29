library(httr)
library(dplyr)
library(data.table)
library(terra)
library(sf)
library(dplyr)
library(ggplot2)
library(ggspatial)
library(scales)
library(exactextractr)

export.dir = path.expand('~/temp/snodas/')

my_simplify <-
  function(x, 
           keep = 0.05,
           quiet = TRUE){
    temp.js <- tempfile(fileext = ".geojson")
    
    sf::write_sf(x, temp.js)
    
    call <- rmapshaper:::make_simplify_call(keep = keep, method = NULL, 
                                            weighting = 0.7, keep_shapes = FALSE, no_repair = FALSE, 
                                            snap = TRUE, explode = FALSE, drop_null_geometries = TRUE, 
                                            snap_interval = NULL)
    command <- paste(rmapshaper:::ms_compact(call), collapse = " ")
    
    sys_mem = getOption("mapshaper.sys_mem", 
                        default = 8)
    
    out.js <- tempfile(fileext = ".geojson")
    
    cmd_args <- c(sys_mem, shQuote(temp.js), command, shQuote(NULL), "-o", shQuote(out.js), if (quiet) "-quiet")
    
    ms_path <- paste0(rmapshaper:::check_sys_mapshaper("mapshaper-xl", verbose = FALSE))
    
    system2(ms_path, cmd_args)
    
    sf::read_sf(out.js)
  }

get_snodas = function(date){
  # Create needed directories if they don't exist
  dir.create(file.path(export.dir, "snodas/raw"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(export.dir, "snodas/processed/swe"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(export.dir, "snodas/processed/snow_depth"), recursive = TRUE, showWarnings = FALSE)
  
  for(d in seq_along(date)){
    tryCatch({
      # Build data URL
      url = paste0("ftp://sidads.colorado.edu/DATASETS/NOAA/G02158/masked/", 
                   format(date[d], "%Y"), "/", format(date[d], "%m_%b"), 
                   "/SNODAS_", format(date[d], "%Y%m%d"), ".tar")
      
      # Define where raw tarball will be stored
      tar.dir = file.path(export.dir, "snodas/raw", paste0("SNODAS_", format(date[d], "%Y%m%d"), ".tar"))
      unzip.dir = file.path(export.dir, "snodas/raw", paste0("SNODAS_", format(date[d], "%Y%m%d")))
      
      # Download zipped data
      httr::GET(url, write_disk(path = tar.dir, overwrite=TRUE))
      
      # Create unzip location
      dir.create(unzip.dir, recursive = TRUE, showWarnings = FALSE)
      
      # Unzip file
      untar(tarfile = tar.dir,  exdir = unzip.dir)
      
      # Get files from unzipped dir
      files = list.files(unzip.dir, full.names = TRUE) 
      
      # NOAA variable IDs
      files_of_interest = c("1034", "1036")
      
      for(i in seq_along(files_of_interest)) {
        file_to_process = files[which(files %like% files_of_interest[i] & files %like% ".dat.gz")]
        
        writeLines(
          "ENVI
samples = 6935
lines   = 3351
bands   = 1
header offset = 0
file type = ENVI Standard
data type = 2
interleave = bsq
byte order = 1", 
          con = gsub(".dat.gz", ".hdr", file_to_process)
        )
        
        R.utils::gunzip(file_to_process, destname = gsub(".gz", "", file_to_process)) 
        
        processed_name = if (files_of_interest[i] == "1034") {
          file.path(export.dir, "snodas/processed/swe", paste0("snodas_swe_conus_", format(date[d], "%Y%m%d"), ".tif"))
        } else {
          file.path(export.dir, "snodas/processed/snow_depth", paste0("snodas_snow_depth_conus_", format(date[d], "%Y%m%d"), ".tif"))
        }
        
        # Different file dimensions after Oct. 2013
        translate_cmd = if(date[d] > as.Date('2013-10-01')){
          "gdal_translate -of GTiff -a_srs '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs' -a_nodata -9999 -a_ullr -124.73333333 52.87500000 -66.94166667 24.95000000"
        } else {
          "gdal_translate -of GTiff -a_srs '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs' -a_nodata -9999 -a_ullr -124.73375000000000 52.87458333333333 -66.94208333333333 24.94958333333333"
        }
        
        system(paste(translate_cmd, gsub(".gz", "", file_to_process), processed_name))
      }
      
      # Clean up raw data
      unlink(c(tar.dir, unzip.dir), recursive = TRUE)
      
    }, error = function(e){
      # Clean up even on error
      unlink(c(tar.dir, unzip.dir), recursive = TRUE)
    })
  }
}

today = Sys.Date()

years = 2004:as.integer(format(today, "%Y"))

dates = as.Date(
  paste0(
    years, "-",
    format(today, "%m-%d")
  )
)

get_snodas(dates)

# ---- Inputs ----
files = base::sort(list.files(glue::glue("{export.dir}/snodas/processed/swe"), full.names = TRUE))

conus = read_sf('https://eric.clst.org/assets/wiki/uploads/Stuff/gz_2010_us_040_00_20m.json') |>
  filter(!NAME %in% c('Hawaii', 'Alaska', 'Virgin Islands', 'Puerto Rico'))

#huc 8 from: https://www.hydroshare.org/resource/b832a6c2f96541808444ec9562c5247e/
huc8 = sf::read_sf("~/Downloads/HUC8_CONUS/HUC8_US.shp")

# Watershed name column in your shapefile
name_col = "NAME"

# Date tag to match your target raster in the filename
target_tag = stringr::str_remove_all(today, '-')

# ---- Baseline vs target ----
baseline_files = files[!grepl(target_tag, files) & grepl(substr(target_tag, 5, 8) , files)]
target_file    = files[ grepl(target_tag, files) & grepl(substr(target_tag, 5, 8) , files)]

target = terra::rast(target_file)

# Project HUC8s to raster CRS (critical for correct extraction)
huc8 = sf::st_transform(huc8, terra::crs(target))

# ---- Build baseline median raster on the SAME grid as target ----
baseline_list = lapply(
  baseline_files,
  function(f) terra::resample(terra::rast(f), target, method = "bilinear")
)
baseline_stack = terra::rast(baseline_list)
baseline_median = terra::app(baseline_stack, median, na.rm = TRUE)

# ---- Clean / guardrails ----
target = terra::clamp(target, lower = 0, values = TRUE)
baseline_median = terra::clamp(baseline_median, lower = 0, values = TRUE)
baseline_median = terra::ifel(baseline_median == 0, NA, baseline_median)

# ---- Basin-wide % of normal (ratio of sums) ----
base_sum = exactextractr::exact_extract(baseline_median, huc8, fun = "sum")
curr_sum = exactextractr::exact_extract(target,         huc8, fun = "sum")

huc8_sf = huc8 |>
  dplyr::mutate(
    base_sum = base_sum,
    curr_sum = curr_sum,
    value = (curr_sum / base_sum) * 100,
    value = dplyr::if_else(is.finite(value), value, NA_real_),
    value = pmax(value, 0)
  )

# ---- Discrete bins (match provided legend) ----
# Bins match the posted legend (0-29,30-49,50-69,70-89,90-109,110-129,130-149,>=150)
bin_breaks = c(0, 30, 50, 70, 90, 110, 130, 150, Inf)
bin_labels = c(
  "≤ 30%",
  "30 to 50%",
  "50 to 70%",
  "70 to 90%",
  "90 to 110%",
  "110 to 130%",
  "130 to 150%",
  "≥ 150%"
)

# ---- Colors matching your legend (top -> bottom: ≥150 -> ... -> 0-29, then no-basin grey) ----
# I picked hex values close to the sample you provided.
palette_named = c(
  "≥ 150%"     = "#2166ac",  # blue
  "130 to 150%"= "#3b8ecf",  # cyan-blue
  "110 to 130%"= "#66c2a5",  # aqua
  "90 to 110%" = "#c7e9c0",  # light mint/green
  "70 to 90%"  = "#ffff99",  # pale yellow
  "50 to 70%"  = "#fdae61",  # orange
  "30 to 50%"  = "#f46d43",  # dark orange
  "≤ 30%"   = "#8b2500"  # brownish/dark-red
)

huc8_sf = huc8_sf |>
  dplyr::mutate(
    pct_bin = cut(value, breaks = bin_breaks, labels = bin_labels, right = FALSE, include.lowest = TRUE),
    pct_bin = as.character(pct_bin),
    pct_bin = ifelse(is.na(value), "No basin value", pct_bin),
    pct_bin = factor(pct_bin, levels = c(rev(bin_labels), "No basin value")) # reverse so >=150 shows first in legend
  ) |>
  dplyr::mutate(
    hex_color = palette_named[as.character(pct_bin)]
  )

# Make sure palette covers all factor levels
palette_named = palette_named[levels(huc8_sf$pct_bin)]

# ---- Labels (watershed name + percent) ----
label_pts = sf::st_point_on_surface(huc8_sf) |>
  sf::st_as_sf() |>
  dplyr::mutate(
    ws_name = as.character(.data[[name_col]]),
    pct_lbl = dplyr::if_else(is.na(value), NA_character_, paste0(round(value), "%")),
    lbl = dplyr::if_else(is.na(value), NA_character_, paste0(ws_name, "\n", pct_lbl))
  )

# ---- Plot (discrete bins) ----
p = ggplot2::ggplot() +
  ggplot2::geom_sf(
    data = huc8_sf,
    ggplot2::aes(fill = pct_bin),
    color = "black",
    linewidth = 0.25
  ) +
  ggplot2::geom_sf(
    data = conus,
    color = "#3b3b3b",
    fill = 'transparent',
    linewidth = 1
  ) +
  ggplot2::scale_fill_manual(
    values = palette_named,
    drop = FALSE,
    na.value = "grey90",
    name = "Basin-wide SWE\nPercent of Normal",
    guide = ggplot2::guide_legend(reverse = FALSE, ncol = 1) # levels already reversed above
  ) +
  ggspatial::annotation_north_arrow(
    location = "bl",
    which_north = "true",
    style = ggspatial::north_arrow_fancy_orienteering,
    pad_x = grid::unit(5, "cm"),
    pad_y = grid::unit(1, "cm")
  ) +
  ggplot2::labs(
    title = "Snow Data Assimilation System (SNODAS) Snow Water Equivalent (SWE)",
    subtitle = "Basin-wide percent of normal (baseline median 2004–2025) — January 26, 2026",
    caption = "SNODAS SWE summary by the Montana Climate Office\n Basin-wide % normal = 100 * (sum current SWE) / (sum baseline median SWE)"
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    legend.position = "right",
    legend.key.height = grid::unit(0.6, "cm"),
    legend.title = ggplot2::element_text(size = 10),
    legend.text = ggplot2::element_text(size = 9),
    plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 18),
    plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 14),
    axis.title = ggplot2::element_blank(),
    axis.text = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    legend.box.spacing = grid::unit(-2, "cm")
  )+
  ggplot2::theme(
    plot.background  = ggplot2::element_rect(fill = "white", color = NA),
    panel.background = ggplot2::element_rect(fill = "white", color = NA)
  )

# ---- Save outputs ----
ggplot2::ggsave(
  filename = "~/temp/snodas_basinwide_percent_normal_discrete.png",
  plot = p,
  width = 11,
  height = 8,
  dpi = 300
)

final_geojson = huc8_sf |>
  select(HUC8, NAME, value, pct_bin, hex_color)

write_sf(my_simplify(final_geojson), '~/temp/current_snodas_summary.fgb')
write_lines(format(today, "%m-%d-%Y"), '~/temp/current_snodas_time.txt')