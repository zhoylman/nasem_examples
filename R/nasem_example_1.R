# --- 1. LIBRARIES & ENVIRONMENT ---
library(reticulate)
library(sf)
library(tidyverse)
library(rgee)
library(slider)
library(purrr)
library(lmomco)
library(viridis)
library(cowplot)

# Use your specific conda environment
use_condaenv("gee_env", required = TRUE)

# Initialize with project ID
ee_Initialize(user = 'zhoylman@gmail.com', drive = TRUE, project = 'ee-zhoylman')
ee$Initialize(project = 'ee-zhoylman')

# Manual session fix for rgee (macOS/Linux bug)
session_file = "/Users/zachary.hoylman/.config/earthengine/rgee_sessioninfo.txt"
info_content = paste('"user"', '"drive_cre"', '"gcs_cre"', '"zhoylman@gmail.com"', 
                     '"~/.config/googledrive/zhoylman@gmail.com"', 'NA', sep = " ")
writeLines(info_content, con = session_file)

# Source external drought functions
source('https://raw.githubusercontent.com/mt-climate-office/mco-drought-indicators/refs/heads/master/processing/ancillary-functions/R/drought-functions.R')

# --- 2. GEE DATA EXTRACTION (3-MONTH JAS) ---
# Define ROI: Blaine County, Idaho
ROI_raw = read_sf('https://eric.clst.org/assets/wiki/uploads/Stuff/gz_2010_us_050_00_20m.json') |>
  filter(NAME == 'Blaine' & STATE == "16") |>
  st_union() 

ROI = ROI_raw |>
  sf_as_ee()

# Extraction parameters
target_bands = c("pr", "pet")
years_ee = ee$List$sequence(1958, 2024)
terra_ic = ee$ImageCollection("IDAHO_EPSCOR/TERRACLIMATE")

# Updated GEE function for Seasonal Summation (Jul-Sep)
get_seasonal_stat = function(year, dataset, band_name) {
  img = dataset$select(band_name)$
    filter(ee$Filter$calendarRange(year, year, "year"))$
    filter(ee$Filter$calendarRange(7, 9, "month"))$
    sum() # Seasonal Total
  
  stat = img$reduceRegion(reducer = ee$Reducer$median(), geometry = ROI, 
                          scale = 4638, maxPixels = 1e9)
  
  return(ee$Feature(NULL, list(year = year, value = stat$get(band_name), variable = band_name)))
}

# Execute extraction
annual_raw = map_dfr(target_bands, function(band) {
  message(paste("Extracting JAS Seasonal:", band))
  mapped = years_ee$map(ee_utils_pyfunc(function(y) get_seasonal_stat(y, terra_ic, band)))
  ee_as_sf(ee$FeatureCollection(mapped), via = "getInfo") |>
    as_tibble() |>
    select(year, value, variable) |>
    mutate(value = if_else(variable == 'pet', value * 0.1, value))
})

# --- 3. CALCULATIONS: DEFICIT & ROLLING NORMALS ---
annual_combined = annual_raw |>
  pivot_wider(names_from = variable, values_from = value) |>
  mutate(def = pr - pet) |>
  pivot_longer(cols = c(pr, pet, def), names_to = "variable")

rolling_data = annual_combined |>
  arrange(variable, year) |>
  group_by(variable) |>
  mutate(rolling_mean_30yr = slide_index_dbl(.x = value, .i = year, .f = mean, 
                                             .before = 29, .complete = TRUE)) |>
  ungroup()

# Calculate D3 Threshold Shift for JAS Deficit
def_subset = annual_combined |> filter(variable == "def")

# First Period (1961-1990)
p1_data = def_subset |> filter(year <= 1990, year > 1960) |> pull(value)
fit_p1 = parglo(pwm2lmom(pwm.ub(p1_data)))
d3_threshold = quaglo(0.049, fit_p1) # Value of the 4th percentile

# Last Period (1991-2020)
p2_data = def_subset |> filter(year <= 2020, year > 1990) |> pull(value)
fit_p2 = parglo(pwm2lmom(pwm.ub(p2_data)))
new_prob = cdfglo(d3_threshold, fit_p2) # New percentile for that same value

# Generate PDF "Slinky" data for Deficit
get_density_curve = function(target_year, df) {
  window_data = df |> filter(year <= target_year, year > (target_year - 30)) |> pull(value)
  if(length(window_data) < 30) return(NULL)
  
  fit_glo = parglo(pwm2lmom(pwm.ub(window_data)))
  x_range = seq(min(df$value) - 100, max(df$value) + 100, length.out = 200) # Wider range for seasonal totals
  return(tibble(x = x_range, y = pdfglo(x_range, fit_glo), year = target_year))
}

slinky_data = map_dfr(unique(rolling_data$year[!is.na(rolling_data$rolling_mean_30yr)]), 
                      ~get_density_curve(.x, def_subset))

time_range = range(slinky_data$year)

# --- 4. INDIVIDUAL PLOT DEFINITIONS ---

# Plot A: The 3-Panel Time Series
plot_a_data = rolling_data |>
  filter(!is.na(rolling_mean_30yr)) |>
  mutate(name_clean = factor(case_when(
    variable == "pr" ~ "Total JAS Precipitation (mm)",
    variable == "pet" ~ "Total JAS Reference ET (mm)",
    variable == "def" ~ "Total JAS Water Deficit (mm)"
  ), levels = c("Total JAS Precipitation (mm)", "Total JAS Reference ET (mm)", "Total JAS Water Deficit (mm)")))

target_breaks = c(1990, 2000, 2010, 2020)
target_labels = c("1961–1990", "1971–2000", "1981–2010", "1991–2020")

p1 = ggplot(plot_a_data, aes(x = year, y = rolling_mean_30yr)) +
  geom_point(alpha = 1, size = 1) +
  geom_smooth(method = "lm", color = "grey60", se = F, linewidth = 0.8) +
  geom_smooth(method = "loess", color = "blue", se = F, span = 0.5, linewidth = 0.8) +
  scale_x_continuous(limits = time_range, breaks = target_breaks, labels = target_labels) +
  facet_wrap(~name_clean, ncol = 1, scales = "free_y") +
  theme_bw() +
  labs(x = "30-Year Normal Period", y = "Parameter Value") +
  theme(strip.background = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 10),
        axis.text.x = element_text(angle = 35, hjust = 1))

# Plot B: Slinky PDF
climate_normal_ends = c(1990, 2000, 2010, 2020)
climate_normal_labels = paste0(climate_normal_ends - 29, "–", climate_normal_ends)

p2 = ggplot(slinky_data, aes(x = x, y = y, group = year, color = year)) +
  geom_line(linewidth = 0.7, alpha = 0.7) +
  # Add the vertical D3 reference line
  geom_vline(xintercept = d3_threshold, linetype = "dashed", color = "black", linewidth = 0.8) +
  # Labels for the threshold shift
  annotate("text", 
           x = -615, 
           y = 0.0125, 
           # Use paste to combine plotmath functions
           # \n does not work with parse = TRUE, so we use atop() for multi-line
           label = paste0("atop(bold(underline('", round(d3_threshold), "mm Deficit')),",
                          "atop('D3 in 1961–1990', 'D1 in 1991–2020'))"), 
           parse = TRUE, 
           hjust = 0, 
           size = 4.5, 
           lineheight = 0.9, 
           color = "black")+
  scale_color_gradientn(
    colors = turbo(10), name = "30-Year Normal Period",
    breaks = climate_normal_ends, labels = climate_normal_labels,
    limits = time_range,
    guide = guide_colorbar(title.position = "top", title.hjust = 0.5, barwidth = 20)
  ) +
  theme_bw() +
  labs(x = "JAS Potential Water Deficit (mm)", y = "Probability Density (PDF)") +
  # Note: Adjusting annotation X positions for seasonal magnitudes
  annotate("segment", x = -400, xend = -500, y = 0, yend = 0, arrow = arrow(length = unit(0.3, "cm"))) +
  annotate("text", x = -450, y = 0.001, label = "Delta*t", parse = TRUE, size = 6) +
  annotate("text", x = min(slinky_data$x)+20, y = 0.002, label = "Drier", fontface = "bold", color = "darkred", size = 5) +
  annotate("text", x = max(slinky_data$x)-20, y = 0.002, label = "Wetter", fontface = "bold", color = "darkblue", size = 5) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

# Map Inset
usa = st_as_sf(maps::map("state", fill = TRUE, plot = FALSE))
idaho = usa |> filter(ID == "idaho")
map_inset = ggplot() +
  geom_sf(data = idaho, fill = "white", color = "black", linewidth = 0.5) +
  geom_sf(data = ROI_raw, fill = "red", color = NA) +
  theme_void()

# --- 5. FINAL ASSEMBLY (COWPLOT) ---

title_gg <- ggdraw() + 
  draw_label("Summertime Nonstationary Drought Trends: Blaine County, ID", fontface = 'bold', size = 18) +
  theme(plot.background = element_rect(fill = "white", color = NA))

subtitle_gg <- ggdraw() + 
  draw_label(paste0('Moving 30-year Climatology: ', time_range[1]-29, '–', time_range[1], ' to ', time_range[2]-29, '–', time_range[2]), 
             fontface = 'italic', size = 12) +
  theme(plot.background = element_rect(fill = "white", color = NA))

# Align Plot B top with Plot A facet contents
plot_row <- plot_grid(
  p1 + theme(plot.background = element_rect(fill = "white", color = NA)), 
  p2 + theme(plot.background = element_rect(fill = "white", color = NA),
             plot.margin = margin(t = 30, r = 10, b = 5, l = 5)), 
  ncol = 2, rel_widths = c(1, 1.8), align = 'v', axis = 'lr'
)

caption_text <- "This analysis extracts 67 years of July–September (JAS) climate data to calculate 30-year rolling precipitation, reference ET, and potential water deficit (P-ETo), fitting\nthe deficit to Generalized Logistic distributions to visualize shifting hydroclimatic conditions in Blaine County. Data Source: TerraClimate (Abatzoglou et al., 2018)"

caption_gg <- ggdraw() + 
  draw_label(caption_text, size = 9, hjust = 0, x = 0.05) +
  theme(plot.background = element_rect(fill = "white", color = NA))

final_plot <- plot_grid(
  title_gg, subtitle_gg, plot_row, caption_gg, 
  ncol = 1, rel_heights = c(0.08, 0.04, 1, 0.12)
) + theme(plot.background = element_rect(fill = "white", color = NA))

# Place Map Inset
final_plot_attributed <- ggdraw(final_plot) +
  draw_plot(map_inset, x = 0.7, y = 0.52, width = 0.4, height = 0.3)

# Save Final Product
ggsave('~/nasem_examples/figs/nasem_example_1_hoylman.png', plot = final_plot_attributed, 
       width = 10, height = 7, bg = "white")

print(final_plot_attributed)

