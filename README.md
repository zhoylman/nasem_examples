# NASEM Nonstationary Drought Figure Examples

This repository contains **three R scripts** developed for the National Academies of Sciences, Engineering, and Medicine (NASEM) study  
*The Future of Drought in the United States*.  
Together, they demonstrate how **nonstationarity and reference period design** influence drought classification, severity, and interpretation.

Two scripts use **observed climate data** (TerraClimate via Google Earth Engine) to illustrate shifting drought thresholds under moving climatologies.  
A third script uses **synthetic data** to isolate and visualize how different forms of nonstationarity affect drought frequency under a percentile-based framework.

---

## Repository contents

### Empirical examples (moving 30-year climatologies)
- **`nasem_example_1_blaine.R`**  
  *Blaine County, Idaho*  
  Demonstrates how a fixed historical drought threshold (e.g., D3) maps to different percentile classes as the 30-year climatology shifts toward drier conditions.

- **`nasem_example_2_palmbeach.R`**  
  *Palm Beach County, Florida*  
  Demonstrates the opposite case, where increasing moisture supply causes a fixed deficit to transition toward wetter percentile classifications (e.g., D1 → D3).

These scripts produce a two-panel figure:
- **Panel A:** Three stacked time series showing moving 30-year means of  
  precipitation (P), reference evapotranspiration (ET₀), and water balance (P − ET₀).  
- **Panel B:** A “slinky” of fitted probability density functions (PDFs) for P − ET₀, colored by the end year of the 30-year normal period, with a historical drought threshold indicated.

---

### Synthetic example (expanding climatology; Figure 6.4)
- **`nie_etal_example_expanded.R`**  

  A synthetic experiment adapted from **Nie et al. (2025)** that isolates how different forms of nonstationarity influence drought classification under a percentile-based framework.

  This script generates daily “drought indicator” time series over 100 years and evaluates drought categories using an **expanding reference climatology**:
  - Minimum baseline: 30 years  
  - After year 30, each year’s drought classes are computed using the **longest available period of record** (e.g., year 40 uses 40 years; year 100 uses 100 years).

  Four idealized scenarios are illustrated:
  - **A:** Stationary climate variability  
  - **B:** Linear drying trend  
  - **C:** Decadal oscillation (low-frequency variability, no trend)  
  - **D:** Increasing variance through time  

  Output figures show:
  - **Left panels:** Daily drought index time series with a seasonal reference signal overlaid  
  - **Right panels:** Annual days in drought by category (D0–D4), plus cumulative percent of time spent in drought over the final 10 years

Nie, W., Kumar, S. V., & Zhao, L. (2025). Anthropogenic influences on the water cycle amplify uncertainty in drought assessments. One Earth, 8(2).

---

## Quick summary

**Empirical scripts**
- **Data:** TerraClimate (`IDAHO_EPSCOR/TERRACLIMATE`) via Google Earth Engine  
- **Variables:**  
  - `pr` — precipitation  
  - `pet` — reference evapotranspiration  
- **Season:** July–September (JAS)  
- **Period:** 1958–2024 (configurable)  
- **Method:**  
  - Trailing 30-year windows  
  - Generalized Logistic distributions fit via L-moments (`lmomco`)  
- **Purpose:** Show how shifting climatologies alter drought percentiles and categories

**Synthetic script**
- **Data:** Fully simulated daily time series  
- **Classification:** Percentile-based thresholds consistent with the U.S. Drought Monitor  
- **Reference period:** Expanding climatology with a 30-year minimum  
- **Purpose:** Isolate effects of trends, oscillations, and variance changes on drought frequency

---

## Requirements

### System
- R (≥ 4.0 recommended)
- Conda / Python environment accessible via `reticulate`

### R packages
```r
install.packages(c(
  "reticulate", "sf", "tidyverse", "slider", "purrr",
  "lmomco", "viridis", "cowplot", "maps", "scales", "shadowtext"
))
remotes::install_github("r-spatial/rgee")
```

### Google Earth Engine
- Earth Engine account and `rgee` authenticated.
- Update `ee_Initialize()` in each script with your `user` and `project`.

## Before you run (edit these in each script)
1. **Conda env name**: `reticulate::use_condaenv("gee_env", required = TRUE)` — change `"gee_env"` if needed.  
2. **GEE user/project**: set `ee_Initialize(user = 'YOUR_EMAIL', drive = TRUE, project = 'YOUR_GEE_PROJECT')`.  
3. **rgee session file path**: scripts write a session file at `~/.config/earthengine/rgee_sessioninfo.txt` — change if your home path differs.  
4. **Output paths**: scripts currently save to `~/nasem_examples/figs/` — change to `output/` for repo relative paths.  
5. **County GeoJSON**: scripts read a remote counties GeoJSON. For reproducibility, add the GeoJSON to `data/` and point scripts to it.  
6. **Units**: scripts multiply `pet` by `0.1` after extraction — keep consistent if you change datasets.

## Run (example)
1. Start R / RStudio and set working directory to repo root.  
2. Edit top-of-script configuration (user/project, env, paths).  
3. Run:
```r
source("nasem_example_1_blaine.R")
# and/or
source("nasem_example_2_palmbeach.R")
```
4. Figures and (optional) CSVs will be saved to the configured output path.

## Output & interpretation
- The vertical dashed line marks a fixed historical drought threshold (e.g., D1 or D3).  
- The slinky PDFs show how the distribution of JAS water-balance shifts across moving 30-yr normals.  
- A single physical deficit may map to different percentiles/categories across windows — illustrating nonstationarity effects.

## License & attribution
- TerraClimate data: Abatzoglou et al., 2018 — cite appropriately.  
- Scripts: CC-BY 4.0

---

Maintainer: Zachary H. Hoylman  
