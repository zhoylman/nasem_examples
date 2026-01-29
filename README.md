
# NASEM Nonstationary Drought Figure Examples

Two example R scripts produce NASEM-style figures that visualize how **moving 30-year climatologies** change the interpretation of summertime drought thresholds.  
Both scripts use TerraClimate via Google Earth Engine (rgee), fit moving-window Generalized Logistic distributions to the JAS (Jul–Aug–Sep) water balance (P − PET), and produce a two-panel figure:
- **A:** 3-panel 30-yr rolling means (pr, PET, P−PET)  
- **B:** “Slinky” of fitted PDFs colored by the 30-yr normal end year with a historical drought threshold marked.

## Repository contents
- `nasem_example_1_blaine.R` — Blaine County, Idaho example (D3 → D1 threshold illustration).  
- `nasem_example_2_palmbeach.R` — Palm Beach County, Florida example (D1 → D3 threshold illustration).  
- `output/` (recommended) — save figures and CSV outputs here.  
- External helper functions are sourced from `mco-drought-indicators` in the scripts.

## Quick summary
- **Input:** TerraClimate (`IDAHO_EPSCOR/TERRACLIMATE`) via GEE  
- **Bands:** `pr` (precipitation) and `pet` (reference ET)  
- **Period:** 1958–2024 (configurable in the scripts)  
- **Window:** trailing 30-year windows for rolling means and PDF fits  
- **Distribution:** Generalized Logistic via `lmomco`  
- **Output:** PNG figures saved with `ggsave()` (paths configurable)

## Requirements

### System
- R (>= 4.0 recommended)
- Conda with a Python environment accessible to `reticulate`

### R packages
```r
install.packages(c(
  "reticulate","sf","tidyverse","slider","purrr","lmomco","viridis","cowplot","maps"
))
# rgee:
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

## Troubleshooting
- **rgee auth issues:** run `rgee::ee_Auth()` and `rgee::ee_Initialize()` per rgee docs.  
- **Slow reduceRegion:** increase `scale` or reduce geometry resolution.  
- **Missing 30-year windows:** early/late years are skipped by `.complete = TRUE` in `slider::slide_index_dbl()`.

## Suggestions (future)
- Add `run_all.R` to drive both examples with a single config file.  
- Vendor county GeoJSON into `data/` and save intermediate CSVs to `output/`.  
- Add a `config.R` or `config.yml` to centralize user/project/path settings.

## License & attribution
- TerraClimate data: Abatzoglou et al., 2018 — cite appropriately.  
- Scripts: add a LICENSE (e.g., MIT) if you want to open-source; otherwise keep internal.

---

Maintainer: Zachary H. Hoylman  
If you want this README shortened, expanded to include diagrams, or converted to multiple docs (`CONTRIBUTING.md`, `run_all.R`), tell me and I’ll update it.
