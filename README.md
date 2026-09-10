# Code for "Long-term effectiveness of Wolbachia introgression in preventing dengue in Rio de Janeiro, Brazil: The potential role of immune dynamics"

This manuscript evaluates the long-term effectiveness of *w*Mel *Wolbachia* releases in Rio de Janeiro, Brazil, over eight years of surveillance using Bayesian spatiotemporal models to relate *w*Mel introgression to dengue incidence across the release area and analyze confounding by pre-intervention spatial risk, and temperature-dependent SEIR-SEI models to explore how population immunity dynamics shape observed intervention effectiveness.

## Repository layout

```
rio-wmel-longterm/
├── R/                            # reusable functions, sourced by scripts
│   ├── packages.R                #   package loading (empirical stack + INLA)
│   ├── theme.R                   #   plotting helpers and project paths
│   ├── inla_helpers.R            #   mesh/SPDE/stack builders, RI extraction, g-computation
│   ├── temperature_functions.R   #   Ae. aegypti thermal-response functions
│   └── seirsei_model.R           #   SEIR-SEI ODE model and simulation runners
├── analysis/                     # pipeline (run in order)
│   ├── 01_prepare_data.R         #   load/wrangle inputs, build grid -> prepared.rds
│   ├── 02_empirical_analysis.R   #   INLA models + Figures 1-4 + supplements/tables
│   ├── 03_run_simulations.R      #   run + cache the 5 SEIR-SEI scenarios
│   └── 04_simulation_figures.R   #   Figure 5 + Supplementary S5-S8
├── data/                         # (not committed, data are restricted)
│   ├── raw/                      #   raw inputs
│   ├── interim/                  #   cached intermediates
├── outputs/                      # (not committed)
├── run_all.R                     # reproduce everything end to end
├── DESCRIPTION                   # dependency list
└── LICENSE
```

## Requirements

- **R ≥ 4.5** (the manuscript used R 4.5.1).
- **INLA**, which is not on CRAN:
  ```r
  install.packages(
    "INLA",
    repos = c(getOption("repos"), INLA = "https://inla.r-inla-download.org/R/stable"),
    dep = TRUE
  )
  ```
- The remaining CRAN packages are loaded via
  [`pacman`](https://cran.r-project.org/package=pacman) at the top of each script
  and are listed in `DESCRIPTION`. Spatial packages (`sf`, `terra`, `raster`,
  `exactextractr`) require system libraries GDAL, GEOS and PROJ.

### Reproducible environment (recommended)

Pin exact package versions with [`renv`](https://rstudio.github.io/renv/)
(run in console after commit)

```r
install.packages("renv")
renv::restore()
```

Commit `renv.lock` once generated; the `renv/library/` folder is git-ignored.

## Reproducing the analysis

From the repository root, with the data in place:

```r
# everything, in order:
Rscript run_all.R

# or run individual steps (each can be run on its own):
Rscript analysis/01_prepare_data.R
Rscript analysis/02_empirical_analysis.R
Rscript analysis/03_run_simulations.R
Rscript analysis/04_simulation_figures.R
```

Notes:

- The **empirical** (01–02) and **simulation** (03–04) halves are independent.
  The simulation scripts do not require INLA or the spatial stack.
- `01_prepare_data.R` writes `data/interim/prepared.rds`; `02` loads it (or reuses
  the objects if run in the same session via `run_all.R`).
- `03_run_simulations.R` **caches** each scenario under `data/interim/`. The ODE
  integrations are the slow step; reruns skip completed scenarios. Delete a
  cache file or set `FORCE_RERUN <- TRUE` to recompute.
- Figures are printed to the active device. `ggsave()` calls are left commented
  (they already point at `outputs/`); uncomment them to write files.

## Where each figure/table is produced

| Manuscript item | Script |
|---|---|
| Figure 1 (descriptive maps, case time series, introgression) | `analysis/02_empirical_analysis.R` |
| Figure 2 (primary cumulative & year-specific models) | `analysis/02_empirical_analysis.R` |
| Figure 3 (pre-intervention risk & future introgression) | `analysis/02_empirical_analysis.R` |
| Figure 4 (effect under alternative spatial-risk adjustments) | `analysis/02_empirical_analysis.R` |
| Figure 5 (simulation results) | `analysis/04_simulation_figures.R` |
| Table S2 (model comparison / sensitivity) | `analysis/02_empirical_analysis.R` |
| Table S3 (adjustment for past cumulative incidence) | `analysis/02_empirical_analysis.R` |
| Figures S2–S4 (introgression tileplot, municipal comparison, ovitrap models) | `analysis/02_empirical_analysis.R` |
| Figures S5–S8 (simulation supplements) | `analysis/04_simulation_figures.R` |

## License

Code is released under the MIT License (see `LICENSE`).