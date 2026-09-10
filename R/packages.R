# -----------------------------------------------------------------------------
# Package loading
# -----------------------------------------------------------------------------
#
# Sourced at the top of every analysis script so the dependency set is declared
# in exactly one place. INLA is not on CRAN; see README.md / DESCRIPTION for its
# install command. Package versions are pinned via renv (renv.lock) once you run
# renv::snapshot() -- see README.md.
# -----------------------------------------------------------------------------

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")

# Empirical spatiotemporal analysis (analysis/01-0x_*.R)
pacman::p_load(
  tidyverse,
  sf,
  terra,
  raster,
  exactextractr,
  fmesher,
  patchwork,
  slider,
  scales,
  geobr,
  ggspatial,
  ggnewscale,
  RColorBrewer,
  zoo,
  nngeo           # st_remove_holes()
)

# Mechanistic simulations (analysis/03-04_*.R)
pacman::p_load(
  deSolve,
  colorspace,
  grid
)

# INLA is loaded explicitly (not via pacman) because it lives in a separate repo.
# install.packages("INLA",
#   repos = c(getOption("repos"), INLA = "https://inla.r-inla-download.org/R/stable"),
#   dep = TRUE)
suppressPackageStartupMessages(library(INLA))

# Optional, only for Supplementary Figure 3 (census tract populations):
#   pacman::p_load(censobr)
