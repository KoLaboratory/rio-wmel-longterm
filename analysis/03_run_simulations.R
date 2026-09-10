# =============================================================================
# 03_run_simulations.R
# -----------------------------------------------------------------------------
# Run temperature-dependent SEIR-SEI simulations for Figure 5 and
# Supplementary Figures S5-S8. Five scenarios are produced over a grid of
# mosquito-to-host ratios (M0 = 1..10) and wMel-driven reductions
# (frac = 0.1..1.0, where frac = 1 means no wMel):
#
#   1. Early introduction, near-fully susceptible population (S_h = 0.999)
#   2. Early introduction, half-susceptible population       (S_h = 0.5)
#   3. Late introduction, near endemic equilibrium           (two-stage burn-in)
#   4. Early introduction with a transient temperature-amplitude shock (year 5)
#   5. Late introduction with a transient temperature-amplitude shock  (year 25)
#
# Each result is cached under data/interim/.
# =============================================================================

# The simulations do not need INLA or the spatial stack, only the ODE solver and
# a little tidyverse, so they load their own light dependency set rather than
# R/packages.R (which pulls in INLA/sf for the empirical analysis).
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(dplyr, readr, tidyr, deSolve)

source(file.path("R", "temperature_functions.R"))
source(file.path("R", "seirsei_model.R"))      # RunSEIRSEI1Sero, run_scenarios, ...

FORCE_RERUN <- FALSE
if (!dir.exists("data/interim")) dir.create("data/interim", recursive = TRUE)

#' Run a scenario only if its cache is missing (or FORCE_RERUN is TRUE).
cache_run <- function(path, fun) {
  if (!FORCE_RERUN && file.exists(path)) {
    message("Cache exists, skipping: ", path)
    return(invisible(readRDS(path)))
  }
  message("Running: ", path)
  result <- fun()
  saveRDS(result, path)
  invisible(result)
}

#' Temperature series with a transient doubling of the annual amplitude during
#' one "shock" year, ramped in and out over 30 days with a cosine taper.
#' Scalar-in, scalar-out (called per time step by the ODE solver).
#'
#' @param shock_year The 1-indexed simulation year in which amplitude doubles.
make_amplitude_temp <- function(shock_year) {
  function(x, mean_temp = 26, ampl_temp = 3) {
    year <- floor(x / 365) + 1
    if (year == shock_year) {
      day_in_year <- x %% 365
      if (day_in_year < 30) {
        scale_factor <- 1 + (1 - cos(pi * day_in_year / 30)) / 2
      } else if (day_in_year > 335) {
        scale_factor <- 1 + (1 + cos(pi * (day_in_year - 335) / 30)) / 2
      } else {
        scale_factor <- 2
      }
      current_ampl <- ampl_temp * scale_factor
    } else {
      current_ampl <- ampl_temp
    }
    tmp <- current_ampl * sin(2 * pi * x / 365) + mean_temp
    names(tmp) <- NULL
    tmp
  }
}

# ---- Shared scenario grid ---------------------------------------------------
results_grid <- expand.grid(
  M0   = seq(1, 10, by = 1),
  frac = seq(0.1, 1, by = 0.1),
  KEEP.OUT.ATTRS = FALSE
)

# Base parameters. `M` is set per run inside the run_* helpers, so it is left
# out here; every other parameter is shared across scenarios except S_h.
base_pars <- c(
  N_h = 1e6,             # host population size
  S_h = 0.999,           # starting proportion susceptible
  I_h = 0.001,           # starting proportion infectious
  mu_h = 1 / (77 * 365), # host mortality rate
  r_h = 1 / (77 * 365),  # host birth rate
  gamma_h = 1 / 5,       # recovery rate (5-day infectious period)
  delta_h = 1 / 5.9,     # intrinsic incubation rate
  trickle = 1e-5         # constant infectious-mosquito influx (proportion of N_h)
)

# ---- Scenario 1: early introduction, near-fully susceptible -----------------
ts_early <- seq_len(365 * 15 + 1)
cache_run(
  "data/interim/seirsei_fig5-early intro suscept.rds",
  function() run_scenarios(results_grid, base_pars, ts_early, temp_fun = GetTemp)
)

# ---- Scenario 2: early introduction, half susceptible -----------------------
pars_half <- base_pars
pars_half["S_h"] <- 0.5
cache_run(
  "data/interim/seirsei_fig5-early intro half suscept.rds",
  function() run_scenarios(results_grid, pars_half, ts_early, temp_fun = GetTemp)
)

# ---- Scenario 3: late introduction, near endemic equilibrium ----------------
n_days_total <- 365 * 41
t_break      <- 365 * 20
times1 <- 1:t_break
times2 <- (t_break + 1):(n_days_total + 1)
cache_run(
  "data/interim/seirsei_fig5-late intro suscept.rds",
  function() run_scenarios_late(results_grid, base_pars, times1, times2, temp_fun = GetTemp)
)

# ---- Scenario 4: early introduction, temperature-amplitude shock (year 5) ----
cache_run(
  "data/interim/seirsei_fig5-amplitude change.rds",
  function() run_scenarios(results_grid, base_pars, ts_early,
                           temp_fun = make_amplitude_temp(shock_year = 5))
)

# ---- Scenario 5: late introduction, temperature-amplitude shock (year 25) ----
cache_run(
  "data/interim/seirsei_fig5-amplitude change late.rds",
  function() run_scenarios_late(results_grid, base_pars, times1, times2,
                                temp_fun = make_amplitude_temp(shock_year = 25))
)

message("All simulation scenarios ready in data/interim/.")
