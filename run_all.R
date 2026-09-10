# =============================================================================
# run_all.R
# -----------------------------------------------------------------------------
# Reproduce the full analysis end to end. Run from the repository root:
#
#   Rscript run_all.R
#
# The empirical (01-02) and simulation (03-04) halves are independent.
# Simulation results are cached, so 03 is fast on reruns.
#
# Requires the raw data files in data/raw/ and the INLA
# package (see README.md) for the empirical steps.
# =============================================================================

# Make sure we are at the project root (the folder containing this file).
stopifnot(dir.exists("R"), dir.exists("analysis"))

message("== 01 prepare data =========================================")
source("analysis/01_prepare_data.R")

message("== 02 empirical analysis ===================================")
source("analysis/02_empirical_analysis.R")

message("== 03 run simulations ======================================")
source("analysis/03_run_simulations.R")

message("== 04 simulation figures ===================================")
source("analysis/04_simulation_figures.R")

message("== done ====================================================")
