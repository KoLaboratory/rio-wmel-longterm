# -----------------------------------------------------------------------------
# Temperature-dependent mosquito life-history and transmission functions
# -----------------------------------------------------------------------------
#
# Temperature-dependent rates for Aedes aegypti, taken from the Mordecai et al.
# thermal-biology papers:
#   - Mordecai et al. 2017, PLoS Negl Trop Dis, doi: 10.1371/journal.pntd.0005568
#   - Tesla et al. 2018 / Mordecai et al. 2019, doi: 10.1371/journal.pntd.0006451
#
# Each rate is a Briere (asymmetric) or quadratic (symmetric) function of
# temperature in degrees Celsius, with hard-coded parameters for Ae. aegypti.
# All functions return *daily* rates and are truncated to zero outside their
# biologically plausible thermal limits.
#
# These functions are sourced by R/seirsei_model.R and by the simulation
# scripts in analysis/. They have no side effects.
# -----------------------------------------------------------------------------

#' Briere function (asymmetric temperature response)
#'
#' @param x Numeric vector of temperatures (deg C).
#' @param a Scaling constant.
#' @param Tmin,Tmax Lower and upper thermal limits (deg C).
#' @return Numeric vector of rates; zero outside `[Tmin, Tmax]`.
Briere <- function(x, a, Tmin, Tmax) {
  result <- rep(0, length(x))
  valid <- (x >= Tmin) & (x <= Tmax) & is.finite(x)
  result[valid] <- a * x[valid] * (x[valid] - Tmin) * sqrt(Tmax - x[valid])
  result
}

#' Quadratic function (symmetric temperature response)
#'
#' @param x Numeric vector of temperatures (deg C).
#' @param a Scaling constant (negative for a downward-opening parabola).
#' @param Tmin,Tmax Lower and upper thermal limits (deg C).
#' @return Numeric vector of rates; zero outside `[Tmin, Tmax]`.
Quad <- function(x, a, Tmin, Tmax) {
  result <- rep(0, length(x))
  valid <- (x >= Tmin) & (x <= Tmax) & is.finite(x)
  result[valid] <- a * (x[valid] - Tmin) * (x[valid] - Tmax)
  result
}

# ---- Individual rate functions ---------------------------------------------

#' Eggs laid per female per day.
EggsLaid <- function(tmp) Briere(tmp, a = 8.56e-3, Tmin = 14.58, Tmax = 34.61)

#' Mosquito egg-to-adult survival probability.
SurvProb <- function(tmp) Quad(tmp, a = -5.99e-3, Tmin = 13.56, Tmax = 38.29)

#' Mosquito egg-to-adult development rate per day.
DevlRate <- function(tmp) Briere(tmp, a = 7.86e-5, Tmin = 11.36, Tmax = 39.17)

#' Mosquito mortality rate per day (inverse of lifespan).
#'
#' The rate is capped at 1 so the minimum mosquito lifespan is at least one day;
#' otherwise it can diverge to infinity near the thermal limits.
MortRate <- function(tmp) {
  m <- 1 / Quad(tmp, a = -1.48e-1, Tmin = 9.16, Tmax = 37.73)
  m[m > 1] <- 1
  m
}

#' Biting rate per day.
BiteRate <- function(tmp) Briere(tmp, a = 2.02e-4, Tmin = 13.35, Tmax = 40.08)

#' Probability of mosquito infection (given a bite on an infectious host).
ProbInfc <- function(tmp) Briere(tmp, a = 4.91e-4, Tmin = 12.22, Tmax = 37.46)

#' Probability of mosquito infectiousness.
ProbInfs <- function(tmp) Briere(tmp, a = 8.49e-4, Tmin = 17.05, Tmax = 35.83)

#' Virus extrinsic incubation rate per day.
IncuRate <- function(tmp) Briere(tmp, a = 6.65e-5, Tmin = 10.68, Tmax = 45.90)

# ---- Carrying capacity ------------------------------------------------------

#' Temperature-dependent mosquito carrying capacity (internal form).
#'
#' Peaks at the reference temperature `Tref` and scales with `N` (typically the
#' host population size times the mosquito-to-host ratio).
#'
#' @param x Numeric vector of temperatures (deg C).
#' @param Tref Reference temperature at which capacity peaks (deg C).
#' @param N Maximum capacity scale.
#' @param E_a Activation-energy-like parameter controlling the width of the peak.
CarrCapFun <- function(x, Tref, N, E_a = 0.05) {
  k <- 8.617e-5  # Boltzmann constant (eV/K)
  term <- EggsLaid(Tref) * SurvProb(Tref) * DevlRate(Tref) / MortRate(Tref)
  (term - MortRate(Tref)) / term *
    N * exp((-E_a * (x - Tref)^2) / (k * (x + 273.15) * (Tref + 273.15)))
}

#' Temperature-dependent mosquito carrying capacity (peaks at 29 deg C).
#'
#' @param tmp Numeric vector of temperatures (deg C).
#' @param N_max Maximum capacity scale.
#' @param E_a Activation-energy-like parameter (default 0.05).
CarrCap <- function(tmp, N_max, E_a = 0.05) {
  CarrCapFun(x = tmp, Tref = 29, N = N_max, E_a = E_a)
}
