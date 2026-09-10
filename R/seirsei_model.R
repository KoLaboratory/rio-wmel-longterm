# -----------------------------------------------------------------------------
# Temperature-dependent SEIR-SEI dengue transmission model
# -----------------------------------------------------------------------------
#
# A deterministic, single-serotype model coupling human SEIR and mosquito SEI
# compartments, with mosquito life-history and viral development rates driven by
# temperature (see R/temperature_functions.R). wMel introduction is represented
# by reducing the mosquito-to-host ratio `M`, which lowers the effective
# carrying capacity of transmission-competent vectors.
#
# Requires: deSolve (for ode() in the run helpers) and the temperature functions.
# Source R/temperature_functions.R before this file.
# -----------------------------------------------------------------------------

if (!exists("Briere")) {
  source(file.path("R", "temperature_functions.R"))
}

# ---- Model definition -------------------------------------------------------

#' Right-hand side of the temperature-dependent SEIR-SEI ODE system.
#'
#' Compartments (state vector `x`, named):
#'   Sh, Eh, Ih, Rh  susceptible / exposed / infectious / recovered humans
#'   Rhc             cumulative recovered humans (a proxy for cumulative incidence)
#'   Sv, Ev, Iv      susceptible / exposed / infectious mosquitoes
#'
#' @param tms Numeric time point(s) at which to evaluate the derivatives.
#' @param x Named numeric vector of current state values.
#' @param pars Named numeric vector of parameters. Must contain:
#'   N_h (host population size), M (mosquito-to-host ratio), mu_h (host mortality
#'   rate), r_h (host birth rate), gamma_h (recovery rate), delta_h (intrinsic
#'   incubation rate), trickle (constant infectious-mosquito influx as a
#'   proportion of N_h).
#' @param input A function returning temperature (deg C) for a given time.
#' @return A list suitable for deSolve::ode(): the derivative vector plus the
#'   human-to-mosquito transmission coefficient and the temperature.
RunSEIRSEI1Sero <- function(tms, x, pars, input) {
  with(as.list(c(pars)), {
    dxdt <- vector("numeric", length = 8)
    names(dxdt) <- paste0("d", c("Sh", "Eh", "Ih", "Rh", "Rhc", "Sv", "Ev", "Iv"))

    # Temperature at this time step and the temperature-dependent rates.
    tmp <- input(tms)
    biterate_a <- BiteRate(tmp)
    probinfs_a <- ProbInfs(tmp)
    probinfc_a <- ProbInfc(tmp)
    mortrate_a <- MortRate(tmp)
    incurate_a <- IncuRate(tmp)
    eggslaid_a <- EggsLaid(tmp)
    survprob_a <- SurvProb(tmp)
    devlrate_a <- DevlRate(tmp)
    carrcapa_a <- CarrCap(tmp = tmp, N_max = N_h * M)

    comp_names <- names(x)
    vs <- grep("v", comp_names)   # vector compartments
    N_v <- sum(x[vs])

    # Force of infection: mosquito -> human.
    hivi <- biterate_a * probinfs_a * x["Iv"] / N_h

    dxdt["dSh"] <- r_h * N_h - (hivi + mu_h) * x["Sh"]
    dxdt["dEh"] <- hivi * x["Sh"] - (delta_h + mu_h) * x["Eh"]
    dxdt["dIh"] <- delta_h * x["Eh"] - (gamma_h + mu_h) * x["Ih"]
    dxdt["dRh"] <- gamma_h * x["Ih"] - mu_h * x["Rh"]
    dxdt["dRhc"] <- gamma_h * x["Ih"]

    # Mosquito recruitment (density-dependent) and force of infection human -> mosquito.
    birth_rate <- eggslaid_a * survprob_a * devlrate_a *
      N_v * (1 - N_v / carrcapa_a) / mortrate_a
    vihi <- biterate_a * probinfc_a * x["Ih"] / N_h

    # Small constant infectious-mosquito influx keeps populations from going
    # extinct (and the solver stable) when temperatures get low.
    in_mosq <- trickle * N_h

    dxdt["dSv"] <- birth_rate - (vihi + mortrate_a) * x["Sv"]
    dxdt["dEv"] <- vihi * x["Sv"] - (incurate_a + mortrate_a) * x["Ev"]
    dxdt["dIv"] <- incurate_a * x["Ev"] - mortrate_a * x["Iv"] + in_mosq

    list(dxdt, beta_h = biterate_a * probinfs_a, temp = tmp)
  })
}

# ---- Helpers ----------------------------------------------------------------

#' Build an initial state vector from a parameter vector.
#'
#' Starts the epidemic with a seed of infectious humans (`I_h`), a susceptible
#' fraction (`S_h`), the remainder recovered, and a fully-susceptible mosquito
#' population sized at `N_h * M`.
#'
#' @param pars_vec Named parameter vector containing N_h, M, S_h, I_h.
#' @return Named numeric state vector (Sh, Eh, Ih, Rh, Rhc, Sv, Ev, Iv).
make_xstart <- function(pars_vec) {
  x0 <- setNames(numeric(8), c("Sh", "Eh", "Ih", "Rh", "Rhc", "Sv", "Ev", "Iv"))
  x0["Sv"] <- pars_vec["N_h"] * pars_vec["M"]
  x0["Sh"] <- pars_vec["S_h"] * pars_vec["N_h"]
  x0["Ih"] <- pars_vec["I_h"] * pars_vec["N_h"]
  x0["Rh"] <- pars_vec["N_h"] - x0["Sh"] - x0["Ih"]
  x0["Eh"] <- 0
  x0["Rhc"] <- 0
  x0["Ev"] <- 0
  x0["Iv"] <- 0
  x0
}

#' Synthetic seasonal temperature series (sine curve).
#'
#' The model queries a temperature function at arbitrary time points; this is a
#' smooth annual cycle. To use real temperatures instead, interpolate a daily
#' series with `approxfun()` and pass that as `input` to the solver.
#'
#' @param x Numeric time point(s) (days).
#' @param mean_temp Mean temperature (deg C).
#' @param ampl_temp Amplitude of the annual cycle (deg C).
GetTemp <- function(x, mean_temp = 26, ampl_temp = 3) {
  tmp <- ampl_temp * sin(2 * pi * x / 365) + mean_temp
  names(tmp) <- NULL
  tmp
}

#' Run one single-stage SEIR-SEI simulation over a wMel scenario grid.
#'
#' For each (M0, frac) combination in `grid`, sets the mosquito-to-host ratio to
#' `M0 * frac` (frac = 1 means no wMel; frac < 1 represents introgression that
#' reduces transmission-competent vectors) and solves the ODE from a fresh
#' initial state over `ts`.
#'
#' @param grid A data.frame with columns M0 and frac (e.g. from expand.grid()).
#' @param pars Base parameter vector (M is overwritten per scenario).
#' @param ts Numeric vector of time points to solve over.
#' @param temp_fun Temperature function passed as `input` to the model.
#' @param rtol Relative tolerance for the solver.
#' @return A long data.frame of trajectories tagged with M0, frac and Mnew.
run_scenarios <- function(grid, pars, ts, temp_fun = GetTemp, rtol = 1e-5) {
  traj_list <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    M0_i <- grid$M0[i]
    frac_i <- grid$frac[i]
    Mnew_i <- M0_i * frac_i

    pars_i <- pars
    pars_i["M"] <- Mnew_i
    xstart_i <- make_xstart(pars_i)

    sim <- deSolve::ode(
      y = xstart_i, times = ts, func = RunSEIRSEI1Sero,
      parms = pars_i, input = temp_fun, rtol = rtol
    )
    a <- as.data.frame(sim)
    traj_list[[i]] <- cbind(
      M0 = M0_i, frac = frac_i, Mnew = Mnew_i,
      a[, c("time", "Sh", "Eh", "Ih", "Rh", "Rhc", "Sv", "Ev", "Iv")]
    )
  }
  dplyr::bind_rows(traj_list)
}

#' Run a two-stage "late introduction" SEIR-SEI simulation over a scenario grid.
#'
#' Stage 1 runs the epidemic to (near) endemic equilibrium at the baseline
#' mosquito-to-host ratio `M0`; wMel is then introduced at `t_break` by scaling
#' M to `M0 * frac` for stage 2. Stage-1 solutions are cached per unique M0 so
#' they are not recomputed for every `frac`.
#'
#' @param grid A data.frame with columns M0 and frac.
#' @param pars Base parameter vector (must NOT already contain M).
#' @param times1,times2 Time vectors for stage 1 (burn-in) and stage 2.
#' @param temp_fun Temperature function passed as `input`.
#' @param rtol Relative tolerance for the solver.
#' @return A long data.frame of full (stage1 + stage2) trajectories tagged with
#'   M0, frac and Mnew.
run_scenarios_late <- function(grid, pars, times1, times2, temp_fun = GetTemp, rtol = 1e-5) {
  state_cols <- c("Sh", "Eh", "Ih", "Rh", "Rhc", "Sv", "Ev", "Iv")
  cache_stage1 <- new.env(parent = emptyenv())
  traj_list <- vector("list", nrow(grid))

  for (i in seq_len(nrow(grid))) {
    M0_i <- grid$M0[i]
    frac_i <- grid$frac[i]
    Mnew_i <- M0_i * frac_i

    # Stage 1: burn-in to (near) equilibrium at baseline M0 (cached by M0).
    key <- as.character(M0_i)
    if (exists(key, envir = cache_stage1)) {
      stage1 <- get(key, envir = cache_stage1)
    } else {
      pars1 <- pars
      pars1["M"] <- M0_i
      sim1 <- deSolve::ode(
        y = make_xstart(pars1), times = times1, func = RunSEIRSEI1Sero,
        parms = pars1, input = temp_fun, rtol = rtol
      )
      stage1 <- as.data.frame(sim1)
      assign(key, stage1, envir = cache_stage1)
    }

    # Stage 2: introduce wMel (scale M) starting from the end of stage 1.
    xstart2 <- as.numeric(stage1[nrow(stage1), state_cols])
    names(xstart2) <- state_cols
    pars2 <- pars
    pars2["M"] <- Mnew_i
    sim2 <- deSolve::ode(
      y = xstart2, times = times2, func = RunSEIRSEI1Sero,
      parms = pars2, input = temp_fun, rtol = rtol
    )
    stage2 <- as.data.frame(sim2)

    full <- rbind(stage1[, c("time", state_cols)], stage2[, c("time", state_cols)])
    traj_list[[i]] <- cbind(M0 = M0_i, frac = frac_i, Mnew = Mnew_i, full)
  }
  dplyr::bind_rows(traj_list)
}
