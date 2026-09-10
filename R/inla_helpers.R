# -----------------------------------------------------------------------------
# Helpers for the INLA spatiotemporal negative-binomial models
# -----------------------------------------------------------------------------
#
# Every model in the empirical analysis shares the same spatial scaffolding:
# a 2D mesh over the cell centroids, an SPDE (Matern) spatial random effect,
# a projector matrix A, and an inla.stack. Only the data frame, the spatial
# range, and the model formula differ between models. These helpers remove that
# boilerplate so each script can focus on the formula and the inla() call, which
# are the parts that actually differ between models.
#
# Requires: INLA (and fmesher, which INLA depends on). Load INLA before use.
# -----------------------------------------------------------------------------

#' Build the spatial mesh, SPDE and stack for an INLA model.
#'
#' Reproduces the scaffolding used throughout the empirical analysis. The mesh
#' is built over the `x`/`y` centroid columns of `data`; the response is stacked
#' under the name `obs` together with an explicit `Intercept` column so models
#' can be written with `-1 + Intercept`.
#'
#' @param data A data frame containing the response, projected centroid columns
#'   `x` and `y`, and any covariates used in the formula.
#' @param response Name of the response column (default "cases").
#' @param spatial_range Range parameter controlling mesh resolution (map units,
#'   here metres). Default 22000, as used in the primary models.
#' @param cutoff Minimum distance between mesh nodes (default 200).
#' @param field_name Name of the SPDE index/effect (default "spatial.field").
#' @return A list with elements `mesh`, `spde`, `s.index`, `A`, `stack`.
build_spatial_stack <- function(data,
                                response = "cases",
                                spatial_range = 22000,
                                cutoff = 200,
                                field_name = "spatial.field") {
  loc <- as.matrix(data[, c("x", "y")])

  mesh <- fmesher::fm_mesh_2d_inla(
    loc = loc,
    max.edge = c(spatial_range / 5, spatial_range),
    cutoff = cutoff,
    offset = c(spatial_range / 5, diff(range(data$x)) / 10)
  )

  spde <- INLA::inla.spde2.matern(mesh = mesh)
  s.index <- INLA::inla.spde.make.index(name = field_name, n.spde = spde$n.spde)
  A.est <- INLA::inla.spde.make.A(mesh = mesh, loc = loc)

  stack <- INLA::inla.stack(
    data = list(obs = data[[response]]),
    A = list(A.est, 1),
    effects = list(
      spatial_field = data.frame(s.index, Intercept = 1),
      data[, names(data) != response]
    ),
    tag = "stdata"
  )

  list(mesh = mesh, spde = spde, s.index = s.index, A = A.est, stack = stack)
}

#' Fit a negative-binomial INLA model on a stack from build_spatial_stack().
#'
#' A thin wrapper around INLA::inla() with the control options used throughout
#' the analysis (INLA factor expansion; config/DIC/WAIC computed).
#'
#' @param formula A model formula referencing `obs`, `Intercept`, `spde`, etc.
#' @param stack An inla.stack (the `stack` element from build_spatial_stack()).
#' @param verbose Passed to inla() (default FALSE).
#' @param ... Additional arguments passed to INLA::inla().
#' @return The fitted inla object.
fit_nbinomial <- function(formula, stack, verbose = FALSE, ...) {
  INLA::inla(
    formula,
    data = INLA::inla.stack.data(stack),
    family = "nbinomial",
    control.predictor = list(A = INLA::inla.stack.A(stack), compute = TRUE),
    control.fixed = list(expand.factor.strategy = "inla"),
    control.compute = list(config = TRUE, dic = TRUE, waic = TRUE),
    verbose = verbose,
    ...
  )
}

#' Ratio of incidences for a +10-percentage-point change in introgression.
#'
#' Transforms a posterior marginal of a coefficient (on the log scale, for the
#' 0-1 introgression proportion) into the incidence ratio associated with a
#' 0.10 increase, and summarises it.
#'
#' @param marginal An INLA fixed-effect marginal, e.g.
#'   `fit$marginals.fixed$mintro`.
#' @param digits Optional rounding for the returned summary (NULL = no rounding).
#' @return A named numeric vector (mean, sd, quantiles) as from inla.zmarginal().
ri_per_10pp <- function(marginal, digits = NULL) {
  m <- INLA::inla.tmarginal(function(x) exp(0.1 * x), marginal)
  out <- unlist(INLA::inla.zmarginal(m, silent = TRUE))
  if (!is.null(digits)) out <- round(out, digits)
  out
}

#' Posterior effect curve of the incidence ratio across the introgression range.
#'
#' Draws posterior samples of a coefficient and returns the mean incidence ratio
#' and a 95% credible band over a grid of introgression values, for plotting the
#' RI vs introgression relationship (Figures 2, 3, 4).
#'
#' @param marginal An INLA fixed-effect marginal for the introgression term.
#' @param n_samples Number of posterior draws (default 1000).
#' @param mintro_grid Grid of introgression values (default seq(0, 1, len = 100)).
#' @param seed Optional RNG seed for reproducibility.
#' @return A data frame with columns mintro, effect, lower, upper.
ri_effect_curve <- function(marginal, n_samples = 1000,
                            mintro_grid = seq(0, 1, length.out = 100),
                            seed = 123) {
  if (!is.null(seed)) set.seed(seed)
  beta_draws <- INLA::inla.rmarginal(n_samples, marginal)
  effect_samples <- sapply(beta_draws, function(b) exp(b * mintro_grid))
  data.frame(
    mintro = mintro_grid,
    effect = apply(effect_samples, 1, mean),
    lower  = apply(effect_samples, 1, quantile, probs = 0.025),
    upper  = apply(effect_samples, 1, quantile, probs = 0.975)
  )
}

#' Model-based (g-computation) standardized ratio of incidences.
#'
#' Draws from the joint posterior and, for each draw, predicts total expected
#' cases under the observed introgression versus a counterfactual with
#' introgression set to zero (subtracting `beta * mintro` from the linear
#' predictor), holding the intercept, spatial field and AR1 fixed. The
#' standardized RI is the ratio of these totals; the prevented fraction is
#' `1 - RI`.
#'
#' @param fit A fitted inla object (must have been fit with config = TRUE).
#' @param stack The inla.stack used to fit `fit` (for the observation index).
#' @param intro_values Observed introgression vector, in the same row order as
#'   the stacked observations (e.g. `data$mintro`).
#' @param coef_name Name of the introgression coefficient in the latent field
#'   (default "mintro"); matched with a leading anchor.
#' @param n_draws Number of posterior draws (default 2000).
#' @param seed Optional RNG seed for reproducibility.
#' @return A list with `ri` and `prevented_fraction`, each a named numeric
#'   vector of mean, median (RI only) and 95% credible limits.
standardized_ri <- function(fit, stack, intro_values,
                            coef_name = "mintro", n_draws = 2000, seed = 1) {
  if (!is.null(seed)) set.seed(seed)
  samp <- INLA::inla.posterior.sample(n_draws, fit)

  idx <- INLA::inla.stack.index(stack, tag = "stdata")$data
  ln <- rownames(samp[[1]]$latent)
  apred_rows <- grep("^APredictor", ln)
  beta_row <- grep(paste0("^", coef_name), ln)

  ri_one <- function(s) {
    eta_obs <- s$latent[apred_rows, 1][idx]
    beta <- s$latent[beta_row, 1]
    mu_obs <- exp(eta_obs)
    mu_0 <- exp(eta_obs - beta * intro_values)
    sum(mu_obs) / sum(mu_0)
  }

  ri_draws <- vapply(samp, ri_one, numeric(1))
  pf <- 1 - ri_draws

  list(
    ri = c(mean = mean(ri_draws), median = median(ri_draws),
           quantile(ri_draws, c(0.025, 0.975))),
    prevented_fraction = c(mean = mean(pf), quantile(pf, c(0.025, 0.975)))
  )
}
