# =============================================================================
# 04_simulation_figures.R
# -----------------------------------------------------------------------------
# Build Figure 5 and Supplementary Figures S5-S8 from cached SEIR-SEI
# simulation outputs produced by analysis/03_run_simulations.R.
#
# Each scenario block below restates its parameters and simulation length for
# reference and then loads the corresponding cached trajectories with read_rds():
#   trajectories_earlyintro  early intro, near-fully susceptible
#   trajectories_midintro    early intro, half susceptible
#   trajectories_lateintro   late intro, near equilibrium
#   trajectories_ampli       early intro, temperature-amplitude shock
#   trajectories_amplilate   late intro, temperature-amplitude shock
# The figures depend only on these cached trajectories plus a few timing
# constants (year_introduction, t_break, results_grid, ...) defined per block.
#
# Requires the caches in data/interim/ (run analysis/03_run_simulations.R first).
# ggsave() lines are left commented; uncomment to write to outputs/.
# =============================================================================

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(dplyr, readr, tidyr, ggplot2, patchwork, colorspace, grid)

source(file.path("R", "temperature_functions.R"))
source(file.path("R", "seirsei_model.R"))
source(file.path("R", "theme.R"))

.needed <- c(
  "data/interim/seirsei_fig5-early intro suscept.rds",
  "data/interim/seirsei_fig5-early intro half suscept.rds",
  "data/interim/seirsei_fig5-late intro suscept.rds",
  "data/interim/seirsei_fig5-amplitude change.rds",
  "data/interim/seirsei_fig5-amplitude change late.rds"
)
if (!all(file.exists(.needed))) {
  stop("Missing simulation caches. Run analysis/03_run_simulations.R first.")
}

# ---- Simulations with early introduction of wMel (near-fully susceptible population) ----

# -- Define model parameters
pars = c(N_h = 1e6,              # Host population size.
         S_h = 0.999,            # Starting proportion of susceptible hosts.
         I_h = 0.001,            # Starting proportion of infected hosts.
         M = 4,                  # Ratio of mosquitoes to hosts.
         mu_h = 1 / (77 * 365),  # Host mortality rate.
         r_h = 1 / (77 * 365),   # Host birth rate.
         gamma_h = 1 / 5,        # Host recovery from infectiousness.
         delta_h = 1 / 5.9,      # Intrinsic incubation period.
         # Constant rate of introduction of mosquitoes, as a proportion of `N_h`.
         trickle = 1e-5)

# -- State variables
xstart = vector('numeric', length = 8)
# Sh, Eh, Ih, Rh: susceptible, exposed, infected, recovered hosts.
# Rhc: cumulative recovered hosts (proxy for incidence).
# Sv, Ev, Iv: susceptible, exposed, infected vectors.

names(xstart) = c('Sh', 'Eh', 'Ih', 'Rh', 'Rhc', 'Sv', 'Ev', 'Iv')

xstart['Sv'] = pars['N_h'] * pars['M']
xstart['Sh'] = pars['S_h'] * pars['N_h']
xstart['Ih'] = pars['I_h'] * pars['N_h']
xstart['Rh'] = pars['N_h'] - xstart['Sh'] - xstart['Ih']

make_xstart <- function(pars_vec) {
  x0 <- setNames(numeric(8), c("Sh","Eh","Ih","Rh","Rhc","Sv","Ev","Iv"))
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

# -- Length of simulation
# Total number of days to simulate.
n = 365 * 15
# Need +1 here or else the last time-step is (always?) returned as NA.
ts = seq_len(n + 1)

# -- Define Ms and wMel % to test
M_grid <- seq(1, 10, by = 1)
reduct_fracs <- seq(0.1, 1, by = 0.1)

results_grid <- expand.grid(
  M0    = M_grid,
  frac  = reduct_fracs,
  KEEP.OUT.ATTRS = FALSE
)

# -- Define temperature function
# Option (i): Arbitrary made-up example time series of temperatures.

#' Get a temperature time series based on a sine curve.
#'
#' @param x :numeric (vector): time point(s) for which to get the temperature.
#' @param mean_temp :numeric: mean temperature in C.
#' @param ampl_temp :numeric: peak-to-peak amplitude of temperature.
GetTemp = function(x, mean_temp = 26, ampl_temp = 3)
{
  tmp = ampl_temp * sin(2 * pi * x / 365) + mean_temp
  names(tmp) = NULL
  return(tmp)
}

# Option (ii): Interpolating an actual temperature time series. Here, I pretend
# I have temperatures using the same function above.
# temps = 5 / 2 * sin(2 * pi * ts / 365) + 26 + rnorm(n, 0, 1)
# GetTemp = approxfun(temps)

# -- Run model
trajectories_earlyintro <- read_rds("data/interim/seirsei_fig5-early intro suscept.rds")

# -- Plot outputs

# -- Time series of infections
trajectories_earlyintro |>
  filter(time/365 < 10) |>
  ggplot(aes(x = time/365, y = Ih, color = factor(100-frac*100))) +
  geom_line() +
  labs(x = "Year", y = "Infections",
       color = "%wMel at intervention") +
  scale_x_continuous(breaks = c(1:10)) +
  facet_wrap(~M0, labeller = labeller(M0 = function(x) paste0("M0 = ", x)), ncol = 5) +
  theme_bw() +
  coord_cartesian(ylim = c(0, 250)) +
  theme(legend.position = "bottom",
        legend.byrow = TRUE)

trajectories_earlyintro |>
  filter((100-frac*100) %in% c(0, 50, 90)) |>
  filter(time/365 < 10) |>
  ggplot(aes(x = time/365, y = Ih, color = factor(100-frac*100))) +
  geom_line() +
  labs(x = "Year", y = "Infections",
       color = "%wMel at intervention") +
  scale_x_continuous(breaks = c(1:10)) +
  facet_wrap(~M0, labeller = labeller(M0 = function(x) paste0("M0 = ", x)), ncol = 5) +
  theme_bw() +
  coord_cartesian(ylim = c(0, 250)) +
  theme(legend.position = "bottom",
        legend.byrow = TRUE)

trajectories_earlyintro |>
  filter((100-frac*100) %in% c(0, 50, 90)) |>
  filter(time/365 < 10) |>
  ggplot(aes(x = time/365, y = Ih, color = factor(100-frac*100))) +
  geom_line() +
  scale_x_continuous(breaks = c(1:10)) +
  labs(x = "Year", y = "Infections",
       color = "%wMel at intervention") +
  facet_wrap(~M0, scales = "free_y", labeller = labeller(M0 = function(x) paste0("M0 = ", x)), ncol = 5) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.byrow = TRUE)

trajectories_earlyintro |>
  filter(time/365 < 10) |>
  ggplot(aes(x = time/365, y = Ih, color = factor(100-frac*100))) +
  geom_line() +
  scale_x_continuous(breaks = c(1:10)) +
  labs(x = "Year", y = "Infections",
       color = "%wMel at intervention") +
  facet_wrap(~M0, scales = "free_y", labeller = labeller(M0 = function(x) paste0("M0 = ", x)), ncol = 5) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.byrow = TRUE)

trajectories_earlyintro |>
  filter(M0 %in% c(1, 6) & (100-frac*100) %in% c(0, 50, 90)) |>
  filter(time/365 < 10) |>
  ggplot(aes(x = time/365, y = Ih, color = factor(100-frac*100))) +
  geom_line(lwd = 0.8) +
  scale_x_continuous(breaks = c(0:10)) +
  scale_color_brewer(palette = "Set1") +
  labs(x = "Year", y = "Infections",
       color = "%wMel at intervention") +
  facet_wrap(~M0, scales = "free_y", labeller = labeller(M0 = function(x) paste0("M0 = ", x)), ncol = 5) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.byrow = TRUE) +
  coord_cartesian(xlim = c(0,3))

# -- Yearly and cumulative incidence
traj_rio_earlyintro <- trajectories_earlyintro |>
  ungroup() |>
  mutate(year = floor((time - 1) / 365)) |>
  group_by(M0, Mnew, frac, year) |>
  summarise(Ih = sum(Ih),
            Rhc = sum(Rhc),
            .groups = "drop") |>
  arrange(M0, Mnew, frac, year) |>
  group_by(M0, year) |>
  arrange(desc(Mnew)) |>
  mutate(Ih_intro0 = Ih[1],
         Rhc_intro0 = Rhc[1],
         Ih_intro10 = Ih[2],
         Rhc_intro10 = Rhc[2]) |>
  ungroup() |>
  rowwise() |>
  mutate(inc_ratio_0 = Ih/Ih_intro0,
         inc_ratio_10 = Ih/Ih_intro10,
         cuminc_ratio_0 = Rhc/Rhc_intro0,
         cuminc_ratio_10 = Rhc/Rhc_intro10) |>
  ungroup() |>
  mutate(years_after_intro = year + 1) |>
  filter(years_after_intro <= 10)

spectral_n <- colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(length(unique(traj_rio_earlyintro$years_after_intro)))

# -- Yearly incidence, facets by M0, wMel in x axis
traj_rio_earlyintro |>
  ggplot(aes(x = (1-frac)*100, y = inc_ratio_0, color = factor(years_after_intro))) +
  geom_point() +
  geom_line() +
  geom_hline(aes(yintercept = 1)) +
  facet_wrap(~M0, labeller = labeller(M0 = function(x) paste0("M0 = ", x)), ncol = 5) +
  scale_color_manual(values = spectral_n) +
  scale_x_continuous(breaks = seq(0, 90, 10)) +
  labs(x = "%wMel", y = "Yearly incidence ratio\n(ref: 0% wMel)",
       color = "Years after wMel\nintroduction") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.byrow = TRUE) +
  guides(color = guide_legend(nrow = 2))

traj_rio_earlyintro |>
  ggplot(aes(x = (1-frac)*100, y = inc_ratio_0, color = factor(years_after_intro))) +
  geom_point() +
  geom_line() +
  geom_hline(aes(yintercept = 1)) +
  facet_wrap(~M0, labeller = labeller(M0 = function(x) paste0("M0 = ", x)), ncol = 5) +
  scale_color_manual(values = spectral_n) +
  scale_x_continuous(breaks = seq(0, 90, 10)) +
  labs(x = "%wMel", y = "Yearly incidence ratio\n(ref: 0% wMel)",
       color = "Years after wMel\nintroduction") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.byrow = TRUE) +
  coord_cartesian(ylim = c(0,1)) +
  guides(color = guide_legend(nrow = 2))

# -- Cumulative incidence, facets by M0, wMel in x axis
traj_rio_earlyintro |>
  ggplot(aes(x = (1-frac)*100, y = cuminc_ratio_0, color = factor(years_after_intro))) +
  geom_point() +
  geom_line() +
  geom_hline(aes(yintercept = 1)) +
  facet_wrap(~M0, labeller = labeller(M0 = function(x) paste0("M0 = ", x)), ncol = 5) +
  scale_color_manual(values = spectral_n) +
  scale_x_continuous(breaks = seq(0, 90, 10)) +
  labs(x = "%wMel", y = "Cumulative incidence ratio\n(ref: 0% wMel)",
       color = "Years after wMel\nintroduction") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.byrow = TRUE) +
  guides(color = guide_legend(nrow = 2))

# -- Yearly incidence, facets by M0, years in x axis
spectral_2 <- colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(length(unique(traj_rio_earlyintro$frac)))

traj_rio_earlyintro |>
  ggplot(aes(x = years_after_intro, y = inc_ratio_0, color = factor((1-frac)*100))) +
  geom_point() +
  geom_line() +
  geom_hline(aes(yintercept = 1)) +
  facet_wrap(~M0, labeller = labeller(M0 = function(x) paste0("M0 = ", x)), ncol = 5) +
  scale_color_manual(values = spectral_2) +
  scale_x_continuous(breaks = c(1:10)) +
  labs(x = "Years after wMel introduction", y = "Yearly incidence ratio\n(ref: 0% wMel)",
       color = "% wMel") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.byrow = TRUE)

traj_rio_earlyintro |>
  ggplot(aes(x = years_after_intro, y = inc_ratio_0, color = factor((1-frac)*100))) +
  geom_point() +
  geom_line() +
  geom_hline(aes(yintercept = 1)) +
  facet_wrap(~M0, labeller = labeller(M0 = function(x) paste0("M0 = ", x)), ncol = 5) +
  scale_color_manual(values = spectral_2) +
  scale_x_continuous(breaks = c(1:10)) +
  labs(x = "Years after wMel introduction", y = "Yearly incidence ratio\n(ref: 0% wMel)",
       color = "% wMel") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.byrow = TRUE) +
  coord_cartesian(ylim = c(0, 1))

# -- Cumulative incidence, facets by M0, years in x axis
traj_rio_earlyintro |>
  ggplot(aes(x = years_after_intro, y = cuminc_ratio_0, color = factor((1-frac)*100))) +
  geom_point() +
  geom_line() +
  geom_hline(aes(yintercept = 1)) +
  facet_wrap(~M0, labeller = labeller(M0 = function(x) paste0("M0 = ", x)), ncol = 5) +
  scale_color_manual(values = spectral_2) +
  scale_x_continuous(breaks = c(1:10)) +
  labs(x = "Years after wMel introduction", y = "Cumulative incidence ratio\n(ref: 0% wMel)",
       color = "% wMel") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.byrow = TRUE)

# -- Only plots for M=1 and M=6
(fig5_leg <- traj_rio_earlyintro |>
  filter(M0 == 1) |>
  ggplot(aes(x = years_after_intro, y = inc_ratio_0, color = factor((1-frac)*100))) +
  geom_point() +
  geom_line() +
  geom_hline(aes(yintercept = 1)) +
  scale_color_manual(values = spectral_2) +
  scale_x_continuous(breaks = c(1:10)) +
  labs(x = "Years after wMel introduction", y = "Yearly incidence ratio\n(ref: 0% wMel)",
       color = "% wMel") +
  theme_bw() +
  theme(
    legend.byrow = TRUE))

(fig5_1 <- traj_rio_earlyintro |>
  filter(M0 == 1) |>
  ggplot(aes(x = years_after_intro, y = inc_ratio_0, color = factor((1-frac)*100))) +
  geom_point() +
  geom_line() +
  geom_hline(aes(yintercept = 1)) +
  scale_color_manual(values = spectral_2) +
  scale_x_continuous(breaks = c(1:10)) +
  labs(x = "Years after wMel introduction", y = "Yearly incidence ratio\n(ref: 0% wMel)",
       color = "% wMel") +
  theme_bw() +
  theme(
    legend.position = "none",
    legend.byrow = TRUE))

traj_rio_earlyintro |>
  filter(M0 == 1 & frac %in% c(0.1, 0.5, 1)) |>
  ggplot(aes(x = years_after_intro, y = inc_ratio_0, color = factor((1-frac)*100))) +
  geom_point() +
  geom_line() +
  geom_hline(aes(yintercept = 1)) +
  scale_color_brewer(palette = "Set1") +
  scale_x_continuous(breaks = c(1:10)) +
  labs(x = "Years after wMel introduction", y = "Yearly incidence ratio\n(ref: 0% wMel)",
       color = "% wMel") +
  theme_bw() +
  theme(
    legend.byrow = TRUE)

traj_rio_earlyintro |>
  filter(M0 == 1) |>
  ggplot(aes(x = years_after_intro, y = inc_ratio_0, color = factor((1-frac)*100))) +
  geom_point() +
  geom_line() +
  geom_hline(aes(yintercept = 1)) +
  scale_color_manual(values = spectral_2) +
  scale_x_continuous(breaks = c(1:10)) +
  labs(x = "Years after wMel introduction", y = "Yearly incidence ratio\n(ref: 0% wMel)",
       color = "% wMel") +
  theme_bw() +
  theme(
    legend.position = "none",
    legend.byrow = TRUE) +
  coord_cartesian(ylim = c(0, 1))

traj_rio_earlyintro |>
  filter(M0 == 1 & frac %in% c(0.1, 0.5, 1)) |>
  ggplot(aes(x = years_after_intro, y = inc_ratio_0, color = factor((1-frac)*100))) +
  geom_point() +
  geom_line() +
  geom_hline(aes(yintercept = 1)) +
  scale_color_brewer(palette = "Set1") +
  scale_x_continuous(breaks = c(1:10)) +
  labs(x = "Years after wMel introduction", y = "Yearly incidence ratio\n(ref: 0% wMel)",
       color = "% wMel") +
  theme_bw() +
  theme(
    legend.byrow = TRUE) +
  coord_cartesian(ylim = c(0, 1))

# ggsave("outputs/fig5_1-feb16.pdf", plot = fig5_1, width = 3, height = 2)

(fig5_2 <- traj_rio_earlyintro |>
  filter(M0 == 6) |>
  ggplot(aes(x = years_after_intro, y = inc_ratio_0, color = factor((1-frac)*100))) +
  geom_point() +
  geom_line() +
  geom_hline(aes(yintercept = 1)) +
  scale_color_manual(values = spectral_2) +
  scale_x_continuous(breaks = c(1:10)) +
  scale_y_continuous(
    ) +
  labs(x = "Years after wMel introduction", y = "Yearly incidence ratio\n(ref: 0% wMel)",
       color = "% wMel") +
  theme_bw() +
  theme(
    legend.position = "none",
    legend.byrow = TRUE))

traj_rio_earlyintro |>
  filter(M0 == 6 & frac %in% c(0.1, 0.5, 1)) |>
  ggplot(aes(x = years_after_intro, y = inc_ratio_0, color = factor((1-frac)*100))) +
  geom_point() +
  geom_line() +
  geom_hline(aes(yintercept = 1)) +
  scale_color_brewer(palette = "Set1") +
  scale_x_continuous(breaks = c(1:10)) +
  scale_y_continuous(
    # limits = c(0, 1),
    ) +
  labs(x = "Years after wMel introduction", y = "Yearly incidence ratio\n(ref: 0% wMel)",
       color = "% wMel") +
  theme_bw() +
  theme(
    legend.byrow = TRUE)

traj_rio_earlyintro |>
  filter(M0 == 6) |>
  ggplot(aes(x = years_after_intro, y = inc_ratio_0, color = factor((1-frac)*100))) +
  geom_point() +
  geom_line() +
  geom_hline(aes(yintercept = 1)) +
  scale_color_manual(values = spectral_2) +
  scale_x_continuous(breaks = c(1:10)) +
  scale_y_continuous(
    ) +
  labs(x = "Years after wMel introduction", y = "Yearly incidence ratio\n(ref: 0% wMel)",
       color = "% wMel") +
  theme_bw() +
  theme(
    legend.position = "none",
    legend.byrow = TRUE) +
  coord_cartesian(ylim = c(0, 1))

traj_rio_earlyintro |>
  filter(M0 == 6 & frac %in% c(0.1, 0.5, 1)) |>
  ggplot(aes(x = years_after_intro, y = inc_ratio_0, color = factor((1-frac)*100))) +
  geom_point() +
  geom_line() +
  geom_hline(aes(yintercept = 1)) +
  scale_color_brewer(palette = "Set1") +
  scale_x_continuous(breaks = c(1:10)) +
  scale_y_continuous(
    ) +
  labs(x = "Years after wMel introduction", y = "Yearly incidence ratio\n(ref: 0% wMel)",
       color = "% wMel") +
  theme_bw() +
  theme(
    legend.byrow = TRUE) +
  coord_cartesian(ylim = c(0, 1))

# ggsave("outputs/fig5_2-feb16.pdf", plot = fig5_2, width = 3, height = 2)

(fig5_3 <- traj_rio_earlyintro |>
  filter(M0 == 1) |>
  ggplot(aes(x = years_after_intro, y = cuminc_ratio_0, color = factor((1-frac)*100))) +
  geom_point() +
  geom_line() +
  geom_hline(aes(yintercept = 1)) +
  scale_color_manual(values = spectral_2) +
  scale_x_continuous(breaks = c(1:10)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Years after wMel introduction", y = "Cumulative incidence ratio\n(ref: 0% wMel)",
       color = "% wMel") +
  theme_bw() +
  theme(legend.position = "none",
        legend.byrow = TRUE))

traj_rio_earlyintro |>
  filter(M0 == 1 & frac %in% c(0.1, 0.5, 1)) |>
  ggplot(aes(x = years_after_intro, y = cuminc_ratio_0, color = factor((1-frac)*100))) +
  geom_point() +
  geom_line() +
  geom_hline(aes(yintercept = 1)) +
  scale_color_brewer(palette = "Set1") +
  scale_x_continuous(breaks = c(1:10)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Years after wMel introduction", y = "Cumulative incidence ratio\n(ref: 0% wMel)",
       color = "% wMel") +
  theme_bw() +
  theme(
      legend.byrow = TRUE)

# ggsave("outputs/fig5_3-feb16.pdf", plot = fig5_3, width = 3, height = 2)

(fig5_4 <- traj_rio_earlyintro |>
  filter(M0 == 6) |>
  ggplot(aes(x = years_after_intro, y = cuminc_ratio_0, color = factor((1-frac)*100))) +
  geom_point() +
  geom_line() +
  geom_hline(aes(yintercept = 1)) +
  scale_color_manual(values = spectral_2) +
  scale_x_continuous(breaks = c(1:10)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Years after wMel introduction", y = "Cumulative incidence ratio\n(ref: 0% wMel)",
       color = "% wMel") +
  theme_bw() +
  theme(legend.position = "none",
        legend.byrow = TRUE))

traj_rio_earlyintro |>
  filter(M0 == 6 & frac %in% c(0.1, 0.5, 1)) |>
  ggplot(aes(x = years_after_intro, y = cuminc_ratio_0, color = factor((1-frac)*100))) +
  geom_point() +
  geom_line() +
  geom_hline(aes(yintercept = 1)) +
  scale_color_brewer(palette = "Set1") +
  scale_x_continuous(breaks = c(1:10)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Years after wMel introduction", y = "Cumulative incidence ratio\n(ref: 0% wMel)",
       color = "% wMel") +
  theme_bw() +
  theme(
      legend.byrow = TRUE)

# ggsave("outputs/fig5_4-feb16.pdf", plot = fig5_4, width = 3, height = 2)

# -- Tileplots comparing to true counterfactual
traj_grid_earlyintro <- trajectories_earlyintro |>
  ungroup() |>
  mutate(year = floor((time - 1) / 365),
         years_after_intro = year+1) |>
  group_by(M0, Mnew, frac, year, years_after_intro) |>
  summarise(Ih = sum(Ih),
            Rhc = sum(Rhc),
            .groups = "drop") |>
  arrange(M0, Mnew, frac, year, years_after_intro) |>
  group_by(M0, Mnew, frac) |>
  mutate(Rhc_post = cumsum(if_else(years_after_intro >= 1, Ih, 0))) |>
  ungroup() |>
  filter(year %in% c(0:10)) |>
  group_by(M0, year) |>
  arrange(desc(Mnew)) |>
  mutate(Ih_intro0 = Ih[1],
         Rhc_intro0 = Rhc_post[1],
         Ih_intro10 = Ih[2],
         Rhc_intro10 = Rhc_post[2]) |>
  ungroup() |>
  rowwise() |>
  mutate(inc_ratio_0 = Ih/Ih_intro0,
         inc_ratio_10 = Ih/Ih_intro10,
         cuminc_ratio_0 = Rhc_post/Rhc_intro0,
         cuminc_ratio_10 = Rhc_post/Rhc_intro10) |>
  ungroup() |>
  filter(years_after_intro < 10) |>
  mutate(t_eff_yearly = (1 - inc_ratio_0) * 100,
         t_eff_cumul = (1 - cuminc_ratio_0) * 100)

traj_grid_earlyintro |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = inc_ratio_0)) +
  geom_tile(data = traj_grid_earlyintro |> filter((1-frac)*100 == 0),
              aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "red", width = 10, lwd = 1) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_viridis_c(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Ratio of mosquitoes to hosts",
       fill = "Incidence relative\nto 0% wMel") +
  theme_bw()

traj_grid_earlyintro |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = inc_ratio_0)) +
  geom_tile(data = traj_grid_earlyintro |> filter((1-frac)*100 == 0),
              aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "red", width = 10, lwd = 1) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_viridis_c(trans = "log1p", breaks = c(0, 1, 100, 800)) +
    scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Ratio of mosquitoes to hosts",
       fill = "Incidence relative\nto 0% wMel") +
  theme_bw()

(n1 <- traj_grid_earlyintro |>
  filter(years_after_intro < 4) |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = inc_ratio_0)) +
  geom_tile(data = traj_grid_earlyintro |> filter((1-frac)*100 == 0 & years_after_intro < 4),
              aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "red", width = 10, lwd = 1) +
  geom_tile(data = traj_grid_earlyintro |> filter((1-frac)*100 == 60 & years_after_intro < 4 & M0 ==3),
            aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "white", width = 10, lwd = 1) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_viridis_c(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "Incidence relative\nto 0% wMel") +
  theme_bw())

# ggsave(plot = n1, "outputs/2026-03-05_fig5pt1.pdf", width = 7.5, height = 2.5)

traj_grid_earlyintro |>
  filter(years_after_intro < 4) |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = t_eff_yearly)) +
  geom_tile(data = traj_grid_earlyintro |> filter((1-frac)*100 == 0 & years_after_intro < 4),
              aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "red", width = 10, lwd = 1) +
  geom_tile(data = traj_grid_earlyintro |> filter((1-frac)*100 == 60 & years_after_intro < 4 & M0 ==3),
            aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "white", width = 10, lwd = 1) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_viridis_c(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
    scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "Yearly\nintervention\neffect\n(%)") +
  theme_bw()

traj_grid_earlyintro |>
  filter((1-frac)*100 == 60 & years_after_intro < 4 & M0 == 3) |>
  dplyr::select(years_after_intro, M0, frac, t_eff_yearly, inc_ratio_0)

(fig_s4c <- traj_grid_earlyintro |>
    ggplot() +
    geom_tile(aes(x = (1-frac)*100, y = M0, fill = cuminc_ratio_0)) +
  geom_tile(data = traj_grid_earlyintro |> filter((1-frac)*100 == 0),
              aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "red", width = 10, lwd = 1) +
    facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
    scale_fill_viridis_c(limits = c(0, NA), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    scale_y_continuous(breaks = c(1:10)) +
    labs(x = "% wMel", y = "Ratio of mosquitoes to hosts",
         fill = "Cumulative incidence\nrelative to 0% wMel") +
    theme_bw())

(n2 <- traj_grid_earlyintro |>
  filter(years_after_intro < 4) |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = cuminc_ratio_0)) +
  geom_tile(data = traj_grid_earlyintro |> filter((1-frac)*100 == 0 & years_after_intro < 4),
              aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "red", width = 10, lwd = 1) +
  geom_tile(data = traj_grid_earlyintro |> filter((1-frac)*100 == 60 & years_after_intro < 4 & M0 ==3),
            aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "white", width = 10, lwd = 1) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_viridis_c(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  # scale_fill_viridis_c() +
    scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "Cumulative incidence\nrelative to 0% wMel") +
  theme_bw())

# ggsave(plot = n2, "outputs/2026-03-05_fig5pt2.pdf", width = 7.5, height = 2.5)

traj_grid_earlyintro |>
  filter(years_after_intro < 4) |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = t_eff_cumul)) +
  geom_tile(data = traj_grid_earlyintro |> filter((1-frac)*100 == 0 & years_after_intro < 4),
              aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "red", width = 10, lwd = 1) +
  geom_tile(data = traj_grid_earlyintro |> filter((1-frac)*100 == 60 & years_after_intro < 4 & M0 ==3),
            aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "white", width = 10, lwd = 1) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_viridis_c(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
    scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "Cumulative\nintervention\neffect\n(%)") +
  theme_bw()

traj_grid_earlyintro |>
  filter((1-frac)*100 == 60 & years_after_intro < 4 & M0 == 3) |>
  dplyr::select(years_after_intro, M0, frac, t_eff_cumul, cuminc_ratio_0)

# -- Tileplots comparing to 0% wMel in M0=8
traj_grid_biased_earlyintro <- trajectories_earlyintro |>
  ungroup() |>
  mutate(year = floor((time - 1) / 365),
         years_after_intro = year+1) |>
  group_by(M0, Mnew, frac, year, years_after_intro) |>
  summarise(Ih = sum(Ih),
            Rhc = sum(Rhc),
            .groups = "drop") |>
  arrange(M0, Mnew, frac, year, years_after_intro) |>
  group_by(M0, Mnew, frac) |>
  mutate(Rhc_post = cumsum(if_else(years_after_intro >= 0, Ih, 0))) |>
  ungroup() |>
  filter(year %in% c(0:10)) |>
  group_by(year) |>
  arrange(desc(M0)) |>
  mutate(Ih_intro0 = Ih[30],
         Rhc_intro0 = Rhc_post[30]) |>
  ungroup() |>
  rowwise() |>
  mutate(inc_ratio_0 = Ih/Ih_intro0,
         cuminc_ratio_0 = Rhc_post/Rhc_intro0) |>
  ungroup() |>
  filter(years_after_intro < 10) |>
  mutate(t_eff_yearly = (1 - inc_ratio_0) * 100,
         t_eff_cumul = (1 - cuminc_ratio_0) * 100)

traj_grid_biased_earlyintro |>
    ggplot() +
    geom_tile(aes(x = (1-frac)*100, y = M0, fill = inc_ratio_0)) +
    geom_tile(data = traj_grid_biased_earlyintro |> filter((1-frac)*100 == 0),
              aes(x = (1-frac)*100, y = 8),
              fill = NA, color = "red", width = 10, lwd = 1) +
    facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
    scale_fill_viridis_c(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    scale_y_continuous(breaks = c(1:10)) +
    labs(x = "% wMel", y = "Ratio of mosquitoes to hosts",
         fill = "Incidence relative\nto 0% wMel for M=8") +
    theme_bw()

traj_grid_biased_earlyintro |>
    ggplot() +
    geom_tile(aes(x = (1-frac)*100, y = M0, fill = inc_ratio_0)) +
    geom_tile(data = traj_grid_biased_earlyintro |> filter((1-frac)*100 == 0),
              aes(x = (1-frac)*100, y = 8),
              fill = NA, color = "red", width = 10, lwd = 1) +
    facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
    scale_fill_viridis_c(trans = "log1p", breaks = c(0, 1, 100, 800)) +
    scale_y_continuous(breaks = c(1:10)) +
    labs(x = "% wMel", y = "Ratio of mosquitoes to hosts",
         fill = "Incidence relative\nto 0% wMel for M=8") +
    theme_bw()

traj_grid_biased_earlyintro |>
    ggplot() +
    geom_tile(aes(x = (1-frac)*100, y = M0, fill = cuminc_ratio_0)) +
    geom_tile(data = traj_grid_biased_earlyintro |> filter((1-frac)*100 == 0),
              aes(x = (1-frac)*100, y = 8),
              fill = NA, color = "red", width = 10, lwd = 1) +
    facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
    scale_fill_viridis_c(limits = c(0, NA), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    scale_y_continuous(breaks = c(1:10)) +
    labs(x = "% wMel", y = "Ratio of mosquitoes to hosts",
         fill = "Cumulative incidence\nrelative to 0% wMel for M=8") +
    theme_bw()

# -- Tileplots comparing to 0% wMel in M0=6
traj_grid_biased_earlyintro2 <- trajectories_earlyintro |>
  ungroup() |>
  mutate(year = floor((time - 1) / 365),
         years_after_intro = year+1) |>
  group_by(M0, Mnew, frac, year, years_after_intro) |>
  summarise(Ih = sum(Ih),
            Rhc = sum(Rhc),
            .groups = "drop") |>
  arrange(M0, Mnew, frac, year, years_after_intro) |>
  group_by(M0, Mnew, frac) |>
  mutate(Rhc_post = cumsum(if_else(years_after_intro >= 0, Ih, 0))) |>
  ungroup() |>
  filter(year %in% c(0:10)) |>
  group_by(year) |>
  arrange(desc(M0)) |>
  mutate(Ih_intro0 = Ih[50],
         Rhc_intro0 = Rhc_post[50]) |>
  ungroup() |>
  rowwise() |>
  mutate(inc_ratio_0 = Ih/Ih_intro0,
         cuminc_ratio_0 = Rhc_post/Rhc_intro0) |>
  ungroup() |>
  filter(years_after_intro < 10) |>
  mutate(t_eff_yearly = (1 - inc_ratio_0) * 100,
         t_eff_cumul = (1 - cuminc_ratio_0) * 100)

traj_grid_biased_earlyintro2 |>
    ggplot() +
    geom_tile(aes(x = (1-frac)*100, y = M0, fill = inc_ratio_0)) +
    geom_tile(data = traj_grid_biased_earlyintro2 |> filter((1-frac)*100 == 0),
              aes(x = (1-frac)*100, y = 6),
              fill = NA, color = "red", width = 10, lwd = 1) +
    facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
    scale_fill_viridis_c(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    scale_y_continuous(breaks = c(1:10)) +
    labs(x = "% wMel", y = "Ratio of mosquitoes to hosts",
         fill = "Incidence relative\nto 0% wMel for M=6") +
    theme_bw()

traj_grid_biased_earlyintro2 |>
    ggplot() +
    geom_tile(aes(x = (1-frac)*100, y = M0, fill = inc_ratio_0)) +
    geom_tile(data = traj_grid_biased_earlyintro2 |> filter((1-frac)*100 == 0),
              aes(x = (1-frac)*100, y = 6),
              fill = NA, color = "red", width = 10, lwd = 1) +
    facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
    scale_fill_viridis_c(trans = "log1p", breaks = c(0, 1, 100, 800)) +
    scale_y_continuous(breaks = c(1:10)) +
    labs(x = "% wMel", y = "Ratio of mosquitoes to hosts",
         fill = "Incidence relative\nto 0% wMel for M=6") +
    theme_bw()

(n3 <- traj_grid_biased_earlyintro2 |>
  filter(years_after_intro < 4) |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = inc_ratio_0)) +
  geom_tile(data = traj_grid_biased_earlyintro2 |> filter((1-frac)*100 == 0 & years_after_intro < 4),
              aes(x = (1-frac)*100, y = 6),
              fill = NA, color = "red", width = 10, lwd = 1) +
  geom_tile(data = traj_grid_biased_earlyintro2 |> filter((1-frac)*100 == 60 & years_after_intro < 4 & M0 == 3),
            aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "white", width = 10, lwd = 1) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_viridis_c(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "Incidence relative\nto 0% wMel") +
  theme_bw())

# ggsave(plot = n3, "outputs/2026-03-05_fig5pt3.pdf", width = 7.5, height = 2.5)

traj_grid_biased_earlyintro2 |>
  filter(years_after_intro < 4) |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = t_eff_yearly)) +
  geom_tile(data = traj_grid_biased_earlyintro2 |> filter((1-frac)*100 == 0 & years_after_intro < 4),
              aes(x = (1-frac)*100, y = 6),
              fill = NA, color = "red", width = 10, lwd = 1) +
  geom_tile(data = traj_grid_biased_earlyintro2 |> filter((1-frac)*100 == 60 & years_after_intro < 4 & M0 == 3),
            aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "white", width = 10, lwd = 1) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_viridis_c(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
    scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "Yearly\nintervention\neffect\n(%)") +
  theme_bw()

traj_grid_biased_earlyintro2 |>
  filter((1-frac)*100 == 60 & years_after_intro < 4 & M0 == 3) |>
  dplyr::select(years_after_intro, M0, frac, t_eff_yearly, inc_ratio_0)

traj_grid_biased_earlyintro2 |>
    ggplot() +
    geom_tile(aes(x = (1-frac)*100, y = M0, fill = cuminc_ratio_0)) +
    geom_tile(data = traj_grid_biased_earlyintro2 |> filter((1-frac)*100 == 0),
              aes(x = (1-frac)*100, y = 6),
              fill = NA, color = "red", width = 10, lwd = 1) +
    facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
    scale_fill_viridis_c(limits = c(0, NA), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    scale_y_continuous(breaks = c(1:10)) +
    labs(x = "% wMel", y = "Ratio of mosquitoes to hosts",
         fill = "Cumulative incidence\nrelative to 0% wMel for M=6") +
    theme_bw()

(n4 <- traj_grid_biased_earlyintro2|>
  filter(years_after_intro < 4) |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = cuminc_ratio_0)) +
  geom_tile(data = traj_grid_biased_earlyintro2 |> filter((1-frac)*100 == 0 & years_after_intro < 4),
              aes(x = (1-frac)*100, y = 6),
              fill = NA, color = "red", width = 10, lwd = 1) +
  geom_tile(data = traj_grid_biased_earlyintro2|> filter((1-frac)*100 == 60 & years_after_intro < 4 & M0 ==3),
            aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "white", width = 10, lwd = 1) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_viridis_c(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "Cumulative incidence\nrelative to 0% wMel") +
  theme_bw())

# ggsave(plot = n4, "outputs/2026-03-05_fig5pt4.pdf", width = 7.5, height = 2.5)

traj_grid_biased_earlyintro2 |>
  filter(years_after_intro < 4) |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = t_eff_cumul)) +
  geom_tile(data = traj_grid_biased_earlyintro2 |> filter((1-frac)*100 == 0 & years_after_intro < 4),
              aes(x = (1-frac)*100, y = 6),
              fill = NA, color = "red", width = 10, lwd = 1) +
  geom_tile(data = traj_grid_biased_earlyintro2 |> filter((1-frac)*100 == 60 & years_after_intro < 4 & M0 ==3),
            aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "white", width = 10, lwd = 1) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_viridis_c(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
    scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "Cumulative\nintervention\neffect\n(%)") +
  theme_bw()

traj_grid_biased_earlyintro2 |>
  filter((1-frac)*100 == 60 & years_after_intro < 4 & M0 == 3) |>
  dplyr::select(years_after_intro, M0, frac, t_eff_cumul, cuminc_ratio_0)

# -- Tileplots comparing true to biased (M0 = 6) counterfactual
comp_earlyintro <- traj_grid_earlyintro |>
  mutate(wmel = (1-frac)*100,
         true_incratio = inc_ratio_0,
         true_cumincratio = cuminc_ratio_0) |>
  dplyr::select(M0, wmel, true_incratio, true_cumincratio, years_after_intro) |>
  left_join(traj_grid_biased_earlyintro2 |>
              mutate(wmel = (1-frac)*100,
                     biased_incratio = inc_ratio_0,
                     biased_cumincratio = cuminc_ratio_0) |>
              dplyr::select(M0, wmel, biased_incratio, biased_cumincratio, years_after_intro),
            by = c("M0", "wmel", "years_after_intro")) |>
  ungroup() |>
  mutate(diff_incratio = true_incratio - biased_incratio,
         diff_cumincratio = true_cumincratio - biased_cumincratio,
         ratio_incratio = true_incratio / biased_incratio,
         ratio_cumincratio = true_cumincratio / biased_cumincratio)

## yearly
comp_earlyintro |>
  ggplot() +
  geom_tile(aes(x = wmel, y = M0, fill = diff_incratio)) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_distiller(palette = "Spectral",
                       limits = c(-max(abs(comp_earlyintro$diff_incratio), na.rm = TRUE),
                                  max(abs(comp_earlyintro$diff_incratio), na.rm = TRUE))) +
  scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Ratio of mosquitoes to hosts",
       fill = "True - biased yearly incidence") +
  theme_bw()

comp_earlyintro |>
  mutate(diff_incratio = cut(diff_incratio, breaks = c(-650, -10, -1, -0.5, 0.5, 1, 10, 650))) |>
  ggplot() +
  geom_tile(aes(x = wmel, y = M0, fill = diff_incratio)) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_brewer(palette = "Spectral") +
  scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Ratio of mosquitoes to hosts",
       fill = "True - biased yearly incidence") +
  theme_bw()

comp_earlyintro |>
  mutate(diff_incratio = factor(case_when(
    diff_incratio > 0 ~ "True incidence higher",
    diff_incratio == 0 ~ "No bias",
    diff_incratio < 0 ~ "True incidence lower"),
    levels = c("True incidence higher", "No bias", "True incidence lower"))) |>
  ggplot() +
  geom_tile(aes(x = wmel, y = M0, fill = diff_incratio)) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_brewer(palette = "Spectral") +
  scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Ratio of mosquitoes to hosts",
       fill = "True - biased yearly incidence") +
  theme_bw()

## cumulative
comp_earlyintro |>
  ggplot() +
  geom_tile(aes(x = wmel, y = M0, fill = diff_cumincratio)) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_distiller(palette = "Spectral",
                       limits = c(-max(abs(comp_earlyintro$diff_cumincratio), na.rm = TRUE),
                                  max(abs(comp_earlyintro$diff_cumincratio), na.rm = TRUE))) +
  scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Ratio of mosquitoes to hosts",
       fill = "True - biased cumulative incidence") +
  theme_bw()

comp_earlyintro |>
  mutate(diff_cumincratio = cut(
    diff_cumincratio,
    breaks = c(-0.2, -0.1, -0.05, -0.01, 0.01, 0.05, 0.1, 0.2))) |>
  ggplot() +
  geom_tile(aes(x = wmel, y = M0, fill = diff_cumincratio)) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_brewer(palette = "Spectral", drop = FALSE) +
  scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Ratio of mosquitoes to hosts",
       fill = "True - biased cumulative incidence") +
  theme_bw()

comp_earlyintro |>
  mutate(diff_cumincratio = factor(case_when(
    diff_cumincratio > 0 ~ "True incidence higher",
    diff_cumincratio == 0 ~ "No bias",
    diff_cumincratio < 0 ~ "True incidence lower"),
    levels = c("True incidence higher", "No bias", "True incidence lower"))) |>
  ggplot() +
  geom_tile(aes(x = wmel, y = M0, fill = diff_cumincratio)) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_brewer(palette = "Spectral") +
  scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Ratio of mosquitoes to hosts",
       fill = "True - biased cumulative incidence") +
  theme_bw()

# ---- Simulations with early introduction of wMel (half susceptible population) ----

# -- Define model parameters
pars = c(N_h = 1e6,              # Host population size.
         S_h = 0.5,            # Starting proportion of susceptible hosts.
         I_h = 0.001,            # Starting proportion of infected hosts.
         M = 4,                  # Ratio of mosquitoes to hosts.
         mu_h = 1 / (77 * 365),  # Host mortality rate.
         r_h = 1 / (77 * 365),   # Host birth rate.
         gamma_h = 1 / 5,        # Host recovery from infectiousness.
         delta_h = 1 / 5.9,      # Intrinsic incubation period.
         # Constant rate of introduction of mosquitoes, as a proportion of `N_h`.
         trickle = 1e-5)

# -- State variables
xstart = vector('numeric', length = 8)
# Sh, Eh, Ih, Rh: susceptible, exposed, infected, recovered hosts.
# Rhc: cumulative recovered hosts (proxy for incidence).
# Sv, Ev, Iv: susceptible, exposed, infected vectors.

names(xstart) = c('Sh', 'Eh', 'Ih', 'Rh', 'Rhc', 'Sv', 'Ev', 'Iv')

xstart['Sv'] = pars['N_h'] * pars['M']
xstart['Sh'] = pars['S_h'] * pars['N_h']
xstart['Ih'] = pars['I_h'] * pars['N_h']
xstart['Rh'] = pars['N_h'] - xstart['Sh'] - xstart['Ih']

make_xstart <- function(pars_vec) {
  x0 <- setNames(numeric(8), c("Sh","Eh","Ih","Rh","Rhc","Sv","Ev","Iv"))
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

# -- Length of simulation
# Total number of days to simulate.
n = 365 * 15
# Need +1 here or else the last time-step is (always?) returned as NA.
ts = seq_len(n + 1)

# -- Define Ms and wMel % to test
M_grid <- seq(1, 10, by = 1)
reduct_fracs <- seq(0.1, 1, by = 0.1)

results_grid <- expand.grid(
  M0    = M_grid,
  frac  = reduct_fracs,
  KEEP.OUT.ATTRS = FALSE
)

# -- Define temperature function
#' Get a temperature time series based on a sine curve.
#'
#' @param x :numeric (vector): time point(s) for which to get the temperature.
#' @param mean_temp :numeric: mean temperature in C.
#' @param ampl_temp :numeric: peak-to-peak amplitude of temperature.
GetTemp = function(x, mean_temp = 26, ampl_temp = 3)
{
  tmp = ampl_temp * sin(2 * pi * x / 365) + mean_temp
  names(tmp) = NULL
  return(tmp)
}

# -- Run model
trajectories_midintro <- read_rds("data/interim/seirsei_fig5-early intro half suscept.rds")

# -- Plot outputs

# -- Time series of infections
trajectories_midintro |>
  filter(time/365 < 10) |>
  ggplot(aes(x = time/365, y = Ih, color = factor(100-frac*100))) +
  geom_line() +
  labs(x = "Year", y = "Infections",
       color = "%wMel at intervention") +
  scale_x_continuous(breaks = c(1:10)) +
  facet_wrap(~M0, labeller = labeller(M0 = function(x) paste0("M0 = ", x)), ncol = 5) +
  theme_bw() +
  coord_cartesian(ylim = c(0, 250)) +
  theme(legend.position = "bottom",
        legend.byrow = TRUE)

trajectories_midintro |>
  filter((100-frac*100) %in% c(0, 50, 90)) |>
  filter(time/365 < 10) |>
  ggplot(aes(x = time/365, y = Ih, color = factor(100-frac*100))) +
  geom_line() +
  labs(x = "Year", y = "Infections",
       color = "%wMel at intervention") +
  scale_x_continuous(breaks = c(1:10)) +
  facet_wrap(~M0, labeller = labeller(M0 = function(x) paste0("M0 = ", x)), ncol = 5) +
  theme_bw() +
  coord_cartesian(ylim = c(0, 250)) +
  theme(legend.position = "bottom",
        legend.byrow = TRUE)

trajectories_midintro |>
  filter((100-frac*100) %in% c(0, 50, 90)) |>
  filter(time/365 < 10) |>
  ggplot(aes(x = time/365, y = Ih, color = factor(100-frac*100))) +
  geom_line() +
  scale_x_continuous(breaks = c(1:10)) +
  labs(x = "Year", y = "Infections",
       color = "%wMel at intervention") +
  facet_wrap(~M0, scales = "free_y", labeller = labeller(M0 = function(x) paste0("M0 = ", x)), ncol = 5) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.byrow = TRUE)

trajectories_midintro |>
  filter(time/365 < 10) |>
  ggplot(aes(x = time/365, y = Ih, color = factor(100-frac*100))) +
  geom_line() +
  scale_x_continuous(breaks = c(1:10)) +
  labs(x = "Year", y = "Infections",
       color = "%wMel at intervention") +
  facet_wrap(~M0, scales = "free_y", labeller = labeller(M0 = function(x) paste0("M0 = ", x)), ncol = 5) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.byrow = TRUE)

trajectories_midintro |>
  filter(M0 %in% c(1, 2) & (100-frac*100) %in% c(0, 50, 90)) |>
  filter(time/365 < 10) |>
  ggplot(aes(x = time/365, y = Ih, color = factor(100-frac*100))) +
  geom_line(lwd = 0.8) +
  scale_x_continuous(breaks = c(0:10)) +
  scale_color_brewer(palette = "Set1") +
  labs(x = "Year", y = "Infections",
       color = "%wMel at intervention") +
  facet_wrap(~M0, scales = "free_y", labeller = labeller(M0 = function(x) paste0("M0 = ", x)), ncol = 5) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.byrow = TRUE) +
  coord_cartesian(xlim = c(0,5))

# -- Yearly and cumulative incidence
traj_rio_midintro <- trajectories_midintro |>
  ungroup() |>
  mutate(year = floor((time - 1) / 365)) |>
  group_by(M0, Mnew, frac, year) |>
  summarise(Ih = sum(Ih),
            Rhc = sum(Rhc),
            .groups = "drop") |>
  arrange(M0, Mnew, frac, year) |>
  group_by(M0, year) |>
  arrange(desc(Mnew)) |>
  mutate(Ih_intro0 = Ih[1],
         Rhc_intro0 = Rhc[1],
         Ih_intro10 = Ih[2],
         Rhc_intro10 = Rhc[2]) |>
  ungroup() |>
  rowwise() |>
  mutate(inc_ratio_0 = Ih/Ih_intro0,
         inc_ratio_10 = Ih/Ih_intro10,
         cuminc_ratio_0 = Rhc/Rhc_intro0,
         cuminc_ratio_10 = Rhc/Rhc_intro10) |>
  ungroup() |>
  mutate(years_after_intro = year + 1) |>
  filter(years_after_intro <= 10)

spectral_n <- colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(length(unique(traj_rio_midintro$years_after_intro)))

# -- Yearly incidence, facets by M0, wMel in x axis
traj_rio_midintro |>
  ggplot(aes(x = (1-frac)*100, y = inc_ratio_0, color = factor(years_after_intro))) +
  geom_point() +
  geom_line() +
  geom_hline(aes(yintercept = 1)) +
  facet_wrap(~M0, labeller = labeller(M0 = function(x) paste0("M0 = ", x)), ncol = 5) +
  scale_color_manual(values = spectral_n) +
  scale_x_continuous(breaks = seq(0, 90, 10)) +
  labs(x = "%wMel", y = "Yearly incidence ratio\n(ref: 0% wMel)",
       color = "Years after wMel\nintroduction") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.byrow = TRUE) +
  guides(color = guide_legend(nrow = 2))

traj_rio_midintro |>
  ggplot(aes(x = (1-frac)*100, y = inc_ratio_0, color = factor(years_after_intro))) +
  geom_point() +
  geom_line() +
  geom_hline(aes(yintercept = 1)) +
  facet_wrap(~M0, labeller = labeller(M0 = function(x) paste0("M0 = ", x)), ncol = 5) +
  scale_color_manual(values = spectral_n) +
  scale_x_continuous(breaks = seq(0, 90, 10)) +
  labs(x = "%wMel", y = "Yearly incidence ratio\n(ref: 0% wMel)",
       color = "Years after wMel\nintroduction") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.byrow = TRUE) +
  coord_cartesian(ylim = c(0,1)) +
  guides(color = guide_legend(nrow = 2))

# -- Cumulative incidence, facets by M0, wMel in x axis
traj_rio_midintro |>
  ggplot(aes(x = (1-frac)*100, y = cuminc_ratio_0, color = factor(years_after_intro))) +
  geom_point() +
  geom_line() +
  geom_hline(aes(yintercept = 1)) +
  facet_wrap(~M0, labeller = labeller(M0 = function(x) paste0("M0 = ", x)), ncol = 5) +
  scale_color_manual(values = spectral_n) +
  scale_x_continuous(breaks = seq(0, 90, 10)) +
  labs(x = "%wMel", y = "Cumulative incidence ratio\n(ref: 0% wMel)",
       color = "Years after wMel\nintroduction") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.byrow = TRUE) +
  guides(color = guide_legend(nrow = 2))

# -- Yearly incidence, facets by M0, years in x axis
spectral_2 <- colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(length(unique(traj_rio_midintro$frac)))

traj_rio_midintro |>
  ggplot(aes(x = years_after_intro, y = inc_ratio_0, color = factor((1-frac)*100))) +
  geom_point() +
  geom_line() +
  geom_hline(aes(yintercept = 1)) +
  facet_wrap(~M0, labeller = labeller(M0 = function(x) paste0("M0 = ", x)), ncol = 5) +
  scale_color_manual(values = spectral_2) +
  scale_x_continuous(breaks = c(1:10)) +
  labs(x = "Years after wMel introduction", y = "Yearly incidence ratio\n(ref: 0% wMel)",
       color = "% wMel") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.byrow = TRUE)

traj_rio_midintro |>
  ggplot(aes(x = years_after_intro, y = inc_ratio_0, color = factor((1-frac)*100))) +
  geom_point() +
  geom_line() +
  geom_hline(aes(yintercept = 1)) +
  facet_wrap(~M0, labeller = labeller(M0 = function(x) paste0("M0 = ", x)), ncol = 5) +
  scale_color_manual(values = spectral_2) +
  scale_x_continuous(breaks = c(1:10)) +
  labs(x = "Years after wMel introduction", y = "Yearly incidence ratio\n(ref: 0% wMel)",
       color = "% wMel") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.byrow = TRUE) +
  coord_cartesian(ylim = c(0, 1))

# -- Cumulative incidence, facets by M0, years in x axis
traj_rio_midintro |>
  ggplot(aes(x = years_after_intro, y = cuminc_ratio_0, color = factor((1-frac)*100))) +
  geom_point() +
  geom_line() +
  geom_hline(aes(yintercept = 1)) +
  facet_wrap(~M0, labeller = labeller(M0 = function(x) paste0("M0 = ", x)), ncol = 5) +
  scale_color_manual(values = spectral_2) +
  scale_x_continuous(breaks = c(1:10)) +
  labs(x = "Years after wMel introduction", y = "Cumulative incidence ratio\n(ref: 0% wMel)",
       color = "% wMel") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.byrow = TRUE)

# -- Only plots for M=1 and M=6
(traj_rio_midintro |>
  filter(M0 == 1) |>
  ggplot(aes(x = years_after_intro, y = inc_ratio_0, color = factor((1-frac)*100))) +
  geom_point() +
  geom_line() +
  geom_hline(aes(yintercept = 1)) +
  scale_color_manual(values = spectral_2) +
  scale_x_continuous(breaks = c(1:10)) +
  labs(x = "Years after wMel introduction", y = "Yearly incidence ratio\n(ref: 0% wMel)",
       color = "% wMel") +
  theme_bw() +
  theme(
    # legend.position = "none",
    legend.byrow = TRUE))

traj_rio_midintro |>
  filter(M0 == 1) |>
  ggplot(aes(x = years_after_intro, y = inc_ratio_0, color = factor((1-frac)*100))) +
  geom_point() +
  geom_line() +
  geom_hline(aes(yintercept = 1)) +
  scale_color_manual(values = spectral_2) +
  scale_x_continuous(breaks = c(1:10)) +
  labs(x = "Years after wMel introduction", y = "Yearly incidence ratio\n(ref: 0% wMel)",
       color = "% wMel") +
  theme_bw() +
  theme(
    legend.position = "none",
    legend.byrow = TRUE)

traj_rio_midintro |>
  filter(M0 == 1 & frac %in% c(0.1, 0.5, 1)) |>
  ggplot(aes(x = years_after_intro, y = inc_ratio_0, color = factor((1-frac)*100))) +
  geom_point() +
  geom_line() +
  geom_hline(aes(yintercept = 1)) +
  scale_color_brewer(palette = "Set1") +
  scale_x_continuous(breaks = c(1:10)) +
  labs(x = "Years after wMel introduction", y = "Yearly incidence ratio\n(ref: 0% wMel)",
       color = "% wMel") +
  theme_bw() +
  theme(
    legend.byrow = TRUE)

traj_rio_midintro |>
  filter(M0 == 1) |>
  ggplot(aes(x = years_after_intro, y = inc_ratio_0, color = factor((1-frac)*100))) +
  geom_point() +
  geom_line() +
  geom_hline(aes(yintercept = 1)) +
  scale_color_manual(values = spectral_2) +
  scale_x_continuous(breaks = c(1:10)) +
  labs(x = "Years after wMel introduction", y = "Yearly incidence ratio\n(ref: 0% wMel)",
       color = "% wMel") +
  theme_bw() +
  theme(
    legend.position = "none",
    legend.byrow = TRUE) +
  coord_cartesian(ylim = c(0, 1))

traj_rio_midintro |>
  filter(M0 == 1 & frac %in% c(0.1, 0.5, 1)) |>
  ggplot(aes(x = years_after_intro, y = inc_ratio_0, color = factor((1-frac)*100))) +
  geom_point() +
  geom_line() +
  geom_hline(aes(yintercept = 1)) +
  scale_color_brewer(palette = "Set1") +
  scale_x_continuous(breaks = c(1:10)) +
  labs(x = "Years after wMel introduction", y = "Yearly incidence ratio\n(ref: 0% wMel)",
       color = "% wMel") +
  theme_bw() +
  theme(
    legend.byrow = TRUE) +
  coord_cartesian(ylim = c(0, 1))

# ggsave("outputs/fig5_1-feb16.pdf", plot = fig5_1, width = 3, height = 2)

traj_rio_midintro |>
  filter(M0 == 6) |>
  ggplot(aes(x = years_after_intro, y = inc_ratio_0, color = factor((1-frac)*100))) +
  geom_point() +
  geom_line() +
  geom_hline(aes(yintercept = 1)) +
  scale_color_manual(values = spectral_2) +
  scale_x_continuous(breaks = c(1:10)) +
  scale_y_continuous(
    ) +
  labs(x = "Years after wMel introduction", y = "Yearly incidence ratio\n(ref: 0% wMel)",
       color = "% wMel") +
  theme_bw() +
  theme(
    legend.position = "none",
    legend.byrow = TRUE)

traj_rio_midintro |>
  filter(M0 == 6 & frac %in% c(0.1, 0.5, 1)) |>
  ggplot(aes(x = years_after_intro, y = inc_ratio_0, color = factor((1-frac)*100))) +
  geom_point() +
  geom_line() +
  geom_hline(aes(yintercept = 1)) +
  scale_color_brewer(palette = "Set1") +
  scale_x_continuous(breaks = c(1:10)) +
  scale_y_continuous(
    ) +
  labs(x = "Years after wMel introduction", y = "Yearly incidence ratio\n(ref: 0% wMel)",
       color = "% wMel") +
  theme_bw() +
  theme(
    legend.byrow = TRUE)

traj_rio_midintro |>
  filter(M0 == 6) |>
  ggplot(aes(x = years_after_intro, y = inc_ratio_0, color = factor((1-frac)*100))) +
  geom_point() +
  geom_line() +
  geom_hline(aes(yintercept = 1)) +
  scale_color_manual(values = spectral_2) +
  scale_x_continuous(breaks = c(1:10)) +
  scale_y_continuous(
    ) +
  labs(x = "Years after wMel introduction", y = "Yearly incidence ratio\n(ref: 0% wMel)",
       color = "% wMel") +
  theme_bw() +
  theme(
    legend.position = "none",
    legend.byrow = TRUE) +
  coord_cartesian(ylim = c(0, 1))

traj_rio_midintro |>
  filter(M0 == 6 & frac %in% c(0.1, 0.5, 1)) |>
  ggplot(aes(x = years_after_intro, y = inc_ratio_0, color = factor((1-frac)*100))) +
  geom_point() +
  geom_line() +
  geom_hline(aes(yintercept = 1)) +
  scale_color_brewer(palette = "Set1") +
  scale_x_continuous(breaks = c(1:10)) +
  scale_y_continuous(
    ) +
  labs(x = "Years after wMel introduction", y = "Yearly incidence ratio\n(ref: 0% wMel)",
       color = "% wMel") +
  theme_bw() +
  theme(
    legend.byrow = TRUE) +
  coord_cartesian(ylim = c(0, 1))

traj_rio_midintro |>
  filter(M0 == 1) |>
  ggplot(aes(x = years_after_intro, y = cuminc_ratio_0, color = factor((1-frac)*100))) +
  geom_point() +
  geom_line() +
  geom_hline(aes(yintercept = 1)) +
  scale_color_manual(values = spectral_2) +
  scale_x_continuous(breaks = c(1:10)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Years after wMel introduction", y = "Cumulative incidence ratio\n(ref: 0% wMel)",
       color = "% wMel") +
  theme_bw() +
  theme(legend.position = "none",
        legend.byrow = TRUE)

traj_rio_midintro |>
  filter(M0 == 1 & frac %in% c(0.1, 0.5, 1)) |>
  ggplot(aes(x = years_after_intro, y = cuminc_ratio_0, color = factor((1-frac)*100))) +
  geom_point() +
  geom_line() +
  geom_hline(aes(yintercept = 1)) +
  scale_color_brewer(palette = "Set1") +
  scale_x_continuous(breaks = c(1:10)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Years after wMel introduction", y = "Cumulative incidence ratio\n(ref: 0% wMel)",
       color = "% wMel") +
  theme_bw() +
  theme(
      legend.byrow = TRUE)

traj_rio_midintro |>
  filter(M0 == 6) |>
  ggplot(aes(x = years_after_intro, y = cuminc_ratio_0, color = factor((1-frac)*100))) +
  geom_point() +
  geom_line() +
  geom_hline(aes(yintercept = 1)) +
  scale_color_manual(values = spectral_2) +
  scale_x_continuous(breaks = c(1:10)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Years after wMel introduction", y = "Cumulative incidence ratio\n(ref: 0% wMel)",
       color = "% wMel") +
  theme_bw() +
  theme(legend.position = "none",
        legend.byrow = TRUE)

traj_rio_midintro |>
  filter(M0 == 6 & frac %in% c(0.1, 0.5, 1)) |>
  ggplot(aes(x = years_after_intro, y = cuminc_ratio_0, color = factor((1-frac)*100))) +
  geom_point() +
  geom_line() +
  geom_hline(aes(yintercept = 1)) +
  scale_color_brewer(palette = "Set1") +
  scale_x_continuous(breaks = c(1:10)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Years after wMel introduction", y = "Cumulative incidence ratio\n(ref: 0% wMel)",
       color = "% wMel") +
  theme_bw() +
  theme(
      legend.byrow = TRUE)

# -- Tileplots comparing to true counterfactual
traj_grid_midintro <- trajectories_midintro |>
  ungroup() |>
  mutate(year = floor((time - 1) / 365),
         years_after_intro = year+1) |>
  group_by(M0, Mnew, frac, year, years_after_intro) |>
  summarise(Ih = sum(Ih),
            Rhc = sum(Rhc),
            .groups = "drop") |>
  arrange(M0, Mnew, frac, year, years_after_intro) |>
  group_by(M0, Mnew, frac) |>
  mutate(Rhc_post = cumsum(if_else(years_after_intro >= 1, Ih, 0))) |>
  ungroup() |>
  filter(year %in% c(0:10)) |>
  group_by(M0, year) |>
  arrange(desc(Mnew)) |>
  mutate(Ih_intro0 = Ih[1],
         Rhc_intro0 = Rhc_post[1],
         Ih_intro10 = Ih[2],
         Rhc_intro10 = Rhc_post[2]) |>
  ungroup() |>
  rowwise() |>
  mutate(inc_ratio_0 = Ih/Ih_intro0,
         inc_ratio_10 = Ih/Ih_intro10,
         cuminc_ratio_0 = Rhc_post/Rhc_intro0,
         cuminc_ratio_10 = Rhc_post/Rhc_intro10) |>
  ungroup() |>
  filter(years_after_intro < 10) |>
  mutate(t_eff_yearly = (1 - inc_ratio_0) * 100,
         t_eff_cumul = (1 - cuminc_ratio_0) * 100)

traj_grid_midintro |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = inc_ratio_0)) +
  geom_tile(data = traj_grid_midintro |> filter((1-frac)*100 == 0),
              aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "red", width = 10, lwd = 1) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_viridis_c(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Ratio of mosquitoes to hosts",
       fill = "Incidence relative\nto 0% wMel") +
  theme_bw()

traj_grid_midintro |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = inc_ratio_0)) +
  geom_tile(data = traj_grid_midintro |> filter((1-frac)*100 == 0),
              aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "red", width = 10, lwd = 1) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_viridis_c(trans = "log1p", breaks = c(0, 1, 100, 800)) +
    scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Ratio of mosquitoes to hosts",
       fill = "Incidence relative\nto 0% wMel") +
  theme_bw()

(q1 <- traj_grid_midintro |>
  filter(years_after_intro < 4) |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = inc_ratio_0)) +
  geom_tile(data = traj_grid_midintro |> filter((1-frac)*100 == 0 & years_after_intro < 4),
              aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "red", width = 10, lwd = 1) +
  geom_tile(data = traj_grid_midintro |> filter((1-frac)*100 == 60 & years_after_intro < 4 & M0 ==3),
            aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "white", width = 10, lwd = 1) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_viridis_c(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "Incidence relative\nto 0% wMel") +
  theme_bw())

# ggsave(plot = n1, "outputs/2026-03-05_fig5pt1.pdf", width = 7.5, height = 2.5)

traj_grid_midintro |>
  filter(years_after_intro < 4) |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = t_eff_yearly)) +
  geom_tile(data = traj_grid_midintro |> filter((1-frac)*100 == 0 & years_after_intro < 4),
              aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "red", width = 10, lwd = 1) +
  geom_tile(data = traj_grid_midintro |> filter((1-frac)*100 == 60 & years_after_intro < 4 & M0 ==3),
            aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "white", width = 10, lwd = 1) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_viridis_c(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
    scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "Yearly\nintervention\neffect\n(%)") +
  theme_bw()

traj_grid_midintro |>
  filter((1-frac)*100 == 60 & years_after_intro < 4 & M0 == 3) |>
  dplyr::select(years_after_intro, M0, frac, t_eff_yearly, inc_ratio_0)

(traj_grid_midintro |>
    ggplot() +
    geom_tile(aes(x = (1-frac)*100, y = M0, fill = cuminc_ratio_0)) +
  geom_tile(data = traj_grid_midintro |> filter((1-frac)*100 == 0),
              aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "red", width = 10, lwd = 1) +
    facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
    scale_fill_viridis_c(limits = c(0, NA), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    scale_y_continuous(breaks = c(1:10)) +
    labs(x = "% wMel", y = "Ratio of mosquitoes to hosts",
         fill = "Cumulative incidence\nrelative to 0% wMel") +
    theme_bw())

(q2 <- traj_grid_midintro |>
  filter(years_after_intro < 4) |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = cuminc_ratio_0)) +
  geom_tile(data = traj_grid_midintro |> filter((1-frac)*100 == 0 & years_after_intro < 4),
              aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "red", width = 10, lwd = 1) +
  geom_tile(data = traj_grid_midintro |> filter((1-frac)*100 == 60 & years_after_intro < 4 & M0 ==3),
            aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "white", width = 10, lwd = 1) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_viridis_c(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "Cumulative incidence\nrelative to 0% wMel") +
  theme_bw())

# ggsave(plot = n2, "outputs/2026-03-05_fig5pt2.pdf", width = 7.5, height = 2.5)

traj_grid_midintro |>
  filter(years_after_intro < 4) |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = t_eff_cumul)) +
  geom_tile(data = traj_grid_midintro |> filter((1-frac)*100 == 0 & years_after_intro < 4),
              aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "red", width = 10, lwd = 1) +
  geom_tile(data = traj_grid_midintro |> filter((1-frac)*100 == 60 & years_after_intro < 4 & M0 ==3),
            aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "white", width = 10, lwd = 1) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_viridis_c(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
    scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "Cumulative\nintervention\neffect\n(%)") +
  theme_bw()

traj_grid_midintro |>
  filter((1-frac)*100 == 60 & years_after_intro < 4 & M0 == 3) |>
  dplyr::select(years_after_intro, M0, frac, t_eff_cumul, cuminc_ratio_0)

# -- Tileplots comparing to 0% wMel in M0=8
traj_grid_biased_midintro <- trajectories_midintro |>
  ungroup() |>
  mutate(year = floor((time - 1) / 365),
         years_after_intro = year+1) |>
  group_by(M0, Mnew, frac, year, years_after_intro) |>
  summarise(Ih = sum(Ih),
            Rhc = sum(Rhc),
            .groups = "drop") |>
  arrange(M0, Mnew, frac, year, years_after_intro) |>
  group_by(M0, Mnew, frac) |>
  mutate(Rhc_post = cumsum(if_else(years_after_intro >= 0, Ih, 0))) |>
  ungroup() |>
  filter(year %in% c(0:10)) |>
  group_by(year) |>
  arrange(desc(M0)) |>
  mutate(Ih_intro0 = Ih[30],
         Rhc_intro0 = Rhc_post[30]) |>
  ungroup() |>
  rowwise() |>
  mutate(inc_ratio_0 = Ih/Ih_intro0,
         cuminc_ratio_0 = Rhc_post/Rhc_intro0) |>
  ungroup() |>
  filter(years_after_intro < 10) |>
  mutate(t_eff_yearly = (1 - inc_ratio_0) * 100,
         t_eff_cumul = (1 - cuminc_ratio_0) * 100)

traj_grid_biased_midintro |>
    ggplot() +
    geom_tile(aes(x = (1-frac)*100, y = M0, fill = inc_ratio_0)) +
    geom_tile(data = traj_grid_biased_midintro |> filter((1-frac)*100 == 0),
              aes(x = (1-frac)*100, y = 8),
              fill = NA, color = "red", width = 10, lwd = 1) +
    facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
    scale_fill_viridis_c(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    scale_y_continuous(breaks = c(1:10)) +
    labs(x = "% wMel", y = "Ratio of mosquitoes to hosts",
         fill = "Incidence relative\nto 0% wMel for M=8") +
    theme_bw()

traj_grid_biased_midintro |>
    ggplot() +
    geom_tile(aes(x = (1-frac)*100, y = M0, fill = inc_ratio_0)) +
    geom_tile(data = traj_grid_biased_midintro |> filter((1-frac)*100 == 0),
              aes(x = (1-frac)*100, y = 8),
              fill = NA, color = "red", width = 10, lwd = 1) +
    facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
    scale_fill_viridis_c(trans = "log1p", breaks = c(0, 1, 100, 800)) +
    scale_y_continuous(breaks = c(1:10)) +
    labs(x = "% wMel", y = "Ratio of mosquitoes to hosts",
         fill = "Incidence relative\nto 0% wMel for M=8") +
    theme_bw()

traj_grid_biased_midintro |>
    ggplot() +
    geom_tile(aes(x = (1-frac)*100, y = M0, fill = cuminc_ratio_0)) +
    geom_tile(data = traj_grid_biased_midintro |> filter((1-frac)*100 == 0),
              aes(x = (1-frac)*100, y = 8),
              fill = NA, color = "red", width = 10, lwd = 1) +
    facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
    scale_fill_viridis_c(limits = c(0, NA), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    scale_y_continuous(breaks = c(1:10)) +
    labs(x = "% wMel", y = "Ratio of mosquitoes to hosts",
         fill = "Cumulative incidence\nrelative to 0% wMel for M=8") +
    theme_bw()

# -- Tileplots comparing to 0% wMel in M0=6
traj_grid_biased_midintro2 <- trajectories_midintro |>
  ungroup() |>
  mutate(year = floor((time - 1) / 365),
         years_after_intro = year+1) |>
  group_by(M0, Mnew, frac, year, years_after_intro) |>
  summarise(Ih = sum(Ih),
            Rhc = sum(Rhc),
            .groups = "drop") |>
  arrange(M0, Mnew, frac, year, years_after_intro) |>
  group_by(M0, Mnew, frac) |>
  mutate(Rhc_post = cumsum(if_else(years_after_intro >= 0, Ih, 0))) |>
  ungroup() |>
  filter(year %in% c(0:10)) |>
  group_by(year) |>
  arrange(desc(M0)) |>
  mutate(Ih_intro0 = Ih[50],
         Rhc_intro0 = Rhc_post[50]) |>
  ungroup() |>
  rowwise() |>
  mutate(inc_ratio_0 = Ih/Ih_intro0,
         cuminc_ratio_0 = Rhc_post/Rhc_intro0) |>
  ungroup() |>
  filter(years_after_intro < 10) |>
  mutate(t_eff_yearly = (1 - inc_ratio_0) * 100,
         t_eff_cumul = (1 - cuminc_ratio_0) * 100)

traj_grid_biased_midintro2 |>
    ggplot() +
    geom_tile(aes(x = (1-frac)*100, y = M0, fill = inc_ratio_0)) +
    geom_tile(data = traj_grid_biased_midintro2 |> filter((1-frac)*100 == 0),
              aes(x = (1-frac)*100, y = 6),
              fill = NA, color = "red", width = 10, lwd = 1) +
    facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
    scale_fill_viridis_c(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    scale_y_continuous(breaks = c(1:10)) +
    labs(x = "% wMel", y = "Ratio of mosquitoes to hosts",
         fill = "Incidence relative\nto 0% wMel for M=6") +
    theme_bw()

traj_grid_biased_midintro2 |>
    ggplot() +
    geom_tile(aes(x = (1-frac)*100, y = M0, fill = inc_ratio_0)) +
    geom_tile(data = traj_grid_biased_midintro2 |> filter((1-frac)*100 == 0),
              aes(x = (1-frac)*100, y = 6),
              fill = NA, color = "red", width = 10, lwd = 1) +
    facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
    scale_fill_viridis_c(trans = "log1p", breaks = c(0, 1, 100, 800)) +
    scale_y_continuous(breaks = c(1:10)) +
    labs(x = "% wMel", y = "Ratio of mosquitoes to hosts",
         fill = "Incidence relative\nto 0% wMel for M=6") +
    theme_bw()

(q3 <- traj_grid_biased_midintro2 |>
  filter(years_after_intro < 4) |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = inc_ratio_0)) +
  geom_tile(data = traj_grid_biased_midintro2 |> filter((1-frac)*100 == 0 & years_after_intro < 4),
              aes(x = (1-frac)*100, y = 6),
              fill = NA, color = "red", width = 10, lwd = 1) +
  geom_tile(data = traj_grid_biased_midintro2 |> filter((1-frac)*100 == 60 & years_after_intro < 4 & M0 ==3),
            aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "white", width = 10, lwd = 1) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_viridis_c(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "Incidence relative\nto 0% wMel") +
  theme_bw())

# ggsave(plot = n3, "outputs/2026-03-05_fig5pt3.pdf", width = 7.5, height = 2.5)

traj_grid_biased_midintro2 |>
  filter(years_after_intro < 4) |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = t_eff_yearly)) +
  geom_tile(data = traj_grid_biased_midintro2 |> filter((1-frac)*100 == 0 & years_after_intro < 4),
              aes(x = (1-frac)*100, y = 6),
              fill = NA, color = "red", width = 10, lwd = 1) +
  geom_tile(data = traj_grid_biased_midintro2 |> filter((1-frac)*100 == 60 & years_after_intro < 4 & M0 ==3),
            aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "white", width = 10, lwd = 1) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_viridis_c(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
    scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "Yearly\nintervention\neffect\n(%)") +
  theme_bw()

traj_grid_biased_midintro2 |>
  filter((1-frac)*100 == 60 & years_after_intro < 4 & M0 == 3) |>
  dplyr::select(years_after_intro, M0, frac, t_eff_yearly, inc_ratio_0)

traj_grid_biased_midintro2 |>
    ggplot() +
    geom_tile(aes(x = (1-frac)*100, y = M0, fill = cuminc_ratio_0)) +
    geom_tile(data = traj_grid_biased_midintro2 |> filter((1-frac)*100 == 0),
              aes(x = (1-frac)*100, y = 6),
              fill = NA, color = "red", width = 10, lwd = 1) +
    facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
    scale_fill_viridis_c(limits = c(0, NA), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    scale_y_continuous(breaks = c(1:10)) +
    labs(x = "% wMel", y = "Ratio of mosquitoes to hosts",
         fill = "Cumulative incidence\nrelative to 0% wMel for M=6") +
    theme_bw()

(q4 <- traj_grid_biased_midintro2|>
  filter(years_after_intro < 4) |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = cuminc_ratio_0)) +
  geom_tile(data = traj_grid_biased_midintro2 |> filter((1-frac)*100 == 0 & years_after_intro < 4),
              aes(x = (1-frac)*100, y = 6),
              fill = NA, color = "red", width = 10, lwd = 1) +
  geom_tile(data = traj_grid_biased_midintro2|> filter((1-frac)*100 == 60 & years_after_intro < 4 & M0 ==3),
            aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "white", width = 10, lwd = 1) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_viridis_c(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "Cumulative incidence\nrelative to 0% wMel") +
  theme_bw())

# ggsave(plot = n4, "outputs/2026-03-05_fig5pt4.pdf", width = 7.5, height = 2.5)

traj_grid_biased_midintro2 |>
  filter(years_after_intro < 4) |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = t_eff_cumul)) +
  geom_tile(data = traj_grid_biased_midintro2 |> filter((1-frac)*100 == 0 & years_after_intro < 4),
              aes(x = (1-frac)*100, y = 6),
              fill = NA, color = "red", width = 10, lwd = 1) +
  geom_tile(data = traj_grid_biased_midintro2 |> filter((1-frac)*100 == 60 & years_after_intro < 4 & M0 ==3),
            aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "white", width = 10, lwd = 1) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_viridis_c(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
    scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "Cumulative\nintervention\neffect\n(%)") +
  theme_bw()

traj_grid_biased_midintro2 |>
  filter((1-frac)*100 == 60 & years_after_intro < 4 & M0 == 3) |>
  dplyr::select(years_after_intro, M0, frac, t_eff_cumul, cuminc_ratio_0)

# -- Tileplots comparing true to biased (M0 = 6) counterfactual
comp_midintro <- traj_grid_midintro |>
  mutate(wmel = (1-frac)*100,
         true_incratio = inc_ratio_0,
         true_cumincratio = cuminc_ratio_0) |>
  dplyr::select(M0, wmel, true_incratio, true_cumincratio, years_after_intro) |>
  left_join(traj_grid_biased_midintro2 |>
              mutate(wmel = (1-frac)*100,
                     biased_incratio = inc_ratio_0,
                     biased_cumincratio = cuminc_ratio_0) |>
              dplyr::select(M0, wmel, biased_incratio, biased_cumincratio, years_after_intro),
            by = c("M0", "wmel", "years_after_intro")) |>
  ungroup() |>
  mutate(diff_incratio = true_incratio - biased_incratio,
         diff_cumincratio = true_cumincratio - biased_cumincratio,
         ratio_incratio = true_incratio / biased_incratio,
         ratio_cumincratio = true_cumincratio / biased_cumincratio)

## yearly
comp_midintro |>
  ggplot() +
  geom_tile(aes(x = wmel, y = M0, fill = diff_incratio)) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_distiller(palette = "Spectral",
                       limits = c(-max(abs(comp_midintro$diff_incratio), na.rm = TRUE),
                                  max(abs(comp_midintro$diff_incratio), na.rm = TRUE))) +
  scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Ratio of mosquitoes to hosts",
       fill = "True - biased yearly incidence") +
  theme_bw()

comp_midintro |>
  mutate(diff_incratio = cut(diff_incratio, breaks = c(-650, -10, -1, -0.5, 0.5, 1, 10, 650))) |>
  ggplot() +
  geom_tile(aes(x = wmel, y = M0, fill = diff_incratio)) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_brewer(palette = "Spectral") +
  scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Ratio of mosquitoes to hosts",
       fill = "True - biased yearly incidence") +
  theme_bw()

comp_midintro |>
  mutate(diff_incratio = factor(case_when(
    diff_incratio > 0 ~ "True incidence higher",
    diff_incratio == 0 ~ "No bias",
    diff_incratio < 0 ~ "True incidence lower"),
    levels = c("True incidence higher", "No bias", "True incidence lower"))) |>
  ggplot() +
  geom_tile(aes(x = wmel, y = M0, fill = diff_incratio)) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_brewer(palette = "Spectral") +
  scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Ratio of mosquitoes to hosts",
       fill = "True - biased yearly incidence") +
  theme_bw()

## cumulative
comp_midintro |>
  ggplot() +
  geom_tile(aes(x = wmel, y = M0, fill = diff_cumincratio)) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_distiller(palette = "Spectral",
                       limits = c(-max(abs(comp_midintro$diff_cumincratio), na.rm = TRUE),
                                  max(abs(comp_midintro$diff_cumincratio), na.rm = TRUE))) +
  scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Ratio of mosquitoes to hosts",
       fill = "True - biased cumulative incidence") +
  theme_bw()

comp_midintro |>
  mutate(diff_cumincratio = cut(
    diff_cumincratio,
    breaks = c(-0.2, -0.1, -0.05, -0.01, 0.01, 0.05, 0.1, 0.2))) |>
  ggplot() +
  geom_tile(aes(x = wmel, y = M0, fill = diff_cumincratio)) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_brewer(palette = "Spectral", drop = FALSE) +
  scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Ratio of mosquitoes to hosts",
       fill = "True - biased cumulative incidence") +
  theme_bw()

comp_midintro |>
  mutate(diff_cumincratio = factor(case_when(
    diff_cumincratio > 0 ~ "True incidence higher",
    diff_cumincratio == 0 ~ "No bias",
    diff_cumincratio < 0 ~ "True incidence lower"),
    levels = c("True incidence higher", "No bias", "True incidence lower"))) |>
  ggplot() +
  geom_tile(aes(x = wmel, y = M0, fill = diff_cumincratio)) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_brewer(palette = "Spectral") +
  scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Ratio of mosquitoes to hosts",
       fill = "True - biased cumulative incidence") +
  theme_bw()

# ---- Simulations with late introduction of wMel (near equilibrium susceptibility) ----

# -- Define model parameters
pars = c(N_h = 1e6,              # Host population size.
         S_h = 0.999,            # Starting proportion of susceptible hosts.
         I_h = 0.001,            # Starting proportion of infected hosts.
         mu_h = 1 / (77 * 365),  # Host mortality rate.
         r_h = 1 / (77 * 365),   # Host birth rate.
         gamma_h = 1 / 5,        # Host recovery from infectiousness.
         delta_h = 1 / 5.9,      # Intrinsic incubation period.
         # Constant rate of introduction of mosquitoes, as a proportion of `N_h`.
         trickle = 1e-5)

# -- State variables
make_xstart <- function(pars_vec) {
  x0 <- setNames(numeric(8), c("Sh","Eh","Ih","Rh","Rhc","Sv","Ev","Iv"))
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

# -- Length of simulation and time of introduction of wMel
n_days_total <- 365 * 41
t_break      <- 365 * 20
times1 <- 1:(t_break)
times2 <- (t_break + 1):(n_days_total + 1)

# -- Define Ms and wMel % to test
M_grid <- seq(1, 10, by = 1)
reduct_fracs <- seq(0.1, 1, by = 0.1)

results_grid <- expand.grid(
  M0    = M_grid,
  frac  = reduct_fracs,
  KEEP.OUT.ATTRS = FALSE
)

# -- Define temperature function
# Option (i): Arbitrary made-up example time series of temperatures.

#' Get a temperature time series based on a sine curve.
#'
#' @param x :numeric (vector): time point(s) for which to get the temperature.
#' @param mean_temp :numeric: mean temperature in C.
#' @param ampl_temp :numeric: peak-to-peak amplitude of temperature.
GetTemp = function(x, mean_temp = 26, ampl_temp = 3)
{
  tmp = ampl_temp * sin(2 * pi * x / 365) + mean_temp
  names(tmp) = NULL
  return(tmp)
}

# -- Run model
trajectories_lateintro <- read_rds("data/interim/seirsei_fig5-late intro suscept.rds")

# -- Plot outputs
year_introduction <- t_break/365
years_incratio <- c((year_introduction):(n_days_total/365))

# -- Time series of infections
trajectories_lateintro |>
  ggplot(aes(x = time/365, y = Ih, color = factor(100-frac*100))) +
  geom_line() +
  geom_vline(xintercept = (t_break + 1)/365, linetype = 2) +
  labs(x = "Year", y = "Infections",
       color = "%wMel at intervention") +
  scale_x_continuous(breaks = c(1:10)) +
  facet_wrap(~M0, labeller = labeller(M0 = function(x) paste0("M0 = ", x)), ncol = 5) +
  theme_bw() +
  coord_cartesian(ylim = c(0, 250)) +
  theme(legend.position = "bottom",
        legend.byrow = TRUE)

trajectories_lateintro |>
  filter((100-frac*100) %in% c(0, 50, 90)) |>
  ggplot(aes(x = time/365, y = Ih, color = factor(100-frac*100))) +
  geom_line() +
  geom_vline(xintercept = (t_break + 1)/365, linetype = 2) +
  labs(x = "Year", y = "Infections",
       color = "%wMel at intervention") +
  facet_wrap(~M0, labeller = labeller(M0 = function(x) paste0("M0 = ", x)), ncol = 5) +
  theme_bw() +
  coord_cartesian(ylim = c(0, 250)) +
  theme(legend.position = "bottom",
        legend.byrow = TRUE)

trajectories_lateintro |>
  filter((100-frac*100) %in% c(0, 50, 90)) |>
  ggplot(aes(x = time/365, y = Ih, color = factor(100-frac*100))) +
  geom_line() +
  geom_vline(xintercept = (t_break + 1)/365, linetype = 2) +
  labs(x = "Year", y = "Infections",
       color = "%wMel at intervention") +
  facet_wrap(~M0, scales = "free_y", labeller = labeller(M0 = function(x) paste0("M0 = ", x)), ncol = 5) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.byrow = TRUE)

trajectories_lateintro |>
  ggplot(aes(x = time/365, y = Ih, color = factor(100-frac*100))) +
  geom_line() +
  geom_vline(xintercept = (t_break + 1)/365, linetype = 2) +
  labs(x = "Year", y = "Infections",
       color = "%wMel at intervention") +
  facet_wrap(~M0, scales = "free_y", labeller = labeller(M0 = function(x) paste0("M0 = ", x)), ncol = 5) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.byrow = TRUE)

trajectories_lateintro |>
  filter(M0 %in% c(1, 6) & (100-frac*100) %in% c(0, 50, 90)) |>
  ggplot(aes(x = time/365, y = Ih, color = factor(100-frac*100))) +
  geom_line(lwd = 1) +
  geom_vline(xintercept = (t_break + 1)/365, linetype = 2) +
  scale_color_brewer(palette = "Set1") +
  labs(x = "Year", y = "Infections",
       color = "%wMel at intervention") +
  facet_wrap(~M0, scales = "free_y", labeller = labeller(M0 = function(x) paste0("M0 = ", x)), ncol = 5) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.byrow = TRUE)

trajectories_lateintro |>
  filter(time/365 > year_introduction) |>
  filter((100-frac*100) %in% c(0, 50, 90)) |>
  ggplot(aes(x = time/365, y = Ih, color = factor(100-frac*100))) +
  geom_line() +
  geom_vline(xintercept = (t_break + 1)/365, linetype = 2) +
  labs(x = "Year", y = "Infections",
       color = "%wMel at intervention") +
  facet_wrap(~M0, scales = "free_y", labeller = labeller(M0 = function(x) paste0("M0 = ", x)), ncol = 5) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.byrow = TRUE)

trajectories_lateintro |>
  filter(time/365 > year_introduction) |>
  ggplot(aes(x = time/365, y = Ih, color = factor(100-frac*100))) +
  geom_line() +
  geom_vline(xintercept = (t_break + 1)/365, linetype = 2) +
  labs(x = "Year", y = "Infections",
       color = "%wMel at intervention") +
  facet_wrap(~M0, scales = "free_y", labeller = labeller(M0 = function(x) paste0("M0 = ", x)), ncol = 5) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.byrow = TRUE)

trajectories_lateintro |>
  filter(time/365 > year_introduction) |>
  filter(M0 %in% c(1, 6) & (100-frac*100) %in% c(0, 50, 90)) |>
  ggplot(aes(x = time/365, y = Ih, color = factor(100-frac*100))) +
  geom_line(lwd = 0.8) +
  geom_vline(xintercept = (t_break + 1)/365, linetype = 2) +
  scale_color_brewer(palette = "Set1") +
  labs(x = "Year", y = "Infections",
       color = "%wMel at intervention") +
  facet_wrap(~M0, scales = "free_y", labeller = labeller(M0 = function(x) paste0("M0 = ", x)), ncol = 5) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.byrow = TRUE)

# -- Yearly and cumulative incidence
traj_rio_lateintro <- trajectories_lateintro |>
  ungroup() |>
  mutate(year = floor((time - 1) / 365),
         years_after_intro = year - year_introduction) |>
  group_by(M0, Mnew, frac, year, years_after_intro) |>
  summarise(Ih = sum(Ih),
            Rhc = sum(Rhc),
            .groups = "drop") |>
  arrange(M0, Mnew, frac, year, years_after_intro) |>
  group_by(M0, Mnew, frac) |>
  mutate(Rhc_post = cumsum(if_else(years_after_intro >= 0, Ih, 0))) |>
  ungroup() |>
  filter(year %in% years_incratio) |>
  group_by(M0, year) |>
  arrange(desc(Mnew)) |>
  mutate(Ih_intro0 = Ih[1],
         Rhc_intro0 = Rhc_post[1],
         Ih_intro10 = Ih[2],
         Rhc_intro10 = Rhc_post[2]) |>
  ungroup() |>
  rowwise() |>
  mutate(inc_ratio_0 = Ih/Ih_intro0,
         inc_ratio_10 = Ih/Ih_intro10,
         cuminc_ratio_0 = Rhc_post/Rhc_intro0,
         cuminc_ratio_10 = Rhc_post/Rhc_intro10) |>
  ungroup() |>
  filter(years_after_intro <= 10)

spectral_n <- colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(length(unique(traj_rio_lateintro$years_after_intro)))

# -- Yearly incidence, facets by M0, wMel in x axis
traj_rio_lateintro |>
  ggplot(aes(x = (1-frac)*100, y = inc_ratio_0, color = factor(years_after_intro))) +
  geom_point() +
  geom_line() +
  geom_hline(aes(yintercept = 1)) +
  facet_wrap(~M0, labeller = labeller(M0 = function(x) paste0("M0 = ", x)), ncol = 5) +
  scale_color_manual(values = spectral_n) +
  scale_x_continuous(breaks = seq(0, 90, 10)) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  labs(x = "%wMel", y = "Yearly incidence ratio\n(ref: 0% wMel)",
       color = "Years after wMel\nintroduction") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.byrow = TRUE) +
  guides(color = guide_legend(nrow = 2)) +
  coord_cartesian(ylim = c(0, 1))

# -- Cumulative incidence, facets by M0, wMel in x axis
traj_rio_lateintro |>
  ggplot(aes(x = (1-frac)*100, y = cuminc_ratio_0, color = factor(years_after_intro))) +
  geom_point() +
  geom_line() +
  geom_hline(aes(yintercept = 1)) +
  facet_wrap(~M0, labeller = labeller(M0 = function(x) paste0("M0 = ", x)), ncol = 5) +
  scale_color_manual(values = spectral_n) +
  scale_x_continuous(breaks = seq(0, 90, 10)) +
  labs(x = "%wMel", y = "Cumulative incidence ratio\n(ref: 0% wMel)",
       color = "Years after wMel\nintroduction") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.byrow = TRUE) +
  guides(color = guide_legend(nrow = 2)) +
  coord_cartesian(ylim = c(0, 1))

# -- Yearly incidence, facets by M0, years in x axis
spectral_2 <- colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(length(unique(traj_rio_lateintro$frac)))

traj_rio_lateintro |>
  ggplot(aes(x = years_after_intro, y = inc_ratio_0, color = factor((1-frac)*100))) +
  geom_point() +
  geom_line() +
  geom_hline(aes(yintercept = 1)) +
  facet_wrap(~M0, labeller = labeller(M0 = function(x) paste0("M0 = ", x)), ncol = 5) +
  scale_color_manual(values = spectral_2) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_x_continuous(breaks = c(1:10)) +
  labs(x = "Years after wMel introduction", y = "Yearly incidence ratio\n(ref: 0% wMel)",
       color = "% wMel") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.byrow = TRUE) +
  coord_cartesian(ylim = c(0, 1))

# -- Cumulative incidence, facets by M0, years in x axis
traj_rio_lateintro |>
  ggplot(aes(x = years_after_intro, y = cuminc_ratio_0, color = factor((1-frac)*100))) +
  geom_point() +
  geom_line() +
  geom_hline(aes(yintercept = 1)) +
  facet_wrap(~M0, labeller = labeller(M0 = function(x) paste0("M0 = ", x)), ncol = 5) +
  scale_color_manual(values = spectral_2) +
  scale_x_continuous(breaks = c(1:10)) +
  labs(x = "Years after wMel introduction", y = "Cumulative incidence ratio\n(ref: 0% wMel)",
       color = "% wMel") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.byrow = TRUE) +
  coord_cartesian(ylim = c(0, 1))

# -- Only plots for M=1 and M=6
(fig5_leg <- traj_rio_lateintro |>
  filter(M0 == 1) |>
  ggplot(aes(x = years_after_intro, y = inc_ratio_0, color = factor((1-frac)*100))) +
  geom_point() +
  geom_line() +
  geom_hline(aes(yintercept = 1)) +
  scale_color_manual(values = spectral_2) +
  scale_x_continuous(breaks = c(1:10)) +
  # scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Years after wMel introduction", y = "Yearly incidence ratio\n(ref: 0% wMel)",
       color = "% wMel") +
  theme_bw() +
  theme(
    legend.byrow = TRUE) +
  coord_cartesian(ylim = c(0, 1)))

traj_rio_lateintro |>
  filter(M0 == 1) |>
  ggplot(aes(x = years_after_intro, y = inc_ratio_0, color = factor((1-frac)*100))) +
  geom_point() +
  geom_line() +
  geom_hline(aes(yintercept = 1)) +
  scale_color_manual(values = spectral_2) +
  scale_x_continuous(breaks = c(1:10)) +
  labs(x = "Years after wMel introduction", y = "Yearly incidence ratio\n(ref: 0% wMel)",
       color = "% wMel") +
  theme_bw() +
  theme(
    legend.position = "none",
    legend.byrow = TRUE) +
  coord_cartesian(ylim = c(0, 1))

traj_rio_lateintro |>
  filter(M0 == 1 & frac %in% c(0.1, 0.5, 1)) |>
  ggplot(aes(x = years_after_intro, y = inc_ratio_0, color = factor((1-frac)*100))) +
  geom_point() +
  geom_line() +
  geom_hline(aes(yintercept = 1)) +
  scale_color_brewer(palette = "Set1") +
  scale_x_continuous(breaks = c(1:10)) +
  labs(x = "Years after wMel introduction", y = "Yearly incidence ratio\n(ref: 0% wMel)",
       color = "% wMel") +
  theme_bw() +
  theme(
    legend.byrow = TRUE) +
  coord_cartesian(ylim = c(0, 1))

# ggsave("outputs/fig5_1-feb16.pdf", plot = fig5_1, width = 3, height = 2)

traj_rio_lateintro |>
  filter(M0 == 6) |>
  ggplot(aes(x = years_after_intro, y = inc_ratio_0, color = factor((1-frac)*100))) +
  geom_point() +
  geom_line() +
  geom_hline(aes(yintercept = 1)) +
  scale_color_manual(values = spectral_2) +
  scale_x_continuous(breaks = c(1:10)) +
  scale_y_continuous(
    ) +
  labs(x = "Years after wMel introduction", y = "Yearly incidence ratio\n(ref: 0% wMel)",
       color = "% wMel") +
  theme_bw() +
  theme(
    legend.position = "none",
    legend.byrow = TRUE) +
  coord_cartesian(ylim = c(0, 1))

traj_rio_lateintro |>
  filter(M0 == 6 & frac %in% c(0.1, 0.5, 1)) |>
  ggplot(aes(x = years_after_intro, y = inc_ratio_0, color = factor((1-frac)*100))) +
  geom_point() +
  geom_line() +
  geom_hline(aes(yintercept = 1)) +
  scale_color_brewer(palette = "Set1") +
  scale_x_continuous(breaks = c(1:10)) +
  scale_y_continuous(
    # limits = c(0, 1),
    ) +
  labs(x = "Years after wMel introduction", y = "Yearly incidence ratio\n(ref: 0% wMel)",
       color = "% wMel") +
  theme_bw() +
  theme(
    # legend.position = "none",
    legend.byrow = TRUE) +
  coord_cartesian(ylim = c(0, 1))

# ggsave("outputs/fig5_2-feb16.pdf", plot = fig5_2, width = 3, height = 2)

(fig5_3 <- traj_rio_lateintro |>
  filter(M0 == 1) |>
  ggplot(aes(x = years_after_intro, y = cuminc_ratio_0, color = factor((1-frac)*100))) +
  geom_point() +
  geom_line() +
  geom_hline(aes(yintercept = 1)) +
  scale_color_manual(values = spectral_2) +
  scale_x_continuous(breaks = c(1:10)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Years after wMel introduction", y = "Cumulative incidence ratio\n(ref: 0% wMel)",
       color = "% wMel") +
  theme_bw() +
  theme(legend.position = "none",
        legend.byrow = TRUE))

traj_rio_lateintro |>
  filter(M0 == 1 & frac %in% c(0.1, 0.5, 1)) |>
  ggplot(aes(x = years_after_intro, y = cuminc_ratio_0, color = factor((1-frac)*100))) +
  geom_point() +
  geom_line() +
  geom_hline(aes(yintercept = 1)) +
  scale_color_brewer(palette = "Set1") +
  scale_x_continuous(breaks = c(1:10)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Years after wMel introduction", y = "Cumulative incidence ratio\n(ref: 0% wMel)",
       color = "% wMel") +
  theme_bw() +
  theme(
    # legend.position = "none",
      legend.byrow = TRUE)

# ggsave("outputs/fig5_3-feb16.pdf", plot = fig5_3, width = 3, height = 2)

(fig5_4 <- traj_rio_lateintro |>
  filter(M0 == 6) |>
  ggplot(aes(x = years_after_intro, y = cuminc_ratio_0, color = factor((1-frac)*100))) +
  geom_point() +
  geom_line() +
  geom_hline(aes(yintercept = 1)) +
  scale_color_manual(values = spectral_2) +
  scale_x_continuous(breaks = c(1:10)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Years after wMel introduction", y = "Cumulative incidence ratio\n(ref: 0% wMel)",
       color = "% wMel") +
  theme_bw() +
  theme(legend.position = "none",
        legend.byrow = TRUE))

traj_rio_lateintro |>
  filter(M0 == 6 & frac %in% c(0.1, 0.5, 1)) |>
  ggplot(aes(x = years_after_intro, y = cuminc_ratio_0, color = factor((1-frac)*100))) +
  geom_point() +
  geom_line() +
  geom_hline(aes(yintercept = 1)) +
  scale_color_brewer(palette = "Set1") +
  scale_x_continuous(breaks = c(1:10)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Years after wMel introduction", y = "Cumulative incidence ratio\n(ref: 0% wMel)",
       color = "% wMel") +
  theme_bw() +
  theme(
    # legend.position = "none",
      legend.byrow = TRUE)

# ggsave("outputs/fig5_4-feb16.pdf", plot = fig5_4, width = 3, height = 2)

# -- Tileplots comparing to true counterfactual
traj_grid_lateintro <- trajectories_lateintro |>
  ungroup() |>
  mutate(year = floor((time - 1) / 365),
         years_after_intro = (year+1) - year_introduction) |>
  group_by(M0, Mnew, frac, year, years_after_intro) |>
  summarise(Ih = sum(Ih),
            Rhc = sum(Rhc),
            .groups = "drop") |>
  arrange(M0, Mnew, frac, year, years_after_intro) |>
  group_by(M0, Mnew, frac) |>
  mutate(Rhc_post = cumsum(if_else(years_after_intro >= 1, Ih, 0))) |>
  ungroup() |>
  filter(year %in% years_incratio) |>
  group_by(M0, year) |>
  arrange(desc(Mnew)) |>
  mutate(Ih_intro0 = Ih[1],
         Rhc_intro0 = Rhc_post[1],
         Ih_intro10 = Ih[2],
         Rhc_intro10 = Rhc_post[2]) |>
  ungroup() |>
  rowwise() |>
  mutate(inc_ratio_0 = Ih/Ih_intro0,
         inc_ratio_10 = Ih/Ih_intro10,
         cuminc_ratio_0 = Rhc_post/Rhc_intro0,
         cuminc_ratio_10 = Rhc_post/Rhc_intro10) |>
  ungroup() |>
  filter(years_after_intro < 10) |>
  mutate(t_eff_yearly = (1 - inc_ratio_0) * 100,
         t_eff_cumul = (1 - cuminc_ratio_0) * 100)

traj_grid_lateintro |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = inc_ratio_0)) +
  geom_tile(data = traj_grid_lateintro |> filter((1-frac)*100 == 0),
              aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "red", width = 10, lwd = 1) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_viridis_c(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  # scale_fill_viridis_c() +
    scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Ratio of mosquitoes to hosts",
       fill = "Incidence relative\nto 0% wMel") +
  theme_bw()

traj_grid_lateintro |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = inc_ratio_0)) +
  geom_tile(data = traj_grid_lateintro |> filter((1-frac)*100 == 0),
              aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "red", width = 10, lwd = 1) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  # scale_fill_viridis_c(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_fill_viridis_c(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1, 1.2)) +
    scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Ratio of mosquitoes to hosts",
       fill = "Incidence relative\nto 0% wMel") +
  theme_bw()

(n5 <- traj_grid_lateintro |>
  filter(years_after_intro < 4) |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = inc_ratio_0)) +
  geom_tile(data = traj_grid_lateintro |> filter((1-frac)*100 == 0 & years_after_intro < 4),
              aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "red", width = 10, lwd = 1) +
  geom_tile(data = traj_grid_lateintro |> filter((1-frac)*100 == 60 & years_after_intro < 4 & M0 ==3),
            aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "white", width = 10, lwd = 1) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_viridis_c(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  # scale_fill_viridis_c() +
    scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "Incidence relative\nto 0% wMel") +
  theme_bw())

# ggsave(plot = n5, "outputs/2026-03-05_fig5pt5.pdf", width = 7.5, height = 2.5)

traj_grid_lateintro |>
  filter(years_after_intro < 4) |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = t_eff_yearly)) +
  geom_tile(data = traj_grid_lateintro |> filter((1-frac)*100 == 0 & years_after_intro < 4),
              aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "red", width = 10, lwd = 1) +
  geom_tile(data = traj_grid_lateintro |> filter((1-frac)*100 == 60 & years_after_intro < 4 & M0 ==3),
            aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "white", width = 10, lwd = 1) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_viridis_c(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
  # scale_fill_viridis_c() +
    scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "Yearly\nintervention\neffect\n(%)") +
  theme_bw()

traj_grid_lateintro |>
  filter((1-frac)*100 == 60 & years_after_intro < 4 & M0 == 3) |>
  dplyr::select(years_after_intro, M0, frac, t_eff_yearly, inc_ratio_0)

(fig_s4c <- traj_grid_lateintro |>
    ggplot() +
    geom_tile(aes(x = (1-frac)*100, y = M0, fill = cuminc_ratio_0)) +
  geom_tile(data = traj_grid_lateintro |> filter((1-frac)*100 == 0),
              aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "red", width = 10, lwd = 1) +
    facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
    scale_fill_viridis_c(limits = c(0, NA), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    scale_y_continuous(breaks = c(1:10)) +
    labs(x = "% wMel", y = "Ratio of mosquitoes to hosts",
         fill = "Cumulative incidence\nrelative to 0% wMel") +
    theme_bw())

(n6 <- traj_grid_lateintro |>
  filter(years_after_intro < 4) |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = cuminc_ratio_0)) +
  geom_tile(data = traj_grid_lateintro |> filter((1-frac)*100 == 0 & years_after_intro < 4),
              aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "red", width = 10, lwd = 1) +
  geom_tile(data = traj_grid_lateintro |> filter((1-frac)*100 == 60 & years_after_intro < 4 & M0 ==3),
            aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "white", width = 10, lwd = 1) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_viridis_c(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  # scale_fill_viridis_c() +
    scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "Cumulative incidence\nrelative to 0% wMel") +
  theme_bw())

# ggsave(plot = n6, "outputs/2026-03-05_fig5pt6.pdf", width = 7.5, height = 2.5)

traj_grid_lateintro |>
  filter(years_after_intro < 4) |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = t_eff_cumul)) +
  geom_tile(data = traj_grid_lateintro |> filter((1-frac)*100 == 0 & years_after_intro < 4),
              aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "red", width = 10, lwd = 1) +
  geom_tile(data = traj_grid_lateintro |> filter((1-frac)*100 == 60 & years_after_intro < 4 & M0 ==3),
            aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "white", width = 10, lwd = 1) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_viridis_c(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
  # scale_fill_viridis_c() +
    scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "Cumulative\nintervention\neffect\n(%)") +
  theme_bw()

traj_grid_lateintro |>
  filter((1-frac)*100 == 60 & years_after_intro < 4 & M0 == 3) |>
  dplyr::select(years_after_intro, M0, frac, t_eff_cumul, cuminc_ratio_0)

# -- Tileplots comparing to 0% wMel in M0=8
traj_grid_biased_lateintro <- trajectories_lateintro |>
  ungroup() |>
  mutate(year = floor((time - 1) / 365),
         years_after_intro = year - year_introduction) |>
  group_by(M0, Mnew, frac, year, years_after_intro) |>
  summarise(Ih = sum(Ih),
            Rhc = sum(Rhc),
            .groups = "drop") |>
  arrange(M0, Mnew, frac, year, years_after_intro) |>
  group_by(M0, Mnew, frac) |>
  mutate(Rhc_post = cumsum(if_else(years_after_intro >= 0, Ih, 0))) |>
  ungroup() |>
  filter(year %in% years_incratio) |>
  filter(years_after_intro < 10) |>
  group_by(year) |>
  arrange(desc(M0)) |>
  mutate(Ih_intro0 = Ih[30],
         Rhc_intro0 = Rhc_post[30]) |>
  ungroup() |>
  rowwise() |>
  mutate(inc_ratio_0 = Ih/Ih_intro0,
         cuminc_ratio_0 = Rhc_post/Rhc_intro0) |>
  ungroup() |>
  mutate(t_eff_yearly = (1 - inc_ratio_0) * 100,
         t_eff_cumul = (1 - cuminc_ratio_0) * 100)

traj_grid_biased_lateintro |>
    ggplot() +
    geom_tile(aes(x = (1-frac)*100, y = M0, fill = inc_ratio_0)) +
    geom_tile(data = traj_grid_biased_lateintro |> filter((1-frac)*100 == 0),
              aes(x = (1-frac)*100, y = 8),
              fill = NA, color = "red", width = 10, lwd = 1) +
    facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
    scale_fill_viridis_c(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    # scale_fill_viridis_c() +
    scale_y_continuous(breaks = c(1:10)) +
    labs(x = "% wMel", y = "Ratio of mosquitoes to hosts",
         fill = "Incidence relative\nto 0% wMel for M=8") +
    theme_bw()

traj_grid_biased_lateintro |>
    ggplot() +
    geom_tile(aes(x = (1-frac)*100, y = M0, fill = inc_ratio_0)) +
    geom_tile(data = traj_grid_biased_lateintro |> filter((1-frac)*100 == 0),
              aes(x = (1-frac)*100, y = 8),
              fill = NA, color = "red", width = 10, lwd = 1) +
    facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
    # scale_fill_viridis_c(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    scale_fill_viridis_c(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1, 1.2)) +
    scale_y_continuous(breaks = c(1:10)) +
    labs(x = "% wMel", y = "Ratio of mosquitoes to hosts",
         fill = "Incidence relative\nto 0% wMel for M=8") +
    theme_bw()

traj_grid_biased_lateintro |>
    ggplot() +
    geom_tile(aes(x = (1-frac)*100, y = M0, fill = cuminc_ratio_0)) +
    geom_tile(data = traj_grid_biased_lateintro |> filter((1-frac)*100 == 0),
              aes(x = (1-frac)*100, y = 8),
              fill = NA, color = "red", width = 10, lwd = 1) +
    facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
    scale_fill_viridis_c(limits = c(0, NA), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    scale_y_continuous(breaks = c(1:10)) +
    labs(x = "% wMel", y = "Ratio of mosquitoes to hosts",
         fill = "Cumulative incidence\nrelative to 0% wMel for M=8") +
    theme_bw()

# -- Tileplots comparing to 0% wMel in M0=6
traj_grid_biased_lateintro2 <- trajectories_lateintro |>
  ungroup() |>
  mutate(year = floor((time - 1) / 365),
         years_after_intro = year - year_introduction) |>
  group_by(M0, Mnew, frac, year, years_after_intro) |>
  summarise(Ih = sum(Ih),
            Rhc = sum(Rhc),
            .groups = "drop") |>
  arrange(M0, Mnew, frac, year, years_after_intro) |>
  group_by(M0, Mnew, frac) |>
  mutate(Rhc_post = cumsum(if_else(years_after_intro >= 0, Ih, 0))) |>
  ungroup() |>
  filter(year %in% years_incratio) |>
  filter(years_after_intro < 10) |>
  group_by(year) |>
  arrange(desc(M0)) |>
  mutate(Ih_intro0 = Ih[50],
         Rhc_intro0 = Rhc_post[50]) |>
  ungroup() |>
  rowwise() |>
  mutate(inc_ratio_0 = Ih/Ih_intro0,
         cuminc_ratio_0 = Rhc_post/Rhc_intro0) |>
  ungroup() |>
  mutate(t_eff_yearly = (1 - inc_ratio_0) * 100,
         t_eff_cumul = (1 - cuminc_ratio_0) * 100)

traj_grid_biased_lateintro2 |>
    ggplot() +
    geom_tile(aes(x = (1-frac)*100, y = M0, fill = inc_ratio_0)) +
    geom_tile(data = traj_grid_biased_lateintro2 |> filter((1-frac)*100 == 0),
              aes(x = (1-frac)*100, y = 6),
              fill = NA, color = "red", width = 10, lwd = 1) +
    facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
    scale_fill_viridis_c(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    # scale_fill_viridis_c() +
    scale_y_continuous(breaks = c(1:10)) +
    labs(x = "% wMel", y = "Ratio of mosquitoes to hosts",
         fill = "Incidence relative\nto 0% wMel for M=8") +
    theme_bw()

traj_grid_biased_lateintro2 |>
    ggplot() +
    geom_tile(aes(x = (1-frac)*100, y = M0, fill = inc_ratio_0)) +
    geom_tile(data = traj_grid_biased_lateintro2 |> filter((1-frac)*100 == 0),
              aes(x = (1-frac)*100, y = 6),
              fill = NA, color = "red", width = 10, lwd = 1) +
    facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
    # scale_fill_viridis_c(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    scale_fill_viridis_c(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1, 1.2)) +
    scale_y_continuous(breaks = c(1:10)) +
    labs(x = "% wMel", y = "Ratio of mosquitoes to hosts",
         fill = "Incidence relative\nto 0% wMel for M=8") +
    theme_bw()

(n7 <- traj_grid_biased_lateintro2 |>
  filter(years_after_intro < 4) |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = inc_ratio_0)) +
  geom_tile(data = traj_grid_biased_lateintro2 |> filter((1-frac)*100 == 0 & years_after_intro < 4),
              aes(x = (1-frac)*100, y = 6),
              fill = NA, color = "red", width = 10, lwd = 1) +
  geom_tile(data = traj_grid_biased_lateintro2 |> filter((1-frac)*100 == 60 & years_after_intro < 4 & M0 ==3),
            aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "white", width = 10, lwd = 1) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_viridis_c(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  # scale_fill_viridis_c() +
    scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "Incidence relative\nto 0% wMel") +
  theme_bw())

# ggsave(plot = n7, "outputs/2026-03-05_fig5pt7.pdf", width = 7.5, height = 2.5)

traj_grid_biased_lateintro2 |>
  filter(years_after_intro < 4) |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = t_eff_yearly)) +
  geom_tile(data = traj_grid_biased_lateintro2 |> filter((1-frac)*100 == 0 & years_after_intro < 4),
              aes(x = (1-frac)*100, y = 6),
              fill = NA, color = "red", width = 10, lwd = 1) +
  geom_tile(data = traj_grid_biased_lateintro2 |> filter((1-frac)*100 == 60 & years_after_intro < 4 & M0 ==3),
            aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "white", width = 10, lwd = 1) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_viridis_c(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
  # scale_fill_viridis_c() +
    scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "Yearly\nintervention\neffect\n(%)") +
  theme_bw()

traj_grid_biased_lateintro2 |>
  filter((1-frac)*100 == 60 & years_after_intro < 4 & M0 == 3) |>
  dplyr::select(years_after_intro, M0, frac, t_eff_yearly, inc_ratio_0)

traj_grid_biased_lateintro2 |>
    ggplot() +
    geom_tile(aes(x = (1-frac)*100, y = M0, fill = cuminc_ratio_0)) +
    geom_tile(data = traj_grid_biased_lateintro2 |> filter((1-frac)*100 == 0),
              aes(x = (1-frac)*100, y = 6),
              fill = NA, color = "red", width = 10, lwd = 1) +
    facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
    scale_fill_viridis_c(limits = c(0, NA), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    scale_y_continuous(breaks = c(1:10)) +
    labs(x = "% wMel", y = "Ratio of mosquitoes to hosts",
         fill = "Cumulative incidence\nrelative to 0% wMel for M=8") +
    theme_bw()

(n8 <- traj_grid_biased_lateintro2|>
  filter(years_after_intro < 4) |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = cuminc_ratio_0)) +
  geom_tile(data = traj_grid_biased_lateintro2 |> filter((1-frac)*100 == 0 & years_after_intro < 4),
              aes(x = (1-frac)*100, y = 6),
              fill = NA, color = "red", width = 10, lwd = 1) +
  geom_tile(data = traj_grid_biased_lateintro2|> filter((1-frac)*100 == 60 & years_after_intro < 4 & M0 ==3),
            aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "white", width = 10, lwd = 1) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_viridis_c(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  # scale_fill_viridis_c() +
    scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "Cumulative incidence\nrelative to 0% wMel") +
  theme_bw())

# ggsave(plot = n8, "outputs/2026-03-05_fig5pt8.pdf", width = 7.5, height = 2.5)

traj_grid_biased_lateintro2 |>
  filter(years_after_intro < 4) |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = t_eff_cumul)) +
  geom_tile(data = traj_grid_biased_lateintro2 |> filter((1-frac)*100 == 0 & years_after_intro < 4),
              aes(x = (1-frac)*100, y = 6),
              fill = NA, color = "red", width = 10, lwd = 1) +
  geom_tile(data = traj_grid_biased_lateintro2 |> filter((1-frac)*100 == 60 & years_after_intro < 4 & M0 ==3),
            aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "white", width = 10, lwd = 1) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_viridis_c(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
  # scale_fill_viridis_c() +
    scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "Cumulative\nintervention\neffect\n(%)") +
  theme_bw()

traj_grid_biased_lateintro2 |>
  filter((1-frac)*100 == 60 & years_after_intro < 4 & M0 == 3) |>
  dplyr::select(years_after_intro, M0, frac, t_eff_cumul, cuminc_ratio_0)

# -- Tileplots comparing true to biased (M0 = 6) counterfactual
comp_lateintro <- traj_grid_lateintro |>
  mutate(wmel = (1-frac)*100,
         true_incratio = inc_ratio_0,
         true_cumincratio = cuminc_ratio_0) |>
  dplyr::select(M0, wmel, true_incratio, true_cumincratio, years_after_intro) |>
  left_join(traj_grid_biased_lateintro2 |>
              mutate(wmel = (1-frac)*100,
                     biased_incratio = inc_ratio_0,
                     biased_cumincratio = cuminc_ratio_0) |>
              dplyr::select(M0, wmel, biased_incratio, biased_cumincratio, years_after_intro),
            by = c("M0", "wmel", "years_after_intro")) |>
  ungroup() |>
  mutate(diff_incratio = true_incratio - biased_incratio,
         diff_cumincratio = true_cumincratio - biased_cumincratio,
         ratio_incratio = true_incratio / biased_incratio,
         ratio_cumincratio = true_cumincratio / biased_cumincratio)

## yearly
comp_lateintro |>
  ggplot() +
  geom_tile(aes(x = wmel, y = M0, fill = diff_incratio)) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_distiller(palette = "Spectral",
                       limits = c(-max(abs(comp_lateintro$diff_incratio), na.rm = TRUE),
                                  max(abs(comp_lateintro$diff_incratio), na.rm = TRUE))) +
  scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Ratio of mosquitoes to hosts",
       fill = "True - biased yearly incidence") +
  theme_bw()

comp_lateintro |>
  mutate(diff_incratio = cut(diff_incratio, breaks = c(-0.4, -0.2, -0.1, -0.05, 0.05, 0.1, 0.2, 0.4))) |>
  ggplot() +
  geom_tile(aes(x = wmel, y = M0, fill = diff_incratio)) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_brewer(palette = "Spectral", drop = FALSE) +
  scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Ratio of mosquitoes to hosts",
       fill = "True - biased yearly incidence") +
  theme_bw()

comp_lateintro |>
   mutate(diff_incratio = factor(case_when(
    diff_incratio > 0 ~ "True incidence higher",
    diff_incratio == 0 ~ "No bias",
    diff_incratio < 0 ~ "True incidence lower"),
    levels = c("True incidence higher", "No bias", "True incidence lower"))) |>
  ggplot() +
  geom_tile(aes(x = wmel, y = M0, fill = diff_incratio)) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_brewer(palette = "Spectral", drop = FALSE) +
  scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Ratio of mosquitoes to hosts",
       fill = "True - biased yearly incidence") +
  theme_bw()

## cumulative
comp_lateintro |>
  ggplot() +
  geom_tile(aes(x = wmel, y = M0, fill = diff_cumincratio)) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_distiller(palette = "Spectral",
                       limits = c(-max(abs(comp_lateintro$diff_cumincratio), na.rm = TRUE),
                                  max(abs(comp_lateintro$diff_cumincratio), na.rm = TRUE))) +
  scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Ratio of mosquitoes to hosts",
       fill = "True - biased cumulative incidence") +
  theme_bw()

comp_lateintro |>
  mutate(diff_cumincratio = cut(diff_cumincratio, breaks = c(-0.4, -0.2, -0.1, -0.05, 0.05, 0.1, 0.2, 0.4))) |>
  ggplot() +
  geom_tile(aes(x = wmel, y = M0, fill = diff_cumincratio)) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_brewer(palette = "Spectral", drop = FALSE) +
  scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Ratio of mosquitoes to hosts",
       fill = "True - biased cumulative incidence") +
  theme_bw()

comp_lateintro |>
  mutate(diff_cumincratio = factor(case_when(
    diff_cumincratio > 0 ~ "True incidence higher",
    diff_cumincratio == 0 ~ "No bias",
    diff_cumincratio < 0 ~ "True incidence lower"),
    levels = c("True incidence higher", "No bias", "True incidence lower")))|>
  ggplot() +
  geom_tile(aes(x = wmel, y = M0, fill = diff_cumincratio)) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_brewer(palette = "Spectral", drop = FALSE) +
  scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Ratio of mosquitoes to hosts",
       fill = "True - biased cumulative incidence") +
  theme_bw()

# ---- Figure 5 ----------------------------------------------------------
## Choose cell to highlight
highlight_wmel <- 50
highlight_m0 <- 3

# -- Relative incidence
## Get relative incidences for highlighted cells
traj_grid_earlyintro |>
  filter((1-frac)*100 == highlight_wmel & years_after_intro < 4 & M0 == highlight_m0) |>
  dplyr::select(years_after_intro, M0, frac, t_eff_yearly, inc_ratio_0, t_eff_cumul, cuminc_ratio_0)

traj_grid_biased_earlyintro2 |>
  filter((1-frac)*100 == highlight_wmel & years_after_intro < 4 & M0 == highlight_m0) |>
  dplyr::select(years_after_intro, M0, frac, t_eff_yearly, inc_ratio_0, t_eff_cumul, cuminc_ratio_0)

traj_grid_lateintro |>
  filter((1-frac)*100 == highlight_wmel & years_after_intro < 4 & M0 == highlight_m0) |>
  dplyr::select(years_after_intro, M0, frac, t_eff_yearly, inc_ratio_0, t_eff_cumul, cuminc_ratio_0)

traj_grid_biased_lateintro2 |>
  filter((1-frac)*100 == highlight_wmel & years_after_intro < 4 & M0 == highlight_m0) |>
  dplyr::select(years_after_intro, M0, frac, t_eff_yearly, inc_ratio_0, t_eff_cumul, cuminc_ratio_0)

n1 <- traj_grid_earlyintro |>
  filter(years_after_intro < 4) |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = inc_ratio_0)) +
  geom_tile(data = traj_grid_earlyintro |> filter((1-frac)*100 == 0 & years_after_intro < 4),
              aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "red", width = 10, lwd = 1) +
  # geom_tile(data = traj_grid_earlyintro |> filter((1-frac)*100 == highlight_wmel & years_after_intro < 4 & M0 == highlight_m0),
  #           aes(x = (1-frac)*100, y = M0),
  #             fill = NA, color = "white", width = 10, lwd = 1) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_viridis_c(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  # scale_fill_viridis_c() +
    scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "Incidence relative     \nto 0% wMel",
       title = "wMel introduced near full susceptibility") +
  theme_bw() +
  theme(legend.position = "left")

n2 <- traj_grid_earlyintro |>
  filter(years_after_intro < 4) |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = cuminc_ratio_0)) +
  geom_tile(data = traj_grid_earlyintro |> filter((1-frac)*100 == 0 & years_after_intro < 4),
              aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "red", width = 10, lwd = 1) +
  # geom_tile(data = traj_grid_earlyintro |> filter((1-frac)*100 == highlight_wmel & years_after_intro < 4 & M0 == highlight_m0),
  #           aes(x = (1-frac)*100, y = M0),
  #             fill = NA, color = "white", width = 10, lwd = 1) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_viridis_c(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  # scale_fill_viridis_c() +
    scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "Cumulative incidence\nrelative to 0% wMel") +
  theme_bw() +
  theme(legend.position = "left")

n3 <- traj_grid_biased_earlyintro2 |>
  filter(years_after_intro < 4) |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = inc_ratio_0)) +
  geom_tile(data = traj_grid_biased_earlyintro2 |> filter((1-frac)*100 == 0 & years_after_intro < 4),
              aes(x = (1-frac)*100, y = 6),
              fill = NA, color = "red", width = 10, lwd = 1) +
  geom_tile(data = traj_grid_biased_earlyintro2 |> filter((1-frac)*100 == highlight_wmel & years_after_intro < 4 & M0 == highlight_m0),
            aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "white", width = 10, lwd = 1) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_viridis_c(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  # scale_fill_viridis_c() +
    scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "Incidence relative\nto 0% wMel") +
  theme_bw() +
  theme(legend.position = "left")

n4 <- traj_grid_biased_earlyintro2|>
  filter(years_after_intro < 4) |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = cuminc_ratio_0)) +
  geom_tile(data = traj_grid_biased_earlyintro2 |> filter((1-frac)*100 == 0 & years_after_intro < 4),
              aes(x = (1-frac)*100, y = 6),
              fill = NA, color = "red", width = 10, lwd = 1) +
  geom_tile(data = traj_grid_biased_earlyintro2|> filter((1-frac)*100 == highlight_wmel & years_after_intro < 4 & M0 == highlight_m0),
            aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "white", width = 10, lwd = 1) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_viridis_c(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  # scale_fill_viridis_c() +
    scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "Cumulative incidence\nrelative to 0% wMel") +
  theme_bw() +
  theme(legend.position = "left")

n5 <- traj_grid_lateintro |>
  filter(years_after_intro < 4) |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = inc_ratio_0)) +
  geom_tile(data = traj_grid_lateintro |> filter((1-frac)*100 == 0 & years_after_intro < 4),
              aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "red", width = 10, lwd = 1) +
  # geom_tile(data = traj_grid_lateintro |> filter((1-frac)*100 == highlight_wmel & years_after_intro < 4 & M0 == highlight_m0),
  #           aes(x = (1-frac)*100, y = M0),
  #             fill = NA, color = "white", width = 10, lwd = 1) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_viridis_c(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  # scale_fill_viridis_c() +
    scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "Incidence relative     \nto 0% wMel",
       title = "wMel introduced near equilibrium susceptibility") +
  theme_bw() +
  theme(legend.position = "left")

n6 <- traj_grid_lateintro |>
  filter(years_after_intro < 4) |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = cuminc_ratio_0)) +
  geom_tile(data = traj_grid_lateintro |> filter((1-frac)*100 == 0 & years_after_intro < 4),
              aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "red", width = 10, lwd = 1) +
  # geom_tile(data = traj_grid_lateintro |> filter((1-frac)*100 == highlight_wmel & years_after_intro < 4 & M0 == highlight_m0),
  #           aes(x = (1-frac)*100, y = M0),
  #             fill = NA, color = "white", width = 10, lwd = 1) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_viridis_c(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  # scale_fill_viridis_c() +
    scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "Cumulative incidence\nrelative to 0% wMel") +
  theme_bw() +
  theme(legend.position = "left")

n7 <- traj_grid_biased_lateintro2 |>
  filter(years_after_intro < 4) |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = inc_ratio_0)) +
  geom_tile(data = traj_grid_biased_lateintro2 |> filter((1-frac)*100 == 0 & years_after_intro < 4),
              aes(x = (1-frac)*100, y = 6),
              fill = NA, color = "red", width = 10, lwd = 1) +
  geom_tile(data = traj_grid_biased_lateintro2 |> filter((1-frac)*100 == highlight_wmel & years_after_intro < 4 & M0 == highlight_m0),
            aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "white", width = 10, lwd = 1) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_viridis_c(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  # scale_fill_viridis_c() +
    scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "Incidence relative\nto 0% wMel") +
  theme_bw() +
  theme(legend.position = "left")

n8 <- traj_grid_biased_lateintro2|>
  filter(years_after_intro < 4) |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = cuminc_ratio_0)) +
  geom_tile(data = traj_grid_biased_lateintro2 |> filter((1-frac)*100 == 0 & years_after_intro < 4),
              aes(x = (1-frac)*100, y = 6),
              fill = NA, color = "red", width = 10, lwd = 1) +
  geom_tile(data = traj_grid_biased_lateintro2|> filter((1-frac)*100 == highlight_wmel & years_after_intro < 4 & M0 == highlight_m0),
            aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "white", width = 10, lwd = 1) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_viridis_c(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  # scale_fill_viridis_c() +
    scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "Cumulative incidence\nrelative to 0% wMel") +
  theme_bw() +
  theme(legend.position = "left")

# -- With gray boxes
n1 <- traj_grid_earlyintro |>
  filter(years_after_intro < 4) |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = inc_ratio_0)) +
  geom_tile(data = traj_grid_earlyintro |> filter((1-frac)*100 == 0 & years_after_intro < 4),
              aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "red", width = 10, lwd = 0.6) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_viridis_c(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "Incidence relative     \nto 0% wMel") +
  theme_bw() +
  theme(legend.position = "left",
        legend.title = element_text(size = 11),
        legend.key.height = unit(5, "mm"))

n2 <- traj_grid_earlyintro |>
  filter(years_after_intro < 4) |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = cuminc_ratio_0)) +
  geom_tile(data = traj_grid_earlyintro |> filter((1-frac)*100 == 0 & years_after_intro < 4),
              aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "red", width = 10, lwd = 0.6) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_viridis_c(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "Cumulative incidence\nrelative to 0% wMel") +
  theme_bw() +
  theme(legend.position = "left",
        legend.title = element_text(size = 11),
        legend.key.height = unit(5, "mm"))

n5 <- traj_grid_lateintro |>
  filter(years_after_intro < 4) |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = inc_ratio_0)) +
  geom_tile(data = traj_grid_lateintro |> filter((1-frac)*100 == 0 & years_after_intro < 4),
              aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "red", width = 10, lwd = 0.6) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_viridis_c(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "Incidence relative     \nto 0% wMel") +
  theme_bw() +
  theme(legend.position = "left",
        legend.title = element_text(size = 11),
        legend.key.height = unit(5, "mm"))

n6 <- traj_grid_lateintro |>
  filter(years_after_intro < 4) |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = cuminc_ratio_0)) +
  geom_tile(data = traj_grid_lateintro |> filter((1-frac)*100 == 0 & years_after_intro < 4),
              aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "red", width = 10, lwd = 0.6) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_viridis_c(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "Cumulative incidence\nrelative to 0% wMel") +
  theme_bw() +
  theme(legend.position = "left",
        legend.title = element_text(size = 11),
        legend.key.height = unit(5, "mm"))

title1 <- wrap_elements(grid::textGrob(expression(bold("wMel introduced near equilibrium susceptibility")))) +
  theme(plot.margin = margin(-10, 0, -10, 0))
title2 <- wrap_elements(grid::textGrob(expression(bold("wMel introduced near full susceptibility")))) +
  theme(plot.margin = margin(-10, 0, -10, 0))

(title1 + title2) / ((n5 + n1) + plot_layout(guides = "collect")) / ((n6 + n2) + plot_layout(guides = "collect")) +
  plot_annotation(tag_levels = list(c("", "", "A", "B", "C", "D"))) + plot_layout(heights = c(0.25, 2, 2))

# ggsave("outputs/2026-03-12_fig5.pdf", width = 12, height = 4.5)

((n5 + n1) + plot_layout(guides = "collect")) / ((n6 + n2) + plot_layout(guides = "collect")) / ((n7 + n3) + plot_layout(guides = "collect")) / ((n8 + n4) + plot_layout(guides = "collect"))   &
  theme(legend.position = "left")

# ggsave("outputs/2026-03-05_fig5_highlight50.pdf", width = 15, height = 10)

# -- With red boxes
(n1 <- traj_grid_earlyintro |>
  filter(years_after_intro < 4) |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = inc_ratio_0)) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_viridis_c(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "Incidence relative     \nto 0% wMel") +
  ggnewscale::new_scale_fill() +
  geom_tile(data = traj_grid_earlyintro |>
              filter(years_after_intro < 4 & inc_ratio_0 > 1) |>
              mutate(inc_ratio_0_bin = case_when(inc_ratio_0 <= 1 ~ NA,
                                                 inc_ratio_0 <= 1.5 ~ "1-1.5",
                                                 inc_ratio_0 > 1.5 ~ ">1.5"),
                     inc_ratio_0_bin = factor(inc_ratio_0_bin, levels = c(">1.5", "1-1.5"))),
              aes(x = (1-frac)*100, y = M0, fill = inc_ratio_0_bin),
            show.legend = TRUE) +
  geom_tile(data = traj_grid_earlyintro |> filter((1-frac)*100 == 0 & years_after_intro < 4),
              aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "black", width = 10, lwd = 0.6) +
  scale_fill_manual(values = c("1-1.5" = "#E8634A", ">1.5" = "#8B1A1A"), na.value = NA,
                    name = "Incidence relative     \nto 0% wMel", drop = FALSE) +
  theme_bw() +
  theme(legend.position = "left",
        legend.title = element_text(size = 11),
        legend.key.height = unit(5, "mm")))

(n2 <- traj_grid_earlyintro |>
  filter(years_after_intro < 4) |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = cuminc_ratio_0)) +
  geom_tile(data = traj_grid_earlyintro |> filter((1-frac)*100 == 0 & years_after_intro < 4),
              aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "black", width = 10, lwd = 0.6) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_viridis_c(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "Cumulative incidence\nrelative to 0% wMel") +
  theme_bw() +
  theme(legend.position = "left",
        legend.title = element_text(size = 11),
        legend.key.height = unit(5, "mm")))

(n5 <- traj_grid_lateintro |>
  filter(years_after_intro < 4) |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = inc_ratio_0)) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_viridis_c(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "Incidence relative     \nto 0% wMel") +
  ggnewscale::new_scale_fill() +
  geom_tile(data = traj_grid_lateintro |>
              filter(years_after_intro < 4 & inc_ratio_0 > 1) |>
              mutate(inc_ratio_0_bin = case_when(inc_ratio_0 <= 1 ~ NA,
                                                 inc_ratio_0 <= 1.5 ~ "1-1.5",
                                                 inc_ratio_0 > 1.5 ~ ">1.5"),
                     inc_ratio_0_bin = factor(inc_ratio_0_bin, levels = c(">1.5", "1-1.5"))),
              aes(x = (1-frac)*100, y = M0, fill = inc_ratio_0_bin),
            show.legend = TRUE) +
  geom_tile(data = traj_grid_lateintro |> filter((1-frac)*100 == 0 & years_after_intro < 4),
              aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "black", width = 10, lwd = 0.6) +
  scale_fill_manual(values = c("1-1.5" = "#E8634A", ">1.5" = "#8B1A1A"), na.value = NA,
                    name = "Incidence relative     \nto 0% wMel", drop = FALSE) +
  theme_bw() +
  theme(legend.position = "left",
        legend.title = element_text(size = 11),
        legend.key.height = unit(5, "mm")))

(n6 <- traj_grid_lateintro |>
  filter(years_after_intro < 4) |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = cuminc_ratio_0)) +
  geom_tile(data = traj_grid_lateintro |> filter((1-frac)*100 == 0 & years_after_intro < 4),
              aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "black", width = 10, lwd = 0.6) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_viridis_c(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "Cumulative incidence\nrelative to 0% wMel") +
  theme_bw() +
  theme(legend.position = "left",
        legend.title = element_text(size = 11),
        legend.key.height = unit(5, "mm")))

title1 <- wrap_elements(grid::textGrob(expression(bold("wMel introduced near equilibrium susceptibility")))) +
  theme(plot.margin = margin(-10, 0, -10, 0))
title2 <- wrap_elements(grid::textGrob(expression(bold("wMel introduced near full susceptibility")))) +
  theme(plot.margin = margin(-10, 0, -10, 0))

(title1 + title2) / ((n5 + n1) + plot_layout(guides = "collect")) / ((n6 + n2) + plot_layout(guides = "collect")) +
  plot_annotation(tag_levels = list(c("", "", "A", "B", "C", "D"))) + plot_layout(heights = c(0.25, 2, 2))

# ggsave("outputs/2026-07-17_fig5.pdf", width = 12, height = 4.5)

# -- Relative incidence with three scenarios
## Get relative incidences for highlighted cells
traj_grid_earlyintro |>
  filter((1-frac)*100 == highlight_wmel & years_after_intro < 4 & M0 == highlight_m0) |>
  dplyr::select(years_after_intro, M0, frac, t_eff_yearly, inc_ratio_0, t_eff_cumul, cuminc_ratio_0)

traj_grid_biased_earlyintro2 |>
  filter((1-frac)*100 == highlight_wmel & years_after_intro < 4 & M0 == highlight_m0) |>
  dplyr::select(years_after_intro, M0, frac, t_eff_yearly, inc_ratio_0, t_eff_cumul, cuminc_ratio_0)

traj_grid_midintro |>
  filter((1-frac)*100 == highlight_wmel & years_after_intro < 4 & M0 == highlight_m0) |>
  dplyr::select(years_after_intro, M0, frac, t_eff_yearly, inc_ratio_0, t_eff_cumul, cuminc_ratio_0)

traj_grid_biased_midintro2 |>
  filter((1-frac)*100 == highlight_wmel & years_after_intro < 4 & M0 == highlight_m0) |>
  dplyr::select(years_after_intro, M0, frac, t_eff_yearly, inc_ratio_0, t_eff_cumul, cuminc_ratio_0)

traj_grid_lateintro |>
  filter((1-frac)*100 == highlight_wmel & years_after_intro < 4 & M0 == highlight_m0) |>
  dplyr::select(years_after_intro, M0, frac, t_eff_yearly, inc_ratio_0, t_eff_cumul, cuminc_ratio_0)

traj_grid_biased_lateintro2 |>
  filter((1-frac)*100 == highlight_wmel & years_after_intro < 4 & M0 == highlight_m0) |>
  dplyr::select(years_after_intro, M0, frac, t_eff_yearly, inc_ratio_0, t_eff_cumul, cuminc_ratio_0)

n1 <- traj_grid_earlyintro |>
  filter(years_after_intro < 4) |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = inc_ratio_0)) +
  geom_tile(data = traj_grid_earlyintro |> filter((1-frac)*100 == 0 & years_after_intro < 4),
              aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "red", width = 10, lwd = 1) +
  geom_tile(data = traj_grid_earlyintro |> filter((1-frac)*100 == highlight_wmel & years_after_intro < 4 & M0 == highlight_m0),
            aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "white", width = 10, lwd = 1) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_viridis_c(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  # scale_fill_viridis_c() +
    scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "Incidence relative\nto 0% wMel") +
  theme_bw() +
  theme(legend.position = "left")

n2 <- traj_grid_earlyintro |>
  filter(years_after_intro < 4) |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = cuminc_ratio_0)) +
  geom_tile(data = traj_grid_earlyintro |> filter((1-frac)*100 == 0 & years_after_intro < 4),
              aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "red", width = 10, lwd = 1) +
  geom_tile(data = traj_grid_earlyintro |> filter((1-frac)*100 == highlight_wmel & years_after_intro < 4 & M0 == highlight_m0),
            aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "white", width = 10, lwd = 1) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_viridis_c(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  # scale_fill_viridis_c() +
    scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "Cumulative incidence\nrelative to 0% wMel") +
  theme_bw() +
  theme(legend.position = "left")

n3 <- traj_grid_biased_earlyintro2 |>
  filter(years_after_intro < 4) |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = inc_ratio_0)) +
  geom_tile(data = traj_grid_biased_earlyintro2 |> filter((1-frac)*100 == 0 & years_after_intro < 4),
              aes(x = (1-frac)*100, y = 6),
              fill = NA, color = "red", width = 10, lwd = 1) +
  geom_tile(data = traj_grid_biased_earlyintro2 |> filter((1-frac)*100 == highlight_wmel & years_after_intro < 4 & M0 == highlight_m0),
            aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "white", width = 10, lwd = 1) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_viridis_c(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  # scale_fill_viridis_c() +
    scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "Incidence relative\nto 0% wMel") +
  theme_bw() +
  theme(legend.position = "left")

n4 <- traj_grid_biased_earlyintro2|>
  filter(years_after_intro < 4) |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = cuminc_ratio_0)) +
  geom_tile(data = traj_grid_biased_earlyintro2 |> filter((1-frac)*100 == 0 & years_after_intro < 4),
              aes(x = (1-frac)*100, y = 6),
              fill = NA, color = "red", width = 10, lwd = 1) +
  geom_tile(data = traj_grid_biased_earlyintro2|> filter((1-frac)*100 == highlight_wmel & years_after_intro < 4 & M0 == highlight_m0),
            aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "white", width = 10, lwd = 1) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_viridis_c(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  # scale_fill_viridis_c() +
    scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "Cumulative incidence\nrelative to 0% wMel") +
  theme_bw() +
  theme(legend.position = "left")

n5 <- traj_grid_lateintro |>
  filter(years_after_intro < 4) |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = inc_ratio_0)) +
  geom_tile(data = traj_grid_lateintro |> filter((1-frac)*100 == 0 & years_after_intro < 4),
              aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "red", width = 10, lwd = 1) +
  geom_tile(data = traj_grid_lateintro |> filter((1-frac)*100 == highlight_wmel & years_after_intro < 4 & M0 == highlight_m0),
            aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "white", width = 10, lwd = 1) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_viridis_c(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  # scale_fill_viridis_c() +
    scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "Incidence relative\nto 0% wMel") +
  theme_bw() +
  theme(legend.position = "left")

n6 <- traj_grid_lateintro |>
  filter(years_after_intro < 4) |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = cuminc_ratio_0)) +
  geom_tile(data = traj_grid_lateintro |> filter((1-frac)*100 == 0 & years_after_intro < 4),
              aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "red", width = 10, lwd = 1) +
  geom_tile(data = traj_grid_lateintro |> filter((1-frac)*100 == highlight_wmel & years_after_intro < 4 & M0 == highlight_m0),
            aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "white", width = 10, lwd = 1) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_viridis_c(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  # scale_fill_viridis_c() +
    scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "Cumulative incidence\nrelative to 0% wMel") +
  theme_bw() +
  theme(legend.position = "left")

n7 <- traj_grid_biased_lateintro2 |>
  filter(years_after_intro < 4) |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = inc_ratio_0)) +
  geom_tile(data = traj_grid_biased_lateintro2 |> filter((1-frac)*100 == 0 & years_after_intro < 4),
              aes(x = (1-frac)*100, y = 6),
              fill = NA, color = "red", width = 10, lwd = 1) +
  geom_tile(data = traj_grid_biased_lateintro2 |> filter((1-frac)*100 == highlight_wmel & years_after_intro < 4 & M0 == highlight_m0),
            aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "white", width = 10, lwd = 1) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_viridis_c(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  # scale_fill_viridis_c() +
    scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "Incidence relative\nto 0% wMel") +
  theme_bw() +
  theme(legend.position = "left")

n8 <- traj_grid_biased_lateintro2|>
  filter(years_after_intro < 4) |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = cuminc_ratio_0)) +
  geom_tile(data = traj_grid_biased_lateintro2 |> filter((1-frac)*100 == 0 & years_after_intro < 4),
              aes(x = (1-frac)*100, y = 6),
              fill = NA, color = "red", width = 10, lwd = 1) +
  geom_tile(data = traj_grid_biased_lateintro2|> filter((1-frac)*100 == highlight_wmel & years_after_intro < 4 & M0 == highlight_m0),
            aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "white", width = 10, lwd = 1) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_viridis_c(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  # scale_fill_viridis_c() +
    scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "Cumulative incidence\nrelative to 0% wMel") +
  theme_bw() +
  theme(legend.position = "left")

(q1 <- traj_grid_midintro |>
  filter(years_after_intro < 4) |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = inc_ratio_0)) +
  geom_tile(data = traj_grid_midintro |> filter((1-frac)*100 == 0 & years_after_intro < 4),
              aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "red", width = 10, lwd = 1) +
  geom_tile(data = traj_grid_midintro |> filter((1-frac)*100 == 60 & years_after_intro < 4 & M0 ==3),
            aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "white", width = 10, lwd = 1) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_viridis_c(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  # scale_fill_viridis_c() +
    scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "Incidence relative\nto 0% wMel") +
  theme_bw())

(q2 <- traj_grid_midintro |>
  filter(years_after_intro < 4) |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = cuminc_ratio_0)) +
  geom_tile(data = traj_grid_midintro |> filter((1-frac)*100 == 0 & years_after_intro < 4),
              aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "red", width = 10, lwd = 1) +
  geom_tile(data = traj_grid_midintro |> filter((1-frac)*100 == 60 & years_after_intro < 4 & M0 ==3),
            aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "white", width = 10, lwd = 1) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_viridis_c(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  # scale_fill_viridis_c() +
    scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "Cumulative incidence\nrelative to 0% wMel") +
  theme_bw())

(q3 <- traj_grid_biased_midintro2 |>
  filter(years_after_intro < 4) |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = inc_ratio_0)) +
  geom_tile(data = traj_grid_biased_midintro2 |> filter((1-frac)*100 == 0 & years_after_intro < 4),
              aes(x = (1-frac)*100, y = 6),
              fill = NA, color = "red", width = 10, lwd = 1) +
  geom_tile(data = traj_grid_biased_midintro2 |> filter((1-frac)*100 == 60 & years_after_intro < 4 & M0 ==3),
            aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "white", width = 10, lwd = 1) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_viridis_c(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  # scale_fill_viridis_c() +
    scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "Incidence relative\nto 0% wMel") +
  theme_bw())

(q4 <- traj_grid_biased_midintro2|>
  filter(years_after_intro < 4) |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = cuminc_ratio_0)) +
  geom_tile(data = traj_grid_biased_midintro2 |> filter((1-frac)*100 == 0 & years_after_intro < 4),
              aes(x = (1-frac)*100, y = 6),
              fill = NA, color = "red", width = 10, lwd = 1) +
  geom_tile(data = traj_grid_biased_midintro2|> filter((1-frac)*100 == 60 & years_after_intro < 4 & M0 ==3),
            aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "white", width = 10, lwd = 1) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_viridis_c(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  # scale_fill_viridis_c() +
    scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "Cumulative incidence\nrelative to 0% wMel") +
  theme_bw())

((n5 + q1 + n1) + plot_layout(guides = "collect")) / ((n6 + q2 + n2) + plot_layout(guides = "collect")) / ((n7 + q3 + n3) + plot_layout(guides = "collect")) / ((n8 + q4 + n4) + plot_layout(guides = "collect"))   &
  theme(legend.position = "left")

# ggsave("outputs/2026-03-09_fig5_highlight50-3scenarios.pdf", width = 22.5, height = 10)

# -- Supp Figure 50% suscept
q1 <- traj_grid_midintro |>
  filter(years_after_intro < 4) |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = inc_ratio_0)) +
  geom_tile(data = traj_grid_midintro |> filter((1-frac)*100 == 0 & years_after_intro < 4),
              aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "black", width = 10, lwd = 0.6) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_viridis_c(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  # scale_fill_viridis_c() +
    scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "Incidence relative\nto 0% wMel") +
  theme_bw() +
  theme(legend.title = element_text(size = 10.5),
        legend.key.height = unit(5, "mm"))

q2 <- traj_grid_midintro |>
  filter(years_after_intro < 4) |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = cuminc_ratio_0)) +
  geom_tile(data = traj_grid_midintro |> filter((1-frac)*100 == 0 & years_after_intro < 4),
              aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "black", width = 10, lwd = 0.6) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_viridis_c(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  # scale_fill_viridis_c() +
    scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "Cumulative incidence\nrelative to 0% wMel") +
  theme_bw() +
  theme(legend.title = element_text(size = 10.5),
        legend.key.height = unit(5, "mm"))

title3 <- wrap_elements(grid::textGrob(expression(bold("wMel introduced near 50% susceptibility")))) +
  theme(plot.margin = margin(-10, 0, -10, 0))

title3 / q1 / q2 +
  plot_annotation(tag_levels = list(c("", "A", "B"))) +
  plot_layout(heights = c(0.25, 2, 2))

# ggsave("outputs/2026-03-12_figS4.pdf", width = 7, height = 5)

q1 <- traj_grid_midintro |>
  filter(years_after_intro < 4) |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = inc_ratio_0)) +
  scale_fill_viridis_c(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1),
                       name = "Incidence relative\nto 0% wMel") +
  ggnewscale::new_scale_fill() +
  geom_tile(data = traj_grid_midintro |>
              filter(years_after_intro < 4 & inc_ratio_0 > 1) |>
              mutate(inc_ratio_0_bin = case_when(inc_ratio_0 <= 1 ~ NA,
                                                 inc_ratio_0 <= 1.5 ~ "1-1.5",
                                                 inc_ratio_0 > 1.5 ~ ">1.5"),
                     inc_ratio_0_bin = factor(inc_ratio_0_bin, levels = c(">1.5", "1-1.5"))),
              aes(x = (1-frac)*100, y = M0, fill = inc_ratio_0_bin),
            show.legend = TRUE) +
  scale_fill_manual(values = c("1-1.5" = "#E8634A", ">1.5" = "#8B1A1A"), na.value = NA,
                    name = "Incidence relative\nto 0% wMel", drop = FALSE) +
  geom_tile(data = traj_grid_midintro |> filter((1-frac)*100 == 0 & years_after_intro < 4),
              aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "black", width = 10, lwd = 0.6) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "Incidence relative\nto 0% wMel") +
  theme_bw() +
  theme(legend.title = element_text(size = 10.5),
        legend.key.height = unit(5, "mm"))

q2 <- traj_grid_midintro |>
  filter(years_after_intro < 4) |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = cuminc_ratio_0)) +
  scale_fill_viridis_c(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  geom_tile(data = traj_grid_midintro |> filter((1-frac)*100 == 0 & years_after_intro < 4),
              aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "black", width = 10, lwd = 0.6) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "Cumulative incidence\nrelative to 0% wMel") +
  theme_bw() +
  theme(legend.title = element_text(size = 10.5),
        legend.key.height = unit(5, "mm"))

title3 <- wrap_elements(grid::textGrob(expression(bold("wMel introduced near 50% susceptibility")))) +
  theme(plot.margin = margin(-10, 0, -10, 0))

title3 / q1 / q2 +
  plot_annotation(tag_levels = list(c("", "A", "B"))) +
  plot_layout(heights = c(0.25, 2, 2))

# ggsave("outputs/2026-07-17_figS4.pdf", width = 7, height = 5)

# -- Treatment effect
m1 <- traj_grid_earlyintro |>
  filter(years_after_intro < 4) |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = t_eff_yearly)) +
  geom_tile(data = traj_grid_earlyintro |> filter((1-frac)*100 == 0 & years_after_intro < 4),
              aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "red", width = 10, lwd = 1) +
  geom_tile(data = traj_grid_earlyintro |> filter((1-frac)*100 == highlight_wmel & years_after_intro < 4 & M0 == highlight_m0),
            aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "white", width = 10, lwd = 1) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_viridis_c(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
  scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "Yearly\nintervention\neffect\n(%)") +
  theme_bw() +
  theme(legend.position = "left")

m2 <- traj_grid_earlyintro |>
  filter(years_after_intro < 4) |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = t_eff_cumul)) +
  geom_tile(data = traj_grid_earlyintro |> filter((1-frac)*100 == 0 & years_after_intro < 4),
              aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "red", width = 10, lwd = 1) +
  geom_tile(data = traj_grid_earlyintro |> filter((1-frac)*100 == highlight_wmel & years_after_intro < 4 & M0 == highlight_m0),
            aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "white", width = 10, lwd = 1) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_viridis_c(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
  scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "Cumulative\nintervention\neffect\n(%)") +
  theme_bw() +
  theme(legend.position = "left")

m3 <- traj_grid_biased_earlyintro2 |>
  filter(years_after_intro < 4) |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = t_eff_yearly)) +
  geom_tile(data = traj_grid_biased_earlyintro2 |> filter((1-frac)*100 == 0 & years_after_intro < 4),
              aes(x = (1-frac)*100, y = 6),
              fill = NA, color = "red", width = 10, lwd = 1) +
  geom_tile(data = traj_grid_biased_earlyintro2 |> filter((1-frac)*100 == highlight_wmel & years_after_intro < 4 & M0 == highlight_m0),
            aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "white", width = 10, lwd = 1) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_viridis_c(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
  scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "Yearly\nintervention\neffect\n(%)") +
  theme_bw() +
  theme(legend.position = "left")

m4 <- traj_grid_biased_earlyintro2|>
  filter(years_after_intro < 4) |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = t_eff_cumul)) +
  geom_tile(data = traj_grid_biased_earlyintro2 |> filter((1-frac)*100 == 0 & years_after_intro < 4),
              aes(x = (1-frac)*100, y = 6),
              fill = NA, color = "red", width = 10, lwd = 1) +
  geom_tile(data = traj_grid_biased_earlyintro2|> filter((1-frac)*100 == highlight_wmel & years_after_intro < 4 & M0 == highlight_m0),
            aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "white", width = 10, lwd = 1) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_viridis_c(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
  scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "Cumulative\nintervention\neffect\n(%)") +
  theme_bw() +
  theme(legend.position = "left")

m5 <- traj_grid_lateintro |>
  filter(years_after_intro < 4) |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = t_eff_yearly)) +
  geom_tile(data = traj_grid_lateintro |> filter((1-frac)*100 == 0 & years_after_intro < 4),
              aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "red", width = 10, lwd = 1) +
  geom_tile(data = traj_grid_lateintro |> filter((1-frac)*100 == highlight_wmel & years_after_intro < 4 & M0 == highlight_m0),
            aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "white", width = 10, lwd = 1) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_viridis_c(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
  scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "Yearly\nintervention\neffect\n(%)") +
  theme_bw() +
  theme(legend.position = "left")

m6 <- traj_grid_lateintro |>
  filter(years_after_intro < 4) |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = t_eff_cumul)) +
  geom_tile(data = traj_grid_lateintro |> filter((1-frac)*100 == 0 & years_after_intro < 4),
              aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "red", width = 10, lwd = 1) +
  geom_tile(data = traj_grid_lateintro |> filter((1-frac)*100 == highlight_wmel & years_after_intro < 4 & M0 == highlight_m0),
            aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "white", width = 10, lwd = 1) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_viridis_c(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
  scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "Cumulative\nintervention\neffect\n(%)") +
  theme_bw() +
  theme(legend.position = "left")

m7 <- traj_grid_biased_lateintro2 |>
  filter(years_after_intro < 4) |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = t_eff_yearly)) +
  geom_tile(data = traj_grid_biased_lateintro2 |> filter((1-frac)*100 == 0 & years_after_intro < 4),
              aes(x = (1-frac)*100, y = 6),
              fill = NA, color = "red", width = 10, lwd = 1) +
  geom_tile(data = traj_grid_biased_lateintro2 |> filter((1-frac)*100 == highlight_wmel & years_after_intro < 4 & M0 == highlight_m0),
            aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "white", width = 10, lwd = 1) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_viridis_c(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
  scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "Yearly\nintervention\neffect\n(%)") +
  theme_bw() +
  theme(legend.position = "left")

m8 <- traj_grid_biased_lateintro2|>
  filter(years_after_intro < 4) |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = t_eff_cumul)) +
  geom_tile(data = traj_grid_biased_lateintro2 |> filter((1-frac)*100 == 0 & years_after_intro < 4),
              aes(x = (1-frac)*100, y = 6),
              fill = NA, color = "red", width = 10, lwd = 1) +
  geom_tile(data = traj_grid_biased_lateintro2|> filter((1-frac)*100 == highlight_wmel & years_after_intro < 4 & M0 == highlight_m0),
            aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "white", width = 10, lwd = 1) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_viridis_c(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
  scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "Cumulative\nintervention\neffect\n(%)") +
  theme_bw() +
  theme(legend.position = "left")

((m5 + m1) + plot_layout(guides = "collect")) /
  ((m6 + m2) + plot_layout(guides = "collect")) / ((m7 + m3) + plot_layout(guides = "collect")) / ((m8 + m4) + plot_layout(guides = "collect"))  &
  theme(legend.position = "left")

# ggsave("outputs/2026-03-05_fig5_highlight50_teff.pdf", width = 15, height = 10)

# -- Comparing plots
lim <- max(
  abs(comp_earlyintro$ratio_incratio[comp_earlyintro$years_after_intro<4]),
  # abs(comp_midintro$ratio_incratio[comp_earlyintro$years_after_intro<4]),
  abs(comp_lateintro$ratio_incratio[comp_lateintro$years_after_intro<4]),
  abs(comp_earlyintro$ratio_cumincratio[comp_earlyintro$years_after_intro<4]),
  # abs(comp_midintro$ratio_incratio[comp_lateintro$years_after_intro<4]),
  abs(comp_lateintro$ratio_cumincratio[comp_lateintro$years_after_intro<4]), na.rm = TRUE)

(c1 <- comp_earlyintro |>
  filter(years_after_intro < 4) |>
  ggplot() +
  geom_tile(aes(x = wmel, y = M0, fill = ratio_incratio)) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_distiller(palette = "Spectral",
                       limits = c(-lim, lim)) +
  scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "True/biased\nyearly incidence   ",
       title = "wMel introduced near full susceptibility") +
  theme_bw())

(c2 <- comp_earlyintro |>
  filter(years_after_intro < 4) |>
  ggplot() +
  geom_tile(aes(x = wmel, y = M0, fill = ratio_cumincratio)) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_distiller(palette = "Spectral",
                       limits = c(-lim, lim)) +
  scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "True/biased\ncumulative incidence") +
  theme_bw())

# comp_midintro |>
#   filter(years_after_intro < 4) |>
#   ggplot() +
#   geom_tile(aes(x = wmel, y = M0, fill = ratio_incratio)) +
#   facet_wrap(~years_after_intro,
#              labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
#   scale_fill_distiller(palette = "Spectral",
#                        limits = c(-lim, lim)) +
#   scale_y_continuous(breaks = c(1:10)) +
#   labs(x = "% wMel", y = "Mosquitoes to hosts",
#        fill = "True/biased yearly incidence",
#        title = "wMel introduced with half susceptibility") +
#   theme_bw()
#
# comp_midintro |>
#   filter(years_after_intro < 4) |>
#   ggplot() +
#   geom_tile(aes(x = wmel, y = M0, fill = ratio_cumincratio)) +
#   facet_wrap(~years_after_intro,
#              labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
#   scale_fill_distiller(palette = "Spectral",
#                        limits = c(-lim, lim)) +
#   scale_y_continuous(breaks = c(1:10)) +
#   labs(x = "% wMel", y = "Mosquitoes to hosts",
#        fill = "True/biased cumulative incidence",
#        title = "wMel introduced near full susceptibility") +
#   theme_bw()

(c3 <- comp_lateintro |>
  filter(years_after_intro < 4) |>
  ggplot() +
  geom_tile(aes(x = wmel, y = M0, fill = ratio_incratio)) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_distiller(palette = "Spectral",
                       limits = c(-lim, lim)) +
  scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "True/biased\nyearly incidence   ",
       title = "wMel introduced near equilibrium susceptibility") +
  theme_bw())

(c4 <- comp_lateintro |>
  filter(years_after_intro < 4) |>
  ggplot() +
  geom_tile(aes(x = wmel, y = M0, fill = ratio_cumincratio)) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_distiller(palette = "Spectral",
                       limits = c(-lim, lim)) +
  scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "True/biased\ncumulative incidence") +
  theme_bw())

((c1 + c3) + plot_layout(guides = "collect")) /
  ((c2 + c4) + plot_layout(guides = "collect"))

min_inc <- min(
  comp_earlyintro$ratio_incratio[comp_earlyintro$years_after_intro<4],
  comp_lateintro$ratio_incratio[comp_lateintro$years_after_intro<4],
  comp_earlyintro$ratio_cumincratio[comp_earlyintro$years_after_intro<4],
  comp_lateintro$ratio_cumincratio[comp_lateintro$years_after_intro<4],
  na.rm = TRUE)

max_inc <- max(
  comp_earlyintro$ratio_incratio[comp_earlyintro$years_after_intro<4],
  comp_lateintro$ratio_incratio[comp_lateintro$years_after_intro<4],
  comp_earlyintro$ratio_cumincratio[comp_earlyintro$years_after_intro<4],
  comp_lateintro$ratio_cumincratio[comp_lateintro$years_after_intro<4],
  na.rm = TRUE)

(c1 <- comp_earlyintro |>
  filter(years_after_intro < 4) |>
  ggplot() +
  geom_tile(aes(x = wmel, y = M0, fill = ratio_incratio)) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  # scale_fill_distiller(palette = "Spectral",
  #                      limits = c(-lim, lim)) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 1,
                       limits = c(min_inc, max_inc), breaks = c(0.2, 0.6, 1, 1.4)) +
  scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "True/biased\nyearly incidence") +
  theme_bw() +
  theme(legend.title = element_text(size = 11),
        legend.key.height = unit(5, "mm")))

(c2 <- comp_earlyintro |>
  filter(years_after_intro < 4) |>
  ggplot() +
  geom_tile(aes(x = wmel, y = M0, fill = ratio_cumincratio)) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  # scale_fill_distiller(palette = "Spectral",
  #                      limits = c(-lim, lim)) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 1,
                       limits = c(min_inc, max_inc), breaks = c(0.2, 0.6, 1, 1.4)) +
  scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "True/biased\ncumulative incidence") +
  theme_bw() +
  theme(legend.title = element_text(size = 11),
        legend.key.height = unit(5, "mm")))

(c3 <- comp_lateintro |>
  filter(years_after_intro < 4) |>
  ggplot() +
  geom_tile(aes(x = wmel, y = M0, fill = ratio_incratio)) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  # scale_fill_distiller(palette = "Spectral",
  #                      limits = c(-lim, lim)) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 1,
                       limits = c(min_inc, max_inc), breaks = c(0.2, 0.6, 1, 1.4)) +
  scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "True/biased\nyearly incidence") +
  theme_bw() +
  theme(legend.title = element_text(size = 11),
        legend.key.height = unit(5, "mm")))

(c4 <- comp_lateintro |>
  filter(years_after_intro < 4) |>
  ggplot() +
  geom_tile(aes(x = wmel, y = M0, fill = ratio_cumincratio)) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  # scale_fill_distiller(palette = "Spectral",
  #                      limits = c(-lim, lim)) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 1,
                       limits = c(min_inc, max_inc), breaks = c(0.2, 0.6, 1, 1.4)) +
  scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "True/biased\ncumulative incidence") +
  theme_bw() +
  theme(legend.title = element_text(size = 11),
        legend.key.height = unit(5, "mm")))

(title1 + title2) / ((c1 + c3) + plot_layout(guides = "collect")) /
  ((c2 + c4) + plot_layout(guides = "collect")) +
  plot_annotation(tag_levels = list(c("", "", "A", "B", "C", "D"))) + plot_layout(heights = c(0.25, 2, 2))

# ggsave("outputs/2026-03-12_figS5.pdf", width = 12, height = 5)

# -- Late measurement of cumulative incidence
traj_grid_earlyintro2 <- trajectories_earlyintro |>
  ungroup() |>
  mutate(year = floor((time - 1) / 365),
         years_after_intro = year+1) |>
  group_by(M0, Mnew, frac, year, years_after_intro) |>
  summarise(Ih = sum(Ih),
            Rhc = sum(Rhc),
            .groups = "drop") |>
  arrange(M0, Mnew, frac, year, years_after_intro) |>
  group_by(M0, Mnew, frac) |>
  mutate(Rhc_post = cumsum(if_else(years_after_intro >= 1, Ih, 0)),
         Rhc_post_late1 = cumsum(if_else(years_after_intro >=2, Ih, 0)),
         Rhc_post_late2 = cumsum(if_else(years_after_intro >=3, Ih, 0)),
         Rhc_post_late3 = cumsum(if_else(years_after_intro >=4, Ih, 0))) |>
  ungroup() |>
  filter(year %in% c(0:10)) |>
  group_by(M0, year) |>
  arrange(desc(Mnew)) |>
  mutate(Ih_intro0 = Ih[1],
         Rhc_intro0 = Rhc_post[1],
         Ih_intro10 = Ih[2],
         Rhc_intro10 = Rhc_post[2],
         Rhc_intro_late1 = Rhc_post_late1[1],
         Rhc_intro_late2 = Rhc_post_late2[1],
         Rhc_intro_late3 = Rhc_post_late3[1]) |>
  ungroup() |>
  rowwise() |>
  mutate(inc_ratio_0 = Ih/Ih_intro0,
         inc_ratio_10 = Ih/Ih_intro10,
         cuminc_ratio_0 = Rhc_post/Rhc_intro0,
         cuminc_ratio_10 = Rhc_post/Rhc_intro10,
         cuminc_ratio_late1 = Rhc_post_late1/Rhc_intro_late1,
         cuminc_ratio_late2 = Rhc_post_late2/Rhc_intro_late2,
         cuminc_ratio_late3 = Rhc_post_late3/Rhc_intro_late3) |>
  ungroup() |>
  filter(years_after_intro < 10)

traj_grid_lateintro2 <-  trajectories_lateintro |>
  ungroup() |>
  mutate(year = floor((time - 1) / 365),
         years_after_intro = (year+1) - year_introduction) |>
  group_by(M0, Mnew, frac, year, years_after_intro) |>
  summarise(Ih = sum(Ih),
            Rhc = sum(Rhc),
            .groups = "drop") |>
  arrange(M0, Mnew, frac, year, years_after_intro) |>
  group_by(M0, Mnew, frac) |>
  mutate(Rhc_post = cumsum(if_else(years_after_intro >= 1, Ih, 0)),
         Rhc_post_late1 = cumsum(if_else(years_after_intro >=2, Ih, 0)),
         Rhc_post_late2 = cumsum(if_else(years_after_intro >=3, Ih, 0)),
         Rhc_post_late3 = cumsum(if_else(years_after_intro >=4, Ih, 0))) |>
  ungroup() |>
  filter(year %in% years_incratio) |>
  group_by(M0, year) |>
  arrange(desc(Mnew)) |>
  mutate(Ih_intro0 = Ih[1],
         Rhc_intro0 = Rhc_post[1],
         Ih_intro10 = Ih[2],
         Rhc_intro10 = Rhc_post[2],
         Rhc_intro_late1 = Rhc_post_late1[1],
         Rhc_intro_late2 = Rhc_post_late2[1],
         Rhc_intro_late3 = Rhc_post_late3[1]) |>
  ungroup() |>
  rowwise() |>
  mutate(inc_ratio_0 = Ih/Ih_intro0,
         inc_ratio_10 = Ih/Ih_intro10,
         cuminc_ratio_0 = Rhc_post/Rhc_intro0,
         cuminc_ratio_10 = Rhc_post/Rhc_intro10,
         cuminc_ratio_late1 = Rhc_post_late1/Rhc_intro_late1,
         cuminc_ratio_late2 = Rhc_post_late2/Rhc_intro_late2,
         cuminc_ratio_late3 = Rhc_post_late3/Rhc_intro_late3) |>
  ungroup() |>
  filter(years_after_intro < 10)

# -- Just years w data
lateplot_a <- traj_grid_lateintro2 |>
  filter(years_after_intro %in% c(2, 3, 4)) |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = cuminc_ratio_late1)) +
  geom_tile(data = traj_grid_lateintro2 |> filter((1-frac)*100 == 0 & years_after_intro %in% c(2, 3, 4)),
              aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "red", width = 10, lwd = 0.6) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_viridis_c(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "Cumulative incidence\nrelative to 0% wMel") +
  theme_bw() +
  theme(legend.position = "left",
        legend.title = element_text(size = 11),
        legend.key.height = unit(5, "mm"),
        plot.background = element_rect(fill = NA))

lateplot_b <- traj_grid_lateintro2 |>
  filter(years_after_intro %in% c(3, 4, 5)) |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = cuminc_ratio_late2)) +
  geom_tile(data = traj_grid_lateintro2 |> filter((1-frac)*100 == 0 & years_after_intro %in% c(3, 4, 5)),
              aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "red", width = 10, lwd = 0.6) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_viridis_c(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "Cumulative incidence\nrelative to 0% wMel") +
  theme_bw() +
  theme(legend.position = "left",
        legend.title = element_text(size = 11),
        legend.key.height = unit(5, "mm"),
        plot.background = element_rect(fill = NA))

lateplot_c <- traj_grid_lateintro2 |>
  filter(years_after_intro %in% c(4, 5, 6)) |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = cuminc_ratio_late3)) +
  geom_tile(data = traj_grid_lateintro2 |> filter((1-frac)*100 == 0 & years_after_intro %in% c(4, 5, 6)),
              aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "red", width = 10, lwd = 0.6) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_viridis_c(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "Cumulative incidence\nrelative to 0% wMel") +
  theme_bw() +
  theme(legend.position = "left",
        legend.title = element_text(size = 11),
        legend.key.height = unit(5, "mm"),
        plot.background = element_rect(fill = NA))

lateplot_d <- traj_grid_earlyintro2 |>
  filter(years_after_intro %in% c(2, 3, 4)) |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = cuminc_ratio_late1)) +
  geom_tile(data = traj_grid_earlyintro2 |> filter((1-frac)*100 == 0 & years_after_intro %in% c(2, 3, 4)),
              aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "red", width = 10, lwd = 0.6) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_viridis_c(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "Cumulative incidence\nrelative to 0% wMel") +
  theme_bw() +
  theme(legend.position = "left",
        legend.title = element_text(size = 11),
        legend.key.height = unit(5, "mm"),
        plot.background = element_rect(fill = NA))

lateplot_e <- traj_grid_earlyintro2 |>
  filter(years_after_intro %in% c(3, 4, 5)) |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = cuminc_ratio_late2)) +
  geom_tile(data = traj_grid_earlyintro2 |> filter((1-frac)*100 == 0 & years_after_intro %in% c(3, 4, 5)),
              aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "red", width = 10, lwd = 0.6) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_viridis_c(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "Cumulative incidence\nrelative to 0% wMel") +
  theme_bw() +
  theme(legend.position = "left",
        legend.title = element_text(size = 11),
        legend.key.height = unit(5, "mm"),
        plot.background = element_rect(fill = NA))

lateplot_f <- traj_grid_earlyintro2 |>
  filter(years_after_intro %in% c(4, 5, 6)) |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = cuminc_ratio_late3)) +
  geom_tile(data = traj_grid_earlyintro2 |> filter((1-frac)*100 == 0 & years_after_intro %in% c(4, 5, 6)),
              aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "red", width = 10, lwd = 0.6) +
  facet_wrap(~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x))) +
  scale_fill_viridis_c(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts",
       fill = "Cumulative incidence\nrelative to 0% wMel") +
  theme_bw() +
  theme(legend.position = "left",
        legend.title = element_text(size = 11),
        legend.key.height = unit(5, "mm"),
        plot.background = element_rect(fill = NA))

title1 <- wrap_elements(grid::textGrob(expression(bold("wMel introduced near equilibrium susceptibility"))))
title2 <- wrap_elements(grid::textGrob(expression(bold("wMel introduced near full susceptibility"))))

(title1 + title2) / (lateplot_a + lateplot_d) / (lateplot_b + lateplot_e) / (lateplot_c + lateplot_f) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = list(c("", "", "A", "B", "C", "D", "E", "F"))) + plot_layout(heights = c(0.25, 2, 2, 2))

# ggsave("outputs/2026-07-16_fig-late-measurement.pdf", width = 12, height = 7)

# -- Grid of years
traj_grid_earlyintro3 <- traj_grid_earlyintro2 |>
  pivot_longer(cols = c(cuminc_ratio_0, cuminc_ratio_late1, cuminc_ratio_late2, cuminc_ratio_late3),
               names_to = "lag_pre", values_to = "cuminc_ratio") |>
  mutate(lag = case_when(lag_pre == "cuminc_ratio_0" ~ "No delay",
                         lag_pre == "cuminc_ratio_late1" ~ "1-year delay",
                         lag_pre == "cuminc_ratio_late2" ~ "2-year delay",
                         lag_pre == "cuminc_ratio_late3" ~ "3-year delay"),
         lag = factor(lag, levels = c("No delay", "1-year delay", "2-year delay", "3-year delay")))

traj_grid_lateintro3 <- traj_grid_lateintro2 |>
  pivot_longer(cols = c(cuminc_ratio_0, cuminc_ratio_late1, cuminc_ratio_late2, cuminc_ratio_late3),
               names_to = "lag_pre", values_to = "cuminc_ratio") |>
  mutate(lag = case_when(lag_pre == "cuminc_ratio_0" ~ "No delay",
                         lag_pre == "cuminc_ratio_late1" ~ "1-year delay",
                         lag_pre == "cuminc_ratio_late2" ~ "2-year delay",
                         lag_pre == "cuminc_ratio_late3" ~ "3-year delay"),
         lag = factor(lag, levels = c("No delay", "1-year delay", "2-year delay", "3-year delay")))

s5_a <- traj_grid_lateintro3 |>
  filter(years_after_intro %in% c(0:6)) |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = cuminc_ratio)) +
  facet_grid(lag~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x)),
             switch = "y") +
  scale_fill_viridis_c(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1),
                       name = "Cumulative incidence\nrelative to 0% wMel") +
  ggnewscale::new_scale_fill() +
  geom_tile(data = traj_grid_lateintro3 |>
              filter(years_after_intro %in% c(0:6) & cuminc_ratio > 1) |>
              mutate(cuminc_ratio_bin = case_when(cuminc_ratio <= 1 ~ NA,
                                                  cuminc_ratio <= 1.5 ~ "1-1.5",
                                                  cuminc_ratio > 1.5 ~ ">1.5",
                                                  .default = NA),
                     cuminc_ratio_bin = factor(cuminc_ratio_bin, levels = c(">1.5", "1-1.5"))),
              aes(x = (1-frac)*100, y = M0, fill = cuminc_ratio_bin),
            show.legend = TRUE) +
  scale_fill_manual(values = c("1-1.5" = "#E8634A", ">1.5" = "#8B1A1A"), na.value = NA,
                    name = "Cumulative incidence\nrelative to 0% wMel", drop = FALSE) +
  geom_tile(data = traj_grid_lateintro3 |>
              filter((1-frac)*100 == 0 & years_after_intro %in% c(0:6) & !is.na(cuminc_ratio)),
              aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "black", width = 10, lwd = 0.6) +
  scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts", title = "wMel introduced near equilibrium susceptibility") +
  theme_bw() +
  theme(legend.title = element_text(size = 11),
        legend.key.height = unit(5, "mm"),
        plot.background = element_rect(fill = NA))

s5_b <- traj_grid_earlyintro3 |>
  filter(years_after_intro %in% c(0:6)) |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = cuminc_ratio)) +
  facet_grid(lag~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x)),
             switch = "y") +
  scale_fill_viridis_c(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1),
                       name = "Cumulative incidence\nrelative to 0% wMel") +
  ggnewscale::new_scale_fill() +
  geom_tile(data = traj_grid_earlyintro3 |>
              filter(years_after_intro %in% c(0:6) & cuminc_ratio > 1) |>
              mutate(cuminc_ratio_bin = case_when(cuminc_ratio <= 1 ~ NA,
                                                  cuminc_ratio <= 1.5 ~ "1-1.5",
                                                  cuminc_ratio > 1.5 ~ ">1.5",
                                                  .default = NA),
                     cuminc_ratio_bin = factor(cuminc_ratio_bin, levels = c(">1.5", "1-1.5"))),
              aes(x = (1-frac)*100, y = M0, fill = cuminc_ratio_bin),
            show.legend = TRUE) +
  scale_fill_manual(values = c("1-1.5" = "#E8634A", ">1.5" = "#8B1A1A"), na.value = NA,
                    name = "Cumulative incidence\nrelative to 0% wMel", drop = FALSE) +
  geom_tile(data = traj_grid_earlyintro3 |>
              filter((1-frac)*100 == 0 & years_after_intro %in% c(0:6) & !is.na(cuminc_ratio)),
              aes(x = (1-frac)*100, y = M0),
              fill = NA, color = "black", width = 10, lwd = 0.6) +
  scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts", title = "wMel introduced near full susceptibility") +
  theme_bw() +
  theme(legend.title = element_text(size = 11),
        legend.key.height = unit(5, "mm"),
        plot.background = element_rect(fill = NA))

s5_a / s5_b +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A")

# ggsave("outputs/2026-07-16_fig-late-measurement.pdf", width = 12, height = 12)

# -- Differences in measured cuminc
traj_grid_lateintro4 <- traj_grid_lateintro2 |>
  mutate(diff_delay1 = cuminc_ratio_late1 - cuminc_ratio_0,
         diff_delay2 = cuminc_ratio_late2 - cuminc_ratio_0,
         diff_delay3 = cuminc_ratio_late3 - cuminc_ratio_0) |>
  pivot_longer(cols = c(diff_delay1, diff_delay2, diff_delay3),
               names_to = "lag_pre", values_to = "cuminc_ratio_diff") |>
  mutate(lag = case_when(lag_pre == "diff_delay1" ~ "1-year delay",
                         lag_pre == "diff_delay2" ~ "2-year delay",
                         lag_pre == "diff_delay3" ~ "3-year delay"),
         lag = factor(lag, levels = c("1-year delay", "2-year delay", "3-year delay")))

traj_grid_earlyintro4 <- traj_grid_earlyintro2 |>
  mutate(diff_delay1 = cuminc_ratio_late1 - cuminc_ratio_0,
         diff_delay2 = cuminc_ratio_late2 - cuminc_ratio_0,
         diff_delay3 = cuminc_ratio_late3 - cuminc_ratio_0) |>
  pivot_longer(cols = c(diff_delay1, diff_delay2, diff_delay3),
               names_to = "lag_pre", values_to = "cuminc_ratio_diff") |>
  mutate(lag = case_when(lag_pre == "diff_delay1" ~ "1-year delay",
                         lag_pre == "diff_delay2" ~ "2-year delay",
                         lag_pre == "diff_delay3" ~ "3-year delay"),
         lag = factor(lag, levels = c("1-year delay", "2-year delay", "3-year delay")))

min_diff <- min(traj_grid_earlyintro4$cuminc_ratio_diff, traj_grid_lateintro4$cuminc_ratio_diff, na.rm = TRUE)
max_diff <- max(traj_grid_earlyintro4$cuminc_ratio_diff, traj_grid_lateintro4$cuminc_ratio_diff, na.rm = TRUE)

(s6_a <- traj_grid_lateintro4 |>
  filter(years_after_intro %in% c(0:6)) |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = cuminc_ratio_diff)) +
  facet_grid(lag~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x)),
             switch = "y") +
  # scale_fill_gradient2(low = "#2166ac", high = "#d6604d", mid = "#d9d9d9", midpoint = 0,
  #                      limits = c(min_diff, max_diff), transform = "log1p",
  #                      breaks = c(-0.8, 0, 5, 50, 500),
  #                      name =
  #                      "Difference between relative\ncumulative incidence measured\nwith delayed start\nand after intervention") +
    scale_fill_continuous_diverging(
      palette   = "Blue-Red 3",
      mid       = 0,
      limits    = c(min_diff, max_diff),
      transform = "log1p",
      breaks    = c(-0.8, 0, 5, 50, 500),
      name = "Difference in relative\ncumulative incidence\n(delayed − no delay)") +
  scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts", title = "wMel introduced near equilibrium susceptibility") +
  theme_bw() +
  theme(legend.title = element_text(size = 11),
        legend.key.height = unit(5, "mm"),
        plot.background = element_rect(fill = NA)))

(s6_b <- traj_grid_earlyintro4 |>
  filter(years_after_intro %in% c(0:6)) |>
  ggplot() +
  geom_tile(aes(x = (1-frac)*100, y = M0, fill = cuminc_ratio_diff)) +
  facet_grid(lag~years_after_intro,
             labeller = labeller(years_after_intro = function(x) paste0("Year ", x)),
             switch = "y") +
  # scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0,
  #                      limits = c(min_diff, max_diff), transform = "log1p",
  #                      breaks = c(-0.8, 0, 5, 50, 500),
  #                      name =
  #                      "Difference between relative\ncumulative incidence measured\nwith delayed start\nand after intervention") +
    scale_fill_continuous_diverging(
      palette   = "Blue-Red 3",  # or "Cyan-Magenta", "Green-Brown"
      mid       = 0,
      limits    = c(min_diff, max_diff),
      transform = "log1p",
      breaks    = c(-0.8, 0, 5, 50, 500),
      name = "Difference in relative\ncumulative incidence\n(delayed − no delay)") +
  scale_y_continuous(breaks = c(1:10)) +
  labs(x = "% wMel", y = "Mosquitoes to hosts", title = "wMel introduced near full susceptibility") +
  theme_bw() +
  theme(legend.title = element_text(size = 11),
        legend.key.height = unit(5, "mm"),
        plot.background = element_rect(fill = NA)))

s6_a / s6_b +
  plot_annotation(tag_levels = "A")

# ggsave("outputs/2026-07-20_fig-late-measurement-diffs.pdf", width = 12, height = 10)

# ---- Simulations with early introduction of wMel (near-fully susceptible population), change in amplitude ----

# -- Define model parameters
pars = c(N_h = 1e6,              # Host population size.
         S_h = 0.999,            # Starting proportion of susceptible hosts.
         I_h = 0.001,            # Starting proportion of infected hosts.
         M = 4,                  # Ratio of mosquitoes to hosts.
         mu_h = 1 / (77 * 365),  # Host mortality rate.
         r_h = 1 / (77 * 365),   # Host birth rate.
         gamma_h = 1 / 5,        # Host recovery from infectiousness.
         delta_h = 1 / 5.9,      # Intrinsic incubation period.
         # Constant rate of introduction of mosquitoes, as a proportion of `N_h`.
         trickle = 1e-5)

# -- State variables
xstart = vector('numeric', length = 8)
# Sh, Eh, Ih, Rh: susceptible, exposed, infected, recovered hosts.
# Rhc: cumulative recovered hosts (proxy for incidence).
# Sv, Ev, Iv: susceptible, exposed, infected vectors.

names(xstart) = c('Sh', 'Eh', 'Ih', 'Rh', 'Rhc', 'Sv', 'Ev', 'Iv')

xstart['Sv'] = pars['N_h'] * pars['M']
xstart['Sh'] = pars['S_h'] * pars['N_h']
xstart['Ih'] = pars['I_h'] * pars['N_h']
xstart['Rh'] = pars['N_h'] - xstart['Sh'] - xstart['Ih']

make_xstart <- function(pars_vec) {
  x0 <- setNames(numeric(8), c("Sh","Eh","Ih","Rh","Rhc","Sv","Ev","Iv"))
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

# -- Length of simulation
# Total number of days to simulate.
n = 365 * 15
# Need +1 here or else the last time-step is (always?) returned as NA.
ts = seq_len(n + 1)

# -- Define Ms and wMel % to test
M_grid <- seq(1, 10, by = 1)
reduct_fracs <- seq(0.1, 1, by = 0.1)

results_grid <- expand.grid(
  M0    = M_grid,
  frac  = reduct_fracs,
  KEEP.OUT.ATTRS = FALSE
)

# -- Define temperature function
GetTemp = function(x, mean_temp = 26, ampl_temp = 3)
{
  # Determine which year we're in (1-indexed)
  year = floor(x / 365) + 1

  # Smooth transition for amplitude in year 5
  if (year == 5) {
    # Create smooth ramp-up at start of year 5 and ramp-down at end
    day_in_year = x %% 365

    # Smooth transition using sigmoid/cosine taper
    # Ramp up over first 30 days, ramp down over last 30 days
    if (day_in_year < 30) {
      # Smooth ramp up from 1x to 2x amplitude
      scale_factor = 1 + (1 - cos(pi * day_in_year / 30)) / 2
    } else if (day_in_year > 335) {
      # Smooth ramp down from 2x to 1x amplitude
      scale_factor = 1 + (1 + cos(pi * (day_in_year - 335) / 30)) / 2
    } else {
      # Full 2x amplitude in middle of year
      scale_factor = 2
    }

    current_ampl = ampl_temp * scale_factor
  } else {
    current_ampl = ampl_temp
  }

  tmp = current_ampl * sin(2 * pi * x / 365) + mean_temp
  names(tmp) = NULL
  return(tmp)
}

days = 1:(365 * 6)
temps = sapply(days, GetTemp)

df_temps <- data.frame(t = days/365, temp = temps)

df_temps |>
  ggplot(aes(x = t, y = temps)) +
  geom_line() +
  theme_bw()

# -- Run model
trajectories_ampli <- read_rds("data/interim/seirsei_fig5-amplitude change.rds")

# -- Plot outputs
trajectories_ampli |>
  filter(M0 %in% c(1, 6) & (100-frac*100) %in% c(0, 50, 90)) |>
  filter(time/365 < 10 & time/365 >= 1) |>
  ggplot(aes(x = time/365, y = Ih, color = factor(100-frac*100))) +
  geom_line(lwd = 0.8) +
  scale_x_continuous(breaks = c(0:10)) +
  scale_color_brewer(palette = "Set1") +
  labs(x = "Year", y = "Infections",
       color = "%wMel at intervention") +
  facet_wrap(~M0, scales = "free_y", labeller = labeller(M0 = function(x) paste0("M0 = ", x)), ncol = 5) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.byrow = TRUE) +
  coord_cartesian(xlim = c(1,7))

trajectories_ampli |>
  filter((100-frac*100) %in% c(0, 50, 90)) |>
  filter(time/365 < 10) |>
  ggplot(aes(x = time/365, y = Ih, color = factor(100-frac*100))) +
  geom_line(lwd = 0.8) +
  scale_x_continuous(breaks = c(0:10)) +
  scale_color_brewer(palette = "Set1") +
  labs(x = "Year", y = "Infections",
       color = "%wMel at intervention") +
  facet_wrap(~M0, scales = "free_y", labeller = labeller(M0 = function(x) paste0("M0 = ", x)), ncol = 5) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.byrow = TRUE) +
  coord_cartesian(xlim = c(0,7))

trajectories_ampli |>
  filter((100-frac*100) %in% c(0, 50, 90)) |>
  filter(time/365 < 10 & time/365 >= 1) |>
  ggplot(aes(x = time/365, y = Ih, color = factor(100-frac*100))) +
  geom_line(lwd = 0.8) +
  scale_x_continuous(breaks = c(0:10)) +
  scale_color_brewer(palette = "Set1") +
  labs(x = "Year", y = "Infections",
       color = "%wMel at intervention") +
  facet_wrap(~M0, scales = "free_y", labeller = labeller(M0 = function(x) paste0("M0 = ", x)), ncol = 5) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.byrow = TRUE) +
  coord_cartesian(xlim = c(1,7))

# ---- Simulations with late introduction of wMel (near equilibrium susceptibility) ----

# -- Length of simulation and time of introduction of wMel
n_days_total <- 365 * 41
t_break      <- 365 * 20
times1 <- 1:(t_break)
times2 <- (t_break + 1):(n_days_total + 1)

# -- Define temperature function
GetTemp = function(x, mean_temp = 26, ampl_temp = 3)
{
  # Determine which year we're in (1-indexed)
  year = floor(x / 365) + 1

  # Smooth transition for amplitude in year 5
  if (year == 25) {
    # Create smooth ramp-up at start of year 5 and ramp-down at end
    day_in_year = x %% 365

    # Smooth transition using sigmoid/cosine taper
    # Ramp up over first 30 days, ramp down over last 30 days
    if (day_in_year < 30) {
      # Smooth ramp up from 1x to 2x amplitude
      scale_factor = 1 + (1 - cos(pi * day_in_year / 30)) / 2
    } else if (day_in_year > 335) {
      # Smooth ramp down from 2x to 1x amplitude
      scale_factor = 1 + (1 + cos(pi * (day_in_year - 335) / 30)) / 2
    } else {
      # Full 2x amplitude in middle of year
      scale_factor = 2
    }

    current_ampl = ampl_temp * scale_factor
  } else {
    current_ampl = ampl_temp
  }

  tmp = current_ampl * sin(2 * pi * x / 365) + mean_temp
  names(tmp) = NULL
  return(tmp)
}

days = 1:(365 * 36)
temps = sapply(days, GetTemp)

df_temps <- data.frame(t = days/365, temp = temps)

df_temps |>
  ggplot(aes(x = t, y = temps)) +
  geom_line() +
  theme_bw()

# -- Run model
trajectories_amplilate <- read_rds("data/interim/seirsei_fig5-amplitude change late.rds")

# -- Plot outputs
year_introduction <- t_break/365
years_incratio <- c((year_introduction+1):(n_days_total/365))

trajectories_amplilate |>
  filter(M0 %in% c(1, 8) & (100-frac*100) %in% c(0, 50, 90)) |>
  filter(time/365 > year_introduction & time/365 <= year_introduction + 6) |>
  ggplot(aes(x = time/365, y = Ih, color = factor(100-frac*100))) +
  geom_line(lwd = 0.8) +
  scale_x_continuous(breaks = c(0:10)+year_introduction) +
  scale_color_brewer(palette = "Set1") +
  labs(x = "Year", y = "Infections",
       color = "%wMel at intervention") +
  facet_wrap(~M0, scales = "free_y", labeller = labeller(M0 = function(x) paste0("M0 = ", x)), ncol = 5) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.byrow = TRUE)

trajectories_amplilate |>
  filter((100-frac*100) %in% c(0, 50, 90)) |>
  filter(time/365 > year_introduction & time/365 <= year_introduction + 6) |>
  ggplot(aes(x = time/365, y = Ih, color = factor(100-frac*100))) +
  geom_line(lwd = 0.8) +
  scale_x_continuous(breaks = c(0:10)+year_introduction) +
  scale_color_brewer(palette = "Set1") +
  labs(x = "Year", y = "Infections",
       color = "%wMel at intervention") +
  facet_wrap(~M0, scales = "free_y", labeller = labeller(M0 = function(x) paste0("M0 = ", x)), ncol = 5) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.byrow = TRUE)
