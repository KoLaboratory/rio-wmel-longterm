# =============================================================================
# 02_empirical_analysis.R
# -----------------------------------------------------------------------------
# Empirical spatiotemporal analysis: INLA negative-binomial models of the effect
# of wMel introgression on dengue incidence, plus all main-text and supplementary
# figures/tables from the empirical arm of the study.
#
# Contents:
#   Figure 1   Descriptive maps, case time series, introgression by release area
#   Figure 2   Primary models  - Model 1  cumulative effect (2017-2024)
#                              - Model 1.5 any-vs-no introgression (cumulative)
#                              - Model 2   year-specific effects
#                              - Model 2.5 any-vs-no introgression (yearly)
#                              - g-computation standardized RI / prevented fraction
#   Table S3   Year-specific effects adjusting for past cumulative incidence
#   Figure 3   Pre-intervention risk & association with future introgression
#   Figure 4   Effect estimates under different historical spatial-risk adjustments
#   Fig S2/S3  Introgression tileplot; Rio municipality vs release-area comparison
#   Figure S4  Ovitrap-based sensitivity models (monthly / yearly / cumulative)
#   Table S2   Model comparison across case windows and exposure definitions
#
# Requires data/interim/prepared.rds (run analysis/01_prepare_data.R first).
# Figures are printed; ggsave() lines are left commented so a full run does not
# overwrite files. Uncomment (they already point at outputs/) to save.
# =============================================================================

source(file.path("R", "packages.R"))
source(file.path("R", "inla_helpers.R"))
source(file.path("R", "theme.R"))

# Load the shared prepared objects into this environment (skip if already present,
# e.g. when sourced straight after 01 in run_all.R).
if (!exists("grid")) {
  if (!file.exists("data/interim/prepared.rds")) {
    stop("data/interim/prepared.rds not found. Run analysis/01_prepare_data.R first.")
  }
  list2env(readRDS("data/interim/prepared.rds"), envir = environment())
}

# Re-read the WorldPop raster (terra rasters are not serialised into prepared.rds).
if (!exists("popdens")) {
  popdens <- rast("data/raw/wp_popdens_2020.tif") |>
    crop(release_areas_raw) |>
    project(crs(release_areas))
}

# ---- Figure 1 ----------------------------------------------------------

# -- Figure 1A (new): Pop density and location
br_map <- geobr::read_country(year = "2022")
rj_map <- geobr::read_municipality(code_muni = "RJ", year = "2022")

popdens_df <- as.data.frame(terra::mask(popdens, vect(release_areas_j)[1,]),
                            xy = TRUE) |>
  st_as_sf(coords = c("x", "y"), crs = st_crs(popdens)) |>
  mutate(geometry = st_buffer(geometry, endCapStyle = "SQUARE", dist = 50))

ext_br <- ext(st_buffer(br_map, dist = 1000))
ext_rj <- ext(st_buffer(rj_map, dist = 500))
ext_rj2 <- ext(st_buffer(rj_map |> filter(name_muni == "Rio de Janeiro"), dist = 500))
ext_rel <- ext(st_buffer(release_areas, dist = 500))

br_1km <- terra::rast("data/raw/bra_pd_2020_1km_UNadj.tif")

br_1km_p <- terra::project(br_1km, release_areas_j)

br_1km_mask <- terra::mask(br_1km_p, vect(release_areas_j))

br_1km_crop <- crop(br_1km_p, release_areas_j)

summary(br_1km_crop)

# get median household income
tracts_sf <- geobr::read_census_tract(
  code_tract = "RJ",
  simplified = FALSE,
  year = 2010,
  showProgress = FALSE
  )

tracts_rj <- tracts_sf |>
  filter(name_muni == "Rio De Janeiro") |>
  filter()

tracts_release <- st_join(release_areas |>
                            st_transform(crs = st_crs(tracts_rj)) |>
                            summarise(geometry = st_union(geometry)) |>
                            st_make_valid(), tracts_rj, join = st_contains)

tracts_release_full <- tracts_rj |>
  filter(code_tract %in% tracts_release$code_tract)

pacman::p_load(censobr)

renda_raw <- censobr::read_tracts(year = 2010, dataset = "DomicilioRenda")

domicilio_raw <- censobr::read_tracts(year = 2010, dataset = "Domicilio")

renda_rj <- renda_raw |>
  filter(code_tract %in% tracts_release_full$code_tract) |>
  collect() |>
  dplyr::select(code_tract, total_income = V002, total_income_pp = V003)

domicilio_rj <- domicilio_raw |>
  filter(code_tract %in% tracts_release_full$code_tract) |>
  collect() |>
  dplyr::select(code_tract, total_households = domicilio01_V001, total_households_pp = domicilio01_V002)

renda_dom <- renda_rj |>
  left_join(domicilio_rj, by = "code_tract") |>
  mutate(mean_hh_income = total_income / total_households,
         mean_hh_income_pp = total_income_pp / total_households_pp,
         mean_hh_income_usd = mean_hh_income * 0.2,
         mean_hh_income_usd_pp = mean_hh_income_pp * 0.2) |>
  left_join(tracts_release_full |> mutate(code_tract = as.character(code_tract)), by = "code_tract") |>
  st_as_sf()

summary(renda_dom$mean_hh_income_usd_pp)

renda_dom |>
  ggplot() +
  geom_sf(aes(fill = mean_hh_income_usd_pp)) +
  scale_fill_viridis_c() +
  theme_void()

(fig1n_a <- br_map |>
  ggplot() +
  geom_sf(fill = "gray30", color = NA) +
  geom_rect(aes(xmin = ext_rj[1],
                xmax = ext_rj[2],
                ymin = ext_rj[3],
                ymax = ext_rj[4]),
            fill = NA, color = "darkred", lwd = 1) +
  coord_sf(xlim = c(ext_br[1], ext_br[2]), ylim = c(ext_br[3], ext_br[4])) +
   scale_y_continuous(expand = c(0, 0)) +
   scale_x_continuous(expand = c(0, 0)) +
  theme_void() +
  theme(panel.background = element_rect(fill = "white"),
        plot.margin = margin(0, -12, 0, 2, "pt")))

(fig1n_b <- rj_map |>
  nngeo::st_remove_holes() |>
  ggplot() +
  geom_sf(fill = "gray30", color = "gray30") +
  geom_sf(data = rj_map |> filter(name_muni == "Rio De Janeiro"),
          fill = "lightgray", color = "lightgray") +
  geom_sf(data = release_areas, fill = "darkred", color = NA) +
  annotate("text", x = -43.53, y = -22.95, label = "Rio de Janeiro\nmunicipality", color = "gray30", fontface = "bold", size = 3) +
  annotate("text", x = -43.38, y = -22.87, label = "Release\narea", color = "darkred", fontface = "bold", size = 3) +
  ggspatial::annotation_scale(location = "bl",
                              pad_x = unit(2, "mm"),
                              text_cex = 1) +
  theme_void() +
  coord_sf(xlim = c(ext_rj2[1], ext_rj2[2]), ylim = c(ext_rj2[3], ext_rj2[4])) +
  theme(text = element_text(size = 16),
        legend.key.height = unit(12, "pt")))

(fig1n_ab <- fig1n_b + inset_element(fig1n_a, -0.15, 0.5, 0.3, 1, align_to = "full", ignore_tag = TRUE) & theme(plot.margin = margin(0, 0, 0, 0)))

(fig1n_c <- popdens_df |>
  ggplot() +
  geom_sf(aes(fill = population), color = NA) +
  scale_fill_distiller(palette = "Blues", trans = "log1p",
                       breaks = c(0, 15, 150, 1500), na.value = "transparent", direction = 1) +
  ggspatial::annotation_scale(location = "bl",
                              pad_x = unit(2, "mm"),
                              text_cex = 1) +
  labs(fill = "Population\n(100m)") +
  theme_void() +
  coord_sf(xlim = c(ext_rel[1], ext_rel[2]), ylim = c(ext_rel[3], ext_rel[4])) +
  theme(text = element_text(size = 16),
        legend.key.height = unit(12, "pt"),
        axis.title = element_blank()))

# -- Figure 1A: Maps of dengue incidence
colors <- c("#25489E", "#5482B4", "#EDC21D", "#E4697A", "#CE263D")
pal <- colorRampPalette(colors)(21)

fig1a_db <- data |>
    filter(pop > 30) |>
  group_by(year, id) |>
  summarise(cases = sum(cases, na.rm = TRUE),
            pop = mean(pop, na.rm = TRUE),
            .groups = "drop") |>
  mutate(inc = cases/pop*1000)

(fig1a <- fig1a_db |>
  mutate(inc = if_else(inc >55, 55, inc)) |>
  ggplot() +
  geom_sf(data = release_areas, fill = "gray", color = NA) +
  geom_sf(aes(fill = inc), color = NA) +
  scale_fill_distiller(palette = "Reds", trans = "log1p", breaks = c(0.5, 5, 50),
                       limits = c(0, NA), na.value = "transparent") +
  facet_wrap(~year, nrow = 2) +
  labs(fill = "Dengue\nincidence\nper 1000") +
  theme_void() +
  theme(text = element_text(size = 16),
        strip.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.position = "inside",
        legend.position.inside = c(0.95, 0.25),
        legend.key.height = unit(12, "pt"),
        axis.title = element_blank(),
        panel.spacing.x = unit(-2, "mm")))

(fig1a_n <- fig1a_db |>
  mutate(inc = if_else(inc >55, 55, inc)) |>
  ggplot() +
  geom_sf(data = release_areas, fill = "gray", color = NA) +
  geom_sf(aes(fill = inc), color = NA) +
  scale_fill_distiller(palette = "Reds", trans = "log1p", breaks = c(0.5, 5, 50),
                       limits = c(0, NA), direction = 1) +
  facet_wrap(~year, nrow = 3) +
  labs(fill = "Dengue\nincidence\nper 1000") +
  theme_void() +
  theme(text = element_text(size = 16),
        strip.text = element_text(size = 14, vjust = -2),
        strip.clip = "off",
        legend.title = element_text(size = 14),
        legend.key.height = unit(12, "pt"),
        axis.title = element_blank(),
        panel.spacing.x = unit(-5, "mm"),
        panel.spacing.y = unit(-3, "mm")))

# -- Figure 1B: Time series of cases
dengues <- dengue_grid |>
  group_by(year, month) |>
  summarise(cases = sum(cases)) |>
  mutate(date = as.Date(paste(year, month, "01", sep = "-"))) |>
  filter(date < as.Date("2024-05-01"))

dengues |>
  group_by(year) |>
  summarise(cases = sum(cases))

releaseses <- release_raw |>
  mutate(monthstart = lubridate::floor_date(date_of_release, "month")) |>
  group_by(monthstart) |>
  summarise(released = sum(average_released))

highcol <- "gray45"
lowcol.hex <- as.hexmode(round(col2rgb(highcol) * 0 + 255 * (1 - 0)))
lowcol <- paste0("#", paste(format(lowcol.hex, width = 2), collapse = ""))

(fig1b <- dengues |>
  ggplot(aes(x = date)) +
  geom_col(data = releaseses, aes(x = monthstart, y = 4700, fill = released/100000)) +
  scale_fill_distiller(palette = "Greys", name = "Mosquitoes\nreleased x10^5", direction = 1) +
  annotate("rect",
           xmin = as.Date("2017-09-01"),
           xmax = as.Date("2021-07-15"),
           ymin = 3200, ymax = 4000, fill = "white") +
  annotate("text", x = mean(release_raw$date_of_release),
           y = 3700, label = "wMel-infected Ae. aegypti releases", size = 5) +
  geom_line(aes(y = cases), size = 1.5) +
  labs(x = "", y = "Monthly\ndengue cases") +
  scale_x_date(date_labels = "%Y",
               breaks = seq.Date(as.Date("2010-01-01"), as.Date("2024-07-01"), by = "12 months"),
               expand = c(0.02, 0.02)) +
  theme_bw() +
  coord_cartesian(ylim = c(0, 4700)) +
  theme(axis.text.x = element_text(),
        axis.title.x = element_blank(),
        text = element_text(size = 16),
        plot.title = element_text(hjust = 0.5),
        legend.position = "inside",
        legend.position.inside = c(0.9, 0.55),
        legend.key.height = unit(2.9, "mm"),
        legend.title = element_text(size = 14),
        legend.background = element_blank()))

# ggsave("outputs/fig1b-may1.pdf", plot = fig1a)

(fig1b_n <- dengues |>
  ggplot(aes(x = date)) +
  geom_col(data = releaseses, aes(x = monthstart, y = 4700, fill = released/100000)) +
  scale_fill_distiller(palette = "Greys", name = "Mosquitoes\nreleased x10^5", direction = 1) +
  annotate("rect",
           xmin = as.Date("2017-09-01"),
           xmax = as.Date("2021-07-15"),
           ymin = 3200, ymax = 4000, fill = "white") +
  annotate("text", x = mean(release_raw$date_of_release),
           y = 3700, label = "wMel-infected Ae. aegypti releases", size = 5) +
  geom_line(aes(y = cases), size = 1.5) +
  labs(x = "", y = "Monthly\ndengue cases") +
  scale_x_date(date_labels = "%Y",
               breaks = seq.Date(as.Date("2010-01-01"), as.Date("2024-07-01"), by = "12 months"),
               expand = c(0.02, 0.02)) +
  theme_bw() +
  coord_cartesian(ylim = c(0, 4700)) +
  theme(axis.text.x = element_text(),
        axis.title.x = element_blank(),
        text = element_text(size = 16),
        plot.title = element_text(hjust = 0.5),
        legend.position = "inside",
        legend.position.inside = c(0.9, 0.55),
        legend.key.height = unit(2.9, "mm"),
        legend.title = element_text(size = 14),
        legend.background = element_blank()))

# -- Figure 1C: Introgression by release area
fig1c_db <- trap_raw |>
  filter(successful == TRUE &
           !is.na(target_species_count) &
           target_species_count > 0 &
           !is.na(screening_wmel_aeg)) |>
  mutate(week_start = floor_date(collected_at, unit = "week"),
         year = year(collected_at),
         intro = screening_wmel_aeg/target_species_count,
         zone = if_else(str_starts(reporting_area, "RJ3"), substr(reporting_area, 1, 5), substr(reporting_area, 1, 3)),
         zone = case_when(zone == "RJ3_1" ~ "RJ3.1",
                          zone == "RJ3_2" ~ "RJ3.2",
                          zone == "RJ3_3" ~ "RJ3.3",
                          .default = zone)) |>
  filter(zone != "TUB") |>
  mutate(study_week = as.numeric(as.Date(week_start) - as.Date("2017-09-24"))/7) |>
  group_by(trap_id) |>
  mutate(n_trap_events = n(),
         sum_tested = sum(target_species_count, na.rm = TRUE)) |>
  mutate(month_start = floor_date(collected_at, unit = "month")) |>
  group_by(zone, month_start, type) |>
  summarise(mean_intro = mean(intro, na.rm = TRUE),
            pos = sum(screening_wmel_aeg, na.rm = TRUE),
            tot = sum(target_species_count, na.rm = TRUE)) |>
  ungroup() |>
  mutate(mean_intro2 = pos/tot,
         trap = if_else(type == "BG", "BG-trap", "Ovitrap"),
         month_start = as.Date(month_start)) |>
  arrange(zone, type, month_start) |>
  group_by(zone, type) |>
  mutate(mean_intro3 = slide_dbl(mean_intro, mean, .before = 1, .after = 1))

(fig1c <- fig1c_db |>
  ggplot() +
  geom_col(data = release_dates, aes(x = monthstart, y = Inf, fill = "Releases"),
           width = 31) +
  scale_fill_manual(values = c("Releases" = "#BCCDDE")) +
  geom_line(aes(x = month_start, y = mean_intro*100, lty = trap), lwd = 1.05) +
  facet_wrap(~zone, ncol = 1, strip.position = "right") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y", expand = c(0.03, 0.03)) +
  scale_y_continuous(breaks = c(0, 50, 100)) +
  labs(x = "", y = "% wMel") +
  theme_bw() +
  theme(text = element_text(size = 16),
        legend.position = "inside",
        legend.position.inside = c(0.12, 0.15),
        legend.margin = margin(-0.8,0,0,0, "cm"),
        legend.spacing.y = unit(18, "mm"),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.frame = element_blank(),
        legend.box.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 12, margin = unit(c(-4, 0, -4, 2), "pt")),
        panel.spacing = unit(2, "mm"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 11),
        plot.margin = unit(c(3, 4, 0, 4), "pt")))

# -- Figure 1D: Introgression in 2019-2020
# trap point data
intro_2020 <- trap_raw |>
  filter(successful == TRUE &
           type == "BG" &
           !is.na(target_species_count) &
           target_species_count > 0 &
           !is.na(screening_wmel_aeg)) |>
  mutate(year = year(collected_at),
         intro = screening_wmel_aeg/target_species_count) |>
  filter(year%in% c(2019:2020)) |>
  group_by(latitude, longitude) |>
  summarise(mean_intro = mean(intro, na.rm = TRUE),
            pos = sum(screening_wmel_aeg, na.rm = TRUE),
            tot = sum(target_species_count, na.rm = TRUE)) |>
  ungroup() |>
  mutate(mean_intro2 = pos/tot) |>
  st_as_sf(coords = c("longitude", "latitude"), crs = "EPSG:4326 - WGS 84") |>
  st_transform(crs = st_crs(release_areas)) |>
  filter(st_intersects(geometry, release_areas |> st_union(), sparse = FALSE)[,1])

buffer_dist <- 500
ext_expanded <- extent(intro_2020)
ext_expanded@xmin <- ext_expanded@xmin - buffer_dist
ext_expanded@xmax <- ext_expanded@xmax + buffer_dist
ext_expanded@ymin <- ext_expanded@ymin - buffer_dist
ext_expanded@ymax <- ext_expanded@ymax + buffer_dist

grid2 <- raster(ext_expanded, resolution = c(500,500), crs = projection(intro_2020))

grid2$mean_intro <- rasterize(intro_2020, grid2, field = "mean_intro", fun = mean, na.rm = TRUE)

grid2$mean_introsm <- focal(grid2$mean_intro, w = matrix(1, nc = 3, nr = 3), fun = mean, na.rm = TRUE)

grid3db <- as.data.frame(grid2, xy = TRUE) |>
  st_as_sf(coords = c("x", "y"), crs = st_crs(release_areas)) |>
  filter(st_intersects(geometry, release_areas |> st_union(), sparse = FALSE)[,1]) |>
  mutate(geometry = st_buffer(geometry, dist = 250, endCapStyle = "SQUARE"))

(fig1d <- grid3db |>
  ggplot() +
  geom_sf(data = release_areas, fill = "lightgray", color = NA) +
  geom_sf(aes(fill = mean_intro*100), color = NA) +
  geom_sf(data = release_areas, fill = NA, color = "black", alpha = 0.5, lwd = 0.7) +
  geom_sf_label(data = dplyr::filter(release_areas, REGIAO_2 != "RJ3.1"),
                aes(label = REGIAO_2),
                linewidth = NA,
                color = "white",
                fill = "black",
                alpha = 0.5,
                size = 4.5,
                label.padding = unit(0.3, "mm"),
                fontface = "bold") +
  geom_sf_label(data = dplyr::filter(release_areas, REGIAO_2 == "RJ3.1"),
                aes(label = REGIAO_2),
                linewidth = NA,
                color = "white",
                fill = "black",
                alpha = 0.5,
                size = 4.5,
                label.padding = unit(0.3, "mm"),
                fontface = "bold",
                nudge_y = 1000) +
  scale_fill_distiller(palette = "Greens", na.value = "transparent", direction = 1) +
  labs(fill = "% wMel") +
  theme_void() +
  theme(text = element_text(size = 16),
        # legend.text = element_text(size = 11),
        legend.key.height = unit(12, "pt"),
        legend.title = element_text(size = 14),
        plot.margin = unit(c(8, 0, 40, -10), "pt")))

# -- Complete figure 1
(fig1 <- free(fig1a) / fig1b / (fig1d + fig1c + plot_layout(widths = c(6.3, 7.7))) +
   plot_annotation(tag_levels = "A") +
   plot_layout(heights = c(4.5, 2.5, 5.3)))

# ggsave("outputs/2026-04-09_fig1.pdf", plot = fig1, width = 14, height = 10.5)
# ggsave("outputs/2026-04-09_fig1.png", plot = fig1, width = 14, height = 10.5)

# -- Complete figure 1 with pop density
(fig1_new <- (((fig1n_ab / fig1n_c) | fig1a_n) + plot_layout(widths = c(2.9, 6.8))) /
   fig1b_n /
   (fig1d + fig1c + plot_layout(widths = c(6.3, 7.7))) +
   plot_layout(heights = c(2.5, 1, 2.5)) +
   plot_annotation(tag_levels = "A") &
   theme(plot.margin = margin(0, 0, 0, 0)))

# ggsave("outputs/2026-04-13_fig1-popdens.pdf", plot = fig1_new, width = 14, height = 10.5)

# ---- Figure 2 ----------------------------------------------------------

# -- Model 1: cumulative effect of introgression 2017-2020
dengue_grid_mod_long <- dengue |>
  st_drop_geometry() |>
  filter(!is.na(year) & !is.na(month)) |>
  group_by(year, month, id) |>
  summarise(cases = n()) |>
  ungroup() |>
  right_join(grid |> st_drop_geometry(), by = "id") |>
  tidyr::complete(year, month, id, fill = list(cases = 0)) |>
  filter(!is.na(year) & !is.na(month))

trap_long <- trap_raw |>
  filter(successful == "TRUE" & !is.na(target_species_count) & target_species_count > 0) |>
  mutate(intro = screening_wmel_aeg/target_species_count,
         year = year(collected_at),
         month = month(collected_at)) |>
  group_by(trap_id, type, latitude, longitude, year, month) |>
  summarise(mean_intro = mean(intro, na.rm = TRUE),
            .groups = "drop") |>
  st_as_sf(coords = c("longitude", "latitude"), crs = "EPSG:4326 - WGS 84", remove = FALSE) |>
  st_transform(crs = 31983)

trap_long2 <- st_join(grid, trap_long, join = st_contains)

trap_long_grid <- trap_long2 |>
  st_drop_geometry() |>
  filter(!is.na(trap_id)) |>
  group_by(id, year, month) |>
  summarise(mintro = mean(mean_intro, na.rm = TRUE),
            .groups = "drop")

data_mod_long <- grid |>
  full_join(dengue_grid_mod_long, by = "id") |>
  left_join(popdens_grid, by = "id") |>
  left_join(trap_long_grid, by = c("id", "year", "month")) |>
  left_join(release_grid, by = "id") |>
  rename(pop2 = sum) |>
  mutate(pop = if_else(pop2 < 50, 0, pop2),
         inc = cases/pop*1000,
         date_month = as.Date(paste(year, month, "01", sep = "-")),
         time_since_rel = date_month - month_first_release,
         mintro = if_else(is.na(mintro) & time_since_rel < 0 & year < 2018, 0, mintro),
         mintro = if_else(year < 2017, 0, mintro)) |>
  filter(st_intersects(geometry, release_areas_j, sparse = FALSE)[,1]) |>
  mutate(geometry = st_centroid(geometry),
         x = st_coordinates(geometry)[,1],
         y = st_coordinates(geometry)[,2]) |>
  st_drop_geometry() |>
  filter(pop != 0) |>
  mutate(logpop = log(pop),
         year = factor(year)) |>
  arrange(year) |>
  group_by(year) |>
  mutate(nyear = cur_group_id()) |>
  ungroup() |>
  filter(!is.na(mintro))

data_mod_long2 <- data_mod_long |>
  mutate(year = as.numeric(as.character(year))) |>
  filter(year >= 2017) |>
  mutate(year = factor(year)) |>
  arrange(year) |>
  group_by(year) |>
  mutate(nyear = cur_group_id()) |>
  ungroup() |>
  filter(!is.na(mintro))

sum(data_mod_long$cases[data_mod_long$year == 2017])
sum(data_mod_long$cases[data_mod_long$year == 2018])
sum(data_mod_long$cases[data_mod_long$year == 2019])
sum(data_mod_long$cases[data_mod_long$year == 2020])
sum(data_mod_long$cases[data_mod_long$year == 2021])
sum(data_mod_long$cases[data_mod_long$year == 2022])
sum(data_mod_long$cases[data_mod_long$year == 2023])
sum(data_mod_long$cases[data_mod_long$year == 2024])

dengue_grid_mod_long |>
  group_by(year) |>
  summarise(cases = sum(cases))

intro_year_summ <- trap_raw |>
  filter(successful == "TRUE" & !is.na(target_species_count) & target_species_count > 0) |>
  mutate(intro = screening_wmel_aeg/target_species_count,
         year = year(collected_at),
         month = month(collected_at),
         n_collected = target_species_count) |>
  group_by(trap_id, year) |>
  summarise(intro = weighted.mean(intro, w = n_collected),
            n_collected = sum(n_collected),
            .groups = "drop") |>
  group_by(year) |>
  summarise(mean_intro = mean(intro, na.rm = TRUE),
            weighted_mean = weighted.mean(intro, w = n_collected, na.rm = TRUE),
            median_intro = median(intro, na.rm = TRUE),
            q25 = quantile(intro, 0.25, na.rm = TRUE),
            q75 = quantile(intro, 0.75, na.rm = TRUE),
            .groups = "drop")

intro_year_summ

table(data_mod_long2$year)

spatial.range <- 22000

mesh <- fmesher::fm_mesh_2d_inla(
  loc = as.matrix(data_mod_long2[,c("x","y")]),
  max.edge = c(spatial.range / 5, spatial.range),
  cutoff = 200,
  offset = c(spatial.range / 5, diff(range(data_mod_long2[ ,c("x")])) / 10))

spde <- inla.spde2.matern(mesh = mesh)

s.index <- inla.spde.make.index(name = "spatial.field",
                                n.spde = spde$n.spde)

A.est = inla.spde.make.A(mesh = mesh,
                         loc = as.matrix(data_mod_long2[ ,c("x","y")]))

stack = inla.stack(data = list(obs = data_mod_long2$cases),
                   A = list(A.est, 1),
                   effects = list(spatial_field = data.frame(s.index, Intercept = 1),
                                  data_mod_long2[ ,c(names(data_mod_long2) != "cases")]),
                   tag = 'stdata')

modlong2b <- obs ~ -1 + Intercept + f(nyear, model = "ar1") + mintro + f(spatial.field, model = spde) + offset(logpop)

modlong2b_out <- inla(
  modlong2b,
  data = inla.stack.data(stack),
  family = "nbinomial",
  control.predictor = list(A = inla.stack.A(stack), compute = TRUE),
  control.fixed = list(expand.factor.strategy = 'inla'),
  control.compute = list(config = TRUE, dic = TRUE, waic = TRUE),
  verbose = TRUE
)

exp(as.data.frame(inla.zmarginal(modlong2b_out$marginals.fixed$mintro)))

# for 10 percentage point change

# marginal of the mintro coefficient (beta_1), on the log scale
m_beta <- modlong2b_out$marginals.fixed$mintro

# RI for a +10 percentage-point (0.10 on the 0-1 proportion scale) increase
m_RI10 <- inla.tmarginal(function(x) exp(0.1 * x), m_beta)

# posterior mean, sd, and quantiles (incl. 95% CrI)
inla.zmarginal(m_RI10)

mintro_marg <- modlong2b_out$marginals.fixed$mintro

set.seed(123)
beta_samples <- inla.rmarginal(1000, mintro_marg)

mintro_values <- seq(0, 1, length.out = 100)

effect_samples <- sapply(beta_samples, function(b) exp(b * mintro_values))

modlong2b_df <- data.frame(
  mintro = mintro_values,
  effect = apply(effect_samples, 1, mean),
  lower  = apply(effect_samples, 1, quantile, 0.025),
  upper  = apply(effect_samples, 1, quantile, 0.975)
)

modlong2b_df |>
  ggplot(aes(x = mintro, y = effect)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_hline(aes(yintercept = 1), lty = 2, color = "red") +
  labs(x = "% wMel introgression", y = "Ratio of incidences") +
  theme_bw() +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1")) +
  coord_trans(y = "log1p") +
  theme(legend.position = "bottom",
        text = element_text(size = 13))

## spatial field
coords_modlong2b <- spde$mesh$loc

xrange_modlong2b <- seq(min(coords_modlong2b[, 1]), max(coords_modlong2b[, 1]), length.out = 100)
yrange_modlong2b <- seq(min(coords_modlong2b[, 2]), max(coords_modlong2b[, 2]), length.out = 100)

grid_modlong2b <- expand.grid(x = xrange_modlong2b, y = yrange_modlong2b)
sp::coordinates(grid_modlong2b) <- ~x + y
proj_grid_modlong2b <- inla.mesh.projector(spde$mesh, xlim = range(xrange_modlong2b), ylim = range(yrange_modlong2b), dims = c(100, 100))

spatial_mean_modlong2b <- inla.mesh.project(proj_grid_modlong2b, modlong2b_out$summary.random$spatial.field$mean)
spatial_lb_modlong2b <- inla.mesh.project(proj_grid_modlong2b, modlong2b_out$summary.random$spatial.field$`0.025quant`)
spatial_ub_modlong2b <- inla.mesh.project(proj_grid_modlong2b, modlong2b_out$summary.random$spatial.field$`0.975quant`)

grid_modlong2b$mean <- as.vector(spatial_mean_modlong2b)
grid_modlong2b$lb <- as.vector(spatial_lb_modlong2b)
grid_modlong2b$ub <- as.vector(spatial_ub_modlong2b)

grid_sf_modlong2b <- as.data.frame(grid_modlong2b) |>
  st_as_sf(coords = c("x", "y"), remove = FALSE, crs = 31983) |>
  mutate(sig = case_when(lb > 0 & ub > 0 ~ 1,
                         lb < 0 & ub < 0 ~ 1,
                         .default = 0)) |>
  st_intersection(release_areas_j)

# -- g-computation nfor population-standardized ratio of total predicted cases (realized impact given incomplete, spatially heterogeneous introgression that actually occurred)
set.seed(1)
S    <- 2000
samp <- inla.posterior.sample(S, modlong2b_out)

idx  <- inla.stack.index(stack, tag = "stdata")$data
mintro_obs <- data_mod_long2$mintro

ln         <- rownames(samp[[1]]$latent)
apred_rows <- grep("^APredictor", ln)
beta_row   <- grep("^mintro", ln)

ri_one <- function(s) {
  eta_obs <- s$latent[apred_rows, 1][idx]
  beta    <- s$latent[beta_row, 1]
  mu_obs  <- exp(eta_obs)
  mu_0    <- exp(eta_obs - beta * mintro_obs)
  sum(mu_obs) / sum(mu_0)
}

ri_draws <- vapply(samp, ri_one, numeric(1))

# standardized cumulative RI + 95% CrI
c(mean = mean(ri_draws), median = median(ri_draws),
  quantile(ri_draws, c(0.025, 0.975)))

# prevented fraction (cases averted)
pf <- 1 - ri_draws
c(mean = mean(pf), quantile(pf, c(0.025, 0.975)))

# -- Model 1.5: yes/no introgression
data_mod_long2b <- data_mod_long2 |>
  mutate(mintro = if_else(mintro > 0, 1, 0))

spatial.range <- 22000

mesh <- fmesher::fm_mesh_2d_inla(
  loc = as.matrix(data_mod_long2b[,c("x","y")]),
  max.edge = c(spatial.range / 5, spatial.range),
  cutoff = 200,
  offset = c(spatial.range / 5, diff(range(data_mod_long2b[ ,c("x")])) / 10))

spde <- inla.spde2.matern(mesh = mesh)

s.index <- inla.spde.make.index(name = "spatial.field",
                                n.spde = spde$n.spde)

A.est = inla.spde.make.A(mesh = mesh,
                         loc = as.matrix(data_mod_long2b[ ,c("x","y")]))

modlong2bb <- obs ~ -1 + Intercept + f(nyear, model = "ar1") + mintro + f(spatial.field, model = spde) + offset(logpop)

stack = inla.stack(data = list(obs = data_mod_long2b$cases),
                   A = list(A.est, 1),
                   effects = list(spatial_field = data.frame(s.index, Intercept = 1),
                                  data_mod_long2b[ ,c(names(data_mod_long2b) != "cases")]),
                   tag = 'stdata')

modlong2bb_out <- inla(
  modlong2bb,
  data = inla.stack.data(stack),
  family = "nbinomial",
  control.predictor = list(A = inla.stack.A(stack), compute = TRUE),
  control.fixed = list(expand.factor.strategy = 'inla'),
  control.compute = list(config = TRUE, dic = TRUE, waic = TRUE),
  verbose = TRUE
)

exp(as.data.frame(inla.zmarginal(modlong2bb_out$marginals.fixed$mintro)))

# -- Model 2: yearly effect of introgression 2017-2020
spatial.range <- 22000

mesh <- inla.mesh.2d(loc = as.matrix(data_mod_long2[,c("x","y")]),
                     max.edge = c(spatial.range / 5, spatial.range),
                     cutoff = 200,
                     offset = c(spatial.range / 5, diff(range(data_mod_long2[ ,c("x")])) / 10))

spde <- inla.spde2.matern(mesh = mesh)

s.index <- inla.spde.make.index(name = "spatial.field",
                                n.spde = spde$n.spde)

A.est = inla.spde.make.A(mesh = mesh,
                         loc = as.matrix(data_mod_long2[ ,c("x","y")]))

stack = inla.stack(data = list(obs = data_mod_long2$cases),
                   A = list(A.est, 1),
                   effects = list(spatial_field = data.frame(s.index, Intercept = 1),
                                  data_mod_long2[ ,c(names(data_mod_long2) != "cases")]),
                   tag = 'stdata')

modlong2 <- obs ~ -1 + Intercept + f(nyear, model = "ar1") + year:mintro + f(spatial.field, model = spde) + offset(logpop)

modlong2_out <- inla(
  modlong2,
  data = inla.stack.data(stack),
  family = "nbinomial",
  control.predictor = list(A = inla.stack.A(stack), compute = TRUE),
  control.fixed = list(expand.factor.strategy = 'inla'),
  control.compute = list(config = TRUE, dic = TRUE, waic = TRUE),
  verbose = TRUE
)

round(exp(as.data.frame(inla.zmarginal(modlong2_out$marginals.fixed$`year2017:mintro`))), 2)
round(exp(as.data.frame(inla.zmarginal(modlong2_out$marginals.fixed$`year2018:mintro`))), 2)
round(exp(as.data.frame(inla.zmarginal(modlong2_out$marginals.fixed$`year2019:mintro`))), 2)
round(exp(as.data.frame(inla.zmarginal(modlong2_out$marginals.fixed$`year2020:mintro`))), 2)
round(exp(as.data.frame(inla.zmarginal(modlong2_out$marginals.fixed$`year2021:mintro`))), 2)
round(exp(as.data.frame(inla.zmarginal(modlong2_out$marginals.fixed$`year2022:mintro`))), 2)
round(exp(as.data.frame(inla.zmarginal(modlong2_out$marginals.fixed$`year2023:mintro`))), 2)
round(exp(as.data.frame(inla.zmarginal(modlong2_out$marginals.fixed$`year2024:mintro`))), 2)

# for 10 percentage point change

# per-10-percentage-point RI for one year
per10 <- function(marg) {
  m <- inla.tmarginal(function(x) exp(0.1 * x), marg)
  round(unlist(inla.zmarginal(m, silent = TRUE)), 3)
}

years <- 2017:2024
res10 <- t(sapply(years, function(y) {
  per10(modlong2_out$marginals.fixed[[paste0("year", y, ":mintro")]])
}))
rownames(res10) <- years
res10

mintro_int_margs_modlong2 <- list(
  `2017` = modlong2_out$marginals.fixed$`year2017:mintro`,
  `2018` = modlong2_out$marginals.fixed$`year2018:mintro`,
  `2019` = modlong2_out$marginals.fixed$`year2019:mintro`,
  `2020` = modlong2_out$marginals.fixed$`year2020:mintro`,
  `2021` = modlong2_out$marginals.fixed$`year2021:mintro`,
  `2022` = modlong2_out$marginals.fixed$`year2022:mintro`,
  `2023` = modlong2_out$marginals.fixed$`year2023:mintro`,
  `2024` = modlong2_out$marginals.fixed$`year2024:mintro`
)

mintro_values <- seq(0, 1, length.out = 100)

set.seed(123)
modlong2_df <- bind_rows(lapply(names(mintro_int_margs_modlong2), function(yr) {
  beta_year_draws <- inla.rmarginal(1000, mintro_int_margs_modlong2[[yr]])
  effect_samples <- sapply(beta_year_draws, function(b) exp(b * mintro_values))

  data.frame(
    year   = as.integer(yr),
    mintro = mintro_values,
    effect = apply(effect_samples, 1, mean),
    lower  = apply(effect_samples, 1, quantile, probs = 0.025),
    upper  = apply(effect_samples, 1, quantile, probs = 0.975)
  )
}))

modlong2_df |>
  ggplot(aes(x = mintro, y = effect)) +
  geom_line(lwd = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_hline(aes(yintercept = 1), lty = 2, color = "red") +
  labs(x = "% wMel introgression", y = "Ratio of incidences", color = "Year", fill = "Year") +
  theme_bw() +
  facet_wrap(~year, ncol = 4) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1")) +
  coord_transform(y = "log1p") +
  theme(legend.position = "bottom",
        text = element_text(size = 15))

## spatial field
coords_modlong2 <- spde$mesh$loc

xrange_modlong2 <- seq(min(coords_modlong2[, 1]), max(coords_modlong2[, 1]), length.out = 100)
yrange_modlong2 <- seq(min(coords_modlong2[, 2]), max(coords_modlong2[, 2]), length.out = 100)

grid_modlong2 <- expand.grid(x = xrange_modlong2, y = yrange_modlong2)
sp::coordinates(grid_modlong2) <- ~x + y
proj_grid_modlong2 <- inla.mesh.projector(spde$mesh, xlim = range(xrange_modlong2), ylim = range(yrange_modlong2), dims = c(100, 100))

spatial_mean_modlong2 <- inla.mesh.project(proj_grid_modlong2, modlong2_out$summary.random$spatial.field$mean)
spatial_lb_modlong2 <- inla.mesh.project(proj_grid_modlong2, modlong2_out$summary.random$spatial.field$`0.025quant`)
spatial_ub_modlong2 <- inla.mesh.project(proj_grid_modlong2, modlong2_out$summary.random$spatial.field$`0.975quant`)

grid_modlong2$mean <- as.vector(spatial_mean_modlong2)
grid_modlong2$lb <- as.vector(spatial_lb_modlong2)
grid_modlong2$ub <- as.vector(spatial_ub_modlong2)

grid_sf_modlong2 <- as.data.frame(grid_modlong2) |>
  st_as_sf(coords = c("x", "y"), remove = FALSE, crs = 31983) |>
  mutate(sig = case_when(lb > 0 & ub > 0 ~ 1,
                         lb < 0 & ub < 0 ~ 1,
                         .default = 0)) |>
  st_intersection(release_areas_j)

# -- Model 2.5: yes/no introgression
spatial.range <- 22000

mesh <- inla.mesh.2d(loc = as.matrix(data_mod_long2b[,c("x","y")]),
                     max.edge = c(spatial.range / 5, spatial.range),
                     cutoff = 200,
                     offset = c(spatial.range / 5, diff(range(data_mod_long2b[ ,c("x")])) / 10))

spde <- inla.spde2.matern(mesh = mesh)

s.index <- inla.spde.make.index(name = "spatial.field",
                                n.spde = spde$n.spde)

A.est = inla.spde.make.A(mesh = mesh,
                         loc = as.matrix(data_mod_long2b[ ,c("x","y")]))

stack = inla.stack(data = list(obs = data_mod_long2b$cases),
                   A = list(A.est, 1),
                   effects = list(spatial_field = data.frame(s.index, Intercept = 1),
                                  data_mod_long2b[ ,c(names(data_mod_long2b) != "cases")]),
                   tag = 'stdata')

modlong2b <- obs ~ -1 + Intercept + f(nyear, model = "ar1") + year:mintro + f(spatial.field, model = spde) + offset(logpop)

modlong25_out <- inla(
  modlong2b,
  data = inla.stack.data(stack),
  family = "nbinomial",
  control.predictor = list(A = inla.stack.A(stack), compute = TRUE),
  control.fixed = list(expand.factor.strategy = 'inla'),
  control.compute = list(config = TRUE, dic = TRUE, waic = TRUE),
  verbose = TRUE
)

round(exp(as.data.frame(inla.zmarginal(modlong25_out$marginals.fixed$`year2017:mintro`))), 2)
round(exp(as.data.frame(inla.zmarginal(modlong25_out$marginals.fixed$`year2018:mintro`))), 2)
round(exp(as.data.frame(inla.zmarginal(modlong25_out$marginals.fixed$`year2019:mintro`))), 2)
round(exp(as.data.frame(inla.zmarginal(modlong25_out$marginals.fixed$`year2020:mintro`))), 2)
round(exp(as.data.frame(inla.zmarginal(modlong25_out$marginals.fixed$`year2021:mintro`))), 2)
round(exp(as.data.frame(inla.zmarginal(modlong25_out$marginals.fixed$`year2022:mintro`))), 2)
round(exp(as.data.frame(inla.zmarginal(modlong25_out$marginals.fixed$`year2023:mintro`))), 2)
round(exp(as.data.frame(inla.zmarginal(modlong25_out$marginals.fixed$`year2024:mintro`))), 2)

# -- Complete figure 2
(fig2a <- modlong2b_df |>
  mutate(year = "Cumulative") |>
  ggplot(aes(x = mintro, y = effect)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_hline(aes(yintercept = 1), lty = 2, color = "red") +
  labs(x = "% wMel introgression", y = "Ratio of incidences") +
  theme_bw() +
  facet_wrap(~year, ncol = 1) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1")) +
  scale_y_continuous(limits = c(0, 2.46)) +
  coord_trans(y = "log1p") +
  theme(text = element_text(size = 15),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()))

(fig2b <- modlong2_df |>
  ggplot(aes(x = mintro, y = effect)) +
  geom_line(lwd = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_hline(aes(yintercept = 1), lty = 2, color = "red") +
  labs(x = "% wMel introgression", y = "Ratio of incidences") +
  theme_bw() +
  facet_wrap(~year, nrow = 1) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1")) +
  coord_trans(y = "log1p") +
  theme(text = element_text(size = 15))) +
  coord_cartesian(ylim = c(0, 2.46))

(fig2 <- fig2b + fig2a + plot_layout(widths = c(8, 1), axis_titles = "collect")  +
   plot_annotation(tag_levels = "A",
                  theme = theme(
                    plot.margin = margin(1, 3, 1, 3),
                    plot.tag.position = c(0, 0))))

# ggsave("outputs/2026-07-06_fig2.png", plot = fig2, height = 3, width = 14)
# ggsave("outputs/2026-07-06_fig2.pdf", plot = fig2, height = 3, width = 14)

# -- Model X: introgression-year interaction adjusting for past incidence
dengue_grid_mod11 <- dengue |>
  st_drop_geometry() |>
  filter(!is.na(year)) |>
  group_by(year, id) |>
  summarise(cases = n()) |>
  ungroup() |>
  right_join(grid |> st_drop_geometry(), by = "id") |>
  tidyr::complete(year, id, fill = list(cases = 0)) |>
  filter(!is.na(year))

data_pastinc <- grid |>
  right_join(dengue_grid_mod_long, by = "id") |>
  left_join(popdens_grid, by = "id") |>
  left_join(trap_long_grid, by = c("id", "year", "month")) |>
  left_join(release_grid, by = "id") |>
  rename(pop2 = sum) |>
  mutate(pop = if_else(pop2 < 50, 0, pop2),
         inc = cases/pop*1000,
         date_month = as.Date(paste(year, month, "01", sep = "-")),
         time_since_rel = date_month - month_first_release,
         mintro = if_else(is.na(mintro) & time_since_rel < 0 & year < 2018, 0, mintro),
         mintro = if_else(year < 2017, 0, mintro)) |>
  filter(st_intersects(geometry, release_areas_j, sparse = FALSE)[,1]) |>
  mutate(geometry = st_centroid(geometry),
         x = st_coordinates(geometry)[,1],
         y = st_coordinates(geometry)[,2]) |>
  st_drop_geometry() |>
  filter(pop != 0) |>
  arrange(id, year, month) |>
  group_by(id) |>
  mutate(across(cases,
                list(inc_1y = \(x) slide_sum(x, before = 1*12, after = 0, complete = FALSE) - x,
                     inc_2y = \(x) slide_sum(x, before = 2*12, after = 0, complete = FALSE) - x,
                     inc_3y = \(x) slide_sum(x, before = 3*12, after = 0, complete = FALSE) - x,
                     inc_4y = \(x) slide_sum(x, before = 4*12, after = 0, complete = FALSE) - x,
                     inc_5y = \(x) slide_sum(x, before = 5*12, after = 0, complete = FALSE) - x,
                     inc_6y = \(x) slide_sum(x, before = 6*12, after = 0, complete = FALSE) - x,
                     inc_7y = \(x) slide_sum(x, before = 7*12, after = 0, complete = FALSE) - x),
                .names = "{.fn}")) |>
  mutate(across(inc_1y:inc_7y, \(x) x / pop * 1000))

data_mod11 <- data_pastinc |>
  ungroup() |>
  filter(year >= 2017) |>
  mutate(logpop = log(pop),
         year = factor(year)) |>
  arrange(year) |>
  group_by(year) |>
  mutate(nyear = cur_group_id()) |>
  ungroup() |>
  filter(!is.na(mintro)) |>
  mutate(inc_1y_binned = inla.group(inc_1y, n = 5, method = "quantile"),
         inc_2y_binned = inla.group(inc_2y, n = 5, method = "quantile"),
         inc_3y_binned = inla.group(inc_3y, n = 5, method = "quantile"),
         inc_4y_binned = inla.group(inc_4y, n = 5, method = "quantile"),
         inc_5y_binned = inla.group(inc_5y, n = 5, method = "quantile"),
         inc_6y_binned = inla.group(inc_6y, n = 5, method = "quantile"),
         inc_7y_binned = inla.group(inc_7y, n = 5, method = "quantile"))

spatial.range <- 22000

mesh <- inla.mesh.2d(loc = as.matrix(data_mod11[,c("x","y")]),
                     max.edge = c(spatial.range / 5, spatial.range),
                     cutoff = 200,
                     offset = c(spatial.range / 5, diff(range(data_mod11[ ,c("x")])) / 10))

spde <- inla.spde2.matern(mesh = mesh)

s.index <- inla.spde.make.index(name = "spatial.field",
                                n.spde = spde$n.spde)

A.est = inla.spde.make.A(mesh = mesh,
                         loc = as.matrix(data_mod11[ ,c("x","y")]))

stack = inla.stack(data = list(obs = data_mod11$cases),
                   A = list(A.est, 1),
                   effects = list(spatial_field = data.frame(s.index, Intercept = 1),
                                  data_mod11[ ,c(names(data_mod11) != "cases")]),
                   tag = 'stdata')

mod29 <- obs ~ -1 + Intercept + f(nyear, model = "ar1") + year:mintro + f(spatial.field, model = spde) + inc_1y + offset(logpop)

mod11 <- obs ~ -1 + Intercept + f(nyear, model = "ar1") + year:mintro + f(spatial.field, model = spde) + inc_2y + offset(logpop)

mod12 <- obs ~ -1 + Intercept + f(nyear, model = "ar1") + year:mintro + f(spatial.field, model = spde) + inc_3y + offset(logpop)

mod13 <- obs ~ -1 + Intercept + f(nyear, model = "ar1") + year:mintro + f(spatial.field, model = spde) + inc_4y + offset(logpop)

mod14 <- obs ~ -1 + Intercept + f(nyear, model = "ar1") + year:mintro + f(spatial.field, model = spde) + inc_5y + offset(logpop)

mod15 <- obs ~ -1 + Intercept + f(nyear, model = "ar1") + year:mintro + f(spatial.field, model = spde) + inc_6y + offset(logpop)

mod16 <- obs ~ -1 + Intercept + f(nyear, model = "ar1") + year:mintro + f(spatial.field, model = spde) + inc_7y + offset(logpop)

mod29_out <- inla(
  mod29,
  data = inla.stack.data(stack),
  family = "nbinomial",
  control.predictor = list(A = inla.stack.A(stack), compute = TRUE),
  control.fixed = list(expand.factor.strategy = 'inla'),
  control.compute = list(config = TRUE, dic = TRUE, waic = TRUE),
  verbose = TRUE
)

mod11_out <- inla(
  mod11,
  data = inla.stack.data(stack),
  family = "nbinomial",
  control.predictor = list(A = inla.stack.A(stack), compute = TRUE),
  control.fixed = list(expand.factor.strategy = 'inla'),
  control.compute = list(config = TRUE, dic = TRUE, waic = TRUE),
  verbose = TRUE
)

mod12_out <- inla(
  mod12,
  data = inla.stack.data(stack),
  family = "nbinomial",
  control.predictor = list(A = inla.stack.A(stack), compute = TRUE),
  control.fixed = list(expand.factor.strategy = 'inla'),
  control.compute = list(config = TRUE, dic = TRUE, waic = TRUE),
  verbose = TRUE
)

mod13_out <- inla(
  mod13,
  data = inla.stack.data(stack),
  family = "nbinomial",
  control.predictor = list(A = inla.stack.A(stack), compute = TRUE),
  control.fixed = list(expand.factor.strategy = 'inla'),
  control.compute = list(config = TRUE, dic = TRUE, waic = TRUE),
  verbose = TRUE
)

mod14_out <- inla(
  mod14,
  data = inla.stack.data(stack),
  family = "nbinomial",
  control.predictor = list(A = inla.stack.A(stack), compute = TRUE),
  control.fixed = list(expand.factor.strategy = 'inla'),
  control.compute = list(config = TRUE, dic = TRUE, waic = TRUE),
  verbose = TRUE
)

mod15_out <- inla(
  mod15,
  data = inla.stack.data(stack),
  family = "nbinomial",
  control.predictor = list(A = inla.stack.A(stack), compute = TRUE),
  control.fixed = list(expand.factor.strategy = 'inla'),
  control.compute = list(config = TRUE, dic = TRUE, waic = TRUE),
  verbose = TRUE
)

mod16_out <- inla(
  mod16,
  data = inla.stack.data(stack),
  family = "nbinomial",
  control.predictor = list(A = inla.stack.A(stack), compute = TRUE),
  control.fixed = list(expand.factor.strategy = 'inla'),
  control.compute = list(config = TRUE, dic = TRUE, waic = TRUE),
  verbose = TRUE
)

n_levels_mod30 <- nrow(data_mod11 |> dplyr::select(inc_1y_binned) |> distinct())
extraconstr_mod30 <- list(A = matrix(0, ncol = n_levels_mod30, nrow = 1),
                          e = matrix(0, ncol = 1))
extraconstr_mod30$A[1, 1] <- 1

n_levels_mod17 <- nrow(data_mod11 |> dplyr::select(inc_2y_binned) |> distinct())
extraconstr_mod17 <- list(A = matrix(0, ncol = n_levels_mod17, nrow = 1),
                          e = matrix(0, ncol = 1))
extraconstr_mod17$A[1, 1] <- 1

n_levels_mod18 <- nrow(data_mod11 |> dplyr::select(inc_3y_binned) |> distinct())
extraconstr_mod18 <- list(A = matrix(0, ncol = n_levels_mod18, nrow = 1),
                          e = matrix(0, ncol = 1))
extraconstr_mod18$A[1, 1] <- 1

n_levels_mod19 <- nrow(data_mod11 |> dplyr::select(inc_4y_binned) |> distinct())
extraconstr_mod19 <- list(A = matrix(0, ncol = n_levels_mod19, nrow = 1),
                          e = matrix(0, ncol = 1))
extraconstr_mod19$A[1, 1] <- 1

n_levels_mod20 <- nrow(data_mod11 |> dplyr::select(inc_5y_binned) |> distinct())
extraconstr_mod20 <- list(A = matrix(0, ncol = n_levels_mod20, nrow = 1),
                          e = matrix(0, ncol = 1))
extraconstr_mod20$A[1, 1] <- 1

n_levels_mod21 <- nrow(data_mod11 |> dplyr::select(inc_6y_binned) |> distinct())
extraconstr_mod21 <- list(A = matrix(0, ncol = n_levels_mod21, nrow = 1),
                          e = matrix(0, ncol = 1))
extraconstr_mod21$A[1, 1] <- 1

n_levels_mod22 <- nrow(data_mod11 |> dplyr::select(inc_7y_binned) |> distinct())
extraconstr_mod22 <- list(A = matrix(0, ncol = n_levels_mod22, nrow = 1),
                          e = matrix(0, ncol = 1))
extraconstr_mod22$A[1, 1] <- 1

mod30 <- obs ~ -1 + Intercept + f(nyear, model = "ar1") + year:mintro + f(spatial.field, model = spde) +
  f(inc_1y_binned, model = "rw1", constr = FALSE, extraconstr = extraconstr_mod30) + offset(logpop)

mod17 <- obs ~ -1 + Intercept + f(nyear, model = "ar1") + year:mintro + f(spatial.field, model = spde) +
  f(inc_2y_binned, model = "rw1", constr = FALSE, extraconstr = extraconstr_mod17) + offset(logpop)

mod18 <- obs ~ -1 + Intercept + f(nyear, model = "ar1") + year:mintro + f(spatial.field, model = spde) +
  f(inc_3y_binned, model = "rw1", constr = FALSE, extraconstr = extraconstr_mod18) + offset(logpop)

mod19 <- obs ~ -1 + Intercept + f(nyear, model = "ar1") + year:mintro + f(spatial.field, model = spde) +
  f(inc_4y_binned, model = "rw1", constr = FALSE, extraconstr = extraconstr_mod19) + offset(logpop)

mod20 <- obs ~ -1 + Intercept + f(nyear, model = "ar1") + year:mintro + f(spatial.field, model = spde) +
  f(inc_5y_binned, model = "rw1", constr = FALSE, extraconstr = extraconstr_mod20) + offset(logpop)

mod21 <- obs ~ -1 + Intercept + f(nyear, model = "ar1") + year:mintro + f(spatial.field, model = spde) +
  f(inc_6y_binned, model = "rw1", constr = FALSE, extraconstr = extraconstr_mod21) + offset(logpop)

mod22 <- obs ~ -1 + Intercept + f(nyear, model = "ar1") + year:mintro + f(spatial.field, model = spde) +
  f(inc_7y_binned, model = "rw1", constr = FALSE, extraconstr = extraconstr_mod22) + offset(logpop)

mod30_out <- inla(
  mod30,
  data = inla.stack.data(stack),
  family = "nbinomial",
  control.predictor = list(A = inla.stack.A(stack), compute = TRUE),
  control.fixed = list(expand.factor.strategy = 'inla'),
  control.compute = list(config = TRUE, dic = TRUE, waic = TRUE),
  verbose = TRUE
)

mod17_out <- inla(
  mod17,
  data = inla.stack.data(stack),
  family = "nbinomial",
  control.predictor = list(A = inla.stack.A(stack), compute = TRUE),
  control.fixed = list(expand.factor.strategy = 'inla'),
  control.compute = list(config = TRUE, dic = TRUE, waic = TRUE),
  verbose = TRUE
)

mod18_out <- inla(
  mod18,
  data = inla.stack.data(stack),
  family = "nbinomial",
  control.predictor = list(A = inla.stack.A(stack), compute = TRUE),
  control.fixed = list(expand.factor.strategy = 'inla'),
  control.compute = list(config = TRUE, dic = TRUE, waic = TRUE),
  verbose = TRUE
)

mod19_out <- inla(
  mod19,
  data = inla.stack.data(stack),
  family = "nbinomial",
  control.predictor = list(A = inla.stack.A(stack), compute = TRUE),
  control.fixed = list(expand.factor.strategy = 'inla'),
  control.compute = list(config = TRUE, dic = TRUE, waic = TRUE),
  verbose = TRUE
)

mod20_out <- inla(
  mod20,
  data = inla.stack.data(stack),
  family = "nbinomial",
  control.predictor = list(A = inla.stack.A(stack), compute = TRUE),
  control.fixed = list(expand.factor.strategy = 'inla'),
  control.compute = list(config = TRUE, dic = TRUE, waic = TRUE),
  verbose = TRUE
)

mod21_out <- inla(
  mod21,
  data = inla.stack.data(stack),
  family = "nbinomial",
  control.predictor = list(A = inla.stack.A(stack), compute = TRUE),
  control.fixed = list(expand.factor.strategy = 'inla'),
  control.compute = list(config = TRUE, dic = TRUE, waic = TRUE),
  verbose = TRUE
)

mod22_out <- inla(
  mod22,
  data = inla.stack.data(stack),
  family = "nbinomial",
  control.predictor = list(A = inla.stack.A(stack), compute = TRUE),
  control.fixed = list(expand.factor.strategy = 'inla'),
  control.compute = list(config = TRUE, dic = TRUE, waic = TRUE),
  verbose = TRUE
)

mod31 <- obs ~ -1 + Intercept + f(nyear, model = "ar1") + year:mintro + f(spatial.field, model = spde) +
  f(inc_1y_binned, model = "iid", constr = FALSE, extraconstr = extraconstr_mod30) + offset(logpop)

mod23 <- obs ~ -1 + Intercept + f(nyear, model = "ar1") + year:mintro + f(spatial.field, model = spde) +
  f(inc_2y_binned, model = "iid", constr = FALSE, extraconstr = extraconstr_mod17) + offset(logpop)

mod24 <- obs ~ -1 + Intercept + f(nyear, model = "ar1") + year:mintro + f(spatial.field, model = spde) +
  f(inc_3y_binned, model = "iid", constr = FALSE, extraconstr = extraconstr_mod18) + offset(logpop)

mod25 <- obs ~ -1 + Intercept + f(nyear, model = "ar1") + year:mintro + f(spatial.field, model = spde) +
  f(inc_4y_binned, model = "iid", constr = FALSE, extraconstr = extraconstr_mod19) + offset(logpop)

mod26 <- obs ~ -1 + Intercept + f(nyear, model = "ar1") + year:mintro + f(spatial.field, model = spde) +
  f(inc_5y_binned, model = "iid", constr = FALSE, extraconstr = extraconstr_mod20) + offset(logpop)

mod27 <- obs ~ -1 + Intercept + f(nyear, model = "ar1") + year:mintro + f(spatial.field, model = spde) +
  f(inc_6y_binned, model = "iid", constr = FALSE, extraconstr = extraconstr_mod21) + offset(logpop)

mod28 <- obs ~ -1 + Intercept + f(nyear, model = "ar1") + year:mintro + f(spatial.field, model = spde) +
  f(inc_7y_binned, model = "iid", constr = FALSE, extraconstr = extraconstr_mod22) + offset(logpop)

mod31_out <- inla(
  mod31,
  data = inla.stack.data(stack),
  family = "nbinomial",
  control.predictor = list(A = inla.stack.A(stack), compute = TRUE),
  control.fixed = list(expand.factor.strategy = 'inla'),
  control.compute = list(config = TRUE, dic = TRUE, waic = TRUE),
  verbose = TRUE
)

mod23_out <- inla(
  mod23,
  data = inla.stack.data(stack),
  family = "nbinomial",
  control.predictor = list(A = inla.stack.A(stack), compute = TRUE),
  control.fixed = list(expand.factor.strategy = 'inla'),
  control.compute = list(config = TRUE, dic = TRUE, waic = TRUE),
  verbose = TRUE
)

mod24_out <- inla(
  mod24,
  data = inla.stack.data(stack),
  family = "nbinomial",
  control.predictor = list(A = inla.stack.A(stack), compute = TRUE),
  control.fixed = list(expand.factor.strategy = 'inla'),
  control.compute = list(config = TRUE, dic = TRUE, waic = TRUE),
  verbose = TRUE
)

mod25_out <- inla(
  mod25,
  data = inla.stack.data(stack),
  family = "nbinomial",
  control.predictor = list(A = inla.stack.A(stack), compute = TRUE),
  control.fixed = list(expand.factor.strategy = 'inla'),
  control.compute = list(config = TRUE, dic = TRUE, waic = TRUE),
  verbose = TRUE
)

mod26_out <- inla(
  mod26,
  data = inla.stack.data(stack),
  family = "nbinomial",
  control.predictor = list(A = inla.stack.A(stack), compute = TRUE),
  control.fixed = list(expand.factor.strategy = 'inla'),
  control.compute = list(config = TRUE, dic = TRUE, waic = TRUE),
  verbose = TRUE
)

mod27_out <- inla(
  mod27,
  data = inla.stack.data(stack),
  family = "nbinomial",
  control.predictor = list(A = inla.stack.A(stack), compute = TRUE),
  control.fixed = list(expand.factor.strategy = 'inla'),
  control.compute = list(config = TRUE, dic = TRUE, waic = TRUE),
  verbose = TRUE
)

mod28_out <- inla(
  mod28,
  data = inla.stack.data(stack),
  family = "nbinomial",
  control.predictor = list(A = inla.stack.A(stack), compute = TRUE),
  control.fixed = list(expand.factor.strategy = 'inla'),
  control.compute = list(config = TRUE, dic = TRUE, waic = TRUE),
  verbose = TRUE
)

extract_marginals <- function(model_out, model_label, years = 2017:2024) {
  map_dfr(years, \(yr) {
    term <- paste0("year", yr, ":mintro")
    zmarg <- inla.zmarginal(model_out$marginals.fixed[[term]], silent = TRUE)
    data.frame(
      model = model_label,
      year  = yr,
      mean  = exp(zmarg$mean),
      lb    = exp(zmarg$quant0.025),
      ub    = exp(zmarg$quant0.975)
    )
  })
}

cuminc_mod_comparison <- bind_rows(
  extract_marginals(modlong2_out, "No adjustment for cumulative incidence"),
  extract_marginals(mod29_out, "1-year cumulative incidence"),
  extract_marginals(mod11_out, "2-year cumulative incidence"),
  extract_marginals(mod12_out, "3-year cumulative incidence"),
  extract_marginals(mod13_out, "4-year cumulative incidence"),
  extract_marginals(mod14_out, "5-year cumulative incidence"),
  extract_marginals(mod15_out, "6-year cumulative incidence"),
  extract_marginals(mod16_out, "7-year cumulative incidence"),
  extract_marginals(mod30_out, "1-year cumulative incidence - NL"),
  extract_marginals(mod17_out, "2-year cumulative incidence - NL"),
  extract_marginals(mod18_out, "3-year cumulative incidence - NL"),
  extract_marginals(mod19_out, "4-year cumulative incidence - NL"),
  extract_marginals(mod20_out, "5-year cumulative incidence - NL"),
  extract_marginals(mod21_out, "6-year cumulative incidence - NL"),
  extract_marginals(mod22_out, "7-year cumulative incidence - NL")
)

cuminc_mod_comparison_summ <- cuminc_mod_comparison |>
  mutate(mean = round(mean, 2),
         lb = round(lb, 2),
         ub = round(ub, 2),
         summ = paste0(mean, " (", lb, "-", ub, ")")) |>
  dplyr::select(-mean, -lb, -ub) |>
  pivot_wider(names_from = year, values_from = "summ")

n_models <- n_distinct(cuminc_mod_comparison$model)
offsets   <- seq(-0.4, 0.4, length.out = n_models)
names(offsets) <- unique(cuminc_mod_comparison$model)

cuminc_mod_comparison |>
  mutate(year_offset = year + offsets[model],
         model = factor(model, levels = c(
           "No adjustment for cumulative incidence", "1-year cumulative incidence",
           "2-year cumulative incidence", "3-year cumulative incidence",
           "4-year cumulative incidence", "5-year cumulative incidence",
           "6-year cumulative incidence", "7-year cumulative incidence",
           "1-year cumulative incidence - NL", "2-year cumulative incidence - NL",
           "3-year cumulative incidence - NL", "4-year cumulative incidence - NL",
           "5-year cumulative incidence - NL", "6-year cumulative incidence - NL",
           "7-year cumulative incidence - NL"
         ))) |>
  ggplot(aes(x = year_offset, y = mean, color = model)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  geom_errorbar(aes(ymin = lb, ymax = ub), width = 0.1, linewidth = 0.7) +
  geom_point(size = 2.5) +
  scale_x_continuous(breaks = 2017:2024, labels = 2017:2024) +
  labs(x = "Year",
       y = "Ratio of incidences",
       color = "Adjustment") +
  theme_bw() +
  theme(legend.position = "right")

cuminc_mod_comparison2 <- bind_rows(
  extract_marginals(modlong2_out, "No adjustment for cumulative incidence"),
  extract_marginals(mod29_out, "1-year cumulative incidence"),
  extract_marginals(mod11_out, "2-year cumulative incidence"),
  extract_marginals(mod12_out, "3-year cumulative incidence"),
  extract_marginals(mod13_out, "4-year cumulative incidence"),
  extract_marginals(mod14_out, "5-year cumulative incidence"),
  extract_marginals(mod15_out, "6-year cumulative incidence"),
  extract_marginals(mod16_out, "7-year cumulative incidence"),
  extract_marginals(mod31_out, "1-year cumulative incidence - NL iid"),
  extract_marginals(mod23_out, "2-year cumulative incidence - NL iid"),
  extract_marginals(mod24_out, "3-year cumulative incidence - NL iid"),
  extract_marginals(mod25_out, "4-year cumulative incidence - NL iid"),
  extract_marginals(mod26_out, "5-year cumulative incidence - NL iid"),
  extract_marginals(mod27_out, "6-year cumulative incidence - NL iid"),
  extract_marginals(mod28_out, "7-year cumulative incidence - NL iid")
)

n_models <- n_distinct(cuminc_mod_comparison2$model)
offsets   <- seq(-0.4, 0.4, length.out = n_models)
names(offsets) <- unique(cuminc_mod_comparison2$model)

cuminc_mod_comparison2 |>
  mutate(year_offset = year + offsets[model],
         model = factor(model, levels = c(
           "No adjustment for cumulative incidence", "1-year cumulative incidence",
           "2-year cumulative incidence", "3-year cumulative incidence",
           "4-year cumulative incidence", "5-year cumulative incidence",
           "6-year cumulative incidence", "7-year cumulative incidence",
           "1-year cumulative incidence - NL iid", "2-year cumulative incidence - NL iid",
           "3-year cumulative incidence - NL iid", "4-year cumulative incidence - NL iid",
           "5-year cumulative incidence - NL iid", "6-year cumulative incidence - NL iid",
           "7-year cumulative incidence - NL iid"
         ))) |>
  ggplot(aes(x = year_offset, y = mean, color = model)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  geom_errorbar(aes(ymin = lb, ymax = ub), width = 0.1, linewidth = 0.7) +
  geom_point(size = 2.5) +
  scale_x_continuous(breaks = 2017:2024, labels = 2017:2024) +
  labs(x = "Year",
       y = "Ratio of incidences",
       color = "Adjustment") +
  theme_bw() +
  theme(legend.position = "right")

cuminc_mod_comparison_summ2 <- cuminc_mod_comparison2 |>
  mutate(mean = round(mean, 2),
         lb = round(lb, 2),
         ub = round(ub, 2),
         summ = paste0(mean, " (", lb, "-", ub, ")")) |>
  dplyr::select(-mean, -lb, -ub) |>
  pivot_wider(names_from = year, values_from = "summ")

extract_fixed <- function(model_out, term_name, model_label) {
  term <- term_name
  zmarg <- inla.zmarginal(model_out$marginals.fixed[[term]], silent = TRUE)
  data.frame(
      model = model_label,
      effect  = exp(zmarg$mean),
      lb    = exp(zmarg$quant0.025),
      ub    = exp(zmarg$quant0.975)
    )
}

extract_rw1 <- function(model_out, model_label, inc_var, inc_var_binned_name, inc_var_binned) {
  rw1 <- model_out$summary.random[[inc_var_binned_name]]
  midpoints <- tapply(inc_var, cut(inc_var, breaks = length(unique(inc_var_binned))),
                      mean, na.rm = TRUE)
  rw1 |>
    mutate(model    = model_label,
           inc_orig = ID,
           effect   = exp(mean),
           lb       = exp(`0.025quant`),
           ub       = exp(`0.975quant`))
}

inc_fixed <- bind_rows(
  extract_fixed(mod29_out, "inc_1y", "1y linear"),
  extract_fixed(mod11_out, "inc_2y", "2y linear"),
  extract_fixed(mod12_out, "inc_3y", "3y linear"),
  extract_fixed(mod13_out, "inc_4y", "4y linear"),
  extract_fixed(mod14_out, "inc_5y", "5y linear"),
  extract_fixed(mod15_out, "inc_6y", "6y linear"),
  extract_fixed(mod16_out, "inc_7y", "7y linear")
)

inc_rw1 <- bind_rows(
  extract_rw1(mod30_out, "1y nonlinear", data_mod11$inc_1y, "inc_1y_binned", data_mod11$inc_1y_binned),
  extract_rw1(mod17_out, "2y nonlinear", data_mod11$inc_2y, "inc_2y_binned", data_mod11$inc_2y_binned),
  extract_rw1(mod18_out, "3y nonlinear", data_mod11$inc_3y, "inc_3y_binned", data_mod11$inc_3y_binned),
  extract_rw1(mod19_out, "4y nonlinear", data_mod11$inc_4y, "inc_4y_binned", data_mod11$inc_4y_binned),
  extract_rw1(mod20_out, "5y nonlinear", data_mod11$inc_5y, "inc_5y_binned", data_mod11$inc_5y_binned),
  extract_rw1(mod21_out, "6y nonlinear", data_mod11$inc_6y, "inc_6y_binned", data_mod11$inc_6y_binned),
  extract_rw1(mod22_out, "7y nonlinear", data_mod11$inc_7y, "inc_7y_binned", data_mod11$inc_7y_binned)
)

inc_iid <- bind_rows(
  extract_rw1(mod31_out, "1y nonlinear iid", data_mod11$inc_1y, "inc_1y_binned", data_mod11$inc_1y_binned),
  extract_rw1(mod23_out, "2y nonlinear iid", data_mod11$inc_2y, "inc_2y_binned", data_mod11$inc_2y_binned),
  extract_rw1(mod24_out, "3y nonlinear iid", data_mod11$inc_3y, "inc_3y_binned", data_mod11$inc_3y_binned),
  extract_rw1(mod25_out, "4y nonlinear iid", data_mod11$inc_4y, "inc_4y_binned", data_mod11$inc_4y_binned),
  extract_rw1(mod26_out, "5y nonlinear iid", data_mod11$inc_5y, "inc_5y_binned", data_mod11$inc_5y_binned),
  extract_rw1(mod27_out, "6y nonlinear iid", data_mod11$inc_6y, "inc_6y_binned", data_mod11$inc_6y_binned),
  extract_rw1(mod28_out, "7y nonlinear iid", data_mod11$inc_7y, "inc_7y_binned", data_mod11$inc_7y_binned)
)

inc_fixed_summ <- inc_fixed |>
  mutate(effect = round(effect, 2),
         lb = round(lb, 2),
         ub = round(ub, 2),
         summ = paste0(effect, " (", lb, "-", ub, ")"))

inc_rw1_summ <- inc_rw1 |>
  mutate(effect = round(effect, 2),
         lb = round(lb, 2),
         ub = round(ub, 2),
         summ = paste0(effect, " (", lb, "-", ub, ")"),
         inc = round(inc_orig, 2)
         )

inc_iid_summ <- inc_iid |>
  mutate(effect = round(effect, 2),
         lb = round(lb, 2),
         ub = round(ub, 2),
         summ = paste0(effect, " (", lb, "-", ub, ")"),
         inc = round(inc_orig, 2)
         )

test <- data_mod11 |>
  mutate(test = cut(inc_1y,
                    breaks = unique(quantile(data_mod11$inc_1y, probs = seq(0, 1, 0.2), na.rm = TRUE)),
                    include.lowest = TRUE,
                    labels = FALSE)) |>
  group_by(test) |>
  mutate(median = median(inc_1y, na.rm = TRUE)) |>
  ungroup()

inc_fixed |>
  ggplot(aes(x = 1, color = model)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  geom_errorbar(aes(ymin = lb, ymax = ub), width = 0.1, linewidth = 0.7, position = position_dodge(width = 0.3)) +
  geom_point(aes(y = effect), size = 2.5, position = position_dodge(width = 0.3)) +
  labs(x = "",
       y = "Ratio of incidences",
       color = "Adjustment",
       subtitle = "Effect of cumulative incidence (linear)") +
  theme_bw() +
  theme(legend.position = "right",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

inc_rw1 |>
  ggplot(aes(x = inc_orig, color = model)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  geom_errorbar(aes(ymin = lb, ymax = ub), width = 0.1, linewidth = 0.7) +
  geom_point(aes(y = effect), size = 2.5) +
  labs(x = "Cumulative incidence",
       y = "Ratio of incidences",
       color = "Adjustment",
       subtitle = "Effect of cumulative incidence (nonlinear - rw1)") +
  theme_bw() +
  theme(legend.position = "right")

inc_rw1 |>
  ggplot(aes(x = inc_orig, color = model)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  geom_errorbar(aes(ymin = lb, ymax = ub), width = 0.2, linewidth = 0.7) +
  geom_point(aes(y = effect), size = 2.5) +
  labs(x = "Cumulative incidence",
       y = "Ratio of incidences",
       color = "Adjustment",
       subtitle = "Effect of cumulative incidence (nonlinear - rw1)") +
  # scale_x_continuous(breaks = seq(0, 11, by = 1)) +
  theme_bw() +
  facet_wrap(~model, ncol = 1) +
  theme(legend.position = "none")

inc_rw1 |>
  ggplot(aes(x = inc_orig, color = model)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  geom_errorbar(aes(ymin = lb, ymax = ub), width = 0.2, linewidth = 0.7) +
  geom_point(aes(y = effect), size = 2.5) +
  labs(x = "Cumulative incidence",
       y = "Ratio of incidences",
       color = "Adjustment",
       subtitle = "Effect of cumulative incidence (nonlinear - rw1)") +
  # scale_x_continuous(breaks = seq(0, 11, by = 1)) +
  theme_bw() +
  facet_wrap(~model, ncol = 1, scales = "free_y") +
  theme(legend.position = "none")

inc_iid |>
  ggplot(aes(x = inc_orig, color = model)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  geom_errorbar(aes(ymin = lb, ymax = ub), width = 0.2, linewidth = 0.7) +
  geom_point(aes(y = effect), size = 2.5) +
  labs(x = "Cumulative incidence",
       y = "Ratio of incidences",
       color = "Adjustment",
       subtitle = "Effect of cumulative incidence (nonlinear - iid)") +
  # scale_x_continuous(breaks = seq(0, 11, by = 1)) +
  theme_bw() +
  facet_wrap(~model, ncol = 1) +
  theme(legend.position = "none")

inc_iid |>
  ggplot(aes(x = inc_orig, color = model)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  geom_errorbar(aes(ymin = lb, ymax = ub), width = 0.2, linewidth = 0.7) +
  geom_point(aes(y = effect), size = 2.5) +
  labs(x = "Cumulative incidence",
       y = "Ratio of incidences",
       color = "Adjustment",
       subtitle = "Effect of cumulative incidence (nonlinear - iid)") +
  # scale_x_continuous(breaks = seq(0, 11, by = 1)) +
  theme_bw() +
  facet_wrap(~model, ncol = 1, scales = "free_y") +
  theme(legend.position = "none")

get_waic <- function(model_out, model_label) {
  waic <- model_out$waic$waic
  dic <- model_out$dic$dic
  data.frame(
      model = model_label,
      waic = waic,
      dic = dic
    )
}

compare_waic <- bind_rows(
  get_waic(modlong2_out, "No adjustment for cumulative incidence"),
  get_waic(mod29_out, "1-year cumulative incidence"),
  get_waic(mod11_out, "2-year cumulative incidence"),
  get_waic(mod12_out, "3-year cumulative incidence"),
  get_waic(mod13_out, "4-year cumulative incidence"),
  get_waic(mod14_out, "5-year cumulative incidence"),
  get_waic(mod15_out, "6-year cumulative incidence"),
  get_waic(mod16_out, "7-year cumulative incidence"),
  get_waic(mod30_out, "1-year cumulative incidence - rw1"),
  get_waic(mod17_out, "2-year cumulative incidence - rw1"),
  get_waic(mod18_out, "3-year cumulative incidence - rw1"),
  get_waic(mod19_out, "4-year cumulative incidence - rw1"),
  get_waic(mod20_out, "5-year cumulative incidence - rw1"),
  get_waic(mod21_out, "6-year cumulative incidence - rw1"),
  get_waic(mod22_out, "7-year cumulative incidence - rw1"),
  get_waic(mod31_out, "1-year cumulative incidence - iid"),
  get_waic(mod23_out, "2-year cumulative incidence - iid"),
  get_waic(mod24_out, "3-year cumulative incidence - iid"),
  get_waic(mod25_out, "4-year cumulative incidence - iid"),
  get_waic(mod26_out, "5-year cumulative incidence - iid"),
  get_waic(mod27_out, "6-year cumulative incidence - iid"),
  get_waic(mod28_out, "7-year cumulative incidence - iid")
)

mod32 <- obs ~ -1 + Intercept + f(nyear, model = "ar1") + year:mintro + f(spatial.field, model = spde) + inc_3y + inc_5y + offset(logpop)

mod33 <- obs ~ -1 + Intercept + f(nyear, model = "ar1") + year:mintro + f(spatial.field, model = spde) +
  f(inc_3y_binned, model = "rw1", constr = FALSE, extraconstr = extraconstr_mod18) +
  f(inc_5y_binned, model = "rw1", constr = FALSE, extraconstr = extraconstr_mod20) + offset(logpop)

mod34 <- obs ~ -1 + Intercept + f(nyear, model = "ar1") + year:mintro + f(spatial.field, model = spde) +
  f(inc_3y_binned, model = "iid", constr = FALSE, extraconstr = extraconstr_mod18) +
  f(inc_5y_binned, model = "iid", constr = FALSE, extraconstr = extraconstr_mod20) + offset(logpop)

mod32_out <- inla(
  mod32,
  data = inla.stack.data(stack),
  family = "nbinomial",
  control.predictor = list(A = inla.stack.A(stack), compute = TRUE),
  control.fixed = list(expand.factor.strategy = 'inla'),
  control.compute = list(config = TRUE, dic = TRUE, waic = TRUE)
)

mod33_out <- inla(
  mod33,
  data = inla.stack.data(stack),
  family = "nbinomial",
  control.predictor = list(A = inla.stack.A(stack), compute = TRUE),
  control.fixed = list(expand.factor.strategy = 'inla'),
  control.compute = list(config = TRUE, dic = TRUE, waic = TRUE)
)

mod34_out <- inla(
  mod34,
  data = inla.stack.data(stack),
  family = "nbinomial",
  control.predictor = list(A = inla.stack.A(stack), compute = TRUE),
  control.fixed = list(expand.factor.strategy = 'inla'),
  control.compute = list(config = TRUE, dic = TRUE, waic = TRUE)
)

cuminc_two_mod_comparison <- bind_rows(
  extract_marginals(modlong2_out, "No adjustment"),
  extract_marginals(mod32_out, "3-year + 5-year linear"),
  extract_marginals(mod33_out, "3-year + 5-year RW1"),
  extract_marginals(mod34_out, "3-year + 5-year IID")
)

n_models <- n_distinct(cuminc_two_mod_comparison$model)
offsets   <- seq(-0.4, 0.4, length.out = n_models)
names(offsets) <- unique(cuminc_two_mod_comparison$model)

cuminc_two_mod_comparison |>
  mutate(year_offset = year + offsets[model],
         model = factor(model, levels = c(
           "No adjustment", "3-year + 5-year linear",
           "3-year + 5-year RW1", "3-year + 5-year IID"
         ))) |>
  ggplot(aes(x = year_offset, y = mean, color = model)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  geom_errorbar(aes(ymin = lb, ymax = ub), width = 0.1, linewidth = 0.7) +
  geom_point(size = 2.5) +
  scale_x_continuous(breaks = 2017:2024, labels = 2017:2024) +
  labs(x = "Year",
       y = "Ratio of incidences",
       color = "Adjustment") +
  theme_bw() +
  theme(legend.position = "right")

inc_fixed_two <- bind_rows(
  extract_fixed(mod32_out, "inc_3y", "3y linear"),
  extract_fixed(mod32_out, "inc_5y", "5y linear")
)

inc_fixed_two |>
  ggplot(aes(x = 1, color = model)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  geom_errorbar(aes(ymin = lb, ymax = ub), width = 0.1, linewidth = 0.7, position = position_dodge(width = 0.3)) +
  geom_point(aes(y = effect), size = 2.5, position = position_dodge(width = 0.3)) +
  labs(x = "",
       y = "Ratio of incidences",
       color = "Adjustment",
       subtitle = "Effect of cumulative incidence (linear)") +
  theme_bw() +
  theme(legend.position = "right",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

inc_rw1_two <- bind_rows(
  extract_rw1(mod33_out, "3y nonlinear", data_mod11$inc_3y, "inc_3y_binned", data_mod11$inc_3y_binned),
  extract_rw1(mod33_out, "5y nonlinear", data_mod11$inc_5y, "inc_5y_binned", data_mod11$inc_5y_binned)
)

inc_iid_two <- bind_rows(
  extract_rw1(mod34_out, "3y nonlinear", data_mod11$inc_3y, "inc_3y_binned", data_mod11$inc_3y_binned),
  extract_rw1(mod34_out, "5y nonlinear", data_mod11$inc_5y, "inc_5y_binned", data_mod11$inc_5y_binned)
)

inc_rw1_two |>
  ggplot(aes(x = inc_orig, color = model)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  geom_errorbar(aes(ymin = lb, ymax = ub), width = 0.2, linewidth = 0.7) +
  geom_point(aes(y = effect), size = 2.5) +
  labs(x = "Cumulative incidence",
       y = "Ratio of incidences",
       color = "Adjustment",
       subtitle = "Effect of cumulative incidence (nonlinear - rw1)") +
  # scale_x_continuous(breaks = seq(0, 11, by = 1)) +
  theme_bw() +
  facet_wrap(~model, ncol = 1, scales = "free_y") +
  theme(legend.position = "none")

inc_iid_two |>
  ggplot(aes(x = inc_orig, color = model)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  geom_errorbar(aes(ymin = lb, ymax = ub), width = 0.2, linewidth = 0.7) +
  geom_point(aes(y = effect), size = 2.5) +
  labs(x = "Cumulative incidence",
       y = "Ratio of incidences",
       color = "Adjustment",
       subtitle = "Effect of cumulative incidence (nonlinear - iid)") +
  # scale_x_continuous(breaks = seq(0, 11, by = 1)) +
  theme_bw() +
  facet_wrap(~model, ncol = 1, scales = "free_y") +
  theme(legend.position = "none")

compare_waic2 <- bind_rows(
  get_waic(modlong2_out, "No adjustment"),
  get_waic(mod32_out, "3 + 5 year linear"),
  get_waic(mod33_out, "3 + 5 year RW1"),
  get_waic(mod34_out, "3 + 5 year IID")
)

# ---- Figure 3 ----------------------------------------------------------

# -- Figure 3A: Time series between 2010 and 2017
(fig3a <- dengues |>
  filter(date < as.Date("2017-01-02")) |>
  ggplot(aes(x = date)) +
  geom_line(aes(y = cases), size = 1.5) +
  labs(x = "", y = "Monthly\ndengue cases") +
  scale_x_date(date_labels = "%Y", breaks = seq.Date(as.Date("2010-01-01"), as.Date("2017-01-01"), by = "12 months"),
               expand = c(0.01, 0.01)) +
  theme_bw() +
  theme(axis.text.x = element_text(),
        axis.title.x = element_blank(),
        text = element_text(size = 16),
        plot.title = element_text(hjust = 0.5)))

# ggsave("outputs/fig1b-may1.pdf", plot = fig1a)

# -- Figure 3B: Future introgression model
data_mod3 <- grid |>
  right_join(dengue_grid, by = "id") |>
  left_join(popdens_grid, by = "id") |>
  left_join(intro_1920_grid, by = "id") |>
  rename(pop2 = sum) |>
  mutate(pop = if_else(pop2 < 50, 0, pop2),
         inc = cases/pop*1000) |>
  filter(st_intersects(geometry, release_areas_j, sparse = FALSE)[,1]) |>
  mutate(geometry = st_centroid(geometry),
         x = st_coordinates(geometry)[,1],
         y = st_coordinates(geometry)[,2]) |>
  st_drop_geometry() |>
  filter(pop != 0 & year <= 2016) |>
  mutate(logpop = log(pop),
         year = factor(year)) |>
  arrange(year) |>
  group_by(year) |>
  mutate(nyear = cur_group_id()) |>
  ungroup() |>
  filter(!is.na(mintro))

sum(data_mod3$cases[data_mod3$year == 2010])
sum(data_mod3$cases[data_mod3$year == 2011])
sum(data_mod3$cases[data_mod3$year == 2012])
sum(data_mod3$cases[data_mod3$year == 2013])
sum(data_mod3$cases[data_mod3$year == 2014])
sum(data_mod3$cases[data_mod3$year == 2015])
sum(data_mod3$cases[data_mod3$year == 2016])

spatial.range <- 22000

mesh_mod3 <- inla.mesh.2d(loc = as.matrix(data_mod3[, c("x", "y")]),
                        max.edge = c(spatial.range / 5, spatial.range),
                        cutoff = 200,
                        offset = c(spatial.range / 5, diff(range(data_mod3[, "x"])) / 10))

spde_mod3 <- inla.spde2.matern(mesh = mesh_mod3)

s.index_mod3 <- inla.spde.make.index(name = "spatial.field",
                                  n.spde = spde_mod3$n.spde)

A.est_mod3 <- inla.spde.make.A(mesh = mesh_mod3,
                            loc = as.matrix(data_mod3[, c("x", "y")]))

stack_mod3 <- inla.stack(
  data = list(obs = data_mod3$cases),
  A = list(A.est_mod3, 1),
  effects = list(spatial_field = data.frame(s.index_mod3, Intercept = 1),
                 data_mod3[, c(names(data_mod3) != "cases")]),
  tag = "stdata"
)

f_mod3 <- obs ~ -1 + Intercept +
  year:mintro +
  f(nyear, model = "ar1") +
  offset(logpop)

mod3_out <- inla(
  f_mod3,
  data = inla.stack.data(stack_mod3),
  family = "nbinomial",
  control.predictor = list(A = inla.stack.A(stack_mod3)),
  control.fixed = list(expand.factor.strategy = 'inla'),
  control.compute = list(config = TRUE, dic = TRUE, waic = TRUE),
  verbose = TRUE
  )

exp(as.data.frame(inla.zmarginal(mod3_out$marginals.fixed$`year2010:mintro`)))
exp(as.data.frame(inla.zmarginal(mod3_out$marginals.fixed$`year2011:mintro`)))
exp(as.data.frame(inla.zmarginal(mod3_out$marginals.fixed$`year2012:mintro`)))
exp(as.data.frame(inla.zmarginal(mod3_out$marginals.fixed$`year2013:mintro`)))
exp(as.data.frame(inla.zmarginal(mod3_out$marginals.fixed$`year2014:mintro`)))
exp(as.data.frame(inla.zmarginal(mod3_out$marginals.fixed$`year2015:mintro`)))
exp(as.data.frame(inla.zmarginal(mod3_out$marginals.fixed$`year2016:mintro`)))

mintro_int_margs_mod3 <- list(
  `2010` = mod3_out$marginals.fixed$`year2010:mintro`,
  `2011` = mod3_out$marginals.fixed$`year2011:mintro`,
  `2012` = mod3_out$marginals.fixed$`year2012:mintro`,
  `2013` = mod3_out$marginals.fixed$`year2013:mintro`,
  `2014` = mod3_out$marginals.fixed$`year2014:mintro`,
  `2015` = mod3_out$marginals.fixed$`year2015:mintro`,
  `2016` = mod3_out$marginals.fixed$`year2016:mintro`
)

set.seed(123)
mod3_df <- bind_rows(lapply(names(mintro_int_margs_mod3), function(yr) {
  beta_year_draws <- inla.rmarginal(1000, mintro_int_margs_mod3[[yr]])
  effect_samples <- sapply(beta_year_draws, function(b) exp(b * mintro_values))

  data.frame(
    year   = as.integer(yr),
    mintro = mintro_values,
    effect = apply(effect_samples, 1, mean),
    lower  = apply(effect_samples, 1, quantile, probs = 0.025),
    upper  = apply(effect_samples, 1, quantile, probs = 0.975)
  )
}))

(fig3b <- mod3_df |>
  ggplot(aes(x = mintro, y = effect)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_hline(aes(yintercept = 1), lty = 2, color = "red") +
  labs(x = "Future % wMel introgression in 2019-2020", y = "Ratio of incidences") +
  theme_bw() +
  facet_wrap(~year, ncol = 7) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1")) +
  scale_y_continuous(breaks = c(0, 3, 6, 9)) +
  coord_trans(y = "log1p") +
  theme(text = element_text(size = 16)))

# -- Figure 3C: Time-varying spatial field
data_mod4 <- data_mod3

spatial.range <- 22000

mesh_mod4 <- inla.mesh.2d(loc = as.matrix(unique(data_mod4[, c("x", "y")])),
                          max.edge = c(spatial.range / 5, spatial.range),
                          cutoff = 200,
                          offset = c(spatial.range / 5, diff(range(data_mod4[, "x"])) / 10))

spde_mod4 <- inla.spde2.matern(mesh = mesh_mod4)

s.index_mod4 <- inla.spde.make.index(name = "spatial.field",
                                     n.spde = spde_mod4$n.spde,
                                     n.group = length(unique(data_mod4$year)))

A.est_mod4 <- inla.spde.make.A(mesh = mesh_mod4,
                               loc = as.matrix(data_mod4[, c("x", "y")]),
                               group = data_mod4$nyear,
                               n.group = length(unique(data_mod4$nyear)))

stack_mod4 <- inla.stack(
  data = list(obs = data_mod4$cases),
  A = list(A.est_mod4, 1),
  effects = list(spatial_field = data.frame(s.index_mod4, Intercept = 1),
                 data_mod4[, c(names(data_mod4) != "cases")]),
  tag = "stdata"
)

n_spde <- spde_mod4$n.spde
n_years <- length(unique(data_mod4$nyear))

A_constr <- matrix(data = 0, nrow = n_years, ncol = n_spde)

for (i in 1:n_years) {
  A_constr[i, ] <- 1 / n_spde
}

e_constr <- rep(0, n_years)

f_mod4 <- obs ~ -1 + Intercept +
  f(nyear, model = "ar1") +
  f(spatial.field, model = spde_mod4, group = spatial.field.group,
    control.group = list(model = "ar1"),
    extraconstr = list(A = A_constr, e = e_constr)) +
  offset(logpop)

mod4_out <- inla(
  f_mod4,
  data = inla.stack.data(stack_mod4),
  family = "nbinomial",
  control.predictor = list(A = inla.stack.A(stack_mod4)),
  control.fixed = list(expand.factor.strategy = 'inla'),
  control.compute = list(config = TRUE, dic = TRUE, waic = TRUE),
  verbose = TRUE
  )

sf_mod4 <- mod4_out$summary.random$spatial.field |>
  mutate(nyear = rep(unique(data_mod4$nyear), each = length(mesh_mod4$loc[,1])),
         year = rep(unique(data_mod4$year), each = length(mesh_mod4$loc[,1])))

unique_years <- unique(data_mod4$year)

coords_mod4 <- spde_mod4$mesh$loc

xrange_mod4 <- seq(min(coords_mod4[, 1]), max(coords_mod4[, 1]), length.out = 100)
yrange_mod4 <- seq(min(coords_mod4[, 2]), max(coords_mod4[, 2]), length.out = 100)
grid_mod4 <- expand.grid(x = xrange_mod4, y = yrange_mod4)
sp::coordinates(grid_mod4) <- ~x + y
proj_grid_mod4 <- inla.mesh.projector(spde_mod4$mesh, xlim = range(xrange_mod4), ylim = range(yrange_mod4), dims = c(100, 100))

sf_year_mod4 <- lapply(unique_years, function(yr) {
  spatial_mean <- inla.mesh.project(proj_grid_mod4, sf_mod4$mean[sf_mod4$year == yr])
  spatial_lb <- inla.mesh.project(proj_grid_mod4, sf_mod4$`0.025quant`[sf_mod4$year == yr])
  spatial_ub <- inla.mesh.project(proj_grid_mod4, sf_mod4$`0.975quant`[sf_mod4$year == yr])
  data.frame(
    x = grid_mod4$x,
    y = grid_mod4$y,
    year = yr,
    mean = as.vector(spatial_mean),
    lb = as.vector(spatial_lb),
    ub = as.vector(spatial_ub)
  )
})

sf_year_mod4 <- do.call(rbind, sf_year_mod4)

grid_sf_mod4 <- sf_year_mod4|>
  st_as_sf(coords = c("x", "y"), remove = FALSE, crs = 31983) |>
  mutate(sig = case_when(lb > 0 & ub > 0 ~ 1,
                         lb < 0 & ub < 0 ~ 1,
                         .default = 0)) |>
  st_intersection(release_areas_j)

(fig3c <- grid_sf_mod4 |>
  ggplot() +
  geom_tile(aes(x = x, y = y, fill = exp(mean))) +
  scale_fill_gradient2(low = "#25489E", mid = "white", high = "#CE263D", midpoint = 1) +
  geom_sf(data = release_areas_j, fill = NA, color = "black") +
  theme_void() +
  labs(fill = "Ratio of incidences", y = "Gaussian field") +
  facet_wrap(~year, ncol = 7) +
  theme(text = element_text(size = 16),
        axis.title.y = element_text(angle = 90),
        axis.text.y = element_text(color = "transparent", size = 8),
        axis.ticks.y = element_line(color = "transparent"),
        strip.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.position = "bottom",
        legend.key.height = unit(12, "pt"),
        axis.title = element_blank(),
        panel.spacing.x = unit(-2, "mm")))

# -- Complete figure 3
(fig3 <- fig3a / fig3c / fig3b +
   plot_annotation(tag_levels = "A") +
   plot_layout(heights = c(2, 3, 3)))

# ggsave("outputs/2026-04-13_fig3.pdf", plot = fig3, width = 13, height = 7)
# ggsave("outputs/fig3-oct08.png", plot = fig3, width = 13, height = 7)

# ---- Figure 4 ----------------------------------------------------------

# -- Get low transmission spatial field (2015-2016)
data_mod5 <- grid |>
  right_join(dengue_grid, by = "id") |>
  left_join(popdens_grid, by = "id") |>
  rename(pop2 = sum) |>
  mutate(pop = if_else(pop2 < 50, 0, pop2),
         inc = cases/pop*1000) |>
  filter(st_intersects(geometry, release_areas_j, sparse = FALSE)[,1]) |>
  mutate(geometry = st_centroid(geometry),
         x = st_coordinates(geometry)[,1],
         y = st_coordinates(geometry)[,2]) |>
  st_drop_geometry() |>
  filter(pop != 0 & year %in% c(2015, 2016)) |>
  mutate(logpop = log(pop),
         year = factor(year)) |>
  arrange(year) |>
  group_by(year) |>
  mutate(nyear = cur_group_id()) |>
  ungroup()

spatial.range <- 22000

mesh_mod5 <- inla.mesh.2d(loc = as.matrix(data_mod5[, c("x", "y")]),
                          max.edge = c(spatial.range / 5, spatial.range),
                          cutoff = 200,
                          offset = c(spatial.range / 5, diff(range(data_mod5[, "x"])) / 10))

spde_mod5 <- inla.spde2.matern(mesh = mesh_mod5)

s.index_mod5 <- inla.spde.make.index(name = "spatial.field",
                                     n.spde = spde_mod5$n.spde)

A.est_mod5 <- inla.spde.make.A(mesh = mesh_mod5,
                               loc = as.matrix(data_mod5[, c("x", "y")]))

stack_mod5 <- inla.stack(data = list(obs = data_mod5$cases),
                         A = list(A.est_mod5, 1),
                         effects = list(spatial_field = data.frame(s.index_mod5,
                                                                   Intercept = 1),
                                        data_mod5[, c(names(data_mod5) != "cases")]),
                         tag = "stdata")

f_mod5 <- obs ~ -1 + Intercept +
  f(nyear, model = "ar1") +
  f(spatial.field, model = spde_mod5) +
  offset(logpop)

mod5_out <- inla(f_mod5,
                 data = inla.stack.data(stack_mod5),
                 family = "nbinomial",
                 control.predictor = list(A = inla.stack.A(stack_mod5)),
                 control.fixed = list(expand.factor.strategy = 'inla'),
                 control.compute = list(config = TRUE, dic = TRUE, waic = TRUE),
                 verbose = TRUE)

sf_mod5 <- mod5_out$summary.random$spatial.field

spest_mod5 <- as.numeric(A.est_mod5 %*% sf_mod5$mean)

spatfield_low <- data_mod5 |>
  mutate(spde = spest_mod5) |>
  distinct(id, spde)

coords_mod5 <- spde_mod5$mesh$loc

xrange_mod5 <- seq(min(coords_mod5[, 1]), max(coords_mod5[, 1]), length.out = 100)
yrange_mod5 <- seq(min(coords_mod5[, 2]), max(coords_mod5[, 2]), length.out = 100)

grid_mod5 <- expand.grid(x = xrange_mod5, y = yrange_mod5)
sp::coordinates(grid_mod5) <- ~x + y
proj_grid_mod5 <- inla.mesh.projector(spde_mod5$mesh, xlim = range(xrange_mod5), ylim = range(yrange_mod5), dims = c(100, 100))

spatial_mean_mod5 <- inla.mesh.project(proj_grid_mod5, mod5_out$summary.random$spatial.field$mean)
spatial_lb_mod5<- inla.mesh.project(proj_grid_mod5, mod5_out$summary.random$spatial.field$`0.025quant`)
spatial_ub_mod5 <- inla.mesh.project(proj_grid_mod5, mod5_out$summary.random$spatial.field$`0.975quant`)

grid_mod5$mean <- as.vector(spatial_mean_mod5)
grid_mod5$lb <- as.vector(spatial_lb_mod5)
grid_mod5$ub <- as.vector(spatial_ub_mod5)

grid_sf_mod5 <- as.data.frame(grid_mod5) |>
  st_as_sf(coords = c("x", "y"), remove = FALSE, crs = 31983) |>
  mutate(sig = case_when(lb > 0 & ub > 0 ~ 1,
                         lb < 0 & ub < 0 ~ 1,
                         .default = 0)) |>
  st_intersection(release_areas_j)

(fig4_lowinset <- grid_sf_mod5 |>
  ggplot() +
  geom_tile(aes(x = x, y = y, fill = exp(mean)), show.legend = FALSE) +
  scale_fill_gradient2(low = "#25489E", mid = "white", high = "#CE263D",
                       midpoint = 1, limits = c(0.27, 6.6)) +
  geom_sf(data = release_areas_j, fill = NA, color = "black") +
  theme_void() +
  labs(fill = "Ratio of\nincidences") +
  theme(text = element_text(size = 15),
        plot.margin = unit(c(-5, -5, -5, -5), "mm")))

# -- Get high transmission spatial field (2011, 2012, 2013)
data_mod6 <- grid |>
  right_join(dengue_grid, by = "id") |>
  left_join(popdens_grid, by = "id") |>
  rename(pop2 = sum) |>
  mutate(pop = if_else(pop2 < 50, 0, pop2),
         inc = cases/pop*1000) |>
  filter(st_intersects(geometry, release_areas_j, sparse = FALSE)[,1]) |>
  mutate(geometry = st_centroid(geometry),
         x = st_coordinates(geometry)[,1],
         y = st_coordinates(geometry)[,2]) |>
  st_drop_geometry() |>
  filter(pop != 0 & year %in% c(2011:2013)) |>
  mutate(logpop = log(pop),
         year = factor(year)) |>
  arrange(year) |>
  group_by(year) |>
  mutate(nyear = cur_group_id()) |>
  ungroup()

spatial.range <- 22000

mesh_mod6 <- inla.mesh.2d(loc = as.matrix(data_mod6[, c("x", "y")]),
                          max.edge = c(spatial.range / 5, spatial.range),
                          cutoff = 200,
                          offset = c(spatial.range / 5, diff(range(data_mod6[, "x"])) / 10))

spde_mod6 <- inla.spde2.matern(mesh = mesh_mod6)

s.index_mod6 <- inla.spde.make.index(name = "spatial.field",
                                     n.spde = spde_mod6$n.spde)

A.est_mod6 <- inla.spde.make.A(mesh = mesh_mod6,
                               loc = as.matrix(data_mod6[, c("x", "y")]))

stack_mod6 <- inla.stack(data = list(obs = data_mod6$cases),
                         A = list(A.est_mod6, 1),
                         effects = list(spatial_field = data.frame(s.index_mod6,
                                                                   Intercept = 1),
                                        data_mod6[, c(names(data_mod6) != "cases")]),
                         tag = "stdata")

f_mod6 <- obs ~ -1 + Intercept +
  f(nyear, model = "ar1") +
  f(spatial.field, model = spde_mod6) +
  offset(logpop)

mod6_out <- inla(f_mod6,
                 data = inla.stack.data(stack_mod6),
                 family = "nbinomial",
                 control.predictor = list(A = inla.stack.A(stack_mod6)),
                 control.fixed = list(expand.factor.strategy = 'inla'),
                 control.compute = list(config = TRUE, dic = TRUE, waic = TRUE),
                 verbose = TRUE)

sf_mod6 <- mod6_out$summary.random$spatial.field

spest_mod6 <- as.numeric(A.est_mod6 %*% sf_mod6$mean)

spatfield_high <- data_mod6 |>
  mutate(spde = spest_mod6) |>
  distinct(id, spde)

coords_mod6 <- spde_mod6$mesh$loc

xrange_mod6 <- seq(min(coords_mod6[, 1]), max(coords_mod6[, 1]), length.out = 100)
yrange_mod6 <- seq(min(coords_mod6[, 2]), max(coords_mod6[, 2]), length.out = 100)

grid_mod6 <- expand.grid(x = xrange_mod6, y = yrange_mod6)
sp::coordinates(grid_mod6) <- ~x + y
proj_grid_mod6 <- inla.mesh.projector(spde_mod6$mesh, xlim = range(xrange_mod6), ylim = range(yrange_mod6), dims = c(100, 100))

spatial_mean_mod6 <- inla.mesh.project(proj_grid_mod6, mod6_out$summary.random$spatial.field$mean)
spatial_lb_mod6<- inla.mesh.project(proj_grid_mod6, mod6_out$summary.random$spatial.field$`0.025quant`)
spatial_ub_mod6 <- inla.mesh.project(proj_grid_mod6, mod6_out$summary.random$spatial.field$`0.975quant`)

grid_mod6$mean <- as.vector(spatial_mean_mod6)
grid_mod6$lb <- as.vector(spatial_lb_mod6)
grid_mod6$ub <- as.vector(spatial_ub_mod6)

grid_sf_mod6 <- as.data.frame(grid_mod6) |>
  st_as_sf(coords = c("x", "y"), remove = FALSE, crs = 31983) |>
  mutate(sig = case_when(lb > 0 & ub > 0 ~ 1,
                         lb < 0 & ub < 0 ~ 1,
                         .default = 0)) |>
  st_intersection(release_areas_j)

(fig4_highinset <- grid_sf_mod6 |>
  ggplot() +
  geom_tile(aes(x = x, y = y, fill = exp(mean))) +
  scale_fill_gradient2(low = "#25489E", mid = "white", high = "#CE263D",
                       midpoint = 1, limits = c(0.27, 6.6)) +
  geom_sf(data = release_areas_j, fill = NA, color = "black") +
  theme_void() +
  labs(fill = "Ratio of\nincidences") +
  theme(text = element_text(size = 15),
        plot.margin = unit(c(-5, -5, -5, -5), "mm"),
        legend.position = "inside",
        legend.position.inside = c(0,0),
        legend.direction = "horizontal"))

# -- Figure 4A: 2023 cases with low transmission field (2015-2016)
(fig4a1 <- dengues |>
  ggplot(aes(x = date)) +
  geom_rect(aes(xmin = as.Date("2015-01-01"), xmax = as.Date("2016-12-31"),
                ymin = -Inf, ymax = Inf), fill = "lightgray") +
  geom_rect(aes(xmin = as.Date("2023-01-01"), xmax = as.Date("2023-11-30"),
                ymin = -Inf, ymax = Inf), fill = "lightgray") +
  geom_line(aes(y = cases), size = 1.5) +
  annotate("text", x = mean(c(as.Date("2015-01-01"), as.Date("2016-12-31"))),
           y = 4000, label = "Spatial field", color = "gray10") +
  annotate("text", x = as.Date("2023-12-01"), y = 4000,
           label = "Cases\nin model", color = "gray10", lineheight = 0.85, hjust = 1) +
  labs(x = "", y = "Cases") +
  scale_x_date(date_labels = "%Y",
               breaks = seq.Date(as.Date("2010-01-01"), as.Date("2024-01-01"), by = "24 months"),
               expand = c(0.01, 0.01)) +
  scale_y_continuous(n.breaks = 3) +
  theme_bw() +
  theme(axis.text.x = element_text(),
        axis.title.x = element_blank(),
        text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5),
        plot.margin = unit(c(3, 4, 0, 4), "pt")) +
  coord_cartesian(ylim = c(0, 5000)))

data_mod7 <- data_mod_long2 |>
  mutate(year = as.numeric(as.character(year))) |>
  filter(year == 2023) |>
  mutate(logpop = log(pop)) |>
  ungroup() |>
  left_join(spatfield_low, by = "id") |>
  dplyr::select(-year)

spatial.range <- 22000

mesh_mod7 <- inla.mesh.2d(loc = as.matrix(data_mod7[, c("x", "y")]),
                          max.edge = c(spatial.range / 5, spatial.range),
                          cutoff = 200,
                          offset = c(spatial.range / 5, diff(range(data_mod7[, "x"])) / 10))

spde_mod7 <- inla.spde2.matern(mesh = mesh_mod7)

s.index_mod7 <- inla.spde.make.index(name = "spatial.field",
                                     n.spde = spde_mod7$n.spde)

A.est_mod7 <- inla.spde.make.A(mesh = mesh_mod7,
                               loc = as.matrix(data_mod7[, c("x", "y")]))

stack_mod7 <- inla.stack(data = list(obs = data_mod7$cases),
                         A = list(A.est_mod7, 1),
                         effects = list(spatial_field = data.frame(s.index_mod7,
                                                                   Intercept = 1),
                                        data_mod7[, c(names(data_mod7) != "cases")]),
                         tag = "stdata")

f_mod7 <- obs ~ -1 + Intercept +
  mintro +
  offset(spde) +
  offset(logpop)

mod7_out <- inla(f_mod7,
                 data = inla.stack.data(stack_mod7),
                 family = "nbinomial",
                 control.predictor = list(A = inla.stack.A(stack_mod7)),
                 control.fixed = list(expand.factor.strategy = 'inla'),
                 control.compute = list(config = TRUE, dic = TRUE, waic = TRUE),
                 verbose = TRUE)

mintro_marg_mod7 <- mod7_out$marginals.fixed$mintro

set.seed(123)
beta_samples_mod7 <- inla.rmarginal(1000, mintro_marg_mod7)

effect_samples_mod7 <- sapply(beta_samples_mod7, function(b) exp(b * mintro_values))

mod7_df <- data.frame(
  mintro = mintro_values,
  effect = apply(effect_samples_mod7, 1, mean),
  lower  = apply(effect_samples_mod7, 1, quantile, 0.025),
  upper  = apply(effect_samples_mod7, 1, quantile, 0.975)
)

(fig4a2 <- mod7_df |>
  ggplot(aes(x = mintro, y = effect)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_hline(aes(yintercept = 1), lty = 2, color = "red") +
  labs(x = "% wMel introgression", y = "\nRatio of incidences in 2023") +
  theme_bw() +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1")) +
  coord_trans(y = "log1p") +
  theme(legend.position = "bottom",
        text = element_text(size = 15),
        plot.margin = unit(c(3, 4, 0, 4), "pt")) +
  coord_cartesian(ylim = c(0.1, 1.75)))

# -- Figure 4B: 2023 cases with high transmission field (2011-2013)
(fig4b1 <- dengues |>
  ggplot(aes(x = date)) +
  geom_rect(aes(xmin = as.Date("2011-01-01"), xmax = as.Date("2013-12-31"),
                ymin = -Inf, ymax = Inf), fill = "lightgray") +
  geom_rect(aes(xmin = as.Date("2023-01-01"), xmax = as.Date("2023-11-30"),
                ymin = -Inf, ymax = Inf), fill = "lightgray") +
  geom_line(aes(y = cases), size = 1.5) +
  annotate("text", x = mean(c(as.Date("2011-01-01"), as.Date("2013-12-31"))),
           y = 4000, label = "Spatial field", color = "gray10") +
  annotate("text", x = as.Date("2023-12-01"), y = 4000,
           label = "Cases\nin model", color = "gray10", lineheight = 0.85, hjust = 1) +
  labs(x = "", y = "Cases") +
  scale_x_date(date_labels = "%Y",
               breaks = seq.Date(as.Date("2010-01-01"), as.Date("2024-01-01"), by = "24 months"),
               expand = c(0.01, 0.01)) +
  scale_y_continuous(n.breaks = 3) +
  theme_bw() +
  theme(axis.text.x = element_text(),
        axis.title.x = element_blank(),
        text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5),
        plot.margin = unit(c(3, 4, 0, 4), "pt")) +
  coord_cartesian(ylim = c(0, 5000)))

data_mod8 <- data_mod_long2 |>
  mutate(year = as.numeric(as.character(year))) |>
  filter(year == 2023) |>
  mutate(logpop = log(pop)) |>
  ungroup() |>
  filter(!is.na(mintro)) |>
  left_join(spatfield_high, by = "id") |>
  dplyr::select(-year)

spatial.range <- 22000

mesh_mod8 <- inla.mesh.2d(loc = as.matrix(data_mod8[, c("x", "y")]),
                          max.edge = c(spatial.range / 5, spatial.range),
                          cutoff = 200,
                          offset = c(spatial.range / 5, diff(range(data_mod8[, "x"])) / 10))

spde_mod8 <- inla.spde2.matern(mesh = mesh_mod8)

s.index_mod8 <- inla.spde.make.index(name = "spatial.field",
                                     n.spde = spde_mod8$n.spde)

A.est_mod8 <- inla.spde.make.A(mesh = mesh_mod8,
                               loc = as.matrix(data_mod8[, c("x", "y")]))

stack_mod8 <- inla.stack(data = list(obs = data_mod8$cases),
                         A = list(A.est_mod8, 1),
                         effects = list(spatial_field = data.frame(s.index_mod8,
                                                                   Intercept = 1),
                                        data_mod8[, c(names(data_mod8) != "cases")]),
                         tag = "stdata")

f_mod8 <- obs ~ -1 + Intercept +
  mintro +
  offset(spde) +
  offset(logpop)

mod8_out <- inla(f_mod8,
                 data = inla.stack.data(stack_mod8),
                 family = "nbinomial",
                 control.predictor = list(A = inla.stack.A(stack_mod8)),
                 control.fixed = list(expand.factor.strategy = 'inla'),
                 control.compute = list(config = TRUE, dic = TRUE, waic = TRUE),
                 verbose = TRUE)

mintro_marg_mod8 <- mod8_out$marginals.fixed$mintro

set.seed(123)
beta_samples_mod8 <- inla.rmarginal(1000, mintro_marg_mod8)

effect_samples_mod8 <- sapply(beta_samples_mod8, function(b) exp(b * mintro_values))

mod8_df <- data.frame(
  mintro = mintro_values,
  effect = apply(effect_samples_mod8, 1, mean),
  lower  = apply(effect_samples_mod8, 1, quantile, 0.025),
  upper  = apply(effect_samples_mod8, 1, quantile, 0.975)
)

(fig4b2 <- mod8_df |>
  ggplot(aes(x = mintro, y = effect)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_hline(aes(yintercept = 1), lty = 2, color = "red") +
  labs(x = "% wMel introgression", y = "\nRatio of incidences in 2023") +
  theme_bw() +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1")) +
  coord_trans(y = "log1p") +
  theme(legend.position = "bottom",
        text = element_text(size = 15),
        plot.margin = unit(c(3, 4, 0, 4), "pt")) +
  coord_cartesian(ylim = c(0.1, 1.75)))

# -- Figure 4C: 2024 cases with low transmission field (2015-2016)
(fig4c1 <- dengues |>
  ggplot(aes(x = date)) +
  geom_rect(aes(xmin = as.Date("2015-01-01"), xmax = as.Date("2016-12-31"),
                ymin = -Inf, ymax = Inf), fill = "lightgray") +
  geom_rect(aes(xmin = as.Date("2024-01-01"), xmax = as.Date("2024-05-01"),
                ymin = -Inf, ymax = Inf), fill = "lightgray") +
  geom_line(aes(y = cases), size = 1.5) +
  annotate("text", x = mean(c(as.Date("2015-01-01"), as.Date("2016-12-31"))),
           y = 4000, label = "Spatial field", color = "gray10") +
  annotate("text", x = as.Date("2023-12-01"), y = 4000,
           label = "Cases\nin model", color = "gray10", lineheight = 0.85, hjust = 1) +
  labs(x = "", y = "Cases") +
  scale_x_date(date_labels = "%Y",
               breaks = seq.Date(as.Date("2010-01-01"), as.Date("2024-01-01"), by = "24 months"),
               expand = c(0.01, 0.01)) +
  scale_y_continuous(n.breaks = 3) +
  theme_bw() +
  theme(axis.text.x = element_text(),
        axis.title.x = element_blank(),
        text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5),
        plot.margin = unit(c(3, 4, 0, 4), "pt")) +
  coord_cartesian(ylim = c(0, 5000)))

data_mod9 <- data_mod_long2 |>
  mutate(year = as.numeric(as.character(year))) |>
  filter(year == 2024) |>
  mutate(logpop = log(pop)) |>
  ungroup() |>
  filter(!is.na(mintro)) |>
  left_join(spatfield_low, by = "id") |>
  dplyr::select(-year)

spatial.range <- 22000

mesh_mod9 <- inla.mesh.2d(loc = as.matrix(data_mod9[, c("x", "y")]),
                          max.edge = c(spatial.range / 5, spatial.range),
                          cutoff = 200,
                          offset = c(spatial.range / 5, diff(range(data_mod9[, "x"])) / 10))

spde_mod9 <- inla.spde2.matern(mesh = mesh_mod9)

s.index_mod9 <- inla.spde.make.index(name = "spatial.field",
                                     n.spde = spde_mod9$n.spde)

A.est_mod9 <- inla.spde.make.A(mesh = mesh_mod9,
                               loc = as.matrix(data_mod9[, c("x", "y")]))

stack_mod9 <- inla.stack(data = list(obs = data_mod9$cases),
                         A = list(A.est_mod9, 1),
                         effects = list(spatial_field = data.frame(s.index_mod9,
                                                                   Intercept = 1),
                                        data_mod9[, c(names(data_mod9) != "cases")]),
                         tag = "stdata")

f_mod9 <- obs ~ -1 + Intercept +
  mintro +
  offset(spde) +
  offset(logpop)

mod9_out <- inla(
  f_mod9,
  data = inla.stack.data(stack_mod9),
  family = "nbinomial",
  control.predictor = list(A = inla.stack.A(stack_mod9)),
  control.fixed = list(expand.factor.strategy = 'inla'),
  control.compute = list(config = TRUE, dic = TRUE, waic = TRUE),
  verbose = TRUE
  )

mintro_marg_mod9 <- mod9_out$marginals.fixed$mintro

set.seed(123)
beta_samples_mod9 <- inla.rmarginal(1000, mintro_marg_mod9)

effect_samples_mod9 <- sapply(beta_samples_mod9, function(b) exp(b * mintro_values))

mod9_df <- data.frame(
  mintro = mintro_values,
  effect = apply(effect_samples_mod9, 1, mean),
  lower  = apply(effect_samples_mod9, 1, quantile, 0.025),
  upper  = apply(effect_samples_mod9, 1, quantile, 0.975)
)

(fig4c2 <- mod9_df |>
  ggplot(aes(x = mintro, y = effect)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_hline(aes(yintercept = 1), lty = 2, color = "red") +
  labs(x = "% wMel introgression", y = "\nRatio of incidences in 2024") +
  theme_bw() +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1")) +
  coord_trans(y = "log1p") +
  theme(legend.position = "bottom",
        text = element_text(size = 15),
        plot.margin = unit(c(3, 4, 0, 4), "pt")) +
  coord_cartesian(ylim = c(0.1, 1.75)))

# -- Figure 4D: 2024 cases with high transmission field (2011-2013)
(fig4d1 <- dengues |>
  ggplot(aes(x = date)) +
  geom_rect(aes(xmin = as.Date("2011-01-01"), xmax = as.Date("2013-12-31"),
                ymin = -Inf, ymax = Inf), fill = "lightgray") +
  geom_rect(aes(xmin = as.Date("2024-01-01"), xmax = as.Date("2024-05-01"),
                ymin = -Inf, ymax = Inf), fill = "lightgray") +
  geom_line(aes(y = cases), size = 1.5) +
  annotate("text", x = mean(c(as.Date("2011-01-01"), as.Date("2013-12-31"))),
           y = 4000, label = "Spatial field", color = "gray10") +
  annotate("text", x = as.Date("2023-12-01"), y = 4000,
           label = "Cases\nin model", color = "gray10", lineheight = 0.85, hjust = 1) +
  labs(x = "", y = "Cases") +
  scale_x_date(date_labels = "%Y",
               breaks = seq.Date(as.Date("2010-01-01"), as.Date("2024-01-01"), by = "24 months"),
               expand = c(0.01, 0.01)) +
  scale_y_continuous(n.breaks = 3) +
  theme_bw() +
  theme(axis.text.x = element_text(),
        axis.title.x = element_blank(),
        text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5),
        plot.margin = unit(c(3, 4, 0, 4), "pt")) +
  coord_cartesian(ylim = c(0, 5000)))

data_mod10 <- data_mod_long2 |>
  mutate(year = as.numeric(as.character(year))) |>
  filter(year == 2024) |>
  mutate(logpop = log(pop)) |>
  ungroup() |>
  filter(!is.na(mintro)) |>
  left_join(spatfield_high, by = "id") |>
  dplyr::select(-year)

spatial.range <- 22000

mesh_mod10 <- inla.mesh.2d(loc = as.matrix(data_mod10[, c("x", "y")]),
                           max.edge = c(spatial.range / 5, spatial.range),
                           cutoff = 200,
                           offset = c(spatial.range / 5, diff(range(data_mod10[, "x"])) / 10))

spde_mod10 <- inla.spde2.matern(mesh = mesh_mod10)

s.index_mod10 <- inla.spde.make.index(name = "spatial.field",
                                      n.spde = spde_mod10$n.spde)

A.est_mod10 <- inla.spde.make.A(mesh = mesh_mod10,
                                loc = as.matrix(data_mod10[, c("x", "y")]))

stack_mod10 <- inla.stack(data = list(obs = data_mod10$cases),
                          A = list(A.est_mod10, 1),
                          effects = list(spatial_field = data.frame(s.index_mod10,
                                                                    Intercept = 1),
                                         data_mod10[, c(names(data_mod10) != "cases")]),

                          tag = "stdata")

f_mod10 <- obs ~ -1 + Intercept +
  mintro +
  offset(spde) +
  offset(logpop)

mod10_out <- inla(f_mod10,
                  data = inla.stack.data(stack_mod10),
                  family = "nbinomial",
                  control.predictor = list(A = inla.stack.A(stack_mod10)),
                  control.fixed = list(expand.factor.strategy = 'inla'),
                  control.compute = list(config = TRUE, dic = TRUE, waic = TRUE),
                  verbose = TRUE)

mintro_marg_mod10 <- mod10_out$marginals.fixed$mintro

set.seed(123)
beta_samples_mod10 <- inla.rmarginal(1000, mintro_marg_mod10)

effect_samples_mod10 <- sapply(beta_samples_mod10, function(b) exp(b * mintro_values))

mod10_df <- data.frame(
  mintro = mintro_values,
  effect = apply(effect_samples_mod10, 1, mean),
  lower  = apply(effect_samples_mod10, 1, quantile, 0.025),
  upper  = apply(effect_samples_mod10, 1, quantile, 0.975)
)

(fig4d2 <- mod10_df |>
  ggplot(aes(x = mintro, y = effect)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_hline(aes(yintercept = 1), lty = 2, color = "red") +
  labs(x = "% wMel introgression", y = "\nRatio of incidences in 2024") +
  theme_bw() +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1")) +
  coord_trans(y = "log1p") +
  theme(legend.position = "bottom",
        text = element_text(size = 15),
        plot.margin = unit(c(3, 4, 0, 4), "pt")) +
  coord_cartesian(ylim = c(0.1, 1.75)))

# -- Time series
(fig4a_v2 <- dengues |>
  ggplot(aes(x = date)) +
  geom_rect(aes(xmin = as.Date("2015-01-01"), xmax = as.Date("2016-12-31"),
                ymin = -Inf, ymax = Inf), fill = "lightgray") +
  geom_line(aes(y = cases), size = 1.5) +
  labs(x = "", y = "Monthly\ndengue cases") +
  scale_x_date(date_labels = "%Y",
               breaks = seq.Date(as.Date("2010-01-01"), as.Date("2024-01-01"), by = "24 months"),
               expand = c(0.01, 0.01)) +
  scale_y_continuous(n.breaks = 3) +
  theme_bw() +
  theme(axis.text.x = element_text(),
        axis.title.x = element_blank(),
        text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5),
        plot.margin = unit(c(3, 4, 0, 4), "pt")) +
  coord_cartesian(ylim = c(0, 5000)))

(fig4b_v2 <- dengues |>
  ggplot(aes(x = date)) +
  geom_rect(aes(xmin = as.Date("2011-01-01"), xmax = as.Date("2013-12-31"),
                ymin = -Inf, ymax = Inf), fill = "lightgray") +
  geom_line(aes(y = cases), size = 1.5) +
  labs(x = "", y = "Monthly\ndengue cases") +
  scale_x_date(date_labels = "%Y",
               breaks = seq.Date(as.Date("2010-01-01"), as.Date("2024-01-01"), by = "24 months"),
               expand = c(0.01, 0.01)) +
  scale_y_continuous(n.breaks = 3) +
  theme_bw() +
  theme(axis.text.x = element_text(),
        axis.title.x = element_blank(),
        text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5),
        plot.margin = unit(c(3, 4, 0, 4), "pt")) +
  coord_cartesian(ylim = c(0, 5000)))

# -- Complete figure 4
fig4_layout <- "
AAAAAABBBBBB
AAAAAABBBBBB
CCCCCCDDDDDD
EEEEEEFFFFFF
EEEEEEFFFFFF
EEEEEEFFFFFF
GGGGGGHHHHHH
GGGGGGHHHHHH
GGGGGGHHHHHH
"

(fig4 <- (free(fig4_lowinset) + free(fig4_highinset)) +
    fig4a_v2 + fig4b_v2 +
    fig4a2 + fig4b2 +
    fig4c2 + fig4d2 +
    plot_layout(design = fig4_layout) +
    plot_annotation(tag_levels = list(c("A", "B", "", "", "C", "D", "E", "F"))))

ggsave("outputs/2026-07-07_fig4.pdf", plot = fig4, width = 10, height = 11)
# ggsave("outputs/2026-07-07_fig4.png", plot = fig4, width = 10, height = 11)

# ---- SX: Tileplot of introgression -------------------------------------
intersections_sep <- st_intersection(grid, release_areas) |>
  mutate(intersection_area = st_area(geometry)) |>
  left_join(grid_prelim |> st_drop_geometry(), by = "id")

grid_sep <- intersections_sep |>
  mutate(perc_intersection = as.numeric(intersection_area / cell_area)) |>
  st_drop_geometry() |>
  arrange(id, desc(perc_intersection)) |>
  group_by(id) |>
  summarise(release_area = REGIAO_2[1])

trap_long_grid_rj <- trap_long2 |>
  st_drop_geometry() |>
  filter(!is.na(trap_id)) |>
  group_by(id, year, month) |>
  summarise(mintro = mean(mean_intro, na.rm = TRUE),
            .groups = "drop") |>
  left_join(grid_sep, by = "id") |>
  mutate(monthyear = as.Date(paste(year, month, "01", sep = "-")))

trap_long_grid_rj |>
  ggplot(aes(x = monthyear, y = factor(id))) +
  geom_tile(aes(fill = mintro*100), width = 31) +
  scale_fill_distiller(palette = "Greens", na.value = "transparent", direction = 1) +
  facet_grid(release_area~., scales = "free_y", space = "free_y") +
  theme_bw() +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(x = "", fill = "% wMel") +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.background = element_rect(fill = "lightgray"))

xr    <- range(as.Date(paste(data_mod_long2$year, data_mod_long2$month, "01", sep = "-")))
xlim  <- c(xr[1] - 20, xr[2] + 20)

periods <- data.frame(
  xmin  = as.Date(c("2017-02-01", "2020-05-12")),
  xmax  = as.Date(c("2020-12-04", "2024-05-03")),
  row   = c(1, 0),                       # BG on top row, ovitraps on lower row
  label = c("BG-traps", "Ovitraps")
)
periods$xmid <- periods$xmin + as.numeric(periods$xmax - periods$xmin) / 2

top <- ggplot(periods) +
  geom_segment(aes(x = xmin, xend = xmax, y = row, yend = row), linewidth = 0.5) +
  geom_segment(aes(x = xmin, xend = xmin, y = row, yend = row - 0.15), linewidth = 0.5) +
  geom_segment(aes(x = xmax, xend = xmax, y = row, yend = row - 0.15), linewidth = 0.5) +
  geom_text(aes(x = xmid, y = row + 0.15, label = label), vjust = 0, size = 3.2) +
  scale_x_date(limits = xlim, expand = expansion(mult = 0.05)) +
  scale_y_continuous(limits = c(-0.4, 1.6)) +
  theme_void()

bottom <- data_mod_long2 |>
  mutate(monthyear = as.Date(paste(year, month, "01", sep = "-"))) |>
  left_join(grid_sep, by = "id") |>
  ggplot(aes(x = monthyear, y = factor(id))) +
  geom_tile(aes(fill = mintro * 100), width = 31) +
  scale_fill_distiller(palette = "Greens", na.value = "transparent", direction = 1) +
  facet_grid(release_area ~ ., scales = "free_y", space = "free_y") +
  theme_bw() +
  scale_x_date(limits = xlim, date_breaks = "1 year", date_labels = "%Y",
               expand = expansion(mult = 0.05)) +
  labs(x = "", fill = "% wMel") +
  theme(axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.background   = element_rect(fill = "gray90"),
        legend.position = "bottom")

top / bottom + plot_layout(heights = c(1.4, 20))

ggsave("outputs/2026-09-07_tileplot-introgression.png", width = 10, height = 8)
ggsave("outputs/2026-09-07_tileplot-introgression.pdf", width = 10, height = 8)

# ---- Figure S2 ---------------------------------------------------------
# Get 3-mo rolling average of ovi-estimated introgression
trapsss <- trap_raw |>
  filter(type == "Oviwmel" & successful == "TRUE" & !is.na(target_species_count) & target_species_count > 0) |>
  mutate(intro = screening_wmel_aeg/target_species_count,
         year = year(collected_at),
         month = month(collected_at)) |>
  group_by(trap_id, latitude, longitude, year, month) |>
  summarise(mean_intro = mean(intro, na.rm = TRUE),
            wmel = mean(screening_wmel_aeg, na.rm = TRUE),
            aedes = mean(target_species_count, na.rm = TRUE)) |>
  st_as_sf(coords = c("longitude", "latitude"), crs = "EPSG:4326 - WGS 84", remove = FALSE) |>
  st_transform(crs = 31983)

trap <- st_join(grid, trapsss, join = st_contains)

trap_grid <- trap |>
  st_drop_geometry() |>
  filter(!is.na(trap_id)) |>
  group_by(id, year, month) |>
  group_by(id, year, month) |>
  summarise(weighted_intro = sum(mean_intro * aedes, na.rm = TRUE) / sum(aedes, na.rm = TRUE)) |>
  ungroup() |>
  arrange(id, year, month) |>
  group_by(id) |>
  mutate(rolling_intro = zoo::rollapply(data = weighted_intro,
                                   width = 3,
                                   FUN = mean,
                                   align = "center",
                                   fill = NA,
                                   na.rm = TRUE)) |>
  ungroup()

# Get BG-trap estimate
bg_raw <- trap_raw |>
  filter(type == "BG" & successful == "TRUE" & collected_at > as.POSIXct("2018-12-31") & !is.na(target_species_count) & target_species_count > 0) |>
  mutate(intro = screening_wmel_aeg/target_species_count) |>
  group_by(trap_id, latitude, longitude) |>
  summarise(wmean_intro = weighted.mean(intro, w = total_aegypti_caught, na.rm = TRUE),
            mean_intro = mean(intro, na.rm = TRUE)) |>
  st_as_sf(coords = c("longitude", "latitude"), crs = "EPSG:4326 - WGS 84", remove = FALSE) |>
  st_transform(crs = 31983)

bg <- st_join(grid, bg_raw, join = st_contains)

bg_grid <- bg |>
  st_drop_geometry() |>
  filter(!is.na(trap_id)) |>
  group_by(id) |>
  summarise(mintro_bg = mean(mean_intro, na.rm = TRUE),
            wmintro_bg = mean(wmean_intro, na.rm = TRUE))

compare_traps <- trap_grid |>
  left_join(bg_grid, by = "id") |>
  na.omit() |>
  mutate(mintro_bg_quartile = ntile(mintro_bg, 4)) |>
  mutate(mintro_bg_quartile = factor(mintro_bg_quartile,
                                     levels = 1:4,
                                     labels = c("Q1", "Q2", "Q3", "Q4"))) |>
  group_by(year) |>
  mutate(rolling_intro_quartile = ntile(rolling_intro, 4)) |>
  ungroup() |>
  mutate(rolling_intro_quartile = factor(rolling_intro_quartile,
                                         levels = 1:4,
                                         labels = c("Q1", "Q2", "Q3", "Q4")))

(fig_s2 <- compare_traps |>
  ggplot(aes(x = rolling_intro, y = mintro_bg)) +
  geom_point(alpha = 0.2) +
  geom_smooth(method = "lm") +
  facet_wrap(~year) +
  labs(x = "Introgression estimated from ovitraps\n(3-month rolling mean)", y = "Introgression estimated from BG-traps\n(mean 2019-2020)") +
  theme_bw() +
  scale_x_continuous(labels = c("0", "0.25", "0.50", "0.75", "1")) +
  scale_y_continuous(labels = c("0", "0.25", "0.50", "0.75", "1")) +
  theme(text = element_text(size = 16)) +
  coord_equal())

# ggsave("outputs/figs1-sep29.pdf", plot = fig_s2, width = 7, height = 5)
# ggsave("outputs/figs1-sep29.png", plot = fig_s2, width = 7, height = 5)

cor_labels <- compare_traps |>
  group_by(year) |>
  summarise(r = cor(rolling_intro, mintro_bg, method = "pearson", use = "complete.obs"),
            p = cor.test(rolling_intro, mintro_bg, method = "pearson")$p.value) |>
  mutate(p = if_else(p < 0.001, "<0.001", as.character(signif(p, 2))),
         label = paste0("r = ", round(r, 2),
                        "\np = ", p))

(fig_s2 <- compare_traps |>
  ggplot(aes(x = rolling_intro, y = mintro_bg)) +
  geom_point(alpha = 0.2) +
  geom_smooth(method = "lm") +
  facet_wrap(~year) +
  geom_text(
    data = cor_labels,
    aes(x = 0.05, y = 0.95, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 5
  ) +
  labs(
    x = "Introgression estimated from ovitraps\n(3-month rolling mean)",
    y = "Introgression estimated from BG-traps\n(mean 2019-2020)"
  ) +
  theme_bw() +
  scale_x_continuous(labels = c("0", "0.25", "0.50", "0.75", "1")) +
  scale_y_continuous(labels = c("0", "0.25", "0.50", "0.75", "1")) +
  theme(text = element_text(size = 16)) +
  coord_equal())

# ggsave("outputs/2026-04-13_figS1.pdf", plot = fig_s2, width = 7, height = 5)
# ggsave("outputs/2026-04-13_figS1.png", plot = fig_s2, width = 7, height = 5)

# ---- Figure S3 ---------------------------------------------------------
(fig_s3a <- grid_sf_modlong2 |>
  ggplot() +
  #geom_sf(data = release_area, fill = "gray25", color = "gray25") +
  geom_tile(aes(x = x, y = y, fill = exp(mean))) +
  scale_fill_gradient2(low = "#25489E", mid = "white", high = "#CE263D", midpoint = 1, limits = c(0.3, 3.7)) +
  geom_sf(data = release_areas_j, fill = NA, color = "black") +
  theme_void() +
  labs(fill = "Ratio of\nincidences") +
  theme(text = element_text(size = 16)))

(fig_s3b <- grid_sf_modlong2b |>
  ggplot() +
  #geom_sf(data = release_area, fill = "gray25", color = "gray25") +
  geom_tile(aes(x = x, y = y, fill = exp(mean))) +
  scale_fill_gradient2(low = "#25489E", mid = "white", high = "#CE263D", midpoint = 1, limits = c(0.3, 3.7)) +
  geom_sf(data = release_areas_j, fill = NA, color = "black") +
  theme_void() +
  labs(fill = "Ratio of\nincidences") +
  theme(text = element_text(size = 16)))

(fig_s3 <- fig_s3a + fig_s3b + plot_layout(widths = c(1, 1))  +
   plot_annotation(tag_levels = "A",
                  theme = theme(
                    plot.margin = margin(1, 3, 1, 3),
                    plot.tag.position = c(0, 0))) +
   plot_layout(guides = "collect"))

# ggsave("outputs/2026-07-06_figs3.png", plot = fig_s3, height = 2.5, width = 8)
# ggsave("outputs/2026-07-06_figs3.pdf", plot = fig_s3, height = 2.5, width = 8)

# ---- Figure S4: Ovi models ---------------------------------------------

# -- Monthly
trap_ovi_month <- trap_raw |>
  filter(type == "Oviwmel" & successful == "TRUE" & !is.na(target_species_count) & target_species_count > 0) |>
  mutate(intro = screening_wmel_aeg/target_species_count,
         year = year(collected_at),
         month = month(collected_at)) |>
  group_by(trap_id, latitude, longitude, year, month) |>
  summarise(mean_intro = mean(intro, na.rm = TRUE),
            wmel = mean(screening_wmel_aeg, na.rm = TRUE),
            aedes = mean(target_species_count, na.rm = TRUE)) |>
  st_as_sf(coords = c("longitude", "latitude"), crs = "EPSG:4326 - WGS 84", remove = FALSE) |>
  st_transform(crs = 31983)

trap_ovi_m <- st_join(grid, trap_ovi_month, join = st_contains)

trap_ovi_grid_m <- trap_ovi_m |>
  st_drop_geometry() |>
  filter(!is.na(trap_id)) |>
  group_by(id, year, month) |>
  summarise(mintro = mean(mean_intro, na.rm = TRUE),
            wmel = sum(wmel, na.rm = TRUE),
            aedes = sum(aedes, na.rm = TRUE))

dengue_grid_s3a <- st_join(grid, dengue_raw1, join = st_contains) |>
  filter(year >= 2020)

dengue_gridmod_s3a <- dengue_grid_s3a |>
  st_drop_geometry() |>
  filter(!is.na(year)) |>
  group_by(year, month, id) |>
  summarise(cases = n()) |>
  ungroup() |>
  right_join(grid |> st_drop_geometry(), by = "id") |>
  tidyr::complete(year, month, id, fill = list(cases = 0)) |>
  filter(year != 0)

data_s3a <- grid |>
  right_join(dengue_gridmod_s3a, by = "id") |>
  left_join(popdens_grid, by = "id") |>
  left_join(trap_ovi_grid_m, by = c("id", "year", "month")) |>
  filter(year >= 2020) |>
  rename(pop2 = sum) |>
  mutate(pop = if_else(pop2 < 50, 0, pop2),
         inc = cases/pop*1000) |>
  filter(st_intersects(geometry, release_areas_j, sparse = FALSE)[,1]) |>
  filter(!is.na(mintro)) |>
  mutate(geometry = st_centroid(geometry),
         x = st_coordinates(geometry)[,1],
         y = st_coordinates(geometry)[,2]) |>
  st_drop_geometry() |>
  filter(pop != 0) |>
  mutate(logpop = log(pop),
         year = factor(year)) |>
  group_by(year) |>
  mutate(nyear = cur_group_id()) |>
  ungroup()

table(data_s3a$year)

spatial.range <- 22000

mesh_s3a <- inla.mesh.2d(loc = as.matrix(data_s3a[, c("x", "y")]),
                        max.edge = c(spatial.range / 5, spatial.range),
                        cutoff = 200,
                        offset = c(spatial.range / 5, diff(range(data_s3a[, "x"])) / 10))

spde_s3a <- inla.spde2.matern(mesh = mesh_s3a)

s.index_s3a <- inla.spde.make.index(name = "spatial.field",
                                  n.spde = spde_s3a$n.spde)

A.est_s3a <- inla.spde.make.A(mesh = mesh_s3a,
                            loc = as.matrix(data_s3a[, c("x", "y")]))

stack_s3a <- inla.stack(
  data = list(obs = data_s3a$cases),
  A = list(A.est_s3a, 1),
  effects = list(spatial_field = data.frame(s.index_s3a, Intercept = 1),
                 data_s3a[, c(names(data_s3a) != "cases")]),
  tag = "stdata"
)

f_s3a <- obs ~ -1 + Intercept +
  year:mintro +
  f(nyear, model = "ar1") +
  f(spatial.field, model = spde_s3a) +
  offset(logpop)

s3a_out <- inla(
  f_s3a,
  data = inla.stack.data(stack_s3a),
  family = "nbinomial",
  control.predictor = list(A = inla.stack.A(stack_s3a)),
  control.fixed = list(expand.factor.strategy = 'inla'),
  control.compute = list(config = TRUE, dic = TRUE, waic = TRUE),
  verbose = TRUE
  )

summary(s3a_out)
exp(as.data.frame(inla.zmarginal(s3a_out$marginals.fixed$`year2020:mintro`)))
exp(as.data.frame(inla.zmarginal(s3a_out$marginals.fixed$`year2021:mintro`)))
exp(as.data.frame(inla.zmarginal(s3a_out$marginals.fixed$`year2022:mintro`)))
exp(as.data.frame(inla.zmarginal(s3a_out$marginals.fixed$`year2023:mintro`)))
exp(as.data.frame(inla.zmarginal(s3a_out$marginals.fixed$`year2024:mintro`)))

mintro_int_margs_s3a <- list(
  `2020` = s3a_out$marginals.fixed$`year2020:mintro`,
  `2021` = s3a_out$marginals.fixed$`year2021:mintro`,
  `2022` = s3a_out$marginals.fixed$`year2022:mintro`,
  `2023` = s3a_out$marginals.fixed$`year2023:mintro`,
  `2024` = s3a_out$marginals.fixed$`year2024:mintro`
)

set.seed(123)
s3a_df <- bind_rows(lapply(names(mintro_int_margs_s3a), function(yr) {
  beta_year_draws <- inla.rmarginal(1000, mintro_int_margs_s3a[[yr]])
  effect_samples <- sapply(beta_year_draws, function(b) exp(b * mintro_values))

  data.frame(
    year   = as.integer(yr),
    mintro = mintro_values,
    effect = apply(effect_samples, 1, mean),
    lower  = apply(effect_samples, 1, quantile, probs = 0.025),
    upper  = apply(effect_samples, 1, quantile, probs = 0.975)
  )
}))

(s3a_plot <- s3a_df |>
  ggplot(aes(x = mintro, y = effect)) +
  geom_line(lwd = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_hline(aes(yintercept = 1), lty = 2, color = "red") +
  labs(x = "% wMel introgression", y = "Ratio of incidences",
       title = "B. Monthly cases") +
  theme_bw() +
  facet_wrap(~year, ncol = 5) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1")) +
  coord_trans(y = "log1p") +
  theme(legend.position = "bottom",
        text = element_text(size = 15)))

# -- Yearly
trap_ovi_y <- trap_raw |>
  filter(type == "Oviwmel" & successful == "TRUE" & !is.na(target_species_count) & target_species_count > 0) |>
  mutate(intro = screening_wmel_aeg/target_species_count,
         year = year(collected_at)) |>
  group_by(trap_id, latitude, longitude, year) |>
  summarise(mean_intro = mean(intro, na.rm = TRUE),
            wmel = mean(screening_wmel_aeg, na.rm = TRUE),
            aedes = mean(target_species_count, na.rm = TRUE)) |>
  st_as_sf(coords = c("longitude", "latitude"), crs = "EPSG:4326 - WGS 84", remove = FALSE) |>
  st_transform(crs = 31983)

trap_ovi_y <- st_join(grid, trap_ovi_y, join = st_contains)

trap_ovi_grid_y <- trap_ovi_y |>
  st_drop_geometry() |>
  filter(!is.na(trap_id)) |>
  group_by(id, year) |>
  summarise(mintro = mean(mean_intro, na.rm = TRUE),
            wmel = sum(wmel, na.rm = TRUE),
            aedes = sum(aedes, na.rm = TRUE))

dengue_gridmod_s3b <- dengue_grid_s3a |>
  st_drop_geometry() |>
  filter(!is.na(year)) |>
  group_by(year, id) |>
  summarise(cases = n()) |>
  ungroup() |>
  right_join(grid |> st_drop_geometry(), by = "id") |>
  tidyr::complete(year, id, fill = list(cases = 0)) |>
  filter(year != 0)

data_s3b <- grid |>
  right_join(dengue_gridmod_s3b, by = "id") |>
  left_join(popdens_grid, by = "id") |>
  left_join(trap_ovi_grid_y, by = c("id", "year")) |>
  filter(year >= 2020) |>
  rename(pop2 = sum) |>
  mutate(pop = if_else(pop2 < 50, 0, pop2),
         inc = cases/pop*1000) |>
  filter(st_intersects(geometry, release_areas_j, sparse = FALSE)[,1]) |>
  filter(!is.na(mintro)) |>
  mutate(geometry = st_centroid(geometry),
         x = st_coordinates(geometry)[,1],
         y = st_coordinates(geometry)[,2]) |>
  st_drop_geometry() |>
  filter(pop != 0) |>
  mutate(logpop = log(pop),
         year = factor(year)) |>
  group_by(year) |>
  mutate(nyear = cur_group_id()) |>
  ungroup()

spatial.range <- 22000

mesh_s3b <- inla.mesh.2d(loc = as.matrix(data_s3b[, c("x", "y")]),
                        max.edge = c(spatial.range / 5, spatial.range),
                        cutoff = 200,
                        offset = c(spatial.range / 5, diff(range(data_s3b[, "x"])) / 10))

spde_s3b <- inla.spde2.matern(mesh = mesh_s3b)

s.index_s3b <- inla.spde.make.index(name = "spatial.field",
                                  n.spde = spde_s3b$n.spde)

A.est_s3b <- inla.spde.make.A(mesh = mesh_s3b,
                            loc = as.matrix(data_s3b[, c("x", "y")]))

stack_s3b <- inla.stack(
  data = list(obs = data_s3b$cases),
  A = list(A.est_s3b, 1),
  effects = list(spatial_field = data.frame(s.index_s3b, Intercept = 1),
                 data_s3b[, c(names(data_s3b) != "cases")]),
  tag = "stdata"
)

f_s3b <- obs ~ -1 + Intercept +
  year:mintro +
  f(nyear, model = "ar1") +
  f(spatial.field, model = spde_s3b) +
  offset(logpop)

s3b_out <- inla(
  f_s3b,
  data = inla.stack.data(stack_s3b),
  family = "nbinomial",
  control.predictor = list(A = inla.stack.A(stack_s3b)),
  control.fixed = list(expand.factor.strategy = 'inla'),
  control.compute = list(config = TRUE, dic = TRUE, waic = TRUE),
  verbose = TRUE
  )

summary(s3b_out)
exp(as.data.frame(inla.zmarginal(s3b_out$marginals.fixed$`year2020:mintro`)))
exp(as.data.frame(inla.zmarginal(s3b_out$marginals.fixed$`year2021:mintro`)))
exp(as.data.frame(inla.zmarginal(s3b_out$marginals.fixed$`year2022:mintro`)))
exp(as.data.frame(inla.zmarginal(s3b_out$marginals.fixed$`year2023:mintro`)))
exp(as.data.frame(inla.zmarginal(s3b_out$marginals.fixed$`year2024:mintro`)))

mintro_int_margs_s3b <- list(
  `2020` = s3b_out$marginals.fixed$`year2020:mintro`,
  `2021` = s3b_out$marginals.fixed$`year2021:mintro`,
  `2022` = s3b_out$marginals.fixed$`year2022:mintro`,
  `2023` = s3b_out$marginals.fixed$`year2023:mintro`,
  `2024` = s3b_out$marginals.fixed$`year2024:mintro`
)

set.seed(123)
s3b_df <- bind_rows(lapply(names(mintro_int_margs_s3b), function(yr) {
  beta_year_draws <- inla.rmarginal(1000, mintro_int_margs_s3b[[yr]])
  effect_samples <- sapply(beta_year_draws, function(b) exp(b * mintro_values))

  data.frame(
    year   = as.integer(yr),
    mintro = mintro_values,
    effect = apply(effect_samples, 1, mean),
    lower  = apply(effect_samples, 1, quantile, probs = 0.025),
    upper  = apply(effect_samples, 1, quantile, probs = 0.975)
  )
}))

(s3b_plot <- s3b_df |>
  ggplot(aes(x = mintro, y = effect)) +
  geom_line(lwd = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_hline(aes(yintercept = 1), lty = 2, color = "red") +
  labs(x = "% wMel introgression", y = "Ratio of incidences",
       title = "A. Yearly cases") +
  theme_bw() +
  facet_wrap(~year, ncol = 5) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1")) +
  coord_trans(y = "log1p") +
  theme(legend.position = "bottom",
        text = element_text(size = 15)))

# -- Cumulative
trap_ovi_f <- trap_raw |>
  filter(type == "Oviwmel" & successful == "TRUE" & !is.na(target_species_count) & target_species_count > 0) |>
  mutate(intro = screening_wmel_aeg/target_species_count,
         year = year(collected_at)) |>
  group_by(trap_id, latitude, longitude) |>
  summarise(mean_intro = mean(intro, na.rm = TRUE),
            wmel = mean(screening_wmel_aeg, na.rm = TRUE),
            aedes = mean(target_species_count, na.rm = TRUE)) |>
  st_as_sf(coords = c("longitude", "latitude"), crs = "EPSG:4326 - WGS 84", remove = FALSE) |>
  st_transform(crs = 31983)

trap_ovi_f <- st_join(grid, trap_ovi_f, join = st_contains)

trap_ovi_grid_f <- trap_ovi_f |>
  st_drop_geometry() |>
  filter(!is.na(trap_id)) |>
  group_by(id) |>
  summarise(mintro = mean(mean_intro, na.rm = TRUE),
            wmel = sum(wmel, na.rm = TRUE),
            aedes = sum(aedes, na.rm = TRUE))

dengue_gridmod_s3c <- dengue_grid_s3a |>
  st_drop_geometry() |>
  filter(!is.na(year) & year >= 2020) |>
  group_by(id) |>
  summarise(cases = n()) |>
  ungroup() |>
  right_join(grid |> st_drop_geometry(), by = "id") |>
  tidyr::complete(id, fill = list(cases = 0))

data_s3c <- grid |>
  right_join(dengue_gridmod_s3c, by = "id") |>
  left_join(popdens_grid, by = "id") |>
  left_join(trap_ovi_grid_f, by = c("id")) |>
  rename(pop2 = sum) |>
  mutate(pop = if_else(pop2 < 50, 0, pop2),
         inc = cases/pop*1000) |>
  filter(st_intersects(geometry, release_areas_j, sparse = FALSE)[,1]) |>
  filter(!is.na(mintro)) |>
  mutate(geometry = st_centroid(geometry),
         x = st_coordinates(geometry)[,1],
         y = st_coordinates(geometry)[,2]) |>
  st_drop_geometry() |>
  filter(pop != 0) |>
  mutate(logpop = log(pop))

spatial.range <- 22000

mesh_s3c <- inla.mesh.2d(loc = as.matrix(data_s3c[, c("x", "y")]),
                        max.edge = c(spatial.range / 5, spatial.range),
                        cutoff = 200,
                        offset = c(spatial.range / 5, diff(range(data_s3c[, "x"])) / 10))

spde_s3c <- inla.spde2.matern(mesh = mesh_s3c)

s.index_s3c <- inla.spde.make.index(name = "spatial.field",
                                  n.spde = spde_s3c$n.spde)

A.est_s3c <- inla.spde.make.A(mesh = mesh_s3c,
                            loc = as.matrix(data_s3c[, c("x", "y")]))

stack_s3c <- inla.stack(
  data = list(obs = data_s3c$cases),
  A = list(A.est_s3c, 1),
  effects = list(spatial_field = data.frame(s.index_s3c, Intercept = 1),
                 data_s3c[, c(names(data_s3c) != "cases")]),
  tag = "stdata"
)

f_s3c <- obs ~ -1 + Intercept +
  mintro +
  f(spatial.field, model = spde_s3c) +
  offset(logpop)

s3c_out <- inla(
  f_s3c,
  data = inla.stack.data(stack_s3c),
  family = "nbinomial",
  control.predictor = list(A = inla.stack.A(stack_s3c)),
  control.fixed = list(expand.factor.strategy = 'inla'),
  control.compute = list(config = TRUE, dic = TRUE, waic = TRUE),
  verbose = TRUE
  )

summary(s3c_out)
exp(as.data.frame(inla.zmarginal(s3c_out$marginals.fixed$mintro)))

mintro_marg_s3c <- s3c_out$marginals.fixed$mintro

set.seed(123)
beta_samples <- inla.rmarginal(1000, mintro_marg_s3c)

effect_samples <- sapply(beta_samples, function(b) exp(b * mintro_values))

s3c_df <- data.frame(
  mintro = mintro_values,
  effect = apply(effect_samples, 1, mean),
  lower  = apply(effect_samples, 1, quantile, 0.025),
  upper  = apply(effect_samples, 1, quantile, 0.975)
)

(s3c_plot <- s3c_df |>
  ggplot(aes(x = mintro, y = effect)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_hline(aes(yintercept = 1), lty = 2, color = "red") +
  labs(x = "% wMel introgression", y = "Ratio of incidences",
       title = "C. Cumulative cases") +
  theme_bw() +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1")) +
  coord_trans(y = "log1p") +
  theme(legend.position = "bottom",
        text = element_text(size = 15)))

# -- Full figure
figs3_layout <- "
AAAAAAAAAA
AAAAAAAAAA
BBBBBBBBBB
BBBBBBBBBB
CC########
CC########
"

(figs3 <- s3b_plot + s3a_plot + s3c_plot +
    plot_layout(design = figs3_layout))

# ggsave("outputs/figs3-jan28.png", plot = figs3, height = 9, width = 10)
# ggsave("outputs/2026-03-13_figS3.pdf", plot = figs3, height = 9, width = 10)

# ---- Table SX: Model comparison ----------------------------------------

# -- Model from 2019-2020

# -- Cumulative
intro_1920 <- trap_raw |>
  filter(type == "BG" &
           successful == "TRUE" &
           collected_at > as.POSIXct("2018-12-31") &
           !is.na(target_species_count) &
           target_species_count > 0) |>
  mutate(intro = screening_wmel_aeg/target_species_count) |>
  group_by(trap_id, latitude, longitude) |>
  summarise(wmean_intro = weighted.mean(intro, w = total_aegypti_caught, na.rm = TRUE),
            mean_intro = mean(intro, na.rm = TRUE)) |>
  st_as_sf(coords = c("longitude", "latitude"), crs = "EPSG:4326 - WGS 84", remove = FALSE) |>
  st_transform(crs = 31983)

intro_1920_gridpre <- st_join(grid, intro_1920, join = st_contains)

intro_1920_grid <- intro_1920_gridpre |>
  st_drop_geometry() |>
  filter(!is.na(trap_id)) |>
  group_by(id) |>
  summarise(mintro = mean(mean_intro, na.rm = TRUE),
            wmintro = mean(wmean_intro, na.rm = TRUE))

dengue_gridpre <- st_join(grid, dengue_raw1, join = st_contains) |>
  filter(year >= 2020)

dengue_gridmod <- dengue_gridpre |>
  st_drop_geometry() |>
  filter(!is.na(year)) |>
  group_by(year, id) |>
  summarise(cases = n()) |>
  ungroup() |>
  right_join(grid |> st_drop_geometry(), by = "id") |>
  tidyr::complete(year, id, fill = list(cases = 0)) |>
  filter(year != 0)

data_mod1 <- grid |>
  right_join(dengue_gridmod, by = "id") |>
  left_join(popdens_grid, by = "id") |>
  left_join(intro_1920_grid, by = "id") |>
  rename(pop2 = sum) |>
  mutate(pop = if_else(pop2 < 50, 0, pop2),
         inc = cases/pop*1000) |>
  filter(st_intersects(geometry, release_areas_j, sparse = FALSE)[,1]) |>
  mutate(geometry = st_centroid(geometry),
         x = st_coordinates(geometry)[,1],
         y = st_coordinates(geometry)[,2]) |>
  st_drop_geometry() |>
  filter(pop != 0 & year >= 2020) |>
  mutate(logpop = log(pop),
         year = factor(year)) |>
  arrange(year) |>
  group_by(year) |>
  mutate(nyear = cur_group_id()) |>
  ungroup() |>
  filter(!is.na(mintro))

spatial.range <- 22000

mesh <- inla.mesh.2d(loc = as.matrix(data_mod1[,c("x","y")]),
                     max.edge = c(spatial.range / 5, spatial.range),
                     cutoff = 200,
                     offset = c(spatial.range / 5, diff(range(data_mod1[ ,c("x")])) / 10))

spde <- inla.spde2.matern(mesh = mesh)

s.modcum <- inla.spde.make.index(name = "spatial.field",
                                 n.spde = spde$n.spde)

A.modcum  = inla.spde.make.A(mesh = mesh,
                             loc = as.matrix(data_mod1[ ,c("x","y")]))

stack.modcum = inla.stack(data = list(obs = data_mod1$cases),
                          A = list(A.modcum, 1),
                          effects = list(spatial_field = data.frame(s.modcum, Intercept = 1),
                                         data_mod1[ ,c(names(data_mod1) != "cases")]),
                          tag = 'stdata')

modcum <- obs ~ -1 + Intercept + mintro + f(spatial.field, model = spde) + offset(logpop)

modcum_out <- inla(
  modcum,
  data = inla.stack.data(stack.modcum),
  family = "nbinomial",
  control.predictor = list(A = inla.stack.A(stack.modcum), compute = TRUE),
  control.fixed = list(expand.factor.strategy = 'inla'),
  control.compute = list(config = TRUE, dic = TRUE, waic = TRUE),
  verbose = TRUE
)

summary(modcum_out)
exp(as.data.frame(inla.zmarginal(modcum_out$marginals.fixed$mintro)))

mintro_marg <- modcum_out$marginals.fixed$mintro

set.seed(123)
beta_samples <- inla.rmarginal(1000, mintro_marg)

mintro_values <- seq(0, 1, length.out = 100)

effect_samples <- sapply(beta_samples, function(b) exp(b * mintro_values))

modcum_df <- data.frame(
  mintro = mintro_values,
  effect = apply(effect_samples, 1, mean),
  lower  = apply(effect_samples, 1, quantile, 0.025),
  upper  = apply(effect_samples, 1, quantile, 0.975)
)

modcum_df |>
  ggplot(aes(x = mintro, y = effect)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_hline(aes(yintercept = 1), lty = 2, color = "red") +
  labs(x = "% wMel introgression in 2019-2020", y = "Ratio of incidences") +
  theme_bw() +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1")) +
  coord_trans(y = "log1p") +
  theme(legend.position = "bottom",
        text = element_text(size = 13))

## spatial field
coords_mod1 <- spde$mesh$loc

xrange_mod1 <- seq(min(coords_mod1[, 1]), max(coords_mod1[, 1]), length.out = 100)
yrange_mod1 <- seq(min(coords_mod1[, 2]), max(coords_mod1[, 2]), length.out = 100)

grid_modcum <- expand.grid(x = xrange_mod1, y = yrange_mod1)
sp::coordinates(grid_modcum) <- ~x + y
proj_grid_modcum <- inla.mesh.projector(spde$mesh, xlim = range(xrange_mod1), ylim = range(yrange_mod1), dims = c(100, 100))

spatial_mean_modcum <- inla.mesh.project(proj_grid_modcum, modcum_out$summary.random$spatial.field$mean)
spatial_lb_modcum <- inla.mesh.project(proj_grid_modcum, modcum_out$summary.random$spatial.field$`0.025quant`)
spatial_ub_modcum <- inla.mesh.project(proj_grid_modcum, modcum_out$summary.random$spatial.field$`0.975quant`)

grid_modcum$mean <- as.vector(spatial_mean_modcum)
grid_modcum$lb <- as.vector(spatial_lb_modcum)
grid_modcum$ub <- as.vector(spatial_ub_modcum)

grid_sf_modcum <- as.data.frame(grid_modcum) |>
  st_as_sf(coords = c("x", "y"), remove = FALSE, crs = 31983) |>
  mutate(sig = case_when(lb > 0 & ub > 0 ~ 1,
                         lb < 0 & ub < 0 ~ 1,
                         .default = 0)) |>
  st_intersection(release_areas_j)

# -- Yearly
data_mod2 <- grid |>
  right_join(dengue_gridmod, by = "id") |>
  left_join(popdens_grid, by = "id") |>
  left_join(intro_1920_grid, by = "id") |>
  rename(pop2 = sum) |>
  mutate(pop = if_else(pop2 < 50, 0, pop2),
         inc = cases/pop*1000) |>
  filter(st_intersects(geometry, release_areas_j, sparse = FALSE)[,1]) |>
  mutate(geometry = st_centroid(geometry),
         x = st_coordinates(geometry)[,1],
         y = st_coordinates(geometry)[,2]) |>
  st_drop_geometry() |>
  filter(pop != 0 & year >= 2020) |>
  mutate(logpop = log(pop),
         year = factor(year)) |>
  arrange(year) |>
  group_by(year) |>
  mutate(nyear = cur_group_id()) |>
  ungroup() |>
  filter(!is.na(mintro))

sum(data_mod2$cases[data_mod2$year == 2020])
sum(data_mod2$cases[data_mod2$year == 2021])
sum(data_mod2$cases[data_mod2$year == 2022])
sum(data_mod2$cases[data_mod2$year == 2023])
sum(data_mod2$cases[data_mod2$year == 2024])

table(data_mod2$year)

spatial.range <- 22000

mesh <- inla.mesh.2d(loc = as.matrix(data_mod2[,c("x","y")]),
                     max.edge = c(spatial.range / 5, spatial.range),
                     cutoff = 200,
                     offset = c(spatial.range / 5, diff(range(data_mod2[ ,c("x")])) / 10))

spde <- inla.spde2.matern(mesh = mesh)

s.index <- inla.spde.make.index(name = "spatial.field",
                                n.spde = spde$n.spde)

A.est = inla.spde.make.A(mesh = mesh,
                         loc = as.matrix(data_mod2[ ,c("x","y")]))

stack = inla.stack(data = list(obs = data_mod2$cases),
                   A = list(A.est, 1),
                   effects = list(spatial_field = data.frame(s.index, Intercept = 1),
                                  data_mod2[ ,c(names(data_mod2) != "cases")]),
                   tag = 'stdata')

mod2 <- obs ~ -1 + Intercept + f(nyear, model = "ar1") + year:mintro + f(spatial.field, model = spde) + offset(logpop)

mod2_out <- inla(
  mod2,
  data = inla.stack.data(stack),
  family = "nbinomial",
  control.predictor = list(A = inla.stack.A(stack), compute = TRUE),
  control.fixed = list(expand.factor.strategy = 'inla'),
  control.compute = list(config = TRUE, dic = TRUE, waic = TRUE),
  verbose = TRUE
)

summary(mod2_out)
exp(as.data.frame(inla.zmarginal(mod2_out$marginals.fixed$`year2020:mintro`)))
exp(as.data.frame(inla.zmarginal(mod2_out$marginals.fixed$`year2021:mintro`)))
exp(as.data.frame(inla.zmarginal(mod2_out$marginals.fixed$`year2022:mintro`)))
exp(as.data.frame(inla.zmarginal(mod2_out$marginals.fixed$`year2023:mintro`)))
exp(as.data.frame(inla.zmarginal(mod2_out$marginals.fixed$`year2024:mintro`)))

mintro_int_margs_mod2 <- list(
  `2020` = mod2_out$marginals.fixed$`year2020:mintro`,
  `2021` = mod2_out$marginals.fixed$`year2021:mintro`,
  `2022` = mod2_out$marginals.fixed$`year2022:mintro`,
  `2023` = mod2_out$marginals.fixed$`year2023:mintro`,
  `2024` = mod2_out$marginals.fixed$`year2024:mintro`
)

set.seed(123)
mod2_df <- bind_rows(lapply(names(mintro_int_margs_mod2), function(yr) {
  beta_year_draws <- inla.rmarginal(1000, mintro_int_margs_mod2[[yr]])
  effect_samples <- sapply(beta_year_draws, function(b) exp(b * mintro_values))

  data.frame(
    year   = as.integer(yr),
    mintro = mintro_values,
    effect = apply(effect_samples, 1, mean),
    lower  = apply(effect_samples, 1, quantile, probs = 0.025),
    upper  = apply(effect_samples, 1, quantile, probs = 0.975)
  )
}))

mod2_df |>
  ggplot(aes(x = mintro, y = effect)) +
  geom_line(lwd = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_hline(aes(yintercept = 1), lty = 2, color = "red") +
  labs(x = "% wMel introgression in 2019-2020", y = "Ratio of incidences", color = "Year", fill = "Year") +
  theme_bw() +
  facet_wrap(~year, ncol = 5) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1")) +
  coord_trans(y = "log1p") +
  theme(legend.position = "bottom",
        text = element_text(size = 15))

## spatial field
coords_mod2 <- spde$mesh$loc

xrange_mod2 <- seq(min(coords_mod2[, 1]), max(coords_mod2[, 1]), length.out = 100)
yrange_mod2 <- seq(min(coords_mod2[, 2]), max(coords_mod2[, 2]), length.out = 100)

grid_mod2 <- expand.grid(x = xrange_mod2, y = yrange_mod2)
sp::coordinates(grid_mod2) <- ~x + y
proj_grid_mod2 <- inla.mesh.projector(spde$mesh, xlim = range(xrange_mod2), ylim = range(yrange_mod2), dims = c(100, 100))

spatial_mean_mod2 <- inla.mesh.project(proj_grid_mod2, mod2_out$summary.random$spatial.field$mean)
spatial_lb_mod2 <- inla.mesh.project(proj_grid_mod2, mod2_out$summary.random$spatial.field$`0.025quant`)
spatial_ub_mod2 <- inla.mesh.project(proj_grid_mod2, mod2_out$summary.random$spatial.field$`0.975quant`)

grid_mod2$mean <- as.vector(spatial_mean_mod2)
grid_mod2$lb <- as.vector(spatial_lb_mod2)
grid_mod2$ub <- as.vector(spatial_ub_mod2)

grid_sf_mod2 <- as.data.frame(grid_mod2) |>
  st_as_sf(coords = c("x", "y"), remove = FALSE, crs = 31983) |>
  mutate(sig = case_when(lb > 0 & ub > 0 ~ 1,
                         lb < 0 & ub < 0 ~ 1,
                         .default = 0)) |>
  st_intersection(release_areas_j)

# -- Model from 2010-2024
dengue_grid_mod_long <- dengue |>
  st_drop_geometry() |>
  filter(!is.na(year) & !is.na(month)) |>
  group_by(year, month, id) |>
  summarise(cases = n()) |>
  ungroup() |>
  right_join(grid |> st_drop_geometry(), by = "id") |>
  tidyr::complete(year, month, id, fill = list(cases = 0)) |>
  filter(!is.na(year) & !is.na(month))

trap_long <- trap_raw |>
  filter(successful == "TRUE" & !is.na(target_species_count) & target_species_count > 0) |>
  mutate(intro = screening_wmel_aeg/target_species_count,
         year = year(collected_at),
         month = month(collected_at)) |>
  group_by(trap_id, type, latitude, longitude, year, month) |>
  summarise(mean_intro = mean(intro, na.rm = TRUE),
            .groups = "drop") |>
  st_as_sf(coords = c("longitude", "latitude"), crs = "EPSG:4326 - WGS 84", remove = FALSE) |>
  st_transform(crs = 31983)

trap_long2 <- st_join(grid, trap_long, join = st_contains)

trap_long_grid <- trap_long2 |>
  st_drop_geometry() |>
  filter(!is.na(trap_id)) |>
  group_by(id, year, month) |>
  summarise(mintro = mean(mean_intro, na.rm = TRUE),
            .groups = "drop")

data_mod_long <- grid |>
  full_join(dengue_grid_mod_long, by = "id") |>
  left_join(popdens_grid, by = "id") |>
  left_join(trap_long_grid, by = c("id", "year", "month")) |>
  left_join(release_grid, by = "id") |>
  rename(pop2 = sum) |>
  mutate(pop = if_else(pop2 < 50, 0, pop2),
         inc = cases/pop*1000,
         date_month = as.Date(paste(year, month, "01", sep = "-")),
         time_since_rel = date_month - month_first_release,
         mintro = if_else(is.na(mintro) & time_since_rel < 0 & year < 2018, 0, mintro),
         mintro = if_else(year < 2017, 0, mintro)) |>
  filter(st_intersects(geometry, release_areas_j, sparse = FALSE)[,1]) |>
  mutate(geometry = st_centroid(geometry),
         x = st_coordinates(geometry)[,1],
         y = st_coordinates(geometry)[,2]) |>
  st_drop_geometry() |>
  filter(pop != 0) |>
  mutate(logpop = log(pop),
         year = factor(year)) |>
  arrange(year) |>
  group_by(year) |>
  mutate(nyear = cur_group_id()) |>
  ungroup() |>
  filter(!is.na(mintro))

table(data_mod_long$year)

spatial.range <- 22000

mesh <- inla.mesh.2d(loc = as.matrix(data_mod_long[,c("x","y")]),
                     max.edge = c(spatial.range / 5, spatial.range),
                     cutoff = 200,
                     offset = c(spatial.range / 5, diff(range(data_mod_long[ ,c("x")])) / 10))

spde <- inla.spde2.matern(mesh = mesh)

s.index <- inla.spde.make.index(name = "spatial.field",
                                n.spde = spde$n.spde)

A.est = inla.spde.make.A(mesh = mesh,
                         loc = as.matrix(data_mod_long[ ,c("x","y")]))

stack = inla.stack(data = list(obs = data_mod_long$cases),
                   A = list(A.est, 1),
                   effects = list(spatial_field = data.frame(s.index, Intercept = 1),
                                  data_mod_long[ ,c(names(data_mod_long) != "cases")]),
                   tag = 'stdata')

modlong1 <- obs ~ -1 + Intercept + f(nyear, model = "ar1") + year:mintro + f(spatial.field, model = spde) + offset(logpop)

modlong1_out <- inla(
  modlong1,
  data = inla.stack.data(stack),
  family = "nbinomial",
  control.predictor = list(A = inla.stack.A(stack), compute = TRUE),
  control.fixed = list(expand.factor.strategy = 'inla'),
  control.compute = list(config = TRUE, dic = TRUE, waic = TRUE),
  verbose = TRUE
)

summary(modlong1_out)

round(exp(as.data.frame(inla.zmarginal(modlong1_out$marginals.fixed$`year2017:mintro`))), 2)
round(exp(as.data.frame(inla.zmarginal(modlong1_out$marginals.fixed$`year2018:mintro`))), 2)
round(exp(as.data.frame(inla.zmarginal(modlong1_out$marginals.fixed$`year2019:mintro`))), 2)
round(exp(as.data.frame(inla.zmarginal(modlong1_out$marginals.fixed$`year2020:mintro`))), 2)
round(exp(as.data.frame(inla.zmarginal(modlong1_out$marginals.fixed$`year2021:mintro`))), 2)
round(exp(as.data.frame(inla.zmarginal(modlong1_out$marginals.fixed$`year2022:mintro`))), 2)
round(exp(as.data.frame(inla.zmarginal(modlong1_out$marginals.fixed$`year2023:mintro`))), 2)
round(exp(as.data.frame(inla.zmarginal(modlong1_out$marginals.fixed$`year2024:mintro`))), 2)

mintro_int_margs_modlong1 <- list(
  `2010` = modlong1_out$marginals.fixed$`year2010:mintro`,
  `2011` = modlong1_out$marginals.fixed$`year2011:mintro`,
  `2012` = modlong1_out$marginals.fixed$`year2012:mintro`,
  `2013` = modlong1_out$marginals.fixed$`year2013:mintro`,
  `2014` = modlong1_out$marginals.fixed$`year2014:mintro`,
  `2015` = modlong1_out$marginals.fixed$`year2015:mintro`,
  `2016` = modlong1_out$marginals.fixed$`year2016:mintro`,
  `2017` = modlong1_out$marginals.fixed$`year2017:mintro`,
  `2018` = modlong1_out$marginals.fixed$`year2018:mintro`,
  `2019` = modlong1_out$marginals.fixed$`year2019:mintro`,
  `2020` = modlong1_out$marginals.fixed$`year2020:mintro`,
  `2021` = modlong1_out$marginals.fixed$`year2021:mintro`,
  `2022` = modlong1_out$marginals.fixed$`year2022:mintro`,
  `2023` = modlong1_out$marginals.fixed$`year2023:mintro`,
  `2024` = modlong1_out$marginals.fixed$`year2024:mintro`
)

mintro_values <- seq(0, 1, length.out = 100)

set.seed(123)
modlong1_df <- bind_rows(lapply(names(mintro_int_margs_modlong1), function(yr) {
  beta_year_draws <- inla.rmarginal(1000, mintro_int_margs_modlong1[[yr]])
  effect_samples <- sapply(beta_year_draws, function(b) exp(b * mintro_values))

  data.frame(
    year   = as.integer(yr),
    mintro = mintro_values,
    effect = apply(effect_samples, 1, mean),
    lower  = apply(effect_samples, 1, quantile, probs = 0.025),
    upper  = apply(effect_samples, 1, quantile, probs = 0.975)
  )
}))

modlong1_df |>
  # filter(year >= 2017) |>
  mutate(effect = if_else(year >= 2017, effect, 1),
         lower = if_else(year >= 2017, lower, -Inf),
         upper = if_else(year >= 2017, upper, Inf)) |>
  ggplot(aes(x = mintro, y = effect)) +
  geom_line(lwd = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_hline(aes(yintercept = 1), lty = 2, color = "red") +
  labs(x = "% wMel introgression", y = "Ratio of incidences", color = "Year", fill = "Year") +
  theme_bw() +
  facet_wrap(~year, ncol = 4) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1")) +
  coord_transform(y = "log1p") +
  theme(legend.position = "bottom",
        text = element_text(size = 15))

modlong1_df |>
  filter(year >= 2017) |>
  ggplot(aes(x = mintro, y = effect)) +
  geom_line(lwd = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_hline(aes(yintercept = 1), lty = 2, color = "red") +
  labs(x = "% wMel introgression", y = "Ratio of incidences", color = "Year", fill = "Year") +
  theme_bw() +
  facet_wrap(~year, ncol = 4) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1")) +
  coord_transform(y = "log1p") +
  theme(legend.position = "bottom",
        text = element_text(size = 15))

mesh <- inla.mesh.2d(loc = as.matrix(data_mod_long[,c("x","y")]),
                     max.edge = c(spatial.range / 5, spatial.range),
                     cutoff = 200,
                     offset = c(spatial.range / 5, diff(range(data_mod_long[ ,c("x")])) / 10))

spde <- inla.spde2.matern(mesh = mesh)

s.index <- inla.spde.make.index(name = "spatial.field",
                                n.spde = spde$n.spde)

A.est = inla.spde.make.A(mesh = mesh,
                         loc = as.matrix(data_mod_long[ ,c("x","y")]))

stack = inla.stack(data = list(obs = data_mod_long$cases),
                   A = list(A.est, 1),
                   effects = list(spatial_field = data.frame(s.index, Intercept = 1),
                                  data_mod_long[ ,c(names(data_mod_long) != "cases")]),
                   tag = 'stdata')

modlong1b <- obs ~ -1 + Intercept + f(nyear, model = "ar1") + mintro + f(spatial.field, model = spde) + offset(logpop)

modlong1b_out <- inla(
  modlong1b,
  data = inla.stack.data(stack),
  family = "nbinomial",
  control.predictor = list(A = inla.stack.A(stack), compute = TRUE),
  control.fixed = list(expand.factor.strategy = 'inla'),
  control.compute = list(config = TRUE, dic = TRUE, waic = TRUE),
  verbose = TRUE
)

summary(modlong1b_out)
exp(as.data.frame(inla.zmarginal(modlong1b_out$marginals.fixed$mintro)))

mintro_marg <- modlong1b_out$marginals.fixed$mintro

set.seed(123)
beta_samples <- inla.rmarginal(1000, mintro_marg)

mintro_values <- seq(0, 1, length.out = 100)

effect_samples <- sapply(beta_samples, function(b) exp(b * mintro_values))

modlong1b_df <- data.frame(
  mintro = mintro_values,
  effect = apply(effect_samples, 1, mean),
  lower  = apply(effect_samples, 1, quantile, 0.025),
  upper  = apply(effect_samples, 1, quantile, 0.975)
)

modlong1b_df |>
  ggplot(aes(x = mintro, y = effect)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_hline(aes(yintercept = 1), lty = 2, color = "red") +
  labs(x = "% wMel introgression", y = "Ratio of incidences") +
  theme_bw() +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1")) +
  coord_trans(y = "log1p") +
  theme(legend.position = "bottom",
        text = element_text(size = 13))

# -- Model from 2017-2024
data_mod_long2 <- data_mod_long |>
  mutate(year = as.numeric(as.character(year))) |>
  filter(year >= 2017) |>
  mutate(year = factor(year)) |>
  arrange(year) |>
  group_by(year) |>
  mutate(nyear = cur_group_id()) |>
  ungroup() |>
  filter(!is.na(mintro))

spatial.range <- 22000

mesh <- inla.mesh.2d(loc = as.matrix(data_mod_long2[,c("x","y")]),
                     max.edge = c(spatial.range / 5, spatial.range),
                     cutoff = 200,
                     offset = c(spatial.range / 5, diff(range(data_mod_long2[ ,c("x")])) / 10))

spde <- inla.spde2.matern(mesh = mesh)

s.index <- inla.spde.make.index(name = "spatial.field",
                                n.spde = spde$n.spde)

A.est = inla.spde.make.A(mesh = mesh,
                         loc = as.matrix(data_mod_long2[ ,c("x","y")]))

stack = inla.stack(data = list(obs = data_mod_long2$cases),
                   A = list(A.est, 1),
                   effects = list(spatial_field = data.frame(s.index, Intercept = 1),
                                  data_mod_long2[ ,c(names(data_mod_long2) != "cases")]),
                   tag = 'stdata')

modlong2 <- obs ~ -1 + Intercept + f(nyear, model = "ar1") + year:mintro + f(spatial.field, model = spde) + offset(logpop)

modlong2_out <- inla(
  modlong2,
  data = inla.stack.data(stack),
  family = "nbinomial",
  control.predictor = list(A = inla.stack.A(stack), compute = TRUE),
  control.fixed = list(expand.factor.strategy = 'inla'),
  control.compute = list(config = TRUE, dic = TRUE, waic = TRUE),
  verbose = TRUE
)

round(exp(as.data.frame(inla.zmarginal(modlong2_out$marginals.fixed$`year2017:mintro`))), 2)
round(exp(as.data.frame(inla.zmarginal(modlong2_out$marginals.fixed$`year2018:mintro`))), 2)
round(exp(as.data.frame(inla.zmarginal(modlong2_out$marginals.fixed$`year2019:mintro`))), 2)
round(exp(as.data.frame(inla.zmarginal(modlong2_out$marginals.fixed$`year2020:mintro`))), 2)
round(exp(as.data.frame(inla.zmarginal(modlong2_out$marginals.fixed$`year2021:mintro`))), 2)
round(exp(as.data.frame(inla.zmarginal(modlong2_out$marginals.fixed$`year2022:mintro`))), 2)
round(exp(as.data.frame(inla.zmarginal(modlong2_out$marginals.fixed$`year2023:mintro`))), 2)
round(exp(as.data.frame(inla.zmarginal(modlong2_out$marginals.fixed$`year2024:mintro`))), 2)

mintro_int_margs_modlong2 <- list(
  `2017` = modlong2_out$marginals.fixed$`year2017:mintro`,
  `2018` = modlong2_out$marginals.fixed$`year2018:mintro`,
  `2019` = modlong2_out$marginals.fixed$`year2019:mintro`,
  `2020` = modlong2_out$marginals.fixed$`year2020:mintro`,
  `2021` = modlong2_out$marginals.fixed$`year2021:mintro`,
  `2022` = modlong2_out$marginals.fixed$`year2022:mintro`,
  `2023` = modlong2_out$marginals.fixed$`year2023:mintro`,
  `2024` = modlong2_out$marginals.fixed$`year2024:mintro`
)

mintro_values <- seq(0, 1, length.out = 100)

set.seed(123)
modlong2_df <- bind_rows(lapply(names(mintro_int_margs_modlong2), function(yr) {
  beta_year_draws <- inla.rmarginal(1000, mintro_int_margs_modlong2[[yr]])
  effect_samples <- sapply(beta_year_draws, function(b) exp(b * mintro_values))

  data.frame(
    year   = as.integer(yr),
    mintro = mintro_values,
    effect = apply(effect_samples, 1, mean),
    lower  = apply(effect_samples, 1, quantile, probs = 0.025),
    upper  = apply(effect_samples, 1, quantile, probs = 0.975)
  )
}))

modlong2_df |>
  ggplot(aes(x = mintro, y = effect)) +
  geom_line(lwd = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_hline(aes(yintercept = 1), lty = 2, color = "red") +
  labs(x = "% wMel introgression", y = "Ratio of incidences", color = "Year", fill = "Year") +
  theme_bw() +
  facet_wrap(~year, ncol = 4) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1")) +
  coord_transform(y = "log1p") +
  theme(legend.position = "bottom",
        text = element_text(size = 15))

spatial.range <- 22000

mesh <- inla.mesh.2d(loc = as.matrix(data_mod_long2[,c("x","y")]),
                     max.edge = c(spatial.range / 5, spatial.range),
                     cutoff = 200,
                     offset = c(spatial.range / 5, diff(range(data_mod_long2[ ,c("x")])) / 10))

spde <- inla.spde2.matern(mesh = mesh)

s.index <- inla.spde.make.index(name = "spatial.field",
                                n.spde = spde$n.spde)

A.est = inla.spde.make.A(mesh = mesh,
                         loc = as.matrix(data_mod_long2[ ,c("x","y")]))

stack = inla.stack(data = list(obs = data_mod_long2$cases),
                   A = list(A.est, 1),
                   effects = list(spatial_field = data.frame(s.index, Intercept = 1),
                                  data_mod_long2[ ,c(names(data_mod_long2) != "cases")]),
                   tag = 'stdata')

modlong2b <- obs ~ -1 + Intercept + f(nyear, model = "ar1") + mintro + f(spatial.field, model = spde) + offset(logpop)

modlong2b_out <- inla(
  modlong2b,
  data = inla.stack.data(stack),
  family = "nbinomial",
  control.predictor = list(A = inla.stack.A(stack), compute = TRUE),
  control.fixed = list(expand.factor.strategy = 'inla'),
  control.compute = list(config = TRUE, dic = TRUE, waic = TRUE),
  verbose = TRUE
)

exp(as.data.frame(inla.zmarginal(modlong2b_out$marginals.fixed$mintro)))

mintro_marg <- modlong2b_out$marginals.fixed$mintro

set.seed(123)
beta_samples <- inla.rmarginal(1000, mintro_marg)

mintro_values <- seq(0, 1, length.out = 100)

effect_samples <- sapply(beta_samples, function(b) exp(b * mintro_values))

modlong2b_df <- data.frame(
  mintro = mintro_values,
  effect = apply(effect_samples, 1, mean),
  lower  = apply(effect_samples, 1, quantile, 0.025),
  upper  = apply(effect_samples, 1, quantile, 0.975)
)

modlong2b_df |>
  ggplot(aes(x = mintro, y = effect)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_hline(aes(yintercept = 1), lty = 2, color = "red") +
  labs(x = "% wMel introgression", y = "Ratio of incidences") +
  theme_bw() +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1")) +
  coord_trans(y = "log1p") +
  theme(legend.position = "bottom",
        text = element_text(size = 13))

# ---- Compare Rio municipality vs release area --------------------------
dengue_releasearea <- st_join(release_areas |>
                            summarise(geometry = st_union(geometry),
                                      .groups = "drop"), dengue_raw, join = st_contains)

pop_releasearea <- exact_extract(popdens, release_areas |>
                                   summarise(geometry = st_union(geometry),
                                             .groups = "drop"), fun = "sum", progress = FALSE)

d_release <- dengue_releasearea |>
  st_drop_geometry() |>
  group_by(year) |>
  summarise(cases = n(),
            .groups = "drop") |>
  mutate(pop = pop_releasearea,
         inc = cases/pop*1000)

rj_joined <- rj_map |>
  filter(name_muni == "Rio de Janeiro") |>
  st_transform(crs = st_crs(dengue_raw)) |>
  dplyr::select(-year)

pop_rio <- exact_extract(popdens, rj_joined, fun = "sum", progress = FALSE)

dengue_riomuni <- st_join(rj_joined, dengue_raw, join = st_contains)

d_riomuni <- dengue_riomuni |>
  st_drop_geometry() |>
  group_by(year) |>
  summarise(cases = n(),
            .groups = "drop") |>
  mutate(pop = pop_rio,
         inc = cases/pop*1000) |>
  left_join(d_release |>
              rename(cases_release = cases,
                     pop_release = pop,
                     inc_release = inc) |>
              dplyr::select(year, cases_release, pop_release, inc_release),
            by = "year") |>
  mutate(cases_norel = cases-cases_release,
         pop_norel = pop-pop_release,
         inc_norel = cases_norel/pop_norel*1000) |>
  pivot_longer(cols = c("inc", "inc_norel", "inc_release"),
               names_to = "type",
               values_to = "incidence") |>
  mutate(type_title = case_when(type == "inc_release" ~ "Release areas",
                                type == "inc" ~ "All of Rio Municipality",
                                type == "inc_norel" ~ "Rest of Rio Municipality"),
         type_title = factor(type_title,
                             levels = c("Release areas", "All of Rio Municipality", "Rest of Rio Municipality")))

# -- Dot plot rio muni vs release area
d_riomuni |>
  filter(type != "inc") |>
  ggplot(aes(x = year)) +
  geom_point(aes(y = incidence, color = type_title), size = 2) +
  scale_color_brewer(palette = "Set1") +
  scale_x_continuous(breaks = seq(2010, 2024, by = 1)) +
  labs(color = NULL, x = NULL, y = "Annual reported dengue incidence per 1,000") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.margin = margin(),
        legend.spacing.y = unit(0, "mm"))

d_riomuni|>
  filter(type != "inc_norel") |>
  ggplot(aes(x = year)) +
  geom_point(aes(y = incidence, color = type_title), size = 2) +
  scale_color_brewer(palette = "Set1") +
  scale_x_continuous(breaks = seq(2010, 2024, by = 1)) +
  labs(color = NULL, x = NULL, y = "Annual reported dengue incidence per 1,000") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.margin = margin(),
        legend.spacing.y = unit(0, "mm"))

pacman::p_load(censobr)
# V0001 is total population by census tract
tracts_df22 <- read_tracts(year = 2022,
                          dataset = "Preliminares",
                          showProgress = FALSE) |>
               filter(name_muni == 'Rio de Janeiro') |>
               collect()

tracts_2022 <- geobr::read_census_tract(code_tract = "RJ",
                                        year = 2022,
                                        simplified = FALSE,
                                        showProgress = FALSE) |>
  filter(name_muni == "Rio de Janeiro")

pop_census <- tracts_2022 |>
  mutate(code_tract = as.character(code_tract)) |>
  dplyr::select(geometry, code_tract) |>
  left_join(tracts_df22, by = "code_tract") |>
  rename(pop = V0001)

subdistrict_shp_state <- read_sf("data/raw/RJ_subdistritos_CD2022.shp")

subdistrict_shp <- subdistrict_shp_state |>
  st_transform(31983) |>
  filter(NM_MUN == "Rio de Janeiro")

pop_subdistrict <- tracts_df22 |>
  filter(!is.na(code_subdistrict)) |>
  group_by(code_subdistrict) |>
  summarise(pop = sum(V0001, na.rm = TRUE),
            .groups = "drop") |>
  left_join(subdistrict_shp, by = c("code_subdistrict" = "CD_SUBDIST")) |>
  st_as_sf()

pop_subdistrict |>
  ggplot() +
  geom_sf(aes(fill = pop)) +
  theme_void()

dengue_rio_sd <- st_join(subdistrict_shp, dengue_raw, join = st_contains)

d_riosd <- dengue_rio_sd |>
  st_drop_geometry() |>
  group_by(year, CD_SUBDIST) |>
  summarise(cases = n(),
            .groups = "drop") |>
  complete(year, CD_SUBDIST, fill = list(cases = 0)) |>
  filter(!is.na(year) & !is.na(CD_SUBDIST)) |>
  left_join(pop_subdistrict |> st_drop_geometry() |> dplyr::select(code_subdistrict, pop),
            by = c("CD_SUBDIST" = "code_subdistrict")) |>
  mutate(inc = cases/pop*1000,
         release_site = if_else(CD_SUBDIST %in% c(33045570515, 33045570516, 33045570525,
                                                  33045570535, 33045570536, 33045570539),
                                "Within release area", "Outside release area")) |>
  left_join(subdistrict_shp |> dplyr::select(CD_SUBDIST, NM_SUBDIST, geometry), by = "CD_SUBDIST") |>
  st_as_sf()

# -- Map of incidence by district
d_riosd |>
  ggplot() +
  geom_sf(aes(fill = inc), color = NA) +
  scale_fill_viridis_c(option = "B") +
  facet_wrap(~year) +
  labs(fill = "Annual reported dengue\nincidence per 1,000") +
  theme_void()

d_riosd |>
  filter(CD_SUBDIST != 33045570526) |>
  ggplot() +
  geom_sf(aes(fill = inc), color = NA) +
  scale_fill_viridis_c(option = "B") +
  facet_wrap(~year) +
  labs(fill = "Annual reported dengue\nincidence per 1,000") +
  theme_void()

pop_subdistrict |>
  ggplot() +
  geom_sf(aes(fill = NM_SUBDIST)) +
  geom_sf(data = release_areas_j, fill = NA, color = "red", lwd = 0.6) +
  theme_void()

pop_subdistrict |>
  filter(NM_SUBDIST %in% c("Penha", "Ilha do Governador", "Maré", "Ramos",
                           "Complexo do Alemão", "Vigário Geral")) |>
  ggplot() +
  geom_sf(aes(fill = NM_SUBDIST)) +
  geom_sf(data = release_areas_j, fill = NA, color = "red", lwd = 0.6) +
  theme_void()

# -- Map of incidence by subdistrict, different legends
years <- sort(unique(d_riosd$year))

popdens_plot <- pop_subdistrict |>
  filter(code_subdistrict != 33045570526) |>
  ggplot() +
  geom_sf(aes(fill = pop), color = NA) +
  geom_sf(data = release_areas_j, fill = NA, color = "red", lwd = 0.6) +
  scale_fill_viridis_c(
    option = "E",
    name = "Population") +
  labs(title = "Population") +
  theme_void() +
  theme(
    plot.title = element_text(size = 10, hjust = 0.5),
    legend.title  = element_text(size = 8),
    legend.text   = element_text(size = 7),
    legend.key.size = unit(0.4, "cm")
    )

plots <- map(years, \(yr) {
  d_riosd |>
    filter(year == yr & CD_SUBDIST != 33045570526) |>
    ggplot() +
    geom_sf(aes(fill = inc), color = NA) +
    geom_sf(data = release_areas_j, fill = NA, color = "red", lwd = 0.6) +
    scale_fill_viridis_c(
      option = "B",
      name   = "Incidence\nper 1000"
    ) +
    labs(title = yr) +
    theme_void() +
    theme(
      plot.title = element_text(size = 10, hjust = 0.5),
      legend.title  = element_text(size = 8),
      legend.text   = element_text(size = 7),
      legend.key.size = unit(0.4, "cm")
    )
})

rio_plot_a <- wrap_plots(c(popdens_plot, plots))

# -- Dotplot of incidence by subdistrict
(rio_plot_b <- d_riosd |>
  filter(release_site == "Outside release area" & CD_SUBDIST != 33045570526) |>
  ggplot(aes(x = year, y = inc, color = release_site)) +
  geom_point(position = position_jitter(width = 0.3), alpha = 0.5, size = 1.5) +
  geom_point(data = d_riosd |> filter(release_site == "Within release area", CD_SUBDIST != 33045570526),
             position = position_jitter(width = 0.3), alpha = 0.7, size = 2) +
  scale_color_brewer(palette = "Set1") +
  labs(color = NULL, x = NULL,
       y = "Dengue incidence\nper 1000") +
  scale_x_continuous(breaks = seq(2010, 2024, by = 1)) +
  theme_bw() +
  theme(legend.position = "inside",
        legend.position.inside = c(0.15, 0.87),
        legend.background = element_rect(color = "black"),
        legend.key = element_rect(fill = NA),
        legend.direction = "horizontal",
        legend.margin = margin(0, 6, 0, 0),
        legend.spacing.y = unit(0, "mm")))

# ggsave("outputs/2026-07-15_compare to rio muni.png", height = 5, width = 7)
# ggsave("outputs/2026-07-15_compare to rio muni.pdf", height = 5, width = 7)

# -- Composite
rio_plot_a / rio_plot_b +
  plot_layout(heights = c(6, 2)) +
  plot_annotation(tag_levels = list(c("A", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "B"))) &
   theme(plot.margin = margin(0, 0, 0, 0))

# ggsave("outputs/2026-07-15_compare to rio muni.pdf", height = 8, width = 12)

d_riosd |>
  filter(release_site == "Outside release area" & CD_SUBDIST != 33045570526) |>
  ggplot(aes(x = year, y = inc)) +
  geom_line(aes(group = NM_SUBDIST, color = release_site), alpha = 0.5) +
  geom_line(data = d_riosd |> filter(release_site == "Within release area", CD_SUBDIST != 33045570526),
             aes(group = NM_SUBDIST, color = NM_SUBDIST), lwd = 1) +
  scale_color_brewer(palette = "Set1") +
  labs(color = NULL, x = NULL,
       y = "Annual reported dengue incidence per 1,000,\n(by subdistricts in Rio de Janeiro municipality)") +
  scale_x_continuous(breaks = seq(2010, 2024, by = 1)) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.margin = margin(),
        legend.spacing.y = unit(0, "mm"))

d_riosd |>
  filter(release_site == "Outside release area" & CD_SUBDIST != 33045570526) |>
  ggplot(aes(x = year, y = log1p(inc))) +
  geom_line(aes(group = NM_SUBDIST, color = release_site), alpha = 0.5) +
  geom_line(data = d_riosd |> filter(release_site == "Within release area", CD_SUBDIST != 33045570526),
             aes(group = NM_SUBDIST, color = NM_SUBDIST), lwd = 1) +
  scale_color_brewer(palette = "Set1") +
  labs(color = NULL, x = NULL,
       y = "Annual reported dengue log-incidence per 1,000,\n(by subdistricts in Rio de Janeiro municipality)") +
  scale_x_continuous(breaks = seq(2010, 2024, by = 1)) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.margin = margin(),
        legend.spacing.y = unit(0, "mm"))
