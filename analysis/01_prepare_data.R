# =============================================================================
# 01_prepare_data.R
# -----------------------------------------------------------------------------
# Load and harmonise the raw inputs (SINAN dengue cases, entomological trap
# monitoring, wMel releases, WorldPop population density, WMP release areas),
# build the 0.25 km2 (500 m) analysis grid, and summarise every input to
# cell-months / cell-years.
#
# Output: data/interim/prepared.rds  (a named list of shared objects)
#         data/interim/grid.rds, data/interim/grid.shp
#
# Raster inputs (WorldPop) are NOT serialised into prepared.rds because terra
# rasters cannot be saved with saveRDS(); downstream scripts re-read them from
# data/raw/ where needed. See data/README.md for the required raw files.
# =============================================================================

source(file.path("R", "packages.R"))

# ---- Release areas and population density -----------------------------------

# WMP wMel release areas (restricted to the Rio de Janeiro municipality).
release_areas_raw <- read_sf("data/raw/release-areas.shp") |>
  filter(NOME_DO_MU == "RIO DE JANEIRO")

release_areas <- release_areas_raw |>
  st_zm() |>
  st_transform(crs = 31983) |>
  group_by(REGIAO_2) |>
  summarise(geometry = st_union(geometry)) |>
  mutate(REGIAO_2 = if_else(REGIAO_2 == "Maré (RJ3.3)", "RJ3.3", REGIAO_2))

# Single dissolved release-area polygon (holes removed) used for clipping.
release_areas_j <- release_areas |>
  st_union() |>
  nngeo::st_remove_holes()

# WorldPop 2020 population density, cropped and reprojected to the release area.
popdens <- rast("data/raw/wp_popdens_2020.tif") |>
  crop(release_areas_raw) |>
  project(crs(release_areas))

# ---- Geocoded dengue cases, 2010 - 2024 -------------------------------------
# 2019-2024 come from the Arbogeo extract (filtered to dengue, ICD-10 A90);
# 2010-2018 come from one shapefile per year. All are cleaned identically:
# drop discarded/inconclusive classifications (CLASSI_FIN 5, 8), keep onset date.

dengue_raw1 <- read_sf("data/raw/Arbogeo_2019_2024_04_30.shp") |>
  filter(ID_AGRAVO == "A90" & !is.na(DT_SIN_PRI) & !CLASSI_FIN %in% c(5, 8)) |>
  dplyr::select(date_onset = DT_SIN_PRI) |>
  mutate(year = year(date_onset), month = month(date_onset)) |>
  st_transform(crs = 31983)

read_dengue_year <- function(yr) {
  read_sf(sprintf("data/raw/Dengue%d.shp", yr)) |>
    mutate(date_onset = as.Date(DT_SIN_PRI)) |>
    filter(!CLASSI_FIN %in% c(5, 8)) |>
    dplyr::select(date_onset) |>
    mutate(year = year(date_onset), month = month(date_onset)) |>
    st_transform(crs = 31983)
}

dengue_layers <- lapply(2010:2018, read_dengue_year)
dengue_raw <- do.call(rbind, c(list(dengue_raw1), dengue_layers))

# ---- Entomological trap monitoring ------------------------------------------
trap_raw <- readr::read_delim(
  "data/raw/monitoring_RJO.csv",
  delim = ";", locale = readr::locale(encoding = "ISO-8859-1")
)

# BG-trap wMel prevalence per trap x cell-month (primary exposure source).
trap_year <- trap_raw |>
  filter(type == "BG" & successful == "TRUE" &
           !is.na(target_species_count) & target_species_count > 0) |>
  mutate(intro = screening_wmel_aeg / target_species_count,
         year = year(collected_at),
         month = month(collected_at)) |>
  group_by(trap_id, latitude, longitude, year, month) |>
  summarise(mean_intro = mean(intro, na.rm = TRUE)) |>
  st_as_sf(coords = c("longitude", "latitude"),
           crs = "EPSG:4326 - WGS 84", remove = FALSE) |>
  st_transform(crs = 31983)

# ---- Releases ----------------------------------------------------------------
release_raw <- read_csv("data/raw/new-release-data.csv")

release <- release_raw |>
  group_by(latitude, longitude) |>
  summarise(first_release = min(date_of_release)) |>
  st_as_sf(crs = "EPSG:4326 - WGS 84", coords = c("longitude", "latitude")) |>
  st_transform(crs = 31983)

release_dates <- release_raw |>
  mutate(monthstart = lubridate::floor_date(date_of_release, "month"),
         zone = if_else(str_starts(identifier, "RJ3"),
                        substr(identifier, 1, 6), substr(identifier, 1, 3)),
         zone = case_when(zone == "RJ3.B1" ~ "RJ3.1",
                          zone == "RJ3.V1" ~ "RJ3.1",
                          zone == "RJ3.B2" ~ "RJ3.2",
                          zone == "RJ3.V2" ~ "RJ3.2",
                          zone == "RJ3.B3" ~ "RJ3.3",
                          .default = zone)) |>
  filter(zone != "RJ4") |>
  group_by(monthstart, zone) |>
  summarise(released = sum(average_released))

# ---- Build the 0.25 km2 (500 m) analysis grid -------------------------------
# Keep grid cells at least 50% inside the dissolved release area.
grid_prelim <- st_make_grid(release_areas_j, cellsize = 500) |>
  st_as_sf() |>
  rename(geometry = x) |>
  mutate(id = row_number(), cell_area = st_area(geometry))

intersections <- st_intersection(grid_prelim, release_areas_j) |>
  mutate(intersection_area = st_area(geometry))

grid <- intersections |>
  mutate(perc_intersection = as.numeric(intersection_area / cell_area)) |>
  filter(perc_intersection >= 0.5) |>
  dplyr::select(id, geometry)

if (!dir.exists("data/interim")) dir.create("data/interim", recursive = TRUE)
readr::write_rds(grid, "data/interim/grid.rds")
write_sf(grid, "data/interim/grid.shp", delete_dsn = TRUE)

# ---- Summarise inputs to cell-months ----------------------------------------
# Cases per cell (spatial join of point cases into grid cells).
dengue <- st_join(grid, dengue_raw, join = st_contains)

# Monthly dengue counts per cell (zero-filled).
dengue_grid <- dengue |>
  st_drop_geometry() |>
  filter(!is.na(year)) |>
  group_by(year, month, id) |>
  summarise(cases = n()) |>
  ungroup() |>
  right_join(grid |> st_drop_geometry(), by = "id") |>
  tidyr::complete(year, month, id, fill = list(cases = 0)) |>
  filter(!is.na(year) & !is.na(month))

# Population per cell (zonal sum of WorldPop density).
popdens_grid <- exact_extract(popdens, grid, fun = "sum",
                              append_cols = "id", progress = FALSE)

# BG-trap introgression per cell-month.
trap <- st_join(grid, trap_year, join = st_contains)
trap_grid <- trap |>
  st_drop_geometry() |>
  filter(!is.na(trap_id)) |>
  group_by(id, year, month) |>
  summarise(mintro = mean(mean_intro, na.rm = TRUE))

# First release month per cell.
release_pre <- st_join(grid, release, join = st_contains)
release_grid <- release_pre |>
  st_drop_geometry() |>
  filter(!is.na(first_release)) |>
  group_by(id) |>
  summarise(first_release = min(first_release), .groups = "drop") |>
  mutate(month_first_release = floor_date(first_release, unit = "month")) |>
  dplyr::select(-first_release)

# Combined cell-month table (descriptive; Figure 1).
data <- grid |>
  right_join(dengue_grid, by = "id") |>
  left_join(popdens_grid, by = "id") |>
  left_join(trap_grid, by = c("id", "year", "month")) |>
  left_join(release_grid, by = "id") |>
  filter(!(is.na(month_first_release) & is.na(mintro))) |>
  rename(pop2 = sum) |>
  mutate(pop = if_else(pop2 == 0, NA, pop2),
         inc = cases / pop * 1000,
         month_start = as.Date(paste(year, month, "01", sep = "-")),
         time_since_release = month_start - month_first_release,
         pre_post = factor(if_else(time_since_release < 0, "pre", "post"))) |>
  filter(st_intersects(geometry, release_areas_j, sparse = FALSE)[, 1])

# ---- Reusable modelling building blocks -------------------------------------
# These are consumed by several downstream modelling scripts. Computing them
# here once (a) removes the duplicate definitions that were scattered across the
# original notebook and (b) fixes an ordering bug: intro_1920_grid was used by
# the Figure 3 models but only defined much later in the file.

# Full monthly case grid (all years, zero-filled) used by the primary models.
dengue_grid_mod_long <- dengue |>
  st_drop_geometry() |>
  filter(!is.na(year) & !is.na(month)) |>
  group_by(year, month, id) |>
  summarise(cases = n()) |>
  ungroup() |>
  right_join(grid |> st_drop_geometry(), by = "id") |>
  tidyr::complete(year, month, id, fill = list(cases = 0)) |>
  filter(!is.na(year) & !is.na(month))

# BG + ovitrap introgression per cell-month (both trap types together).
trap_long <- trap_raw |>
  filter(successful == "TRUE" &
           !is.na(target_species_count) & target_species_count > 0) |>
  mutate(intro = screening_wmel_aeg / target_species_count,
         year = year(collected_at), month = month(collected_at)) |>
  group_by(trap_id, type, latitude, longitude, year, month) |>
  summarise(mean_intro = mean(intro, na.rm = TRUE), .groups = "drop") |>
  st_as_sf(coords = c("longitude", "latitude"),
           crs = "EPSG:4326 - WGS 84", remove = FALSE) |>
  st_transform(crs = 31983)

trap_long_grid <- st_join(grid, trap_long, join = st_contains) |>
  st_drop_geometry() |>
  filter(!is.na(trap_id)) |>
  group_by(id, year, month) |>
  summarise(mintro = mean(mean_intro, na.rm = TRUE), .groups = "drop")

# Time-invariant BG-trap introgression from 2019-2020 collections; used as the
# exposure in the Figure 3 "future introgression" models and Model D.
intro_1920 <- trap_raw |>
  filter(type == "BG" & successful == "TRUE" &
           collected_at > as.POSIXct("2018-12-31") &
           !is.na(target_species_count) & target_species_count > 0) |>
  mutate(intro = screening_wmel_aeg / target_species_count) |>
  group_by(trap_id, latitude, longitude) |>
  summarise(wmean_intro = weighted.mean(intro, w = total_aegypti_caught, na.rm = TRUE),
            mean_intro = mean(intro, na.rm = TRUE)) |>
  st_as_sf(coords = c("longitude", "latitude"),
           crs = "EPSG:4326 - WGS 84", remove = FALSE) |>
  st_transform(crs = 31983)

intro_1920_grid <- st_join(grid, intro_1920, join = st_contains) |>
  st_drop_geometry() |>
  filter(!is.na(trap_id)) |>
  group_by(id) |>
  summarise(mintro = mean(mean_intro, na.rm = TRUE),
            wmintro = mean(wmean_intro, na.rm = TRUE))

# ---- Save shared objects ----------------------------------------------------
# Everything except the terra raster (which cannot be serialised with saveRDS).
.shared <- setdiff(
  ls(),
  c("popdens", "read_dengue_year", "dengue_layers", ".shared")
)
saveRDS(mget(.shared), "data/interim/prepared.rds")
message("Wrote data/interim/prepared.rds with ", length(.shared), " objects.")
