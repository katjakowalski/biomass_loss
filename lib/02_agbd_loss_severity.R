# ---
# Title:        Estimating severity of AGBD loss
# Description:  Derive AGBD from map and calculate AGBD loss severity, output already generated in repo, additional input needed to run the code: https://zenodo.org/records/8154445
# Author:       Katja Kowalski
# Affiliation:  Technical University of Munich, School of Life Sciences
# Date Created: 2025-08-01
# Last Updated: 2025-08-04
# ---

library(future)
library(future.apply)
library(terra)
library(sf)
library(tidyverse)
library(data.table)
library(terra)


# 01 Derive AGBD time series from Planet AGBD map -------------------------

res <- 20 # set spatial resolution of grid

grid <- read_sf(
  sprintf("data/gis/referencegrid/grid_%skm_3035.shp", res)
)

process_grids <- function(gr, res) {
  grid <- read_sf(
    sprintf("data/gis/referencegrid/grid_%skm_3035.shp", res)
  )

  disturbances <- rast("data/latest_disturbance_eu_v211_3035.tif")

  gridcell <- terra::vect(grid[gr, ])

  gridid <- gridcell$gridid

  planet <- terra::rast("data/planet/planet_agb_30m_v0.1.tif")
  names(planet) <- "planet_agb_2019"

  planet_grid <- terra::crop(planet, gridcell)
  rm(planet)

  planet_grid <- (planet_grid * (10000 / 900)) / 1000

  disturbance_grid <- terra::crop(disturbances, gridcell)
  rm(disturbances)

  # AGBD after tree cover loss
  planet_grid_post <- planet_grid
  planet_grid_post[is.na(planet_grid_post[])] <- 0

  planet_grid_post <- resample(
    planet_grid_post,
    disturbance_grid,
    method = "near"
  )

  # AGBD before tree cover loss
  planet_grid_pre <- planet_grid

  planet_grid_pre <- resample(
    planet_grid_pre,
    disturbance_grid,
    method = "near"
  )

  disturbance_years <- unique(disturbance_grid)
  disturbance_years <- disturbance_years$latest_disturbance_eu_v211_3035

  disturbance_years_pre <- disturbance_years[disturbance_years > 2019]
  disturbance_years_post <- disturbance_years[disturbance_years < 2019]

  df_planet_grid <- data.frame()

  # pre-disturbance
  if (length(disturbance_years_pre) != 0) {
    for (disturbance_year_pre in disturbance_years_pre) {
      dist_mask <- terra::ifel(disturbance_grid != disturbance_year_pre, NA, 1)
      planet_grid_pre_masked <- terra::mask(
        x = planet_grid_pre,
        mask = dist_mask
      )

      # mean per grid cell and year
      df_planet_grid_yr_pre <- global(planet_grid_pre_masked, 'mean', na.rm = T)

      # standard deviation and error per grid cell and year
      planet_agbd_sd = global(planet_grid_pre_masked, 'sd', na.rm = T)
      n_px <- global(planet_grid_pre_masked, fun = "notNA")
      planet_agbd_se = planet_agbd_sd$sd / sqrt(n_px$notNA)

      time_since_disturbance <- 2019 - disturbance_year_pre

      df_planet_grid_yr_pre$sd <- planet_agbd_sd$sd
      df_planet_grid_yr_pre$se <- planet_agbd_se
      df_planet_grid_yr_pre$time_since_disturbance <- time_since_disturbance
      df_planet_grid_yr_pre$nr_dist_px <- n_px$notNA
      df_planet_grid_yr_pre$gridid <- gridid
      rownames(df_planet_grid_yr_pre) <- NULL

      df_planet_grid <- rbind(df_planet_grid, df_planet_grid_yr_pre)
    }
  }

  # post disturbance AGBD

  if (length(disturbance_years_post) != 0) {
    for (disturbance_year_post in disturbance_years_post) {
      dist_mask <- terra::ifel(disturbance_grid != disturbance_year_post, NA, 1)
      planet_grid_post_masked <- terra::mask(
        x = planet_grid_post,
        mask = dist_mask
      )

      # mean per grid cell and year
      df_planet_grid_yr_post <- global(
        planet_grid_post_masked,
        'mean',
        na.rm = T
      )

      # standard deviation and error per grid cell and year
      planet_agbd_sd = global(planet_grid_post_masked, 'sd', na.rm = T)
      n_px <- global(planet_grid_post_masked, fun = "notNA")
      planet_agbd_se = planet_agbd_sd$sd / sqrt(n_px$notNA)

      time_since_disturbance <- 2019 - disturbance_year_post

      df_planet_grid_yr_post$sd <- planet_agbd_sd$sd
      df_planet_grid_yr_post$se <- planet_agbd_se
      df_planet_grid_yr_post$time_since_disturbance <- time_since_disturbance
      df_planet_grid_yr_post$nr_dist_px <- n_px$notNA
      df_planet_grid_yr_post$gridid <- gridid
      rownames(df_planet_grid_yr_post) <- NULL

      df_planet_grid <- rbind(df_planet_grid, df_planet_grid_yr_post)
    }
  }

  if (length(df_planet_grid) != 0) {
    outfile <- sprintf(
      "data/04_aggregation_planet_agbd/agbd_%skm/grid_%s.csv",
      res,
      gridid
    )

    write.csv(df_planet_grid, outfile, row.names = F)
  }
  invisible(gc())
}


num_cores <- 4 # adjust based on available CPU (resource-intensive!)
plan(multisession, workers = num_cores)

grid_ids <- 1:nrow(grid)

future_lapply(grid_ids, function(gr) {
  process_grids(gr, res)
})

future:::ClusterRegistry("stop")

csv_list <- list.files(
  sprintf("data/04_aggregation_planet_agbd/agbd_%skm", res),
  glob2rx("*.csv$"),
  full.names = TRUE
)

df <- plyr::ldply(csv_list, fread)

fwrite(df, sprintf("data/04_aggregation_planet_agbd/planet_agbd_%skm.csv", res))


# 02 Severity of AGBD loss ------------------------------------------------

df_agbd_ts <- fread(sprintf(
  "data/04_aggregation_planet_agbd/planet_agbd_%skm.csv",
  res
))

grid <- read_sf(
  sprintf("data/gis/referencegrid/grid_%skm_3035.shp", res)
)


df_agbd_pre <- df_agbd_ts %>%
  filter(time_since_disturbance < 0) %>%
  rename(planet_agbd = mean) %>%
  group_by(
    gridid
  ) %>%
  summarize(
    planet_agbd_before_mean = mean(planet_agbd, na.rm = TRUE),
    planet_agbd_before_sd = sd(planet_agbd, na.rm = TRUE),
    planet_agbd_before_se = planet_agbd_before_sd / sqrt(n()),
    n_pre = n()
  ) %>%
  ungroup()

df_agbd_post <- df_agbd_ts %>%
  filter(time_since_disturbance > 0 & time_since_disturbance <= 5) %>%
  rename(planet_agbd = mean) %>%
  group_by(
    gridid
  ) %>%
  summarize(
    planet_agbd_after_mean = mean(planet_agbd, na.rm = TRUE),
    planet_agbd_after_sd = sd(planet_agbd, na.rm = TRUE),
    planet_agbd_after_se = planet_agbd_after_sd / sqrt(n()),
    n_post = n()
  ) %>%
  ungroup()


df_agbd <- df_agbd_pre %>%
  full_join(df_agbd_post)


# AGBD loss
df_agbd <- df_agbd %>%
  group_by(gridid) %>%
  mutate(
    gridid = gridid,

    # Planet AGBD
    planet_agbd_before_mean = mean(rnorm(
      10000,
      mean = planet_agbd_before_mean,
      sd = planet_agbd_before_se
    )),
    planet_agbd_after_min = min(rnorm(
      10000,
      mean = planet_agbd_after_mean,
      sd = planet_agbd_after_se
    )),

    planet_agbd_after_q05 = quantile(
      rnorm(10000, mean = planet_agbd_after_mean, sd = planet_agbd_after_se),
      0.05,
      na.rm = T
    ),

    planet_agbd_after_med = median(rnorm(
      10000,
      mean = planet_agbd_after_mean,
      sd = planet_agbd_after_se
    )),

    planet_agbd_after_mean = mean(rnorm(
      10000,
      mean = planet_agbd_after_mean,
      sd = planet_agbd_after_se
    )),

    planet_severity_abs_min = planet_agbd_before_mean - planet_agbd_after_min,
    planet_severity_abs_q05 = planet_agbd_before_mean - planet_agbd_after_q05,
    planet_severity_abs_med = planet_agbd_before_mean - planet_agbd_after_med,
    planet_severity_abs_mean = planet_agbd_before_mean - planet_agbd_after_mean,

    planet_severity_rel_min = planet_severity_abs_min / planet_agbd_before_mean,
    planet_severity_rel_q05 = planet_severity_abs_q05 / planet_agbd_before_mean,
    planet_severity_rel_med = planet_severity_abs_med / planet_agbd_before_mean,
    planet_severity_rel_mean = planet_severity_abs_mean /
      planet_agbd_before_mean
  )


fwrite(
  df_agbd,
  sprintf("data/04_aggregation_planet_agbd/planet_agbd_loss_%skm.csv", res),
  row.names = F
)
