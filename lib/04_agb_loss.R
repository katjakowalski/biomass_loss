# ---
# Title:        Total AGB losses from tree cover loss
# Description:  Calculates total AGB losses based on area corrected maps of tree cover loss, output already created in repo
# Author:       Katja Kowalski
# Affiliation:  Technical University of Munich, School of Life Sciences
# Date Created: 2025-08-01
# Last Updated: 2025-08-04
# ---

library(tidyverse)
library(sf)
library(data.table)
library(terra)
library(ggpubr)
library(ggalluvial)
library(lemon)
library(broom)
library(ggrepel)
library(ggplot2)
library(cowplot)
library(ggh4x)
library(colorspace)
library(future)
library(future.apply)


# 01 Map area correction --------------------------------------------------

all_disturbances <- read_csv(
  "data/03_aggregation_area/disturbance_area_20km.csv"
)

forest <- read_csv("data/03_aggregation_area/forest_area_20km.csv")

all_disturbances <- all_disturbances %>%
  left_join(forest) %>%
  mutate(
    dist_ha_sum = dist_px * 0.09,
    forest_ha_sum = forest_px_sum * 0.09
  ) %>%
  rename(
    dist_px_sum = dist_px
  )

# commission errors
errors <- read_csv("data/disturbancemaps_temporal_validation_errors_v2.csv")
errors <- errors[, c(1, 4:6, 9:11)]

colnames(errors) <- c(
  "year",
  "commission_error_undisturbed",
  "se_ce_undisturbed",
  "se_ce_undisturbed_boot",
  "commission_error_disturbed",
  "se_ce_disturbed",
  "se_ce_disturbed_boot"
)

errors <- errors %>%
  mutate(
    commission_error_undisturbed = commission_error_undisturbed / 100,
    se_ce_undisturbed = se_ce_undisturbed / 100,
    se_ce_undisturbed_boot = se_ce_undisturbed_boot / 100,

    commission_error_disturbed = commission_error_disturbed / 100,
    se_ce_disturbed = se_ce_disturbed / 100,
    se_ce_disturbed_boot = se_ce_disturbed_boot / 100
  )


mean_errors <- errors %>%
  filter(year >= 2014) %>%
  summarise_all(mean) %>%
  uncount(length(c(2019:2023))) %>%
  mutate(year = c(2019:2023))

errors <- rbind(errors, mean_errors)

all_disturbances <- dplyr::left_join(all_disturbances, errors, by = "year")

dist_area_orig <- all_disturbances

# write for creating figure
write_csv(
  dist_area_orig,
  "data/03_aggregation_area/disturbance_area_map_orig.csv"
)


area_correction <- function(
  n_sim,
  px_undisturbed,
  px_disturbed,
  com_error_undisturbed,
  se_undisturbed,
  com_error_disturbed,
  se_disturbed
) {
  area_corrected <- c()
  for (i in 1:n_sim) {
    success_rate_0 <- 1 -
      truncnorm::rtruncnorm(1, 0, 1, com_error_undisturbed, se_undisturbed) # commission error and standard error for undisturbed

    success_rate_1 <- 1 -
      truncnorm::rtruncnorm(1, 0, 1, com_error_disturbed, se_disturbed) # comission error and standard error for disturbed

    sim_0 <- px_undisturbed -
      px_disturbed -
      rbinom(1, px_undisturbed - px_disturbed, (success_rate_0)) # px forest - px disturbed

    sim_1 <- rbinom(1, px_disturbed, (success_rate_1)) # nr. px disturbed

    area_corrected <- c(area_corrected, sum(c(sim_1, sim_0)))
  }

  return(area_corrected)
}


process_grids <- function(row_nr) {
  grid_disturbance <- all_disturbances[row_nr, ]

  grid_disturbance <- grid_disturbance %>%
    reframe(
      year = year,
      gridid = gridid,
      forest_px_sum = forest_px_sum,
      forest_ha_sum = forest_ha_sum,
      dist_px_sum = dist_px_sum,
      dist_px_sum_corrected = area_correction(
        n_sim = 100,
        com_error_undisturbed = commission_error_undisturbed,
        se_undisturbed = se_ce_undisturbed_boot,
        com_error_disturbed = commission_error_disturbed,
        se_disturbed = se_ce_disturbed_boot,
        px_undisturbed = as.integer(forest_px_sum),
        px_disturbed = as.integer(dist_px_sum)
      )
    )

  area_correction(
    n_sim = 100,
    com_error_undisturbed = grid_disturbance$commission_error_undisturbed,
    se_undisturbed = grid_disturbance$se_ce_undisturbed_boot,
    com_error_disturbed = grid_disturbance$commission_error_disturbed,
    se_disturbed = grid_disturbance$se_ce_disturbed_boot,
    px_undisturbed = as.integer(grid_disturbance$forest_px_sum),
    px_disturbed = as.integer(grid_disturbance$dist_px_sum)
  )

  out_file <- sprintf("disturbance_area_100_draws_rownr_%s.csv", row_nr)
  write.csv(
    grid_disturbance,
    file.path(
      "data/05_disturbance_area_correction/dist_area_corrected_20km_upd/",
      out_file
    ),
    row.names = F
  )

  invisible(gc())
}

num_cores <- 1 # adjust according to available CPU (resource-intensive!)
plan(multisession, workers = num_cores)

rows <- 1:nrow(all_disturbances)

future_lapply(
  rows,
  function(row_nr) {
    process_grids(row_nr)
  },
  future.seed = TRUE
)

future:::ClusterRegistry("stop")


# 02 Calculate AGB loss ---------------------------------------------------

# Disturbance area

csv_list <- list.files(
  "data/05_disturbance_area_correction/dist_area_corrected_20km",
  glob2rx("*.csv$"),
  full.names = TRUE
)

all_disturbances <- plyr::ldply(csv_list, fread)

dist_area <- all_disturbances %>%
  group_by(gridid, year) %>%
  summarize(
    dist_px_sum_corrected_mean = mean(dist_px_sum_corrected, na.rm = T),
    dist_px_sum_corrected_sd = sd(dist_px_sum_corrected, na.rm = T)
  ) %>%
  mutate(
    dist_ha_mean = (dist_px_sum_corrected_mean * 900) / 10000,
    dist_ha_sd = (dist_px_sum_corrected_sd * 900) / 10000
  )

# agents
agents <- fread("data/03_aggregation_area/disturbance_agents_20km.csv")

agents <- agents %>%
  mutate(
    agent = ifelse(agent == 3, "harvest", "natural") # 1 bb/wind, 2 fire, 3 harvest
  ) %>%
  group_by(gridid, year) %>%
  mutate(
    dist_px_sum_gridcell = sum(dist_px),
    dist_percent = dist_px * 100 / dist_px_sum_gridcell
  ) %>%
  dplyr::select(
    year,
    gridid,
    dist_percent,
    agent
  ) %>%
  mutate(
    dist_percent = ifelse(is.na(dist_percent), 0, dist_percent)
  ) %>%
  group_by(gridid, year, agent) %>%
  summarize(
    dist_percent = sum(dist_percent)
  ) %>%
  pivot_wider(names_from = agent, values_from = dist_percent) %>%
  mutate(
    year = as.numeric(year)
  )

# AGBD loss

df_severity_planet <- read.csv(
  "data/04_aggregation_planet_agbd/planet_agbd_loss_20km.csv"
)

grid <- read_sf(
  "data/gis/referencegrid/grid_20km_3035.shp"
)

# add forest area
forest <- read.csv("data/03_aggregation_area/forest_area_20km.csv")

# filter severity grid by forest share
df_severity_planet <- df_severity_planet %>%
  left_join(forest) %>%
  mutate(
    forest_ha_sum = forest_px_sum * 0.09,
    forest_share = forest_ha_sum / 40000
  ) %>%
  filter(forest_share >= 0.1)

# add country

countries <- read_sf("data/gis/borders/countries_europe_map_3035.shp")

country_intersection <- st_join(
  grid,
  countries,
  join = st_intersects,
  largest = TRUE
)

country_intersection <- country_intersection %>%
  st_drop_geometry() %>%
  group_by(
    gridid
  ) %>%
  summarise(
    country = unique(name0)[1],
    country_iso = unique(isocode)[1]
  )

df_severity_planet <- df_severity_planet %>%
  left_join(country_intersection)


# add ecoregions
ecoregions <- read_sf(
  paste0("data/gis/ecoregions/terrestrial_ecoregions_olson_europe_3035.shp")
)

eco_intersection <- st_join(
  grid,
  ecoregions,
  join = st_intersects,
  largest = TRUE
)

eco_intersection <- eco_intersection %>%
  st_drop_geometry() %>%
  group_by(
    gridid
  ) %>%
  summarise(
    biome = unique(BIOME)[1]
  ) %>%
  mutate(
    biome_name = ifelse(biome == 4, "Temperate", NA),
    biome_name = ifelse(biome == 12, "Mediterranean", biome_name),
    biome_name = ifelse(biome == 8, "Pontic Steppe", biome_name),
    biome_name = ifelse(biome == 5, "Temperate", biome_name),
    biome_name = ifelse(biome == 6, "Boreal", biome_name),
    biome_name = ifelse(biome == 11, "Tundra", biome_name)
  ) %>%
  select(gridid, biome_name)

# filter ecoregions
df_severity_planet <- df_severity_planet %>%
  left_join(eco_intersection) %>%
  filter(biome_name != "Tundra" & biome_name != "Pontic Steppe")

# filter area estimates with grid cells from severity estimates
gridids <- unique(df_severity_planet$gridid)
dist_area <- dist_area[dist_area$gridid %in% gridids, ]


df <- dist_area %>%
  left_join(df_severity_planet, by = c("gridid")) %>%
  mutate(
    # set negative values to 0
    across(starts_with("planet_severity"), ~ ifelse(. < 0, 0, .))
  ) %>%
  group_by(year, gridid) %>%
  mutate(
    planet_bm_loss_min = planet_severity_abs_min * dist_ha_mean,
    planet_bm_loss_q05 = planet_severity_abs_q05 * dist_ha_mean,
    planet_bm_loss_med = planet_severity_abs_med * dist_ha_mean,
    planet_bm_loss_mean = planet_severity_abs_mean * dist_ha_mean,

    planet_bm_loss_min_sd = planet_severity_abs_min * dist_ha_sd,
    planet_bm_loss_q05_sd = planet_severity_abs_q05 * dist_ha_sd,
    planet_bm_loss_med_sd = planet_severity_abs_med * dist_ha_sd,
    planet_bm_loss_mean_sd = planet_severity_abs_mean * dist_ha_sd
  )

fwrite(df, "data/06_agb_loss/agb_loss_planet_20km.csv")

df <- fread("data/06_agb_loss/agb_loss_planet_20km.csv")

df <- df %>%
  dplyr::select(
    gridid,
    year,
    dist_ha_mean,
    dist_ha_sd,
    forest_ha_sum,
    forest_share,
    country,
    country_iso,
    biome_name,
    planet_agbd_before_mean,
    planet_severity_abs_q05,
    planet_severity_abs_mean,
    planet_severity_rel_q05,
    planet_severity_rel_mean,
    planet_bm_loss_q05,
    planet_bm_loss_mean,
    planet_bm_loss_q05_sd,
    planet_bm_loss_mean_sd
  ) %>%
  mutate(
    planet_bm_loss = rowMeans(
      dplyr::select(., planet_bm_loss_q05, planet_bm_loss_mean),
      na.rm = TRUE
    ),
    planet_bm_loss_sd = rowMeans(
      dplyr::select(., planet_bm_loss_q05_sd, planet_bm_loss_mean_sd),
      na.rm = TRUE
    ),

    planet_mean_sev_abs = rowMeans(
      dplyr::select(., planet_severity_abs_q05, planet_severity_abs_mean),
      na.rm = TRUE
    ),
    planet_mean_sev_rel = rowMeans(
      dplyr::select(., planet_severity_rel_q05, planet_severity_rel_mean),
      na.rm = TRUE
    )
  )


df_full <- df

# remove overestimation in Sweden and Finland
df <- df %>%
  filter(!(biome_name == "Boreal" & (year == 2018 | year == 2023)))


fwrite(df, "data/06_agb_loss/agb_loss_planet_20km_filtered.csv")
fwrite(df_full, "data/06_agb_loss/agb_loss_planet_20km.csv")
