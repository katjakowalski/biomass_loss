# ---
# Title:        Area extraction for tree cover loss, forest land use, disturbance agents
# Description:  Prepare most important area datasets, output already created in repo, additional input needed to run the code: https://zenodo.org/records/13333034
# Author:       Katja Kowalski
# Affiliation:  Technical University of Munich, School of Life Sciences
# Date Created: 2025-08-01
# Last Updated: 2025-08-04
# ---

library(gdalUtilities)
library(future)
library(future.apply)
library(sf)
library(terra)
library(tidyverse)
library(data.table)


# 01 Reproject disturbance area ------------------------------------------

reproject_years <- function(yr) {
  in_file = sprintf(
    "data/mosaics_annual_disturb_v211/%s_disturb_mosaic_v211_22.tif",
    yr
  )

  ras_vrt = sprintf(
    "data/mosaics_annual_disturb_v211/%s_disturb_mosaic_v211_22_epsg3035.vrt",
    yr
  )

  out_file = sprintf(
    "data/mosaics_annual_disturb_v211/%s_disturb_mosaic_v211_22_epsg3035.tif",
    yr
  )

  gdalwarp(
    srcfile = in_file,
    dstfile = ras_vrt,
    r = "near",
    ot = "Byte",
    t_srs = "EPSG:3035",
    srcnodata = 0,
    dstnodata = 0,
    overwrite = TRUE,
    tr = c(30, 30)
  )

  gdal_translate(
    src_dataset = ras_vrt,
    dst_dataset = out_file,
    ot = "Byte",
    co = "COMPRESS=LZW"
  )

  invisible(gc())
}

num_cores <- 1

plan(multisession, workers = num_cores)

years <- c(1995:2022)

future_lapply(years, function(yr) {
  reproject_years(yr)
})

future:::ClusterRegistry("stop")


# 02 Reproject agents ----------------------------------------------------

reproject_years <- function(yr) {
  in_file = sprintf(
    "data/02_mosaics_agents/v211/%s_disturb_agent_v211_reclass_compv211.tif",
    yr
  )

  ras_vrt = sprintf(
    "data/02_mosaics_agents/v211/%s_disturb_agent_v211_reclass_compv211_epsg3035.vrt",
    yr
  )

  out_file = sprintf(
    "data/02_mosaics_agents/v211/%s_disturb_agent_v211_reclass_compv211_epsg3035.tif",
    yr
  )

  gdalwarp(
    srcfile = in_file,
    dstfile = ras_vrt,
    r = "near",
    ot = "Byte",
    t_srs = "EPSG:3035",
    srcnodata = 0,
    dstnodata = 0,
    overwrite = TRUE,
    tr = c(30, 30)
  )

  gdal_translate(
    src_dataset = ras_vrt,
    dst_dataset = out_file,
    ot = "Byte",
    co = "COMPRESS=LZW"
  )

  invisible(gc())
}

num_cores <- 1

plan(multisession, workers = num_cores)

years <- c(1985:2023)

future_lapply(years, function(yr) {
  reproject_years(yr)
})

future:::ClusterRegistry("stop")


# 03 Create regular grid -------------------------------------------------

countries <- read_sf(
  "data/gis/borders/europe_outline_3035.shp"
)

for (res in c(10, 20, 30, 40, 50)) {
  grid <- st_make_grid(
    x = st_as_sfc(round(st_bbox(countries), 0), crs = st_crs(countries)),
    cellsize = c(res * 1000, res * 1000)
  )

  grid <- st_as_sf(grid)
  grid$gridid <- 1:nrow(grid)
  grid <- grid[countries, ]

  st_write(
    grid,
    paste0("data/gis/referencegrid/grid_", res, "km_3035.shp"),
    append = FALSE
  )
}


# 04 Extract disturbance area per grid cell ----------------------

res <- 20
grid <- st_read(
  sprintf("data/gis/referencegrid/grid_%skm_3035.shp", res)
)

grid_ids <- unique(grid$gridid)

process_grids <- function(gr) {
  res <- 20
  grid <- read_sf(
    paste0("data/gis/referencegrid/grid_", res, "km_3035.shp")
  )

  distdat <- list.files(
    "data/01_mosaics_annual_disturb_v211/",
    pattern = "epsg3035.tif$",
    full.names = TRUE
  )

  disturbances <- terra::rast(distdat)

  gridcell <- terra::vect(grid[gr, ])

  gridid <- gridcell$gridid

  disturbance_gridcell <- terra::crop(disturbances, gridcell)

  dist_px <- as.data.frame(global(disturbance_gridcell, "sum", na.rm = TRUE))
  colnames(dist_px) <- "dist_px"
  dist_px$layer_name <- rownames(dist_px)
  dist_px$year <- str_split_i(dist_px$layer_name, "_", 1)
  dist_px$gridid <- gridid

  outfile <- sprintf(
    "data/03_aggregation_area/disturbance_area_%skm/grid_%s.csv",
    res,
    gridid
  )

  write.csv(dist_px, outfile, row.names = F)

  invisible(gc())
}


num_cores <- 15
plan(multisession, workers = num_cores)

grid_ids <- 1:nrow(grid)

future_lapply(grid_ids, function(gr) {
  process_grids(gr)
})

future:::ClusterRegistry("stop")


res <- 10
dist_area_cells <- list.files(
  sprintf("data/03_aggregation_area/disturbance_agent_%skm", res),
  pattern = ".csv$",
  full.names = TRUE
)

dist_area_cells <- do.call(rbind, lapply(dist_area_cells, fread))


write.csv(
  dist_area_cells,
  sprintf("data/03_aggregation_area/disturbance_area_%skm.csv", res)
)


# 05 Extract forest land use per grid cell --------------------------

res <- 20
grid <- read_sf(paste0("data/gis/referencegrid/grid_", res, "km_3035.shp"))

forest <- rast("data/forestlanduse_mask_EUmosaic_3035.tif")

forest_extract <- exactextractr::exact_extract(
  x = forest,
  y = grid,
  fun = function(values, coverage_fraction, ...) {
    sum(values[coverage_fraction > 0.5], na.rm = T)
  },
  max_cells_in_memory = 3e+10
)

dat_forest <- data.frame(gridid = grid$gridid, forest_px_sum = forest_extract)

write_csv(
  dat_forest,
  sprintf("data/03_aggregation_area/forest_area_%skm.csv", res)
)

# 06 Extract disturbance agents pixel counts per grid cell ---------------

res <- 20
grid <- st_read(sprintf("data/gis/referencegrid/grid_%skm_3035.shp", res))

grid_ids <- unique(grid$gridid)


process_grids <- function(gr) {
  res <- 30
  grid <- read_sf(paste0("data/gis/referencegrid/grid_", res, "km_3035.shp"))

  distdat <- list.files(
    "data/02_mosaics_agents/v211",
    pattern = "*epsg3035.tif$",
    full.names = TRUE
  )

  disturbances <- terra::rast(distdat)

  gridcell <- terra::vect(grid[gr, ])

  gridid <- gridcell$gridid

  disturbance_gridcell <- terra::crop(disturbances, gridcell) #ext(gridcell))

  dist_px <- data.frame()

  for (agent in 1:3) {
    # 1 = bark-beetle/wind, 2 = fire, 3 = harvest

    disturbance_gridcell_agent <- ifel(disturbance_gridcell == agent, 1, 0)

    dist_px_agent <- as.data.frame(global(
      disturbance_gridcell_agent,
      "sum",
      na.rm = TRUE
    ))
    colnames(dist_px_agent) <- "dist_px"
    dist_px_agent$layer_name <- rownames(dist_px_agent)
    dist_px_agent$year <- str_split_i(dist_px_agent$layer_name, "_", 1)
    dist_px_agent$gridid <- gridid
    dist_px_agent$agent <- agent

    dist_px <- rbind(dist_px, dist_px_agent)
  }

  outfile <- sprintf(
    "data/03_aggregation_area/disturbance_agent_%skm/grid_%s.csv",
    res,
    gridid
  )

  write.csv(dist_px, outfile, row.names = F)

  invisible(gc())
}


num_cores <- 35
plan(multisession, workers = num_cores)

grid_ids <- 1:nrow(grid)

future_lapply(grid_ids, function(gr) {
  process_grids(gr)
})

future:::ClusterRegistry("stop")


dist_area_agent_cells <- list.files(
  sprintf("data/03_aggregation_area/disturbance_agent_%skm", res),
  pattern = ".csv$",
  full.names = TRUE
)

dist_area_agent_cells <- do.call(rbind, lapply(dist_area_agent_cells, fread))


write.csv(
  dist_area_agent_cells,
  sprintf("data/03_aggregation_area/disturbance_agents_%skm.csv", res)
)
