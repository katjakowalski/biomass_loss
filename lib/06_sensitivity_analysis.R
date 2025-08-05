# ---
# Title:        Analysis on sensitivity of AGB loss to area of tree cover loss
# Description:  Creates models per ecoregion and per agent, figures 3, S10-S16, and table S3. Output already created in repo.
# Author:       Katja Kowalski
# Affiliation:  Technical University of Munich, School of Life Sciences
# Date Created: 2025-08-01
# Last Updated: 2025-08-04
# ---

library(tidyverse)
library(sdmTMB)
library(spdep)
library(sf)
library(ggplot2)
library(glmmTMB)
library(data.table)
library(RhpcBLASctl)
library(sp)
library(ggpubr)
library(TMB)
library(segmented)

# 01 Prepare data ---------------------------------------------------------

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


# Disturbance agents
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


df <- df %>%
  left_join(agents, by = c("gridid", "year")) %>%
  mutate(
    # set negative values to 0
    across(starts_with("planet_severity"), ~ ifelse(. < 0, 0, .))
  ) %>%
  group_by(year, gridid) %>%
  mutate(
    bm_loss_harvest = planet_bm_loss * harvest / 100,
    bm_loss_harvest_sd = planet_bm_loss_sd * harvest / 100,

    bm_loss_natural = planet_bm_loss * natural / 100,
    bm_loss_natural_sd = planet_bm_loss_sd * natural / 100,

    bm_loss_harvest_share = bm_loss_harvest / planet_bm_loss,
    bm_loss_natural_share = bm_loss_natural / planet_bm_loss,

    dist_area_harvest = dist_ha_mean * harvest / 100,
    dist_area_natural = dist_ha_mean * natural / 100
  )


df_full <- df

# remove overestimation in Sweden and Finland
df <- df %>%
  filter(!(biome_name == "Boreal" & (year == 2018 | year == 2023)))


grid <- read_sf(
  "data/gis/referencegrid/grid_20km_3035.shp"
)

df_mod <- grid %>%
  st_centroid() %>%
  mutate(X = st_coordinates(.)[, 1], Y = st_coordinates(.)[, 2]) %>%
  left_join(df) %>%
  mutate(
    bm_loss_natural = bm_loss_natural / forest_ha_sum,
    bm_loss_harvest = bm_loss_harvest / forest_ha_sum,

    dist_rate_natural = dist_area_natural / forest_ha_sum,
    dist_rate_harvest = dist_area_harvest / forest_ha_sum
  ) %>%
  drop_na(
    year,
    bm_loss_natural,
    bm_loss_harvest,
    dist_rate_natural,
    dist_rate_harvest
  ) %>%
  sf::st_drop_geometry() %>%
  mutate(
    x_scaled = X / 1000,
    y_scaled = Y / 1000,
    year_num = as.integer(year),
    year_fct = factor(year)
  )


# 02 Models temperate -----------------------------------------------------

blas_set_num_threads(10) # set according to available CPUs (resource-intensive)
blas_get_num_procs() # verify

# natural
df_temperate_natural <- df_mod %>%
  filter(biome_name == "Temperate") %>%
  filter(bm_loss_natural > 0)


mesh_temperate <- make_mesh(
  df_temperate_natural,
  xy_cols = c("x_scaled", "y_scaled"),
  cutoff = 40
)

fit_natural_temperate <- sdmTMB(
  bm_loss_natural ~ 0 + dist_rate_natural + (0 + dist_rate_natural | year_fct),
  data = df_temperate_natural,
  mesh = mesh_temperate,
  family = gaussian(),
  spatial = "on",
  spatiotemporal = "ar1",
  time = "year_num",
  silent = FALSE
)

sanity(fit_natural_temperate)

effects_natural_temperate <- data.frame(ranef(fit_natural_temperate))
effects_natural_temperate$year <- as.numeric(rownames(
  effects_natural_temperate
))

fix_natural_temperate <- as.numeric(fixef(fit_natural_temperate))
effects_natural_temperate$slope <- fix_natural_temperate +
  effects_natural_temperate$dist_rate_natural


# harvest
df_temperate_harvest <- df_temperate %>%
  filter(bm_loss_harvest > 0)

mesh_temperate <- make_mesh(
  df_temperate_harvest,
  xy_cols = c("x_scaled", "y_scaled"),
  cutoff = 60
)

fit_harvest_temperate <- sdmTMB(
  bm_loss_harvest ~ 0 + dist_rate_harvest + (0 + dist_rate_harvest | year_fct),
  data = df_temperate_harvest,
  mesh = mesh_temperate,
  family = gaussian(),
  spatial = "on",
  spatiotemporal = "iid",
  time = "year_num",
  silent = FALSE
)

sanity(fit_harvest_temperate)

effects_harvest_temperate <- data.frame(ranef(fit_harvest_temperate))
effects_harvest_temperate$year <- as.numeric(rownames(
  effects_harvest_temperate
))

fix_harvest_temperate <- as.numeric(fixef(fit_harvest_temperate))
effects_harvest_temperate$slope <- fix_harvest_temperate +
  effects_harvest_temperate$dist_rate_harvest


# 03 Models boreal --------------------------------------------------------

# natural
df_boreal_natural <- df_mod %>%
  filter(biome_name == "Boreal") %>%
  filter(bm_loss_natural > 0)

mesh_boreal <- make_mesh(
  df_boreal_natural,
  xy_cols = c("x_scaled", "y_scaled"),
  cutoff = 40
)

fit_natural_boreal <- sdmTMB(
  bm_loss_natural ~ 0 + dist_rate_natural + (0 + dist_rate_natural | year_fct),
  data = df_boreal_natural,
  mesh = mesh_boreal,
  family = gaussian(),
  spatial = "on",
  spatiotemporal = "ar1",
  time = "year_num",
  extra_time = c(2018, 2023),
  silent = FALSE
)

sanity(fit_natural_boreal)

effects_natural_boreal <- data.frame(ranef(fit_natural_boreal))
effects_natural_boreal$year <- as.numeric(rownames(effects_natural_boreal))

fix_natural_boreal <- as.numeric(fixef(fit_natural_boreal))
effects_natural_boreal$slope <- fix_natural_boreal +
  effects_natural_boreal$dist_rate_natural


# harvest
df_boreal_harvest <- df_mod %>%
  filter(biome_name == "Boreal") %>%
  filter(bm_loss_harvest > 0)

mesh_boreal <- make_mesh(
  df_boreal_harvest,
  xy_cols = c("x_scaled", "y_scaled"),
  cutoff = 40
)

fit_harvest_boreal <- sdmTMB(
  bm_loss_harvest ~ 0 + dist_rate_harvest + (0 + dist_rate_harvest | year_fct),
  data = df_boreal_harvest,
  mesh = mesh_boreal,
  family = gaussian(),
  spatial = "on",
  spatiotemporal = "iid",
  time = "year_num",
  extra_time = c(2018, 2023),
  silent = FALSE
)

sanity(fit_harvest_boreal)

effects_harvest_boreal <- data.frame(ranef(fit_harvest_boreal))
effects_harvest_boreal$year <- as.numeric(rownames(effects_harvest_boreal))

fix_harvest_boreal <- as.numeric(fixef(fit_harvest_boreal))
effects_harvest_boreal$slope <- fix_harvest_boreal +
  effects_harvest_boreal$dist_rate_harvest


# 04 Models Mediterranean -------------------------------------------------

# natural
df_med_natural <- df_mod %>%
  filter(biome_name == "Mediterranean") %>%
  filter(bm_loss_natural > 0)

mesh_med <- make_mesh(
  df_med_natural,
  xy_cols = c("x_scaled", "y_scaled"),
  cutoff = 40
)

fit_natural_med <- sdmTMB(
  bm_loss_natural ~ 0 + dist_rate_natural + (0 + dist_rate_natural | year_fct),
  data = df_med_natural,
  mesh = mesh_med,
  family = gaussian(),
  spatial = "on",
  spatiotemporal = "ar1",
  time = "year_num",
  silent = FALSE
)

sanity(fit_natural_med)

effects_natural_med <- data.frame(ranef(fit_natural_med))
effects_natural_med$year <- as.numeric(rownames(effects_natural_med))

fix_natural_med <- as.numeric(fixef(fit_natural_med))
effects_natural_med$slope <- fix_natural_med +
  effects_natural_med$dist_rate_natural

# harvest
df_med_harvest <- df_mod %>%
  filter(biome_name == "Mediterranean") %>%
  filter(bm_loss_harvest > 0)

mesh_med <- make_mesh(
  df_med_harvest,
  xy_cols = c("x_scaled", "y_scaled"),
  cutoff = 40
)

fit_harvest_med <- sdmTMB(
  bm_loss_harvest ~ 0 + dist_rate_harvest + (0 + dist_rate_harvest | year_fct),
  data = df_med_harvest,
  mesh = mesh_med,
  family = gaussian(),
  spatial = "on",
  spatiotemporal = "iid",
  time = "year_num",
  silent = FALSE
)

sanity(fit_harvest_med)

effects_harvest_med <- data.frame(ranef(fit_harvest_med))
effects_harvest_med$year <- as.numeric(rownames(effects_harvest_med))

fix_harvest_med <- as.numeric(fixef(fit_harvest_med))
effects_harvest_med$slope <- fix_harvest_med +
  effects_harvest_med$dist_rate_harvest


# 05 Extract sensitivity (random effects) ---------------------------------

# Temperate Harvest

sdr_harvest_temperate <- TMB::sdreport(fit_harvest_temperate$tmb_obj)
re_harvest_temperate <- as.data.frame(summary(sdr_harvest_temperate, "random"))
re_harvest_temperate$type <- rownames(re_harvest_temperate)

re_harvest_temperate <- re_harvest_temperate %>%
  filter(grepl('re_b_pars', type)) %>%
  mutate(
    idx = if_else(
      type == "re_b_pars",
      0L,
      as.integer(sub("re_b_pars\\.", "", type))
    ),
    year = c(1985:2023)
  )

effect_distrate_harvest_temperate <- tidy(fit_harvest_temperate)$estimate
se_distrate_harvest_temperate <- tidy(fit_harvest_temperate)$std.error

re_harvest_temperate <- re_harvest_temperate %>%
  mutate(
    slope_year = effect_distrate_harvest_temperate + Estimate,
    combined_se = sqrt(se_distrate_harvest_temperate^2 + `Std. Error`^2),
    lower_95 = slope_year - 1.96 * combined_se,
    upper_95 = slope_year + 1.96 * combined_se
  ) %>%
  mutate(
    biome_name = "Temperate",
    agent = "Harvest"
  ) %>%
  dplyr::select(
    year,
    slope_year,
    combined_se,
    lower_95,
    upper_95,
    biome_name,
    agent
  )


# Temperate Natural

sdr_natural_temperate <- TMB::sdreport(fit_natural_temperate$tmb_obj)
re_natural_temperate <- as.data.frame(summary(sdr_natural_temperate, "random"))
re_natural_temperate$type <- rownames(re_natural_temperate)

re_natural_temperate <- re_natural_temperate %>%
  filter(grepl('re_b_pars', type)) %>%
  mutate(
    idx = if_else(
      type == "re_b_pars",
      0L,
      as.integer(sub("re_b_pars\\.", "", type))
    ),
    year = c(1985:2023)
  )

effect_distrate_natural_temperate <- tidy(fit_natural_temperate)$estimate
se_distrate_natural_temperate <- tidy(fit_natural_temperate)$std.error

re_natural_temperate <- re_natural_temperate %>%
  mutate(
    slope_year = effect_distrate_natural_temperate + Estimate,
    combined_se = sqrt(se_distrate_natural_temperate^2 + `Std. Error`^2),
    lower_95 = slope_year - 1.96 * combined_se,
    upper_95 = slope_year + 1.96 * combined_se
  ) %>%
  mutate(
    biome_name = "Temperate",
    agent = "Natural"
  ) %>%
  dplyr::select(
    year,
    slope_year,
    combined_se,
    lower_95,
    upper_95,
    biome_name,
    agent
  )


# Boreal Harvest

sdr_harvest_boreal <- TMB::sdreport(fit_harvest_boreal$tmb_obj)
re_harvest_boreal <- as.data.frame(summary(sdr_harvest_boreal, "random"))
re_harvest_boreal$type <- rownames(re_harvest_boreal)

re_harvest_boreal <- re_harvest_boreal %>%
  filter(grepl('re_b_pars', type)) %>%
  mutate(
    idx = if_else(
      type == "re_b_pars",
      0L,
      as.integer(sub("re_b_pars\\.", "", type))
    ),
    year = c(1985:2017, 2019:2022)
  )

effect_distrate_harvest_boreal <- tidy(fit_harvest_boreal)$estimate
se_distrate_harvest_boreal <- tidy(fit_harvest_boreal)$std.error

re_harvest_boreal <- re_harvest_boreal %>%
  mutate(
    slope_year = effect_distrate_harvest_boreal + Estimate,
    combined_se = sqrt(se_distrate_harvest_boreal^2 + `Std. Error`^2),
    lower_95 = slope_year - 1.96 * combined_se,
    upper_95 = slope_year + 1.96 * combined_se
  ) %>%
  mutate(
    biome_name = "Boreal",
    agent = "Harvest"
  ) %>%
  dplyr::select(
    year,
    slope_year,
    combined_se,
    lower_95,
    upper_95,
    biome_name,
    agent
  )


# Boreal Natural

sdr_natural_boreal <- TMB::sdreport(fit_natural_boreal$tmb_obj)
re_natural_boreal <- as.data.frame(summary(sdr_natural_boreal, "random"))
re_natural_boreal$type <- rownames(re_natural_boreal)

re_natural_boreal <- re_natural_boreal %>%
  filter(grepl('re_b_pars', type)) %>%
  mutate(
    idx = if_else(
      type == "re_b_pars",
      0L,
      as.integer(sub("re_b_pars\\.", "", type))
    ),
    year = c(1985:2017, 2019:2022)
  )

effect_distrate_natural_boreal <- tidy(fit_natural_boreal)$estimate
se_distrate_natural_boreal <- tidy(fit_natural_boreal)$std.error

re_natural_boreal <- re_natural_boreal %>%
  mutate(
    slope_year = effect_distrate_natural_boreal + Estimate,
    combined_se = sqrt(se_distrate_natural_boreal^2 + `Std. Error`^2),
    lower_95 = slope_year - 1.96 * combined_se,
    upper_95 = slope_year + 1.96 * combined_se
  ) %>%
  mutate(
    biome_name = "Boreal",
    agent = "Natural"
  ) %>%
  dplyr::select(
    year,
    slope_year,
    combined_se,
    lower_95,
    upper_95,
    biome_name,
    agent
  )


# Mediterranean Harvest

sdr_harvest_med <- TMB::sdreport(fit_harvest_med$tmb_obj)
re_harvest_med <- as.data.frame(summary(sdr_harvest_med, "random"))
re_harvest_med$type <- rownames(re_harvest_med)

re_harvest_med <- re_harvest_med %>%
  filter(grepl('re_b_pars', type)) %>%
  mutate(
    idx = if_else(
      type == "re_b_pars",
      0L,
      as.integer(sub("re_b_pars\\.", "", type))
    ),
    year = c(1985:2023)
  )

effect_distrate_harvest_med <- tidy(fit_harvest_med)$estimate
se_distrate_harvest_med <- tidy(fit_harvest_med)$std.error

re_harvest_med <- re_harvest_med %>%
  mutate(
    slope_year = effect_distrate_harvest_med + Estimate,
    combined_se = sqrt(se_distrate_harvest_med^2 + `Std. Error`^2),
    lower_95 = slope_year - 1.96 * combined_se,
    upper_95 = slope_year + 1.96 * combined_se
  ) %>%
  mutate(
    biome_name = "Mediterranean",
    agent = "Harvest"
  ) %>%
  dplyr::select(
    year,
    slope_year,
    combined_se,
    lower_95,
    upper_95,
    biome_name,
    agent
  )

# Mediterranean Natural

sdr_natural_med <- TMB::sdreport(fit_natural_med$tmb_obj)
re_natural_med <- as.data.frame(summary(sdr_natural_med, "random"))
re_natural_med$type <- rownames(re_natural_med)

re_natural_med <- re_natural_med %>%
  filter(grepl('re_b_pars', type)) %>%
  mutate(
    idx = if_else(
      type == "re_b_pars",
      0L,
      as.integer(sub("re_b_pars\\.", "", type))
    ),
    year = c(1985:2023)
  )

effect_distrate_natural_med <- tidy(fit_natural_med)$estimate
se_distrate_natural_med <- tidy(fit_natural_med)$std.error

re_natural_med <- re_natural_med %>%
  mutate(
    slope_year = effect_distrate_natural_med + Estimate,
    combined_se = sqrt(se_distrate_natural_med^2 + `Std. Error`^2),
    lower_95 = slope_year - 1.96 * combined_se,
    upper_95 = slope_year + 1.96 * combined_se
  ) %>%
  mutate(
    biome_name = "Mediterranean",
    agent = "Natural"
  ) %>%
  dplyr::select(
    year,
    slope_year,
    combined_se,
    lower_95,
    upper_95,
    biome_name,
    agent
  )


# combine random effects from all models

re <- bind_rows(
  re_harvest_boreal,
  re_harvest_temperate,
  re_harvest_med,
  re_natural_boreal,
  re_natural_temperate,
  re_natural_med
)

write_csv(re, "data/07_results/models/annual_slope_biome_agent_sdmtmb_40km.csv")


# 06 Segmentation of sensitivity time series ------------------------------

re <- read_csv(
  "data/07_results/models/annual_slope_biome_agent_sdmtmb_40km.csv"
)

df_re <- re %>%
  ungroup() %>%
  mutate(
    year_nr = year - 1984
  )


# Test linear vs. segmented model

biomes <- c("Temperate", "Boreal", "Mediterranean")
agents <- c("Natural", "Harvest")


df_res <- data.frame()

for (biome in biomes) {
  for (agent_sel in agents) {
    df_biome <- df_re %>%
      filter(biome_name == biome & agent == agent_sel) %>%
      dplyr::select(
        year,
        slope_year,
        lower_95,
        upper_95,
        combined_se,
        biome_name,
        agent,
        year_nr
      )

    lm_fit <- lm(slope_year ~ year_nr, data = df_biome)

    lm_aic <- AIC(lm_fit)

    seg_fit <- segmented(lm_fit, seg.Z = ~year_nr, psi = 20)

    seg_aic <- AIC(seg_fit)

    if ((lm_aic - seg_aic) > 2) {
      # breakpoints
      bp_summary <- summary(seg_fit)$psi
      bp <- bp_summary[1, "Est."]
      bp_se <- bp_summary[1, "St.Err"]
      bp_ci_low <- bp - 1.96 * bp_se
      bp_ci_high <- bp + 1.96 * bp_se

      # fit
      slopes <- slope(seg_fit)$year_nr
      slope_est <- slopes[, "Est."]
      slope_se <- slopes[, "St.Err."]

      intercept_seg <- as.numeric(coef(seg_fit)[1])

      df_out <- data.frame(
        "biome_name" = unique(df_biome$biome_name),
        "agent" = unique(df_biome$agent),
        "year_nr" = df_biome$year_nr,
        "fitted" = seg_fit$fitted.values,
        "intercept" = intercept_seg,
        "model_type" = "segmented",
        "AIC" = seg_aic,
        "breakpoint" = bp,
        "breakpoint_se" = bp_se,
        "bp_ci_low" = bp_ci_low,
        "bp_ci_high" = bp_ci_high
      )

      df_out <- df_out %>%
        mutate(
          slope = ifelse(year_nr <= breakpoint, slope_est[1], slope_est[2]),
          slope_se = ifelse(year_nr <= breakpoint, slope_se[1], slope_se[2]),

          slope_ci_low = slope - 1.96 * slope_se,
          slope_ci_high = slope + 1.96 * slope_se
        )
    } else {
      lm_sum <- summary(lm_fit)$coefficients

      slope_est <- lm_sum[2, 1]
      slope_se <- lm_sum[2, 2]

      intercept_lin <- as.numeric(coef(lm_fit)[1])
      slope_ci_low <- slope_est - 1.96 * slope_se
      slope_ci_high <- slope_est + 1.96 * slope_se

      df_out <- data.frame(
        "biome_name" = unique(df_biome$biome_name),
        "agent" = unique(df_biome$agent),
        "year_nr" = df_biome$year_nr,
        "fitted" = lm_fit$fitted.values,
        "intercept" = intercept_lin,
        "model_type" = "linear",
        "AIC" = lm_aic,
        "breakpoint" = NA,
        "breakpoint_se" = NA,
        "bp_ci_low" = NA,
        "bp_ci_high" = NA,
        "slope" = slope_est,
        "slope_se" = slope_se,
        "slope_ci_low" = slope_ci_low,
        "slope_ci_high" = slope_ci_high
      )
    }
    df_res <- rbind(df_res, df_out)
  }
}


df_res <- df_res %>%
  mutate(year = year_nr + 1984) %>%
  left_join(df_re) %>%
  mutate(
    biome_name = factor(
      biome_name,
      levels = c("Boreal", "Temperate", "Mediterranean")
    ),
    fitted_upper = intercept + slope_ci_low * year_nr,
    fitted_lower = intercept + slope_ci_high * year_nr
  )

write_csv(df_res, "data/07_results/models/model_results_slopes_breakpoints.csv")

df_res <- read_csv(
  "data/07_results/models/model_results_slopes_breakpoints.csv"
)

df_ci_bp <- df_res[!duplicated(df_res$biome_name), ]
df_ci_bp <- df_ci_bp %>%
  mutate(
    biome_name = factor(
      biome_name,
      levels = c("Boreal", "Temperate", "Mediterranean")
    )
  )


df_slopes <- df_res %>%
  distinct(
    biome_name,
    agent,
    model_type,
    breakpoint,
    slope,
    slope_ci_low,
    slope_ci_high
  ) %>%
  mutate(
    significant = ifelse(slope_ci_low > 0 | slope_ci_high <= 0, TRUE, FALSE),
    slope_label = ifelse(significant, "Significant", "Not significant")
  )

# fixed slopes
effects_distrate_biomes <- data.frame(
  slope = c(
    effect_distrate_harvest_boreal,
    effect_distrate_harvest_temperate,
    effect_distrate_harvest_med,
    effect_distrate_natural_boreal,
    effect_distrate_natural_temperate,
    effect_distrate_natural_med
  ),
  biome_name = factor(
    rep(c("Boreal", "Temperate", "Mediterranean"), 2),
    levels = c("Boreal", "Temperate", "Mediterranean")
  ),
  agent = c(rep("Harvest", 3), rep("Natural", 3))
)


# Fig. 3 Time series of AGB loss sensitivity  -----------------------------

p_fig_3 <- df_res %>%
  ggplot(aes(x = year)) +
  geom_hline(
    data = effects_distrate_biomes,
    aes(yintercept = slope, color = biome_name),
    linetype = "dashed",
    linewidth = 0.5
  ) +
  geom_rect(
    data = df_ci_bp,
    aes(
      ymin = -Inf,
      ymax = Inf,
      xmin = bp_ci_low + 1984,
      xmax = bp_ci_high + 1984,
      fill = biome_name
    ),
    color = NA,
    alpha = .1
  ) +
  geom_line(aes(y = fitted, color = biome_name)) +
  geom_errorbar(
    aes(ymin = lower_95, ymax = upper_95),
    width = 0.3,
    color = "darkgrey",
    linewidth = 0.1
  ) +
  geom_point(aes(y = slope_year, color = biome_name), size = 0.9, alpha = .8) +
  geom_vline(
    aes(xintercept = breakpoint + 1984, color = biome_name),
    linetype = "dotted",
    alpha = 0.5
  ) +
  facet_rep_grid(agent ~ biome_name) +
  scale_color_manual(values = c("#648fff", "#785ef0", "#dc267f")) +
  scale_fill_manual(values = c("#648fff", "#785ef0", "#dc267f")) +
  scale_x_continuous(
    breaks = seq(1990, 2025, 10),
    minor_breaks = seq(1985, 2025, 5),
    guide = "axis_minor"
  ) +
  scale_y_continuous(labels = seq(0.3, 1.2, 0.3), breaks = seq(30, 120, 30)) +
  labs(
    y = expression(paste(
      "Sensitivity (Mg ",
      ha[F]^{
        -1
      },
      ")"
    ))
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(size = 8),
    legend.position = "None",
    strip.background = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(
      linewidth = 0.1,
      linetype = "solid",
      colour = "black"
    ),
    axis.line.y = element_line(
      linewidth = 0.1,
      linetype = "solid",
      colour = "black"
    ),
    axis.title.x = element_blank(),
    panel.spacing = unit(2, "mm"),
    axis.ticks = element_line(linewidth = 0.1),
    axis.title.y = element_text(size = 7),
    axis.text.x = element_text(color = "black", size = 7), #angle = 45, hjust = 1),
    axis.text.y = element_text(color = "black", size = 7),
    ggh4x.axis.ticks.length.minor = rel(1)
  )

ragg::agg_png(
  "figures/fig_3_sensitivity.jpg",
  width = 18,
  height = 10,
  units = "cm",
  res = 500
)
p_fig_3
dev.off()


# Tab. S3 Breakpoint identification sensitivity  --------------------------

tab_s3 <- df_res %>%
  group_by(biome_name, agent, model_type) %>%
  summarise(
    breakpoint = unique(breakpoint),
    bp_ci_low = unique(bp_ci_low),
    bp_ci_high = unique(bp_ci_high),

    slope1 = ifelse(
      first(model_type) == "segmented",
      unique(slope[year_nr <= breakpoint]),
      unique(slope)
    ),
    slope1_se = ifelse(
      first(model_type) == "segmented",
      unique(slope_se[year_nr <= breakpoint]),
      unique(slope_se)
    ),
    slope1_ci_low = ifelse(
      first(model_type) == "segmented",
      unique(slope_ci_low[year_nr <= breakpoint]),
      unique(slope_ci_low)
    ),
    slope1_ci_high = ifelse(
      first(model_type) == "segmented",
      unique(slope_ci_high[year_nr <= breakpoint]),
      unique(slope_ci_high)
    ),

    slope2 = ifelse(
      first(model_type) == "segmented",
      unique(slope[year_nr > breakpoint]),
      NA
    ),
    slope2_se = ifelse(
      first(model_type) == "segmented",
      unique(slope_se[year_nr > breakpoint]),
      NA
    ),
    slope2_ci_low = ifelse(
      first(model_type) == "segmented",
      unique(slope_ci_low[year_nr > breakpoint]),
      NA
    ),
    slope2_ci_high = ifelse(
      first(model_type) == "segmented",
      unique(slope_ci_high[year_nr > breakpoint]),
      NA
    ),
    .groups = "drop"
  ) %>%
  mutate(
    breakpoint_year = breakpoint + 1984,
    bp_ci_low_year = bp_ci_low + 1984,
    bp_ci_high_year = bp_ci_high + 1984
  ) %>%
  # Create formatted columns
  mutate(
    slope1_CI = paste0(
      round(slope1, 2),
      " [",
      round(slope1_ci_low, 2),
      ", ",
      round(slope1_ci_high, 2),
      "]"
    ),
    slope2_CI = ifelse(
      !is.na(slope2),
      paste0(
        round(slope2, 2),
        " [",
        round(slope2_ci_low, 2),
        ", ",
        round(slope2_ci_high, 2),
        "]"
      ),
      NA
    ),
    breakpoint_CI = ifelse(
      !is.na(breakpoint_year),
      paste0(
        round(breakpoint_year, 0),
        " [",
        round(bp_ci_low_year, 0),
        ", ",
        round(bp_ci_high_year, 0),
        "]"
      ),
      NA
    )
  ) %>%
  dplyr::select(
    biome_name,
    agent,
    model_type,
    slope1_CI,
    slope2_CI,
    breakpoint_CI
  )

readr::write_excel_csv(
  tab_s3,
  "tables/tab_s3_sensitivity_breakpoints.csv",
  quote = "none"
)


# 07 Posterior predictive checks ------------------------------------------

# Boreal

sims_boreal_natural <- predict(fit_natural_boreal, nsim = 100)
sims_boreal_harvest <- predict(fit_harvest_boreal, nsim = 100)

df_sims_boreal_harvest <- as.data.frame(sims_boreal_harvest)
df_sims_boreal_harvest$year <- as.numeric(rownames(sims_boreal_harvest))

df_sims_boreal_harvest <- df_sims_boreal_harvest %>%
  pivot_longer(-year)

write_csv(
  df_sims_boreal_harvest,
  "data/07_results/models/simulations/sim_100_boreal_harvest.csv"
)

df_sims_boreal_natural <- as.data.frame(sims_boreal_natural)
df_sims_boreal_natural$year <- as.numeric(rownames(sims_boreal_natural))

df_sims_boreal_natural <- df_sims_boreal_natural %>%
  pivot_longer(-year)

write_csv(
  df_sims_boreal_natural,
  "data/07_results/models/simulations/sim_100_boreal_natural.csv"
)


# Fig. S10 ----------------------------------------------------------------

p_fig_s10 <- ggplot(df_sims_boreal_harvest, aes(x = value)) +
  geom_line(
    stat = "density",
    aes(group = name),
    color = "black",
    linewidth = 0.1,
    alpha = 0.3
  ) +
  geom_line(
    data = df_boreal_harvest,
    stat = "density",
    aes(x = bm_loss_harvest),
    color = "red",
    linewidth = 0.3
  ) +
  facet_rep_wrap(~year) +
  labs(
    x = expression(paste(
      "AGB loss (Mg ",
      ha[F]^{
        -1
      },
      ")"
    )),
    y = "Density",
    title = "Boreal - Harvest"
  ) +
  scale_x_continuous(
    limits = c(0, quantile(df_boreal_harvest$bm_loss_harvest, 0.99)),
    labels = c(0, 1, 2, 3),
    breaks = c(0, 1, 2, 3)
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 8),
    strip.background = element_rect(fill = NA, color = NA),
    axis.text = element_text(size = 7, color = "black"),
    axis.title = element_text(size = 8),
    panel.grid = element_blank(),
    axis.ticks = element_line(linewidth = 0.1),
    panel.border = element_blank(),
    axis.line.x = element_line(
      linewidth = 0.1,
      linetype = "solid",
      colour = "black"
    ),
    axis.line.y = element_line(
      linewidth = 0.1,
      linetype = "solid",
      colour = "black"
    ),
    plot.title = element_text(size = 9)
  )


ragg::agg_png(
  "figures/fig_s10_sim_boreal_harvest.jpg",
  width = 16,
  height = 18,
  units = "cm",
  res = 500
)
p_fig_s10
dev.off()

# Fig. S11 ----------------------------------------------------------------

p_fig_s11 <- ggplot(df_sims_boreal_natural, aes(x = value)) +
  geom_line(
    stat = "density",
    aes(group = name),
    color = "black",
    linewidth = 0.1,
    alpha = 0.3
  ) +
  geom_line(
    data = df_boreal_natural,
    stat = "density",
    aes(x = bm_loss_natural),
    color = "red",
    linewidth = 0.3
  ) +
  facet_rep_wrap(~year) +
  labs(
    x = expression(paste(
      "AGB loss (Mg ",
      ha[F]^{
        -1
      },
      ")"
    )),
    y = "Density",
    title = "Boreal - Natural"
  ) +
  scale_x_continuous(
    limits = c(0, quantile(df_boreal_natural$bm_loss_natural, 0.99)),
    labels = c(0, 1, 2, 3),
    breaks = c(0, 1, 2, 3)
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 8),
    strip.background = element_rect(fill = NA, color = NA),
    axis.text = element_text(size = 7, color = "black"),
    axis.title = element_text(size = 8),
    panel.grid = element_blank(),
    axis.ticks = element_line(linewidth = 0.1),
    panel.border = element_blank(),
    axis.line.x = element_line(
      linewidth = 0.1,
      linetype = "solid",
      colour = "black"
    ),
    axis.line.y = element_line(
      linewidth = 0.1,
      linetype = "solid",
      colour = "black"
    ),
    plot.title = element_text(size = 9)
  )

ragg::agg_png(
  "figures/fig_s10_sim_boreal_natural.jpg",
  width = 16,
  height = 18,
  units = "cm",
  res = 500
)
p_fig_s11
dev.off()


# Temperate

sims_temperate_natural <- predict(fit_natural_temperate, nsim = 100)
sims_temperate_harvest <- predict(fit_harvest_temperate, nsim = 100)

df_sims_temperate_harvest <- as.data.frame(sims_temperate_harvest)
df_sims_temperate_harvest$year <- as.numeric(rownames(sims_temperate_harvest))

df_sims_temperate_harvest <- df_sims_temperate_harvest %>%
  pivot_longer(-year)

write_csv(
  df_sims_temperate_harvest,
  "data/07_results/models/simulations/sim_100_temperate_harvest.csv"
)

df_sims_temperate_natural <- as.data.frame(sims_temperate_natural)
df_sims_temperate_natural$year <- as.numeric(rownames(sims_temperate_natural))

df_sims_temperate_natural <- df_sims_temperate_natural %>%
  pivot_longer(-year)

write_csv(
  df_sims_temperate_natural,
  "data/07_results/models/simulations/sim_100_temperate_natural.csv"
)

# Fig. S12 ----------------------------------------------------------------

p_fig_s12 <- ggplot(df_sims_temperate_harvest, aes(x = value)) +
  geom_line(
    stat = "density",
    aes(group = name),
    color = "black",
    linewidth = 0.1,
    alpha = 0.3
  ) +
  geom_line(
    data = df_temperate_harvest,
    stat = "density",
    aes(x = bm_loss_harvest),
    color = "red",
    linewidth = 0.3
  ) +
  facet_rep_wrap(~year) +
  labs(
    x = expression(paste(
      "AGB loss (Mg ",
      ha[F]^{
        -1
      },
      ")"
    )),
    y = "Density",
    title = "Temperate - Harvest"
  ) +
  scale_x_continuous(
    limits = c(0, quantile(df_temperate_harvest$bm_loss_harvest, 0.99)),
    labels = c(0, 1, 2, 3),
    breaks = c(0, 1, 2, 3)
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 8),
    strip.background = element_rect(fill = NA, color = NA),
    axis.text = element_text(size = 7, color = "black"),
    axis.title = element_text(size = 8),
    panel.grid = element_blank(),
    axis.ticks = element_line(linewidth = 0.1),
    panel.border = element_blank(),
    axis.line.x = element_line(
      linewidth = 0.1,
      linetype = "solid",
      colour = "black"
    ),
    axis.line.y = element_line(
      linewidth = 0.1,
      linetype = "solid",
      colour = "black"
    ),
    plot.title = element_text(size = 9)
  )


ragg::agg_png(
  "figures/fig_s12_sim_temperate_harvest.jpg",
  width = 16,
  height = 18,
  units = "cm",
  res = 500
)
p_fig_s12
dev.off()

# Fig. S13 ----------------------------------------------------------------

p_fig_s13 <- ggplot(df_sims_temperate_natural, aes(x = value)) +
  geom_line(
    stat = "density",
    aes(group = name),
    color = "black",
    linewidth = 0.1,
    alpha = 0.3
  ) +
  geom_line(
    data = df_temperate_natural,
    stat = "density",
    aes(x = bm_loss_natural),
    color = "red",
    linewidth = 0.3
  ) +
  facet_rep_wrap(~year) +
  labs(
    x = expression(paste(
      "AGB loss (Mg ",
      ha[F]^{
        -1
      },
      ")"
    )),
    y = "Density",
    title = "Temperate - Natural"
  ) +
  scale_x_continuous(
    limits = c(0, quantile(df_temperate_natural$bm_loss_natural, 0.99)),
    labels = c(0, 1, 2, 3),
    breaks = c(0, 1, 2, 3)
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 8),
    strip.background = element_rect(fill = NA, color = NA),
    axis.text = element_text(size = 7, color = "black"),
    axis.title = element_text(size = 8),
    panel.grid = element_blank(),
    axis.ticks = element_line(linewidth = 0.1),
    panel.border = element_blank(),
    axis.line.x = element_line(
      linewidth = 0.1,
      linetype = "solid",
      colour = "black"
    ),
    axis.line.y = element_line(
      linewidth = 0.1,
      linetype = "solid",
      colour = "black"
    ),
    plot.title = element_text(size = 9)
  )

ragg::agg_png(
  "figures/fig_s12_sim_temperate_natural.jpg",
  width = 16,
  height = 18,
  units = "cm",
  res = 500
)
p_fig_s13
dev.off()


# Mediterranean
sims_med_natural <- predict(fit_natural_med, nsim = 100)
sims_med_harvest <- predict(fit_harvest_med, nsim = 100)

df_sims_med_harvest <- as.data.frame(sims_med_harvest)
df_sims_med_harvest$year <- as.numeric(rownames(sims_med_harvest))

df_sims_med_harvest <- df_sims_med_harvest %>%
  pivot_longer(-year)

write_csv(
  df_sims_med_harvest,
  "data/07_results/models/simulations/sim_100_med_harvest.csv"
)

df_sims_med_natural <- as.data.frame(sims_med_natural)
df_sims_med_natural$year <- as.numeric(rownames(sims_med_natural))

df_sims_med_natural <- df_sims_med_natural %>%
  pivot_longer(-year)

write_csv(
  df_sims_med_natural,
  "data/07_results/models/simulations/sim_100_med_natural.csv"
)

# Fig. S14 ----------------------------------------------------------------

p_fig_s14 <- ggplot(df_sims_med_harvest, aes(x = value)) +
  geom_line(
    stat = "density",
    aes(group = name),
    color = "black",
    linewidth = 0.1,
    alpha = 0.3
  ) +
  geom_line(
    data = df_med_harvest,
    stat = "density",
    aes(x = bm_loss_harvest),
    color = "red",
    linewidth = 0.3
  ) +
  facet_rep_wrap(~year) +
  labs(
    x = expression(paste(
      "AGB loss (Mg ",
      ha[F]^{
        -1
      },
      ")"
    )),
    y = "Density",
    title = "Mediterranean - Harvest"
  ) +
  scale_x_continuous(
    limits = c(0, quantile(df_med_harvest$bm_loss_harvest, 0.99)),
    labels = c(0, 1, 2, 3),
    breaks = c(0, 1, 2, 3)
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 8),
    strip.background = element_rect(fill = NA, color = NA),
    axis.text = element_text(size = 7, color = "black"),
    axis.title = element_text(size = 8),
    panel.grid = element_blank(),
    axis.ticks = element_line(linewidth = 0.1),
    panel.border = element_blank(),
    axis.line.x = element_line(
      linewidth = 0.1,
      linetype = "solid",
      colour = "black"
    ),
    axis.line.y = element_line(
      linewidth = 0.1,
      linetype = "solid",
      colour = "black"
    ),
    plot.title = element_text(size = 9)
  )


ragg::agg_png(
  "figures/fig_s14_sim_med_harvest.jpg",
  width = 16,
  height = 18,
  units = "cm",
  res = 500
)
p_fig_s14
dev.off()

# Fig. S15 ----------------------------------------------------------------

p_fig_s15 <- ggplot(df_sims_med_natural, aes(x = value)) +
  geom_line(
    stat = "density",
    aes(group = name),
    color = "black",
    linewidth = 0.1,
    alpha = 0.3
  ) +
  geom_line(
    data = df_med_natural,
    stat = "density",
    aes(x = bm_loss_natural),
    color = "red",
    linewidth = 0.3
  ) +
  facet_rep_wrap(~year) +
  labs(
    x = expression(paste(
      "AGB loss (Mg ",
      ha[F]^{
        -1
      },
      ")"
    )),
    y = "Density",
    title = "Mediterranean - Natural"
  ) +
  scale_x_continuous(
    limits = c(0, quantile(df_med_natural$bm_loss_natural, 0.99)),
    labels = c(0, 1, 2, 3),
    breaks = c(0, 1, 2, 3)
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 8),
    strip.background = element_rect(fill = NA, color = NA),
    axis.text = element_text(size = 7, color = "black"),
    axis.title = element_text(size = 8),
    panel.grid = element_blank(),
    axis.ticks = element_line(linewidth = 0.1),
    panel.border = element_blank(),
    axis.line.x = element_line(
      linewidth = 0.1,
      linetype = "solid",
      colour = "black"
    ),
    axis.line.y = element_line(
      linewidth = 0.1,
      linetype = "solid",
      colour = "black"
    ),
    plot.title = element_text(size = 9)
  )

ragg::agg_png(
  "figures/fig_s15_sim_med_natural.jpg",
  width = 16,
  height = 18,
  units = "cm",
  res = 500
)
p_fig_s15
dev.off()


# 08 Spatial autocorrelation tests ----------------------------------------

# get residuals for each model
df_temperate_harvest$residuals <- residuals(
  fit_harvest_temperate,
  type = "pearson"
)
df_temperate_natural$residuals <- residuals(
  fit_natural_temperate,
  type = "pearson"
)

df_boreal_harvest$residuals <- residuals(fit_harvest_boreal, type = "pearson")
df_boreal_natural$residuals <- residuals(fit_natural_boreal, type = "pearson")

df_med_harvest$residuals <- residuals(fit_harvest_med, type = "pearson")
df_med_natural$residuals <- residuals(fit_natural_med, type = "pearson")

coordinates <- c("X", "Y")

distances <- c(21, 41, 81, 151, 251) * 1000

# Temperate

mi_temperate_natural <- map_df(distances, function(dist) {
  df_temperate_natural %>%
    filter(is.finite(residuals)) %>%
    group_by(year) %>%
    group_split() %>%
    map_df(function(df_year) {
      coords <- as.matrix(df_year[, coordinates])
      nb <- dnearneigh(coords, 0, dist)
      lw <- nb2listw(nb, style = "W", zero.policy = TRUE)

      moran <- moran.test(df_year$residuals, lw, zero.policy = TRUE)

      se <- sqrt(moran$estimate[["Variance"]])
      ci_low <- moran$estimate[["Moran I statistic"]] - 1.96 * se
      ci_high <- moran$estimate[["Moran I statistic"]] + 1.96 * se

      tibble(
        year = unique(df_year$year),
        distance = dist,
        morans_i = moran$estimate[["Moran I statistic"]],
        ci_low = ci_low,
        ci_high = ci_high,
        p_value = moran$p.value,
        agent = "Natural",
        biome_name = "temperate"
      )
    })
})

mi_temperate_natural <- mi_temperate_natural %>%
  mutate(
    significance = ifelse(p_value < 0.05, "p < 0.05", "p >= 0.05"),
  )

write_csv(
  mi_temperate_natural,
  "data/07_results/models/temperate/morans_i_residuals_temperate_natural.csv"
)


mi_temperate_harvest <- map_df(distances, function(dist) {
  df_temperate_harvest %>%
    filter(is.finite(residuals)) %>%
    group_by(year) %>%
    group_split() %>%
    map_df(function(df_year) {
      coords <- as.matrix(df_year[, coordinates])
      nb <- dnearneigh(coords, 0, dist)
      lw <- nb2listw(nb, style = "W", zero.policy = TRUE)

      moran <- moran.test(df_year$residuals, lw, zero.policy = TRUE)

      se <- sqrt(moran$estimate[["Variance"]])
      ci_low <- moran$estimate[["Moran I statistic"]] - 1.96 * se
      ci_high <- moran$estimate[["Moran I statistic"]] + 1.96 * se

      tibble(
        year = unique(df_year$year),
        distance = dist,
        morans_i = moran$estimate[["Moran I statistic"]],
        ci_low = ci_low,
        ci_high = ci_high,
        p_value = moran$p.value,
        agent = "harvest",
        biome_name = "temperate"
      )
    })
})

mi_temperate_harvest <- mi_temperate_harvest %>%
  mutate(
    significance = ifelse(p_value < 0.05, "p < 0.05", "p >= 0.05"),
  )

write_csv(
  mi_temperate_harvest,
  "data/07_results/models/temperate/morans_i_residuals_temperate_harvest.csv"
)


# Boreal

mi_boreal_natural <- map_df(distances, function(dist) {
  df_boreal_natural %>%
    filter(is.finite(residuals)) %>%
    group_by(year) %>%
    group_split() %>%
    map_df(function(df_year) {
      coords <- as.matrix(df_year[, coordinates])
      nb <- dnearneigh(coords, 0, dist)
      lw <- nb2listw(nb, style = "W", zero.policy = TRUE)

      moran <- moran.test(df_year$residuals, lw, zero.policy = TRUE)

      se <- sqrt(moran$estimate[["Variance"]])
      ci_low <- moran$estimate[["Moran I statistic"]] - 1.96 * se
      ci_high <- moran$estimate[["Moran I statistic"]] + 1.96 * se

      tibble(
        year = unique(df_year$year),
        distance = dist,
        morans_i = moran$estimate[["Moran I statistic"]],
        ci_low = ci_low,
        ci_high = ci_high,
        p_value = moran$p.value,
        agent = "Natural",
        biome_name = "Boreal"
      )
    })
})

mi_boreal_natural <- mi_boreal_natural %>%
  mutate(
    significance = ifelse(p_value < 0.05, "p < 0.05", "p >= 0.05"),
  )

write_csv(
  mi_boreal_natural,
  "data/07_results/models/boreal/morans_i_residuals_boreal_natural.csv"
)


mi_boreal_harvest <- map_df(distances, function(dist) {
  df_boreal_harvest %>%
    filter(is.finite(residuals)) %>%
    group_by(year) %>%
    group_split() %>%
    map_df(function(df_year) {
      coords <- as.matrix(df_year[, coordinates])
      nb <- dnearneigh(coords, 0, dist)
      lw <- nb2listw(nb, style = "W", zero.policy = TRUE)

      moran <- moran.test(df_year$residuals, lw, zero.policy = TRUE)

      se <- sqrt(moran$estimate[["Variance"]])
      ci_low <- moran$estimate[["Moran I statistic"]] - 1.96 * se
      ci_high <- moran$estimate[["Moran I statistic"]] + 1.96 * se

      tibble(
        year = unique(df_year$year),
        distance = dist,
        morans_i = moran$estimate[["Moran I statistic"]],
        ci_low = ci_low,
        ci_high = ci_high,
        p_value = moran$p.value,
        agent = "harvest",
        biome_name = "Boreal"
      )
    })
})

mi_boreal_harvest <- mi_boreal_harvest %>%
  mutate(
    significance = ifelse(p_value < 0.05, "p < 0.05", "p >= 0.05"),
  )

write_csv(
  mi_boreal_harvest,
  "data/07_results/models/boreal/morans_i_residuals_boreal_harvest.csv"
)


# Mediterranean
mi_med_natural <- map_df(distances, function(dist) {
  df_med_natural %>%
    filter(is.finite(residuals)) %>%
    group_by(year) %>%
    group_split() %>%
    map_df(function(df_year) {
      coords <- as.matrix(df_year[, coordinates])
      nb <- dnearneigh(coords, 0, dist)
      lw <- nb2listw(nb, style = "W", zero.policy = TRUE)

      moran <- moran.test(df_year$residuals, lw, zero.policy = TRUE)

      se <- sqrt(moran$estimate[["Variance"]])
      ci_low <- moran$estimate[["Moran I statistic"]] - 1.96 * se
      ci_high <- moran$estimate[["Moran I statistic"]] + 1.96 * se

      tibble(
        year = unique(df_year$year),
        distance = dist,
        morans_i = moran$estimate[["Moran I statistic"]],
        ci_low = ci_low,
        ci_high = ci_high,
        p_value = moran$p.value,
        agent = "Natural",
        biome_name = "med"
      )
    })
})

mi_med_natural <- mi_med_natural %>%
  mutate(
    significance = ifelse(p_value < 0.05, "p < 0.05", "p >= 0.05"),
  )

write_csv(
  mi_med_natural,
  "data/07_results/models/med/morans_i_residuals_med_natural.csv"
)


mi_med_harvest <- map_df(distances, function(dist) {
  df_med_harvest %>%
    filter(is.finite(residuals)) %>%
    group_by(year) %>%
    group_split() %>%
    map_df(function(df_year) {
      coords <- as.matrix(df_year[, coordinates])
      nb <- dnearneigh(coords, 0, dist)
      lw <- nb2listw(nb, style = "W", zero.policy = TRUE)

      moran <- moran.test(df_year$residuals, lw, zero.policy = TRUE)

      se <- sqrt(moran$estimate[["Variance"]])
      ci_low <- moran$estimate[["Moran I statistic"]] - 1.96 * se
      ci_high <- moran$estimate[["Moran I statistic"]] + 1.96 * se

      tibble(
        year = unique(df_year$year),
        distance = dist,
        morans_i = moran$estimate[["Moran I statistic"]],
        ci_low = ci_low,
        ci_high = ci_high,
        p_value = moran$p.value,
        agent = "harvest",
        biome_name = "med"
      )
    })
})

mi_med_harvest <- mi_med_harvest %>%
  mutate(
    significance = ifelse(p_value < 0.05, "p < 0.05", "p >= 0.05"),
  )

write_csv(
  mi_med_harvest,
  "data/07_results/models/mediterranean/morans_i_residuals_med_harvest.csv"
)

# combine all results
mi_europe_agents <- bind_rows(
  mi_med_harvest,
  mi_med_natural,
  mi_temperate_harvest,
  mi_temperate_natural,
  mi_boreal_harvest,
  mi_boreal_natural
) %>%
  mutate(
    significance = ifelse(p_value < 0.05, "p < 0.05", "p >= 0.05"),
    biome_name = ifelse(
      biome_name == "med",
      "Mediterranean",
      ifelse(biome_name == "temperate", "Temperate", biome_name)
    ),
    biome_name = factor(
      biome_name,
      levels = c("Boreal", "Temperate", "Mediterranean")
    )
  )


# Fig. S16 Spatial autocorrelation of model residuals ---------------------

p_fig_s16 <- mi_europe_agents %>%
  filter(distance > 21000) %>%
  mutate(
    distance = factor(
      distance,
      levels = c(41000, 81000, 151000, 251000),
      labels = c(41, 81, 151, 251)
    ),
    agent = ifelse(agent == "harvest", "Harvest", agent)
  ) %>%
  ggplot() +
  geom_point(
    aes(x = year, y = morans_i, color = distance, shape = significance),
    size = 1.2,
    alpha = 0.8
  ) +
  scale_shape_manual(
    values = c(4, 16),
    name = "",
    labels = c("p < 0.05", expression(p >= 0.05))
  ) +
  scale_color_manual(
    values = RColorBrewer::brewer.pal(n = 8, "YlGnBu")[4:8],
    name = "Distance (km)"
  ) +
  facet_rep_grid(agent ~ biome_name) +
  labs(y = "Moran's I") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom",
    axis.title.x = element_blank(),
    axis.title = element_text(size = 7),
    axis.text.y = element_text(size = 7, color = "black"),
    axis.text.x = element_text(size = 7, color = "black"),
    strip.background = element_rect(fill = NA, color = NA),
    strip.text = element_text(size = 8),
    axis.ticks = element_line(linewidth = 0.1),
    panel.border = element_blank(),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 7),
    legend.key.height = unit(3, "mm"),
    legend.key.width = unit(3, "mm"),
    axis.line.x = element_line(
      linewidth = 0.1,
      linetype = "solid",
      colour = "black"
    ),
    axis.line.y = element_line(
      linewidth = 0.1,
      linetype = "solid",
      colour = "black"
    ),
    plot.tag = element_text(size = 7)
  )

ragg::agg_png(
  "figures/fig_s16_spatial_autocorrelation.jpg",
  width = 16,
  height = 14,
  units = "cm",
  res = 500
)
p_fig_s16
dev.off()
