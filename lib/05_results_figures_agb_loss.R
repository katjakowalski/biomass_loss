# ---
# Title:        Figures and tables for total AGB loss
# Description:  Creates and saves figures 2, S2, S4, S9, and tables S1, S2
# Author:       Katja Kowalski
# Affiliation:  Technical University of Munich, School of Life Sciences
# Date Created: 2025-08-01
# Last Updated: 2025-08-04
# ---

library(tidyverse)
library(ragg)
library(ggpubr)
library(sf)

# Prepare data ------------------------------------------------------------

df <- read_csv("data/06_agb_loss/agb_loss_planet_20km_filtered.csv")

df_full <- read_csv("data/06_agb_loss/agb_loss_planet_20km.csv")

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


# 01 AGB loss total Europe ------------------------------------------------

bm_loss_europe <- df %>%
  ungroup() %>%
  summarize(
    loss_sum = sum(planet_bm_loss, na.rm = T),
    loss_sd = sum(planet_bm_loss_sd, na.rm = T)
  ) %>%
  pivot_longer(cols = everything()) %>%
  mutate(type = str_split_i(name, "_", 2)) %>%
  dplyr::select(-name) %>%
  pivot_wider(names_from = type, values_from = value)

# biomass loss estimate +/- sd
round(bm_loss_europe$sum * 1e-09, 1)
round(bm_loss_europe$sd * 1e-09, 1)

# comparison to Pugh et al. 2019
df %>%
  filter(year >= 2001 & year <= 2014) %>%
  summarize(
    loss_sum = round(
      sum(planet_bm_loss, na.rm = T) * 1e-09 / length(unique(year)),
      2
    ),
    loss_sd = round(
      sum(planet_bm_loss_sd, na.rm = T) * 1e-09 / length(unique(year)),
      2
    )
  ) %>%
  mutate(
    percent_loss_global = loss_sum / 2 * 100,
    percent_loss_global_sd = loss_sd / 2 * 100
  )


# 02 AGB loss per ecoregion -----------------------------------------------

tab_europe <- df %>%
  summarize(
    agbd_before_mean = round(mean(planet_agbd_before_mean, na.rm = T), 1),
    agbd_before_sd = round(sd(planet_agbd_before_mean, na.rm = T), 1),

    agbd_loss_mean = round(mean(planet_mean_sev_abs, na.rm = T), 1),
    agbd_loss_sd = round(sd(planet_mean_sev_abs, na.rm = T), 1),

    agbd_loss_rel_mean = round(mean(planet_mean_sev_rel * 100, na.rm = T), 1),
    agbd_loss_rel_sd = round(sd(planet_mean_sev_rel * 100, na.rm = T), 1),

    agb_loss_total = round(sum(planet_bm_loss, na.rm = T) * 1e-09, 1),
    agb_loss_total_sd = round(sum(planet_bm_loss_sd, na.rm = T) * 1e-09, 1),

    agb_loss_haf = round(
      sum(planet_bm_loss, na.rm = T) / sum(forest_ha_sum),
      1
    ),
    agb_loss_haf_sd = round(
      sum(planet_bm_loss_sd, na.rm = T) / sum(forest_ha_sum),
      1
    )
  )

biome_order <- c("Boreal", "Temperate", "Mediterranean", "Europe")
tab1 <- df %>%
  group_by(biome_name) %>%
  summarize(
    agbd_before_mean = round(mean(planet_agbd_before_mean, na.rm = T), 1),
    agbd_before_sd = round(sd(planet_agbd_before_mean, na.rm = T), 1),

    agbd_loss_mean = round(mean(planet_mean_sev_abs, na.rm = T), 1),
    agbd_loss_sd = round(sd(planet_mean_sev_abs, na.rm = T), 1),

    agbd_loss_rel_mean = round(mean(planet_mean_sev_rel * 100, na.rm = T), 1),
    agbd_loss_rel_sd = round(sd(planet_mean_sev_rel * 100, na.rm = T), 1),

    agb_loss_total = round(sum(planet_bm_loss, na.rm = T) * 1e-09, 1),
    agb_loss_total_sd = round(sum(planet_bm_loss_sd, na.rm = T) * 1e-09, 1),

    agb_loss_haf = round(
      sum(planet_bm_loss, na.rm = T) / sum(forest_ha_sum),
      1
    ),
    agb_loss_haf_sd = round(
      sum(planet_bm_loss_sd, na.rm = T) / sum(forest_ha_sum),
      1
    )
  ) %>%
  drop_na() %>%
  add_row(biome_name = "Europe", tab_europe) %>%
  slice(match(biome_order, biome_name)) %>%
  mutate(
    agbd_before_mean = paste(agbd_before_mean, "\u00B1", agbd_before_sd),
    agbd_loss_mean = paste(agbd_loss_mean, "\u00B1", agbd_loss_sd),
    agbd_loss_rel_mean = paste(agbd_loss_rel_mean, "\u00B1", agbd_loss_rel_sd),
    agb_loss_total = paste(agb_loss_total, "\u00B1", agb_loss_total_sd),
    agb_loss_haf = paste(agb_loss_haf, "\u00B1", agb_loss_haf_sd)
  ) %>%
  dplyr::select(
    biome_name,
    agbd_before_mean,
    agbd_loss_mean,
    agbd_loss_rel_mean,
    agb_loss_total,
    agb_loss_haf
  )

write.csv(
  tab1,
  "tables/tab1_agbd_loss_ecoregions.csv",
  row.names = F,
  quote = F
)


# 03 AGB loss per ecoregion and agent -------------------------------------

bm_loss_agent_biome <- df %>%
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

# estimate share of natural disturbances relative to harvest AGB loss (comparison Schelhaas et al. 2003)

bm_loss_natural <- bm_loss_agent_biome %>%
  ungroup() %>%
  summarize(
    bm_loss_harvest = sum(bm_loss_harvest, na.rm = T),
    bm_loss_natural = sum(bm_loss_natural, na.rm = T)
  ) %>%
  mutate(
    percent_natural = round(bm_loss_natural / bm_loss_harvest * 100, 2)
  )

# estimate of total biomass loss per for harvest and natural dist.
bm_loss_agent_total <- bm_loss_agent_biome %>%
  ungroup() %>%
  summarize(
    bm_loss_harvest = sum(bm_loss_harvest, na.rm = T),
    bm_loss_harvest_sd = sum(bm_loss_harvest_sd, na.rm = T),

    bm_loss_natural = sum(bm_loss_natural, na.rm = T),
    bm_loss_natural_sd = sum(bm_loss_natural_sd, na.rm = T),

    dist_area_harvest = sum(dist_area_harvest, na.rm = T),
    dist_area_natural = sum(dist_area_natural, na.rm = T)
  )

share_harvest_area <- bm_loss_agent_total$dist_area_harvest /
  (bm_loss_agent_total$dist_area_harvest +
    bm_loss_agent_total$dist_area_natural)
round(share_harvest_area * 100, 1)

# biomass loss estimate +/- sd
round(bm_loss_agent_total$bm_loss_harvest * 1e-09, 1)
round(bm_loss_agent_total$bm_loss_harvest_sd * 1e-09, 1)

round(bm_loss_agent_total$bm_loss_natural * 1e-09, 1)
round(bm_loss_agent_total$bm_loss_natural_sd * 1e-09, 1)

# percentages

round(
  (bm_loss_agent_total$bm_loss_harvest) /
    (bm_loss_agent_total$bm_loss_harvest +
      bm_loss_agent_total$bm_loss_natural) *
    100,
  1
)
round(
  (bm_loss_agent_total$bm_loss_natural) /
    (bm_loss_agent_total$bm_loss_harvest +
      bm_loss_agent_total$bm_loss_natural) *
    100,
  1
)


# percentages of total AGB loss for ecoregions and agents

bm_loss_agent_biome_percent <- bm_loss_agent_biome %>%
  group_by(biome_name) %>%
  summarize(
    bm_loss_harvest = sum(bm_loss_harvest, na.rm = T),
    bm_loss_harvest_sd = sum(bm_loss_harvest_sd, na.rm = T),

    bm_loss_natural = sum(bm_loss_natural, na.rm = T),
    bm_loss_natural_sd = sum(bm_loss_natural_sd, na.rm = T),

    dist_area_harvest = sum(dist_area_harvest, na.rm = T),
    dist_area_natural = sum(dist_area_natural, na.rm = T)
  ) %>%
  mutate(
    percentage_harvest = round(
      bm_loss_harvest / (bm_loss_harvest + bm_loss_natural) * 100,
      2
    ),
    percentage_natural = round(
      bm_loss_natural / (bm_loss_harvest + bm_loss_natural) * 100,
      2
    )
  ) %>%
  dplyr::select(
    biome_name,
    percentage_harvest,
    percentage_natural
  )


# Fig. 2 AGB losses from natural disturbances and harvests 1985-2023 -------

# Fig. 2a Time series of AGB loss -----------------------------------------

bm_loss_biome <- df %>%
  group_by(year, biome_name) %>%
  summarize(
    loss_sum = sum(planet_bm_loss, na.rm = T) / sum(forest_ha_sum),
    loss_sd = sum(planet_bm_loss_sd, na.rm = T) / sum(forest_ha_sum),
    year = unique(year)
  )

bm_loss_biome <- bm_loss_biome %>%
  pivot_longer(-c(year, biome_name)) %>%
  mutate(type = str_split_i(name, "_", 2)) %>%
  dplyr::select(-name) %>%
  pivot_wider(names_from = type, values_from = value)

# get data for boreal 2018 2023
bm_loss_boreal_18_23 <- df_full %>%
  group_by(year, biome_name) %>%
  summarize(
    loss_sum = sum(planet_bm_loss, na.rm = T) / sum(forest_ha_sum),
    loss_sd = sum(planet_bm_loss_sd, na.rm = T) / sum(forest_ha_sum),
    year = unique(year)
  ) %>%
  filter(biome_name == "Boreal" & (year == 2018 | year == 2023))

bm_loss_boreal_18_23 <- bm_loss_boreal_18_23 %>%
  pivot_longer(-c(year, biome_name)) %>%
  mutate(type = str_split_i(name, "_", 2)) %>%
  dplyr::select(-name) %>%
  pivot_wider(names_from = type, values_from = value) %>%
  mutate(
    biome_name = factor(
      biome_name,
      levels = c("Boreal", "Temperate", "Mediterranean")
    )
  )

# add mean bm loss per time period
mean_loss_biome_tp <- bm_loss_biome %>%
  mutate(
    biome_name = factor(
      biome_name,
      levels = c("Boreal", "Temperate", "Mediterranean")
    )
  ) %>%
  ungroup() %>%
  mutate(
    time_period = ifelse(
      year < 2000,
      "1985-2000",
      ifelse(year < 2018, "2000-2018", "2018-2023")
    )
  ) %>%
  group_by(time_period, biome_name) %>%
  summarize(
    loss_sum = mean(sum)
  ) %>%
  ungroup() %>%
  mutate(
    start_yr = as.numeric(str_split_i(time_period, "-", 1)),
    end_yr = as.numeric(str_split_i(time_period, "-", 2))
  )


p_loss_biome_ts <- bm_loss_biome %>%
  mutate(
    biome_name = factor(
      biome_name,
      levels = c("Boreal", "Temperate", "Mediterranean")
    )
  ) %>%
  ggplot(aes(x = year, y = sum)) +
  geom_point(alpha = 1, aes(color = biome_name, shape = biome_name), size = 1) +
  geom_point(
    data = bm_loss_boreal_18_23,
    aes(color = biome_name),
    size = 1,
    shape = 4
  ) +
  geom_line(aes(color = biome_name), alpha = .3, linewidth = 0.3) +
  geom_segment(
    data = mean_loss_biome_tp,
    aes(
      x = start_yr,
      xend = end_yr,
      y = loss_sum,
      color = factor(biome_name)
    ),
    linewidth = 0.5
  ) +
  scale_color_manual(values = c("#1192e8", "#6929c4", "#ee538b")) +
  scale_x_continuous(labels = seq(1985, 2023, 5), breaks = seq(1985, 2023, 5)) +
  labs(
    y = expression(paste(
      "AGB loss (Mg ",
      ha[F]^{
        -1
      },
      yr.^{
        -1
      },
      ")"
    ))
  ) +
  coord_cartesian(clip = 'off', xlim = c(1985, 2023)) +
  theme_bw() +
  theme(
    plot.margin = margin(2, 2, -5, 2, unit = "mm"),
    panel.grid = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 7),
    axis.text = element_text(size = 7, color = "black"),
    strip.background = element_rect(fill = NA, color = NA),
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
    legend.position = c(0.15, 0.9),
    legend.title = element_blank(),
    legend.text = element_text(size = 7, color = "black"),
    legend.direction = "vertical",
    legend.key.spacing.y = unit(0.07, "cm"),
    legend.key.height = unit(0.1, "cm"),
    legend.key.width = unit(0.25, "cm"),
    axis.ticks = element_line(linewidth = 0.1),
    plot.tag = element_text(size = 7)
  )


# Fig 2b ------------------------------------------------------------------

df_plot_agent_biome <- bm_loss_agent_biome %>%
  group_by(year, biome_name) %>%
  summarize(
    bm_loss_harvest = sum(bm_loss_harvest, na.rm = T) /
      sum(forest_ha_sum, na.rm = T),
    bm_loss_natural = sum(bm_loss_natural, na.rm = T) /
      sum(forest_ha_sum, na.rm = T),

    dist_area_harvest = sum(dist_area_harvest, na.rm = T) /
      sum(forest_ha_sum, na.rm = T),
    dist_area_natural = sum(dist_area_natural, na.rm = T) /
      sum(forest_ha_sum, na.rm = T)
  ) %>%
  group_by(biome_name) %>%
  summarize(
    bm_loss_harvest = sum(bm_loss_harvest, na.rm = T),
    bm_loss_natural = sum(bm_loss_natural, na.rm = T),

    dist_area_harvest = sum(dist_area_harvest, na.rm = T),
    dist_area_natural = sum(dist_area_natural, na.rm = T)
  ) %>%
  dplyr::select(biome_name, bm_loss_harvest, bm_loss_natural) %>%
  pivot_longer(-biome_name)


bm_loss_europe <- df %>%
  group_by(biome_name, year) %>%
  summarize(
    loss_sum = sum(planet_bm_loss, na.rm = T),
    loss_sd = sum(planet_bm_loss_sd, na.rm = T),
    forest_ha_sum = sum(forest_ha_sum, na.rm = T)
  ) %>%
  group_by(biome_name) %>%
  summarize(
    loss_sum = sum(loss_sum) / mean(forest_ha_sum),
    loss_sd = sum(loss_sd) / mean(forest_ha_sum)
  ) %>%
  pivot_longer(-biome_name) %>%
  mutate(type = str_split_i(name, "_", 2)) %>%
  dplyr::select(-name) %>%
  pivot_wider(names_from = type, values_from = value) %>%
  rename("loss_sum" = "sum", "loss_sd" = "sd")

top_labels <- data.frame(
  "biome_name" = c("Temperate", "Boreal", "Mediterranean"),
  "loss_sum" = rep(43, 3)
)

p_biomes <- bm_loss_europe %>%
  left_join(df_plot_agent_biome) %>%
  pivot_wider(names_from = name, values_from = value) %>%
  mutate(bm_total = bm_loss_harvest + bm_loss_natural) %>%
  pivot_longer(-c(biome_name, loss_sum, loss_sd, bm_total)) %>%
  mutate(
    biome_name = factor(
      biome_name,
      levels = c("Temperate", "Boreal", "Mediterranean")
    ),
    name = factor(
      name,
      levels = c("bm_loss_harvest", "bm_loss_natural"),
      labels = c("Harvest", "Natural")
    )
  ) %>%
  ggplot(aes(x = biome_name, y = value)) +
  geom_bar(stat = "identity", width = 0.5, aes(fill = name)) + #
  geom_text(
    data = top_labels,
    aes(x = biome_name, y = loss_sum, label = biome_name),
    size = 2.5
  ) +
  geom_errorbar(
    aes(
      x = biome_name,
      ymin = (bm_total - loss_sd),
      ymax = (bm_total + loss_sd)
    ),
    width = 0.1,
    linewidth = 0.1
  ) +
  scale_fill_manual(values = c("grey70", "grey30")) +
  labs(
    y = expression(paste(
      "AGB loss (Mg ",
      ha[F]^{
        -1
      },
      ")"
    ))
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "None",
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 7),
    axis.text.y = element_text(size = 7, color = "black"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.background = element_rect(fill = NA, color = NA),
    strip.text = element_text(size = 7),
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
    )
  )


# Fig. 2c AGB loss maps ---------------------------------------------------

europe <- read_sf("data/gis/borders/europe_outline_3035.shp")

europe <- st_transform(
  europe,
  crs(grid)
)

grid <- read_sf(
  "data/gis/referencegrid/grid_20km_3035.shp"
)


grid_plt <- grid %>%
  st_intersection(europe) %>%
  left_join(df) %>%
  mutate(planet_bm_loss = planet_bm_loss / forest_ha_sum)

# 1985-1999
grid_plt_1985 <- grid_plt %>%
  filter(year <= 1999) %>%
  group_by(gridid) %>%
  summarize(
    planet_bm_loss = mean(planet_bm_loss, na.rm = T),
    forest_share = mean(forest_share, na.rm = T)
  )

p_1985 <- ggplot(grid_plt_1985) +
  geom_sf(data = europe, fill = "grey75", color = NA) +
  geom_sf(aes(fill = planet_bm_loss), color = NA) +
  scale_fill_gradientn(
    colors = viridisLite::plasma(10),
    limits = c(0, 4),
    oob = scales::squish,
    na.value = NA
  ) +
  annotate(
    geom = 'text',
    label = '1985-1999',
    x = -Inf,
    y = Inf,
    hjust = 0,
    vjust = 5,
    size = 2.5
  ) +
  theme_void() +
  theme(
    panel.border = element_rect(color = NA, fill = NA),
    legend.position = "None"
  ) +
  guides(
    color = guide_colourbar(
      theme = theme(
        legend.key.width = unit(3, "lines"),
        legend.key.height = unit(0.5, "lines")
      ),
      legend.ticks.length = unit(.05, "cm"),
      title.position = "top"
    )
  )


# 2000-2017
grid_plt_2000 <- grid_plt %>%
  filter(year <= 2017 & year >= 2000) %>%
  group_by(gridid) %>%
  summarize(
    planet_bm_loss = mean(planet_bm_loss, na.rm = T),
    forest_share = mean(forest_share, na.rm = T)
  )

p_2000 <- ggplot(grid_plt_2000) +
  geom_sf(data = europe, fill = "grey75", color = NA) +
  geom_sf(aes(fill = planet_bm_loss), color = NA) +
  scale_fill_gradientn(
    colors = viridisLite::plasma(10),
    limits = c(0, 4),
    oob = scales::squish,
    na.value = NA
  ) +
  annotate(
    geom = 'text',
    label = '2000-2017',
    x = -Inf,
    y = Inf,
    hjust = 0,
    vjust = 5,
    size = 2.5
  ) +
  theme_void() +
  theme(
    panel.border = element_rect(color = NA, fill = NA),
    legend.position = "None",
    legend.direction = "horizontal",
    legend.spacing.y = unit(0.1, "cm"),
    legend.text = element_text(size = 6, vjust = 3),
    legend.title = element_text(size = 6, vjust = -1)
  ) +
  guides(
    color = guide_colourbar(
      theme = theme(
        legend.key.width = unit(3, "lines"),
        legend.key.height = unit(0.5, "lines")
      ),
      legend.ticks.length = unit(.05, "cm"),
      title.position = "top"
    )
  )

# 2018-2023
grid_plt_2018 <- grid_plt %>%
  filter(year >= 2018) %>%
  group_by(gridid) %>%
  summarize(
    planet_bm_loss = mean(planet_bm_loss, na.rm = T),
    forest_share = mean(forest_share, na.rm = T)
  )


p_2018 <- ggplot(grid_plt_2018) +
  geom_sf(data = europe, fill = "grey75", color = NA) +
  geom_sf(aes(fill = planet_bm_loss), color = NA) +
  scale_fill_gradientn(
    colors = viridisLite::plasma(10),
    limits = c(0, 4),
    labels = c(0, 1, 2, 3, ">4"),
    breaks = c(0, 1, 2, 3, 4),
    oob = scales::squish,
    name = expression(paste(
      "AGB loss (Mg ",
      ha[F]^{
        -1
      },
      yr.^{
        -1
      },
      ")"
    )),
    na.value = NA
  ) +
  annotate(
    geom = 'text',
    label = '2018-2023',
    x = -Inf,
    y = Inf,
    hjust = 0,
    vjust = 5,
    size = 2.5
  ) +
  theme_void() +
  theme(
    plot.margin = margin(0, 0, 0, 1, unit = "mm"),
    panel.border = element_rect(color = NA, fill = NA),
    legend.position = c(0.2, -0.01),
    legend.direction = "horizontal",
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 7, margin = margin(r = 5, unit = "mm"))
  ) +
  guides(
    fill = guide_colourbar(
      theme = theme(
        legend.key.width = unit(1, "cm"),
        legend.key.height = unit(0.15, "cm")
      ),
      legend.ticks.length = unit(.05, "cm"),
      title.position = "left"
    )
  )


# Fig. 2d -----------------------------------------------------------------

bm_loss_biome_tp <- df %>%
  left_join(agents, by = c("gridid", "year")) %>%
  group_by(year, biome_name) %>%
  summarize(
    loss_sum_harvest = sum(planet_bm_loss, na.rm = T) *
      mean(harvest / 100) /
      sum(forest_ha_sum),
    loss_sum_natural = sum(planet_bm_loss, na.rm = T) *
      mean(natural / 100) /
      sum(forest_ha_sum),
    loss_sd = sum(planet_bm_loss_sd, na.rm = T) / sum(forest_ha_sum),
    year = unique(year)
  ) %>%
  ungroup() %>%
  mutate(
    time_period = ifelse(
      year < 2000,
      "1985-1999",
      ifelse(year < 2018, "2000-2017", "2018-2023")
    )
  ) %>%
  group_by(time_period, biome_name) %>%
  summarize(
    loss_sum_harvest = mean(loss_sum_harvest),
    loss_sum_natural = mean(loss_sum_natural),
    loss_sum_total = loss_sum_harvest + loss_sum_natural,
    loss_sd = mean(loss_sd)
  )


top_labels <- data.frame(
  "biome_name" = factor(
    c("Temperate", "Boreal", "Mediterranean"),
    levels = c("Temperate", "Boreal", "Mediterranean")
  ),
  "loss_sum_total" = rep(1.4, 3),
  "time_period" = rep("2000-2017", 3)
)

p_biomes_tp <- bm_loss_biome_tp %>%
  mutate(
    biome_name = factor(
      biome_name,
      levels = c("Temperate", "Boreal", "Mediterranean")
    )
  ) %>%
  pivot_longer(-c(time_period, biome_name, loss_sum_total, loss_sd)) %>%
  mutate(
    name = factor(
      name,
      levels = c("loss_sum_harvest", "loss_sum_natural"),
      labels = c("Harvest", "Natural")
    )
  ) %>%
  ggplot() +
  geom_bar(aes(x = time_period, y = value, fill = name), stat = "identity") +
  geom_errorbar(
    aes(
      x = time_period,
      ymin = (loss_sum_total - loss_sd),
      ymax = (loss_sum_total + loss_sd)
    ),
    width = 0.2,
    linewidth = 0.1
  ) +
  geom_text(
    data = top_labels,
    aes(x = time_period, y = loss_sum_total, label = biome_name),
    size = 2.5
  ) +
  scale_fill_manual(values = c("grey70", "grey30")) +
  coord_cartesian(clip = "off") +
  facet_wrap(~biome_name) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = c(0.85, 0.7),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 7),
    axis.text.y = element_text(size = 7, color = "black"),
    axis.text.x = element_text(
      size = 7,
      angle = 50,
      hjust = 1,
      color = "black"
    ),
    strip.background = element_rect(fill = NA, color = NA),
    strip.text = element_blank(),
    axis.ticks = element_line(linewidth = 0.1),
    panel.border = element_blank(),
    legend.text = element_text(size = 7),
    legend.title = element_blank(),
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
  ) +
  labs(
    y = expression(paste(
      "AGB loss (Mg ",
      ha[F]^{
        -1
      },
      "yr."^{
        -1
      },
      ")"
    ))
  )


# combine plots
p_hist <- ggarrange(
  p_biomes,
  p_biomes_tp,
  ncol = 1,
  nrow = 2,
  labels = c("b", "d"),
  font.label = list(size = 8, face = "bold"),
  label.y = c(1, 1.05)
)

p_map <- ggarrange(p_1985, p_2000, p_2018, ncol = 3, nrow = 1)

p_ts <- ggarrange(
  p_loss_biome_ts,
  p_map,
  ncol = 1,
  nrow = 2,
  heights = c(0.8, 1),
  labels = c("a", "c"),
  font.label = list(size = 8, face = "bold")
)

p_fig_2 <- ggarrange(p_ts, p_hist, ncol = 2, nrow = 1, widths = c(2, 1))

ragg::agg_png(
  "figures/fig_2_agb_loss_time_series_maps.jpg",
  width = 18,
  height = 9,
  units = "cm",
  res = 500
)
p_fig_2
dev.off()


# Tab. S1. Country summary for severity of AGBD loss and total AGB --------

tab_s1 <- df %>%
  dplyr::group_by(country) %>%
  summarize(
    agbd_before_mean = round(mean(planet_agbd_before_mean, na.rm = T), 1),
    agbd_before_sd = round(sd(planet_agbd_before_mean, na.rm = T), 1),

    agbd_loss_mean = round(mean(planet_mean_sev_abs, na.rm = T), 1),
    agbd_loss_sd = round(sd(planet_mean_sev_abs, na.rm = T), 1),

    agbd_loss_rel_mean = round(mean(planet_mean_sev_rel * 100, na.rm = T), 1),
    agbd_loss_rel_sd = round(sd(planet_mean_sev_rel * 100, na.rm = T), 1),

    agb_loss_total = round(sum(planet_bm_loss, na.rm = T) * 1e-09, 3),
    agb_loss_total_sd = round(sum(planet_bm_loss_sd, na.rm = T) * 1e-09, 3),

    agb_loss_haf = round(
      sum(planet_bm_loss, na.rm = T) / sum(forest_ha_sum),
      1
    ),
    agb_loss_haf_sd = round(
      sum(planet_bm_loss_sd, na.rm = T) / sum(forest_ha_sum),
      1
    ),

    country_iso = unique(country_iso)
  ) %>%
  filter(country != "") %>%
  mutate(
    country = ifelse(
      country == "U.K. of Great Britain and Northern Ireland",
      "United Kingdom",
      country
    ),
    country = ifelse(
      country == "The former Yugoslav Republic of Macedonia",
      "North Macedonia",
      country
    ),
    country = ifelse(country == "Moldova, Republic of", "Moldova", country)
  ) %>%
  mutate(
    agbd_before_mean = paste(agbd_before_mean, "\u00B1", agbd_before_sd),
    agbd_loss_mean = paste(agbd_loss_mean, "\u00B1", agbd_loss_sd),
    agbd_loss_rel_mean = paste(agbd_loss_rel_mean, "\u00B1", agbd_loss_rel_sd),
    agb_loss_total = paste(agb_loss_total, "\u00B1", agb_loss_total_sd),
    agb_loss_haf = paste(agb_loss_haf, "\u00B1", agb_loss_haf_sd),
  ) %>%
  filter(
    country != "Andorra"
  ) %>%
  select(
    country,
    country_iso,
    agbd_before_mean,
    agbd_loss_mean,
    agbd_loss_rel_mean,
    agb_loss_total,
    agb_loss_haf
  )

tab_s1 <- tab_s1 %>% arrange(country)

readr::write_excel_csv(
  tab_s1,
  "tables/tab_s1_country_summary.csv",
  quote = "none"
)


# Tab. S2. Percent change AGB loss from 2018-2023 and 1985-2017 -----------

bm_loss_biome_tp <- df %>%
  left_join(agents, by = c("gridid", "year")) %>%
  group_by(year, biome_name) %>%
  summarize(
    loss_sum_harvest = sum(planet_bm_loss, na.rm = T) *
      mean(harvest / 100) /
      sum(forest_ha_sum),
    loss_sum_natural = sum(planet_bm_loss, na.rm = T) *
      mean(natural / 100) /
      sum(forest_ha_sum),
    loss_sd = sum(planet_bm_loss_sd, na.rm = T) / sum(forest_ha_sum),
    year = unique(year)
  ) %>%
  ungroup() %>%
  mutate(time_period = ifelse(year < 2017, "1985-2017", "2018-2023")) %>%
  group_by(time_period, biome_name) %>%
  summarize(
    loss_sum_harvest = mean(loss_sum_harvest),
    loss_sum_natural = mean(loss_sum_natural),
    loss_sum_total = loss_sum_harvest + loss_sum_natural
  )

bm_loss_change_biome_tp <- bm_loss_biome_tp %>%
  pivot_wider(
    names_from = time_period,
    values_from = c(loss_sum_harvest, loss_sum_natural, loss_sum_total)
  ) %>%
  mutate(
    percent_change_harvest = (`loss_sum_harvest_2018-2023` -
      `loss_sum_harvest_1985-2017`) /
      `loss_sum_harvest_1985-2017` *
      100,
    percent_change_natural = (`loss_sum_natural_2018-2023` -
      `loss_sum_natural_1985-2017`) /
      `loss_sum_natural_1985-2017` *
      100,
    percent_change_total = (`loss_sum_total_2018-2023` -
      `loss_sum_total_1985-2017`) /
      `loss_sum_total_1985-2017` *
      100
  )

bm_loss_tp <- df %>%
  left_join(agents, by = c("gridid", "year")) %>%
  group_by(year) %>%
  summarize(
    loss_sum_harvest = sum(planet_bm_loss, na.rm = T) *
      mean(harvest / 100) /
      sum(forest_ha_sum),
    loss_sum_natural = sum(planet_bm_loss, na.rm = T) *
      mean(natural / 100) /
      sum(forest_ha_sum),
    loss_sd = sum(planet_bm_loss_sd, na.rm = T) / sum(forest_ha_sum),
    year = unique(year)
  ) %>%
  ungroup() %>%
  mutate(time_period = ifelse(year < 2017, "1985-2017", "2018-2023")) %>%
  group_by(time_period) %>%
  summarize(
    loss_sum_harvest = mean(loss_sum_harvest),
    loss_sum_natural = mean(loss_sum_natural),
    loss_sum_total = loss_sum_harvest + loss_sum_natural
  )

bm_loss_change_tp <- bm_loss_tp %>%
  pivot_wider(
    names_from = time_period,
    values_from = c(loss_sum_harvest, loss_sum_natural, loss_sum_total)
  ) %>%
  mutate(
    percent_change_harvest = (`loss_sum_harvest_2018-2023` -
      `loss_sum_harvest_1985-2017`) /
      `loss_sum_harvest_1985-2017` *
      100,
    percent_change_natural = (`loss_sum_natural_2018-2023` -
      `loss_sum_natural_1985-2017`) /
      `loss_sum_natural_1985-2017` *
      100,
    percent_change_total = (`loss_sum_total_2018-2023` -
      `loss_sum_total_1985-2017`) /
      `loss_sum_total_1985-2017` *
      100,
    biome_name = "Europe"
  )


bm_loss_change_eur_tp <- rbind(bm_loss_change_biome_tp, bm_loss_change_tp)

tab_s2 <- bm_loss_change_eur_tp %>%
  dplyr::select(
    biome_name,
    percent_change_harvest,
    percent_change_natural,
    percent_change_total
  ) %>%
  mutate(across(2:4, \(x) round(x, 1))) %>%
  arrange(match(
    biome_name,
    c("Boreal", "Temperate", "Mediterranean", "Europe")
  ))


readr::write_excel_csv(
  tab_s2,
  "tables/tab_s2_bm_loss_change_tp.csv",
  quote = "none"
)

# Fig. S4 -----------------------------------------------------------------

bm_loss_agent_biome <- df %>%
  left_join(agents, by = c("gridid", "year")) %>%
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

bm_area_loss_natural <- bm_loss_agent_biome %>%
  group_by(biome_name) %>%
  summarize(
    bm_loss_natural = sum(bm_loss_natural, na.rm = T),
    dist_area_natural = sum(dist_area_natural, na.rm = T)
  ) %>%
  mutate(
    area_total_natural = sum(dist_area_natural),
    bm_loss_total_natural = sum(bm_loss_natural),
    area_share_natural = dist_area_natural / area_total_natural,
    bm_share_natural = bm_loss_natural / bm_loss_total_natural,
    diff_area_bm = round((bm_share_natural - area_share_natural) * 100, 2),
    diff_area_bm = ifelse(
      diff_area_bm > 0,
      paste("+", as.character(diff_area_bm), sep = ""),
      as.character(diff_area_bm)
    )
  ) %>%
  ungroup() %>%
  dplyr::select(
    biome_name,
    area_share_natural,
    bm_share_natural,
    diff_area_bm
  ) %>%
  mutate(
    ratios = paste(
      "1:",
      round(bm_share_natural / area_share_natural, 2),
      sep = ""
    )
  ) %>%
  pivot_longer(-c(biome_name, diff_area_bm, ratios)) %>%
  mutate(
    biome_name = factor(
      biome_name,
      levels = c("Boreal", "Temperate", "Mediterranean")
    ),
    name = ifelse(name == "area_share_natural", "Area", "AGB loss"),
    name = factor(name, levels = c("Area", "AGB loss")),
    agent = "Natural disturbances"
  )


# harvest
bm_area_loss_harvest <- bm_loss_agent_biome %>%
  group_by(biome_name) %>%
  summarize(
    bm_loss_harvest = sum(bm_loss_harvest, na.rm = T),
    dist_area_harvest = sum(dist_area_harvest, na.rm = T)
  ) %>%
  mutate(
    area_total_harvest = sum(dist_area_harvest),

    bm_loss_total_harvest = sum(bm_loss_harvest),

    area_share_harvest = dist_area_harvest / area_total_harvest,
    bm_share_harvest = bm_loss_harvest / bm_loss_total_harvest,

    diff_area_bm = round((bm_share_harvest - area_share_harvest) * 100, 2),
    diff_area_bm = ifelse(
      diff_area_bm > 0,
      paste("+", as.character(diff_area_bm), sep = ""),
      as.character(diff_area_bm)
    )
  ) %>%
  ungroup() %>%
  dplyr::select(
    biome_name,
    area_share_harvest,
    bm_share_harvest,
    diff_area_bm
  ) %>%
  mutate(
    ratios = paste(
      "1:",
      round(bm_share_harvest / area_share_harvest, 2),
      sep = ""
    )
  ) %>%
  pivot_longer(-c(biome_name, diff_area_bm, ratios)) %>%
  mutate(
    biome_name = factor(
      biome_name,
      levels = c("Boreal", "Temperate", "Mediterranean")
    ),
    name = ifelse(name == "area_share_harvest", "Area", "AGB loss"),
    name = factor(name, levels = c("Area", "AGB loss")),
    agent = "Harvest"
  )

df_plt <- rbind(bm_area_loss_harvest, bm_area_loss_natural)

df_plt <- df_plt %>%
  mutate(agent = factor(agent, levels = c("Harvest", "Natural disturbances")))

ratio_labs <- df_plt %>%
  dplyr::select(-c(value, name, diff_area_bm)) %>%
  group_by(biome_name, agent) %>%
  slice(1)

p_fig_s4 <-
  ggplot(
    df_plt,
    aes(x = name, alluvium = biome_name, y = value, stratum = biome_name)
  ) +
  geom_stratum(
    aes(fill = biome_name, colour = biome_name),
    linewidth = 0.2,
    alpha = .6,
    width = 0.2
  ) +
  geom_alluvium(
    aes(fill = biome_name, color = biome_name),
    linewidth = 0.2,
    curve_type = "cubic",
    width = 0.2,
    alpha = .2
  ) +
  scale_fill_manual(values = c("#1192e8", "#6929c4", "#ee538b"), name = "") +
  scale_color_manual(values = c("#1192e8", "#6929c4", "#ee538b", name = "")) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = c(0, 1),
    labels = c("0%", "100%")
  ) +
  scale_x_discrete(expand = c(0.55, 0.05)) +
  facet_rep_wrap(
    ~agent,
    scales = "free_x",
    nrow = 2,
    ncol = 1,
    labeller = label_wrap_gen(multi_line = TRUE),
    repeat.tick.labels = F
  ) +
  geom_text_repel(
    aes(
      label = ifelse(name == "Area", round(value * 100, 1), NA),
      color = biome_name
    ),
    stat = "stratum",
    size = 2,
    direction = "y",
    nudge_x = -.3,
    min.segment.length = unit(5, "cm"),
    hjust = 1,
  ) +
  geom_text_repel(
    aes(
      label = ifelse(name == "AGB loss", round(value * 100, 1), NA),
      color = biome_name
    ),
    stat = "stratum",
    size = 2,
    direction = "y",
    nudge_x = .3,
    min.segment.length = unit(5, "cm")
  ) +
  geom_text_repel(
    aes(label = ifelse(name == "AGB loss", ratios, NA), color = biome_name),
    stat = "flow",
    size = 2,
    nudge_x = -.5,
    min.segment.length = unit(5, "cm")
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.direction = "vertical",
    legend.justification = "left",
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text = element_text(size = 7, color = "black"),
    axis.ticks.x = element_line(colour = "black", linewidth = 0.2),
    axis.ticks.y = element_line(linewidth = 0.1),
    axis.ticks.length = unit(.05, "cm"),
    plot.title = element_text(hjust = 0.5, size = 7),
    legend.text = element_text(
      size = 6,
      margin = margin(r = 0, l = 0.5, unit = "mm")
    ),
    legend.key.width = unit(2, "mm"),
    legend.key.height = unit(1, "mm"),
    legend.spacing.x = unit(-5, "mm"),
    legend.margin = margin(-10, 0, 0, -5, unit = "mm"),
    panel.border = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 7, margin = margin(b = 0))
  ) +
  annotate("segment", x = -Inf, xend = -Inf, y = 0, yend = 1, linewidth = 0.2) +
  guides(color = "none")

ragg::agg_png(
  "figures/fig_s4_agb_loss_area_relation.jpg",
  width = 80,
  height = 120,
  units = "mm",
  res = 500
)
p_fig_s4
dev.off()


# Fig. S9 Corrected annual area of tree cover loss ------------------------

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

p_fig_s9 <- dist_area_orig %>%
  dplyr::select(gridid, year, dist_ha_sum) %>%
  rename("dist_ha_sum_orig" = "dist_ha_sum") %>%
  left_join(dist_area) %>%
  group_by(year) %>%
  summarize(
    dist_ha_mean = sum(dist_ha_mean, na.rm = T),
    dist_ha_sd = sum(dist_ha_sd, na.rm = T),
    dist_ha_sum_orig = sum(dist_ha_sum_orig, na.rm = T)
  ) %>%
  pivot_longer(-c(year, dist_ha_sd)) %>%
  ggplot() +
  geom_bar(
    aes(x = year, y = value * 0.01, fill = name),
    stat = "identity",
    position = "dodge"
  ) +
  scale_fill_manual(
    values = c("grey40", "grey85"),
    labels = c("Corrected area estimate", "Map area"),
    name = ""
  ) +
  scale_x_continuous(
    minor_breaks = seq(1985, 2023, by = 1),
    breaks = seq(1985, 2023, 2),
    labels = seq(1985, 2023, 2)
  ) +
  guides(
    x = guide_axis(minor.ticks = TRUE)
  ) +
  theme_bw() +
  theme(
    legend.key.size = unit(4, "mm"),
    panel.grid = element_blank(),
    axis.text = element_text(
      size = 7,
      color = "black"
    ),
    axis.title = element_text(size = 8),
    legend.position = "bottom",
    legend.text = element_text(size = 7),
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
    )
  ) +
  labs(
    x = "Year",
    y = expression(paste(
      "Tree cover loss (",
      km^{
        2
      },
      ")"
    ))
  )

ragg::agg_png(
  "figures/fig_s9_area_tree_cover_loss.jpg",
  width = 180,
  height = 80,
  units = "mm",
  res = 500
)
p_fig_s9
dev.off()
