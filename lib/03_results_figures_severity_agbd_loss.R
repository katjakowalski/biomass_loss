# ---
# Title:        Figures for severity of AGBD loss
# Description:  Creates and saves figures 1, S1, S2, and S6
# Author:       Katja Kowalski
# Affiliation:  Technical University of Munich, School of Life Sciences
# Date Created: 2025-08-01
# Last Updated: 2025-08-04
# ---

library(ggplot2)
library(ggrepel)
library(ggpubr)
library(patchwork)
library(ggpmisc)
library(lemon)
library(sf)
library(tidyverse)
library(data.table)
library(ragg)

# Prepare data ------------------------------------------------------------

df_agbd <- read.csv("data/04_aggregation_planet_agbd/planet_agbd_loss_20km.csv")

grid <- read_sf(
  "data/gis/referencegrid/grid_20km_3035.shp"
)

# add forest area
forest <- read.csv("data/03_aggregation_area/forest_area_20km.csv")
europe <- read_sf("data/gis/borders/europe_outline_3035.shp")

df_agbd <- df_agbd %>%
  left_join(forest) %>%
  mutate(
    forest_ha_sum = forest_px_sum * 0.09,
    forest_share = forest_ha_sum / 40000
  ) %>%
  filter(forest_share >= 0.1)


# add ecoregions
ecoregions <- read_sf(
  "data/gis/ecoregions/terrestrial_ecoregions_olson_europe_3035.shp"
)

eco_intersection <- st_join(
  grid,
  ecoregions,
  join = st_intersects,
  largest = TRUE
)

# aggregate ecoregions
eco_intersection_europe <- st_intersection(europe, ecoregions)

eco_intersection_europe <- eco_intersection_europe %>%
  filter(BIOME == 4 | BIOME == 5 | BIOME == 6 | BIOME == 12) %>%
  mutate(
    biome_name = ifelse(BIOME == 4, "Temperate", NA),
    biome_name = ifelse(BIOME == 12, "Mediterranean", biome_name),
    biome_name = ifelse(BIOME == 5, "Temperate", biome_name),
    biome_name = ifelse(BIOME == 6, "Boreal", biome_name)
  ) %>%
  group_by(biome_name) %>%
  summarize()

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
  dplyr::select(gridid, biome_name)

# filter ecoregions
df_agbd <- df_agbd %>%
  left_join(eco_intersection) %>%
  filter(biome_name != "Tundra" & biome_name != "Pontic Steppe")

# add countries
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
    country = unique(isocode)[1]
  )

df_agbd <- df_agbd %>%
  left_join(country_intersection)

# add agents

agents <- fread("data/03_aggregation_area/disturbance_agents_20km.csv")

agents <- agents %>%
  filter(
    year >= 2014
  ) %>%
  mutate(
    agent = ifelse(agent == 3, "harvest", "natural") # 1 bb/wind, 2 fire, 3 harvest
  ) %>%
  group_by(gridid, agent) %>%
  summarize(dist_px_sum_gridcell = sum(dist_px, na.rm = T)) %>%
  group_by(gridid) %>%
  mutate(
    dist_percent = dist_px_sum_gridcell * 100 / sum(dist_px_sum_gridcell)
  ) %>%
  ungroup() %>%
  mutate(dist_percent = ifelse(is.na(dist_percent), 0, dist_percent)) %>%
  dplyr::select(-dist_px_sum_gridcell) %>%
  pivot_wider(names_from = agent, values_from = dist_percent) %>%
  mutate(
    agent_type = ifelse(natural >= 50, "natural", "harvest")
  ) %>%
  dplyr::select(
    gridid,
    agent_type
  )

df_agbd <- df_agbd %>%
  left_join(agents)

# add grid
grid_agbd <- grid %>%
  st_centroid() %>%
  left_join(
    df_agbd,
    by = "gridid"
  ) %>%
  dplyr::select(
    -c(
      planet_severity_abs_q05,
      planet_severity_abs_mean,
      planet_severity_rel_q05,
      planet_severity_rel_mean
    )
  ) %>%
  group_by(gridid) %>%
  mutate(
    across(starts_with("planet_severity"), ~ ifelse(. < 0, 0, .)),
    across(starts_with("planet_severity_rel"), ~ ifelse(. > 1, 1, .))
  ) %>%
  mutate(
    planet_mean_sev_abs = mean(
      c_across(starts_with("planet_severity_abs")),
      na.rm = TRUE
    ),
    planet_mean_sev_rel = mean(
      c_across(starts_with("planet_severity_rel")),
      na.rm = T
    )
  )

df_agbd <- df_agbd %>%
  dplyr::select(
    -c(
      planet_severity_abs_q05,
      planet_severity_abs_mean,
      planet_severity_rel_q05,
      planet_severity_rel_mean
    )
  ) %>%
  group_by(gridid) %>%
  mutate(
    across(starts_with("planet_severity"), ~ ifelse(. < 0, 0, .)),
    across(starts_with("planet_severity_rel"), ~ ifelse(. > 1, 1, .))
  ) %>%
  mutate(
    planet_mean_sev_abs = mean(
      c_across(starts_with("planet_severity_abs")),
      na.rm = TRUE
    ),
    planet_mean_sev_rel = mean(
      c_across(starts_with("planet_severity_rel")),
      na.rm = T
    )
  )


# Fig. 1. Average severity of AGBD losses due to tree cover loss. ---------

xlim <- c(2750000, 6350000)
ylim <- c(1600000, 5250000)


p1 <- grid_agbd %>%
  filter(!is.na(planet_mean_sev_abs)) %>%
  ggplot() +
  geom_sf(data = europe, fill = "grey85", color = NA) +
  geom_sf(aes(col = planet_mean_sev_abs, size = forest_share), stroke = 0) +
  geom_sf(
    data = eco_intersection_europe,
    color = "grey30",
    fill = NA,
    linewidth = 0.1
  ) +
  scale_size("size_area", range = c(0.2, 0.7), guide = 'none') +
  coord_sf(xlim = xlim, ylim = ylim) +
  scale_color_gradientn(
    colours = viridisLite::plasma(10),
    limits = c(0, 200),
    breaks = c(0, 100, 200),
    labels = c(0, 100, ">200"),
    oob = scales::squish,
    name = expression(paste(
      "Absolute severity (Mg ",
      ha[Dist]^{
        -1
      },
      ")"
    )),
    na.value = NA
  ) +
  theme_void() +
  theme(
    panel.border = element_rect(color = NA, fill = NA),
    legend.position = c(0.26, 0.94), #c(0.2,0.8),
    plot.margin = margin(0, 0, 0, -3, unit = "mm"),
    legend.direction = "horizontal",
    legend.text = element_text(size = 7, margin = margin(t = 0.5, unit = "mm")),
    legend.title = element_text(size = 7, margin = margin(b = 1, unit = "mm"))
  ) +
  guides(
    color = guide_colourbar(
      theme = theme(
        legend.key.width = unit(3, "lines"),
        legend.key.height = unit(0.4, "lines")
      ),
      title.position = "top",
      legend.ticks.length = unit(.05, "cm")
    )
  )


p2 <- grid_agbd %>%
  filter(!is.na(planet_mean_sev_rel)) %>%
  ggplot() +
  geom_sf(data = europe, fill = "grey85", color = NA) +
  geom_sf(aes(col = planet_mean_sev_rel, size = forest_share), stroke = 0) +
  geom_sf(
    data = eco_intersection_europe,
    color = "grey10",
    fill = NA,
    linewidth = 0.05
  ) +
  scale_size("size_area", range = c(0.2, 0.7), guide = 'none') +
  coord_sf(xlim = xlim, ylim = ylim) +
  scale_color_gradientn(
    colours = viridisLite::plasma(10),
    limits = c(0, 1),
    breaks = c(0, 0.5, 1),
    labels = c(0, 50, 100),
    oob = scales::squish,
    name = "Relative severity (%)",
    na.value = NA
  ) +
  theme_void() +
  theme(
    panel.border = element_rect(color = NA, fill = NA),
    legend.position = c(0.28, 0.94), #c(0.2,0.8),
    plot.margin = margin(0, 0, 0, -3, unit = "mm"),
    legend.direction = "horizontal",
    legend.text = element_text(size = 7, margin = margin(t = 0.5, unit = "mm")),
    legend.title = element_text(size = 7, margin = margin(b = 1, unit = "mm"))
  ) +
  guides(
    color = guide_colourbar(
      theme = theme(
        legend.key.width = unit(3, "lines"),
        legend.key.height = unit(0.4, "lines")
      ),
      title.position = "top",
      legend.ticks.length = unit(.05, "cm")
    )
  )

p_fig_1 <- ggarrange(
  p1,
  p2,
  labels = c("a", "b"),
  font.label = list(size = 7, face = "bold")
)

ragg::agg_png(
  "figures/fig_1_abs_rel_severity.jpg",
  width = 180,
  height = 90,
  units = "mm",
  res = 500
)
p_fig_1
dev.off()


# Fig. S1. Distribution of absolute and relative severity of AGBD  --------

grid_agbd_plt <- grid_agbd %>%
  dplyr::select(
    planet_mean_sev_abs,
    planet_mean_sev_rel,
    agent_type,
    forest_share,
    gridid,
    biome_name
  ) %>%
  pivot_longer(-c(gridid, agent_type, forest_share, geometry, biome_name))


p_boreal_abs <- grid_agbd_plt %>%
  filter(name == "planet_mean_sev_abs" & biome_name == "Boreal") %>%
  ggplot() +
  geom_histogram(aes(x = value), fill = "#648fff", alpha = 0.7) + #
  scale_x_continuous(
    limits = c(-10, 220),
    labels = c(0, 100, 200),
    breaks = c(0, 100, 200)
  ) +
  facet_wrap(~biome_name, scales = "free_y") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 7),
    axis.text.y = element_text(size = 7, color = "black"),
    axis.text.x = element_blank(),
    legend.position = "None",
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
    axis.ticks = element_line(colour = "black", linewidth = 0.2),
    axis.ticks.length = unit(.05, "cm"),
    strip.background = element_rect(fill = NA, color = NA),
    strip.text = element_text(size = 7),
    plot.margin = unit(c(0, 1, 0, 0.5), "mm")
  )


p_temperate_abs <- grid_agbd_plt %>%
  filter(name == "planet_mean_sev_abs" & biome_name == "Temperate") %>%
  ggplot() +
  geom_histogram(aes(x = value), fill = "#785ef0", alpha = 0.7) +
  scale_x_continuous(
    limits = c(-10, 220),
    labels = c(0, 100, 200),
    breaks = c(0, 100, 200)
  ) +
  facet_wrap(~biome_name, scales = "free_y") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.title.y = element_text(size = 7),
    axis.text.y = element_text(size = 7, color = "black"),
    axis.text.x = element_blank(),
    legend.position = "None",
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
    axis.ticks = element_line(colour = "black", linewidth = 0.2),
    axis.ticks.length = unit(.05, "cm"),
    strip.background = element_rect(fill = NA, color = NA),
    strip.text = element_text(size = 7),
    plot.margin = unit(c(0, 1, 0, 0.5), "mm")
  )


p_med_abs <- grid_agbd_plt %>%
  filter(name == "planet_mean_sev_abs" & biome_name == "Mediterranean") %>%
  ggplot() +
  geom_histogram(aes(x = value), fill = "#dc267f", alpha = 0.7) +
  scale_x_continuous(
    limits = c(-10, 220),
    labels = c(0, 100, 200),
    breaks = c(0, 100, 200)
  ) +
  scale_y_continuous(labels = seq(0, 300, 150), breaks = seq(0, 300, 150)) +
  facet_wrap(~biome_name, scales = "free_y") +
  theme_bw() +
  labs(
    x = expression(paste(
      "Absolute severity (Mg ",
      ha[Dist]^{
        -1
      },
      ")"
    ))
  ) +
  theme(
    panel.grid = element_blank(),
    #axis.title = element_blank(),
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 7, color = "black"),
    legend.position = "None",
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
    axis.ticks = element_line(colour = "black", linewidth = 0.2),
    axis.ticks.length = unit(.05, "cm"),
    strip.background = element_rect(fill = NA, color = NA),
    strip.text = element_text(size = 7),
    plot.margin = unit(c(0, 1, 0, 0.5), "mm")
  )


p_boreal_rel <- grid_agbd_plt %>%
  filter(name == "planet_mean_sev_rel" & biome_name == "Boreal") %>%
  ggplot() +
  geom_histogram(aes(x = value), fill = "#648fff", alpha = 0.7) +
  scale_x_continuous(
    limits = c(-0.02, 1.05),
    breaks = c(0, 0.5, 1),
    labels = c(0, 50, 100)
  ) +
  scale_y_continuous(labels = seq(0, 200, 100), breaks = seq(0, 200, 100)) +
  facet_wrap(~biome_name, scales = "free_y") + # labeller = as_labeller(strip_names))+
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 7),
    axis.text.y = element_text(size = 7, color = "black"),
    axis.text.x = element_blank(),
    legend.position = "None",
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
    axis.ticks = element_line(colour = "black", linewidth = 0.2),
    axis.ticks.length = unit(.05, "cm"),
    strip.background = element_rect(fill = NA, color = NA),
    strip.text = element_text(size = 7),
    plot.margin = unit(c(0, 0.5, 0, 0), "mm")
  )


p_temperate_rel <- grid_agbd_plt %>%
  filter(
    name == "planet_mean_sev_rel" &
      biome_name == "Temperate"
  ) %>%
  ggplot() +
  geom_histogram(aes(x = value), fill = "#785ef0", alpha = 0.7) +
  scale_x_continuous(
    limits = c(-0.02, 1.05),
    breaks = c(0, 0.5, 1),
    labels = c(0, 50, 100)
  ) +
  scale_y_continuous(labels = seq(0, 800, 400), breaks = seq(0, 800, 400)) +
  facet_wrap(~biome_name, scales = "free_y") + # labeller = as_labeller(strip_names))+
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 7),
    axis.text.y = element_text(size = 7, color = "black"),
    axis.text.x = element_blank(),
    legend.position = "None",
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
    axis.ticks = element_line(colour = "black", linewidth = 0.2),
    axis.ticks.length = unit(.05, "cm"),
    strip.background = element_rect(fill = NA, color = NA),
    strip.text = element_text(size = 7),
    plot.margin = unit(c(0, 0.5, 0, 0), "mm")
  )


p_med_rel <- grid_agbd_plt %>%
  filter(name == "planet_mean_sev_rel" & biome_name == "Mediterranean") %>%
  ggplot() +
  geom_histogram(aes(x = value), fill = "#dc267f", alpha = 0.7) +
  scale_x_continuous(
    limits = c(-0.02, 1.05),
    breaks = c(0, 0.5, 1),
    labels = c(0, 50, 100)
  ) +
  scale_y_continuous(labels = c(0, 50, 100), breaks = c(0, 50, 100)) +
  facet_wrap(~biome_name, scales = "free_y") +
  labs(x = "Relative severity (%)") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 7, color = "black"),
    legend.position = "None",
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
    axis.ticks = element_line(colour = "black", linewidth = 0.2),
    axis.ticks.length = unit(.05, "cm"),
    strip.background = element_rect(fill = NA, color = NA),
    strip.text = element_text(size = 7),
    plot.margin = unit(c(0, 0.5, 0, 0), "mm")
  )


p_hist_abs <- ggarrange(
  p_boreal_abs,
  p_temperate_abs,
  p_med_abs,
  ncol = 1,
  nrow = 3
)
p_hist_rel <- ggarrange(
  p_boreal_rel,
  p_temperate_rel,
  p_med_rel,
  ncol = 1,
  nrow = 3
)
p_fig_s1 <- ggarrange(
  p_hist_abs,
  p_hist_rel,
  ncol = 2,
  labels = c("a", "b"),
  font.label = list(size = 7, face = "bold")
)

ragg::agg_png(
  "figures/fig_s1_abs_rel_histograms.jpg",
  width = 120,
  height = 100,
  units = "mm",
  res = 500
)

p_fig_s1
dev.off()


# Fig. S2. Average absolute and relative severity of AGBD loss ------------

df_country <- grid_agbd %>%
  drop_na(country) %>%
  group_by(country) %>%
  summarise(
    abs = mean(planet_mean_sev_abs, na.rm = T),
    rel = mean(planet_mean_sev_rel, na.rm = T),
    abs_sd = sd(planet_mean_sev_abs, na.rm = T),
    rel_sd = sd(planet_mean_sev_rel, na.rm = T),
    pre_dist = mean(planet_agbd_before_mean, na.rm = T),
    pre_dist_sd = sd(planet_agbd_before_mean, na.rm = T),
    biome_name = paste(unique(biome_name), collapse = "_"),
    .groups = 'drop'
  ) %>%
  mutate(
    biome_name = ifelse(
      grepl("Temperate_Mediterranean|Mediterranean_Temperate", biome_name),
      'Mediterranean/Temperate',
      biome_name
    ),
    biome_name = ifelse(
      grepl("Temperate_Boreal|Boreal_Temperate", biome_name),
      'Boreal/Temperate',
      biome_name
    ),
    biome_name = factor(
      biome_name,
      levels = c('Boreal/Temperate', 'Temperate', 'Mediterranean/Temperate')
    )
  ) %>%
  filter(
    country != "AD"
  )

df_country$rel <- df_country$rel * 100
df_country$rel_sd <- df_country$rel_sd * 100

df_country_plt <- df_country %>%
  sf::st_drop_geometry() %>%
  dplyr::select(country, biome_name, abs, rel, abs_sd, rel_sd) %>%
  pivot_longer(-c(country, biome_name)) %>%
  mutate(
    type = ifelse(name == "abs" | name == "rel", "val", "sd"),
    name = str_split_i(name, "_", 1)
  ) %>%
  pivot_wider(names_from = type, values_from = value)

p_country_abs <- ggplot(
  df_country,
  aes(x = reorder(country, -abs), y = abs, col = biome_name)
) +
  geom_point(stat = "identity") +
  geom_errorbar(
    aes(ymin = abs - abs_sd, ymax = abs + abs_sd),
    linewidth = 0.3,
    width = 0.3,
    alpha = .5
  ) +
  scale_color_manual(values = c("#648fff", "#785ef0", "#dc267f"), name = "") +
  scale_y_continuous(limits = c(0, 175)) +
  labs(y = paste("Absolute\nseverity")) + #(Mg ", ha[Dist]^{-1}, ")")))+
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(linewidth = 0.1, color = "grey70"),
    panel.grid.minor.y = element_blank(),
    axis.title.y = element_text(size = 7, color = "black"),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 7, color = "black"),
    legend.position = "None",
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
    axis.ticks = element_line(colour = "black", linewidth = 0.2),
    axis.ticks.length = unit(.05, "cm"),
    plot.background = element_rect(fill = "transparent", colour = NA)
  )

p_country_rel <- ggplot(
  df_country,
  aes(x = reorder(country, -rel), y = rel, col = biome_name)
) +
  geom_point(stat = "identity") +
  geom_errorbar(
    aes(ymin = rel - rel_sd, ymax = rel + rel_sd),
    linewidth = 0.3,
    width = 0.3,
    alpha = .5
  ) +
  scale_color_manual(values = c("#648fff", "#785ef0", "#dc267f"), name = "") +
  scale_y_continuous(limits = c(0, 100)) +
  labs(y = paste("Relative\nseverity")) + #(%)")+
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(linewidth = 0.1, color = "grey70"),
    panel.grid.minor.y = element_blank(),
    axis.title.y = element_text(size = 7, color = "black"),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 7, color = "black"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.key.spacing.x = unit(5, "mm"),
    panel.border = element_blank(),
    legend.background = element_rect(fill = NA),
    legend.key.height = unit(0.1, "cm"),
    legend.key.width = unit(0.1, "cm"),
    legend.text = element_text(size = 7),
    legend.justification = "left",
    legend.margin = margin(-2, 0, 0, 0, unit = "mm"),
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
    axis.ticks = element_line(colour = "black", linewidth = 0.2),
    axis.ticks.length = unit(.05, "cm"),
    plot.background = element_rect(fill = "transparent", colour = NA)
  )


p_fig_s2 <- ggarrange(
  p_country_abs,
  p_country_rel,
  ncol = 1,
  labels = c("a", "b"),
  font.label = list(size = 7, face = "bold")
)

ragg::agg_png(
  "figures/fig_s2_abs_rel_country.jpg",
  width = 180,
  height = 90,
  units = "mm",
  res = 500
)

p_fig_s2
dev.off()


# Fig. S6. Ecoregions of Europe. ------------------------------------------

xlim <- c(2750000, 6350000)
ylim <- c(1600000, 5250000)

p_fig_s6 <- eco_intersection_europe %>%
  mutate(
    biome_name = factor(
      biome_name,
      levels = c("Boreal", "Temperate", "Mediterranean")
    )
  ) %>%
  ggplot() +
  geom_sf(data = europe, color = NA, fill = "grey85") +
  geom_sf(
    aes(fill = biome_name),
    col = "black",
    alpha = 0.8,
    linewidth = 0.08
  ) +
  scale_fill_manual(values = c("#648fff", "#785ef0", "#dc267f"), name = "") +
  theme_void() +
  theme(
    legend.position = c(0.2, 0.9),
    panel.border = element_blank(),
    legend.text = element_text(size = 8),
    legend.key.height = unit(5, "mm"),
    legend.key.width = unit(5, "mm")
  )

ragg::agg_png(
  "figures/fig_s6_ecoregions.jpg",
  width = 89,
  height = 89,
  units = "mm",
  res = 500
)
p_fig_s6
dev.off()
