### mapping and land-use change descriptives

rm(list = ls())
library(terra)
library(sf)
library(tidyverse)
library(cowplot)
library(gtable)
library(ggspatial)

# read input data
captures <- read_csv("/Users/merlin/Documents/Senckenberg/WildLive/Mammals/Code/Data/species_data.csv")
camtraps <- read_csv("/Users/merlin/Documents/Senckenberg/WildLive/Mammals/Code/Data/camtraps_clean.csv")
covariates <- read_csv("/Users/merlin/Documents/Senckenberg/WildLive/Mammals/Code/Data/forest_covariates.csv")

# Station buffers (1500m radius)
station_buffer_1500 <- camtraps %>%
  distinct(Station, lat, long, .keep_all = TRUE) %>%
  st_as_sf(coords = c("long", "lat"), crs = 4326) %>%
  st_transform(crs = 32720) %>%
  st_buffer(dist = 1500)

# Create convex hull bounding box and extend by 10km
research_area_outline <- station_buffer_1500 %>%
  st_bbox() %>%
  st_as_sfc() %>%
  st_buffer(10000)

# NDVI rasters as SpatRaster list (use terra only)
ndvi_files <- list.files("/Users/merlin/Documents/Senckenberg/WildLive/Mammals/Code/Remote Sensing/NDVI (Annual means)", pattern = "ndvi_20\\d{2}\\.grd$", full.names = TRUE)
ndvi_list <- lapply(ndvi_files, rast)
names(ndvi_list) <- str_extract(ndvi_files, "20\\d{2}")

# Covariates in long format
cov_long <- covariates %>%
  dplyr::select(-patch_density_forest, -shape_index_forest) %>%
  pivot_longer(c(treecover, edge_density_forest),
               names_to = "Metric", values_to = "Value") %>%
  mutate(
    year = factor(year),
    Metric = dplyr::recode(Metric,
                    treecover = "Tree cover",
                    edge_density_forest = "Edge Density")
  )

cov_long %>%
  group_by(Metric) %>%
  mutate(z_value = scale(Value)[, 1]) %>%
  ungroup() %>%
  ggplot(aes(x = year, y = z_value, color = Metric, group = interaction(station, Metric))) +
  geom_line(size = 0.8) +
  facet_wrap(~ station, scales = "free_y", ncol = 3, nrow = 7) +
  scale_color_manual(values = c("forestgreen", "goldenrod"),
                     labels = c("Forest cover", "Edge density"),
                     name = "Metric") +
  labs(x = NULL, y = "Z-transformed value") +
  theme_few(base_size = 12) +
  theme(
    strip.text = element_text(size = 12),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(hjust = 1, size = 10, angle = 45),  
    legend.position = "top",
    plot.margin = margin(t = 5, r = 20, b = 5, l = 10),
    panel.background = element_blank(),
    plot.background = element_blank(),
    strip.background = element_blank()
  )

ggsave("/Users/merlin/Documents/Senckenberg/WildLive/Mammals/Code/Output/Figures/supp_4.png",
       last_plot(),
       width = 200,
       height = 250,
       units = "mm",
       dpi = 600,
       bg = "white")

# Covariate changes 2017–2023
covariate_changes <- covariates %>%
  filter(year %in% c(2017, 2023)) %>%
  dplyr::select(station, year, treecover, edge_density_forest) %>%
  pivot_wider(names_from = year, values_from = c(treecover, edge_density_forest)) %>%
  mutate(
    delta_treecover = treecover_2023 - treecover_2017,
    delta_fragmentation = edge_density_forest_2023 - edge_density_forest_2017
  )

# Summary stats
covariate_changes %>%
  summarise(
    mean = mean(delta_treecover, na.rm = TRUE),
    sd = sd(delta_treecover, na.rm = TRUE),
    loss = sum(delta_treecover < 0, na.rm = TRUE),
    gain = sum(delta_treecover > 0, na.rm = TRUE),
    unchanged = sum(delta_treecover == 0, na.rm = TRUE))

covariate_changes %>%
  summarise(
    mean = mean(delta_fragmentation, na.rm = TRUE),
    sd = sd(delta_fragmentation, na.rm = TRUE),
    loss = sum(delta_fragmentation < 0, na.rm = TRUE),
    gain = sum(delta_fragmentation > 0, na.rm = TRUE),
    unchanged = sum(delta_fragmentation == 0, na.rm = TRUE))

# Station points
stations <- st_transform(station_buffer_1500, crs = 32720)




prep_ndvi_for_plot <- function(ndvi_raster) {
  reclass_matrix <- matrix(c(-Inf, 0.64, 0, 0.64, Inf, 1), ncol = 3, byrow = TRUE)
  rc <- classify(ndvi_raster, rcl = reclass_matrix)
  
  # Extract x, y and values
  pts <- as.data.frame(rc, xy = TRUE, na.rm = TRUE)
  
  # Assign labels
  pts$landcover <- factor(pts$mean,
                          levels = c(0, 1),
                          labels = c("No forest", "Forest"))
  return(pts)
}

basemap_2017_df <- prep_ndvi_for_plot(ndvi_list[["2017"]])
basemap_2023_df <- prep_ndvi_for_plot(ndvi_list[["2023"]])


p1 <- basemap_2017_df %>%
  ggplot() +
  geom_raster(aes(x = x, y = y, fill = landcover)) +
  scale_fill_manual(values = c("oldlace", "darkgreen"), name = NULL) +
  geom_sf(data = station_buffer_1500, color = "black", fill = NA, size = 0.5) +
  coord_sf(crs = 32720, expand = FALSE) +
  theme_void(base_size = 11) +
  labs(title = "2017") +
  theme(legend.position = "bottom")

p2 <- basemap_2023_df %>%
  ggplot() +
  geom_raster(aes(x = x, y = y, fill = landcover)) +
  scale_fill_manual(values = c("oldlace", "darkgreen"), name = NULL) +
  geom_sf(data = station_buffer_1500, color = "black", fill = NA, size = 0.5) +
  coord_sf(crs = 32720, expand = FALSE) +
  theme_void(base_size = 11) +
  labs(title = "2023") +
  theme(legend.position = "bottom")


final_map <- (p1 + plot_spacer() + p2) +
  plot_layout(ncol = 3, widths = c(1, 0.05, 1), guides = "collect") &
  theme(legend.position = "bottom")

# Figure 2
ggsave("/Users/merlin/Documents/Senckenberg/WildLive/Mammals/Code/Output/Figures/fig1.png",
       last_plot(),
       width = 180,
       height = 80,
       units = "mm",
       dpi = 600)

