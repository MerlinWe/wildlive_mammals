### mapping and land-use change descriptives

rm(list = ls())
library(terra)
library(sf)
library(ggthemes)
library(tidyverse)
library(cowplot)
library(gtable)
library(stringr)
library(patchwork)
library(ggspatial)
library(rnaturalearth)
library(rnaturalearthdata)


# read input data
captures <- read_csv("Data/species_data.csv")
camtraps <- read_csv("Data/camtraps_clean.csv")
covariates <- read_csv("Data/forest_covariates.csv")

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

# NDVI rasters as SpatRaster list 
ndvi_files <- list.files("Remote Sensing/NDVI (Annual means)", pattern = "ndvi_20\\d{2}\\.grd$", full.names = TRUE)
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

ggsave("Output/Figures/supp_4.png",
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
    min = min(delta_treecover, na.rm = TRUE),
    max = max(delta_treecover, na.rm = TRUE),
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

# Figure 1
ggsave("Output/Figures/fig1.png",
       last_plot(),
       width = 180,
       height = 80,
       units = "mm",
       dpi = 600)

## PANEL A

research_area_wgs84 <- st_transform(research_area_outline, 4326)
study_bbox_poly <- st_as_sfc(st_bbox(research_area_wgs84), crs = 4326)
study_centroid  <- st_centroid(research_area_wgs84)

bolivia <- ne_countries(scale = "medium", returnclass = "sf") %>%
  filter(name == "Bolivia")
south_america <- ne_countries(continent = "South America", 
                              scale = "medium", returnclass = "sf")
main_map <- ggplot() +
  geom_sf(data = bolivia, fill = "grey95", color = "grey40", linewidth = 0.4) +
  geom_sf(data = study_bbox_poly, fill = NA, color = "black", linewidth = 1.1) +
  geom_sf(data = study_centroid, color = "black", size = 2) +
  coord_sf(expand = FALSE) +
  theme_void(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0, face = "bold", size = 12),
    plot.margin = margin(6, 6, 6, 6)
  ) +
  annotation_scale(location = "br", text_cex = 0.5, line_width = 0.4, height = unit(0.2, "cm")) +
  annotation_north_arrow(location = "tr", which_north = "true",
                         style = north_arrow_fancy_orienteering(text_size = 6))

panel_A <- ggdraw() +
  draw_plot(main_map)
panel_A
ggsave("Output/Figures/fig1/panel_A.png", panel_A,
       width = 95, height = 75, units = "mm", dpi = 300)

## PANEL B
plot_extent <- st_buffer(research_area_outline, 1)

#  crop/mask & convert to df 
ndvi_to_df <- function(ndvi_rast, clip_poly) {
  ndvi_clip <- terra::mask(terra::crop(ndvi_rast, vect(clip_poly)), vect(clip_poly))
  df <- as.data.frame(ndvi_clip, xy = TRUE, na.rm = TRUE)
  names(df)[3] <- "ndvi"
  df
}

years_for_maps <- c("2017","2023")
ndvi_vals <- do.call(rbind, lapply(years_for_maps, function(y) {
  df <- ndvi_to_df(ndvi_list[[y]], plot_extent); df$year <- y; df
}))
ndvi_min <- quantile(ndvi_vals$ndvi, 0.01, na.rm = TRUE)
ndvi_max <- quantile(ndvi_vals$ndvi, 0.99, na.rm = TRUE)

# year for Panel B ?
backdrop_year <- "2019"  
ndvi_bg <- ndvi_to_df(ndvi_list[[backdrop_year]], plot_extent)

buffer_edge  <- "black"
buffer_fill  <- scales::alpha(buffer_edge, 0.12)

panel_B <- ggplot() +
  # NDVI backdrop as continuous raster, low opacity
  geom_raster(data = ndvi_bg,
              aes(x = x, y = y, fill = ndvi),
              alpha = 0.35) + 
  scale_fill_viridis_c(option = "viridis",
                       limits = c(ndvi_min, ndvi_max),
                       name = "NDVI") +
  
  # Buffers and stations on top
  geom_sf(data = station_buffer_1500, fill = buffer_fill,
          color = buffer_edge, linewidth = 0.5) +
  coord_sf(crs = 32720,
           xlim = st_bbox(plot_extent)[c("xmin","xmax")],
           ylim = st_bbox(plot_extent)[c("ymin","ymax")],
           expand = FALSE) +
  
  geom_sf(data = st_as_sf(plot_extent), fill = NA, color = "black", linewidth = 0.6) +
  coord_sf(crs = 32720,
           xlim = st_bbox(plot_extent)[c("xmin","xmax")],
           ylim = st_bbox(plot_extent)[c("ymin","ymax")],
           expand = FALSE) +
  
  theme_void(base_size = 11) +
  theme(plot.title = element_text(hjust = 0, face = "bold", size = 12),
        plot.subtitle = element_text(hjust = 0, size = 10),
        plot.margin = margin(6, 6, 6, 6),
        legend.position = "none")               

panel_B

ggsave("Output/Figures/fig1/panel_B.png", panel_B,
       width = 95, height = 75, units = "mm", dpi = 300)

## PANLES C & D 

make_ndvi_panel <- function(year, show_legend = TRUE, title = year) {
  df <- ndvi_to_df(ndvi_list[[year]], plot_extent)
  p <- ggplot() +
    geom_raster(data = df, aes(x = x, y = y, fill = ndvi), alpha = 1) +
    
    geom_sf(data = station_buffer_1500, fill = NA,
            color = buffer_edge, linewidth = 0.5) +
    
    scale_fill_viridis_c(option = "viridis",
                         limits = c(ndvi_min, ndvi_max),
                         name = "NDVI") +
    coord_sf(crs = 32720,
             xlim = st_bbox(plot_extent)[c("xmin","xmax")],
             ylim = st_bbox(plot_extent)[c("ymin","ymax")],
             expand = FALSE) +
    labs(title = title) +
    theme_void(base_size = 11) +
    theme(plot.title = element_text(hjust = 0, face = "bold", size = 12),
          plot.margin = margin(6, 6, 6, 6),
          legend.position = if (show_legend) "right" else "none")
  p
}

panel_C <- make_ndvi_panel("2017", show_legend = FALSE)
panel_D <- make_ndvi_panel("2023", show_legend = FALSE)
panel_legend <- make_ndvi_panel("2023", show_legend = TRUE)

ggsave("Output/Figures/fig1/panel_C.png", panel_C, width = 95, height = 75, units = "mm", dpi = 300)
ggsave("Output/Figures/fig1/panel_D.png", panel_D, width = 95, height = 75, units = "mm", dpi = 300)
ggsave("Output/Figures/fig1/panel_legend.png", panel_legend, width = 95, height = 75, units = "mm", dpi = 300)
