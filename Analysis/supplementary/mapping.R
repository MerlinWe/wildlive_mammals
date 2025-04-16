### mapping and land-use change descriptives

pathname <- "/Users/serpent/Library/CloudStorage/GoogleDrive-merlin.s.weiss@gmail.com/My Drive/WildLive! /Report"

library(cowplot)
library(raster) 
library(terra)
library(sf)        
library(sp)
library(grid)
library(gtable)
library(gridExtra)
library(gridGraphics)
library(ggspatial)
library(tidyverse) 

camtraps <- read_csv(paste0(pathname, "/camtraps.csv")) 

station_buffer_1500 <- camtraps %>% 
	# Keep only relevant columns
	dplyr::select(Station, lat, long) %>% 
	# Keep only one row per Station 
	distinct(Station, .keep_all = TRUE) %>%
	# Transform to spatial object 
	st_as_sf(coords = c("long", "lat"), crs = 4326) %>%
	# Generate buffer zones with 1500m radius 
	st_buffer(dist = 1500) %>%
	# Match format and projections
	st_as_sf() %>%
	st_transform(crs = st_crs("+proj=utm +zone=20 +south +datum=WGS84 +units=m +no_defs"))


research_area_outline <- station_buffer_1500 %>%
	# define hull convex as bounding box 
	st_bbox() %>%
	st_make_grid(cellsize = c(
		st_bbox(station_buffer_1500)[3]-st_bbox(station_buffer_1500)[1],
		st_bbox(station_buffer_1500)[4]-st_bbox(station_buffer_1500)[2])) %>%
	# Extend my 10km 
	st_buffer(dist = 10000) %>%
	st_as_sf() %>%
	st_transform(crs = st_crs("+proj=utm +zone=20 +south +datum=WGS84 +units=m +no_defs")) %>%
	# Simply extend 
	st_as_sfc(st_bbox())

ndvi_2017 <- raster("/Users/serpent/Documents/BSc/Thesis/Data/Sat/NDVI (Annual means)/ndvi_2017.grd")
ndvi_2018 <- raster("/Users/serpent/Documents/BSc/Thesis/Data/Sat/NDVI (Annual means)/ndvi_2018.grd")
ndvi_2019 <- crop(raster("/Users/serpent/Documents/BSc/Thesis/Data/Sat/NDVI (Annual means)/ndvi_2019.grd"), station_buffer_1500)
ndvi_2020 <- raster("/Users/serpent/Documents/BSc/Thesis/Data/Sat/NDVI (Annual means)/ndvi_2020.grd")
ndvi_2021 <- raster("/Users/serpent/Documents/BSc/Thesis/Data/Sat/NDVI (Annual means)/ndvi_2021.grd")
ndvi_2022 <- raster("/Users/serpent/Documents/BSc/Thesis/Data/Sat/NDVI (Annual means)/ndvi_2022.grd")
ndvi_2023 <- raster("/Users/serpent/Documents/BSc/Thesis/Data/Sat/NDVI (Annual means)/ndvi_2023.grd")

par(mfrow = c(2, 4))
stations <- as(as(station_buffer_1500, "Spatial"), "SpatialPolygonsDataFrame")
ext <- extent(ndvi_2017)


covariates <- read_csv("/Users/serpent/Desktop/Paper/covariates.csv")

# Prepare data in long format
cov_long <- covariates %>%
	filter(year != 2023, station != "G-21", station %in% camtraps$Station) %>%
	pivot_longer(cols = c(treecover, aggregation_forest), names_to = "Metric", values_to = "Value") %>%
	mutate(
		year = factor(year),
		Metric = recode(Metric,
										"treecover" = "Forest cover",
										"aggregation_forest" = "Fragmentation")
	)

# Build plot
landuse_plot <- ggplot(cov_long, aes(x = year, y = Value, group = Metric, color = Metric)) +
	# Add agriculture background shading
	geom_tile(
		data = covariates %>%
			filter(year != 2023, station != "G-21", station %in% camtraps$Station) %>%
			mutate(year = factor(year)) %>%
			mutate(Agriculture = factor(agriculture)),
		aes(x = year, y = 0, fill = Agriculture),
		width = 1, height = Inf, alpha = 0.2, inherit.aes = FALSE
	) +
	geom_line() +
	geom_point(size = 1) +
	scale_color_manual(values = c("Forest cover" = "forestgreen", "Fragmentation" = "black")) +
	scale_fill_manual(values = c("1" = "orange", "0" = "white"), labels = c("No", "Yes")) +
	facet_wrap(~ station) +
	theme_bw() +
	theme(
		axis.text.x = element_text(angle = 45, vjust = 0.5),
		legend.position = "top",
		legend.box = "horizontal"
	) +
	labs(
		x = NULL, y = "Percentage (%)",
		fill = "Agricultural activity", color = "Landscape metric"
	)

ggsave("/Users/serpent/Desktop/Paper/landuse_plot.png", 
			 landuse_plot, 
			 width = 220,
			 height = 200, 
			 units = "mm",
			 dpi = 600)

covariate_changes <- covariates %>%
	filter(year %in% c(2017, 2023)) %>%
	select(station, year, treecover, aggregation_forest) %>%
	pivot_wider(names_from = year, values_from = c(treecover, aggregation_forest)) %>%
	mutate(
		delta_treecover = treecover_2023 - treecover_2017,
		delta_fragmentation = aggregation_forest_2023 - aggregation_forest_2017
	)

print(covariate_changes, n = 21)

summary_stats <- covariate_changes %>%
	summarise(
		mean_treecover_change = mean(delta_treecover, na.rm = TRUE),
		sd_treecover_change = sd(delta_treecover, na.rm = TRUE),
		stations_loss = sum(delta_treecover < 0, na.rm = TRUE),
		stations_gain = sum(delta_treecover > 0, na.rm = TRUE),
		stations_nochange = sum(delta_treecover == 0, na.rm = TRUE)
	)



## mapping
basemap <- raster("/Users/serpent/Documents/BSc/Thesis/Data/Sat/NDVI (Annual means)/ndvi_2017.grd")

## Set correct spatial projections 
stations <- camtraps %>% dplyr::select(Station, lat, long) %>% distinct(Station, .keep_all = TRUE)
stations <- st_as_sf(stations, coords = c("long", "lat"), crs = 4326) # sf object
stations <- st_transform(stations, crs = st_crs("+proj=utm +zone=20 +south +datum=WGS84 +units=m +no_defs")) 

## Build buffer zone. Use 1.5km radius (2 would be better, see Semper-Pascual et al. (2023), 
## but that would create too much overlap between buffer zones). 

station_buffer_1500 <- st_buffer(stations, dist = 1500) # buffer zone 
station_buffer_1500 <- st_as_sf(station_buffer_1500) # set format & match projection
#station_buffer_1500 <- st_transform(station_buffer_1500, crs = st_crs("+proj=utm +zone=20 +south +datum=WGS84 +units=m +no_defs")) 

# plot(basemap, main = "basemap")
# plot(station_buffer_1500, add = TRUE) 
# plot(stations, colour = "red", add = TRUE)

# OK 

basemap_2017 <- raster("/Users/serpent/Documents/BSc/Thesis/Data/Sat/NDVI (Annual means)/ndvi_2017.grd")

basemap_2017 <- reclassify(basemap_2017, 
													 rcl = as.matrix(data.frame(
													 	from = c(-Inf,.15,.4,.6), 
													 	to = c(.15,.4, .6, Inf), 
													 	becomes = c(1,2,3,4))), 
													 right = FALSE)

basemap_2017 <- as_tibble(rasterToPoints(basemap_2017))
colnames(basemap_2017) <- c("x", "y", "landuse")

basemap_2017$landuse <- factor(basemap_2017$landuse,
															 levels = c(1, 2, 3, 4),
															 labels = c("Minimal vegetation", 
															 										"Low vegetation (e.g. agriculture)", 
															 										"Intermediate vegetation", 
															 										"High vegetation (dense canopy)"))

map_start <- ggplot() +
	# Basemap
	geom_raster(data = basemap_2017, aes(x = x, y = y, fill = landuse)) +
	
	scale_fill_manual(
		values = c("grey95", "wheat1", "darkseagreen", "darkgreen"),
		name = "Landuse") +
	
	# Buffer zones
	#geom_sf(data = station_buffer_1500, color = "grey60", linetype = "dashed", size = 0.3) +
	
	# Station points
	geom_sf(data = stations, shape = 21, fill = "white", color = "black", size = 2.5, stroke = 1) +
	
	# CRS and annotation
	coord_sf(crs = 32720, expand = FALSE) +
	# Theme
	theme_bw() +
	theme(text = element_text(size = 11),
				legend.position = "none") +
	labs(x = NULL, y = NULL, title = "2017")




basemap <- raster("/Users/serpent/Documents/BSc/Thesis/Data/Sat/NDVI (Annual means)/ndvi_2022.grd")

basemap <- reclassify(basemap, # Reclassify raster to show roads (1), agriculture (2), and forest (3/4)
											rcl = as.matrix(data.frame(
												from = c(-Inf,.15,.4,.6), to = c(.15,.4, .6, Inf), 
												becomes = c(1,2,3,4))), right = FALSE) 
basemap <- as_tibble(rasterToPoints(basemap)) # tranfrom to data frame
colnames(basemap) <- c("x", "y", "landuse")

basemap$landuse <- factor(basemap$landuse,
													levels = c(1, 2, 3, 4),
													labels = c("Minimal vegetation", 
																		 "Low vegetation (e.g. agriculture)", 
																		 "Intermediate vegetation", 
																		 "High vegetation (dense canopy)"))

map_end <- ggplot() +
	# Basemap
	geom_raster(data = basemap, aes(x = x, y = y, fill = landuse)) +
	
	scale_fill_manual(
		values = c("grey95", "wheat1", "darkseagreen", "darkgreen"),
		name = NULL) +
	
	# Buffer zones
	#geom_sf(data = station_buffer_1500, fill = NA, color = "grey60", linetype = "dashed", size = 0.3) +
	
	# Station points
	geom_sf(data = stations, shape = 21, fill = "white", color = "black", size = 2.5, stroke = 1) +
	
	# CRS and annotation
	coord_sf(crs = 32720, expand = FALSE) +
	ggspatial::annotation_scale(
		location = "br",
		bar_cols = c("black", "white"),
		pad_y = unit(.45, "cm"),
		text_col = "white") +
	
	# Theme
	theme_bw() +
	theme(text = element_text(size = 11),
				legend.position = "none") +
	labs(x = NULL, y = NULL, title = "2022")

combined_maps <- plot_grid(
	map_start,
	map_end,
	ncol = 1,  # stacked vertically
	align = "v")

legend_plot <- ggplot() +
	geom_raster(data = basemap_2017, aes(x = x, y = y, fill = landuse)) +
	scale_fill_manual(
		values = c("grey95", "wheat1", "darkseagreen", "darkgreen"),
		name = NULL
	) +
	theme_minimal(base_size = 10) +
	theme(
		legend.position = "top"
	)

g <- ggplotGrob(legend_plot)
legend <- gtable_filter(g, "guide-box")


final_map_plot <- plot_grid(
	legend,
	combined_maps,
	ncol = 1,
	rel_heights = c(0.1, 1))

ggsave("/Users/serpent/Desktop/Paper/map.png", 
			 final_map_plot, 
			 bg = "white",
			 width = 200,
			 height = 220, 
			 units = "mm",
			 dpi = 600)
