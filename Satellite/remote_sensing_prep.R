###########################################################################
##### WildLive - Trait-Filtered Responses, Covariate Data Preparation ##### 
###########################################################################

rm(list=ls()); gc()

library(cowplot)
library(raster) 
library(terra)
library(mapview)
library(sf)        
library(sp)
library(betareg) 
library(car)   
library(landscapemetrics)
library(mgcv)
library(ggeffects)
library(GGally)
library(tidyverse) 

camtraps <- read_csv("/Users/serpent/Documents/Senckenberg/WildLive/Mammals/Code/Data/camtraps_clean.csv")
species <- read_csv("/Users/serpent/Documents/Senckenberg/WildLive/Mammals/Code/Data/species_data.csv")

## Build buffer zones
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

mapview(station_buffer_1500, map.types = "Esri.WorldImagery", legend = FALSE, lwd = 3, color = "white")

# Get research area as a bounding box polygon to use as a clipping feature 
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

## =================== Forest Cover =================== 

ndvi_2017 <- crop(raster("/Users/serpent/Documents/Senckenberg/WildLive/Mammals/Code/Satellite/NDVI (Annual means)/ndvi_2017.grd"), station_buffer_1500)
ndvi_2018 <- crop(raster("/Users/serpent/Documents/Senckenberg/WildLive/Mammals/Code/Satellite/NDVI (Annual means)/ndvi_2018.grd"), station_buffer_1500)
ndvi_2019 <- crop(raster("/Users/serpent/Documents/Senckenberg/WildLive/Mammals/Code/Satellite/NDVI (Annual means)/ndvi_2019.grd"), station_buffer_1500)
ndvi_2020 <- crop(raster("/Users/serpent/Documents/Senckenberg/WildLive/Mammals/Code/Satellite/NDVI (Annual means)/ndvi_2020.grd"), station_buffer_1500)
ndvi_2021 <- crop(raster("/Users/serpent/Documents/Senckenberg/WildLive/Mammals/Code/Satellite/NDVI (Annual means)/ndvi_2021.grd"), station_buffer_1500)
ndvi_2022 <- crop(raster("/Users/serpent/Documents/Senckenberg/WildLive/Mammals/Code/Satellite/NDVI (Annual means)/ndvi_2022.grd"), station_buffer_1500)
ndvi_2023 <- crop(raster("/Users/serpent/Documents/Senckenberg/WildLive/Mammals/Code/Satellite/NDVI (Annual means)/ndvi_2023.grd"), station_buffer_1500)

par(mfrow = c(2, 4))
stations <- as(as(station_buffer_1500, "Spatial"), "SpatialPolygonsDataFrame")
ext <- extent(ndvi_2017)

plot(ndvi_2017, main = "Mean NDVI 2017 (n=3)", ylim = c(ext@ymin, ext@ymax), mar = c(1, 1, 1, 1)); plot(stations, add = TRUE); 
plot(ndvi_2018, main = "Mean NDVI 2018 (n=3)", ylim = c(ext@ymin, ext@ymax), mar = c(1, 1, 1, 1)); plot(stations, add = TRUE); 
plot(ndvi_2019, main = "Mean NDVI 2019 (n=3)", ylim = c(ext@ymin, ext@ymax), mar = c(1, 1, 1, 1)); plot(stations, add = TRUE);
plot(ndvi_2020, main = "Mean NDVI 2020 (n=4)", ylim = c(ext@ymin, ext@ymax), mar = c(1, 1, 1, 1)); plot(stations, add = TRUE);
plot(ndvi_2020, main = "Mean NDVI 2021 (n=4)", ylim = c(ext@ymin, ext@ymax), mar = c(1, 1, 1, 1)); plot(stations, add = TRUE);
plot(ndvi_2021, main = "Mean NDVI 2022 (n=3)", ylim = c(ext@ymin, ext@ymax), mar = c(1, 1, 1, 1)); plot(stations, add = TRUE);
plot(ndvi_2022, main = "Mean NDVI 2023 (n=2)", ylim = c(ext@ymin, ext@ymax), mar = c(1, 1, 1, 1)); plot(stations, add = TRUE)

# Match projection of research area from UTM to lat/long
mask <- station_buffer_1500 %>%
	st_transform(crs = st_crs("+proj=longlat +datum=WGS84 +no_defs")) %>%
	st_as_sf()

# Crop the GFW data to research extent
trees <- raster::raster(raster::crop(
	# Read GFW tree cover
	terra::rast("/Users/serpent/Documents/Senckenberg/WildLive/Mammals/Code/Satellite/10S_070W.tif"),
	# crop to mask
	mask))

dev.off()
plot(trees, main = "Tree cover in 2020 (%)"); plot(mask, add = TRUE, colour = NA)

# Extract tree cover data inside buffer zones from GFW data 
trees_df <- raster::extract(trees, mask) %>%
	# Calculate mean and SD
	map(~ tibble(
		mean_treecover = mean(.x),
		sd_treecover = sd(.x)
	)) %>%
	# Build databse
	bind_rows() %>%
	# Get stations
	mutate(station = mask$Station) %>%
	# Set formats
	mutate(across(c(mean_treecover, sd_treecover), as.numeric))

# Extract tree cover data inside buffer zones from NDVI data 
ndvi_2020_df <- raster::extract(ndvi_2020, stations) %>%
	map(~ tibble(
		mean_ndvi = mean(.x),
		sd_ndvi = sd(.x)
	)) %>%
	bind_rows() %>%
	mutate(station = mask$Station) %>%
	mutate(across(c(mean_ndvi, sd_ndvi), as.numeric))

# Build a data set to train our regression model 
training_data <- trees_df %>% 
	left_join(ndvi_2020_df, by = "station") %>%
	dplyr::select(station, mean_treecover, mean_ndvi, sd_treecover, sd_ndvi) %>%
	mutate(mean_treecover = round(mean_treecover/100, digits = 3))

# Build compound figure
plot_grid(
	
	ggplot(training_data, aes(x=mean_treecover)) + 
		geom_histogram(bins = 30, colour = "black", fill = "forestgreen", alpha = .75) +
		theme_bw()+
		labs(x="Mean forest cover (%)", y="Frequency"),
	
	ggplot(training_data, aes(y=mean_treecover, x = mean_ndvi)) + 
		geom_point(colour = "grey60", alpha = .75) +
		geom_smooth(method = "lm", formula = "y~x", colour = "darkgreen", fill = "lightgreen")+
		theme_bw()+
		labs(x="Mean NDVI", y="Mean forest cover (%)"),
	
	ggplot(training_data, aes(y=mean_treecover, x = sd_ndvi)) + 
		geom_point(colour = "grey60", alpha = .75) +
		geom_smooth(method = "lm", formula = "y~x", colour = "forestgreen", fill = "lightgreen")+
		theme_bw()+
		labs(x="Standard deviation NDVI", y="Mean forest cover (%)"),
	
	ncol = 3, nrow = 1
)

# Start with the most complex (fully factorial) model 
mod_beta1 <- betareg(mean_treecover ~ mean_ndvi*sd_ndvi ,data = training_data)

# Try a more simple, additive model 
mod_beta2 <- betareg(mean_treecover ~ mean_ndvi + sd_ndvi, data = training_data)

# Try a simple model with just the mean 
mod_beta3 <- betareg(mean_treecover ~ mean_ndvi, data = training_data)

# Try a simple model with just the sd
mod_beta4 <- betareg(mean_treecover ~ sd_ndvi, data = training_data)

# Check model performance 
AIC(mod_beta1, mod_beta2, mod_beta3, mod_beta4)

summary(mod_beta1)
par(mfrow = c(1, 4))
plot(mod_beta1)

# To check for multicollinearity, we calculate the correlation coefficient and the VIF
cor(training_data$mean_ndvi, training_data$sd_ndvi) # not great 
vif(mod_beta1) # full factorial model

summary(mod_beta2)
par(mfrow = c(1, 4))
plot(mod_beta2)

vif(mod_beta2) # additive model 

# Predict tree cover of our training data set 
training_data$predicted_mean_treecover <- predict(mod_beta2, newdata = training_data)

# Linear R-squared : Indicates the proportion of variance explained by the model
rsquared <- cor(training_data$predicted_mean_treecover, training_data$mean_treecover)^2 

# Root Mean Squared Error (RMSE): Measures the average difference between predicted values and actual values.
rmse <- sqrt(mean((training_data$predicted_mean_treecover - training_data$mean_treecover)^2)) 

# Mean Absolute Error (MAE): Measures the average absolute differences between predicted and actual values.
mae <- mean(abs(training_data$predicted_mean_treecover - training_data$mean_treecover))

# mod2 is our final model :) let's call it tree_mod
tree_mod <- betareg(mean_treecover ~ mean_ndvi + sd_ndvi, data = training_data)

# Summary of metrics
cat("R-squared:", rsquared, "\n")
cat("RMSE:", rmse, "\n")
cat("MAE:", mae, "\n")

rm(mod_beta1, mod_beta2, mod_beta3, mod_beta4, mae, rmse)

# Next, we get mean and SD NDVI values for the years 2017 - 2023 and predict forest cover using our model. 

# Extract 2017 NDVI data from station buffer polygons
ndvi_2017_df <- raster::extract(ndvi_2017, stations) %>%
	# Extract mean and SD
	map(~ tibble(
		mean_ndvi = mean(.x),
		sd_ndvi = sd(.x)
	)) %>%
	# Build database
	bind_rows() %>%
	# Get station names
	mutate(station = mask$Station) %>%
	# Set formats
	mutate(across(c(mean_ndvi, sd_ndvi), as.numeric)) %>%
	# Predict true tree cover
	mutate(predicted_mean_treecover = predict(tree_mod, newdata = .)) %>%
	# Add year variable 
	mutate(year = "2017")

# 2018
ndvi_2018_df <- raster::extract(ndvi_2018, stations) %>%
	map(~ tibble(
		mean_ndvi = mean(.x),
		sd_ndvi = sd(.x)
	)) %>%
	bind_rows() %>%
	mutate(station = mask$Station) %>%
	mutate(across(c(mean_ndvi, sd_ndvi), as.numeric)) %>%
	mutate(predicted_mean_treecover = predict(tree_mod, newdata = .)) %>%
	mutate(year = "2018")

# 2019
ndvi_2019_df <- raster::extract(ndvi_2019, stations) %>%
	map(~ tibble(
		mean_ndvi = mean(.x),
		sd_ndvi = sd(.x)
	)) %>%
	bind_rows() %>%
	mutate(station = mask$Station) %>%
	mutate(across(c(mean_ndvi, sd_ndvi), as.numeric)) %>%
	mutate(predicted_mean_treecover = predict(tree_mod, newdata = .)) %>%
	mutate(year = "2019")

# 2020 (this is already in our training data)
ndvi_2020_df <- training_data %>% 
	dplyr::select(mean_ndvi, sd_ndvi, station, predicted_mean_treecover) %>%
	mutate(year = "2020")

# 2021
ndvi_2021_df <- raster::extract(ndvi_2021, stations) %>%
	map(~ tibble(
		mean_ndvi = mean(.x),
		sd_ndvi = sd(.x)
	)) %>%
	bind_rows() %>%
	mutate(station = mask$Station) %>%
	mutate(across(c(mean_ndvi, sd_ndvi), as.numeric)) %>%
	mutate(predicted_mean_treecover = predict(tree_mod, newdata = .)) %>%
	mutate(year = "2021")

# 2022
ndvi_2022_df <- raster::extract(ndvi_2022, stations) %>%
	map(~ tibble(
		mean_ndvi = mean(.x),
		sd_ndvi = sd(.x)
	)) %>%
	bind_rows() %>%
	mutate(station = mask$Station) %>%
	mutate(across(c(mean_ndvi, sd_ndvi), as.numeric)) %>%
	mutate(predicted_mean_treecover = predict(tree_mod, newdata = .)) %>%
	mutate(year = "2022")

# 2023
ndvi_2023_df <- raster::extract(ndvi_2023, stations) %>%
	map(~ tibble(
		mean_ndvi = mean(.x),
		sd_ndvi = sd(.x)
	)) %>%
	bind_rows() %>%
	mutate(station = mask$Station) %>%
	mutate(across(c(mean_ndvi, sd_ndvi), as.numeric)) %>%
	mutate(predicted_mean_treecover = predict(tree_mod, newdata = .)) %>%
	mutate(year = "2023")

# Build database
treecover <- bind_rows(ndvi_2017_df, ndvi_2018_df, ndvi_2019_df, ndvi_2020_df, ndvi_2021_df, ndvi_2022_df, ndvi_2023_df) %>%
	as_tibble() %>%
	dplyr::select(station, year, predicted_mean_treecover) %>%
	# Format treecover as actual percentage (1-100)
	mutate(treecover = predicted_mean_treecover *100)

rm(ndvi_2017_df, ndvi_2018_df, ndvi_2019_df, ndvi_2020_df, ndvi_2021_df, ndvi_2022_df, ndvi_2023_df, research_area_outline, station_buffer_1500, training_data, tree_mod, rsquared)


## =================== Fragmentation =================== 


# Define forest thresholds
forest_reclass <- data.frame(from = c(-Inf, .65), to = c(.65, Inf), becomes = c(0, 1))
# Write function to reclassify listed raster layers 
reclassify_raster <- function(raster_layer, reclass_matrix) {
	reclassified <- reclassify(raster_layer, rcl = as.matrix(reclass_matrix), right = FALSE)
	return(reclassified)
}

# Write function to extract fragmentation metrics
extract_fragmentation_metrics <- function(ndvi_layer, stations, year, reclass_matrix) {
	masked <- map(1:nrow(stations), ~ mask(ndvi_layer, stations[.x,])) %>%
		set_names(stations$Station)
	
	reclassified <- map(masked, reclassify_raster, reclass_matrix = reclass_matrix)
	
	metrics_to_calc <- list(
		ed = lsm_c_ed,
		pd = lsm_c_pd,
		shape = lsm_c_shape_mn
	)
	
	frag_metrics <- map(metrics_to_calc, function(metric_fun) {
		map(reclassified, metric_fun) %>%
			bind_rows(.id = "Station") %>%
			filter(class == 1) %>%
			dplyr::select(Station, value)
	})
	
	out <- reduce(names(frag_metrics), function(df, metric) {
		df %>%
			left_join(
				frag_metrics[[metric]] %>% rename(!!metric := value),
				by = "Station"
			)
	}, .init = tibble(Station = stations$Station)) %>%
		mutate(year = year) %>%
		rename(
			edge_density_forest = ed,
			patch_density_forest = pd,
			shape_index_forest = shape
		)
	
	return(out)
}

frag_2017 <- extract_fragmentation_metrics(ndvi_2017, stations, 2017, forest_reclass)
frag_2018 <- extract_fragmentation_metrics(ndvi_2018, stations, 2018, forest_reclass)
frag_2019 <- extract_fragmentation_metrics(ndvi_2019, stations, 2019, forest_reclass)
frag_2020 <- extract_fragmentation_metrics(ndvi_2020, stations, 2020, forest_reclass)
frag_2021 <- extract_fragmentation_metrics(ndvi_2021, stations, 2021, forest_reclass)
frag_2022 <- extract_fragmentation_metrics(ndvi_2022, stations, 2022, forest_reclass)
frag_2023 <- extract_fragmentation_metrics(ndvi_2023, stations, 2023, forest_reclass)

fragmentation_metrics <- bind_rows(frag_2017, frag_2018, frag_2019, frag_2020, frag_2021, frag_2022, frag_2023) %>%
	rename(station = Station) %>%
	mutate(year = as.character(year))

rm(buffer, forest_reclass, frag_2017, frag_2018, frag_2019, frag_2020, frag_2021, frag_2022, frag_2023, mask, ndvi_2017, ndvi_2018, ndvi_2019, ndvi_2020, ndvi_2021, ndvi_2022, ndvi_2023, stations, i, reclassify_raster)

# Find station-year combinations where cattle were found
cattle <- species %>%
	# Keep only cattle
	filter(Trivial == "cattle") %>%
	# Keep only grid stations
	filter(Station %in% unique(treecover$station)) %>%
	# Get year variable
	mutate(year = year(DateTimeOriginal)) %>%
	filter(year %in% unique(treecover$year)) %>%
	mutate(year = as.character(year)) %>%
	# Group by Station and year
	group_by(Trivial, Station, year) %>%
	# Summarise cattle presence/absence
	summarise(agriculture = ifelse(n() > 1, 1, 0)) %>%
	ungroup() %>%
	# Fix capitalization
	rename(station = Station)

# Build databae
covariates <- treecover %>% 
	left_join(fragmentation_metrics, by = c("station", "year")) %>%
	dplyr::select(-predicted_mean_treecover)

write_csv(covariates, file = "/Users/serpent/Documents/Senckenberg/WildLive/Mammals/Code/Data/forest_covariates.csv")
rm(cattle, fragmentation, treecover, camtraps)

## Examine covariates

covariates %>%
	dplyr::select(treecover, edge_density_forest, patch_density_forest, shape_index_forest) %>%
	ggpairs()


# Standardize each metric by year and Plot standardized values
covariates %>%
	filter(year != "2023") %>%
	mutate(year = factor(year, levels = 2017:2022)) %>%
	pivot_longer(cols = -c(station, year),
							 names_to = "Metric", 
							 values_to = "value") %>%
	group_by(Metric) %>%
	mutate(value_z = scale(value)[,1]) %>%
	ungroup() %>%
	
	ggplot(aes(x = year, y = value_z, group = Metric, color = Metric)) +
	geom_line() +
	geom_point(size = 1) +
	scale_color_manual(values = c("forestgreen", "black", "red", "blue")) +
	facet_wrap(~ station) +
	theme_bw() +
	theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1)) +
	labs(x = NULL, y = "Standardized Value (Z-score)")






