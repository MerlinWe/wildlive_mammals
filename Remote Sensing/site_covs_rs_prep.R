###########################################################################
##### WildLive - Trait-Filtered Responses, Covariate Data Preparation ##### 
###########################################################################

rm(list = ls()); gc()

# Load packages
library(tidyverse)
library(terra)
library(landscapemetrics)
library(GGally)
library(sf)
library(fs)
library(car)
library(pROC)
library(betareg)
library(cowplot)

# Set working directory
ndvi_base <- file.path("Remote Sensing", "NDVI")
output_dir <- file.path("Remote Sensing", "NDVI (Annual Means)")
# dir_tree()

camtraps <- read_csv("Data/camtraps_clean.csv")
species <- read_csv("Data/species_data.csv")

# Construct spatial mask from station buffers
station_buffer_1500 <- camtraps %>%
	select(Station, lat, long) %>%
	distinct(Station, .keep_all = TRUE) %>%
	st_as_sf(coords = c("long", "lat"), crs = 4326) %>%
	st_transform("EPSG:32720") %>%
	st_buffer(dist = 1500)

# Outline
crop_vect <- st_buffer(st_as_sfc(st_bbox(station_buffer_1500)), dist = 10000) %>%
	st_as_sf() %>%
	vect()

# Stations
mask <- station_buffer_1500 %>%
	st_transform(crs = st_crs("+proj=longlat +datum=WGS84 +no_defs")) %>%
	st_as_sf()

years <- 2017:2023
ndvi_base <- "NDVI"
output_dir <- "NDVI (Annual Means)"

ndvi_metadata <- list()
ndvi_means <- map(years, function(year) {
  
  subdirs <- dir_ls(file.path(project_root, ndvi_base, year), type = "directory")
  
  ndvi_rasters <- list()
  metadata_rows <- list()
  
  for (subdir in subdirs) {
    
    b04 <- dir_ls(subdir, regexp = "B04.*\\.tiff$")
    b08 <- dir_ls(subdir, regexp = "B08.*\\.tiff$")
    
    if (length(b04) == 1 && length(b08) == 1) {
      red <- rast(b04)
      nir <- rast(b08)
      
      if (!as.character(crs(red)) == as.character(crs(crop_vect))) {
        red <- project(red, crs(crop_vect))
        nir <- project(nir, crs(crop_vect))
      }
      
      ndvi <- (nir - red) / (nir + red)
      ndvi <- mask(crop(ndvi, crop_vect), crop_vect)
      
      ndvi_rasters <- append(ndvi_rasters, list(ndvi))
      
      # Extract sensing date (assumed in folder name or filename)
      sensing_date <- NA
      file_name_04 <- basename(b04)
      file_name_08 <- basename(b08)
      date_match <- str_extract(file_name_04, "\\d{4}-\\d{2}-\\d{2}")
      sensing_date <- if (!is.na(date_match)) lubridate::ymd(date_match) else NA
      
      if (!is.na(date_match)) {
        sensing_date <- lubridate::ymd(date_match)
      } else {
        sensing_date <- NA
      }
      
      metadata_rows[[length(metadata_rows) + 1]] <- tibble(
        year = year,
        subdir = subdir,
        b04_file = file_name_04,
        b08_file = file_name_08,
        sensing_date = sensing_date
      )
    } else {
      warning(paste("Missing or multiple B04/B08 files in:", subdir))
    }
  }
  
  if (length(ndvi_rasters) == 0) {
    warning(paste("No valid NDVI rasters found for year", year))
    return(NULL)
  }
  
  ndvi_stack <- rast(ndvi_rasters)
  mean_ndvi <- mean(ndvi_stack, na.rm = TRUE)
  
  out_file <- file.path(project_root, output_dir, paste0("ndvi_", year, ".grd"))
  writeRaster(mean_ndvi, out_file, overwrite = TRUE)
  
  ndvi_metadata[[as.character(year)]] <<- bind_rows(metadata_rows)
  
  return(mean_ndvi)
})

# Combine metadata for all years
ndvi_metadata_table <- bind_rows(ndvi_metadata)
print(ndvi_metadata_table)

names(ndvi_means) <- paste0("ndvi_", years)
list2env(ndvi_means, envir = .GlobalEnv)

# List NDVI rasters and titles
ndvi_years <- 2017:2023
ndvi_rasters <- mget(paste0("ndvi_", ndvi_years))
titles <- paste("Mean NDVI", ndvi_years, "(n=3)")

par(mfrow = c(2, 4), mar = c(2, 2, 2, 1))
ext <- ext(ndvi_rasters[[1]])
stations <- vect(station_buffer_1500)

for (i in seq_along(ndvi_rasters)) {
	plot(ndvi_rasters[[i]],
			 main = titles[i],
			 ylim = c(ext[3], ext[4]))
	plot(stations, add = TRUE)
}

par(mfrow = c(1, 1))

## ---------- Predict true tree cover ----------

# Crop the GFW data to research extent
trees <- treecover <- rast('Remote Sensing/10S_070W.tif')
crop_vect_latlon <- project(crop_vect, crs(treecover))

trees <- crop(treecover, crop_vect_latlon)
trees <- project(trees, crs(mask))

plot(trees, main = "Tree cover in 2020 (%)"); plot(mask, add = TRUE, colour = NA)

# Define the color gradient
tree_col <- colorRampPalette(c("white", "forestgreen"))
plot(trees, col = tree_col(100), main = "Tree cover in 2020 (%)"); plot(mask, add = TRUE, colour = NA)

# Convert mask to terra vector (required for terra::extract)
stations <- vect(mask)

# Extract tree cover data inside buffer zones using terra
trees_df <- terra::extract(trees, stations) %>%
	as_tibble() %>%
	group_by(ID) %>%
	rename(value = `10S_070W`) %>%
	dplyr::summarise(
		mean_treecover = mean(value, na.rm = TRUE),
		sd_treecover   = sd(value, na.rm = TRUE),
		.groups = "drop"
	) %>%
	mutate(station = stations$Station[ID]) %>%
	select(station, mean_treecover, sd_treecover)

# Extract NDVI data from station buffer polygons for 2020
ndvi_2020_df <- terra::extract(ndvi_2020, stations) %>%
	as_tibble() %>%
	group_by(ID) %>%
	summarise(
		mean_ndvi = mean(mean, na.rm = TRUE),
		sd_ndvi = sd(mean, na.rm = TRUE)
	) %>%
	mutate(station = stations$Station) %>%
	select(station, mean_ndvi, sd_ndvi)

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

# Define the years you want to extract (excluding 2020 for now)
years <- c(2017:2019, 2021:2023)

# Extract NDVI rasters from environment
ndvi_rasters <- mget(paste0("ndvi_", years), inherits = TRUE)

# Loop through years and rasters to extract and predict
treecover_list <- map2(ndvi_rasters, years, function(ndvi, yr) {
	df <- terra::extract(ndvi, stations) %>%
		as_tibble() %>%
		group_by(ID) %>%
		summarise(
			mean_ndvi = mean(mean, na.rm = TRUE),
			sd_ndvi   = sd(mean, na.rm = TRUE),
			.groups = "drop"
		) %>%
		mutate(
			station = stations$Station,
			predicted_mean_treecover = predict(tree_mod, newdata = .),
			year = as.character(yr)
		) %>%
		select(station, year, predicted_mean_treecover)
	
	return(df)
})

# Add 2020 manually
ndvi_2020_df <- training_data %>%
	select(station, predicted_mean_treecover) %>%
	mutate(year = "2020")

treecover <- treecover_list %>%
	unname() %>%
	bind_rows() %>%
	bind_rows(ndvi_2020_df) %>%
	mutate(treecover = predicted_mean_treecover * 100)

## =================== Fragmentation =================== 

# Stack and align NDVI and tree cover raster
r_stack <- c(ndvi_2020, 'Remote Sensing/10S_070W.tif' |> 
		rast() |> 
		project(ndvi_2020) |>  
		terra::crop(ndvi_2020))
names(r_stack) <- c("ndvi", "treecover")

# Sample ~10,000 pixels and create binary forest classification
sample_df <- spatSample(r_stack, size = 10000, na.rm = TRUE, xy = FALSE) %>%
	as_tibble() %>%
	mutate(forest_true = ifelse(treecover > 85, 1, 0))

# ROC curve analysis
par(mfrow = c(1, 1))
roc_obj <- roc(sample_df$forest_true, sample_df$ndvi)
plot(roc_obj, col = "darkgreen", lwd = 2, main = "NDVI vs. GFW Tree Cover")
cat("AUC:", auc(roc_obj), "\n")

# Get best threshold from ROC
opt_thresh <- coords(roc_obj, "best", ret = "threshold") %>% pluck("threshold")
cat("Optimal NDVI threshold:", round(opt_thresh, 3), "\n")

# Add binary NDVI classes using default and optimized thresholds
sample_df <- sample_df %>%
	mutate(
		ndvi_class_opt = as.integer(ndvi > opt_thresh))

# confusion matrices
cat("Confusion (opt):\n"); print(table(sample_df$ndvi_class_opt, sample_df$forest_true))

# Build reclassification matrix from optimal threshold
reclass_matrix <- matrix(c(-Inf, opt_thresh, 0,
													 opt_thresh, Inf, 1), ncol = 3, byrow = TRUE)

# Reclassification helper
stations <- terra::project(stations, crs(ndvi_2020))
reclassify_raster <- function(r, rcl) terra::classify(r, rcl = rcl, right = FALSE)

# Fragmentation extractor
extract_fragmentation_metrics <- function(ndvi, stations, year, rcl) {
	map_dfr(1:nrow(stations), function(i) {
		station_id <- stations$Station[i]
		masked <- terra::mask(ndvi, stations[i, ])
		reclassed <- reclassify_raster(masked, rcl)
		
		# Convert to class table (can be inspected)
		tab <- terra::freq(reclassed)
		if (!any(tab$value == 1, na.rm = TRUE)) {
			message(glue::glue("No forest in station {station_id}, year {year}"))
			return(tibble(
				station = station_id,
				year = year,
				edge_density_forest = NA,
				patch_density_forest = NA,
				shape_index_forest = NA
			))
		}
		
		tibble(
			station = station_id,
			year = year,
			edge_density_forest  = tryCatch(lsm_c_ed(reclassed)        %>% filter(class == 1) %>% pull(value), error = function(e) NA),
			patch_density_forest = tryCatch(lsm_c_pd(reclassed)        %>% filter(class == 1) %>% pull(value), error = function(e) NA),
			shape_index_forest   = tryCatch(lsm_c_shape_mn(reclassed)  %>% filter(class == 1) %>% pull(value), error = function(e) NA)
		)
	})
}

# Loop over all years and extract fragmentation
fragmentation_metrics <- map2_dfr(
	mget(paste0("ndvi_", 2017:2023), inherits = TRUE),
	2017:2023,
	~ extract_fragmentation_metrics(.x, stations, .y, reclass_matrix)
)

covariates <- treecover %>%
	mutate(year = as.integer(year)) %>%
	left_join(fragmentation_metrics, by = c("station", "year")) %>%
	select(-predicted_mean_treecover)

glimpse(covariates)
write_csv(covariates, file = "Data/forest_covariates.csv")

## Examine covariates
covariates %>%
	dplyr::select(treecover, edge_density_forest, patch_density_forest, shape_index_forest) %>%
	ggpairs()

# Treecover
covariates %>%
  pivot_longer(
    cols = c(treecover, edge_density_forest, patch_density_forest, shape_index_forest),
    names_to = "variable",
    values_to = "value") %>%
	filter(variable == "treecover") %>%
  
  ggplot(aes(x = year, y = value)) +
  geom_line() +
  geom_point(size = 1.5) +
  facet_wrap(~station, scales = "fixed") +
  labs(
    x = "Year",
    y = "Value",
    title = "Treecover"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.spacing = unit(0.5, "lines"))

# Edge Density
covariates %>%
	pivot_longer(
		cols = c(treecover, edge_density_forest, patch_density_forest, shape_index_forest),
		names_to = "variable",
		values_to = "value") %>%
	filter(variable == "edge_density_forest") %>%
	
	ggplot(aes(x = year, y = value)) +
	geom_line() +
	geom_point(size = 1.5) +
	facet_wrap(~station, scales = "fixed") +
	labs(
		x = "Year",
		y = "Value",
		title = "Edge Density"
	) +
	theme_bw() +
	theme(
		strip.text = element_text(size = 10),
		axis.text.x = element_text(angle = 45, hjust = 1),
		panel.spacing = unit(0.5, "lines"))

# Patch Density
covariates %>%
	pivot_longer(
		cols = c(treecover, edge_density_forest, patch_density_forest, shape_index_forest),
		names_to = "variable",
		values_to = "value") %>%
	filter(variable == "patch_density_forest") %>%
	
	ggplot(aes(x = year, y = value)) +
	geom_line() +
	geom_point(size = 1.5) +
	facet_wrap(~station, scales = "fixed") +
	labs(
		x = "Year",
		y = "Value",
		title = "Edge Density"
	) +
	theme_bw() +
	theme(
		strip.text = element_text(size = 10),
		axis.text.x = element_text(angle = 45, hjust = 1),
		panel.spacing = unit(0.5, "lines"))

# Shape Index
covariates %>%
	pivot_longer(
		cols = c(treecover, edge_density_forest, patch_density_forest, shape_index_forest),
		names_to = "variable",
		values_to = "value") %>%
	filter(variable == "shape_index_forest") %>%
	
	ggplot(aes(x = year, y = value)) +
	geom_line() +
	geom_point(size = 1.5) +
	facet_wrap(~station, scales = "fixed") +
	labs(
		x = "Year",
		y = "Value",
		title = "Shape Index"
	) +
	theme_bw() +
	theme(
		strip.text = element_text(size = 10),
		axis.text.x = element_text(angle = 45, hjust = 1),
		panel.spacing = unit(0.5, "lines"))





