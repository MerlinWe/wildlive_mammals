library(abind)
library(dplyr)
library(tidyr)

read_csv("/Users/merlin/Documents/Senckenberg/WildLive/Mammals/Code/Data/forest_covariates.csv") %>%
  mutate(
    treecover_z = as.numeric(scale(treecover)),
    edge_density_z = as.numeric(scale(edge_density_forest)),
    patch_density_z = as.numeric(scale(patch_density_forest)),
    shape_index_z = as.numeric(scale(shape_index_forest))
  ) %>%
	ggplot(aes(x = year, y = edge_density_z)) +
	geom_line(color = "forestgreen", linewidth = 1) +
	geom_point(color = "darkgreen") +
	facet_wrap(~ station, scales = "fixed") +
	theme_bw() +
	labs(
		title = "Edge Density over Time",
		x = "Year",
		y = "Edge Density (m/ha)"
	)

read_csv("/Users/merlin/Documents/Senckenberg/WildLive/Mammals/Code/Data/forest_covariates.csv") %>%
  mutate(
    treecover_z = as.numeric(scale(treecover)),
    edge_density_z = as.numeric(scale(edge_density_forest)),
    patch_density_z = as.numeric(scale(patch_density_forest)),
    shape_index_z = as.numeric(scale(shape_index_forest))
  ) %>%
	ggplot(aes(x = year, y = treecover_z)) +
	geom_line(color = "sienna", linewidth = 1) +
	geom_point(color = "darkred") +
	facet_wrap(~ station, scales = "fixed") +
	theme_bw() +
	labs(
		title = "Predicted Tree Cover over Time",
		x = "Year",
		y = "Tree Cover (%)"
	)

# Collapse months to presence/absence per station-year
site_year_species <- apply(detection_array, c(1, 2), function(x) as.integer(any(x > 0)))

# Convert to data frame
site_year_df <- as.data.frame(site_year_species)
rownames(site_year_df) <- dimnames(detection_array)[[1]]
colnames(site_year_df) <- dimnames(detection_array)[[2]]

# Long format for merging
site_year_long <- site_year_df %>%
	tibble::rownames_to_column("species") %>%
	pivot_longer(-species, names_to = "station_year", values_to = "detection") %>%
	filter(detection == 1)

site_trait_matrix <- site_year_long %>%
	left_join(species_traits, by = c("species" = "accepted_bin")) %>%
	filter(!is.na(home_range_km2))  # Or any trait filter to ensure trait data coverage

library(FD)

# Create a site-year × species matrix again, but this time we keep species with traits
site_species_matrix <- site_trait_matrix %>%
	select(station_year, species) %>%
	mutate(presence = 1) %>%
	pivot_wider(names_from = species, values_from = presence, values_fill = 0)

# Match trait data for only these species
trait_matrix <- species_traits %>%
	filter(accepted_bin %in% colnames(site_species_matrix)[-1]) %>%
	column_to_rownames("accepted_bin")  # Must match species names

# Extract species names
species_in_traits <- rownames(trait_matrix)
species_in_matrix <- colnames(site_species_matrix)[-1]

# Ensure matching names
shared_species <- intersect(species_in_traits, species_in_matrix)

# Filter and reorder both
trait_matrix_aligned <- trait_matrix[shared_species, ] %>%
	.[order(rownames(.)), ]

site_species_matrix_aligned <- site_species_matrix %>%
	select(station_year, all_of(shared_species)) %>%
	select(station_year, sort(shared_species))

# Final check (optional)
stopifnot(identical(rownames(trait_matrix_aligned), colnames(site_species_matrix_aligned)[-1]))

# Compute functional composition
fd_results <- dbFD(
	x = trait_matrix_aligned,
	a = site_species_matrix_aligned[,-1],
	calc.CWM = TRUE, 
	corr = "cailliez"  # Use correction if traits are mixed or non-Euclidean
)




forest_change <- read_csv("/Users/merlin/Documents/Senckenberg/WildLive/Mammals/Code/Data/forest_covariates.csv") %>%
	group_by(station) %>%
	arrange(station, year) %>%
	filter(year != 2023) %>%
	summarise(change = last(treecover) - first(treecover)) %>%
	mutate(change_class = case_when(
		change > 5 ~ "Gain",
		change < -5 ~ "Loss",
		TRUE ~ "Stable"
	))

fdis_df <- data.frame(
	station_year = site_species_matrix$station_year,
	FDis = fd_results$FDis
) %>%
	separate(station_year, into = c("station", "year"), sep = "_", convert = TRUE) %>%
	left_join(forest_change, by = "station")

ggplot(fdis_df, aes(x = year, y = FDis, group = station, color = change_class)) +
	geom_line(alpha = 0.5) +
	theme_minimal() +
	labs(x = "Year", y = "Functional Dispersion", color = "Forest change")

glimpse(fdis_df)

ggplot(fdis_df, aes(x = change_class, y = FDis, fill = change_class)) +
	geom_boxplot(alpha = 0.7) +
	geom_jitter(width = 0.2, alpha = 0.4) +
	theme_minimal() +
	labs(x = "Forest Cover Change", y = "Functional Dispersion")


library(broom)
station_slopes <- fdis_df %>%
	group_by(station) %>%
	filter(n() >= 3) %>%
	do(tidy(lm(FDis ~ year, data = .))) %>%
	filter(term == "year") %>%
	left_join(forest_change, by = "station")

ggplot(station_slopes, aes(x = change_class, y = estimate, fill = change_class)) +
	geom_boxplot() +
	labs(y = "FDis temporal trend (slope)", x = "Forest change class")




library(vegan)
library(tidyverse)

# Step 1: Trait groups — e.g. small vs large home range
trait_groups <- species_traits %>%
  mutate(group = case_when(
    home_range_km2 <= quantile(home_range_km2, 0.25, na.rm = TRUE) ~ "SmallHR",
    home_range_km2 >= quantile(home_range_km2, 0.75, na.rm = TRUE) ~ "LargeHR",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(group))

# Step 2: Build binary presence–absence matrix per site-year
site_pa_long <- apply(detection_array, c(1, 2), function(x) as.integer(any(x > 0))) %>%
  as.data.frame() %>%
  rownames_to_column("species") %>%
  pivot_longer(-species, names_to = "station_year", values_to = "presence") %>%
  filter(presence == 1)

# Join trait group
site_trait_long <- site_pa_long %>%
  inner_join(trait_groups, by = c("species" = "accepted_bin")) %>%
  separate(station_year, into = c("station", "year"), sep = "_", convert = TRUE)

# Step 3: Calculate turnover for each group and station
# Step 3: Calculate turnover for each group and station
trait_turnover <- site_trait_long %>%
  group_by(group, station, year) %>%
  summarise(species = list(unique(species)), .groups = "drop") %>%
  arrange(group, station, year) %>%
  group_by(group, station) %>%
  mutate(prev_species = lag(species)) %>%
  filter(!is.na(prev_species)) %>%
  rowwise() %>%
  mutate(
    # Get only species in this trait group
    group_species_pool = list(trait_groups$accepted_bin[trait_groups$group == group]),
    jaccard = vegan::vegdist(rbind(
      as.integer(unlist(group_species_pool) %in% species),
      as.integer(unlist(group_species_pool) %in% prev_species)
    ), method = "jaccard")[1]
  ) %>%
  ungroup()

# Step 4: Add forest change class
turnover_forest <- trait_turnover %>%
  left_join(forest_change, by = "station")

# Step 5: Visualize
ggplot(turnover_forest, aes(x = change_class, y = jaccard, fill = change_class)) +
  geom_boxplot(alpha = 0.7) +
  facet_wrap(~ group) +
  labs(
    y = "Jaccard Turnover (year-to-year)",
    x = "Forest Cover Change",
    title = "Trait-group-specific turnover across forest change gradients"
  ) +
  theme_minimal()











