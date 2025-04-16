rm(list=ls()); gc()

library(patchwork)
library(ggh4x)
library(tidyverse)
library(spOccupancy)
library(broom)
library(scales)


species <- read_csv()
camtraps <- read_csv("/Users/serpent/Desktop/Paper/camtraps_clean.csv")
camop <- as.matrix(read_csv("/Users/serpent/Desktop/Paper/camop_problem.csv"))
covariates <- read_csv("/Users/serpent/Desktop/Paper/covariates.csv")

species <- species %>%
	filter(Station %in% camtraps$Station) %>%
	filter(Category %in% c("artiodactyla", "carnivora", "marsupialia", "perissodactyla", "primates", "rodentia", "xenarthra")) %>%
	filter(DateTimeOriginal >= "2017-01-10 17:01:26") %>%
	filter(DateTimeOriginal <= "2023-10-18 06:32:56")

# Integrate pantheria traits, see which traits predict which response
pantheria <- read.table(
	file = "/Users/serpent/Documents/MSc/Computational Methods in Biodiversity Research (BIOS0002)/Data Science (Lectures)/PanTHERIA (W2)/PanTHERIA_1-0_WR05_Aug2008.txt",
	header = TRUE, sep = "\t", na.strings = c("-999", "-999.00")) %>%
	as_tibble() %>%
	select(
		binomial = MSW05_Binomial,
		order = MSW05_Order,
		adult_body_mass_g = X5.1_AdultBodyMass_g,
		home_range_km2 = X22.1_HomeRange_km2,
		activity_cycle = X1.1_ActivityCycle,
		diet_breadth = X6.1_DietBreadth,
		habitat_breadth = X12.1_HabitatBreadth,
		terrestriality = X12.2_Terrestriality,
		trophic_level = X6.2_TrophicLevel)

# update outdated species names 
species <- species %>%
	mutate(accepted_bin = case_when(
		accepted_bin == "dicotyles_tajacu" ~ "pecari_tajacu",
		accepted_bin == "herpaliurus_yagouaroundi" ~ "puma_yagouaroundi",
		accepted_bin == "sapajus_apella" ~ "cebus_apella",
		accepted_bin == "myrmecophaga_tridaytyla" ~ "myrmecophaga_tridactyla",
		TRUE ~ accepted_bin
	))

# Clean and join PanTHERIA traits with your species list
species_traits <- pantheria %>%
	mutate(accepted_bin = str_to_lower(str_replace_all(binomial, " ", "_"))) %>%
	filter(accepted_bin %in% unique(species$accepted_bin)) %>%
	select(accepted_bin, order, adult_body_mass_g, home_range_km2,
				 activity_cycle, diet_breadth, habitat_breadth,
				 terrestriality, trophic_level)


# Analysis 

create_detection_array <- function(species_data, species_min = 10, K = 12) {
	library(dplyr)
	library(tidyr)
	library(lubridate)
	
	species_data <- species_data %>%
		mutate(Date = as.Date(DateTimeOriginal),
					 year = year(Date),
					 site_year = paste(Station, year, sep = "_"))
	
	# Step 1: Count total detections per species
	species_counts <- species_data %>%
		count(accepted_bin) %>%
		dplyr::filter(n >= 10)
	
	species_list <- species_counts$accepted_bin
	
	# Step 2: Filter to valid species
	species_data <- species_data %>%
		filter(accepted_bin %in% species_list)
	
	sites <- sort(unique(species_data$site_year))
	S <- length(species_list)
	J <- length(sites)
	
	y_array <- array(NA, dim = c(J, K, S))
	dimnames(y_array)[[1]] <- sites
	dimnames(y_array)[[3]] <- species_list
	
	for (s in seq_along(species_list)) {
		sp <- species_list[s]
		sp_data <- species_data %>%
			filter(accepted_bin == sp) %>%
			mutate(month = month(Date)) %>%
			count(site_year, month) %>%
			mutate(detected = 1)
		
		sp_sites <- unique(sp_data$site_year)
		
		for (i in seq_along(sites)) {
			site <- sites[i]
			if (!(site %in% sp_sites)) next
			for (k in 1:K) {
				val <- sp_data %>%
					filter(site_year == site, month == k) %>%
					pull(detected)
				y_array[i, k, s] <- ifelse(length(val) > 0, 1, 0)
			}
		}
	}
	
	return(y_array)
}

detection_array <- create_detection_array(species, camop)

detection_array <- aperm(detection_array, perm = c(3, 1, 2))
dim(detection_array)      

site_years <- dimnames(detection_array)[[2]]

site_covs <- covariates %>%
	mutate(site_year = paste(station, year, sep = "_")) %>%
	filter(site_year %in% site_years) %>%
	mutate(
		treecover_z = as.numeric(scale(treecover)),
		aggregation_z = as.numeric(scale(aggregation_forest))
	) %>%
	arrange(match(site_year, site_years)) %>%
	select(treecover_z, aggregation_z, agriculture) %>%
	as.data.frame()

# FIXED: match rownames to the new site_years
rownames(site_covs) <- site_years

# Optional checks (should return TRUE)
nrow(site_covs) == dim(detection_array)[2]
all.equal(rownames(site_covs), dimnames(detection_array)[[2]])

# Clean up detection array
detection_array[is.na(detection_array)] <- 0
any(is.na(detection_array))  # Should return FALSE

# name months dimension 
dimnames(detection_array)[[3]] <- paste0("Month", 1:dim(detection_array)[3])

# Fit the model
ms_model <- msPGOcc(
	occ.formula = ~ treecover_z + aggregation_z + agriculture,
	det.formula = ~ 1,
	data = list(
		y = detection_array,
		occ.covs = site_covs
	),
	n.samples = 15000,
	n.burn = 2000,
	n.thin = 10,
	n.chains = 3,
	verbose = TRUE
)

summary(ms_model, level = 'community')
summary(ms_model, level = 'species')
summary(ms_model, level = "community")$beta.comm
summary(ms_model, level = "species")$beta
psi <- fitted(ms_model)

# Extract posterior summary stats for coefficients
beta_samples <- ms_model$beta.samples
beta_summary <- as.data.frame(summary(beta_samples)$statistics)
beta_ci <- as.data.frame(summary(beta_samples)$quantiles)

# Prepare posterior summaries and join with trait data
coef_traits <- ms_model$beta.samples %>%
	summary() %>%
	{
		bind_cols(
			as.data.frame(.$statistics) %>% select(Mean = Mean),
			as.data.frame(.$quantiles) %>% select(`2.5%`, `97.5%`)
		)
	} %>%
	rownames_to_column("term") %>%
	separate(term, into = c("covariate", "species"), sep = "-", extra = "merge") %>%
	filter(covariate %in% c("treecover_z", "aggregation_z", "agriculture")) %>%
	left_join(species_traits, by = c("species" = "accepted_bin")) %>%
	mutate(
		log_mass = log(adult_body_mass_g),
		log_range = log(home_range_km2),
		diet_breadth = as.numeric(diet_breadth),
		# Recode activity cycle: 1 = nocturnal, 2 = diurnal, 3 = cathemeral
		activity_cycle_coded = recode(activity_cycle, `1` = -1, `3` = 0, `2` = 1),
		# Recode habitat breadth: 1 = specialist, 2 = intermediate, 3 = generalist (rescaled to -1 to 1)
		habitat_breadth_coded = rescale(habitat_breadth, to = c(-1, 1))
	) %>%
	select(covariate, species, Mean, `2.5%`, `97.5%`,
				 log_mass, log_range, diet_breadth,
				 activity_cycle_coded, habitat_breadth_coded)

# Define traits to test
trait_vars <- c("log_mass", "log_range", "diet_breadth",
								"activity_cycle_coded", "habitat_breadth_coded")

# Run trait ~ coefficient regressions across covariates
results <- map_dfr(unique(coef_traits$covariate), function(cov) {
	dat <- coef_traits %>% filter(covariate == cov)
	
	map_dfr(trait_vars, function(trait) {
		formula <- as.formula(paste0("Mean ~ ", trait))
		mod <- lm(formula, data = dat)
		tidy(mod) %>%
			filter(term != "(Intercept)") %>%
			mutate(covariate = cov, trait = trait)
	})
})

# Clean up for plotting
results <- results %>%
	distinct(covariate, trait, .keep_all = TRUE) %>%
	mutate(trait = recode(trait,
												"log_mass" = "Body mass",
												"log_range" = "Home range size",
												"diet_breadth" = "Diet breadth",
												"activity_cycle_coded" = "Diurnality",
												"habitat_breadth_coded" = "Habitat breadth"))

results$covariate <- factor(results$covariate, levels = c("agriculture", "treecover_z", "aggregation_z"))


beta_comm <- ms_model$beta.comm.samples %>%
	summary() %>%
	{
		bind_cols(
			as.data.frame(.$statistics) %>% select(Mean = Mean),
			as.data.frame(.$quantiles) %>% select(`2.5%`, `97.5%`)
		)
	} %>%
	rownames_to_column("term") %>%
	filter(term != "(Intercept)") %>%
	mutate(term = recode(term,
											 "aggregation_z" = "Fragmentation",
											 "agriculture" = "Agriculture",
											 "treecover_z" = "Forestcover"))

beta_comm$term <- factor(beta_comm$term, levels = c("Agriculture", "Forestcover", "Fragmentation"))

# Community level plot 
comm_plot <- ggplot(beta_comm, aes(x = term, y = Mean, ymin = `2.5%`, ymax = `97.5%`)) +
	geom_pointrange(aes(colour = term), size = 1.2) +
	geom_hline(yintercept = 0, linetype = "dashed") +
	labs(title = "(a.) Community-level responses",
			 x = NULL, y = "Effect size (logit scale)") +
	scale_colour_viridis_d(option = "magma", end = .8) +
	theme_bw() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1),
				text = element_text(size = 10),
				legend.position = "none")

# Individual trait plots per covariate
traits_treecover <- results %>%
	filter(covariate == "treecover_z") %>%
	ggplot(aes(x = trait, y = estimate, 
						 ymin = estimate - std.error, ymax = estimate + std.error,
						 color = trait)) +
	geom_hline(yintercept = 0, linetype = "dashed") +
	geom_pointrange(position = position_dodge(width = 0.5), size = 0.8) +
	scale_color_viridis_d() +
	labs(title = "(b.) Trait responses - Forestcover", 
			 x = NULL, y = "Effect size (slope estimate ± SE)") +
	theme_bw() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1),
				text = element_text(size = 10),
				legend.position = "none")

traits_fragmentation <- results %>%
	filter(covariate == "aggregation_z") %>%
	ggplot(aes(x = trait, y = estimate, 
						 ymin = estimate - std.error, ymax = estimate + std.error,
						 color = trait)) +
	geom_hline(yintercept = 0, linetype = "dashed") +
	geom_pointrange(position = position_dodge(width = 0.5), size = 0.8) +
	scale_color_viridis_d() +
	labs(title = "(c.) Trait responses - Fragmentation", 
			 x = NULL, y = "Effect size (slope estimate ± SE)") +
	theme_bw() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1),
				text = element_text(size = 10),
				legend.position = "none")

traits_agriculture <- results %>%
	filter(covariate == "agriculture") %>%
	ggplot(aes(x = trait, y = estimate, 
						 ymin = estimate - std.error, ymax = estimate + std.error,
						 color = trait)) +
	geom_hline(yintercept = 0, linetype = "dashed") +
	geom_pointrange(position = position_dodge(width = 0.5), size = 0.8) +
	scale_color_viridis_d() +
	labs(title = "(d.) Trait responses - Agriculture", 
			 x = NULL, y = "Effect size (slope estimate ± SE)") +
	theme_bw() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1),
				text = element_text(size = 10),
				legend.position = "none")

effects_plot <- (comm_plot | traits_treecover) /
	(traits_fragmentation | traits_agriculture)

ggsave("/Users/serpent/Desktop/Paper/traits_effects.png", 
			 effects_plot, 
			 width = 180,
			 height = 160, 
			 units = "mm",
			 dpi = 600)



## PDP analyis

# Define covariate sequences
tree_seq <- seq(-2, 2, length.out = 100)
frag_seq <- seq(-2, 2, length.out = 100)
agri_seq <- c(0, 1)

# Extract species names from coefficient column names
species_names <- gsub("treecover_z-", "", grep("^treecover_z-", colnames(ms_model$beta.samples), value = TRUE))

# Define focal species and get their indices
focal_species <- c("panthera_onca", "tayassu_pecari", "mazama_gouazoubira", "cerdocyon_thous")
focal_ids <- match(focal_species, species_names)
stopifnot(!any(is.na(focal_ids)))  # safety check

# Define a robust PDP function
get_pdp <- function(predictor, predictor_seq, species_ids, species_names, ms_model) {
	n_vals <- length(predictor_seq)
	
	cov_data <- data.frame(
		treecover_z   = if (predictor == "treecover_z") predictor_seq else rep(0, n_vals),
		aggregation_z = if (predictor == "aggregation_z") predictor_seq else rep(0, n_vals),
		agriculture   = if (predictor == "agriculture") predictor_seq else rep(0, n_vals)
	)
	
	X.0 <- model.matrix(~ treecover_z + aggregation_z + agriculture, data = cov_data)
	pred <- predict(ms_model, X.0 = X.0, type = "occupancy", ignore.RE = FALSE)
	
	pdp_df <- map_dfr(seq_along(species_ids), function(i) {
		psi_samples <- pred$psi.0.samples[, species_ids[i], ]
		data.frame(
			species = species_names[species_ids[i]],
			predictor = predictor,
			value = predictor_seq,
			mean = apply(psi_samples, 2, mean),
			lower = apply(psi_samples, 2, quantile, probs = 0.025),
			upper = apply(psi_samples, 2, quantile, probs = 0.975)
		)
	})
	
	return(pdp_df)
}

# Generate PDPs for focal species
pdp_tree_focal <- get_pdp("treecover_z", tree_seq, focal_ids, species_names, ms_model)
pdp_frag_focal <- get_pdp("aggregation_z", frag_seq, focal_ids, species_names, ms_model)
pdp_agri_focal <- get_pdp("agriculture", agri_seq, focal_ids, species_names, ms_model)

pdp_focal_all <- bind_rows(pdp_tree_focal, pdp_frag_focal, pdp_agri_focal)

# Plot focal PDPs

# Mapping scientific to common names
species_labels <- c(
	panthera_onca = "Jaguar",
	tayassu_pecari = "White-lipped Peccary",
	mazama_gouazoubira = "Brown Brocket",
	cerdocyon_thous = "Crab-eating Fox"
)
# Order predictors for facets
pdp_focal_all$predictor <- factor(
	pdp_focal_all$predictor,
	levels = c("agriculture", "treecover_z", "aggregation_z")
)

# Labels
predictor_labels <- c(
	agriculture = "Agriculture (binary)",
	treecover_z = "Forest cover (z-score)",
	aggregation_z = "Fragmentation (z-score)"
)

# Recode species column
pdp_focal_all$species <- recode(pdp_focal_all$species, !!!species_labels)

pdp_focal_plot <- ggplot(pdp_focal_all, aes(x = value, y = mean, ymin = lower, ymax = upper, fill = species)) +
	geom_ribbon(alpha = 0.2, colour = NA) +
	geom_line(aes(color = species), linewidth = 1) +
	ggh4x::facet_nested(
		rows = vars(species),
		cols = vars(predictor),
		scales = "free",
		labeller = labeller(predictor = predictor_labels)
	) +
	scale_color_viridis_d(end = .85) +
	scale_fill_viridis_d(end = .85) +
	labs(
		x = "Covariate value", y = "Predicted occupancy probability") +
	theme_bw() +
	theme(
		legend.position = "none",
		strip.text = element_text(size = 10),
		text = element_text(size = 11)
	)

ggsave("/Users/serpent/Desktop/Paper/pdp_focal.png", 
			 pdp_focal_plot, 
			 width = 180,
			 height = 180, 
			 units = "mm",
			 dpi = 600)



# Generate PDPs for ALL species (for supplementary)
all_ids <- seq_along(species_names)
pdp_all_tree <- get_pdp("treecover_z", tree_seq, all_ids, species_names, ms_model)
pdp_all_frag <- get_pdp("aggregation_z", frag_seq, all_ids, species_names, ms_model)
pdp_all_agri <- get_pdp("agriculture", agri_seq, all_ids, species_names, ms_model)

pdp_all <- bind_rows(pdp_all_tree, pdp_all_frag, pdp_all_agri)

# Rename predictors for display
pdp_all <- pdp_all %>%
	mutate(
		species = tools::toTitleCase(gsub("_", " ", species)),
		predictor = recode(predictor,
											 treecover_z   = "Forest cover",
											 aggregation_z = "Fragmentation",
											 agriculture   = "Agriculture"
		),
		facet_label = paste0(species_clean, "\n(", predictor, ")")
	)

# Create combined facet label
pdp_all <- pdp_all %>%
	mutate(facet_label = paste(species, "\n(", predictor, ")", sep = ""))

# Plot with facet_wrap in a wide layout
pdp_all_plot <- ggplot(pdp_all, aes(x = value, y = mean, ymin = lower, ymax = upper, fill = predictor, colour = predictor)) +
	geom_ribbon(alpha = 0.2) +
	geom_line() +
	scale_fill_viridis_d(end = .8) +
	scale_colour_viridis_d(end = .8) +
	facet_wrap(~ facet_label, nrow = 9, ncol = 7, scales = "free") +  # Adjust nrow to fit your figure size
	labs(
		x = "Covariate value", y = "Predicted occupancy probability"
	) +
	theme_bw() +
	theme(
		text = element_text(size = 8),
		strip.text = element_text(size = 8),
		axis.text = element_text(size = 7),
		panel.spacing = unit(0.2, "lines"),
		legend.position = "none")

ggsave("/Users/serpent/Desktop/Paper/pdp_all.png", 
			 pdp_all_plot, 
			 width = 270,
			 height = 320, 
			 units = "mm",
			 dpi = 600)


# Clean column names
trait_table <- results %>%
	select(Trait = trait,
				 Covariate = covariate,
				 Estimate = estimate,
				 SE = std.error,
				 p = p.value) %>%
	mutate(across(where(is.numeric), round, 3)) %>%
	mutate(
		Covariate = recode(Covariate,
											 treecover_z = "Forest cover",
											 aggregation_z = "Fragmentation",
											 agriculture = "Agriculture")
	) %>%
	mutate(sig = ifelse(p < 0.05, TRUE, FALSE))

# Combine Estimate ± SE into one column
trait_table <- trait_table %>%
	mutate(`Estimate (± SE)` = sprintf("%.2f ± %.2f", Estimate, SE)) %>%
	select(Covariate, Trait, `Estimate (± SE)`, p, sig)

# Pivot wider to have covariates as columns
trait_table_wide <- trait_table %>%
	select(-sig) %>%  # Drop sig if not used in export
	pivot_wider(
		names_from = Covariate,
		values_from = c(`Estimate (± SE)`, p),
		names_sep = ": "
	) %>%
	relocate(Trait)

write_csv(trait_table_wide, "/Users/serpent/Desktop/Paper/trait_table.csv")


##### More tables #####

trait_table <- coef_traits %>%
	select(species, log_mass, log_range, diet_breadth, 
				 activity_cycle_coded, habitat_breadth_coded) %>%
	distinct() %>%
	mutate(
		species = tools::toTitleCase(gsub("_", " ", species))
	) %>%
	rename(
		`Species`             = species,
		`Log Body Mass`       = log_mass,
		`Log Range Size`      = log_range,
		`Diet Breadth`        = diet_breadth,
		`Activity Cycle`      = activity_cycle_coded,
		`Habitat Breadth`     = habitat_breadth_coded
	)

write_csv(trait_table, "/Users/serpent/Desktop/Paper/trait_table.csv")



# diagnostics
# 1. Extract Rhat and ESS as named vectors
rhat_vec <- ms_model$rhat$beta
ess_vec  <- ms_model$ESS$beta

length(rhat_vec)
names(rhat_vec)[1:5]

# 2. Build a dataframe from the beta samples summary
beta_df <- as.data.frame(summary(beta_samples)$statistics) %>%
	rownames_to_column("param") %>%
	filter(str_detect(param, "^(treecover_z|aggregation_z|agriculture)-")) %>%
	mutate(
		covariate   = str_extract(param, "^[^\\-]+"),
		species_raw = str_extract(param, "(?<=-).*"),
		species     = tools::toTitleCase(gsub("_", " ", species_raw))
	)

# 3. Turn diagnostics into a dataframe with param names
param_names <- beta_df$param

rhat_filtered <- rhat_vec[param_names]
ess_filtered  <- ess_vec[param_names]

length(rhat_filtered) == nrow(beta_df)  # should now be TRUE

diagnostics_table <- beta_df %>%
	mutate(
		Rhat = rhat_filtered,
		ESS  = ess_filtered,
		converged = Rhat < 1.1,
		ESS_gt_1000 = ESS > 1000
	) %>%
	select(species, covariate, Mean, SD, Rhat, ESS, converged, ESS_gt_1000)

all(names(rhat_filtered) == beta_df$param)  # should be TRUE

diagnostics_clean <- diagnostics_table %>%
	rename(
		Species     = species,
		Covariate   = covariate,
		Mean        = Mean,
		SD          = SD,
		Rhat        = Rhat,
		ESS         = ESS,
		Converged   = converged,
		ESS_gt_1000 = ESS_gt_1000
	) %>%
	mutate(
		Covariate = recode(Covariate,
											 "treecover_z"   = "Forest cover",
											 "aggregation_z" = "Fragmentation",
											 "agriculture"   = "Agricultural land"
		),
		Mean = round(Mean, 3),
		SD   = round(SD, 3),
		Rhat = round(Rhat, 3),
		ESS  = round(ESS, 0)
	) %>%
	select(-Converged, -ESS_gt_1000)

write_csv(diagnostics_clean, "/Users/serpent/Desktop/Paper/diagnostics_table.csv")
