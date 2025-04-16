library(dplyr)
library(tidyr)
library(stringr)

# Latin names already in snake_case; convert to title case with underscores replaced
common_names <- c(
	"cerdocyon_thous"           = "Crab-eating Fox",
	"cuniculus_paca"            = "Paca",
	"dasyprocta_punctata"       = "Central American Agouti",
	"dasypus_novemcinctus"      = "Nine-banded Armadillo",
	"pecari_tajacu"             = "Collared Peccary",
	"didelphis_marsupialis"     = "Common Opossum",
	"eira_barbara"              = "Tayra",
	"euphractus_sexcinctus"     = "Six-banded Armadillo",
	"puma_yagouaroundi"         = "Jaguarundi",
	"hydrochoerus_hydrochaeris"= "Capybara",
	"leopardus_pardalis"        = "Ocelot",
	"leopardus_wiedii"          = "Margay",
	"mazama_gouazoubira"        = "Brown Brocket",
	"nasua_nasua"               = "South American Coati",
	"panthera_onca"             = "Jaguar",
	"procyon_cancrivorus"       = "Crab-eating Raccoon",
	"puma_concolor"             = "Puma",
	"cebus_apella"              = "Brown Capuchin",
	"sciurus_ignitus"           = "Bolivian Squirrel",
	"tamandua_tetradactyla"     = "Southern Tamandua",
	"tapirus_terrestris"        = "Lowland Tapir",
	"tayassu_pecari"            = "White-lipped Peccary",
	"myrmecophaga_tridactyla"   = "Giant Anteater",
	"coendou_prehensilis"       = "Brazilian Porcupine")

species_metadata <- tibble(
	latin_name = names(common_names),
	common_name = unname(common_names))

# Assuming capture_data has columns: species (latin name), station, datetime
detections <- species %>%
	group_by(accepted_bin) %>%
	summarise(
		n_captures = n(),
		n_stations = n_distinct(Station),
		capture_rate = n() / n_distinct(Station), 
		.groups = "drop")

# Get posterior summaries for psi and p
psi_summ <- apply(ms_model$psi.samples, 2, mean)  # mean occupancy per species
naive_psi <- rowMeans(ms_model$y > 0)             # naïve occupancy = % sites detected

# Average number of detections per occasion, per species (naive p)
det_raw <- ms_model$y   # [species x site x occasion]
det_naive <- apply(det_raw, 1, function(x) mean(x > 0, na.rm = TRUE))

# You can combine this with the occupancy estimates:
psi_means <- apply(ms_model$psi.samples, 2, mean)  # mean ψ per species


# Extract treecover slopes
beta_tree <- as.data.frame(ms_model$beta.samples) %>%
	select(starts_with("treecover_z-"))

# Clean species names
species_names <- gsub("treecover_z-", "", names(beta_tree))


# Flag whether 95% CI excludes zero (i.e., significant effect)
add_sig_flag <- function(df, mean_col, lower_col, upper_col, new_col = "sig") {
	df %>%
		mutate(
			{{new_col}} := case_when(
				!!sym(lower_col) > 0 ~ "+",
				!!sym(upper_col) < 0 ~ "−",
				TRUE ~ ""
			)
		)
}

# Calculate posterior summaries for each species
slope_summary <- map_dfr(seq_along(species_names), function(i) {
	vals <- beta_tree[[i]]
	data.frame(
		species = species_names[i],
		mean = mean(vals),
		lower = quantile(vals, 0.025),
		upper = quantile(vals, 0.975)
	)
})

# Optional: prettify names
slope_summary$species <- gsub("_", " ", slope_summary$species)

# First: prepare slope summary nicely
slope_summary_clean <- slope_summary %>%
	rename(
		latin_name = species,
		beta_mean = mean,
		beta_lower = lower,
		beta_upper = upper
	) %>%
	mutate(latin_name = gsub(" ", "_", tolower(latin_name)))  # match snake_case

# Add to each slope summary
slope_summary_clean <- add_sig_flag(slope_summary_clean, "beta_mean", "beta_lower", "beta_upper", "sig_forest")

# Combine all data
summary_table <- detections %>%
	rename(latin_name = accepted_bin) %>%
	left_join(species_metadata, by = "latin_name") %>%
	left_join(
		tibble(latin_name = ms_model$sp.names,
					 mean_psi = psi_means,
					 naive_psi = naive_psi,
					 naive_p = det_naive),
		by = "latin_name"
	) %>%
	left_join(slope_summary_clean, by = "latin_name") %>%
	select(common_name, latin_name, n_captures, n_stations, capture_rate,
				 naive_psi, mean_psi, naive_p, beta_mean, beta_lower, beta_upper, sig_forest)

# Round for display
summary_table_fmt <- summary_table %>%
	mutate(across(where(is.numeric), ~ round(., 3)))  %>%
	mutate(
		latin_name = gsub("_", " ", tools::toTitleCase(latin_name))
	)






# Extract posterior slopes for fragmentation and agriculture
beta_frag <- as.data.frame(ms_model$beta.samples) %>%
	select(starts_with("aggregation_z-"))
beta_agri <- as.data.frame(ms_model$beta.samples) %>%
	select(starts_with("agriculture-"))

# Summarize posterior stats 
summarise_posterior <- function(beta_df) {
	species <- gsub(".*?-", "", names(beta_df))
	
	beta_summary <- map_dfr(seq_along(species), function(i) {
		post_samples <- beta_df[[i]]
		tibble(
			latin_name = tools::toTitleCase(gsub("_", " ", species[i])),  # Match summary_table_fmt format
			mean = mean(post_samples),
			lower = quantile(post_samples, 0.025),
			upper = quantile(post_samples, 0.975)
		)
	})
	
	return(beta_summary)
}

# Summarize
frag_summary <- summarise_posterior(beta_frag) %>%
	rename(beta_frag_mean = mean, beta_frag_lower = lower, beta_frag_upper = upper) 

agri_summary <- summarise_posterior(beta_agri) %>%
	rename(beta_agri_mean = mean, beta_agri_lower = lower, beta_agri_upper = upper) 

fix_species_case <- function(name) {
	parts <- strsplit(name, " ")[[1]]
	if (length(parts) == 2) {
		paste(parts[1], tolower(parts[2]))
	} else {
		name
	}
}

# Apply to both summaries
frag_summary <- frag_summary %>%
	mutate(latin_name = sapply(latin_name, fix_species_case))

agri_summary <- agri_summary %>%
	mutate(latin_name = sapply(latin_name, fix_species_case))

frag_summary        <- add_sig_flag(frag_summary, "beta_frag_mean", "beta_frag_lower", "beta_frag_upper", "sig_frag")
agri_summary        <- add_sig_flag(agri_summary, "beta_agri_mean", "beta_agri_lower", "beta_agri_upper", "sig_agri")


summary_table_full <- summary_table_fmt %>%
	left_join(frag_summary, by = "latin_name") %>%
	left_join(agri_summary, by = "latin_name")

summary_table_full


summary_table_final <- summary_table_full %>%
	mutate(
		`β Forest cover`     = format_estimate(beta_mean, beta_lower, beta_upper, sig_forest),
		`β Fragmentation`    = format_estimate(beta_frag_mean, beta_frag_lower, beta_frag_upper, sig_frag),
		`β Agriculture`      = format_estimate(beta_agri_mean, beta_agri_lower, beta_agri_upper, sig_agri)
	) %>%
	select(
		`Common name` = common_name,
		`Latin name` = latin_name,
		`Detections` = n_captures,
		`Stations` = n_stations,
		`Capture rate` = capture_rate,
		`Naïve Ψ` = naive_psi,
		`Modelled Ψ` = mean_psi,
		`Naïve p` = naive_p,
		`β Forest cover`,
		`β Fragmentation`,
		`β Agriculture`
	)

# Format estimate with 95% CI and significance flag
format_estimate <- function(mean, lower, upper, sig) {
	ifelse(
		is.na(mean), NA,
		sprintf("%.2f [%.2f, %.2f]%s", mean, lower, upper, sig)
	)
}

# Assemble final display table
summary_table_final <- summary_table_full %>%
	mutate(
		`β Forest cover` = format_estimate(beta_mean, beta_lower, beta_upper, sig_forest),
		`β Fragmentation` = format_estimate(beta_frag_mean, beta_frag_lower, beta_frag_upper, sig_frag),
		`β Agriculture` = format_estimate(beta_agri_mean, beta_agri_lower, beta_agri_upper, sig_agri)
	) %>%
	select(
		`Common name` = common_name,
		`Latin name` = latin_name,
		`Detections` = n_captures,
		`Stations` = n_stations,
		`Capture rate` = capture_rate,
		`Naïve Ψ` = naive_psi,
		`Modelled Ψ` = mean_psi,
		`Naïve p` = naive_p,
		`β Forest cover`,
		`β Fragmentation`,
		`β Agriculture`
	)

write_csv(summary_table_final, file = "/Users/serpent/Desktop/Paper/summary_table.csv")
