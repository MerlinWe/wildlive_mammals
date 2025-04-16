
# Load required libraries
library(dplyr)
library(tidyr)
library(stringr)

# 1. Define species metadata
common_names <- c(
  "cerdocyon_thous" = "Crab-eating Fox", "cuniculus_paca" = "Paca",
  "dasyprocta_punctata" = "Central American Agouti", "dasypus_novemcinctus" = "Nine-banded Armadillo",
  "pecari_tajacu" = "Collared Peccary", "didelphis_marsupialis" = "Common Opossum",
  "eira_barbara" = "Tayra", "euphractus_sexcinctus" = "Six-banded Armadillo",
  "puma_yagouaroundi" = "Jaguarundi", "hydrochoerus_hydrochaeris" = "Capybara",
  "leopardus_pardalis" = "Ocelot", "leopardus_wiedii" = "Margay",
  "mazama_gouazoubira" = "Brown Brocket", "nasua_nasua" = "South American Coati",
  "panthera_onca" = "Jaguar", "procyon_cancrivorus" = "Crab-eating Raccoon",
  "puma_concolor" = "Puma", "cebus_apella" = "Brown Capuchin",
  "sciurus_ignitus" = "Bolivian Squirrel", "tamandua_tetradactyla" = "Southern Tamandua",
  "tapirus_terrestris" = "Lowland Tapir", "tayassu_pecari" = "White-lipped Peccary",
  "myrmecophaga_tridactyla" = "Giant Anteater", "coendou_prehensilis" = "Brazilian Porcupine"
)
species_metadata <- tibble(latin_name = names(common_names), common_name = unname(common_names))

# 2. Summarise detections
detections <- species %>%
  group_by(accepted_bin) %>%
  summarise(
    n_captures = n(), n_stations = n_distinct(Station),
    capture_rate = n() / n_distinct(Station), .groups = "drop"
  )

# 3. Occupancy & detection summaries
naive_psi <- rowMeans(ms_model$y > 0)
psi_means <- apply(ms_model$psi.samples, 2, mean)
det_naive <- apply(ms_model$y, 1, function(x) mean(x > 0, na.rm = TRUE))

# 4. Clean slope summaries for forest cover
slope_summary_clean <- slope_summary %>%
  rename(latin_name = species, beta_mean = mean, beta_lower = lower, beta_upper = upper) %>%
  mutate(latin_name = gsub(" ", "_", tolower(latin_name)))

# 5. Merge all metadata
summary_table_fmt <- detections %>%
  rename(latin_name = accepted_bin) %>%
  left_join(species_metadata, by = "latin_name") %>%
  left_join(tibble(latin_name = ms_model$sp.names, mean_psi = psi_means, naive_psi = naive_psi, naive_p = det_naive), by = "latin_name") %>%
  left_join(slope_summary_clean, by = "latin_name") %>%
  mutate(across(where(is.numeric), round, 3)) %>%
  mutate(latin_name = gsub("_", " ", tools::toTitleCase(latin_name)))

# 6. Slope summaries for fragmentation and agriculture
summarise_posterior <- function(beta_df, prefix) {
  species <- gsub(".*?-", "", names(beta_df))
  map_dfr(seq_along(species), function(i) {
    post <- beta_df[[i]]
    tibble(
      latin_name = tools::toTitleCase(gsub("_", " ", species[i])),
      !!paste0(prefix, "_mean") := mean(post),
      !!paste0(prefix, "_lower") := quantile(post, 0.025),
      !!paste0(prefix, "_upper") := quantile(post, 0.975)
    )
  })
}

frag_summary <- summarise_posterior(as.data.frame(ms_model$beta.samples) %>% select(starts_with("aggregation_z-")), "beta_frag")
agri_summary <- summarise_posterior(as.data.frame(ms_model$beta.samples) %>% select(starts_with("agriculture-")), "beta_agri")

# 7. Harmonize species case
fix_species_case <- function(name) {
  parts <- strsplit(name, " ")[[1]]
  if (length(parts) == 2) paste(parts[1], tolower(parts[2])) else name
}
frag_summary <- frag_summary %>% mutate(latin_name = sapply(latin_name, fix_species_case))
agri_summary <- agri_summary %>% mutate(latin_name = sapply(latin_name, fix_species_case))

# 8. Join all and add significance flags
add_sig_flag <- function(df, mean_col, lower_col, upper_col, new_col) {
  df %>%
    mutate(
      !!new_col := case_when(
        !!sym(lower_col) > 0 ~ "+",
        !!sym(upper_col) < 0 ~ "−",
        TRUE ~ ""
      )
    )
}
summary_table_full <- summary_table_fmt %>%
  left_join(add_sig_flag(slope_summary_clean, "beta_mean", "beta_lower", "beta_upper", "sig_forest"), by = "latin_name") %>%
  left_join(add_sig_flag(frag_summary, "beta_frag_mean", "beta_frag_lower", "beta_frag_upper", "sig_frag"), by = "latin_name") %>%
  left_join(add_sig_flag(agri_summary, "beta_agri_mean", "beta_agri_lower", "beta_agri_upper", "sig_agri"), by = "latin_name")

# 9. Format final table
format_estimate <- function(mean, lower, upper, sig) {
  ifelse(is.na(mean), NA, sprintf("%.2f [%.2f, %.2f]%s", mean, lower, upper, sig))
}
summary_table_final <- summary_table_full %>%
  mutate(
    `β Forest cover`  = format_estimate(beta_mean, beta_lower, beta_upper, sig_forest),
    `β Fragmentation` = format_estimate(beta_frag_mean, beta_frag_lower, beta_frag_upper, sig_frag),
    `β Agriculture`   = format_estimate(beta_agri_mean, beta_agri_lower, beta_agri_upper, sig_agri)
  ) %>%
  select(
    `Common name` = common_name, `Latin name` = latin_name,
    `Detections` = n_captures, `Stations` = n_stations, `Capture rate` = capture_rate,
    `Naïve Ψ` = naive_psi, `Modelled Ψ` = mean_psi, `Naïve p` = naive_p,
    `β Forest cover`, `β Fragmentation`, `β Agriculture`
  )

# Optional export for Word/Google Docs
# write.csv(summary_table_final, "summary_table.csv", row.names = FALSE)
# knitr::kable(summary_table_final, format = "markdown")
