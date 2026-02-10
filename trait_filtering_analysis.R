#########################################################################################################
##### Trait-Filtered Responses of Mammal Communities to Land-Use Change in a Neotropical Dry Forest #####
#########################################################################################################

## Code for the main analysis in Weiss et al. (2026) - Trait-filtered responses of mammal communities to land-use change in a neotropical dry forest

rm(list = ls()); gc()

library(patchwork)
library(ggh4x)
library(spOccupancy)
library(broom)
library(scales)
library(scico)
library(brms)
library(hrbrthemes)
library(ggthemes)
library(geometry)
library(ggridges)
library(tidybayes)
library(tidyverse)

select <- dplyr::select
recode <- dplyr::recode

export = FALSE # Export plots?
model_all = FALSE # fit all models?
clean = FALSE # clean output in between? 

# read input data
captures <- read_csv("Data/species_data.csv")
camtraps <- read_csv("Data/camtraps_clean.csv")
camop <- as.matrix(read_csv("Data/camop_problem.csv"))
covariates <- read_csv("Data/forest_covariates.csv")

# ========= Data Preparation ==========

# clean species data for analysis
species <- captures %>%
  filter(Station %in% camtraps$Station) %>%
  filter(Category %in% c("artiodactyla", "carnivora", "marsupialia", "perissodactyla", "primates", "rodentia", "xenarthra")) %>%
  filter(DateTimeOriginal >= "2017-01-10 17:01:26") %>%
  filter(DateTimeOriginal <= "2023-10-23 23:59:59") %>%
  mutate(
    Station = as.character(Station),
    accepted_bin = paste(Genus, Species, sep = "_"),
    accepted_bin = str_to_lower(str_replace_all(accepted_bin, " ", "_")),
    accepted_bin = str_replace_all(accepted_bin, "\\.", "_"),
    accepted_bin = str_replace_all(accepted_bin, "\\-", "_"),
    accepted_bin = str_replace_all(accepted_bin, "\\(", ""),
    accepted_bin = str_replace_all(accepted_bin, "\\)", ""),
    accepted_bin = str_replace_all(accepted_bin, "\\/", "_")
  ) %>%
  # Make sure all data is identified to the species level
  filter(complete.cases(Species)) %>%
  # Assess temporal independence - Arrange by station, name and time
  arrange(Station, accepted_bin, DateTimeOriginal) %>%
  # Calculate the time difference between consecutive captures of the same species at the same station
  group_by(Station, accepted_bin) %>%
  mutate(delta_time = difftime(DateTimeOriginal, lag(DateTimeOriginal), units = "mins")) %>%
  # Keep only the captures that are at least 30 minutes apart from the last capture of the same species
  filter(is.na(delta_time) | delta_time >= 30) %>%
  ungroup() %>%
  # Remove the temporary delta_time column
  select(-delta_time)

# Integrate pantheria traits, see which traits predict which response
pantheria <- read.table(
  file = "Data/PanTHERIA.txt",
  header = TRUE, sep = "\t", na.strings = c("-999", "-999.00")
) %>%
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
    trophic_level = X6.2_TrophicLevel
  )

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
  select(
    accepted_bin, order, adult_body_mass_g, home_range_km2,
    activity_cycle, diet_breadth, habitat_breadth,
    terrestriality, trophic_level
  )

## Home range size based on added literature records, originally not represented in PanTheria version used...
species_traits[5,4] <- 0.02 # agouti
species_traits[8,4] <- 0.03 # paca
species_traits[13,4] <- 1.17 # deer
species_traits[18,4] <- 1.11 # racoon
species_traits[22,4] <- 8.31 # tapir

# ========= Occupancy Analysis =========

create_detection_array <- function(species_data, species_min = 10, K = 12) {
  library(dplyr)
  library(tidyr)
  library(lubridate)
  
  species_data <- species_data %>%
    mutate(
      Date = as.Date(DateTimeOriginal),
      year = year(Date),
      site_year = paste(Station, year, sep = "_")
    )
  
  # Count total detections per species
  species_counts <- species_data %>%
    count(accepted_bin) %>%
    dplyr::filter(n >= 10)
  
  species_list <- species_counts$accepted_bin
  
  # Filter to valid species
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
  filter(site_year %in% dimnames(detection_array)[[2]]) %>%
  group_by(year) %>%
  mutate(
    treecover_z = as.numeric(scale(treecover)),
    edge_density_z = as.numeric(scale(edge_density_forest)),
    patch_density_z = as.numeric(scale(patch_density_forest)),
    shape_index_z = as.numeric(scale(shape_index_forest))
  ) %>%
  ungroup() %>%
  arrange(match(site_year, dimnames(detection_array)[[2]])) %>%
  as.data.frame()

rownames(site_covs) <- dimnames(detection_array)[[2]]

# make sure things match - should return TRUE
nrow(site_covs) == dim(detection_array)[2]
all.equal(rownames(site_covs), dimnames(detection_array)[[2]])

# Clean up detection array
detection_array[is.na(detection_array)] <- 0
any(is.na(detection_array)) # Should return FALSE

# name months dimension
dimnames(detection_array)[[3]] <- paste0("Month", 1:dim(detection_array)[3])

# assess correlations
r_edge  <- cor(site_covs$treecover_z, site_covs$edge_density_z,  use = "complete.obs")
r_patch <- cor(site_covs$treecover_z, site_covs$patch_density_z, use = "complete.obs")
r_shape <- cor(site_covs$treecover_z, site_covs$shape_index_z,   use = "complete.obs")

# model meta
model_meta <- tibble::tribble(
  ~model,            ~covariates,                         ~uses_edge, ~uses_patch, ~uses_shape, ~interaction,
  "Edge",            "treecover + edge_density",          TRUE,       FALSE,       FALSE,       FALSE,
  "Patch",           "treecover + patch_density",         FALSE,      TRUE,        FALSE,       FALSE,
  "Shape",           "treecover + shape_index",           FALSE,      FALSE,       TRUE,        FALSE,
  "Edge + Shape",    "treecover + edge_density + shape",  TRUE,       FALSE,       TRUE,        FALSE,
  "Edge Interaction","treecover * edge_density",          TRUE,       FALSE,       FALSE,       TRUE
)

# Define model framework 
fit_model <- function(formula) {
  msPGOcc(
    occ.formula = formula,
    det.formula = ~1,
    data = list(y = detection_array, occ.covs = site_covs),
    n.samples = 15000,
    n.burn = 2000,
    n.thin = 10,
    n.chains = 3,
    verbose = TRUE
  )
}

if (model_all) {
  
  # Find best fragmentation metric with additive models, then see if we need an interaction
  models <- list(
    "Edge" = fit_model(~ treecover_z + edge_density_z),
    "Patch" = fit_model(~ treecover_z + patch_density_z),
    "Shape" = fit_model(~ treecover_z + shape_index_z),
    "Edge + Shape" = fit_model(~ treecover_z + edge_density_z + shape_index_z),
    "Edge Interaction" = fit_model(~ treecover_z * edge_density_z)
  )
  
  waics <- purrr::map(models, spOccupancy::waicOcc)
  waic_tbl <- tibble(
    model = names(waics),
    WAIC = sapply(waics, `[[`, "WAIC")
  ) %>% arrange(WAIC)
  
  waic_tbl %>%
    mutate(delta_WAIC = WAIC - min(WAIC)) %>%
    ggplot(aes(reorder(model, delta_WAIC), delta_WAIC)) +
    geom_col() +
    coord_flip() +
    ylab("ΔWAIC") + xlab("Model") +
    theme_minimal()
  
  best_model <- waic_tbl$model[1]
  ms_model <- models[[best_model]]
  cat("Best model selected:", best_model, "with WAIC =", waic_tbl$WAIC[1], "\n")
  
  model_selection_tbl <- waic_tbl %>%
    dplyr::left_join(model_meta, by = "model") %>%
    dplyr::mutate(
      delta_WAIC = WAIC - min(WAIC, na.rm = TRUE),
      weight_raw = exp(-0.5 * delta_WAIC),
      weight     = weight_raw / sum(weight_raw, na.rm = TRUE),
      r_treecover_edge  = dplyr::if_else(uses_edge,  r_edge,  NA_real_),
      r_treecover_patch = dplyr::if_else(uses_patch, r_patch, NA_real_),
      r_treecover_shape = dplyr::if_else(uses_shape, r_shape, NA_real_)
    ) %>%
    dplyr::select(model, covariates, WAIC, delta_WAIC, weight,
                  r_treecover_edge, r_treecover_patch, r_treecover_shape, interaction) %>%
    dplyr::arrange(WAIC) %>%
    dplyr::mutate(
      dplyr::across(c(WAIC, delta_WAIC, weight,
                      r_treecover_edge, r_treecover_patch, r_treecover_shape),
                    ~ round(.x, 3))
    )
  
  print(model_selection_tbl)
  
} else {
  # just fit the interaction model
  ms_model <- fit_model(~ treecover_z * edge_density_z)
}


# ----- Check Model output -----

summary(ms_model, level = "community")
summary(ms_model, level = "species")
psi <- fitted(ms_model)
summary(psi)

# Extract posterior summary stats for coefficients
beta_samples <- ms_model$beta.samples
beta_summary <- as.data.frame(summary(beta_samples)$statistics)
beta_ci <- as.data.frame(summary(beta_samples)$quantiles)

# ----- Assess Model Diagnostics -----

# Extract Rhat and ESS as named vectors
beta_stats <- summary(beta_samples)$statistics
names(ms_model$rhat$beta) <- rownames(beta_stats)
names(ms_model$ESS$beta)  <- rownames(beta_stats)

rhat_vec <- ms_model$rhat$beta
ess_vec  <- ms_model$ESS$beta

# Build a tibble from the beta samples summary
beta_df <- as.data.frame(summary(beta_samples)$statistics) %>%
  rownames_to_column("param") %>%
  filter(str_detect(param, "^(treecover_z|edge_density|treecover_z:edge_density_z)-")) %>%
  mutate(
    covariate   = str_extract(param, "^[^\\-]+"),
    species_raw = str_extract(param, "(?<=-).*"),
    species     = tools::toTitleCase(gsub("_", " ", species_raw))
  )

param_names <- beta_df$param

# Subset by name
rhat_filtered <- rhat_vec[param_names]
ess_filtered  <- ess_vec[param_names]

# Sanity check
stopifnot(all(names(rhat_filtered) == param_names))

# Build a tibble from the beta samples summary
beta_df <- as.data.frame(summary(beta_samples)$statistics) %>%
  rownames_to_column("param") %>%
  filter(str_detect(param, "^(treecover_z|edge_density_z|treecover_z:edge_density_z)-")) %>%
  mutate(
    covariate   = str_extract(param, "^[^\\-]+"),
    species_raw = str_extract(param, "(?<=-).*"),
    species     = tools::toTitleCase(gsub("_", " ", species_raw))
  )

# Turn diagnostics into a df with param names
param_names <- beta_df$param
rhat_filtered <- rhat_vec[param_names]
ess_filtered <- ess_vec[param_names]

length(rhat_filtered) == nrow(beta_df) # should  be TRUE

diagnostics_table <- beta_df %>%
  mutate(
    Rhat = rhat_filtered,
    ESS = ess_filtered,
    converged = Rhat < 1.1,
    ESS_gt_1000 = ESS > 1000
  ) %>%
  select(species, covariate, Mean, SD, Rhat, ESS, converged, ESS_gt_1000)

all(names(rhat_filtered) == beta_df$param) # should be TRUE

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
    Covariate = dplyr::recode(Covariate,
                              "treecover_z"   = "Forest cover",
                              "edge_density_z" = "Edge Density",
                              "treecover_z:edge_density_z"   = "Interaction"
    ),
    Mean = round(Mean, 3),
    SD = round(SD, 3),
    Rhat = round(Rhat, 3),
    ESS = round(ESS, 0)
  ) %>%
  select(-Converged, -ESS_gt_1000)

if (export) {
  write_csv(diagnostics_clean, "Output/Tables/model_diagnostics.csv")
}

# ----- Community level effects 

# Update the community-level effect plot as well
beta_comm <- ms_model$beta.comm.samples %>%
  as_tibble() %>%
  pivot_longer(names_to = "term", values_to = "value", cols = everything()) %>%
  filter(term != "(Intercept)") %>%
  mutate(term = dplyr::recode(term,
                              "edge_density_z" = "Edge Density",
                              "treecover_z" = "Forest Cover",
                              "treecover_z:edge_density_z" = "Interaction"
  ),
  term = factor(term, levels = c("Forest Cover", "Edge Density", "Interaction")))

p_comm <- ggplot(beta_comm, aes(x = value, y = fct_rev(term), fill = term)) +
  stat_halfeye(.width = c(0.5, 0.8, 0.95), 
               point_interval = median_qi, 
               slab_alpha = 0.5,  # softer fill
               fill_type = "gradient") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  scale_fill_manual(values = c("darkgreen", "blue4", "cyan4")) +
  labs(x = "Community-level effect size (posterior median ± 95% credible interval)", y = NULL) +
  theme_few(base_size = 8) +
  theme(
    strip.text = element_text(face = "italic", size = 8, hjust = 0),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_text(hjust = 1, size = 8),  
    legend.position = "none",
    panel.background = element_blank(),
    plot.background = element_blank(),
    strip.background = element_blank()
  )


# Check among species variance 
beta_draws_df <- as.data.frame(ms_model$beta.samples) %>%
  rowid_to_column("iteration") %>%
  pivot_longer(-iteration, names_to = "term", values_to = "value") %>%
  separate(term, into = c("covariate", "species"), sep = "-", extra = "merge")

beta_var_df <- beta_draws_df %>%
  filter(covariate != "(Intercept)") %>%
  group_by(iteration, covariate) %>%
  summarise(among_species_var = var(value), .groups = "drop")

beta_var_df %>%
  group_by(covariate) %>%
  summarise(
    mean = mean(among_species_var),
    `2.5%` = quantile(among_species_var, 0.025),
    `97.5%` = quantile(among_species_var, 0.975)
  )

beta_var_df %>%
  mutate(covariate = recode(covariate,
                            "treecover_z" = "Forest Cover",
                            "edge_density_z" = "Edge Density",
                            "treecover_z:edge_density_z" = "Interaction"),
         covariate = factor(covariate, levels = c("Forest Cover", "Edge Density", "Interaction"))) %>%
  ggplot(aes(x = among_species_var, y = fct_rev(covariate), fill = covariate)) +
  stat_halfeye(.width = c(0.66, 0.95), slab_alpha = 0.6) +
  scale_fill_manual(values = c("darkgreen", "blue4", "cyan4")) +
  labs(
    title = "Among-species variance in responses",
    x = expression("Posterior variance of " * beta), y = NULL, fill = NULL
  ) +
  theme_linedraw() +
  theme(
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(hjust = 1),
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0),
    plot.margin = margin(t = 5, r = 10, b = 0, l = 10),
    panel.background = element_blank(),
    plot.background = element_blank(),
    strip.background = element_blank(),
    panel.grid = element_blank())

p_species <- as_tibble(ms_model$beta.samples) %>%
  pivot_longer(cols = everything(),
               names_to = "param",
               values_to = "value") %>%
  separate(param, into = c("parameter", "species"), sep = "-", extra = "merge") %>%
  # Filter to focus species
  filter(species %in% c("cerdocyon_thous",
                        "dasyprocta_punctata",
                        "eira_barbara",
                        "leopardus_pardalis",
                        "panthera_onca",
                        "pecari_tajacu",
                        "puma_concolor",
                        "tapirus_terrestris")) %>%
  filter(parameter != "(Intercept)") %>%
  mutate(
    parameter = recode(parameter,
                       treecover_z = "Forest\nCover",
                       edge_density_z = "Edge\nDensity",
                       `treecover_z:edge_density_z` = "Interaction")
  ) %>%
  mutate(
    species = str_replace_all(species, "_", " "),
    species = str_to_lower(species),
    species = str_replace(species, "^(\\w+)", ~ str_to_title(.x))
  ) %>%
  mutate(parameter = factor(parameter, levels = c("Forest\nCover", "Edge\nDensity", "Interaction"))) %>%
  mutate(species = recode(species,
                          "Cerdocyon thous" = "Cerdocyon\nthous",
                          "Dasyprocta punctata" = "Dasyprocta\npunctata",
                          "Eira barbara" = "Eira barbara",
                          "Leopardus pardalis" = "Leopardus\npardalis",
                          "Panthera onca" = "Panthera onca",
                          "Pecari tajacu" = "Pecari tajacu",
                          "Puma concolor" = "Puma concolor",
                          "Tapirus terrestris" = "Tapirus\nterrestris"
  )) %>%
  
  ggplot(aes(x = value, y = fct_rev(parameter), fill = parameter)) +
  stat_halfeye(.width = c(0.66, 0.95), slab_alpha = 0.6, linewidth = 0.5, point_size = .6) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  facet_wrap(~ species, nrow = 2, ncol = 4) + 
  scale_fill_manual(values = c("darkgreen", "blue4", "cyan4")) +
  labs(
    x = "Species-level effect size (posterior median ± 95% credible interval)", y = NULL, fill = NULL) +
  theme_few(base_size = 8) +
  theme(
    strip.text = element_text(face = "italic", size = 8, hjust = 0),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_text(hjust = 1, size = 8),  
    legend.position = "none",
    panel.background = element_blank(),
    plot.background = element_blank(),
    strip.background = element_blank()
  )

# ----- Trait Filtering ----- 

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
  filter(covariate %in% c("treecover_z", "edge_density_z", "treecover_z:edge_density_z")) %>%
  left_join(species_traits, by = c("species" = "accepted_bin")) %>%
  dplyr::mutate(
    log_mass = scale(log(adult_body_mass_g)),
    log_range = scale(log(home_range_km2)),
    diet_breadth = scale(as.numeric(diet_breadth)),
    activity_cycle_coded = dplyr::recode(activity_cycle, `1` = -1, `3` = 0, `2` = 1),
    habitat_breadth_coded = scale(habitat_breadth)
  ) %>%
  dplyr::select(
    covariate, species, Mean, `2.5%`, `97.5%`,
    log_mass, log_range, diet_breadth,
    activity_cycle_coded, habitat_breadth_coded
  )

# Add standard deviations to the coeff dtata
beta_sds <- beta_summary %>%
  as.data.frame() %>%
  rownames_to_column("term") %>%
  separate(term, into = c("covariate", "species"), sep = "-", extra = "merge") %>%
  select(covariate, species, SD = SD)

coef_traits <- coef_traits %>%
  left_join(beta_sds, by = c("covariate", "species"))

# Create a function for running brms per trait
run_bayes_trait_model <- function(df, covariate_name, trait_name) {
  subdat <- df %>%
    filter(covariate == covariate_name, !is.na(.data[[trait_name]]), !is.na(SD))
  
  formula <- as.formula(paste0("Mean | se(SD, sigma = TRUE) ~ ", trait_name))
  
  brm(
    formula = formula,
    data = subdat,
    family = gaussian(),
    prior = c(
      prior(normal(0, 1), class = "b"),
      prior(normal(0, 1), class = "Intercept")
    ),
    chains = 4, iter = 4000, cores = 4,
    silent = TRUE, refresh = 0
  )
}

trait_vars <- c("log_mass", "log_range", "diet_breadth", "activity_cycle_coded", "habitat_breadth_coded")
covariates <- unique(coef_traits$covariate)

# loop Bayesian models through traits
bayes_trait_models <- crossing(cov = covariates, trait = trait_vars) %>%
  mutate(model = map2(cov, trait, ~ run_bayes_trait_model(coef_traits, .x, .y)))

# Extract posterior draws from the models
bayes_results <- pmap_dfr(
  list(bayes_trait_models$model, bayes_trait_models$trait, bayes_trait_models$cov),
  function(mod, tr, cv) {
    parname <- paste0("b_", tr)
    mod %>%
      tidybayes::gather_draws(!!sym(parname), b_Intercept) %>%
      mutate(trait = tr, cov = cv)
  }
) %>%
  mutate(
    trait_label = dplyr::recode(trait,
                                log_mass = "Body mass",
                                log_range = "Home range size",
                                diet_breadth = "Diet breadth",
                                activity_cycle_coded = "Diurnality",
                                habitat_breadth_coded = "Habitat breadth"),
    covariate = dplyr::recode(cov,
                              treecover_z = "Forest Cover",
                              edge_density_z = "Edge Density",
                              `treecover_z:edge_density_z` = "Interaction")
  )

# Get Bayesian r squared values 
bayes_rsq <- bayes_trait_models %>%
  mutate(r2 = map_dbl(model, ~ brms::bayes_R2(.x)[1])) %>%
  select(cov, trait, r2) %>%
  mutate(r2_label = paste0("Bayes R² = ", round(r2, 2))) %>%
  transmute(
    trait = trait,
    covariate = case_when(
      cov == "treecover_z" ~ "Forest Cover",
      cov == "edge_density_z" ~ "Edge Density",
      cov == "treecover_z:edge_density_z" ~ "Interaction"
    ),
    r2_label = paste0("R² = ", round(r2, 2))
  ) %>%
  mutate(
    covariate_label = factor(covariate, levels = c("Forest Cover", "Edge Density", "Interaction")),
    trait_label = dplyr::recode(trait,
                                log_mass = "Body\nmass",
                                log_range = "Homerange\nsize",
                                diet_breadth = "Diet\nbreadth",
                                activity_cycle_coded = "Diurnality",
                                habitat_breadth_coded = "Habitat\nbreadth"
    ),
    trait_label_r2 = paste0(trait_label, "\n", r2_label)
  )

trait_draws <- bayes_results %>%
  filter(.variable != "b_Intercept") %>%
  mutate(
    covariate = fct_relevel(covariate, "Forest Cover", "Edge Density", "Interaction"),
    trait_label = dplyr::recode(trait,
                                log_mass = "Body\nmass",
                                log_range = "Homerange\nsize",
                                diet_breadth = "Diet\nbreadth",
                                activity_cycle_coded = "Diurnality",
                                habitat_breadth_coded = "Habitat\nbreadth"),
    covariate_label = recode(covariate,
                             "Forest Cover" = "Forest Cover",
                             "Edge Density" = "Edge Density",
                             "Interaction" = "Interaction")
  ) %>%
  mutate(
    covariate_label = factor(covariate_label, levels = c("Forest Cover", "Edge Density", "Interaction")),
    trait_label = factor(trait_label, levels = c("Body\nmass", "Homerange\nsize", "Diet\nbreadth", "Diurnality", "Habitat\nbreadth"))
  ) %>%
  left_join(bayes_rsq %>% 
              select(trait_label, covariate_label, trait_label_r2),
            by = c("trait_label", "covariate_label")) 

r2_annot <- bayes_rsq %>%
  left_join(
    trait_draws %>%
      group_by(trait_label, covariate_label) %>%
      summarise(median_x = median(.value), .groups = "drop"),
    by = c("trait_label", "covariate_label"))

# Trait filtering
p_traits <- ggplot(trait_draws, aes(x = .value, y = fct_rev(trait_label), fill = trait_label)) +
  stat_halfeye(.width = c(0.66, 0.95), slab_alpha = 0.6, fill_type = "gradient") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  facet_wrap(~ covariate_label, nrow = 1, strip.position = "top") +
  scale_fill_viridis_d(option = "viridis", end = 0.95) +
  labs(
    x = "Trait–Environment slope (posterior median ± 95% credible interval)",
    y = NULL
  ) +
  # R² annotation
  geom_label(
    data = r2_annot,
    aes(x = median_x, y = fct_rev(trait_label), label = r2_label),
    inherit.aes = FALSE,
    vjust = 1.5,           
    size = 2.5,
    label.size = NA,       
    fill = "white",       
    alpha = .9,           
    label.padding = unit(0.15, "lines"), 
    fontface = "italic"
  ) +
  theme_few(base_size = 8) +
  theme(
    strip.text = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_text(hjust = 1, size = 8),  
    legend.position = "none",
    plot.margin = margin(t = 5, r = 20, b = 5, l = 10),
    panel.background = element_blank(),
    plot.background = element_blank(),
    strip.background = element_blank()
  )

# Build fig2 as compound figure with a. labeling
fig2 <- (((p_comm / p_species) + plot_layout(heights = c(1, 1.6))) | p_traits) +
  plot_layout(widths = c(1, 1.2)) +
  plot_annotation(tag_levels = 'a', tag_prefix = "", tag_suffix = ".") &
  theme(plot.tag = element_text(size = 13))

if (export) {
  
  fig2 # check plot
  
  # Figure 2
  ggsave("Output/Figures/fig2.png",
         fig2,
         width = 250,
         height = 170,
         units = "mm",
         dpi = 600)
}

if (clean) {
  
  #Clean environment
  rm(bayes_rsq, beta_ci, beta_comm, beta_df, beta_draws_df, beta_sds, beta_summary, beta_stats,
     beta_var_df, diagnostics_clean, diagnostics_table, psi, r2_annot, trait_draws,
     beta_samples, ess_filtered, ess_vec, param_names, rhat_filtered, rhat_vec,
     p_comm, p_species, p_traits)
}

## ========== PDP Analysis ==========

# Standardize function
scale_cov <- function(x, mean_x, sd_x) (x - mean_x) / sd_x

# Means and SDs
mean_treecover <- mean(site_covs$treecover, na.rm = TRUE)
sd_treecover <- sd(site_covs$treecover, na.rm = TRUE)
mean_edgedens <- mean(site_covs$edge_density_forest, na.rm = TRUE)
sd_edgedens <- sd(site_covs$edge_density_forest, na.rm = TRUE)

# Sequences in raw units
tree_seq <- seq(min(site_covs$treecover, na.rm = TRUE), max(site_covs$treecover, na.rm = TRUE), length.out = 50)
edge_seq <- seq(min(site_covs$edge_density_forest, na.rm = TRUE), max(site_covs$edge_density_forest, na.rm = TRUE), length.out = 50)

# ----- Community surface -----

pdp_grid <- expand.grid(treecover = tree_seq, edge_density_forest = edge_seq) %>%
  mutate(
    treecover_z = scale_cov(treecover, mean_treecover, sd_treecover),
    edge_density_z = scale_cov(edge_density_forest, mean_edgedens, sd_edgedens)
  )

X.0 <- model.matrix(~ treecover_z * edge_density_z, data = pdp_grid)
pred <- predict(ms_model, X.0 = X.0, type = "occupancy", ignore.RE = FALSE)

comm_occ <- map_dfr(seq_len(nrow(pdp_grid)), function(i) {
  psi_samples <- pred$psi.0.samples[, , i]
  row_mean <- rowMeans(psi_samples)
  tibble(
    treecover = pdp_grid$treecover[i],
    edge_density_forest = pdp_grid$edge_density_forest[i],
    mean = mean(row_mean)
  )
})

# Mask extrapolated values
hull_pts <- site_covs %>%
  select(treecover, edge_density_forest) %>%
  drop_na() %>%
  as.matrix() %>%
  convhulln(return.non.triangulated.facets = TRUE)

inside_hull <- inhulln(hull_pts, as.matrix(comm_occ[, c("treecover", "edge_density_forest")]))
comm_occ$inside_hull <- inside_hull

p_comm <- ggplot(comm_occ %>% filter(inside_hull), 
                 aes(x = treecover, y = edge_density_forest, fill = mean)) +
  geom_raster(alpha = 0.85) +
  scale_fill_scico(name = "Est. Community-level\nOccupancy", palette = "cork", midpoint = median(comm_occ$mean)) +
  labs(x = "Forest cover (%)", y = "Edge density (m/ha)") +
  theme_few(base_size = 6) +
  theme(
    strip.text = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    axis.title.x = element_text(size = 6),  
    legend.position = "bottom",
    plot.margin = margin(t = 5, r = 20, b = 5, l = 10),
    panel.background = element_blank(),
    plot.background = element_blank(),
    strip.background = element_blank()
  )

# ----- Home Range Size Filtering -----

trait_panel <- function(species_group, label, var_name) {
  pdp_data <- tibble(
    treecover = tree_seq,
    edge_density_forest = mean(site_covs$edge_density_forest, na.rm = TRUE)
  ) %>%
    mutate(
      treecover_z = scale_cov(treecover, mean_treecover, sd_treecover),
      edge_density_z = scale_cov(edge_density_forest, mean_edgedens, sd_edgedens)
    )
  
  X.0 <- model.matrix(~ treecover_z * edge_density_z, data = pdp_data)
  pred <- predict(ms_model, X.0 = X.0, type = "occupancy", ignore.RE = FALSE)
  
  species_names <- gsub("treecover_z-", "", grep("^treecover_z-", colnames(ms_model$beta.samples), value = TRUE))
  species_ids <- setNames(seq_along(species_names), species_names)
  sel <- species_ids[names(species_ids) %in% species_group]
  
  map_dfr(seq_len(nrow(pdp_data)), function(i) {
    psi_samples <- pred$psi.0.samples[, sel, i]
    row_mean <- rowMeans(psi_samples)
    tibble(
      treecover = pdp_data$treecover[i],
      group = label,
      mean = mean(row_mean),
      lower = quantile(row_mean, 0.025),
      upper = quantile(row_mean, 0.975)
    )
  })
}

# Panel B: Home range
hr_quant <- quantile(species_traits$home_range_km2, c(0.25, 0.75), na.rm = TRUE)
hr_occ <- bind_rows(
  trait_panel(species_traits$accepted_bin[species_traits$home_range_km2 <= hr_quant[1]], "Small home range"),
  trait_panel(species_traits$accepted_bin[species_traits$home_range_km2 >= hr_quant[2]], "Large home range")
)

p_hr <- ggplot(hr_occ, aes(x = treecover, y = mean, color = group, fill = group)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA) +
  geom_line(size = 1) +
  scale_color_manual(values = c("darkgreen", "goldenrod"), name = "Trait-Filtering") +
  scale_fill_manual(values = c("darkgreen", "goldenrod"), name = "Trait-Filtering") +
  labs(x = "Forest cover (%)", y = "Predicted occupancy") +
  theme_few(base_size = 6) +
  theme(
    strip.text = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    axis.title.x = element_text(size = 6),  
    legend.position = "bottom",
    plot.margin = margin(t = 5, r = 20, b = 5, l = 10),
    panel.background = element_blank(),
    plot.background = element_blank(),
    strip.background = element_blank()
  )

fig3 <- (p_comm | p_hr) +
  plot_annotation(tag_levels = 'a', tag_prefix = "", tag_suffix = ".") &
  theme(plot.tag = element_text(size = 10))

if (export) {
  
  fig3 # check plot
  
  # Figure 2
  ggsave("Output/Figures/fig3.tiff",
         fig3,
         width = 180,
         height = 100,
         units = "mm",
         dpi = 600)
}
