#########################################################################################################
##### Trait-Filtered Responses of Mammal Communities to Land-Use Change in a Neotropical Dry Forest #####
#########################################################################################################

rm(list = ls()); gc()

library(patchwork)
library(ggh4x)
library(spOccupancy)
library(broom)
library(scales)
library(scico)
library(brms)
library(geometry)
library(ggridges)
library(tidybayes)
library(tidyverse)

select <- dplyr::select
recode <- dplyr::recode

export = FALSE # Export plots?
model_all = FALSE # fit all models?

# read input data
captures <- read_csv("/Users/merlin/Documents/Senckenberg/WildLive/Mammals/Code/Data/species_data.csv")
camtraps <- read_csv("/Users/merlin/Documents/Senckenberg/WildLive/Mammals/Code/Data/camtraps_clean.csv")
camop <- as.matrix(read_csv("/Users/merlin/Documents/Senckenberg/WildLive/Mammals/Code/Data/camop_problem.csv"))
covariates <- read_csv("/Users/merlin/Documents/Senckenberg/WildLive/Mammals/Code/Data/forest_covariates.csv")

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
  select(-Batch) %>%
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
  file = "/Users/merlin/Documents/Senckenberg/WildLive/Mammals/Code/Data/PanTHERIA.txt",
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
cor(site_covs$patch_density_z, site_covs$edge_density_z)
cor(site_covs$patch_density_z, site_covs$shape_index_z)
cor(site_covs$edge_density_z, site_covs$shape_index_z)

cor(site_covs$treecover_z, site_covs$edge_density_z)
cor(site_covs$treecover_z, site_covs$patch_density_z)
cor(site_covs$treecover_z, site_covs$shape_index_z)

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
  
} else {
  
  # just fit the interaction model
  ms_model <- fit_model(~ treecover_z * edge_density_z)
}

# ----- Check Model output -----

summary(ms_model, level = "community")
summary(ms_model, level = "species")
summary(ms_model, level = "community")$beta.comm
summary(ms_model, level = "species")$beta
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

param_names <- beta_df$param  # from your previous code

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
  write_csv(diagnostics_clean, "/Users/merlin/Documents/Senckenberg/WildLive/Mammals/Code/Output/Tables/model_diagnostics.csv")
}

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
                       treecover_z = "Forestcover",
                       edge_density_z = "Edge Density",
                       `treecover_z:edge_density_z` = "Interaction")
  )


# --- Build plot ---

# Update the community-level effect plot as well
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
  mutate(term = dplyr::recode(term,
                              "edge_density_z" = "Edge Density",
                              "treecover_z" = "Forest Cover",
                              "treecover_z:edge_density_z" = "Interaction"
  )) %>%
  mutate(sig_label = ifelse(`2.5%` > 0 | `97.5%` < 0, "*", ""))

beta_comm$term <- factor(beta_comm$term,
                         levels = c("Forest Cover", "Edge Density", "Interaction")
)

# Summarize community effects (no trait filtering)
comm_plot_data <- beta_comm %>%
  rownames_to_column("covariate") %>%
  filter(covariate != "(Intercept)") %>%
  mutate(term = dplyr::recode(term,
                            "Forest Cover" = "Treecover",
                            "Edge Density" = "Edge Density",
                            "Interaction" = "Interaction"))

theme_set(theme_bw(base_family = "sans"))

comm_plot <- ggplot(comm_plot_data, aes(x = term, y = Mean, ymin = `2.5%`, ymax = `97.5%`)) +
  geom_pointrange(aes(color = covariate), size = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "A. Community-level effects",
       x = NULL, y = "Effect size\n(logit scale)") +
  scale_color_viridis_d(option = "magma", end = 0.8) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 8),
      axis.title = element_text(size = 10),
      plot.title = element_text(size = 12),
      legend.position = "none")

# ----- Trait filtering plots -----

# Only keep the trait slopes, not the intercepts
trait_slopes <- bayes_results %>%
  filter(.variable != "b_Intercept") %>%
  # Summarize posterior draws per trait-covariate
  group_by(covariate, trait_label) %>%
  mean_qi(.value, .width = 0.95) %>%
  ungroup()

traits_interaction <- trait_slopes %>%
  filter(covariate == "Interaction") %>%
  ggplot(aes(x = trait_label, y = .value, ymin = .lower, ymax = .upper, color = trait_label)) +
  geom_pointrange(position = position_dodge(width = 0.5), size = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_viridis_d(option = "viridis", end = 0.95) +
  labs(title = "B3. Interaction",
       x = NULL, y = NULL) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 10),
        plot.title = element_text(size = 10, hjust = 0.5),
        legend.position = "none")

# Tree cover trait plot
traits_treecover <- trait_slopes %>%
  filter(covariate == "Forestcover") %>%
  ggplot(aes(x = trait_label, y = .value, ymin = .lower, ymax = .upper, color = trait_label)) +
  geom_pointrange(position = position_dodge(width = 0.5), size = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_viridis_d(option = "viridis", end = 0.95) +
  labs(title = "B1. Tree Cover",
       x = NULL, y = "Trait–Environment Slope\n(posterior mean ± 95% CI)") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 10),
        plot.title = element_text(size = 10, hjust = 0.5),
        legend.position = "none")

# Edge density trait plot
traits_edge <- trait_slopes %>%
  filter(covariate == "Edge Density") %>%
  ggplot(aes(x = trait_label, y = .value, ymin = .lower, ymax = .upper, color = trait_label)) +
  geom_pointrange(position = position_dodge(width = 0.5), size = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_viridis_d(option = "viridis", end = 0.95) +
  labs(title = "B2. Edge Density",
       x = NULL, y = NULL) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 10),
        plot.title = element_text(size = 10, hjust = 0.5),
        legend.position = "none")

# Adjust all trait plots to same y scale
trait_ylim <- range(trait_slopes$.lower, trait_slopes$.upper)
traits_treecover <- traits_treecover + coord_cartesian(ylim = trait_ylim)
traits_edge <- traits_edge + coord_cartesian(ylim = trait_ylim)
traits_interaction <- traits_interaction + coord_cartesian(ylim = trait_ylim)
trait_ylim <- range(trait_slopes$.lower, trait_slopes$.upper)
traits_treecover <- traits_treecover + coord_cartesian(ylim = trait_ylim)
traits_edge <- traits_edge + coord_cartesian(ylim = trait_ylim)
traits_interaction <- traits_interaction + coord_cartesian(ylim = trait_ylim)

# Assemble final figure with heading above trait plots
trait_heading <-  wrap_elements(
  grid::textGrob("B. Trait filtering", 
                 gp = grid::gpar(fontsize = 12), 
                 just = "left", 
                 x = unit(0, "npc")))


final_plot <- comm_plot / 
  trait_heading / 
  (traits_treecover | traits_edge | traits_interaction) +
  plot_layout(heights = c(1, 0.2, 1))

if (export) {
  
  # Figure 2
  ggsave("/Users/merlin/Documents/Senckenberg/WildLive/Mammals/Code/Output/Figures/fig2.png",
         final_plot,
         width = 150,
         height = 130,
         units = "mm",
         dpi = 600)
}

# Check probabilities of positive effects
bayes_results %>%
  filter(.variable != "b_Intercept") %>%
  group_by(covariate, trait_label) %>%
  summarise(prob_positive = mean(.value > 0)) %>%
  arrange(desc(prob_positive))

# Posteriour Distributions of trait slopes
bayes_results %>%
  filter(.variable != "b_Intercept") %>%
  mutate(
    trait_label = factor(trait_label, levels = c("Body mass", "Home range size", "Diet breadth", "Diurnality", "Habitat breadth"))
  ) %>%
  ggplot(aes(x = .value, y = trait_label, fill = covariate)) +
  geom_density_ridges(alpha = 0.7, scale = 1.2, rel_min_height = 0.01) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_viridis_d(option = "magma", end = 0.8) +
  labs(
    x = "Posterior slope estimate",
    y = NULL,
    title = "Posterior distributions of trait slopes by covariate"
  ) +
  facet_wrap(~covariate) +
  theme_bw() +
  theme(legend.position = "none")

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

# ----- Panel A1: Community surface -----
pdp_grid <- expand.grid(treecover = tree_seq, edge_density_forest = edge_seq) %>%
  mutate(
    treecover_z = scale_cov(treecover, mean_treecover, sd_treecover),
    edge_density_z = scale_cov(edge_density_forest, mean_edgedens, sd_edgedens)
  )

X.0 <- model.matrix(~ treecover_z * edge_density_z, data = pdp_grid)
pred <- predict(ms_model, X.0 = X.0, type = "occupancy", ignore.RE = FALSE)

panel_a1_data <- map_dfr(seq_len(nrow(pdp_grid)), function(i) {
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

inside_hull <- inhulln(hull_pts, as.matrix(panel_a1_data[, c("treecover", "edge_density_forest")]))
panel_a1_data$inside_hull <- inside_hull

panel_a1 <- ggplot(panel_a1_data %>% filter(inside_hull), 
                   aes(x = treecover, y = edge_density_forest, fill = mean)) +
  geom_raster(alpha = 0.85) +
  scale_fill_scico(name = NULL, palette = "cork", midpoint = median(panel_a1_data$mean)) +
  labs(title = "A1. Community occupancy", x = "Forest cover (%)", y = "Edge density (m/ha)") +
  theme_classic(base_size = 13) +
  theme(legend.position = c(0.1, 0.85))

# ----- Panel A2: Interaction PDP -----
edge_q <- quantile(site_covs$edge_density_forest, c(0.25, 0.75), na.rm = TRUE)

panel_a2_data <- map_dfr(edge_q, function(ed) {
  df <- tibble(treecover = tree_seq, edge_density_forest = ed) %>%
    mutate(
      treecover_z = scale_cov(treecover, mean_treecover, sd_treecover),
      edge_density_z = scale_cov(edge_density_forest, mean_edgedens, sd_edgedens)
    )
  X.0 <- model.matrix(~ treecover_z * edge_density_z, data = df)
  pred <- predict(ms_model, X.0 = X.0, type = "occupancy", ignore.RE = FALSE)
  
  map_dfr(seq_len(nrow(df)), function(i) {
    psi_samples <- pred$psi.0.samples[, , i]
    tibble(
      treecover = df$treecover[i],
      edge_density_forest = ed,
      mean = mean(rowMeans(psi_samples)),
      lower = quantile(rowMeans(psi_samples), 0.025),
      upper = quantile(rowMeans(psi_samples), 0.975)
    )
  })
}) %>%
  mutate(edge_class = factor(edge_density_forest, labels = c("Low edge density", "High edge density")))

panel_a2 <- ggplot(panel_a2_data, aes(x = treecover, y = mean, color = edge_class, fill = edge_class)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA) +
  geom_line(size = 1) +
  scale_color_manual(values = c("steelblue4", "firebrick")) +
  scale_fill_manual(values = c("steelblue4", "firebrick")) +
  labs(title = "A2. Tree cover × fragmentation", x = "Forest cover (%)", y = "Community occupancy", color = NULL, fill = NULL) +
  theme_classic(base_size = 13) +
  theme(legend.position = c(.85, .1))

# ----- Trait-filtered panels (generic) -----
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
panel_b_data <- bind_rows(
  trait_panel(species_traits$accepted_bin[species_traits$home_range_km2 <= hr_quant[1]], "Small home range"),
  trait_panel(species_traits$accepted_bin[species_traits$home_range_km2 >= hr_quant[2]], "Large home range")
)

panel_b <- ggplot(panel_b_data, aes(x = treecover, y = mean, color = group, fill = group)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA) +
  geom_line(size = 1) +
  scale_color_manual(values = c("darkgreen", "goldenrod"), name = NULL) +
  scale_fill_manual(values = c("darkgreen", "goldenrod"), name = NULL) +
  labs(title = "B. Trait filtering (Home range)", x = "Forest cover (%)", y = "Predicted occupancy") +
  theme_classic(base_size = 13) +
  theme(legend.position = "top")

# Panel C: Body mass
mass_quant <- quantile(species_traits$adult_body_mass_g, c(0.25, 0.75), na.rm = TRUE)
panel_c_data <- bind_rows(
  trait_panel(species_traits$accepted_bin[species_traits$adult_body_mass_g <= mass_quant[1]], "Small body size"),
  trait_panel(species_traits$accepted_bin[species_traits$adult_body_mass_g >= mass_quant[2]], "Large body size")
)

panel_c <- ggplot(panel_c_data, aes(x = treecover, y = mean, color = group, fill = group)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA) +
  geom_line(size = 1) +
  scale_color_manual(values = c("darkgreen", "goldenrod"), name = NULL) +
  scale_fill_manual(values = c("darkgreen", "goldenrod"), name = NULL) +
  labs(title = "C. Trait filtering (Body size)", x = "Forest cover (%)", y = "Predicted occupancy") +
  theme_classic(base_size = 13) +
  theme(legend.position = "top")

# ----- Combine all panels -----
final_fig <- (panel_a1 | panel_a2) / (panel_b | panel_c)
final_fig


