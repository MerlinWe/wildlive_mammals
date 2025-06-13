#########################################################################################################
##### Trait-Filtered Responses of Mammal Communities to Land-Use Change in a Neotropical Dry Forest #####
#########################################################################################################

rm(list = ls()); gc()

library(patchwork)
library(ggh4x)
library(spOccupancy)
library(broom)
library(scales)
library(tidyverse)
select <- dplyr::select

# read input data
captures <- read_csv("/Users/serpent/Documents/Senckenberg/WildLive/Mammals/Code/Data/species_data.csv")
camtraps <- read_csv("/Users/serpent/Documents/Senckenberg/WildLive/Mammals/Code/Data/camtraps_clean.csv")
camop <- as.matrix(read_csv("/Users/serpent/Documents/Senckenberg/WildLive/Mammals/Code/Data/camop_problem.csv"))
covariates <- read_csv("/Users/serpent/Documents/Senckenberg/WildLive/Mammals/Code/Data/forest_covariates.csv")

# ========= Data Preparation ==========

# clean species data for analysis
species <- captures %>%
  filter(Station %in% camtraps$Station) %>%
  filter(Category %in% c("artiodactyla", "carnivora", "marsupialia", "perissodactyla", "primates", "rodentia", "xenarthra")) %>%
  filter(DateTimeOriginal >= "2017-01-10 17:01:26") %>%
  filter(DateTimeOriginal <= "2023-12-31 23:59:59") %>%
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
  file = "/Users/serpent/Documents/Senckenberg/WildLive/Mammals/Code/Data/PanTHERIA.txt",
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
  mutate(
    treecover_z = scale(treecover)[, 1],
    edge_density_z = scale(edge_density_forest)[, 1],
    patch_density_z = scale(patch_density_forest)[, 1],
    shape_index_z = scale(shape_index_forest)[, 1]
  ) %>%
  arrange(match(site_year, dimnames(detection_array)[[2]])) %>% #
  select(treecover_z, edge_density_z, patch_density_z, shape_index_z) %>%
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

## Find best fragmentation metric with additive models, then see if we need an interaction:

# Fit the model with edge density
ms_model_edge <- msPGOcc(
  occ.formula = ~ treecover_z + edge_density_z,
  det.formula = ~1,
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

# Fit the model with patch density
ms_model_patch <- msPGOcc(
  occ.formula = ~ treecover_z + patch_density_z,
  det.formula = ~1,
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

# Fit the model with shape index
ms_model_shape <- msPGOcc(
  occ.formula = ~ treecover_z + shape_index_z,
  det.formula = ~1,
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

# See which model is best
waic_edge <- spOccupancy::waicOcc(ms_model_edge)
waic_patch <- spOccupancy::waicOcc(ms_model_patch)
waic_shape <- spOccupancy::waicOcc(ms_model_shape)

waic_edge
waic_patch
waic_shape

cor(site_covs$patch_density_z, site_covs$edge_density_z)
cor(site_covs$patch_density_z, site_covs$shape_index_z)
cor(site_covs$edge_density_z, site_covs$shape_index_z)

cor(site_covs$treecover_z, site_covs$edge_density_z)
cor(site_covs$treecover_z, site_covs$patch_density_z)
cor(site_covs$treecover_z, site_covs$shape_index_z)

# the model with edge and patch variables are equally good based on the WAIC, the shape index performs worse.
# because we definitely want to keep tree cover in the model we continue with edge density since this is less
# correlated with tree cover than patch density. We can check if including shape index additively does any good:

ms_model_edge_and_shape <- msPGOcc(
  occ.formula = ~ treecover_z + edge_density_z + shape_index_z,
  det.formula = ~1,
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

waic_edge_shape <- spOccupancy::waicOcc(ms_model_edge_and_shape)
waic_edge_shape
waic_edge

# No, much worse. treecover + edge density is the top model. We now check if we need an interaction term.

ms_model_edge_interaction <- msPGOcc(
  occ.formula = ~ treecover_z * edge_density_z,
  det.formula = ~1,
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

waic_edge_interaction <- spOccupancy::waicOcc(ms_model_edge_interaction)
waic_edge_interaction["WAIC"] - waic_edge["WAIC"]

# The interaction model is worse by over 8 WAIC points than the additive model.
# That’s support to drop the interaction between treecover_z and edge_density_z, we drop it for parsimony.

ms_model <- ms_model_edge
rm(ms_model_edge_interaction, ms_model_patch, ms_model_shape, ms_model_edge_and_shape, ms_model_edge) # remove bad models

summary(ms_model, level = "community")
summary(ms_model, level = "species")
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
  filter(covariate %in% c("treecover_z", "edge_density_z")) %>%
  left_join(species_traits, by = c("species" = "accepted_bin")) %>%
  dplyr::mutate(
    log_mass = log(adult_body_mass_g),
    log_range = log(home_range_km2),
    diet_breadth = as.numeric(diet_breadth),
    activity_cycle_coded = dplyr::recode(activity_cycle, `1` = -1, `3` = 0, `2` = 1),
    habitat_breadth_coded = scales::rescale(habitat_breadth, to = c(-1, 1))
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

# Run trait ~ coefficient regressions across covariates
trait_vars <- c(
  "log_mass", "log_range", "diet_breadth",
  "activity_cycle_coded", "habitat_breadth_coded"
)

results <- map_dfr(unique(coef_traits$covariate), function(cov) {
  dat <- coef_traits %>% filter(covariate == cov)

  map_dfr(trait_vars, function(trait) {
    formula <- as.formula(paste0("Mean ~ ", trait))
    mod <- lm(formula, data = dat)

    # Plot diagnostic plots
    par(mfrow = c(2, 2)) # 2x2 plot layout
    plot(mod, main = paste("Diagnostics:", trait, "on", cov))

    # Return the tidy results
    tidy(mod) %>%
      filter(term != "(Intercept)") %>%
      mutate(covariate = cov, trait = trait)
  })
})

# Clean up for plotting
results <- results %>%
  distinct(covariate, trait, .keep_all = TRUE) %>%
  mutate(trait = dplyr::recode(trait,
    "log_mass" = "Body mass",
    "log_range" = "Home range size",
    "diet_breadth" = "Diet breadth",
    "activity_cycle_coded" = "Diurnality",
    "habitat_breadth_coded" = "Habitat breadth"
  ))

results$covariate <- factor(results$covariate,
  levels = c("treecover_z", "edge_density_z")
)

# Add significance flags for p-values
results <- results %>%
  mutate(sig_label = case_when(
    p.value < 0.001 ~ "***",
    p.value < 0.01 ~ "**",
    p.value < 0.05 ~ "*",
    TRUE ~ ""
  ))

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
    "treecover_z" = "Forestcover"
  ))

beta_comm$term <- factor(beta_comm$term,
  levels = c("Forestcover", "Edge Density")
)

# Community level plot
comm_plot <- ggplot(beta_comm, aes(x = term, y = Mean, ymin = `2.5%`, ymax = `97.5%`)) +
  geom_pointrange(aes(colour = term), size = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "(a.) Community-level responses",
    x = NULL, y = "Effect size (logit scale)"
  ) +
  scale_colour_viridis_d(option = "magma", end = .8) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(size = 10),
    title = element_text(size = 10),
    legend.position = "none"
  )

# Individual trait plots per covariate
traits_treecover <- results %>%
  filter(covariate == "treecover_z") %>%
  ggplot(aes(
    x = trait, y = estimate,
    ymin = estimate - std.error, ymax = estimate + std.error,
    color = trait
  )) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_pointrange(position = position_dodge(width = 0.5), size = 0.8) +
  geom_text(aes(label = sig_label, y = estimate + std.error + 0.02),
    size = 4, vjust = 0, show.legend = FALSE
  ) +
  scale_color_viridis_d(end = .95) +
  labs(
    title = "(b.) Trait responses - Tree Cover",
    x = NULL, y = "Effect size (slope estimate ± SE)"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(size = 10),
    title = element_text(size = 10),
    legend.position = "none"
  )

traits_fragmentation <- results %>%
  filter(covariate == "edge_density_z") %>%
  ggplot(aes(
    x = trait, y = estimate,
    ymin = estimate - std.error, ymax = estimate + std.error,
    color = trait
  )) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_pointrange(position = position_dodge(width = 0.5), size = 0.8) +
  geom_text(aes(label = sig_label, y = estimate + std.error + 0.02),
    size = 4, vjust = 0, show.legend = FALSE
  ) +
  scale_color_viridis_d(end = .95) +
  labs(
    title = "(c.) Trait responses - Fragmentation (Edge Density)",
    x = NULL, y = "Effect size (slope estimate ± SE)"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(size = 10),
    title = element_text(size = 10),
    legend.position = "none"
  )

effects_plot <- (comm_plot | (traits_treecover / traits_fragmentation)) +
  plot_layout(widths = c(1.3, 2))

ggsave("/Users/serpent/Documents/Senckenberg/WildLive/Mammals/Code/Output/Figures/fig2.png",
  effects_plot,
  width = 200,
  height = 140,
  units = "mm",
  dpi = 400
)

## ========== PDP Analysis ==========

# Define pdp sequences
tree_seq <- seq(-2, 2, length.out = 100)
frag_seq <- seq(-2, 2, length.out = 100)

# Define helper function to fit pdp curves
get_pdp <- function(predictor, predictor_seq, species_ids, species_names, ms_model) {
  cov_data <- data.frame(
    treecover_z = if (predictor == "treecover_z") predictor_seq else rep(0, length(predictor_seq)),
    edge_density_z = if (predictor == "edge_density_z") predictor_seq else rep(0, length(predictor_seq))
  )
  X.0 <- model.matrix(~ treecover_z + edge_density_z, data = cov_data)
  pred <- predict(ms_model, X.0 = X.0, type = "occupancy", ignore.RE = FALSE)
  
  map_dfr(seq_along(species_ids), function(i) {
    psi_samples <- pred$psi.0.samples[, species_ids[i], ]
    data.frame(
      species = species_names[species_ids[i]],
      predictor = predictor,
      value = predictor_seq,
      mean = apply(psi_samples, 2, mean),
      lower = apply(psi_samples, 2, quantile, 0.025),
      upper = apply(psi_samples, 2, quantile, 0.975)
    )
  })
}

# Get Species IDs
species_names <- gsub("treecover_z-", "", grep("^treecover_z-", colnames(ms_model$beta.samples), value = TRUE))
species_ids <- seq_along(species_names)

# Get PDPs and prep for plotting (supplementary)
pdp_tree <- get_pdp("treecover_z", tree_seq, species_ids, species_names, ms_model)
pdp_frag <- get_pdp("edge_density_z", frag_seq, species_ids, species_names, ms_model)

pdp_all <- bind_rows(pdp_tree, pdp_frag) %>%
  mutate(
    species = tools::toTitleCase(gsub("_", " ", species)),
    predictor = dplyr::recode(predictor,
                       treecover_z = "Forest cover",
                       edge_density_z = "Edge Density"
    ),
    facet_label = paste0(species, "\n(", predictor, ")")
  )  %>%
  mutate(facet_label = paste(species, "\n(", predictor, ")", sep = ""))

# Plot with facet_wrap in a wide layout
pdp_all_plot <- ggplot(pdp_all, aes(x = value, y = mean, ymin = lower, ymax = upper, fill = predictor, colour = predictor)) +
  geom_ribbon(alpha = 0.2) +
  geom_line() +
  scale_fill_viridis_d(end = .8) +
  scale_colour_viridis_d(end = .8) +
  facet_wrap(~facet_label, nrow = 7, ncol = 6, scales = "free") +
  labs(
    x = "Covariate value", y = "Predicted occupancy probability") +
  theme_bw() +
  theme(
    text = element_text(size = 8),
    strip.text = element_text(size = 8),
    axis.text = element_text(size = 7),
    panel.spacing = unit(0.2, "lines"),
    legend.position = "none")

ggsave("/Users/serpent/Documents/Senckenberg/WildLive/Mammals/Code/Output/Figures/supp6.png",
       pdp_all_plot,
       width = 270,
       height = 320,
       units = "mm",
       dpi = 600)

## Get trait bins and aggregate them to get a smoothes trait pdp structure

# Bin traits
trait_bins <- species_traits %>%
  select(accepted_bin, adult_body_mass_g, home_range_km2, diet_breadth, habitat_breadth, activity_cycle) %>%
  mutate(
    body_mass_bin = case_when(
      adult_body_mass_g <= quantile(adult_body_mass_g, 1/3, na.rm = TRUE) ~ "Small-bodied",
      adult_body_mass_g >= quantile(adult_body_mass_g, 2/3, na.rm = TRUE) ~ "Large-bodied"
    ),
    home_range_bin = case_when(
      home_range_km2 <= quantile(home_range_km2, 1/3, na.rm = TRUE) ~ "Small range",
      home_range_km2 >= quantile(home_range_km2, 2/3, na.rm = TRUE) ~ "Wide range"
    ),
    diet_breadth_bin = case_when(
      diet_breadth <= quantile(diet_breadth, 1/3, na.rm = TRUE) ~ "Specialist",
      diet_breadth >= quantile(diet_breadth, 2/3, na.rm = TRUE) ~ "Generalist"
    ),
    habitat_breadth_bin = case_when(
      habitat_breadth <= quantile(habitat_breadth, 1/3, na.rm = TRUE) ~ "Narrow",
      habitat_breadth >= quantile(habitat_breadth, 2/3, na.rm = TRUE) ~ "Broad"
    ),
    diurnality_bin = case_when(
      activity_cycle == 1 ~ "Nocturnal",
      activity_cycle == 2 ~ "Diurnal"
    )
  ) %>%
  select(accepted_bin, ends_with("_bin")) %>%
  pivot_longer(cols = -accepted_bin, names_to = "trait", values_to = "trait_group") %>%
  filter(!is.na(trait_group))  # drop mid-tertile and cathemeral species

# Join PDP with traits and aggregate across trait groups
pdp_trait_grouped <- pdp_all %>%
  mutate(species_snake = tolower(gsub(" ", "_", species))) %>%
  left_join(trait_bins, by = c("species_snake" = "accepted_bin")) %>%
  filter(!is.na(trait_group), !is.na(trait), !is.na(predictor)) %>%
  group_by(trait, trait_group, predictor, value) %>%
  summarise(
    mean = mean(mean),
    lower = mean(lower),
    upper = mean(upper),
    .groups = "drop"
  ) %>%
  mutate(
    # Human-readable facet labels
    facet_label = case_when(
      trait == "body_mass_bin" & trait_group == "Small-bodied" ~ "(a.) Small-bodied species",
      trait == "body_mass_bin" & trait_group == "Large-bodied" ~ "(b.) Large-bodied species",
      trait == "diet_breadth_bin" & trait_group == "Specialist" ~ "(e.) Dietary specialist",
      trait == "diet_breadth_bin" & trait_group == "Generalist" ~ "(f.) Dietary generalist",
      trait == "habitat_breadth_bin" & trait_group == "Narrow" ~ "(g.) Habitat specialist",
      trait == "habitat_breadth_bin" & trait_group == "Broad" ~ "(h.) Habitat generalist",
      trait == "home_range_bin" & trait_group == "Small range" ~ "(c.) Small home range",
      trait == "home_range_bin" & trait_group == "Wide range" ~ "(d.) Large home range",
      trait == "diurnality_bin" & trait_group == "Nocturnal" ~ "(i.) Nocturnal species",
      trait == "diurnality_bin" & trait_group == "Diurnal" ~ "(k.) Diurnal species"
    ),
    predictor = dplyr::recode(predictor,
                       "Edge Density" = "Fragmentation",
                       "Forest cover" = "Tree cover"),
    predictor = factor(predictor, levels = c("Fragmentation", "Tree cover"))
  )

# set facet layout 
custom_facet_levels <- c(
  "(a.) Small-bodied species", "(b.) Large-bodied species", "(c.) Small home range", "(d.) Large home range",
  "(e.) Dietary specialist", "(f.) Dietary generalist",   "(g.) Habitat specialist", "(h.) Habitat generalist",
  "(i.) Nocturnal species", "(k.) Diurnal species"
)

pdp_trait_grouped <- pdp_trait_grouped %>%
  mutate(facet_label = factor(facet_label, levels = custom_facet_levels))

#  Plot figure
trait_pdp <- ggplot(pdp_trait_grouped, aes(x = value, y = mean, ymin = lower, ymax = upper)) +
  geom_ribbon(aes(fill = predictor), alpha = 0.1, color = NA) +
  geom_line(aes(color = predictor), linewidth = 1) +
  facet_wrap(~facet_label, nrow = 3, ncol = 4) +
  scale_color_manual(values = c(
    "Fragmentation" = "#D95F02",
    "Tree cover" = "#1B9E77"
  )) +
  scale_fill_manual(values = c(
    "Fragmentation" = "#D95F02",
    "Tree cover" = "#1B9E77"
  )) +
  labs(
    x = "Predictor value (z-score)",
    y = "Predicted occupancy probability"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 14),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    strip.background = element_rect(fill = "transparent"),
    panel.spacing = unit(0.4, "lines")
  )

ggsave("/Users/serpent/Documents/Senckenberg/WildLive/Mammals/Code/Output/Figures/fig3.png",
       trait_pdp,
       width = 300,
       height = 220,
       units = "mm",
       dpi = 300)


# Build a summary table describing the data input for this table 
traits_to_summarize <- c("adult_body_mass_g", "home_range_km2", "diet_breadth", "habitat_breadth")

# Function to summarize continuous traits
continuous_trait_summary <- map_dfr(traits_to_summarize, function(trait) {
  lower_cutoff <- quantile(species_traits[[trait]], probs = 1/3, na.rm = TRUE)
  upper_cutoff <- quantile(species_traits[[trait]], probs = 2/3, na.rm = TRUE)
  
  species_traits %>%
    select(accepted_bin, !!sym(trait)) %>%
    filter(!is.na(!!sym(trait))) %>%
    mutate(
      tertile = case_when(
        !!sym(trait) <= lower_cutoff ~ "Lower tertile",
        !!sym(trait) >= upper_cutoff ~ "Upper tertile",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(tertile)) %>%
    group_by(trait = trait, tertile) %>%
    summarise(
      n_species = n(),
      species = paste(sort(accepted_bin), collapse = ", "),
      lower_threshold = lower_cutoff,
      upper_threshold = upper_cutoff,
      .groups = "drop"
    ) %>%
    relocate(trait, tertile, lower_threshold, upper_threshold, n_species, species)
})

# Handle diurnality separately
diurnality_summary <- species_traits %>%
  select(accepted_bin, activity_cycle) %>%
  mutate(
    diurnality = case_when(
      activity_cycle == 1 ~ "Nocturnal",
      activity_cycle == 3 ~ "Diurnal",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(diurnality)) %>%
  group_by(trait = "activity_cycle", tertile = diurnality) %>%
  summarise(
    n_species = n(),
    species = paste(sort(accepted_bin), collapse = ", "),
    lower_threshold = NA,
    upper_threshold = NA,
    .groups = "drop"
  ) %>%
  relocate(trait, tertile, lower_threshold, upper_threshold, n_species, species)

# Combine all trait summaries
trait_summary_all <- bind_rows(continuous_trait_summary, diurnality_summary)
write.csv(trait_summary_all, "/Users/serpent/Documents/Senckenberg/WildLive/Mammals/Code/Output/Tables/supplementary7.csv", row.names = FALSE)

# Clean column names
trait_table <- results %>%
  select(
    Trait = trait,
    Covariate = covariate,
    Estimate = estimate,
    SE = std.error,
    p = p.value
  ) %>%
  mutate(across(where(is.numeric), round, 3)) %>%
  mutate(
    Covariate = dplyr::recode(Covariate,
      treecover_z = "Forest cover",
      aggregation_z = "Fragmentation" 
      )
  ) %>%
  mutate(sig = ifelse(p < 0.05, TRUE, FALSE))

# Combine Estimate ± SE into one column
trait_table <- trait_table %>%
  mutate(`Estimate (± SE)` = sprintf("%.2f ± %.2f", Estimate, SE)) %>%
  select(Covariate, Trait, `Estimate (± SE)`, p, sig)

# Pivot wider to have covariates as columns
trait_table_wide <- trait_table %>%
  select(-sig) %>% # Drop sig if not used in export
  pivot_wider(
    names_from = Covariate,
    values_from = c(`Estimate (± SE)`, p),
    names_sep = ": "
  ) %>%
  relocate(Trait)

write_csv(trait_table_wide, "/Users/serpent/Documents/Senckenberg/WildLive/Mammals/Code/Output/Tables/traits_table2.csv")

trait_table <- coef_traits %>%
  select(
    species, log_mass, log_range, diet_breadth,
    activity_cycle_coded, habitat_breadth_coded
  ) %>%
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

write_csv(trait_table, "/Users/serpent/Documents/Senckenberg/WildLive/Mammals/Code/Output/Tables/traits_summary.csv")


# ========== Model Diagnostics ==========

# Extract Rhat and ESS as named vectors
rhat_vec <- ms_model$rhat$beta
ess_vec <- ms_model$ESS$beta

length(rhat_vec)
names(rhat_vec)[1:5]

# Build a tibble from the beta samples summary
beta_df <- as.data.frame(summary(beta_samples)$statistics) %>%
  rownames_to_column("param") %>%
  filter(str_detect(param, "^(treecover_z|aggregation_z)-")) %>%
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
    Covariate = recode(Covariate,
      "treecover_z"   = "Forest cover",
      "aggregation_z" = "Fragmentation",
      "agriculture"   = "Agricultural land"
    ),
    Mean = round(Mean, 3),
    SD = round(SD, 3),
    Rhat = round(Rhat, 3),
    ESS = round(ESS, 0)
  ) %>%
  select(-Converged, -ESS_gt_1000)

write_csv(diagnostics_clean, "/Users/serpent/Documents/Senckenberg/WildLive/Mammals/Code/Output/Tables/model_diagnostics.csv")

# ===== END =====
