library(tidyverse)
library(patchwork)
library(ggthemes)
library(dplyr)
library(tibble)
library(stringr)
library(tidyr)
library(purrr)


# Helper function to scale covariates
scale_cov <- function(x, mean_x, sd_x) (x - mean_x) / sd_x

# Set up predictor sequences
tree_seq <- seq(min(site_covs$treecover, na.rm = TRUE),
                max(site_covs$treecover, na.rm = TRUE),
                length.out = 50)
mean_treecover <- mean(site_covs$treecover, na.rm = TRUE)
sd_treecover <- sd(site_covs$treecover, na.rm = TRUE)
mean_edgedens <- mean(site_covs$edge_density_forest, na.rm = TRUE)
sd_edgedens <- sd(site_covs$edge_density_forest, na.rm = TRUE)

# General trait PDP function
trait_panel <- function(species_group, label) {
  pdp_data <- tibble(
    treecover = tree_seq,
    edge_density_forest = mean_edgedens
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

# Trait names and quantile splits
trait_list <- c("adult_body_mass_g", "home_range_km2",  "diet_breadth", "habitat_breadth")
trait_titles <- c(
  "adult_body_mass_g" = "Body mass",
  "home_range_km2" = "Home range size",
  "diet_breadth" = "Diet breadth",
  "habitat_breadth" = "Habitat breadth"
)
plots <- list()

for (trait in trait_list) {
  quant <- quantile(species_traits[[trait]], c(0.25, 0.75), na.rm = TRUE)
  low_group <- species_traits$accepted_bin[species_traits[[trait]] <= quant[1]]
  high_group <- species_traits$accepted_bin[species_traits[[trait]] >= quant[2]]
  
  pdp <- bind_rows(
    trait_panel(low_group, "Low"),
    trait_panel(high_group, "High")
  )
  
  p <- ggplot(pdp, aes(x = treecover, y = mean, color = group, fill = group)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA) +
    geom_line(size = 1) +
    scale_color_manual(values = c("darkgreen", "goldenrod")) +
    scale_fill_manual(values = c("darkgreen", "goldenrod")) +
    labs(x = "Forest cover (%)", y = "Predicted occupancy",
         title = trait_titles[[trait]]) +
    theme_few(base_size = 8) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold", hjust = 0.5, size = 10),
      strip.text = element_text(size = 8)
    )
  
  plots[[trait]] <- p
}

# Diurnality: 1 = diurnal, 2 = nocturnal
diurnal_group <- species_traits$accepted_bin[species_traits$activity_cycle == 1]
nocturnal_group <- species_traits$accepted_bin[species_traits$activity_cycle == 2]

diurnality_pdp <- bind_rows(
  trait_panel(diurnal_group, "Diurnal"),
  trait_panel(nocturnal_group, "Nocturnal")
)

p_diurnality <- ggplot(diurnality_pdp, aes(x = treecover, y = mean, color = group, fill = group)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA) +
  geom_line(size = 1) +
  scale_color_manual(values = c("darkgreen", "goldenrod")) +
  scale_fill_manual(values = c("darkgreen", "goldenrod")) +
  labs(x = "Forest cover (%)", y = "Predicted occupancy", title = "Diurnality") +
  theme_few(base_size = 8) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", hjust = 0.5, size = 10),
    strip.text = element_text(size = 8)
  )

# Add to existing plot list
plots[["diurnality"]] <- p_diurnality

# Combine with patchwork
combined_traits_plot <- wrap_plots(plots, ncol = 2, guides = "collect") +
  plot_annotation(title = "Trait-based Partial Dependence Plots (PDPs)",
                  theme = theme(plot.title = element_text(size = 12)))


ggsave("/Users/merlin/Documents/Senckenberg/WildLive/Mammals/Code/Output/Figures/supp_8.png",
       last_plot(),
       width = 220,
       height = 160,
       units = "mm",
       dpi = 400)


# Define traits and labels
trait_info <- list(
  "Adult Body Mass (g)"   = "adult_body_mass_g",
  "Home Range (km2)"      = "home_range_km2",
  "Diet Breadth"          = "diet_breadth",
  "Habitat Breadth"       = "habitat_breadth"
)

# Helper to format species nicely
format_species <- function(x) {
  str_to_lower(x) %>% str_replace_all("_", " ") %>% sort() %>% paste(collapse = ", ")
}

# Continuous traits summary
trait_table <- map_dfr(names(trait_info), function(trait_label) {
  trait_col <- trait_info[[trait_label]]
  vals <- species_traits[[trait_col]]
  
  q <- quantile(vals, c(0.25, 0.75), na.rm = TRUE)
  
  low_spp <- species_traits$accepted_bin[vals <= q[1]]
  high_spp <- species_traits$accepted_bin[vals >= q[2]]
  
  tibble(
    Trait = trait_label,
    Quantile = c("Lower", "Upper"),
    Threshold = round(c(q[1], q[2]), 2),
    `# Species` = c(length(low_spp), length(high_spp)),
    Species = c(format_species(low_spp), format_species(high_spp))
  )
})

# Categorical trait: Activity cycle
activity_map <- c("Diurnal" = 1, "Nocturnal" = 2)

activity_table <- map_dfr(names(activity_map), function(label) {
  spp <- species_traits$accepted_bin[species_traits$activity_cycle == activity_map[[label]]]
  tibble(
    Trait = "Activity Cycle",
    Quantile = label,
    Threshold = "-",
    `# Species` = length(spp),
    Species = format_species(spp)
  )
})

# Combine both
trait_table$Threshold <- as.character(trait_table$Threshold)
final_trait_summary <- bind_rows(trait_table, activity_table)

print(final_trait_summary)
write_csv(final_trait_summary, "/Users/merlin/Documents/Senckenberg/WildLive/Mammals/Code/Output/Tables/trait_quantiles.csv")
