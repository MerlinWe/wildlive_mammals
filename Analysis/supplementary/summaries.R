# Example: Slice at low and high edge density
slice_lowED <- comm_occ %>% filter(inside_hull, edge_density_forest < quantile(edge_density_forest, 0.2))
slice_highED <- comm_occ %>% filter(inside_hull, edge_density_forest > quantile(edge_density_forest, 0.8))

# Linear slope estimates along treecover
lowED_fit <- lm(mean ~ treecover, data = slice_lowED)
highED_fit <- lm(mean ~ treecover, data = slice_highED)

summary(lowED_fit)$coef
summary(highED_fit)$coef

comm_occ %>%
  filter(inside_hull) %>%
  arrange(desc(mean)) %>%
  slice(1)  # maximum occupancy point

# Predict from additive model
X_add <- model.matrix(~ treecover_z + edge_density_z, data = pdp_grid)
pred_add <- predict(ms_model, X.0 = X_add, type = "occupancy", ignore.RE = FALSE)

comm_add <- map_dfr(seq_len(nrow(pdp_grid)), function(i) {
  psi_samples <- pred_add$psi.0.samples[, , i]
  row_mean <- rowMeans(psi_samples)
  tibble(mean_add = mean(row_mean))
})

# Combine and compute interaction effect
comm_occ$mean_add <- comm_add$mean_add
comm_occ$interaction_effect <- comm_occ$mean - comm_occ$mean_add

# Quantify non-additivity
summary(comm_occ$interaction_effect)




# General PDP summary function
pdp_summary_for_trait <- function(trait_df, pdp_fn, trait_col, lower_label = "Low", upper_label = "High", q_low = 0.25, q_high = 0.75) {
  q_vals <- quantile(trait_df[[trait_col]], probs = c(q_low, q_high), na.rm = TRUE)
  
  # Apply the PDP function for both groups
  pdp_low <- pdp_fn(trait_df$accepted_bin[trait_df[[trait_col]] <= q_vals[1]], lower_label)
  pdp_high <- pdp_fn(trait_df$accepted_bin[trait_df[[trait_col]] >= q_vals[2]], upper_label)
  
  # Combine
  pdp_all <- bind_rows(pdp_low, pdp_high)
  
  # Summary extraction
  groups <- unique(pdp_all$group)
  for (g in groups) {
    summarize_pdp(pdp_all, g)
  }
  
  # Gap at high forest cover
  gap <- pdp_all %>%
    filter(treecover == max(treecover)) %>%
    group_by(treecover) %>%
    summarize(gap = diff(mean), .groups = "drop")
  
  cat("\nOccupancy gap at high forest cover (", upper_label, " – ", lower_label, "): ", round(gap$gap, 3), "\n", sep = "")
  
  return(pdp_all)
}

# Run for body mass
pdp_bodymass <- pdp_summary_for_trait(
  trait_df = species_traits,
  pdp_fn = trait_panel,
  trait_col = "body_mass_g",
  lower_label = "Small-bodied",
  upper_label = "Large-bodied"
)

# Run for diet breadth
pdp_diet <- pdp_summary_for_trait(
  trait_df = species_traits,
  pdp_fn = trait_panel,
  trait_col = "diet_breadth",
  lower_label = "Narrow diet",
  upper_label = "Broad diet"
)

# Run for habitat breadth
pdp_habitat <- pdp_summary_for_trait(
  trait_df = species_traits,
  pdp_fn = trait_panel,
  trait_col = "habitat_breadth",
  lower_label = "Habitat specialist",
  upper_label = "Habitat generalist"
)









#Function to extract PDP summaries for a group
summarize_pdp <- function(data, group_name) {
  group_data <- data %>% filter(group == group_name)
  
  start_tree <- min(group_data$treecover)
  end_tree   <- max(group_data$treecover)
  
  start_val <- group_data$mean[group_data$treecover == start_tree]
  end_val   <- group_data$mean[group_data$treecover == end_tree]
  
  delta_psi <- end_val - start_val
  
  start_CI <- group_data %>%
    filter(treecover == start_tree) %>%
    select(lower, upper)
  
  end_CI <- group_data %>%
    filter(treecover == end_tree) %>%
    select(lower, upper)
  
  cat("\n---", group_name, "---\n")
  cat("Start tree cover (%):", start_tree, "\n")
  cat("End tree cover (%):", end_tree, "\n")
  cat("Predicted occupancy at low forest cover:", round(start_val, 3), "\n")
  cat("Predicted occupancy at high forest cover:", round(end_val, 3), "\n")
  cat("Δψ (change in occupancy):", round(delta_psi, 3), "\n")
  cat("95% CI at low forest cover: [", round(start_CI$lower, 3), ",", round(start_CI$upper, 3), "]\n")
  cat("95% CI at high forest cover: [", round(end_CI$lower, 3), ",", round(end_CI$upper, 3), "]\n")
}

# Apply to both groups
summarize_pdp(hr_occ, "Small home range")
summarize_pdp(hr_occ, "Large home range")

# Optional: also print occupancy gap at high forest cover
gap_at_high_cover <- hr_occ %>%
  filter(treecover == max(treecover)) %>%
  group_by(treecover) %>%
  summarize(gap = diff(mean), .groups = "drop")

cat("\nOccupancy gap at high forest cover (Large – Small):", round(gap_at_high_cover$gap, 3), "\n")