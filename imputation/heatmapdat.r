# Generate plotting data for heat maps

# test. In full version use all 999.
comparison_df_list <- list(comparisondf1, comparisondf2)

comparison_df_all <- do.call(rbind, comparison_df_list)

comparison_means <- comparison_df_all %>%
  mutate(delta = abs((true_trait-imputed_trait)/imputed_trait)) %>%
  group_by(trait_id, species_id) %>%
  summarize(delta = mean(delta))

# Use comparison_means to draw heat map.