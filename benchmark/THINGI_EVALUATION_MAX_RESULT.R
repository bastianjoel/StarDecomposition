library(dplyr)
library(ggplot2)
library(forcats)

# Initialize an empty list to store results from each 'amount' iteration
all_amount_results <- list()

# Loop through 'amount' from 3 to 20 to generate data
for (amount_val in 2:20) {
  
  filtered_data <- results %>%
    filter(algorithm %in% c("tet", "boundary-lp")) %>%
    filter(result_components > 1) %>%
    filter((algorithm == "boundary-lp" & result_components <= amount_val) |
             (algorithm != "boundary-lp"))
  
  total_unique_algorithms <- n_distinct(filtered_data$algorithm)
  
  data_for_best_calc <- filtered_data %>%
    group_by(filename) %>%
    filter(n_distinct(algorithm) == total_unique_algorithms) %>%
    ungroup()
  
  if (nrow(data_for_best_calc) == 0) {
    message(paste("No data for amount =", amount_val, "after filtering. Skipping."))
    next
  }
  
  data_for_best_calc$result_components <- as.numeric(data_for_best_calc$result_components)
  
  min_components_per_filename <- data_for_best_calc %>%
    group_by(filename) %>%
    summarise(min_result_components = min(result_components, na.rm = TRUE))
    
  max_components_per_filename <- data_for_best_calc %>%
    group_by(filename) %>%
    summarise(max_result_components = max(result_components, na.rm = TRUE))
  
  best_algorithms <- data_for_best_calc %>%
    inner_join(min_components_per_filename, by = "filename") %>%
    #inner_join(max_components_per_filename, by = "filename") %>%
    filter(result_components == min_result_components) %>%
    #filter(result_components < max_result_components) %>%
    distinct(filename, algorithm)
  
  best_algorithm_counts <- best_algorithms %>%
    group_by(algorithm) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(amount = amount_val)
  
  all_amount_results[[as.character(amount_val)]] <- best_algorithm_counts
}

# Combine all results into a single data frame
combined_results <- bind_rows(all_amount_results)

new_legend_title <- "Algorithm" 

new_algorithm_labels <- c(
  "boundary" = "Dynamic enclosure vertex move",
  "boundary-lp" = "Component normal direction",
  "tet" = "Tetrahedron decompose"
)

ggplot(combined_results, aes(x = amount, y = count, fill = algorithm)) +
  # Use position_dodge to place bars side-by-side
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  coord_flip() + # Makes the bars horizontal
  theme_minimal() +
  labs(
    x = "Amount Filter Value",
    y = "Number of ",
    fill = "Algorithm"
  ) +
  scale_fill_manual(
    values = c("#f8766d", "#00b0f6"), # Use ggplot's default colors
    labels = new_algorithm_labels,
    name = new_legend_title
  ) +
  theme(
    axis.text.y = element_text(angle = 0, hjust = 1), # Keep y-axis labels horizontal
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  ) +
  # Add count labels for each individual bar
  geom_text(
    aes(label = count),
    position = position_dodge(width = 0.9), # Match dodge width with geom_bar
    hjust = -0.1, # Adjust horizontal position for horizontal bars
    size = 3.5
  ) +
  scale_x_continuous(
    breaks = seq(0, max(combined_results$count, na.rm = TRUE) + 1, by = 1), # Show every integer break
    expand = expansion(mult = c(0, 0)) # Add space for labels
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) # Add space for labels

