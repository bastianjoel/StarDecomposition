library(magrittr)
library(dplyr)
library(ggplot2)
library(readr)

amount<-3

filtered_data<-results %>%
  #filter(algorithm %in% c("tet", "boundary-lp")) %>%
  filter(result_components > 1)# %>%
  #filter((algorithm == "boundary-lp" & result_components < 3) |
  #         (algorithm != "boundary-lp"))
  
total_unique_algorithms <- n_distinct(filtered_data$algorithm)
  
data<-filtered_data# %>%
  #group_by(filename) %>%
  #filter(n_distinct(algorithm) == total_unique_algorithms) %>%
  #ungroup()

count(data)

# Ensure 'result_components' is numeric
data$result_components <- as.numeric(data$result_components)

# Find the minimum result_components for each filename
min_components_per_filename <- data %>%
  group_by(filename) %>%
  summarise(min_result_components = min(result_components, na.rm = TRUE))

# Find the maximum result_components for each filename
max_components_per_filename <- data %>%
  group_by(filename) %>%
  summarise(max_result_components = max(result_components, na.rm = TRUE))

count(max_components_per_filename)

# Join back with the original data to find the algorithm(s) that achieved the minimum
best_algorithms <- data %>%
  inner_join(min_components_per_filename, by = "filename") %>%
  filter(result_components == min_result_components) %>%
  distinct(filename, algorithm) # Use distinct to handle cases where multiple algorithms
# might achieve the same minimum for a filename

# Count how many times each algorithm was the best
best_algorithm_counts <- best_algorithms %>%
  group_by(algorithm) %>%
  summarise(count = n()) %>%
  arrange(desc(count))

new_legend_title <- "Algorithm" 

new_algorithm_labels <- c(
  "boundary" = "Dynamic enclosure vertex move",
  "boundary-lp" = "Component normal direction",
  "tet" = "Tetrahedron decompose"
)

# Create the bar graph
ggplot(best_algorithm_counts, aes(x = reorder(algorithm, -count), y = count, fill = algorithm)) +
  geom_bar(stat = "identity", color= "black") +
  theme_minimal() +
  labs(
    x = "",
    y = "Number of Times Best",
    fill = "Algorithm"
  ) +  
  scale_fill_manual(
    values = c("#f8766d", "#00b0f6", "#a39a42"), # Automatically use ggplot's default colors
    labels = new_algorithm_labels,
    name = new_legend_title
  ) +
  theme(axis.text.x = element_blank(),
        legend.position = "bottom") +
  geom_text(aes(label = count), vjust = -0.5, size = 3.5) # Add count labels on bars

# Print the table of best algorithm counts
print(best_algorithm_counts)

