library(magrittr)
library(dplyr)
library(ggplot2)
library(readr)

data<-results

filtered_data <- data

total_unique_algorithms <- n_distinct(filtered_data$algorithm)

boxplot_data<-filtered_data %>%
  filter(feasible == 1) %>%
  group_by(filename) %>%
  filter(n_distinct(algorithm) == total_unique_algorithms) %>%
  ungroup()

new_algorithm_labels <- c(
  "boundary" = "Dynamic enclosure vertex move",
  "boundary-lp" = "Component normal direction",
  "tet" = "Tetrahedron decompose"
)

new_legend_title <- "Algorithm" 

ggplot() +
  # Add boxplots for the non-baseline algorithms
  geom_boxplot(data = boxplot_data, aes(x = filename, y = result_components, fill = algorithm)) +

  scale_fill_manual(
    values = c("#f8766d", "#00b0f6", "#a39a42"), # Automatically use ggplot's default colors
    labels = new_algorithm_labels,
    name = new_legend_title
  ) +
  labs(title = "",
       x = "",
       y = "Amount components") + 
  theme(axis.text.x = element_text(angle = 50, vjust = 1, hjust=1), legend.position = "bottom")

boxplot_data$time <- parse_number(boxplot_data$time)

summarized_data <- boxplot_data %>%
  filter(!is.na(time)) %>% 
  group_by(filename, algorithm) %>%
  summarise(mean_time = mean(time), .groups = 'drop')

ggplot(summarized_data, aes(x = filename, y = mean_time, fill = algorithm)) +
  geom_col(position = position_dodge(width = 0.8), color = "black") + # geom_col uses pre-computed y values, position_dodge separates bars
  scale_fill_manual(
    values = c("#f8766d", "#00b0f6", "#a39a42"), # Automatically use ggplot's default colors
    labels = new_algorithm_labels,
    name = new_legend_title
  ) +
  labs(x = "",
       y = "Mean Time",
       fill = "Algorithm") + # Label for the fill legend
  theme(axis.text.x = element_text(angle = 50, vjust = 1, hjust=1), legend.position = "bottom")
