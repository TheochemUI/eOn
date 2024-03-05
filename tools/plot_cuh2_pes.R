#!/usr/bin/env Rscript

# Load necessary libraries
library(ggplot2)
library(readr)
library(dplyr)

# Function to create and display the plot
plot_cuh2_pes <- function(csv_path) {
  # Read the data from the CSV file
  data <- read_csv(csv_path)
  # Shift by the minimum
  data$Energy <- data$Energy - min(data$Energy)
  # Default clipping range
  clip_min <- 0
  clip_max <- 5

  # Calculate upper and lower limits for the Energy variance
  data <- data %>% mutate(Energy = replace(Energy, Energy > clip_max, clip_max)) %>% mutate(Energy = replace(Energy, Energy < clip_min, clip_min))

plot_cuh2_pes <- ggplot(data, ggplot2::aes(x = HH_Distance, y = HCu_Distance, z = Energy)) + ggplot2::geom_raster(interpolate = T, ggplot2::aes(fill = Energy)) +
    ggplot2::geom_contour(color = "white") +
    ggplot2::scale_fill_gradientn(colors = khroma::color("batlow")(10)) +
    ggplot2::labs(x = "H-H distance", y = "Cu-H2 distance", title = "CuH2 Potential Energy (True)")

  # Define the output file name based on the input CSV file name
  output_file_name <- gsub("csv$", "pdf", basename(csv_path))

  # Save the combined plot to a PDF file
  ggsave(output_file_name, plot_cuh2_pes, device = "pdf", width = 11, height = 8.5)
}

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if a CSV file path is provided
if (length(args) == 0) {
  stop("No CSV file path provided. Usage: Rscript plot_cuh2_pes.R <path_to_csv_file>", call. = FALSE)
}

# Call the function with the CSV file path
plot_cuh2_pes(args[1])
