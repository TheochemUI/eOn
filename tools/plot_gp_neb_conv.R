#!/usr/bin/env Rscript

# Load necessary libraries
library(ggplot2)
library(readr)
library(patchwork)
library(dplyr)
library(scales)

# Function to create and display the plot
plot_gp_neb_conv <- function(csv_path) {
  # Read the data from the CSV file
  data <- read_csv(csv_path)

  # Calculate upper and lower limits for the energy variance
  data <- data %>% mutate(Upper = `MAE Energy` + `Energy Variance`, Lower = `MAE Energy` - `Energy Variance`)

  # Plot for MAE Energy with Energy Variance
  plot_mae_energy <- ggplot(data, aes(x = Iteration)) +
    geom_point(aes(y = `MAE Energy`)) +
    geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.1) +
    geom_hline(yintercept = 0.001, color = "red", linetype = "dashed") +
    scale_y_continuous(breaks = pretty_breaks(n = 10)) +
    scale_x_continuous(breaks = pretty_breaks(n = 10)) +
    labs(y = "MAE Energy", x = "Iteration")

  # Plot for RMSF CI
  plot_rmsf_ci <- ggplot(data, aes(x = Iteration, y = `RMSF CI`)) +
    geom_point() +
    geom_line() +
    geom_hline(yintercept = 0.003, color = "red", linetype = "dashed") +
    scale_y_continuous(breaks = pretty_breaks(n = 10)) +
    scale_x_continuous(breaks = pretty_breaks(n = 10)) +
    labs(y = "RMSF CI", x = "Iteration")

  # Combine plots side by side
  combined_plot <- plot_mae_energy + plot_rmsf_ci + plot_layout(ncol = 2)

  # Add a combined title
  combined_plot <- combined_plot + plot_annotation(title = "GP NEB CI Convergence Plots", theme = theme(plot.title = element_text(hjust = 0.5)))

  # Define the output file name based on the input CSV file name
  output_file_name <- gsub("csv$", "pdf", basename(csv_path))

  # Save the combined plot to a PDF file
  ggsave(output_file_name, combined_plot, device = "pdf", width = 11, height = 8.5)
}

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if a CSV file path is provided
if (length(args) == 0) {
  stop("No CSV file path provided. Usage: Rscript plot_gp_neb_conv.R <path_to_csv_file>", call. = FALSE)
}

# Call the function with the CSV file path
plot_gp_neb_conv(args[1])
