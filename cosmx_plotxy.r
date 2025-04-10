# Load necessary libraries
library(Seurat)
library(ggplot2)
library(dplyr)

# Assuming 'combined_seurat' is your Seurat object already loaded
# Ensure that 'x_slide_mm', 'y_slide_mm', 'Celltypes', 'ROI', and 'Genotype' are present in metadata

# Step 1: Extract the metadata
metadata <- combined_seurat@meta.data

# Step 2: Check if the required columns are present
required_columns <- c("x_slide_mm", "y_slide_mm", "Celltypes", "Section", "Slide")
if (!all(required_columns %in% colnames(metadata))) {
  stop("One or more required columns are missing from the metadata.")
}

# Step 3: Create a data frame for plotting
plot_data <- data.frame(
  x = metadata$x_slide_mm,
  y = metadata$y_slide_mm,
  Celltype = metadata$Celltypes,
  Section = metadata$Section,
  Slide = metadata$Slide
)

# Remove cells with NA Celltypes if necessary
plot_data <- plot_data[!is.na(plot_data$Celltype), ]

# Step 4: Define custom colors for cell types
celltype_colors <- c(
  "Excitatory Neurons" = "#90EE90",    # Light green
  "Inhibitory Neurons" = "#FF00FF",    # Magenta
  "Astrocytes" = "#FFA500",            # Orange
  "Oligodendrocytes" = "#D8BFD8",      # Light purple
  "OPCs" = "#9932CC",                  # Bright purple
  "Microglia" = "#0000FF",             # Bright blue
  "Macrophages" = "#FFD700",           # Gold
  "Endothelial cells" = "#808080"      # Gray
)

# Ensure all cell types have assigned colors
unique_celltypes <- unique(plot_data$Celltype)
missing_colors <- setdiff(unique_celltypes, names(celltype_colors))
if (length(missing_colors) > 0) {
  # Assign default colors to missing cell types
  default_colors <- rainbow(length(missing_colors))
  names(default_colors) <- missing_colors
  celltype_colors <- c(celltype_colors, default_colors)
}

# Step 5: Get unique combinations of ROI and Genotype
combinations <- unique(plot_data[, c("Section", "Slide")])

# Step 6: Loop over each combination and generate plots
for (i in 1:nrow(combinations)) {
  current_section <- combinations$Section[i]
  current_slide <- combinations$Slide[i]
  
  # Filter data for the current combination
  subset_data <- plot_data %>%
    filter(Section == current_section & Slide == current_slide)
  
  # Check if there is data to plot
  if (nrow(subset_data) == 0) {
    next  # Skip if no data for this combination
  }
  
  # Create the scatter plot
  p <- ggplot(subset_data, aes(x = x, y = y, color = Celltype)) +
    geom_point(size = 1.5, alpha = 0.8) +  # Adjust point size as needed
    theme_minimal() +
    labs(
      title = paste("Cortex:\nSpatial Distribution of Cell Types\nSection:", current_section, "- Slide:", current_slide),
      x = "X Coordinate (mm)",
      y = "Y Coordinate (mm)"
    ) +
    scale_color_manual(values = celltype_colors) +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5)
    )
  
  # Define the filename for saving
  filename <- paste0("Celltypes_Spatial_Plot_Ctx_Section_", current_section, "_Slide_", current_slide, ".png")
  
  # Save the plot with increased size and resolution
  ggsave(
    filename = filename,
    plot = p,
    width = 10,        # Adjust width as needed
    height = 8,        # Adjust height as needed
    dpi = 300          # Increase dpi for higher resolution
  )
  
  # Optionally, print the plot to the screen
  print(p)
}


# Load necessary libraries
library(Seurat)
library(ggplot2)
library(dplyr)

# Assuming 'combined_seurat' is your Seurat object already loaded
# Ensure that 'x_slide_mm', 'y_slide_mm', 'Celltypes', 'ROI', and 'Genotype' are present in metadata

# Step 1: Extract the metadata
metadata <- combined_seurat@meta.data

# Step 2: Check if the required columns are present
required_columns <- c("x_slide_mm", "y_slide_mm", "Celltypes", "Section", "Slide")
if (!all(required_columns %in% colnames(metadata))) {
  stop("One or more required columns are missing from the metadata.")
}

# Step 3: Create a data frame for plotting
plot_data <- data.frame(
  x = metadata$x_slide_mm,
  y = metadata$y_slide_mm,
  Celltype = metadata$Celltypes,
  Section = metadata$Section,
  Slide = metadata$Slide
)

# Remove cells with NA Celltypes if necessary
plot_data <- plot_data[!is.na(plot_data$Celltype), ]

# Step 4: Get unique combinations of ROI and Genotype
combinations <- unique(plot_data[, c("Section", "Slide")])

# Step 5: Loop over each combination and generate plots
for (i in 1:nrow(combinations)) {
  current_section <- combinations$Section[i]
  current_slide <- combinations$Slide[i]
  
  # Filter data for the current combination - all cells
  subset_all_cells <- plot_data %>%
    filter(Section == current_section & Slide == current_slide)
  
  # Check if there is data to plot
  if (nrow(subset_all_cells) == 0) {
    next  # Skip if no data for this combination
  }
  
  # Step 5a: Create a new column 'Highlight' to distinguish microglia
  subset_all_cells <- subset_all_cells %>%
    mutate(Highlight = ifelse(Celltype == "Microglia", "Microglia", "Other"))
  
  # Step 5b: Reorder the data so that 'Other' cells are plotted first and 'Microglia' are plotted last
  subset_all_cells <- subset_all_cells %>%
    arrange(factor(Highlight, levels = c("Other", "Microglia")))
  
  # Step 5c: Define colors
  highlight_colors <- c("Microglia" = "#0000FF", "Other" = "gray80")
  
  # Step 5d: Define alpha values
  highlight_alpha <- c("Microglia" = 1, "Other" = 0.5)
  
  # Step 5e: Determine the axis limits from all cells
  x_limits <- range(subset_all_cells$x, na.rm = TRUE)
  y_limits <- range(subset_all_cells$y, na.rm = TRUE)
  
  # Step 5f: Create the scatter plot
  p <- ggplot(subset_all_cells, aes(x = x, y = y)) +
    geom_point(aes(color = Highlight, alpha = Highlight), size = 1.5) +
    theme_minimal() +
    labs(
      title = paste("Cortex:\nMicroglia Highlighted Among All Cell Types\nSection:", current_section, "- Slide:", current_slide),
      x = "X Coordinate (mm)",
      y = "Y Coordinate (mm)"
    ) +
    scale_color_manual(values = highlight_colors) +
    scale_alpha_manual(values = highlight_alpha) +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5)
    ) +
    xlim(x_limits) +
    ylim(y_limits)
  
  # Define the filename for saving the plot
  filename <- paste0("Microglia_Highlighted_Plot_Section_", current_section, "_Slide_", current_slide, ".png")
  
  # Save the plot with increased size and resolution
  ggsave(
    filename = filename,
    plot = p,
    width = 10,        # Adjust width as needed
    height = 8,        # Adjust height as needed
    dpi = 300          # Increase dpi for higher resolution
  )
  
  # Optionally, print the plot to the screen
  print(p)
}

