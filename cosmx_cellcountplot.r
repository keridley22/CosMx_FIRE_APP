# Load necessary libraries
library(Seurat)
library(dplyr)
library(tidyr)
library(geometry)  # For calculating convex hull area

# Assuming 'combined_seurat' is your Seurat object already loaded

# Step 1: Check if required columns are present
required_columns <- c("x_slide_mm", "y_slide_mm", "Celltypes", "ROI", "Slide", "Condition")
missing_columns <- setdiff(required_columns, colnames(combined_seurat@meta.data))
if (length(missing_columns) > 0) {
  stop(paste("Missing required metadata columns:", paste(missing_columns, collapse = ", ")))
}

# Step 2: Extract and prepare metadata
metadata <- combined_seurat@meta.data

# Remove rows with missing coordinates
metadata <- metadata %>%
  filter(!is.na(x_slide_mm) & !is.na(y_slide_mm))

# Step 3: Calculate ROI area for each combination of ROI and Genotype
roi_area_df <- metadata %>%
  group_by(ROI, Slide, Condition) %>%
  do({
    # Extract unique coordinates in mm
    coords <- data.frame(x = .$x_slide_mm, y = .$y_slide_mm)
    coords <- coords[!duplicated(coords), ]
    
    # Check if there are at least 3 unique points to form a polygon
    if (nrow(coords) >= 3) {
      # Calculate convex hull
      hull_indices <- chull(coords)
      hull_coords <- coords[hull_indices, ]
      # Calculate area of convex hull polygon in mm²
      area_mm2 <- abs(sum((hull_coords$x * c(hull_coords$y[-1], hull_coords$y[1])) -
                          (hull_coords$y * c(hull_coords$x[-1], hull_coords$x[1])))) / 2
    } else {
      area_mm2 <- 0  # Area is zero if less than 3 points
    }
    data.frame(ROI_Area_mm2 = area_mm2)
  }) %>%
  ungroup()

# Step 4: Add ROI_Area_mm2 to the metadata
metadata <- metadata %>%
  left_join(roi_area_df, by = c("ROI", "Slide", "Condition"))

# Update the Seurat object's metadata
combined_seurat@meta.data <- metadata

# Step 5: Count the number of cells for each cell type in 'Celltypes'
celltype_counts_df <- metadata %>%
  filter(!is.na(Celltypes)) %>%
  group_by(ROI, Slide, Condition, Celltype = Celltypes) %>%
  summarise(Cell_Count = n(), .groups = "drop")  # Add .groups = "drop"

# Step 6: Pivot the data to have one row per ROI and Genotype combination
celltype_counts_pivot <- celltype_counts_df %>%
  pivot_wider(names_from = Celltype, values_from = Cell_Count, values_fill = 0)

# Step 7: Merge counts with ROI area
final_df <- celltype_counts_pivot %>%
  left_join(roi_area_df, by = c("ROI", "Slide", "Condition"))

# Step 8: Calculate normalized counts (cells per mm²)
# Get the list of cell type columns
celltype_columns <- setdiff(colnames(final_df), c("ROI", "Slide", "Condition", "ROI_Area_mm2"))

# For each cell type, calculate the normalized count
for (celltype in celltype_columns) {
  norm_col_name <- paste(celltype, "count normalized to ROI area mm2", sep = " ")
  final_df[[norm_col_name]] <- final_df[[celltype]] / final_df$ROI_Area_mm2
}

# Step 9: Rearrange columns
ordered_columns <- c(
  "ROI", "Slide", "Condition", "ROI_Area_mm2",
  celltype_columns,
  paste(celltype_columns, "count normalized to ROI area mm2", sep = " ")
)

final_df <- final_df[, ordered_columns]

# Step 10: Save the final dataframe to a CSV file
write.csv(final_df, file = "celltype_counts_normalized_mm2.csv", row.names = FALSE)

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(wesanderson)

# Step 1: Load the data
data <- read.csv("celltype_counts_normalized_mm2.csv")

print(colnames(data))
# Step 2: Prepare the data for plotting
# Adjust the grep pattern based on your column names
normalized_cols <- grep("\\.count\\.normalized\\.to\\.ROI\\.area\\.mm2$", colnames(data), value = TRUE)


# Check that normalized_cols is not empty
if (length(normalized_cols) == 0) {
  stop("No normalized count columns found. Please check the column names in your data.")
}

plot_data <- data %>%
  select(ROI, Slide, Condition, all_of(normalized_cols)) %>%
  pivot_longer(
    cols = all_of(normalized_cols),
    names_to = "Celltype",
    values_to = "Normalized_Count"
  ) %>%
  mutate(Celltype = gsub("\\.count\\.normalized\\.to\\.ROI\\.area\\.mm2$", "", Celltype), Celltype = gsub("\\.", " ", Celltype)
)


# Step 3: Create the box plot with larger text and increased width
p <- ggplot(plot_data, aes(x = Celltype, y = Normalized_Count, fill = Condition)) +
  geom_boxplot(color = "black", outlier.shape = NA) +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  theme_minimal(base_size = 14) +  # Increase base text size
  theme(
    panel.border = element_rect(color = "black", fill = NA),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.position = "right"
  ) +
  labs(
    title = "Normalized Cell Type Counts per ROI Area (mm2)",
    x = "Cell Type",
    y = "Normalized Count (cells per mm²)"
  ) +
  scale_fill_manual(values = wes_palette("Royal2", n = length(unique(plot_data$Condition)), type = "discrete"))

# Step 4: Save the plot with increased width and height
ggsave("normalized_celltype_counts_boxplot.png", plot = p, width = 12, height = 8, dpi = 300)

# Optionally, display the plot
print(p)


plot_data <- data %>%
  select(ROI, Slide, Condition, all_of(normalized_cols)) %>%
  pivot_longer(
    cols = all_of(normalized_cols),
    names_to = "Celltype",
    values_to = "Normalized_Count"
  ) %>%
  mutate(
    # Extract the cell type name from the column names
    Celltype = gsub("\\.count\\.normalized\\.to\\.ROI\\.area\\.mm2$", "", Celltype),
    # Replace periods with spaces for readability
    Celltype = gsub("\\.", " ", Celltype)
  )

# Remove any infinite or NaN values resulting from division by zero
plot_data <- plot_data %>%
  filter(is.finite(Normalized_Count))

# Set the order of the x-axis (Condition)
plot_data$Condition <- factor(plot_data$Condition, levels = c("WT", "FIRE", 'APP-WT', 'APP-FIRE'))

# Step 3: Create the faceted box plot
p <- ggplot(plot_data, aes(x = Condition, y = Normalized_Count, fill = Condition)) +
  geom_boxplot(color = "black", outlier.shape = NA) +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  facet_wrap(~ Celltype, scales = "free_y") +
  theme_minimal(base_size = 14) +  # Increase base text size
  theme(
    panel.border = element_rect(color = "black", fill = NA),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.position = "none"
  ) +
  labs(
    title = "Normalized Cell Type Counts per ROI Area (mm2)",
    x = "Condition",
    y = "Normalized Count (cells per mm2)"
  ) +
  scale_fill_manual(values = wes_palette("Royal2", n = length(unique(plot_data$Condition)), type = "discrete"))

# Step 4: Save the plot with increased width and height
ggsave("normalized_celltype_counts_facet_boxplot.png", plot = p, width = 12, height = 8, dpi = 300)

# Optionally, display the plot
print(p)

