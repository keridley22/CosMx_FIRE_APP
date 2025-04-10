# Load required libraries
library(Seurat)
library(dplyr)

# Helper function to read plaque coordinates and compute min distance for each cell
assign_plaque_distance <- function(seurat_obj,
                                   plaque_dir_a2 = "/Users/katherineridley/Projects/CosMx/APP/A2_rois/",
                                   plaque_dir_b1 = "/Users/katherineridley/Projects/CosMx/APP/B1_rois/") {
  # Ensure we have the columns we need
  required_cols <- c("Slide", "Section", "x_slide_mm", "y_slide_mm")
  if (!all(required_cols %in% colnames(seurat_obj@meta.data))) {
    stop("Seurat metadata must contain 'Slide', 'Section', 'x_slide_mm', 'y_slide_mm' columns.")
  }
  
  # Initialize a new metadata column for plaque distance (NA by default)
  seurat_obj$Plaque_Distance <- NA_real_
  
  # Extract metadata into a data frame for easier manipulation
  meta_df <- seurat_obj@meta.data
  
  # Identify unique (Slide, Section) pairs
  unique_slide_sections <- meta_df %>% distinct(Slide, Section)
  
  # Loop over each unique combination of Slide and Section
  for (i in seq_len(nrow(unique_slide_sections))) {
    current_slide   <- unique_slide_sections$Slide[i]
    current_section <- unique_slide_sections$Section[i]
    
    # Construct CSV file path for plaque coordinates
    # e.g. "A2_S1_ROIs_transformed.csv"
    # If Slide is "A2", read from plaque_dir_a2; if "B1", read from plaque_dir_b1
    if (current_slide == "A2") {
      csv_file <- file.path(plaque_dir_a2,
                            paste0(current_slide, "_S", current_section, "_ROIs_transformed.csv"))
    } else if (current_slide == "B1") {
      csv_file <- file.path(plaque_dir_b1,
                            paste0(current_slide, "_S", current_section, "_ROIs_transformed.csv"))
    } else {
      # If there are other Slides not matching A2 or B1, skip or handle accordingly
      message("Slide ", current_slide, " does not match A2 or B1. Skipping.")
      next
    }
    
    if (!file.exists(csv_file)) {
      message("File not found: ", csv_file, " - skipping Slide ", current_slide,
              ", Section ", current_section)
      next
    }
    
    # Read plaque centroid coordinates
    plaque_coords <- read.csv(csv_file)
    # Expect columns: plaque_coords$x_mm, plaque_coords$y_mm
    
    # Get the indices of cells that match this Slide + Section
    cell_indices <- which(meta_df$Slide == current_slide & meta_df$Section == current_section)
    
    # For each cell in this subset, calculate min distance to any plaque centroid
    for (cell_idx in cell_indices) {
      x_cell <- meta_df$x_slide_mm[cell_idx]
      y_cell <- meta_df$y_slide_mm[cell_idx]
      
      # Compute Euclidean distances to each plaque centroid
      # Note that plaque_coords$x_mm and plaque_coords$y_mm are in mm
      distances <- sqrt((x_cell - plaque_coords$x_mm)^2 + (y_cell - plaque_coords$y_mm)^2)
      
      # Minimum distance
      meta_df$Plaque_Distance[cell_idx] <- min(distances)
    }
  }
  
  # Write updated metadata back to the Seurat object
  seurat_obj@meta.data <- meta_df
  
  return(seurat_obj)
}

seurat_cortex <- readRDS('/Users/katherineridley/Projects/CosMx/APP/Cortex Results/combined_seurat_C_uncompressed.RDS')
seurat_hippocampus <- readRDS('/Users/katherineridley/Projects/CosMx/APP/Hpc Results/combined_seurat_H_uncompressed.RDS')


seurat_cortex <- assign_plaque_distance(seurat_cortex)
seurat_hippocampus <- assign_plaque_distance(seurat_hippocampus)


