# Load required libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(DESeq2)
library(ggrepel)
setwd("/Users/katherineridley/Projects/CosMx/APP/Cortex Results")
### -----------------------------------------------------------
### 1. Analysis of Cell-Type Composition as a Function of Plaque Distance
### -----------------------------------------------------------
seurat_cortex <- readRDS('combined_seurat_plaquedistance_cleanup.RDS')
meta_df <- seurat_cortex@meta.data


# Filter out cells with NA Celltypes and only keep the two conditions of interest
composition_df <- meta_df %>%
  filter(!is.na(Celltypes) & Condition %in% c("APP-WT", "APP-FIRE")) %>%
  # Exclude Microglia & Macrophages for APP-FIRE condition
  filter(!(Condition == "APP-FIRE" & Celltypes == "Microglia & Macrophages"))

  # Set Condition as a factor with APP-WT first, then APP-FIRE
composition_df$Condition <- factor(composition_df$Condition, levels = c("APP-WT", "APP-FIRE"))

# Option A: Stacked Density Plot using continuous x-axis (0-100 µm)
p_density <- ggplot(composition_df, aes(x = Plaque_Distance * 1000, fill = Celltypes_r)) +
  geom_density(alpha = 0.7, position = "fill") +
  facet_wrap(~ Condition, ncol = 1) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) +
  xlab("Distance to Plaque (um)") +
  ylab("Proportion of cells") +
  ggtitle("Cell-type Composition by Continuous Plaque Distance and Condition") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = "Celltype_Composition_by_Continuous_Distance_and_Condition_density.png",
       plot = p_density, width = 10, height = 8)

# Option B: Stacked Histogram Plot with bins and continuous x-axis (0-100 µm)
# Define the number of bins; here we set breaks to have 5µm increments if desired,
# but since we want the x-axis to go to 100 µm, you can adjust the bin count accordingly.
p_hist <- ggplot(composition_df, aes(x = Plaque_Distance * 1000, fill = Celltypes_r)) +
  geom_histogram(bins = 50, position = "fill", color = "black") +
  facet_wrap(~ Condition, ncol = 1) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) +
  xlab("Distance to Plaque (um)") +
  ylab("Proportion of cells") +
  ggtitle("Cell-type Composition by Plaque Distance and Condition (Histogram)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = "Celltype_Composition_by_Continuous_Distance_and_Condition_hist.png",
       plot = p_hist, width = 10, height = 8)

print("Continuous-distance composition plots (0-100 µm) saved.")


# Option C: Individual Cell-type Plot Overlaid by Condition
# For each cell type, plot a smoothed line showing the fraction of cells at each plaque distance,
# with both conditions overlaid.
# We'll compute a kernel density estimate per condition for each cell type and then plot.
# One approach is to use geom_density with position="fill" but that gives stack proportions. 
# Instead, here we compute densities separately and then plot them on the same axis for direct comparison.
library(ggplot2)

# Get list of unique cell types from the filtered data
unique_ct <- unique(composition_df$Celltypes)

# Loop over each cell type and produce an overlaid density plot by Condition.
for (ct in unique_ct) {
  comp_ct <- composition_df %>%
    filter(Celltypes == ct)
  
  p_ct <- ggplot(comp_ct, aes(x = Plaque_Distance * 1000, color = Condition)) +
    geom_density(size = 1.2, adjust = 1.5) +
    scale_x_continuous(limits = c(0, 50), breaks = seq(0, 50, by = 5)) +
    xlab("Distance to Plaque (um)") +
    ylab("Density") +
    ggtitle(paste("Density Plot for", ct, "by Condition")) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  out_file <- paste0("Celltype_Density_", gsub(" ", "_", ct), ".png")
  ggsave(filename = out_file, plot = p_ct, width = 8, height = 6)
}

print("All cell-type composition plots (continuous x-axis) saved.")




