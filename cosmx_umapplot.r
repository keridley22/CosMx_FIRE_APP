setwd("/Users/katherineridley/Projects/CosMx/APP/Cortex Results")

# Load required libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(SingleR)
library(SingleCellExperiment)
library(celldex)
library(ggrepel)
### ---------------------
### UMAP Plotting
### ---------------------

# Overall UMAP colored by Celltypes with labels

combined_seurat <- readRDS('combined_seurat_C_plaquedistance.RDS')
umap_overall <- DimPlot(combined_seurat, reduction = "umap", group.by = "Celltypes", label = TRUE, repel = TRUE) +
  ggtitle("Overall UMAP: Celltypes")
ggsave(filename = "UMAP_Celltypes_overall.png", plot = umap_overall, width = 8, height = 6)

# UMAP by Condition colored by Celltypes
# Get unique conditions
conditions <- sort(unique(combined_seurat$Condition))
for (cond in conditions) {
  subset_cond <- subset(combined_seurat, subset = Condition == cond)
  p <- DimPlot(subset_cond, reduction = "umap", group.by = "Celltypes", label = FALSE, repel = TRUE) +
    ggtitle(paste("UMAP - Condition:", cond))
  out_file <- paste0("UMAP_Celltypes_", cond, ".png")
  ggsave(filename = out_file, plot = p, width = 8, height = 6)
}