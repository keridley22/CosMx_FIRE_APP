setwd("/Users/katherineridley/Projects/CosMx/APP")

# Load required libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(SingleR)
library(SingleCellExperiment)
library(celldex)
library(ggrepel)

# Start of the script
print("Starting the Seurat analysis script...")

# Load the two Seurat objects
print("Loading Seurat objects for genotypes A2 and B1...")
seurat_A2 <- readRDS("/Users/katherineridley/Projects/CosMx/APP/seuratObject_with_ROI_A2.RDS")
seurat_B1 <- readRDS("/Users/katherineridley/Projects/CosMx/APP/seuratObject_with_ROI_B1.RDS")
print("Seurat objects loaded successfully.")

# Function to subset Seurat object to include cells where ROI contains 'C' and save
subset_and_save_ROIs_with_C <- function(seurat_object, genotype_label) {
  # Check if the 'ROI' column exists
  if (!"ROI" %in% colnames(seurat_object@meta.data)) {
    stop(paste0("The 'ROI' column was not found in the metadata of genotype ", genotype_label, ". Please check the column names in your metadata."))
  }
  
  print(paste0("Analyzing ROI metadata for genotype ", genotype_label, "..."))
  
  # Print the number of values containing 'C', 'H', or NA
  num_C <- sum(grepl("C", seurat_object@meta.data$ROI))
  num_H <- sum(grepl("H", seurat_object@meta.data$ROI))
  num_NA <- sum(is.na(seurat_object@meta.data$ROI))
  
  print(paste("Number of cells with ROI containing 'C' in genotype", genotype_label, ":", num_C))
  print(paste("Number of cells with ROI containing 'H' in genotype", genotype_label, ":", num_H))
  print(paste("Number of cells with NA ROI in genotype", genotype_label, ":", num_NA))
  
  # Create a new metadata column indicating whether ROI contains 'C'
  seurat_object$ROI_contains_C <- !is.na(seurat_object$ROI) & grepl("C", seurat_object$ROI)
  
  seurat_object$Genotype <- genotype_label
  
  # Remove graphs to avoid mismatches
  seurat_object@graphs <- list()

  # Now subset using valid cell names
  cells_to_keep <- rownames(seurat_object@meta.data)[which(seurat_object@meta.data$ROI_contains_C)]
  cells_to_keep <- intersect(cells_to_keep, Cells(seurat_object))
  subset_object <- subset(seurat_object, cells = cells_to_keep)

  
  # Check if subset_object has any cells
  if (ncol(subset_object) == 0) {
    print(paste0("No cells with ROI containing 'C' found in genotype ", genotype_label, "."))
    return(NULL)
  } else {
    # Create a name for the subset object
    object_name <- paste0("seurat_", genotype_label, "_ROIs_with_C")
    
    # Save the subset Seurat object
    saveRDS(subset_object, file = paste0(object_name, ".RDS"))
    print(paste0("Saved subset Seurat object: ", object_name, ".RDS"))
    
    # Return the subset object
    return(subset_object)
  }
}

# Subset and save Seurat objects for ROIs containing 'C'
print("Subsetting Seurat object A2 for ROIs containing 'C'...")
subset_object_A2 <- subset_and_save_ROIs_with_C(seurat_A2, "A2")
print("Subsetting Seurat object B1 for ROIs containing 'C'...")
subset_object_B1 <- subset_and_save_ROIs_with_C(seurat_B1, "B1")

options(future.globals.maxSize = 20 * 1024^3)  # Allow up to 20 GiB

# Combine the subsetted Seurat objects into one for joint analysis
print("Combining the subsetted Seurat objects into one for joint analysis...")
subset_object_A2 <- readRDS("seurat_A2_ROIs_with_C.RDS")
subset_object_B1 <- readRDS("seurat_B1_ROIs_with_C.RDS")

# Create a list of subsetted Seurat objects (exclude NULLs)
subset_objects_list <- list()
if (!is.null(subset_object_A2)) {
  subset_objects_list[["A2"]] <- subset_object_A2
}
if (!is.null(subset_object_B1)) {
  subset_objects_list[["B1"]] <- subset_object_B1
}

# Check if there are any objects to merge
if (length(subset_objects_list) == 0) {
  stop("No Seurat objects with ROIs containing 'C' found in either genotype. Exiting script.")
} else if (length(subset_objects_list) == 1) {
  combined_seurat <- subset_objects_list[[1]]
  print("Only one genotype has cells with ROI containing 'C'. Proceeding with that genotype only.")
} else {
  # Merge the Seurat objects
  combined_seurat <- merge(
    x = subset_objects_list[[1]],
    y = subset_objects_list[-1],
    add.cell.ids = names(subset_objects_list)
  )
  print("Seurat objects merged successfully.")
}

#save combined seurat object
#saveRDS(combined_seurat, file = "combined_seurat_cortex.RDS")


print("Starting preprocessing with SCTransform...")

# Use SCTransform for normalization, which also performs variable gene selection 
combined_seurat <- SCTransform(combined_seurat, verbose = FALSE)

print("SCTransform normalization complete.")

# Optionally, perform PCA on the SCTransformed data; use more PCs
combined_seurat <- RunPCA(combined_seurat, verbose = FALSE, npcs = 30)
print("PCA complete.")

# Visualize PCA results: 
ElbowPlot(combined_seurat, ndims = 30)

# Find neighbors and clusters using a larger set of dimensions
combined_seurat <- FindNeighbors(combined_seurat, dims = 1:20)
combined_seurat <- FindClusters(combined_seurat, resolution = 0.5)
print("Clustering complete.")

# Run UMAP with tuned parameters.
# Adjust n.neighbors and min.dist to see improved separation.
combined_seurat <- RunUMAP(combined_seurat, dims = 1:20, n.neighbors = 30, min.dist = 0.3)
print("UMAP complete.")

# Plot UMAP colored by clusters (or other metadata)
DimPlot(combined_seurat, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE) +
  ggtitle("UMAP - Clusters (SCTransform, 20 PCs, n.neighbors=30, min.dist=0.3)")

# Optionally, plot feature plots for known markers to assess cell identities
FeaturePlot(combined_seurat, features = c("Gad1", "Pdgfra", "Cx3cr1"), reduction = "umap")


# Convert to SingleCellExperiment (required by SingleR)
sce <- as.SingleCellExperiment(combined_seurat)

# Load a reference dataset from celldex (MouseRNAseqData)
ref <- MouseRNAseqData()  # This is a reference dataset from celldex for mouse brain

# Find common genes between the Seurat object and the reference dataset
common_genes <- intersect(rownames(sce), rownames(ref))
sce <- sce[common_genes, ]
ref <- ref[common_genes, ]

# Run SingleR for cell type prediction
predictions <- SingleR(test = sce, ref = ref, labels = ref$label.main)

# Add the SingleR labels to the Seurat object's metadata
# Add the SingleR labels to the Seurat object's metadata
combined_seurat$SingleR_labels <- predictions$labels

# Check the updated metadata to confirm the new column
head(combined_seurat@meta.data)

# 1. Identify neurons based on SingleR labels
neuron_cells <- combined_seurat$SingleR_labels == "Neurons"

# 2. Initialize 'Celltypes' column in metadata (using a consistent name)
combined_seurat$Celltypes <- NA  # Use "Celltypes" consistently for all cell type assignments

# 3. Label 'Inhibitory Neurons' based on Gad1 expression
gad1_expression <- GetAssayData(combined_seurat, assay = "RNA", slot = "data")["Gad1", ]
inhibitory_neurons <- neuron_cells & (gad1_expression > 3)
combined_seurat$Celltypes[inhibitory_neurons] <- "Inhibitory Neurons"

# 4. Label the remaining neurons as 'Excitatory Neurons'
excitatory_neurons <- neuron_cells & !inhibitory_neurons
combined_seurat$Celltypes[excitatory_neurons] <- "Excitatory Neurons"

# 5. Label 'OPCs' based on Pdgfra expression
pdgfra_expression <- GetAssayData(combined_seurat, assay = "RNA", slot = "data")["Pdgfra", ]
opc_cells <- pdgfra_expression > 4
combined_seurat$Celltypes[opc_cells] <- "OPCs"

# 6. Assign other cell types based on SingleR labels
other_cell_types <- c("Astrocytes", "Oligodendrocytes", "Microglia", "Macrophages", "Endothelial cells")
for (ct in other_cell_types) {
  cells_of_type <- is.na(combined_seurat$Celltypes) & (combined_seurat$SingleR_labels == ct)
  combined_seurat$Celltypes[cells_of_type] <- ct
}

# Define the list of cell types you expect (used later for DE analysis)
cell_types <- c("Excitatory Neurons", "Inhibitory Neurons", "Astrocytes", 
                "Oligodendrocytes", "OPCs", "Microglia", "Macrophages", "Endothelial cells")

# --- Integrate Plaque and Microglia-based Condition Assignment ---
# Build a metadata data frame (meta_df) that includes the desired columns.
meta_df <- data.frame(
  Cell = colnames(combined_seurat),
  x_slide_mm = combined_seurat@meta.data$x_slide_mm,  # use $ to access columns
  y_slide_mm = combined_seurat@meta.data$y_slide_mm,
  Celltypes = combined_seurat@meta.data$Celltypes,     # note: using $ not @
  ROI = combined_seurat@meta.data$ROI,
  Slide = combined_seurat@meta.data$Slide,
  Section = combined_seurat@meta.data$Section,
  Cell_Type = Idents(combined_seurat),
  stringsAsFactors = FALSE
)

# Optionally, inspect the resulting meta_df:
head(meta_df)



library(dplyr)

# Assign PlaqueStatus based on Slide and Section
library(dplyr)

meta_df <- meta_df %>%
  mutate(Condition = case_when(
    Slide == "A2" & Section == "1" ~ "APP-WT",
    Slide == "A2" & Section == "2" ~ "APP-WT",
    Slide == "A2" & Section == "3" ~ "APP-FIRE",
    Slide == "A2" & Section == "4" ~ "WT",
    Slide == "A2" & Section == "5" ~ "WT",
    Slide == "A2" & Section == "6" ~ "FIRE",
    Slide == "B1" & Section == "1" ~ "FIRE",
    Slide == "B1" & Section == "2" ~ "FIRE",
    Slide == "B1" & Section == "3" ~ "APP-FIRE",
    Slide == "B1" & Section == "4" ~ "APP-FIRE",
    Slide == "B1" & Section == "5" ~ "APP-WT",
    Slide == "B1" & Section == "6" ~ "WT",
    TRUE ~ NA_character_
  ))

rownames(meta_df) <- meta_df$Cell
combined_seurat@meta.data <- meta_df

# Now check that the Condition column is present:
print(head(combined_seurat@meta.data$Condition))

# Set working directory and options
setwd("/Users/katherineridley/Projects/CosMx/APP")
options(future.globals.maxSize = 20 * 1024^3)  # allow up to 20 GiB for parallel processing

# Load required libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel)

# Create a folder for results if it doesn't exist
results_folder <- "Cortex Results"
if (!dir.exists(results_folder)) {
  dir.create(results_folder)
}

# Define the list of cell types that you want to examine
# (this list should match the values in the Celltypes column)
cell_types <- c("Excitatory Neurons", "Inhibitory Neurons", "Astrocytes", 
                "Oligodendrocytes", "OPCs", "Microglia", "Macrophages", "Endothelial cells")

# Initialize a list to store DE results
de_results <- list()

# Loop over each cell type
for (cell_type in cell_types) {
  message("Analyzing cell type: ", cell_type)
  
  # Subset combined_seurat to the current cell type
  cells_of_type <- subset(combined_seurat, subset = Celltypes == cell_type)
  
  # Skip if no cells found
  if (ncol(cells_of_type) == 0) {
    message("  Skipping ", cell_type, " because it has no cells.")
    next
  }
  
  # Create a list to store results for the current cell type
  de_results[[cell_type]] <- list()
  
  # Get the condition distribution in this cell type
  cond_table <- table(cells_of_type$Condition)
  message("Condition counts for ", cell_type, ":")
  print(cond_table)
  
  ## Comparison 1: WT vs non-WT
  # Check: at least one cell with condition "WT" and some with other conditions.
  if ("WT" %in% names(cond_table) && sum(cond_table[names(cond_table) != "WT"]) > 0) {
    # Create a temporary grouping: assign "WT" if Condition == "WT", otherwise "nonWT"
    cells_temp <- cells_of_type
    cells_temp$Group <- ifelse(cells_temp$Condition == "WT", "WT", "nonWT")
    Idents(cells_temp) <- cells_temp$Group
    
    markers <- FindMarkers(
      cells_temp,
      ident.1 = "WT",
      ident.2 = "nonWT",
      logfc.threshold = 0.25,
      min.pct = 0.1,
      layer = "data"  # Use the "data" layer instead of the deprecated slot
    )
    
    de_results[[cell_type]][["WT_vs_nonWT"]] <- markers
    
    # Create volcano plot for Comparison 1
    markers$gene <- rownames(markers)
    markers <- markers %>% mutate(Significant = ifelse(p_val_adj < 0.05 & abs(avg_log2FC) >= 0.25, "Yes", "No"))
    
    p <- ggplot(markers, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
      geom_point(aes(color = Significant), alpha = 0.6) +
      scale_color_manual(values = c("No" = "grey", "Yes" = "red")) +
      theme_minimal() +
      ggtitle(paste("Volcano Plot -", cell_type, "- WT vs. non-WT")) +
      xlab("Average Log2 Fold Change") +
      ylab("-log10 Adjusted P-value") +
      theme(plot.title = element_text(hjust = 0.5))
    
    volcano_file <- file.path(results_folder, paste0("volcano_", cell_type, "_WT_vs_nonWT.png"))
    ggsave(filename = volcano_file, plot = p, width = 8, height = 6)
    message("  Saved volcano plot for WT vs. non-WT in ", volcano_file)
  } else {
    message("  Skipping WT vs. non-WT for ", cell_type, " due to insufficient WT or nonWT cells.")
  }
  
  ## Comparison 2: FIRE vs. APP-FIRE
  if (all(c("FIRE", "APP-FIRE") %in% names(cond_table))) {
    Idents(cells_of_type) <- cells_of_type$Condition
    markers <- FindMarkers(
      cells_of_type,
      ident.1 = "FIRE",
      ident.2 = "APP-FIRE",
      logfc.threshold = 0.25,
      min.pct = 0.1,
      layer = "data"  # Specify the layer here
    )
    
    de_results[[cell_type]][["FIRE_vs_APP-FIRE"]] <- markers
    
    markers$gene <- rownames(markers)
    markers <- markers %>% mutate(Significant = ifelse(p_val_adj < 0.05 & abs(avg_log2FC) >= 0.25, "Yes", "No"))
    
    p <- ggplot(markers, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
      geom_point(aes(color = Significant), alpha = 0.6) +
      scale_color_manual(values = c("No" = "grey", "Yes" = "red")) +
      theme_minimal() +
      ggtitle(paste("Volcano Plot -", cell_type, "- FIRE vs. APP-FIRE")) +
      xlab("Average Log2 Fold Change") +
      ylab("-log10 Adjusted P-value") +
      theme(plot.title = element_text(hjust = 0.5))
    
    volcano_file <- file.path(results_folder, paste0("volcano_", cell_type, "_FIRE_vs_APP-FIRE.png"))
    ggsave(filename = volcano_file, plot = p, width = 8, height = 6)
    message("  Saved volcano plot for FIRE vs. APP-FIRE in ", volcano_file)
  } else {
    message("  Skipping FIRE vs. APP-FIRE for ", cell_type, " due to insufficient cells in one or both groups.")
  }
  
  ## Comparison 3: APP-WT vs. APP-FIRE
  if (all(c("APP-WT", "APP-FIRE") %in% names(cond_table))) {
    Idents(cells_of_type) <- cells_of_type$Condition
    markers <- FindMarkers(
      cells_of_type,
      ident.1 = "APP-WT",
      ident.2 = "APP-FIRE",
      logfc.threshold = 0.25,
      min.pct = 0.1,
      layer = "data"  # Using the correct layer argument
    )
    
    de_results[[cell_type]][["APP-WT_vs_APP-FIRE"]] <- markers
    
    markers$gene <- rownames(markers)
    markers <- markers %>% mutate(Significant = ifelse(p_val_adj < 0.05 & abs(avg_log2FC) >= 0.25, "Yes", "No"))
    
    p <- ggplot(markers, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
      geom_point(aes(color = Significant), alpha = 0.6) +
      scale_color_manual(values = c("No" = "grey", "Yes" = "red")) +
      theme_minimal() +
      ggtitle(paste("Volcano Plot -", cell_type, "- APP-WT vs. APP-FIRE")) +
      xlab("Average Log2 Fold Change") +
      ylab("-log10 Adjusted P-value") +
      theme(plot.title = element_text(hjust = 0.5))
    
    volcano_file <- file.path(results_folder, paste0("volcano_", cell_type, "_APP-WT_vs_APP-FIRE.png"))
    ggsave(filename = volcano_file, plot = p, width = 8, height = 6)
    message("  Saved volcano plot for APP-WT vs. APP-FIRE in ", volcano_file)
  } else {
    message("  Skipping APP-WT vs. APP-FIRE for ", cell_type, " due to insufficient cells in one or both groups.")
  }
}

# Optionally, save the DE results as an RDS file and/or CSV files.
saveRDS(de_results, file = file.path(results_folder, "de_results.rds"))

# Save each comparison per cell type as CSV
for (cell_type in names(de_results)) {
  for (comp in names(de_results[[cell_type]])) {
    result_df <- de_results[[cell_type]][[comp]]
    out_file <- file.path(results_folder, paste0("markers_", cell_type, "_", comp, ".csv"))
    write.csv(result_df, file = out_file, row.names = TRUE)
  }
}

message("Differential expression analysis and plotting complete!")

#saverds

saveRDS(combined_seurat, file = "combined_seurat_C_uncompressed.RDS", compress = FALSE)