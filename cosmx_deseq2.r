# Set working directory
setwd("/Users/katherineridley/Projects/CosMx/APP/Cortex Results")

# Load required libraries
library(Seurat)
library(dplyr)
library(DESeq2)
library(ggplot2)
library(ggrepel)

### ---------------------
### DESeq2 Analysis for All Condition Pairs per Cell Type
### ---------------------

print("Performing differential expression analysis using DESeq2...")

# Load combined Seurat object (if not already loaded)
combined_seurat <- readRDS("combined_seurat_C_uncompressed.RDS")

# Get unique cell types from metadata; drop NA.
cell_types <- unique(combined_seurat$Celltypes)
cell_types <- cell_types[!is.na(cell_types)]
print("Cell types to analyze:")
print(cell_types)

# Get unique conditions (assume Condition is a factor or character column)
conditions <- sort(unique(combined_seurat$Condition))
print("Unique conditions:")
print(conditions)

# Generate all unique pairs (combinations) of conditions:
condition_pairs <- combn(conditions, 2, simplify = FALSE)

# Create directory to save DESeq2 results and plots
results_folder <- "deseq2_results_per_celltype"
if (!dir.exists(results_folder)) {
  dir.create(results_folder)
}

# Initialize list to store DESeq2 results
de_results <- list()

# Loop over each cell type
for (current_cell_type in cell_types) {
  message("Analyzing cell type: ", current_cell_type)
  
  # Use WhichCells to subset cells based on Celltypes
  cells_to_use <- WhichCells(combined_seurat, expression = Celltypes == current_cell_type)
  if (length(cells_to_use) == 0) {
    message("  Skipping ", current_cell_type, " because it has no cells.")
    next
  }
  cells_of_type <- subset(combined_seurat, cells = cells_to_use)
  
  # Get the condition distribution for this cell type
  cond_table <- table(cells_of_type$Condition)
  message("Condition counts for ", current_cell_type, ":")
  print(cond_table)
  
  # Initialize sub-list for this cell type
  de_results[[current_cell_type]] <- list()
  
  # Loop over each pair of conditions
  for (pair in condition_pairs) {
    cond1 <- pair[[1]]  # alphabetically smaller
    cond2 <- pair[[2]]
    
    # Check if both conditions are present in this cell-type subset
    if (!(cond1 %in% names(cond_table)) || !(cond2 %in% names(cond_table))) {
      message("  Skipping DE for ", current_cell_type, " for pair ", cond1, " vs. ", cond2, ": one condition is missing.")
      next
    }
    
    # Prepare DESeq2 input
    # Extract counts (for all genes) from the RNA assay. Using the old slot argument.
    counts_mat <- as.matrix(GetAssayData(cells_of_type, assay = "RNA", slot = "counts"))
    
    # Create metadata for DESeq2 from the cell metadata
    # Ensure that Condition is a factor with at least these two levels (order alphabetically)
    metadata_dds <- data.frame(Condition = factor(cells_of_type$Condition, levels = sort(c(cond1, cond2))))
    
    # Subset to columns (cells) belonging only to the current pair:
    pair_cells <- which(metadata_dds$Condition %in% c(cond1, cond2))
    counts_pair <- counts_mat[, pair_cells]
    metadata_pair <- metadata_dds[pair_cells, , drop = FALSE]
    
    # Create DESeq2 dataset and run DESeq2 (using poscounts normalization)
    dds <- DESeqDataSetFromMatrix(countData = counts_pair, colData = metadata_pair, design = ~ Condition)
    dds <- estimateSizeFactors(dds, type = "poscounts")
    dds <- DESeq(dds)
    
    # Contrast: by convention, let's make the result "cond2 vs cond1"
    res <- results(dds, contrast = c("Condition", cond2, cond1))
    res <- as.data.frame(res)
    
    # Add gene names and cell type info
    res$gene <- rownames(res)
    res$cell_type <- current_cell_type
    res$comparison <- paste(cond2, "vs", cond1)
    
    # Store results in the de_results list
    de_results[[current_cell_type]][[paste(cond2, "vs", cond1, sep = "_")]] <- res
    
    # Save the DE results to a CSV file
    out_csv <- file.path(results_folder, paste0("results_", current_cell_type, "_", cond2, "_vs_", cond1, ".csv"))
    write.csv(res, file = out_csv, row.names = FALSE)
    
    # Volcano plot
    res <- res %>% mutate(Significant = ifelse(padj < 0.05 & abs(log2FoldChange) >= 0.25, "Yes", "No"))
    p <- ggplot(res, aes(x = log2FoldChange, y = -log10(padj))) +
      geom_point(aes(color = Significant), alpha = 0.6) +
      scale_color_manual(values = c("No" = "grey", "Yes" = "red")) +
      ggtitle(paste("Volcano Plot -", current_cell_type, "-", cond2, "vs", cond1)) +
      xlab("Log2 Fold Change") +
      ylab("-Log10 Adjusted P-value") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5))
    
    # Label top genes
    top_genes <- res %>% filter(Significant == "Yes") %>% top_n(10, wt = -padj)
    p <- p + geom_text_repel(data = top_genes, aes(label = gene), size = 3, max.overlaps = 10)
    
    # Save the volcano plot as PNG
    volcano_file <- file.path(results_folder, paste0("volcano_", current_cell_type, "_", cond2, "_vs_", cond1, ".png"))
    ggsave(filename = volcano_file, plot = p, width = 8, height = 6)
    message("  Saved volcano plot for ", current_cell_type, " ", cond2, " vs ", cond1, " in ", volcano_file)
  }
}

# Optionally, save the entire de_results object as an RDS file
saveRDS(de_results, file = file.path(results_folder, "de_results.rds"))

print("DESeq2 differential expression analysis completed for all cell types and condition pairs.")


### ---------------------
### UMAP Plotting
### ---------------------

# Overall UMAP colored by Celltypes with labels
umap_overall <- DimPlot(combined_seurat, reduction = "umap", group.by = "Celltypes", label = TRUE, repel = TRUE) +
  ggtitle("Overall UMAP: Celltypes")
ggsave(filename = "UMAP_Celltypes_overall.png", plot = umap_overall, width = 8, height = 6)

# UMAP by Condition colored by Celltypes
# Get unique conditions
conditions <- sort(unique(combined_seurat$Condition))
for (cond in conditions) {
  subset_cond <- subset(combined_seurat, subset = Condition == cond)
  p <- DimPlot(subset_cond, reduction = "umap", group.by = "Celltypes", label = TRUE, repel = TRUE) +
    ggtitle(paste("UMAP - Condition:", cond))
  out_file <- paste0("UMAP_Celltypes_", cond, ".png")
  ggsave(filename = out_file, plot = p, width = 8, height = 6)
}
