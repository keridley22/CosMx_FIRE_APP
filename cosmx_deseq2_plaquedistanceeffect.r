# Set working directory
setwd("/Users/katherineridley/Projects/CosMx/APP/Cortex Results")

# Load required libraries
library(Seurat)
library(dplyr)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(crayon)  # For colored console output

### -----------------------------------------------------------
### 1. DESeq2 Analysis with Condition + Plaque_Distance (Continuous) per Cell Type
### -----------------------------------------------------------

# Define the cell types to analyze (these should match exactly the labels in your metadata)
cell_types <- c("Excitatory Neurons", "Inhibitory Neurons", "Astrocytes", 
                "Oligodendrocytes", "OPCs", "Microglia & Macrophages", "Endothelial cells")

# Create a folder for saving DESeq2 results and MA plots (combined analysis)
deseq_folder <- "deseq2_results_per_celltype_by_distance"
if (!dir.exists(deseq_folder)) {
  dir.create(deseq_folder)
}

# Initialize a list to store DESeq2 results for each cell type
de_results <- list()

for (ct in cell_types) {
  message("Analyzing differential expression for cell type: ", ct)
  
  # Subset cells belonging to current cell type using WhichCells to ensure we get cell names
  cells_to_use <- WhichCells(seurat_cortex, expression = Celltypes == ct)
  if (length(cells_to_use) == 0) {
    message("  Skipping ", ct, " because it has no cells.")
    next
  }
  subset_ct <- subset(seurat_cortex, cells = cells_to_use)
  
  # For DESeq2, only include cells from conditions with plaques (APP-WT and APP-FIRE)
  cells_pair <- WhichCells(subset_ct, expression = Condition %in% c("APP-WT", "APP-FIRE"))
  if (length(cells_pair) == 0) {
    message("  Skipping ", ct, " because no cells in APP-WT or APP-FIRE conditions.")
    next
  }
  subset_ct <- subset(subset_ct, cells = cells_pair)
  
  # Prepare count matrix and metadata for DESeq2
  counts_mat <- as.matrix(GetAssayData(subset_ct, assay = "RNA", slot = "counts"))
  meta_dds <- data.frame(
    Condition = factor(subset_ct$Condition, levels = c("APP-WT", "APP-FIRE")),
    Plaque_Distance = subset_ct$Plaque_Distance,
    row.names = colnames(subset_ct)
  )
  
  # Create DESeq2 dataset with design ~ Condition + Plaque_Distance
  dds <- DESeqDataSetFromMatrix(countData = counts_mat,
                                colData = meta_dds,
                                design = ~ Condition + Plaque_Distance)
  dds <- estimateSizeFactors(dds, type = "poscounts")
  dds <- DESeq(dds)
  
  # Extract the result for the continuous covariate Plaque_Distance
  res <- results(dds, name = "Plaque_Distance")
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  res_df$cell_type <- ct
  res_df$comparison <- "Plaque_Distance (continuous)"
  
  # Save DE results as CSV for this cell type
  out_csv <- file.path(deseq_folder, paste0("results_", gsub(" ", "_", ct), "_Plaque_Distance.csv"))
  write.csv(res_df, file = out_csv, row.names = FALSE)
  message("  Saved DESeq2 CSV for ", ct, " in ", out_csv)
  
  # Mark significance: padj < 0.05 and abs(log2FoldChange) >= 0.25
  res_df <- res_df %>% mutate(Significant = ifelse(padj < 0.05 & abs(log2FoldChange) >= 0.25, "Yes", "No"))
  
  # Print the names of significant genes in red and non-significant in grey
  sig_genes <- res_df %>% filter(Significant == "Yes") %>% arrange(padj)
  if(nrow(sig_genes) > 0) {
    message("Significant genes for ", ct, ":")
    for(i in 1:nrow(sig_genes)) {
      message(red(sig_genes$gene[i]))
    }
  } else {
    message("No significant genes found for ", ct)
  }
  
  # Create an MA plot for the effect of Plaque_Distance
  ma_plot <- ggplot(res_df, aes(x = baseMean, y = log2FoldChange)) +
    geom_point(aes(color = Significant), alpha = 0.6) +
    scale_x_log10() +
    scale_color_manual(values = c("No" = "grey", "Yes" = "red")) +
    ggtitle(paste("MA Plot:", ct, "- Effect of Plaque Distance")) +
    xlab("Mean Expression (log scale)") +
    ylab("Log2 Fold Change per mm") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  ma_file <- file.path(deseq_folder, paste0("MAPlot_", gsub(" ", "_", ct), ".png"))
  ggsave(filename = ma_file, plot = ma_plot, width = 8, height = 6)
  message("  Saved MA plot for ", ct, " in ", ma_file)
  
  # Add result to our list
  de_results[[ct]] <- res_df
}

# Print interpretation of the y-axis for MA plots:
message("Interpretation: In the MA plots, the y-axis (log2FoldChange for Plaque_Distance) represents the change in gene expression per unit (mm) increase in plaque distance. A negative log2FoldChange indicates that gene expression is higher closer to the plaque (i.e., as distance increases, expression decreases), whereas a positive log2FoldChange indicates that expression increases with distance.")

# Optionally, save the entire de_results object as an RDS file for later use
saveRDS(de_results, file = file.path(deseq_folder, "de_results_all_celltypes.rds"))

### -----------------------------------------------------------
### 3. Differential Expression with Plaque_Distance as Continuous Only (By Genotype)
### -----------------------------------------------------------
# Now, for each genotype (APP-WT and APP-FIRE), run a separate DE analysis using design ~ Plaque_Distance only.
genotypes <- c("APP-WT", "APP-FIRE")
deseq_dist_folder <- "deseq2_results_by_distance_only"
if (!dir.exists(deseq_dist_folder)) {
  dir.create(deseq_dist_folder)
}

de_results_distance <- list()

for (gen in genotypes) {
  message("Analyzing distance-only DE for genotype: ", gen)
  
  for (ct in cell_types) {
    # Subset to current cell type and current genotype
    cells_to_use <- WhichCells(seurat_cortex, expression = Celltypes == ct & Condition == gen)
    if (length(cells_to_use) == 0) {
      message("  Skipping ", ct, " for genotype ", gen, " because it has no cells.")
      next
    }
    subset_ct <- subset(seurat_cortex, cells = cells_to_use)
    
    # Prepare DESeq2 input: design ~ Plaque_Distance only
    counts_mat <- as.matrix(GetAssayData(subset_ct, assay = "RNA", slot = "counts"))
    meta_dds <- data.frame(
      Plaque_Distance = subset_ct$Plaque_Distance,
      row.names = colnames(subset_ct)
    )
    
    dds <- DESeqDataSetFromMatrix(countData = counts_mat,
                                  colData = meta_dds,
                                  design = ~ Plaque_Distance)
    dds <- estimateSizeFactors(dds, type = "poscounts")
    dds <- DESeq(dds)
    
    res <- results(dds, name = "Plaque_Distance")
    res_df <- as.data.frame(res)
    res_df$gene <- rownames(res_df)
    res_df$cell_type <- ct
    res_df$genotype <- gen
    res_df$comparison <- "Distance only"
    
    # Save result as CSV for this cell type and genotype
    out_csv <- file.path(deseq_dist_folder, paste0("results_", gsub(" ", "_", ct), "_", gen, "_DistanceOnly.csv"))
    write.csv(res_df, file = out_csv, row.names = FALSE)
    message("  Saved distance-only DE CSV for ", ct, " in genotype ", gen, " at ", out_csv)
    
    # Mark significance
    res_df <- res_df %>% mutate(Significant = ifelse(padj < 0.05 & abs(log2FoldChange) >= 0.25, "Yes", "No"))
    
    # Create an MA plot for distance-only analysis
    ma_plot <- ggplot(res_df, aes(x = baseMean, y = log2FoldChange)) +
      geom_point(aes(color = Significant), alpha = 0.6) +
      scale_x_log10() +
      scale_color_manual(values = c("No" = "grey", "Yes" = "red")) +
      ggtitle(paste("MA Plot:", ct, "-", gen, "- Distance Only")) +
      xlab("Mean Expression (log scale)") +
      ylab("Log2 Fold Change per mm (Distance effect)") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5))
    
    ma_file <- file.path(deseq_dist_folder, paste0("MAPlot_", gsub(" ", "_", ct), "_", gen, "_DistanceOnly.png"))
    ggsave(filename = ma_file, plot = ma_plot, width = 8, height = 6)
    message("  Saved distance-only MA plot for ", ct, " in genotype ", gen, " at ", ma_file)
    
    # Store results in our list
    de_results_distance[[paste(gen, ct, sep = "_")]] <- res_df
  }
}

# Optionally, save the entire de_results_distance object as an RDS file
saveRDS(de_results_distance, file = file.path(deseq_dist_folder, "de_results_distance_only.rds"))

print("Differential expression analysis (Distance-only) completed for each genotype and cell type.")
