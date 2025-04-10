# Set working directory for Cortex Results
setwd("/Users/katherineridley/Projects/CosMx/APP/Cortex Results")

# Load required libraries
library(Seurat)
library(dplyr)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(crayon)  # For colored console output

### -----------------------------------------------------------
### DESeq2 Analysis with Interaction: Condition x Plaque_Distance per Cell Type
### -----------------------------------------------------------

# Define the cell types to analyze (should match exactly the labels in your metadata)
cell_types <- c("Excitatory Neurons", "Inhibitory Neurons", "Astrocytes", 
                "Oligodendrocytes", "OPCs", "Microglia & Macrophages", "Endothelial cells")

# Create folder for saving DESeq2 results and MA plots (interaction analysis)
deseq_folder <- "deseq2_results_per_celltype_interaction"
if (!dir.exists(deseq_folder)) {
  dir.create(deseq_folder)
}

# Initialize list to store DESeq2 results for each cell type
de_results <- list()

for (ct in cell_types) {
  message("Analyzing differential expression for cell type: ", ct)
  
  # Subset cells belonging to current cell type using WhichCells for reliability
  cells_to_use <- WhichCells(seurat_cortex, expression = Celltypes == ct)
  if (length(cells_to_use) == 0) {
    message("  Skipping ", ct, " because it has no cells.")
    next
  }
  subset_ct <- subset(seurat_cortex, cells = cells_to_use)
  
  # For DESeq2, include only cells from the two conditions with plaques: APP-WT and APP-FIRE.
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
  
  # Create DESeq2 dataset with design including interaction: ~ Condition * Plaque_Distance
  dds <- DESeqDataSetFromMatrix(countData = counts_mat,
                                colData = meta_dds,
                                design = ~ Condition * Plaque_Distance)
  dds <- estimateSizeFactors(dds, type = "poscounts")
  dds <- DESeq(dds)
  
  # List available coefficient names (for debugging)
  rn <- resultsNames(dds)
  message("Available coefficient names: ", paste(rn, collapse = ", "))
  
  # Extract results for the interaction term
  interaction_term <- "ConditionAPP.FIRE.Plaque_Distance"
  if (!(interaction_term %in% rn)) {
    stop("Interaction term '", interaction_term, "' not found in the DESeq2 results.")
  }
  res_int <- results(dds, name = interaction_term)
  res_int_df <- as.data.frame(res_int)
  res_int_df$gene <- rownames(res_int_df)
  res_int_df$cell_type <- ct
  res_int_df$comparison <- "Interaction: Condition x Plaque_Distance"
  
  # Save the DE results as CSV for this cell type
  out_csv <- file.path(deseq_folder, paste0("results_", gsub(" ", "_", ct), "_Interaction.csv"))
  write.csv(res_int_df, file = out_csv, row.names = FALSE)
  message("  Saved DESeq2 CSV for ", ct, " in ", out_csv)
  
  # Mark significance: padj < 0.05 and |log2FoldChange| >= 0.25
  res_int_df <- res_int_df %>% 
    mutate(Significant = ifelse(padj < 0.05 & abs(log2FoldChange) >= 0.25, "Yes", "No"))
  
  # Print the names of significant genes in red (if any); others are not printed
  sig_genes <- res_int_df %>% filter(Significant == "Yes") %>% arrange(padj)
  if(nrow(sig_genes) > 0) {
    message("Significant interaction genes for ", ct, ":")
    for(i in 1:nrow(sig_genes)) {
      message(red(sig_genes$gene[i]))
    }
  } else {
    message("No significant interaction genes found for ", ct)
  }
  
  # Create an MA plot for the interaction effect
  ma_plot_int <- ggplot(res_int_df, aes(x = baseMean, y = log2FoldChange)) +
    geom_point(aes(color = Significant), alpha = 0.6) +
    scale_x_log10() +
    scale_color_manual(values = c("No" = "grey", "Yes" = "red")) +
    ggtitle(paste("MA Plot:", ct, "- Interaction (Condition x Plaque_Distance)")) +
    xlab("Mean Expression (log scale)") +
    ylab("Log2 Fold Change per mm (Interaction)") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  ma_file <- file.path(deseq_folder, paste0("MAPlot_", gsub(" ", "_", ct), "_Interaction.png"))
  ggsave(filename = ma_file, plot = ma_plot_int, width = 8, height = 6)
  message("  Saved interaction MA plot for ", ct, " in ", ma_file)
  
  # Save result in our list
  de_results[[ct]] <- res_int_df
}

# Interpretation:
message("Interpretation: The interaction coefficient (on the y-axis of the MA plots) represents how much the effect of plaque distance on gene expression differs between APP-FIRE and APP-WT. For example, a negative value means that the gene’s expression decreases more with increasing plaque distance in APP-FIRE compared to APP-WT (i.e., higher expression close to plaques in APP-FIRE), whereas a positive value indicates the opposite.")

# Optionally, save the entire de_results object as an RDS file for later use
saveRDS(de_results, file = file.path(deseq_folder, "de_results_interaction_all_celltypes.rds"))

