# CosMx Analysis Repository

## Methods Writeup:

Datasets

The datasets comprised 20x magnification immunohistochemistry (IHC) images of amyloid plaques in TIFF format, alongside a CosMx-generated Seurat object containing single-cell RNA transcript counts (mouse RNA panel) and corresponding xy spatial coordinates for individual cells.

Registration

Manual landmark matching was performed to register the IHC images with the spatial coordinates from the CosMx dataset. Landmarks were selected using Fiji (for IHC images) and a Python Dash web application (for visualizing CosMx xy cell coordinates). Subsequently, an affine transformation using the Python skimage library was applied to align the IHC coordinates to the CosMx spatial reference. Transformed images and the corresponding transformation matrices were saved for subsequent analyses.

Region Segmentation

A Python Dash application was developed to visualize CosMx xy coordinates interactively, enabling manual polygon drawing to delineate anatomical regions. Polygons defining cortical and hippocampal regions were manually drawn. The xy coordinates of these polygons were then used to generate convex hulls, which defined region boundaries. Cells whose coordinates fell within these convex hulls were annotated accordingly in the metadata of the Seurat object.

Plaque Region of Interest (ROI) Coordinates

The original IHC images were loaded into Fiji, and plaque centroids were manually identified and selected using Fiji's ROI Manager. Centroid xy coordinates were exported to CSV format and subsequently transformed into CosMx dataset coordinates using previously generated transformation matrices. The Euclidean distance from each CosMx cell to its nearest plaque centroid was calculated, and these distances were recorded as metadata within the Seurat object, expressed in micrometers (µm).

CosMx Pre-processing

Datasets were split into hippocampal and cortical subsets for downstream analyses. Normalization and variable gene selection were performed using SCTransform within the Seurat framework. Principal Component Analysis (PCA) was conducted on normalized data, retaining 30 principal components (PCs). An elbow plot was generated to identify PCs for further analyses. Neighborhood graphs were constructed using the first 20 PCs, and clustering was performed at a resolution of 0.5. Uniform Manifold Approximation and Projection (UMAP) dimensionality reduction was executed, adjusting parameters to 30 neighbors and a minimum distance of 0.3 to optimize cluster separation. Cell clusters were visualized using UMAP, and feature plots for selected markers (Gad1, Pdgfra, Cx3cr1) were produced to assess cluster identity.

Cell Type Identification

The dataset was converted into a SingleCellExperiment (SCE) object for compatibility with downstream analyses. Reference mouse RNA sequencing data from celldex was utilized alongside SingleR for transcriptomic-based cell type prediction. Cell identities were further refined: inhibitory neurons were defined based on Gad1 expression, while oligodendrocyte precursor cells (OPCs) were identified by expression of Pdgfra. Final cell type labels included Excitatory Neurons, Inhibitory Neurons, Astrocytes, Oligodendrocytes, OPCs, Microglia, Macrophages, and Endothelial cells.

Genotypes

Genotype information was imported from external sources and integrated into the Seurat object's metadata to annotate each cell with corresponding genotype labels.

Subclustering and Cell Type Contamination Removal

Due to known transcript contamination issues in CosMx data, particularly in aged samples, subclustering analyses were performed for each major cell type. Subclusters were visualized using UMAP, and clusters were assessed through differential expression gene (DEG) analyses. Clusters displaying inappropriate gene expression patterns relative to their designated cell types, such as neuronal clusters expressing oligodendrocyte markers (e.g., Mbp, Mog), were excluded from further analyses.

Differential Expression Analysis

Differential gene expression analyses were conducted across four distinct genotypes (WT, FIRE, APP-WT, APP-FIRE) using DESeq2 in R. Pairwise genotype comparisons were performed independently, and resulting data, including CSV outputs and volcano plots, were saved for each genotype comparison.

Plaque Distance and Genotype Analysis

Differential expression analyses were conducted to assess transcriptomic changes relative to distance from amyloid plaques. Three analytical approaches were employed: (1) plaque distance alone (genotype effects ignored), (2) plaque distance while accounting for genotype differences (no interaction effect), and (3) the interaction between genotype and plaque distance to reveal condition-specific gene expression changes associated with plaque proximity. 


Datasets

- **Seurat Objects**: Located in the Dropbox folder:
  - `FIRE-APP_data14424/Cosmx_FIREAPP/Results/Final Seurat Objects/`
  - Naming conventions: `C` for cortex, `H` for hippocampus

Registration

- Landmarks manually assigned between IHC images and CosMx data.
- Transformation and registration performed using `transform_image_rois.ipynb`.

Region Segmentation

- `createregionrois.py`: Python Dash application for interactive polygon drawing to define regions of interest (ROIs).
- `processregionrois.ipynb`: Processes manually defined ROIs and assigns them to Seurat object metadata.

Plaque Region of Interest (ROI) Coordinates

- Plaque centroids manually identified and saved using FIJI ROI manager.
- Centroid coordinates exported as ROIs file.
- `cosmx_assignplaquexy.r`: R script that calculates the distance from each CosMx cell to the nearest plaque centroid and stores this information as metadata.

Main Analysis Script

- `cosmx_analysis_1.r`: Performs the following tasks:
  - CosMx data pre-processing
  - Cell type identification
  - Genotype assignment

Subclustering and Cell Type Contamination Removal

- `subclustering_contamination_script.ipynb`: Identifies and removes clusters contaminated with inappropriate transcript expression.

Differential Expression Analysis

- `cosmx_deseq2.r`: R script performing differential expression analysis across genotypes.

Plaque Distance and Genotype Analysis

- `cosmx_deseq2_plaquedistance_x_condition_effect.r`: DESeq2 analysis evaluating the interaction effect of plaque distance and genotype.
- `cosmx_deseq2_plaquedistanceeffect.r`: DESeq2 analysis evaluating the effect of plaque distance alone.

