# CosMx Analysis Repository

## Datasets
- **Seurat Objects**: Located in the Dropbox folder:
  - `FIRE-APP_data14424/Cosmx_FIREAPP/Results/Final Seurat Objects/`
  - Naming conventions: `C` for cortex, `H` for hippocampus

## Registration
- Landmarks manually assigned between IHC images and CosMx data.
- Transformation and registration performed using `transform_image_rois.ipynb`.

## Region Segmentation
- `createregionrois.py`: Python Dash application for interactive polygon drawing to define regions of interest (ROIs).
- `processregionrois.ipynb`: Processes manually defined ROIs and assigns them to Seurat object metadata.

## Plaque Region of Interest (ROI) Coordinates
- Plaque centroids manually identified and saved using FIJI ROI manager.
- Centroid coordinates exported as ROIs file.
- `cosmx_assignplaquexy.r`: R script that calculates the distance from each CosMx cell to the nearest plaque centroid and stores this information as metadata.

## Main Analysis Script
- `cosmx_analysis_1.r`: Performs the following tasks:
  - CosMx data pre-processing
  - Cell type identification
  - Genotype assignment

## Subclustering and Cell Type Contamination Removal
- `subclustering_contamination_script.ipynb`: Identifies and removes clusters contaminated with inappropriate transcript expression.

## Differential Expression Analysis
- `cosmx_deseq2.r`: R script performing differential expression analysis across genotypes.

## Plaque Distance and Genotype Analysis
- `cosmx_deseq2_plaquedistance_x_condition_effect.r`: DESeq2 analysis evaluating the interaction effect of plaque distance and genotype.
- `cosmx_deseq2_plaquedistanceeffect.r`: DESeq2 analysis evaluating the effect of plaque distance alone.

