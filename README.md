# selection-scaling-seurat
This package contains two functions a) `variable_features()` and b) `scaling_seurat()` which correspond to the third (3rd) and fourth (4th) steps of the standard scRNAseq processing workflow.

The `variable_features()` function calculated which genes have high cell-to-cell variation in the dataset and will be used to distinguish cell types. 

The `scaling_seurat()` function removes unwanted sources of variation to ensure that the succeeding clusterings (for dimensionality reduction) are due to actual differences in gene expression and not due to unwanted sources of variations.

## Installation
The package can be installed using
```
devtools::install_github("BTIP2024/selection-scaling-seurat")
```
## Example
The output of this function would be an rds file which when loaded in R, would be a Seurat object.
```
# for the selection of variable feature 
variable_features("normalized_seurat.rds")

# to scale after variable feature selection
scaling_seurat("variable_features.rds")
```
## scRNAseq processing workflow 

The standard scRNAseq processing workflow with the R package Seurat consists of seven (7) steps. The output of this package and function should be used as input for the scRNAseq processing pipeline. 

The following are the repositories of the packages for every step of the pipeline:
1. QC and filtering: [qualitycontrolseurat package](https://github.com/BTIP2024/quality-control-seurat)
2. Normalization: [qualitycontrolseurat package](https://github.com/BTIP2024/quality-control-seurat)
3. Identification of highly variable features: [selectionscalingseurat package](https://github.com/BTIP2024/selection-scaling-seurat)
4. Scaling: [selectionscalingseurat package](https://github.com/BTIP2024/selection-scaling-seurat)
5. Linear Dimensionality Reduction (PCA): [pcaseurat package](https://github.com/BTIP2024/pca-seurat)
6. Clustering: [nonlinearreduction package](https://github.com/BTIP2024/non-linear-reduction)
7. Non-linear dimensionality reduction (t-SNE and UMAP): [nonlinearreduction package](https://github.com/BTIP2024/non-linear-reduction)

An overview of the pipeline and its outputs can be observed below:
![](https://github.com/user-attachments/assets/3e49e900-6c6f-4124-98e7-9fde68c0d4c8)
