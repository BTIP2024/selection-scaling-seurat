#' Variable feature selection and scaling of seurat object
#' 
#' This package contains functions that identifies the variable features in a seurat object after quality filtering and normalization. The second function is scaling which should be done after variable feature identification.
#' 
#' @param input is the file after variable feature selection and scaling of the seurat object
#' @examples
#' variable_features("input.rds")
#' @export
#' 
variable_features <- function(input){
   normalized <- readRDS("input")
   normalized <- Seurat::FindVariableFeatures(normalized, selection.method = "vst", nfeatures = 2000)
   
   top10 <- head(Seurat::VariableFeatures(normalized), 10)
   
   plot1 <- VariableFeaturePlot(normalized)
   LabelPoints(plot = plot1, points = top10, repel = TRUE)
   
   ggplot2::ggsave(filename = variable_features_plot.png)
   saveRDS(top10, file = "top10_variablefeatures.rds")
   saveRDS(normalized, file = "variablefeatures.rds")
}

scaling_seurat <- function(input){
   for_scaling <- readRDS(input)
   all.genes <- rownames(for_scaling)
   scaling <- Seurat::ScaleData(input, features = all.genes)
}