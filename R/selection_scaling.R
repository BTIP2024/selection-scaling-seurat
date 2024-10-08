#' Variable feature selection and scaling of seurat object
#' 
#' This package contains functions that identifies the variable features in a seurat object after quality filtering and normalization. The second function is scaling which should be done after variable feature identification.
#' 
#' @param input is the file after the normalization of the seurat object
#' @examples
#' variable_features("input.rds")
#' @export
variable_features <- function(input){
   if(!(tools::file_ext(input)) == "rds") {
      return("Input file should be an rds file")
   } else if(tools::file_ext(input) == "rds") {
   normalized <- readRDS(input)
   
   if(class(normalized) != "Seurat") {
      return("File is not a seurat object")
   } else {
      normalized <- Seurat::FindVariableFeatures(normalized, selection.method = "vst", nfeatures = 2000)
      
      top10 <- head(Seurat::VariableFeatures(normalized), 10)
      
      plot1 <- Seurat::VariableFeaturePlot(normalized)
      plot2 <- Seurat::LabelPoints(plot = plot1, points = top10, repel = TRUE)
      
      finalplots <- plot1 + plot2
      
      ggplot2::ggsave(finalplots, filename = "variable_features_plot.png", width = 15, height = 8, dpi = 300)
      saveRDS(top10, file = "top10_variablefeatures.rds")
      saveRDS(normalized, file = "variablefeatures.rds")
   }
   }
}

scaling_seurat <- function(input){
   if(!(tools::file_ext(input)) == "rds") {
      return("Input file should be an rds file")
   } else if(tools::file_ext(input) == "rds") {
      for_scaling <- readRDS(input)
      
      if(class(for_scaling) != "Seurat") {
         return("File is not a seurat object")
      } else {
   all.genes <- rownames(for_scaling)
   for_scaling <- Seurat::ScaleData(for_scaling, features = all.genes)
   saveRDS(for_scaling, file = "after_scaling.rds")
      }
   }
   }
