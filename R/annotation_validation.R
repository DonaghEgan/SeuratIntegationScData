BiocManager::install("googleVis")
BiocManager::install("PCAtools")
library(googleVis)
library(dplyr)
library(tidyverse)
library(fgsea)
library(msigdbr)
library(pheatmap)

#' @param batch "string" input of the batch (melanoma/bcc)
#' @param orig.ident "string" input of the originally identity of the cluster being compared. 
#' @param compare list of length 2; contains new identity in string format of each cluster
#' @param counts.matrix counts matrix - should relate to counts from specified batch
#' @param cell.ident.matrix matrix containing batch, new annotation, and old annotation info 


compare_cluster <- function(batch, orig.ident, compare, counts.matrix, cell.ident.matrix){
  
  #filter cells to compare
  cluster_1 <- cell.ident.matrix %>% filter(batch== batch & cluster == orig.ident & annotation == compare[1])
  cluster_2 <- cell.ident.matrix %>% filter(batch== batch & cluster == orig.ident & annotation == compare[2])
  
  # combining data_frames
  cluster_combined <- rbind(cluster_1,cluster_2)
  cluster_combined <-  cluster_combined %>% remove_rownames %>% column_to_rownames(var="Cell_number")
  
  # subsetting counts matrix  
  counts <- counts.matrix[, which(colnames(counts.matrix) %in% rownames(cluster_combined))]
  
  # add data to seurat object for analyses 
  seurat.object <- CreateSeuratObject(counts, meta.data = cluster_combined, project = "comparing_clusters")
  
  #setting response as identity to perform DE
  Idents(seurat.object) <- "annotation"
  
  # Positive values indicate that the feature is more highly expressed in the first group
  markers <- FindMarkers(seurat.object, ident.1 = compare[1], ident.2 = compare[2],  
                         p.adjust.methods = 0.05)
  
  # pathway enrichment analysis
  named_genes_rank <- setNames(markers$avg_log2FC, as.character(rownames(markers)))
  h_gene_sets = msigdbr(species = "Homo sapiens", category = "H")
  msigdbr_list = split(x = h_gene_sets$gene_symbol, f = h_gene_sets$gs_name)
  fgseaRes <- fgsea(pathways = msigdbr_list, 
                    stats    = named_genes_rank,
                    minSize  = 30,
                    maxSize  = 300)
  
  return_list <- list(fgseaRes, markers)
  return(return_list)
  
}


test_compare <- compare_cluster("melanoma","Regulatory_T-cells", list("Chronically act/ex CD8 T cells", "Tregs"), counts_mel, cell.ident.matrix = cell_ident_matrix)


















