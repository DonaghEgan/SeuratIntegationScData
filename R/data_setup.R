# Library
library(readr)
library(Seurat)

# Read in BCC metadata from TISH et al.,
bcc_meta <- read_tsv("/home/degan/SeuratIntegationScData/TISCH_BCC/BCC_GSE123813_aPD1_CellMetainfo_table.tsv")
melanoma_meta <- read_tsv("/home/degan/SeuratIntegationScData/TISCH_mel/SKCM_GSE120575_aPD1aCTLA4_CellMetainfo_table.tsv")

# read in seurat object from integration analysis
immune_combined <- readRDS("/home/degan/integration_2.0_march/saved_objects/immune_combined.Rds")

#read in cell identity matrix

# providing cells with annotations 
immune_combined <- RenameIdents(immune_combined, `0` = "Naive CD4 T cells", `1` = "Chronically act/ex CD8 T cells", `2` = "CD8 Memory T cells", 
                                `3` = "B cells_1", `4` = "CD8 active T cells", `5` = "Monocytes/Macrophages", `6` = "CD8 effector memory T cells", 
                                `7` = "Tregs", `8` = "Plasma Cells", 
                                `9` = "early activated T cells", `10` = "Tumor_1", `11` = "early activated CD8 T cells", 
                                `12` = "Proliferative T cells", 
                                `13` = "Fibroblasts/Melanocytes", `14` = "CD8 exhausted T cells", `15` = "pDCs", `16` = "Tumor_2", `17`= "NK cells" , 
                                `18` = "Mature DCs",`19` = "Endothelial cells", `20` = "B cells_2", `21` = "Intermediate monocyte")

# read in cell identity matrix
cell_ident_matrix <- readRDS("/home/degan/integration_2.0_march/saved_objects/integrated_annotated_2.0.Rds")

# loading in melanoma counts matrix 
counts_mel <- read.csv("/mnt/Data/Vadim/POIAZ/Marina/SCENICdata/whole_set/exprMat.csv", row.names = "X")
colnames(counts_mel) <- gsub("\\.", "-", colnames(counts_mel))


