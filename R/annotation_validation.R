BiocManager::install("googleVis")
library(googleVis)
library(dplyr)
library(tidyverse)

# sankey plots for our annotations and tisch annotations 
plot(getSankey(cell_ident_matrix$cluster, clusters = cell_ident_matrix$annotation))
plot(getSankey(bcc_meta$`Celltype (original)`, bcc_meta$`Celltype (minor-lineage)`))


exhausted_tregs_mel <- cell_ident_matrix %>% filter(batch=="melanoma" & cluster=="Regulatory_T-cells" & annotation == "Chronically act/ex CD8 T cells")
