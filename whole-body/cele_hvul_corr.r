#!/use/bin/env Rscript 

# © EMBL-European Bioinformatics Institute, 2025
# Yuyao Song <ysong@ebi.ac.uk>
# Mar 2025

# this script is to get cross-species correlation heatmap

library(tidyverse)
library(Seurat)
library(pheatmap)
library(scGOclust)
library(stats)

cele_go_obj <- readRDS( "Cele_go_object.rds")

hydra_go_obj <- readRDS( "Hvul_go_object_whole_genome_ct_annots.rds")



cele_ct_go <- getCellTypeGO(go_seurat_obj = cele_go_obj, cell_type_col = 'cell_type_group')

hydra_ct_go <- getCellTypeGO(go_seurat_obj = hydra_go_obj, cell_type_col = 'cluster_no_num')

saveRDS(cele_ct_go, "cele_ct_go.rds")
saveRDS(hydra_ct_go, "hydra_ct_go.rds")

corr = crossSpeciesCellTypeGOCorr(species_1 = 'celegans', 
                                  species_2 = 'hvulgaris', 
                                  cell_type_go_sp1 = cele_ct_go, 
                                  cell_type_go_sp2 = hydra_ct_go, 
                                  corr_method = 'pearson')

saveRDS(corr, "cele_hvul_corr.rds")

plt = plotCellTypeCorrHeatmap(corr, width = 9, height = 10)

pdf("cele_hvul_corr.pdf")
print(plt)
dev.off()


# scale per column or row to see the relative similarity
plt2 = plotCellTypeCorrHeatmap(corr, scale = 'column', width = 9, height = 10)

pdf("cele_hvul_corr_colscaled.pdf")
print(plt2)
dev.off()
