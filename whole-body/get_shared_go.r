#!/use/bin/env Rscript 

# © EMBL-European Bioinformatics Institute, 2025
# Yuyao Song <ysong@ebi.ac.uk>
# Mar 2025

# this script is to get shared co-upregulated GO terms between species

library(tidyverse)
library(Seurat)
library(scGOclust)
library(stats)
library(pheatmap)

cele_go_analyzed = readRDS("../cele_go_analyzed.rds")
hvul_go_analyzed = readRDS("../hydra_go_analyzed.rds")

shared_go = getCellTypeSharedGO(species_1 = 'celegans', species_2 = 'hvulgaris', 
                    cell_type_col_sp1 = 'cell_type_group', 
                    cell_type_col_sp2 = 'cluster_no_num',
                    analyzed_go_seurat_sp1 = cele_go_analyzed, 
                    analyzed_go_seurat_sp2 = hvul_go_analyzed,
                                slot_use = 'data'
                       )


saveRDS(shared_go, "cele_hvul_shared_go.rds")
shared_go$shared_sig_markers %>% write_csv("cele_hvul_shared_go_signif.csv")