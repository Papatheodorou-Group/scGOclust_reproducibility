#!/use/bin/env Rscript 

# © EMBL-European Bioinformatics Institute, 2025
# Yuyao Song <ysong@ebi.ac.uk>
# Mar 2025

# this script is to create the scGOclust object of C.elegans using eggnog-mapper anotated GOs and scRNA-seq data

library(tidyverse)
library(Seurat)
library(scGOclust)

# this is the eggnog-mapper annotation
cele_all = read_tsv("./celegans_WS260_GO_non_electronic_GO_only.tsv", col_names = FALSE)

gene_name_to_go = cele_all %>% separate_rows(X2, sep = ',') %>% filter(X2 != '-') %>% 
rename(external_gene_name = X1, go_id = X2)

# filter for BP processes
bp_terms = read_tsv("../go_ontoloy/biological_process_go_terms.txt", col_names = FALSE)

gene_name_to_go_bp = gene_name_to_go %>% filter(go_id %in% bp_terms$X1)  %>% mutate(go_category = 'biological_process')


gene_name_to_go_bp = merge(gene_name_to_go_bp, bp_terms, by.x = 'go_id', by.y = 'X1')
dim(gene_name_to_go_bp)

gene_name_to_go_bp %>% write_tsv("Cele_eggnog_name_to_go_id.tsv")

# map gene name, which was in the proteome, to WB id, which is in the scRNA-seq data
gene_name_to_wb_id = read_tsv("./c_elegans.PRJNA13758.WS260.protein.name_to_WBid.tsv", col_names = FALSE)

gene_name_and_id_to_go = merge(gene_name_to_go_bp, gene_name_to_wb_id, by.x = 'external_gene_name', by.y = 'X1')

gene_name_and_id_to_go = gene_name_and_id_to_go %>% 
rename(name_1006 = X2.x) %>% 
rename(WB_id = X2.y)

# 11542 peptide names only correspond to 7076 WB gene ids

gene_name_and_id_to_go %>% select(external_gene_name, WB_id) %>% unique() %>% filter(duplicated(WB_id) | duplicated(WB_id, fromLast = T) )

# well, thats alright, I do see lots of differnet peptides correspond to the same gene

head(gene_name_and_id_to_go, n=200)

gene_name_and_id_to_go %>% write_tsv("Cele_eggnog_name_to_WB_id_to_go_id.tsv")

# load seurat count object

obj <- readRDS("./cds_baseline_post_sub_seurat_with_fasta_name.rds")


seurat_obj <- obj
ensembl_to_GO <- gene_name_and_id_to_go
feature_type = 'WB_id'

  message("collect data")
  counts <- as.matrix(seurat_obj@assays$RNA@counts)

head(rownames(counts))

  if(!(feature_type %in% colnames(ensembl_to_GO))){
    stop(paste0(feature_type, " is not in colnames(ensembl_to_GO), please check the var type"))
  }

colnames(gene_name_and_id_to_go)

  ## pivot GO to feature type table to matrix
  go_matrix <- ensembl_to_GO %>%
    dplyr::mutate(placehold = 1) %>%
    tidyr::pivot_wider(id_cols = eval(feature_type), names_from = name_1006, values_from = placehold, values_fill = 0, values_fn = dplyr::first)

# this way of pivoting, if the gene WB id has a GO annot in any peptide, it will be 1

shared <- intersect(go_matrix[[feature_type]], rownames(counts))


# anyway, all genes with GO annotation are in the atlas, so we have to go with this

  go_matrix <- go_matrix %>%
    dplyr::filter(get(feature_type) %in% shared) %>%
    dplyr::arrange(get(feature_type))
  go_matrix <- go_matrix %>% tibble::column_to_rownames(eval(feature_type))


  counts <- t(counts)
  counts <- counts %>% as.data.frame()
  counts <- counts %>% dplyr::select(eval(shared))
  counts <- counts[, match(rownames(go_matrix), colnames(counts))]


  if (!(all(colnames(counts) == rownames(go_matrix)))) {
    stop("error during calculation, please check input format")
  }

  message("compute GO to cell matrix, might take a few secs")
  start_time <- Sys.time()
  go_mtx <- Matrix::Matrix(as.matrix(counts), sparse = TRUE) %*% Matrix::Matrix(as.matrix(go_matrix), sparse = TRUE)
  go_obj <- Seurat::CreateSeuratObject(counts = t(as.matrix(go_mtx)), meta.data = seurat_obj@meta.data)
  end_time <- Sys.time()

  message(paste0("time used: ", round(end_time - start_time, 2), " secs"))

  all_zero_terms = names(rowSums(as.matrix(go_obj@assays$RNA@counts))[which(rowSums(as.matrix(go_obj@assays$RNA@counts)) == 0)])
  message(paste0("removing ", length(all_zero_terms)), " all zero terms")

  go_obj <- go_obj[which(!(rownames(go_obj) %in% all_zero_terms)), ]

  message("returning GO Seurat object")

saveRDS(go_obj, "Cele_go_object.rds")