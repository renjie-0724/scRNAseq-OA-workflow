## 

suppressPackageStartupMessages({
  library(Seurat)
  library(GSVA)
  library(GSEABase)
  library(msigdbr)
  library(dplyr)
  library(pheatmap)
  library(yaml)
})

## 
seu <- readRDS("results/seurat_with_annotation.rds")
expr_sc <- GetAssayData(seu, slot = "data") |> as.matrix()

## 
cfg <- yaml::read_yaml("config/config.yaml")

collections <- cfg$gsea$collections     
min_size    <- cfg$gsea$min_size
max_size    <- cfg$gsea$max_size
method      <- cfg$gsea$score          
parallel    <- cfg$gsea$parallel
ncores      <- cfg$gsea$ncores

## 
get_msig_sets <- function(cat){
  msig <- msigdbr(species = "Homo sapiens", category = cat)

  if(cat == "C2" && length(cfg$gsea$c2_filter) > 0){
    msig <- msig |> filter(gs_subcat %in% cfg$gsea$c2_filter)
  }
  if(cat == "C5" && length(cfg$gsea$c5_filter) > 0){
    msig <- msig |> filter(ont %in% cfg$gsea$c5_filter)
  }

  split(msig$gene_symbol, msig$gs_name) |>
    lapply(unique)
}

## 
dir.create("results/gsva", showWarnings = FALSE, recursive = TRUE)

run_one_collection <- function(cat){

  message(sprintf("\n[INFO] Running %s ...", cat))

  genesets <- get_msig_sets(cat)

  param <- ssgseaParam(expr_sc,
                       genesets,
                       minSize = min_size,
                       maxSize = max_size,
                       normalize = TRUE)

  res <- gsva(param, parallel.sz = if(parallel) ncores else 1)

  ## pathway × cell matrix
  mat <- assay(res)

  ## 
  write.csv(mat,
            sprintf("results/gsva/gsva_%s_byCell.csv", cat),
            row.names = TRUE)
  saveRDS(res,
          sprintf("results/gsva/gsva_%s.rds", cat))

  ##
  clu <- Idents(seu)
  common_cells <- intersect(colnames(mat), names(clu))
  mat <- mat[, common_cells, drop = FALSE]
  clu <- droplevels(clu[common_cells])

  cluster_mat <- t(rowsum(t(mat), group = clu)) / as.vector(table(clu))

  write.csv(cluster_mat,
            sprintf("results/gsva/gsva_%s_byCluster.csv", cat),
            row.names = TRUE)

  ##
  pheatmap(cluster_mat,
           scale = "row",
           clustering_distance_rows = "correlation",
           clustering_distance_cols = "correlation",
           filename = sprintf("results/gsva/gsva_%s_clusterHeatmap.pdf", cat),
           width = 8, height = 10)

  message(sprintf("✅", cat))
}

##
for(cat in collections){
  run_one_collection(cat)
}

message("\n✅ All GSVA collections completed.\n")
