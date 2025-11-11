suppressPackageStartupMessages({
  library(Seurat); library(dplyr)
})
source("scripts/_utils_io.R")
cfg <- read_cfg(opt_config())

raw_dir <- cfg$paths$raw
min_cells <- cfg$qc$min_cells
min_features <- cfg$qc$min_features
mt_cutoff <- cfg$qc$percent_mt

sample_dirs <- list.dirs(raw_dir, recursive = FALSE, full.names = TRUE)
seurat_list <- list()

for (sample_path in sample_dirs) {
  sample_name <- basename(sample_path)
  counts <- Read10X(data.dir = sample_path)
  seu <- CreateSeuratObject(counts, project = sample_name,
                            min.cells = min_cells, min.features = min_features)
  seu$sample <- sample_name
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
  seu <- subset(seu, subset = percent.mt < mt_cutoff &
                          nFeature_RNA > min_features &
                          nFeature_RNA < cfg$qc$max_features)
  seu <- NormalizeData(seu, normalization.method = "LogNormalize",
                       scale.factor = cfg$norm$scale_factor)
  seu <- FindVariableFeatures(seu, selection.method = "vst",
                              nfeatures = cfg$norm$hvg)
  seurat_list[[sample_name]] <- seu
}
combined <- Reduce(function(a,b) merge(a,y=b), seurat_list)
saveRDS_safe(combined, "results/objects/01_combined.rds")
