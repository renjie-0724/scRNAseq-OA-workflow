suppressPackageStartupMessages({ library(Seurat) })
source("scripts/_utils_io.R")
cfg <- read_cfg(opt_config())

combined <- readRDS("results/objects/01_combined.rds")
combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = cfg$dimred$npcs, verbose = FALSE)
combined <- RunUMAP(combined, dims = 1:cfg$clustering$ndims)
combined <- FindNeighbors(combined, dims = 1:cfg$clustering$ndims)
combined <- FindClusters(combined, resolution = cfg$clustering$resolution)

pdf("results/figures/02_umap_by_cluster.pdf", width=6, height=5)
print(DimPlot(combined, group.by="seurat_clusters", label=TRUE))
dev.off()

saveRDS_safe(combined, "results/objects/02_clustered.rds")
