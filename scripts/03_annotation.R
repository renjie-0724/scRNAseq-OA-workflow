suppressPackageStartupMessages({ library(Seurat); library(dplyr) })
source("scripts/_utils_io.R"); cfg <- read_cfg(opt_config())
seu <- readRDS("results/objects/02_clustered.rds")

# 简易标注示例：用你在 config 里给的 reference markers 进行打分/匹配
ref <- cfg$annotation$reference_markers   # list(celltype = c("GENE1","GENE2"))

markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
clusters <- unique(markers$cluster)
ann <- setNames(rep(NA_character_, length(clusters)), clusters)

for (cl in clusters) {
  top_genes <- head(markers |> dplyr::filter(cluster==cl) |> pull(gene), 40)
  scores <- sapply(ref, function(gs) length(intersect(toupper(gs), toupper(top_genes))))
  ann[as.character(cl)] <- names(which.max(scores))
}
seu <- RenameIdents(seu, ann); seu$cell_type <- Idents(seu)

pdf("results/figures/03_umap_by_celltype.pdf", width=6, height=5)
print(DimPlot(seu, label=TRUE) + NoLegend())
dev.off()

saveRDS_safe(seu, "results/objects/03_annotated.rds")
write.csv(ann, "results/tables/03_cluster_to_celltype.csv", quote=FALSE)
