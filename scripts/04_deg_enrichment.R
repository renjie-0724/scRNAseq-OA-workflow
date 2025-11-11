suppressPackageStartupMessages({
  library(Seurat); library(clusterProfiler); library(org.Hs.eg.db); library(dplyr)
})
source("scripts/_utils_io.R"); cfg <- read_cfg(opt_config())
seu <- readRDS("results/objects/03_annotated.rds"); Idents(seu) <- "cell_type"

contrast <- cfg$de$contrast  # list(ident.1="Lining_fibroblasts", ident.2="Sublining_fibroblasts")
deg <- FindMarkers(seu, ident.1 = contrast$ident.1, ident.2 = contrast$ident.2,
                   logfc.threshold = 0, min.pct = 0)
write.csv(deg, "results/tables/04_deg.csv")

sig <- deg %>% filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25)
genes <- rownames(sig)
gene_df <- bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)

ego <- enrichGO(gene_df$ENTREZID, OrgDb=org.Hs.eg.db, ont="BP", readable=TRUE)
ekg <- enrichKEGG(gene_df$ENTREZID, organism="hsa")

openxlsx::write.xlsx(list(GO=as.data.frame(ego), KEGG=as.data.frame(ekg)),
                     file="results/tables/04_enrichment.xlsx")
