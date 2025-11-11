suppressPackageStartupMessages({ library(Seurat); library(GSVA); library(msigdbr); library(ggplot2) })
source("scripts/_utils_io.R"); cfg <- read_cfg(opt_config())
seu <- readRDS("results/objects/03_annotated.rds")
expr <- as.matrix(GetAssayData(seu, slot="data"))

msig_h <- msigdbr(species="Homo sapiens", category="H")
gene_sets <- split(msig_h$gene_symbol, msig_h$gs_name)

gsva_mat <- gsva(expr, gene_sets, method="ssgsea", kcdf="Gaussian")
saveRDS_safe(gsva_mat, "results/objects/06_gsva_ssgsea.rds")
