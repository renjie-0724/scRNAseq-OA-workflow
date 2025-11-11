suppressPackageStartupMessages({ library(Seurat); library(CellChat) })
source("scripts/_utils_io.R"); cfg <- read_cfg(opt_config())
seu <- readRDS("results/objects/03_annotated.rds")
Idents(seu) <- "cell_type"

data.input <- GetAssayData(seu, slot="data")
meta <- seu@meta.data
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "cell_type")
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

dir.create("results/cellchat", showWarnings = FALSE, recursive = TRUE)
pdf("results/cellchat/05_net_circle.pdf", width=7, height=7)
netVisual_circle(cellchat@net$count, vertex.weight=as.numeric(table(cellchat@idents)), weight.scale=TRUE)
dev.off()

saveRDS_safe(cellchat, "results/objects/05_cellchat.rds")
