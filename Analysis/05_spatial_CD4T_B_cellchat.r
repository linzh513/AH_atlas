# Spatial transcriptomics: CD4 T-cell and B-cell communication (CellChat)

# Environment and packages
.libPaths("/home/lzh/miniconda3/envs/singlecell/lib/R/library")
library(reticulate)
use_python("/home/lzh/miniconda3/envs/cellchat_py/bin/python", required = TRUE)
py_config()
py_module_available("umap")
py_module_available("numpy")
py_module_available("pandas")
py_module_available("sklearn")
set.seed(1234)
setwd("/mnt/data_nvme1/lzh/test/AH/atlas/CCC/cellchat_sc_tb")
library(CellChat)
packageVersion("CellChat")
library(patchwork)
library(Seurat)
library(ggplot2)
library(svglite)
library(circlize)
library(dplyr)
library(ggalluvial)
library(ComplexHeatmap)
library(forcats)
library(data.table)
options(stringsAsFactors = FALSE)

# 01. Load annotated CD4 T-cell and B-cell objects
t.all <- readRDS("/mnt/data_nvme1/lzh/test/AH/atlas/recluster/CD4T/CD4T.all.RDS")
b.all <- readRDS("/mnt/data_nvme1/lzh/test/AH/atlas/recluster/B/B.all.RDS")
table(t.all$celltype_t)
table(t.all$TH_meta)
table(b.all$celltype_b)
table(b.all$TH_meta)

# 02. Harmonize cell-type labels and merge objects
t.all$celltype_ccc <- paste0("CD4_T_", t.all$celltype_t)
b.all$celltype_ccc <-  b.all$celltype_b
spatial.all <- merge(t.all, b.all)
spatial.all$Group <- ifelse(spatial.all$TH_meta == "N", "AH","ATH")
table(spatial.all$celltype_ccc)
# Numeric group label derived from TH_meta.

# Align spatial coordinates
patient_col <- "Patient"
get_spatial_df <- function(obj, obj_name, patient_col){
  emb <- Embeddings(obj, reduction = "spatial")
  df <- as.data.frame(emb)
  colnames(df)[1:2] <- c("spatialct_1","spatialct_2")
  df$cell <- rownames(emb)
  df$orig.ident <- obj_name
  df$Patient <- obj@meta.data[df$cell, patient_col, drop=TRUE]
  df
}
t_spatial_df <- get_spatial_df(t.all, "t.all", patient_col)
b_spatial_df <- get_spatial_df(b.all, "b.all", patient_col)
spatial_all_df <- bind_rows(t_spatial_df, b_spatial_df)
print(head(spatial_all_df))
print(table(spatial_all_df$orig.ident))
print(summary(spatial_all_df[, c("spatialct_1","spatialct_2")]))

# Restore spatial coordinates in the merged Seurat object
library(Seurat)
stopifnot(all(c("cell","spatialct_1","spatialct_2") %in% colnames(spatial_all_df)))
cells <- Cells(spatial.all)
df_aligned <- spatial_all_df[match(cells, spatial_all_df$cell), ]
if(any(is.na(df_aligned$spatialct_1))) {
  bad <- sum(is.na(df_aligned$spatialct_1))
  stop(paste0("匹配失败：有 ", bad, " 个 cell 的空间坐标没对上。请检查 cell 列是否一致。"))
}
emb_new <- as.matrix(df_aligned[, c("spatialct_1","spatialct_2")])
rownames(emb_new) <- cells
colnames(emb_new) <- c("spatialct_1","spatialct_2")
spatial.all[["spatial"]] <- CreateDimReducObject(
  embeddings = emb_new,
  key = "spatial_",
  assay = DefaultAssay(spatial.all)
)
print(names(spatial.all@reductions))
print(head(Embeddings(spatial.all, "spatial")))
saveRDS(spatial.all, "/mnt/data_nvme1/lzh/test/AH/atlas/CCC/cellchat_sc_tb/spatial.tb.RDS")

# 03. Split AH and ATH datasets
obj.list <- SplitObject(spatial.all, split.by = "Group")
print(obj.list)
table(obj.list[["AH"]]$celltype_ccc)
table(obj.list[["ATH"]]$celltype_ccc)
celltype_order <- c("CD4_T_Central_memory","CD4_T_Memory","CD4_T_Naive","CD4_T_Proliferative","CD4_T_Tfh","CD4_T_Treg",
                    "GCB_DZ","GCB_LZ","Memory_B","Naive_B","Plasma_cell","Pre_GCB")

# Prepare expression matrices and group metadata
data.input_AH <- GetAssayData(obj.list[["AH"]], slot = 'data')  # Log-normalized expression matrix.
data.input_ATH <- GetAssayData(obj.list[["ATH"]], slot = 'data')  # Log-normalized expression matrix.
meta_AH <- obj.list[["AH"]]@meta.data[,c("Group","celltype_ccc")]
meta_ATH <- obj.list[["ATH"]]@meta.data[,c("Group","celltype_ccc")]
identical(rownames(meta_AH), colnames(data.input_AH))
meta_AH$celltype <- factor(meta_AH$celltype_ccc ,levels = celltype_order)
ordered_indices <- order(meta_AH$celltype)
meta_AH <- meta_AH[ordered_indices, ]
data.input_AH <- data.input_AH[, ordered_indices]
identical(rownames(meta_AH),colnames(data.input_AH))
identical(rownames(meta_ATH), colnames(data.input_ATH))
meta_ATH$celltype <- factor(meta_ATH$celltype_ccc ,levels = celltype_order)
ordered_indices <- order(meta_ATH$celltype)
meta_ATH <- meta_ATH[ordered_indices, ]
data.input_ATH <- data.input_ATH[, ordered_indices]
identical(rownames(meta_ATH),colnames(data.input_ATH))

# 04. Infer CellChat networks for each group
CellChatDB <- CellChatDB.human  # Human ligand-receptor database.
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- CellChatDB

# AH communication model
cellchat_AH <- createCellChat(object = data.input_AH, meta = meta_AH, group.by = "celltype")
cellchat_AH <- addMeta(cellchat_AH, meta = meta_AH)
cellchat_AH <- setIdent(cellchat_AH, ident.use = "celltype")
levels(cellchat_AH@idents)
groupSize_AH <- as.numeric(table(cellchat_AH@idents))
cellchat_AH@DB <- CellChatDB.use
cellchat_AH <- subsetData(cellchat_AH)
cellchat_AH <- identifyOverExpressedGenes(cellchat_AH)
cellchat_AH <- identifyOverExpressedInteractions(cellchat_AH)
cellchat_AH <- computeCommunProb(cellchat_AH, type = "triMean")
cellchat_AH <- filterCommunication(cellchat_AH, min.cells = 10)
cellchat_AH <- computeCommunProbPathway(cellchat_AH)
cellchat_AH <- aggregateNet(cellchat_AH)
cellchat_AH <- netAnalysis_computeCentrality(cellchat_AH, slot.name = "netP")

# ATH communication model
cellchat_ATH <- createCellChat(object = data.input_ATH, meta = meta_ATH, group.by = "celltype")
cellchat_ATH <- addMeta(cellchat_ATH, meta = meta_ATH)
cellchat_ATH <- setIdent(cellchat_ATH, ident.use = "celltype")
levels(cellchat_ATH@idents)
groupSize_ATH <- as.numeric(table(cellchat_ATH@idents))
cellchat_ATH@DB <- CellChatDB.use
cellchat_ATH <- subsetData(cellchat_ATH)
cellchat_ATH <- identifyOverExpressedGenes(cellchat_ATH)
cellchat_ATH <- identifyOverExpressedInteractions(cellchat_ATH)
cellchat_ATH <- computeCommunProb(cellchat_ATH, type = "triMean")
cellchat_ATH <- filterCommunication(cellchat_ATH, min.cells = 10)
cellchat_ATH <- computeCommunProbPathway(cellchat_ATH)
cellchat_ATH <- aggregateNet(cellchat_ATH)
cellchat_ATH <- netAnalysis_computeCentrality(cellchat_ATH, slot.name = "netP")
saveRDS(cellchat_AH, "cellchat_AH.RDS")
saveRDS(cellchat_ATH, "cellchat_ATH.RDS")

# 05. Merge group-specific CellChat networks
object.list <- list(AH = cellchat_AH, ATH = cellchat_ATH)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
cellchat
outdir <- "CellChat_comparison"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Compare overall interaction numbers and strengths
gg1 <- compareInteractions(cellchat, show.legend = FALSE, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = FALSE, group = c(1,2), measure = "weight")
pdf(file.path(outdir, "01_compareInteractions_barplot.pdf"), width = 8, height = 5)
print(gg1 + gg2)
dev.off()

pdf(file.path(outdir, "02_diffInteraction_circle.pdf"), width = 8, height = 10)
netVisual_diffInteraction(cellchat, weight.scale = TRUE)
netVisual_diffInteraction(cellchat, weight.scale = TRUE, measure = "weight")
dev.off()
gg1 <- netVisual_heatmap(cellchat)
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
pdf(file.path(outdir, "03_diffInteraction_heatmap.pdf"), width = 10, height = 5)
print(gg1 + gg2)
dev.off()
weight.max.count <- getMaxWeight(object.list, attribute = c("idents", "count"))
pdf(file.path(outdir, "04_eachDataset_circle_count.pdf"), width = 13, height = 7)
par(mfrow = c(1, 2), xpd = TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(
    object.list[[i]]@net$count,
    weight.scale = TRUE,
    label.edge = FALSE,
    edge.weight.max = weight.max.count[2],
    edge.width.max = 12,
    title.name = paste0("Number of interactions - ", names(object.list)[i])
  )
}
dev.off()
weight.max.weight <- getMaxWeight(object.list, attribute = c("idents", "weight"))
pdf(file.path(outdir, "05_eachDataset_circle_weight.pdf"), width = 13, height = 7)
par(mfrow = c(1, 2), xpd = TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(
    object.list[[i]]@net$weight,
    weight.scale = TRUE,
    label.edge = FALSE,
    edge.weight.max = weight.max.weight[2],
    edge.width.max = 12,
    title.name = paste0("Interaction weights - ", names(object.list)[i])
  )
}
dev.off()

# 06. Compare outgoing and incoming signaling roles
num.link <- sapply(object.list, function(x) {
  rowSums(x@net$count) + colSums(x@net$count) - diag(x@net$count)
})
weight.MinMax <- c(min(num.link), max(num.link))
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(
    object.list[[i]],
    title = names(object.list)[i],
    weight.MinMax = weight.MinMax)
}
pdf(file.path(outdir, "06_signalingRole_scatter.pdf"), width = 10, height = 5)
print(patchwork::wrap_plots(plots = gg))
dev.off()

# 07. Compare pathway information flow
gg1 <- rankNet(cellchat, mode = "comparison", stacked = TRUE, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = FALSE, do.stat = TRUE)
pdf(file.path(outdir, "10_rankNet_comparison.pdf"), width = 12, height = 11)
print(gg1 + gg2)
dev.off()
i <- 1
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)

# Outgoing pathway roles
ht1 <- netAnalysis_signalingRole_heatmap(
  object.list[[i]],
  pattern = "outgoing",
  signaling = pathway.union,
  title = names(object.list)[i],
  width = 12,
  height = 24,
  color.heatmap = "Reds")
ht2 <- netAnalysis_signalingRole_heatmap(
  object.list[[i+1]],
  pattern = "outgoing",
  signaling = pathway.union,
  title = names(object.list)[i+1],
  width = 12,
  height = 24,
  color.heatmap = "Reds")
pdf(file.path(outdir, "11_signalingRole_heatmap_outgoing.pdf"), width = 12, height = 12)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

# Incoming pathway roles
ht1 <- netAnalysis_signalingRole_heatmap(
  object.list[[i]],
  pattern = "incoming",
  signaling = pathway.union,
  title = names(object.list)[i],
  width = 12,
  height = 24,
  color.heatmap = "Reds")
ht2 <- netAnalysis_signalingRole_heatmap(
  object.list[[i+1]],
  pattern = "incoming",
  signaling = pathway.union,
  title = names(object.list)[i+1],
  width = 12,
  height = 24,
  color.heatmap = "Reds")
pdf(file.path(outdir, "12_signalingRole_heatmap_incoming.pdf"), width = 12, height = 12)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

# 08. Compare network similarity

# Functional similarity
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")
pdf(file.path(outdir, "14_embedding_functional.pdf"), width = 8, height = 7.5)
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
dev.off()

pdf(file.path(outdir, "15_embedding_functional_zoomIn.pdf"), width = 12, height = 10)
netVisual_embeddingZoomIn(cellchat, type = "functional", nCol = 2)
dev.off()

# Structural similarity
cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
cellchat <- netClustering(cellchat, type = "structural")
pdf(file.path(outdir, "16_embedding_structural.pdf"), width = 8, height = 7.5)
netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5)
dev.off()

pdf(file.path(outdir, "17_embedding_structural_zoomIn.pdf"), width = 12, height = 10)
netVisual_embeddingPairwiseZoomIn(cellchat, type = "structural", nCol = 2)
dev.off()

pdf(file.path(outdir, "18_rankSimilarity_functional.pdf"), width = 5, height = 7.5)
rankSimilarity(cellchat, type = "functional")
dev.off()

pdf(file.path(outdir, "19_rankSimilarity_structural.pdf"), width = 5, height = 7.5)
rankSimilarity(cellchat, type = "structural")
dev.off()

# 09. Visualize CD45 signaling
sort(cellchat@netP$pathways)
pathways.show <- c("CD45")

# Pathway circle plots
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show)
pdf(file.path(outdir, "20_pathway_circle_CD45.pdf"), width = 10, height = 5)
par(mfrow = c(1, 2), xpd = TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(
    object.list[[i]], signaling = pathways.show,
    layout = "circle",
    edge.weight.max = weight.max[1], edge.width.max = 10,
    signaling.name = paste(pathways.show, names(object.list)[i])
  )
}
dev.off()

# Pathway heatmaps
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(
    object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",
    title.name = paste(pathways.show, "signaling", names(object.list)[i])
  )
}
pdf(file.path(outdir, "21_pathway_heatmap_CD45.pdf"), width = 12, height = 6)
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))
dev.off()

# Pathway chord diagrams
pdf(file.path(outdir, "22_pathway_chord_CD45.pdf"), width = 14, height = 7)
par(mfrow = c(1, 2), xpd = TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]],
    signaling = pathways.show, layout = "chord",
    signaling.name = paste(pathways.show, names(object.list)[i])
  )
}
dev.off()

# 10. Compare T-to-B ligand-receptor interactions
levels(object.list[[1]]@idents)
source.idx <- which(levels(object.list[[1]]@idents) %in% c(
  "CD4_T_Central_memory", "CD4_T_Memory", "CD4_T_Naive", "CD4_T_Proliferative", "CD4_T_Tfh", "CD4_T_Treg"
))
target.idx <- which(levels(object.list[[1]]@idents) %in% c(
  "GCB_DZ","GCB_LZ", "Memory_B", "Naive_B", "Plasma_cell", "Pre_GCB"
))

# Ligand-receptor bubble plot
pdf(file.path(outdir, "23_bubble_T_to_selectedTargets.pdf"), width = 13, height = 4.5)
netVisual_bubble(cellchat,
  sources.use = source.idx, targets.use = target.idx,
  comparison = c(1, 2), angle.x = 45
)
dev.off()

# Increased and decreased interactions
gg1 <- netVisual_bubble(cellchat,
  sources.use = source.idx, targets.use = target.idx,
  comparison = c(1, 2), max.dataset = 2,
  title.name = paste0("Increased signaling in ", names(object.list)[2]),
  angle.x = 90,
  remove.isolate = TRUE)
gg2 <- netVisual_bubble( cellchat,
  sources.use = source.idx, targets.use = target.idx,
  comparison = c(1, 2), max.dataset = 1,
  title.name = paste0("Decreased signaling in ", names(object.list)[2]),
  angle.x = 90,
  remove.isolate = TRUE)
pdf(file.path(outdir, "24_bubble_up_down_T_to_selectedTargets.pdf"), width = 22, height = 6)
print(gg1 + gg2)
dev.off()

# 11. Identify differential ligand-receptor interactions
table(cellchat@meta$datasets)
cellchat@meta$datasets <- factor(cellchat@meta$datasets, levels = c("AH", "ATH"))

# Differential expression analysis
pos.dataset <- "ATH"
features.name <- pos.dataset
# ATH-relative positive direction.
cellchat <- identifyOverExpressedGenes(cellchat,
  group.dataset = "datasets", pos.dataset = pos.dataset,
  features.name = features.name,only.pos = FALSE,
  thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 0.05)
net <- netMappingDEG(cellchat, features.name = features.name)
net.up <- subsetCommunication(cellchat, net = net,datasets = "ATH", ligand.logFC = 0.2, receptor.logFC = NULL)
net.down <- subsetCommunication(cellchat, net = net, datasets = "AH", ligand.logFC = -0.1, receptor.logFC = NULL)
net.de1 <- subsetCommunication(cellchat, net = net, ligand.logFC = 0.2, receptor.logFC = NULL)
net.de2 <- subsetCommunication(cellchat, net = net, ligand.logFC = -0.2, receptor.logFC = NULL)
net.de <- rbind(net.de1, net.de2)
write.csv(net.up, file.path(outdir, "25_net_up_in_ATH.csv"), row.names = FALSE)
write.csv(net.down, file.path(outdir, "26_net_down_in_ATH.csv"), row.names = FALSE)
write.csv(net.de, file.path(outdir, "27_net_de.csv"), row.names = FALSE)
gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)
gene.de <- extractGeneSubsetFromPair(net.de, cellchat)
write.csv(data.frame(gene = gene.up), file.path(outdir, "25_gene_up_in_ATH.csv"), row.names = FALSE)
write.csv(data.frame(gene = gene.down), file.path(outdir, "26_gene_down_in_ATH.csv"), row.names = FALSE)
write.csv(data.frame(gene = gene.de), file.path(outdir, "27_gene_de.csv"), row.names = FALSE)

# Differential ligand-receptor bubble plots
pairLR.use.up <- net.up[, "interaction_name", drop = FALSE]
pairLR.use.down <- net.down[, "interaction_name", drop = FALSE]
pairLR.use.de <- net.de[, "interaction_name", drop = FALSE]
levels(object.list[[1]]@idents)
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up,
  sources.use = c(1:6), targets.use = c(7:12),
  comparison = c(1, 2),angle.x = 90, remove.isolate = TRUE,
  title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down,
  sources.use = c(1:6), targets.use = c(7:12),
  comparison = c(1, 2), angle.x = 90, remove.isolate = TRUE,
  title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
gg3 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.de,
  sources.use = c(1:6), targets.use = c(7:12),
  comparison = c(1, 2), angle.x = 90, remove.isolate = TRUE,
  title.name = "Differential-regulated signaling")
pdf(file.path(outdir, "29_bubble_DEG_based_up_down.pdf"), width = 24, height = 8)
gg1
gg2
gg3
dev.off()

# Differential ligand-receptor chord diagrams
pdf(file.path(outdir, "30_chord_DEG_based_up_down.pdf"), width = 12, height = 10)
netVisual_chord_gene(object.list[[2]],
  sources.use = c(1:6), targets.use = c(7:12),
  slot.name = "net", net = net.up,
  lab.cex = 0.8, small.gap = 3.5,
  title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
netVisual_chord_gene(object.list[[1]],
  sources.use = c(1:6), targets.use = c(7:12),
  slot.name = "net", net = net.down,
  lab.cex = 0.8, small.gap = 3.5,
  title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
dev.off()

# 12. Compare signaling gene expression
pdf(file.path(outdir, "31_geneExpression_CD45.pdf"), width = 8, height = 6)
plotGeneExpression(cellchat, signaling = "CD45", split.by = "datasets", colors.ggplot = TRUE)
dev.off()

pdf(file.path(outdir, "32_geneExpression_CD6.pdf"), width = 8, height = 8)
plotGeneExpression(cellchat, signaling = "CD6", split.by = "datasets", colors.ggplot = TRUE)
dev.off()

pdf(file.path(outdir, "33_geneExpression_CD23.pdf"), width = 8, height = 12)
plotGeneExpression(cellchat, signaling = "CD23", split.by = "datasets", colors.ggplot = TRUE)
dev.off()

pdf(file.path(outdir, "34_geneExpression_SIRP.pdf"), width = 8, height = 12)
plotGeneExpression(cellchat, signaling = "SIRP", split.by = "datasets", colors.ggplot = TRUE)
dev.off()

pdf(file.path(outdir, "36_geneExpression_MHC_I.pdf"), width = 8, height = 12)
plotGeneExpression(cellchat, signaling = "MHC-I", split.by = "datasets", colors.ggplot = TRUE)
dev.off()

# 13. Save merged CellChat results
saveRDS(cellchat, file = file.path(outdir, "cellchat_merged_comparison_AH_vs_ATH.RDS"))
saveRDS(object.list, file = file.path(outdir, "object.list_AH.RDS"))

# 14. Plot spatial cell-type distributions
pal <- c("#d387b5","#BC9DCC", "#f8e6e0", "#f1c6dd","#7FB2D5","#e2edf9", "#CB5C5B","#F5B47F","#8DBA8F", "#e7efc8","#FFD569","#8DD2C5")
theme_global <- theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank()) + theme(aspect.ratio = 1)
pdf("CD4T_B_spatial_byPatient_celltype_ccc.pdf", width = 8, height = 6)
for (p in unique(spatial.all$Patient)) {
  message("Plotting Patient: ", p)
  obj.p <- subset(spatial.all, subset = Patient == p)
  p2 <- DimPlot(obj.p, reduction = "spatial", group.by = "celltype_ccc", pt.size = 1, cols = pal) + ggtitle(p) + theme_global
  print(p2)
}
dev.off()

# 15. Plot spatial distributions of signaling genes
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
lr_gene <- c("CD6", "ALCAM", "PTPRC", "CD22", "SIRPG","CD47")
DefaultAssay(spatial.all) <- "Spatial"
pdf("CD4T_B_spatial_byPatient_signal_gene.pdf", width = 8, height = 6)
for (p in unique(spatial.all$Patient)) {
  message("Plotting Patient: ", p)
  for (g in lr_gene) {
    message("Plotting Gene: ", g)
    obj.p <- subset(spatial.all, subset = Patient == p)
    p2 <- FeaturePlot(obj.p, reduction = "spatial", features = g, pt.size = 1) + ggtitle(paste0(p, "_", g)) + theme_global + scale_color_gradientn(colors = rev(brewer.pal(10, "RdYlBu")))
    print(p2)
  }
}
dev.off()



# 16. Donor-level expression and positivity of selected ligand-receptor genes
library(Seurat)
library(Matrix)
library(dplyr)
library(tidyr)
library(ggplot2)

obj <- spatial.all
assay_use <- "Spatial"
slot_use <- "data"
positive_threshold <- 0
AH_samples <- c("AH7", "AH8", "AH35")
ATH_samples <- c("AH5", "AH15", "AH17")
group_colors <- c("AH" = "#9EC4BE", "ATH-A" = "#FAC795")
targets <- data.frame(
  Interaction = rep(c("PTPRC-CD22", "CD6-ALCAM", "SIRPG-CD47"), each = 2),
  Gene = c("PTPRC", "CD22", "CD6", "ALCAM", "SIRPG", "CD47"),
  CellType = c("CD4_T_Tfh", "Naive_B", "CD4_T_Central_memory",
               "GCB_LZ", "CD4_T_Proliferative", "GCB_LZ"),
  stringsAsFactors = FALSE
)
targets$Pair <- paste0(targets$Gene, " in\n", targets$CellType)

# Extract only the selected genes and convert to a regular matrix.
DefaultAssay(obj) <- assay_use
missing_genes <- setdiff(targets$Gene, rownames(obj[[assay_use]]))
if (length(missing_genes)) {
  stop("Missing genes: ", paste(missing_genes, collapse = ", "))
}
if (!all(c("Patient", "celltype_ccc") %in% colnames(obj@meta.data))) {
  stop("Patient and celltype_ccc are required in the metadata.")
}
if (length(intersect(AH_samples, ATH_samples))) stop("Donor groups overlap.")
expr_mat <- as.matrix(GetAssayData(obj, assay = assay_use, slot = slot_use)[
  unique(targets$Gene), , drop = FALSE
])
stopifnot(all(colnames(expr_mat) %in% rownames(obj@meta.data)))
md <- obj@meta.data[colnames(expr_mat), , drop = FALSE]
stopifnot(identical(rownames(md), colnames(expr_mat)))
md$Patient <- as.character(md$Patient)
md$CellType <- as.character(md$celltype_ccc)
md$Group <- case_when(
  md$Patient %in% AH_samples ~ "AH",
  md$Patient %in% ATH_samples ~ "ATH-A",
  TRUE ~ NA_character_
)
if ("TH_meta" %in% names(md)) {
  expected <- case_when(
    as.character(md$TH_meta) == "N" ~ "AH",
    as.character(md$TH_meta) == "Y" ~ "ATH-A",
    TRUE ~ NA_character_
  )
  bad_group <- !is.na(md$Group) & !is.na(expected) & md$Group != expected
  if (any(bad_group)) stop("The donor lists disagree with TH_meta.")
}
donors <- data.frame(
  Patient = c(AH_samples, ATH_samples),
  Group = rep(c("AH", "ATH-A"), c(length(AH_samples), length(ATH_samples))),
  stringsAsFactors = FALSE
)

# Average all cells within each donor and target subtype, including zeros.
# Donors without captured target cells have missing, not zero, outcomes.
donor_df <- bind_rows(lapply(seq_len(nrow(targets)), function(i) {
  ct <- targets$CellType[i]
  gene <- targets$Gene[i]
  bind_rows(lapply(seq_len(nrow(donors)), function(j) {
    selected <- which(!is.na(md$Patient) & md$Patient == donors$Patient[j] &
                        !is.na(md$CellType) & md$CellType == ct &
                        !is.na(md$Group))
    values <- as.numeric(expr_mat[gene, selected, drop = TRUE])
    values <- values[is.finite(values)]
    n_valid <- length(values)
    n_positive <- if (n_valid) sum(values > positive_threshold) else NA_integer_
    data.frame(
      Patient = donors$Patient[j], Group = donors$Group[j],
      Interaction = targets$Interaction[i], Gene = gene, CellType = ct,
      Pair = targets$Pair[i], N_Cells = length(selected), N_Valid = n_valid,
      N_Positive = n_positive,
      MeanExpression = if (n_valid) mean(values) else NA_real_,
      PositivePct = if (n_valid) 100 * n_positive / n_valid else NA_real_
    )
  }))
}))
donor_df$Group <- factor(donor_df$Group, levels = c("AH", "ATH-A"))
donor_df$Pair <- factor(donor_df$Pair, levels = targets$Pair)
if (any(donor_df$N_Valid == 0)) {
  warning("Some donor-target combinations have no valid cells; see donor CSV.")
}
plot_df <- donor_df %>%
  pivot_longer(c(MeanExpression, PositivePct),
               names_to = "Metric", values_to = "Value")
valid_df <- plot_df %>% dplyr::filter(is.finite(Value))
summary_df <- valid_df %>%
  group_by(Metric, Pair, Group) %>%
  summarise(Mean = mean(Value),
            SEM = if (n() > 1) sd(Value) / sqrt(n()) else NA_real_,
            N_Donors = n(), .groups = "drop")

# Two-sided unpaired Welch t-tests, using donors as independent replicates.
test_df <- plot_df %>%
  group_by(Metric, Pair) %>%
  group_modify(~ {
    x <- .x$Value[.x$Group == "AH" & is.finite(.x$Value)]
    y <- .x$Value[.x$Group == "ATH-A" & is.finite(.x$Value)]
    p <- if (length(x) >= 2 && length(y) >= 2) {
      tryCatch(t.test(x, y, paired = FALSE, var.equal = FALSE,
                      alternative = "two.sided")$p.value,
               error = function(e) NA_real_)
    } else NA_real_
    data.frame(N_AH = length(x), N_ATHA = length(y),
               Mean_AH = if (length(x)) mean(x) else NA_real_,
               Mean_ATHA = if (length(y)) mean(y) else NA_real_,
               P_value = p)
  }) %>%
  ungroup() %>%
  group_by(Metric) %>%
  mutate(Q_BH = p.adjust(P_value, method = "BH", n = nrow(targets))) %>%
  ungroup() %>%
  mutate(P_Label = ifelse(is.na(P_value), "P = NA",
                          paste0("P = ", sprintf("%.3g", P_value))))

# Position P values above the donor points and error bars.
point_top <- valid_df %>%
  group_by(Metric, Pair) %>%
  summarise(PointTop = max(Value), .groups = "drop")
bar_top <- summary_df %>%
  group_by(Metric, Pair) %>%
  summarise(BarTop = max(Mean + coalesce(SEM, 0)), .groups = "drop")
label_df <- test_df %>%
  left_join(point_top, by = c("Metric", "Pair")) %>%
  left_join(bar_top, by = c("Metric", "Pair")) %>%
  mutate(Top = pmax(coalesce(PointTop, 0), coalesce(BarTop, 0)),
         LabelY = Top * 1.15 + ifelse(Metric == "PositivePct", 2, 0.01))

make_donor_plot <- function(metric, y_label, plot_title) {
  ggplot(dplyr::filter(summary_df, Metric == metric),
         aes(x = Group, y = Mean, fill = Group)) +
    geom_col(width = 0.7, alpha = 0.8, color = "black", linewidth = 0.3) +
    geom_errorbar(aes(ymin = pmax(0, Mean - SEM), ymax = Mean + SEM),
                  width = 0.18, linewidth = 0.6, na.rm = TRUE) +
    geom_point(data = dplyr::filter(valid_df, Metric == metric),
               aes(x = Group, y = Value, fill = Group), inherit.aes = FALSE,
               position = position_jitter(width = 0.07, height = 0, seed = 1234),
               color = "#354B52", size = 2, shape = 21) +
    geom_text(data = dplyr::filter(label_df, Metric == metric),
              aes(x = 1.5, y = LabelY, label = P_Label),
              inherit.aes = FALSE, size = 3.8) +
    facet_wrap(~ Pair, scales = "free_y", ncol = 2, drop = FALSE) +
    scale_x_discrete(limits = c("AH", "ATH-A"), drop = FALSE) +
    scale_fill_manual(values = group_colors, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(x = NULL, y = y_label, title = plot_title#, subtitle = "Mean +/- SEM; each point = one donor; unpaired Welch t-test"
         ) +
    theme_classic(base_size = 12) +
    theme(legend.position = "none", strip.background = element_blank(),
          strip.text = element_text(face = "bold", size = 11),
          plot.title = element_text(face = "bold"),
          panel.spacing = grid::unit(1.1, "lines"))
}

# Export donor-level metrics, statistical results, and figures.
write.csv(donor_df, "AH_ATH_LR_donor_metrics.csv", row.names = FALSE)
stats_out <- test_df %>%
  mutate(Pair = as.character(Pair)) %>%
  left_join(targets, by = "Pair") %>%
  dplyr::select(Interaction, Gene, CellType, Metric, N_AH, N_ATHA,
                Mean_AH, Mean_ATHA, P_value, Q_BH)
write.csv(stats_out, "AH_ATH_LR_donor_ttests.csv", row.names = FALSE)

p_expr <- make_donor_plot("MeanExpression", "Mean log-normalized expression", "Ligand-receptor gene expression by donor")
p_pos <- make_donor_plot("PositivePct", "Expression-positive cells (%)", "Ligand-receptor gene positivity by donor")
                         
pdf("AH_ATH_LR_mean_expression.pdf", width = 4, height = 8)
print(p_expr)
dev.off()
pdf("AH_ATH_LR_positive_proportion.pdf", width = 4, height = 8)
print(p_pos)
dev.off()

