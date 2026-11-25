# Spatial transcriptomics: CD4+ T-cell subclustering and downstream analysis
.libPaths("/home/lzh/miniconda3/envs/singlecell/lib/R/library")
setwd("/mnt/data_nvme1/lzh/test/AH/atlas/recluster/CD4T")
rm(list = ls())
set.seed(1234)
library(Seurat)
library(dplyr)
library(limma)
library(Matrix)
library(magrittr)
library(stringr)
library(patchwork)
library(clustree)
library(paletteer)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(forcats)
theme_global <- theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank()) + theme(aspect.ratio = 1)
spatial.all <- readRDS("/mnt/data_nvme1/lzh/test/AH/AH_all_spatial.RDS")
str(spatial.all@meta.data)
table(spatial.all$celltype)

# 01. Subset CD4+ T, Tfh, and Treg cells
objd.t <- subset(spatial.all, celltype == "CD4_T"| celltype == "Tfh"| celltype == "Treg")
DefaultAssay(objd.t) <- "Spatial"

# 02. Normalize and recluster CD4+ T cells
objd.t <- NormalizeData(objd.t, normalization.method = "LogNormalize", scale.factor = 10000) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>%
  ScaleData()
objd.t <- RunPCA(objd.t, npcs = 30)
pdf("pre/ElbowPlot.pdf", width = 8, height = 5)
ElbowPlot(objd.t, ndims = 30)
dev.off()
pc.select <- 1:20
objd.t@meta.data <- objd.t@meta.data %>% select(-starts_with("Spatial_snn_res"))

# Harmony integration
library(harmony)
pca_embed <- Embeddings(objd.t, reduction = "pca")[, pc.select, drop = FALSE]
h_embed <- harmony::HarmonyMatrix(data_mat = pca_embed,
  meta_data = objd.t@meta.data, vars_use = "Patient",  # Batch correction by patient
  do_pca = FALSE, verbose = TRUE)
objd.t[["harmony"]] <- CreateDimReducObject(embeddings = h_embed, key = "harmony_", assay = DefaultAssay(objd.t))
objd.t <- FindNeighbors(objd.t, reduction = "harmony", dims = pc.select)
objd.t <- FindClusters(objd.t, resolution = seq(0.1, 0.8, 0.1))
objd.t <- RunUMAP(objd.t, reduction = "harmony", dims = pc.select)
str(objd.t@meta.data)
table(objd.t$Spatial_snn_res.0.5)
table(objd.t$Spatial_snn_res.0.6)

# Check T-cell identity
DefaultAssay(objd.t)<-"Spatial"
pdf("pre/identity.pdf", width = 6, height = 6)
for (i in c("CD3D", "CD3E", "TRAC", "CD4","IL7R","TRBC2")) {
  try({
    p = FeaturePlot(objd.t, features = i, pt.size = 0.8, cols = c("grey", "red"),
                    raster = F, reduction="umap") + theme_global
    print(p)
  },silent = T)
}
dev.off()

# Set cluster order
library(dplyr)
objd.t$seurat_clusters <- objd.t$Spatial_snn_res.0.6
objd.t$seurat_clusters <- as.character(objd.t$seurat_clusters)
objd.t$seurat_clusters <- factor(objd.t$seurat_clusters, levels = c("0","1","2","3","4","5","6","7","8","9","10","11","12"))
objd.t$seurat_clusters <- as.character(objd.t$seurat_clusters)
objd.t$seurat_clusters <- paste("C", objd.t$seurat_clusters, sep ="")
objd.t$seurat_clusters <- factor(objd.t$seurat_clusters, levels = paste0("C", c(0:12)))
levels <- paste0("C", c(0:12))
objd.t$seurat_clusters <- factor(objd.t$seurat_clusters, levels = levels)
levels(objd.t$seurat_clusters)
table(objd.t$seurat_clusters)
pdf('pre/cluster_umap.pdf', useDingbats=F, width = 6, height = 6)
DimPlot(objd.t, reduction = "umap", label = TRUE, repel = TRUE, raster = FALSE,
    group.by = 'seurat_clusters', pt.size = 0.6)  + theme_global + ggtitle("Cluster")
dev.off()

# 03. CD4+ T-cell subtype annotation
markers<-c("CCR7","TCF7","LEF1","NOSIP","SELL",  # Naive T
           "IL7R","RUNX2","ANK3","KLRB1","LTB",  # Memory T
           "ANXA1","S100A4","S100A10","PCNX1","ZBTB16",  # Central_memory
           "FOXP3","IL2RA","RGS1","IKZF2","FCRL3","CLNK",  # Treg
           "CXCR5","PASK","ICOS","BCL6","PDCD1","TOX2",  # Tfh
           "MKI67", "CDK1", "STMN1",  # Proliferative
           "HAVCR2","CTLA4","LAG3","TIGIT",  # Inhibitory
           "CXCR3","IFNG","TBX21",  # Th1
           "IL4","IL5","GATA3","CCR4",  # Th2
           "CCR6","IL17A","RORC","IL23R",  # Th17
           "CD3D","CD3G","CD3E",
           "MS4A1","CD79A","CD74","IGHM"
)
Idents(objd.t) <- "seurat_clusters"
pdf('pre/markerBubble.pdf',useDingbats=F, width = 15, height = 8)
DotPlot(objd.t, features = markers, assay = 'Spatial') + scale_color_distiller(palette = "RdYlBu") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
dev.off()
pdf("pre/marker.umap.pdf", width = 6, height = 6)
for (i in markers) {
  try({
    p = FeaturePlot(objd.t, features = i, pt.size = 0.8, cols = c("grey", "red"),
                    raster = F, reduction="umap") + theme_global
    print(p)
  },silent = T)
}
dev.off()
cluster2celltype_t <- c("C0"="Tfh", "C1"="Memory",
                         "C2"="Naive", "C3"="Tfh",
                         "C4"="Naive", "C5"="Memory",
                         "C6"="Mix", "C7"="Proliferative",
                         "C8"="Treg", "C9"="Tfh",
                         "C10"="Naive", "C11"="Central_memory",
                         "C12"="Memory")
objd.t[['celltype_t']]<-unname(cluster2celltype_t[objd.t@meta.data$seurat_clusters])
objd.t <- subset(objd.t, celltype_t != "Mix")
dim(objd.t)
objd.t$celltype_t <- factor(objd.t$celltype_t, levels = c("Central_memory","Memory", "Naive", "Proliferative", "Tfh", "Treg"))
pal1 <- c("#984EA3","#BC9DCC", "#8DBA8F","#CB5C5B", "#7FB2D5", "#F5B47F")
pdf('pre/celltype_umap.pdf',useDingbats=F,width = 6, height = 6)
DimPlot(objd.t, reduction = "umap", label = TRUE, repel = TRUE, cols = pal1,
    group.by = 'celltype_t', pt.size = 0.6)  + theme_global + ggtitle("Celltype")
dev.off()

# Visualize subtype markers
library(ggplot2)
library(grid)
markers.show<-list(
  Central_memory = c("ANXA1","S100A4","S100A10"),
  Memory = c("IL7R","ANK3","LTB"),
  Naive = c("CCR7","TCF7","LEF1"),
  Proliferative = c("MKI67", "CDK1", "STMN1"),
  Tfh = c("ICOS","BCL6","PDCD1"),
  Treg = c("FOXP3","IL2RA","IKZF2")
)
markers.show2 <- c("ANXA1","S100A4","S100A10", "IL7R","ANK3","LTB", "CCR7","TCF7","LEF1",
                   "MKI67", "CDK1", "STMN1", "ICOS","BCL6","PDCD1","FOXP3","IL2RA","IKZF2")
Idents(objd.t) <- "celltype_t"
pdf("pre/DotPlot.Ttype.pdf",width = 8, height = 4)
DotPlot(objd.t, features = markers.show2,  cols = c("grey", "red"),
    dot.min = 0, dot.scale = 6, idents = NULL, group.by = "celltype_t",
        split.by = NULL, cluster.idents = FALSE, scale = TRUE, scale.by = "radius", scale.min = NA, scale.max = NA)+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_y_discrete(limits = levels(objd.t$celltype_t)) +
  theme(panel.border = element_rect(color="black"),
        panel.spacing = unit(1, "mm"),
        strip.text = element_text(margin=margin(b=3, unit="mm")),
        strip.placement = 'outlet',
        axis.line = element_blank()) + labs(x="", y="")
DotPlot(objd.t, features = markers.show2,
    dot.min = 0, dot.scale = 6, idents = NULL, group.by = "celltype_t",
        split.by = NULL, cluster.idents = FALSE, scale = TRUE, scale.by = "radius", scale.min = NA, scale.max = NA) +
    scale_color_distiller(palette = "RdYlBu") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_y_discrete(limits = levels(objd.t$celltype_t)) + labs(x="", y="") +
  theme(panel.border = element_rect(color = "black", size = 0.5, fill = NA))
dev.off()

# 04. Differential markers by cluster
levels <- paste0("C", c(0:5, 7:12))
objd.t$seurat_clusters <- factor(objd.t$seurat_clusters, levels = levels)
table(objd.t$seurat_clusters)
DefaultAssay(objd.t) <- "Spatial"
Idents(objd.t) <- "seurat_clusters"
logFCfilter=0.50
adjPvalFilter=1e-5
mrksT <- FindAllMarkers(objd.t, logfc.threshold = 0.25, min.pct = 0.2, only.pos = TRUE)
length(mrksT$gene)
sig.markers<-mrksT[(abs(as.numeric(mrksT$avg_log2FC))>logFCfilter & as.numeric(mrksT$p_val_adj)<adjPvalFilter),]
length(sig.markers$gene)
sig.markers$diff.pct <- sig.markers$pct.1 - sig.markers$pct.2  # Difference in detection fractions
sig.markers$diff.pct.fc <- sig.markers$pct.1/sig.markers$pct.2  # Ratio of detection fractions
sig.markers2 <- sig.markers[which(sig.markers$diff.pct > 0.0),]
sig.markers2 <- sig.markers[which(sig.markers$diff.pct.fc > 1.25),]
length(sig.markers2$gene)
sig.markers <- sig.markers %>%  group_by(cluster) %>% arrange(desc(avg_log2FC), .by_group = TRUE)
write.csv(sig.markers, "sig/sig.markers.csv", row.names = F)
saveRDS(sig.markers, file = "sig/sig.markers.RDS")
sig.markers2 <- sig.markers2 %>%  group_by(cluster) %>% arrange(desc(avg_log2FC), .by_group = TRUE)
write.csv(sig.markers2, "sig/sig.markers2.csv", row.names = F)
saveRDS(sig.markers2, file = "sig/sig.markers2.RDS")

# Cluster quality control
pdf("sig/qc.pdf",width = 15, height = 5)
VlnPlot(objd.t, features=c("nFeature_Spatial","nCount_Spatial","percent.mt"),
        pt.size=0, ncol=3)
dev.off()

# Assign named subclusters
cluster2cellname <- c("C0"="TOX2_Tfh", "C1"="CMTM8_Memory",
                      "C2"="TRABD2A_Naive", "C3"="TNFRSF1B_Tfh",
                      "C4"="PDE3B_Naive", "C5"="AHNAK_Memory",
                      "C7"="MKI67_Proliferative",
                      "C8"="IKZF2_Treg", "C9"="GFOD1_Tfh",
                      "C10"="FCRL3_Naive", "C11"="MYB_Central_memory",
                      "C12"="ITGAE_Memory")
objd.t[['cellname_t']]<-unname(cluster2cellname[objd.t@meta.data$seurat_clusters])
objd.t$cellname_t <- factor(objd.t$cellname_t, levels = c("MYB_Central_memory", "AHNAK_Memory", "CMTM8_Memory", "ITGAE_Memory", "FCRL3_Naive", "PDE3B_Naive", "TRABD2A_Naive",
                            "MKI67_Proliferative", "TNFRSF1B_Tfh","GFOD1_Tfh", "TOX2_Tfh", "IKZF2_Treg"))
table(objd.t$cellname_t)

# Visualize annotated clusters
Idents(objd.t) <- "cellname_t"
pal <-c("#FFD569","#F7D8B7","#F5B47F","#FD8052","#8DD2C5","#8DBA8F","#00A087","#D37295","#B3CDE3","#7FB2D5","#4581B2","#BC9DCC")
pdf('sig/cellname_umap.pdf',useDingbats=F,width = 8, height = 6)
DimPlot(objd.t, reduction = "umap", label = TRUE, repel = TRUE, cols = pal,
    group.by = 'cellname_t', pt.size = 0.8)  + theme_global + ggtitle("Cluster")
dev.off()

# Spatial distribution by patient
pdf("sig/CD4T_spatial_byPatient_cellname.pdf", width = 8, height = 6)
for (p in levels(objd.t$Patient)) {
  message("Plotting Patient: ", p)
  obj.p <- subset(objd.t, subset = Patient == p)
  p2 <- DimPlot(obj.p, reduction = "spatial", group.by = "cellname_t", pt.size = 1, cols = pal) + ggtitle(p) + theme_global
  print(p2)
}
dev.off()

# 05. Differential markers by named subcluster
DefaultAssay(objd.t) <- "Spatial"
Idents(objd.t) <- "cellname_t"
logFCfilter=0.50
adjPvalFilter=1e-5
mrksT <- FindAllMarkers(objd.t, logfc.threshold = 0.25, min.pct = 0.2, only.pos = TRUE)
length(mrksT$gene)
sig.markers<-mrksT[(abs(as.numeric(mrksT$avg_log2FC))>logFCfilter & as.numeric(mrksT$p_val_adj)<adjPvalFilter),]
length(sig.markers$gene)
sig.markers$diff.pct <- sig.markers$pct.1 - sig.markers$pct.2  # Difference in detection fractions
sig.markers$diff.pct.fc <- sig.markers$pct.1/sig.markers$pct.2  # Ratio of detection fractions
sig.markers2 <- sig.markers[which(sig.markers$diff.pct > 0.0),]
sig.markers2 <- sig.markers[which(sig.markers$diff.pct.fc > 1.25),]
length(sig.markers2$gene)
sig.markers <- sig.markers %>%  group_by(cluster) %>% arrange(desc(avg_log2FC), .by_group = TRUE)
write.csv(sig.markers, "sig/sig.markers.name.csv", row.names = F)
saveRDS(sig.markers, file = "sig/sig.markers.name.RDS")
sig.markers2 <- sig.markers2 %>%  group_by(cluster) %>% arrange(desc(avg_log2FC), .by_group = TRUE)
write.csv(sig.markers2, "sig/sig.markers2.name.csv", row.names = F)
saveRDS(sig.markers2, file = "sig/sig.markers2.name.RDS")

# Select top-ranked genes
top50T<-sig.markers2 %>% group_by(cluster) %>% top_n(n=50, wt=avg_log2FC)
top30T<-sig.markers2 %>% group_by(cluster) %>% top_n(n=30, wt=avg_log2FC)
top20T<-sig.markers2 %>% group_by(cluster) %>% top_n(n=20, wt=avg_log2FC)
top10T<-sig.markers2 %>% group_by(cluster) %>% top_n(n=10, wt=avg_log2FC)
top5T<-sig.markers2 %>% group_by(cluster) %>% top_n(n=5, wt=avg_log2FC)
write.csv(top50T, "sig/top50T.csv", row.names = F)
write.csv(top30T, "sig/top30T.csv", row.names = F)
write.csv(top20T, "sig/top20T.csv", row.names = F)
write.csv(top10T, "sig/top10T.csv", row.names = F)
write.csv(top5T, "sig/top5T.csv", row.names = F)

# Marker heatmaps
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
objd.t <- ScaleData(objd.t, verbose = FALSE)
mapal <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256)
pdf(file="sig/top10HeatmapT.pdf", width=16, height=14)
DoHeatmap(objd.t, features = as.character(unique(top10T$gene)),
        slot = "scale.data", raster = F, size = 3,
        group.by = "cellname_t", group.colors = pal) +
        scale_fill_gradient2(low = '#0099CC',mid = 'white',high = '#CC0033',
                       name = 'Z-score')
dev.off()
pdf(file="sig/top5HeatmapT.pdf", width=12, height=10.5)
DoHeatmap(objd.t, features = as.character(unique(top5T$gene)),
        slot = "scale.data", raster = F, size = 3,
        group.by = "cellname_t", group.colors = pal) +
        scale_fill_gradient2(low = '#0099CC',mid = 'white',high = '#CC0033',
                       name = 'Z-score')
dev.off()

# Cluster-averaged expression
mapal <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256)
cluster.averages <- AverageExpression(objd.t,return.seurat = TRUE)
pdf(file="sig/top10HeatmapT_mean.pdf", width=10, height=13)
DoHeatmap(cluster.averages, features = top10T$gene, slot = "scale.data", group.by = "ident",
          group.colors = pal, size = 3, draw.lines = FALSE) + scale_fill_gradientn(colours = rev(mapal))
dev.off()
saveRDS(objd.t, file = '/mnt/data_nvme1/lzh/test/AH/atlas/recluster/CD4T/CD4T.all.RDS')

# 06. Cell-type composition by group
library(ggplot2)
table(objd.t$TH_meta, objd.t$cellname_t)
rt1 <- data.frame(count=c(401, 1275, 3380, 153, 670, 1986, 1818, 676, 1714, 614, 3878, 883),
                   category=levels(objd.t$cellname_t))
rt1$category <- factor(rt1$category, levels=levels(objd.t$cellname_t))
rt1$fraction <- rt1$count / sum(rt1$count)
rt1$ymax <- cumsum(rt1$fraction)
rt1$ymin <- c(0, head(rt1$ymax, n=-1))
rt1$labelPosition <- (rt1$ymax + rt1$ymin) / 2
rt1$label <- paste0(rt1$category, "\n value: ", round(rt1$fraction*100,2), "%")
rt2 <- data.frame(count=c(811, 1067, 4838, 124, 576, 2217, 2821, 1112, 2493, 951, 5125, 901),
                   category=levels(objd.t$cellname_t))
rt2$category <- factor(rt2$category, levels=levels(objd.t$cellname_t))
rt2$fraction <- rt2$count / sum(rt2$count)
rt2$ymax <- cumsum(rt2$fraction)
rt2$ymin <- c(0, head(rt2$ymax, n=-1))
rt2$labelPosition <- (rt2$ymax + rt2$ymin) / 2
rt2$label <- paste0(rt2$category, "\n value: ", round(rt2$fraction*100,2), "%")
pdf("sig/pie.NY.pdf",width = 4, height = 4)
ggplot(rt1, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
  geom_rect() +
  geom_text(x=2, aes(y=labelPosition, label=label, color=category), size=1.8) +
  scale_fill_manual(values=pal) +
  scale_color_manual(values=pal) + coord_polar(theta="y") +
  xlim(c(-1, 4)) + theme_void() + theme(legend.position = "none")
ggplot(rt2, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
  geom_rect() +
  geom_text(x=2, aes(y=labelPosition, label=label, color=category), size=1.8) +
  scale_fill_manual(values=pal) +
  scale_color_manual(values=pal) + coord_polar(theta="y") +
  xlim(c(-1, 4)) + theme_void() + theme(legend.position = "none")
dev.off()

# Subcluster proportions
table(objd.t$TH_meta, objd.t$cellname_t)
num2 <- table(objd.t$TH_meta, objd.t$cellname_t)
write.csv(num2, "sig/num2.csv")
num2 <- as.data.frame(read.csv("sig/num2.csv", header = T, row.names = 1))
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
pct <- num2 %>%
  mutate(across(everything(), ~ .x / sum(.x, na.rm = TRUE) * 100))
dat2 <- pct %>%
  rownames_to_column("Group") %>%
  pivot_longer(cols = -Group, names_to = "Cluster", values_to = "percent")
dat2$Group <- factor(dat2$Group, levels = c("N", "Y"))
dat2$Cluster <- factor(dat2$Cluster, levels = c(
  "MYB_Central_memory", "AHNAK_Memory", "CMTM8_Memory", "ITGAE_Memory",
  "FCRL3_Naive", "PDE3B_Naive", "TRABD2A_Naive", "MKI67_Proliferative",
  "TNFRSF1B_Tfh", "GFOD1_Tfh", "TOX2_Tfh", "IKZF2_Treg"
))
dat2$label <- paste0(sprintf("%.2f", dat2$percent), "%")
pdf("sig/percent.cl.pdf", width = 8, height = 6)
ggplot(dat2, aes(x = Cluster, y = percent, fill = Group)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c(N = "#9EC4BE", Y = "#FAC795")) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  labs(x = "", y = "Percentage", fill = "Group") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    legend.position = "top"
  )
dev.off()

# 07. Tfh pathway enrichment
library(msigdbr)
library(clusterProfiler)
library(enrichplot)
msigdbr_collections()
Hdb <- msigdbr(species = "Homo sapiens", category = "H") %>% dplyr::select(gs_name, gene_symbol)
Kdb <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG_LEGACY") %>% dplyr::select(gs_name, gene_symbol)
Wdb <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:WIKIPATHWAYS") %>% dplyr::select(gs_name, gene_symbol)
BPdb <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP") %>% dplyr::select(gs_name, gene_symbol)
df <- readRDS("/mnt/data_nvme1/lzh/test/AH/atlas/recluster/CD4T/sig/sig.markers2.name.RDS")
df <- read.csv("/mnt/data_nvme1/lzh/test/AH/atlas/recluster/CD4T/sig/top50T.csv")
table(df$cluster)
df2 <- df %>% filter(cluster %in% c("TNFRSF1B_Tfh","GFOD1_Tfh","TOX2_Tfh"))
pal.tfh <-c("#B3CDE3","#7FB2D5","#4581B2")

# Hallmark enrichment
hm1 <- enricher(df2[df2$cluster == "TNFRSF1B_Tfh",]$gene, TERM2GENE=Hdb)
hm2 <- enricher(df2[df2$cluster == "GFOD1_Tfh",]$gene, TERM2GENE=Hdb)
hm3 <- enricher(df2[df2$cluster == "TOX2_Tfh",]$gene, TERM2GENE=Hdb)
hm.all <- list(hm1@result, hm2@result, hm3@result)
names(hm.all) <- c("TNFRSF1B_Tfh","GFOD1_Tfh","TOX2_Tfh")
hm.all.fli <- hm.all
process_df.hk <- function(hm.res) {
  labels <- hm.res[order(hm.res$p.adjust, decreasing = TRUE), "ID"]
  hm.res$ID <- factor(hm.res$ID, levels = labels)
  hm.res$ID <- substr(hm.res$ID, 10, 100)
  hm.res$ID <- gsub("_", " ", hm.res$ID)
  return(hm.res)
}
hm.all.fli <- lapply(hm.all.fli, process_df.hk)
str(hm.all.fli)
process_df <- function(hm.res) {
  df <- hm.res[1:10, ]
  df$LogP <- -log10(df$p.adjust)
  df$labelx <- rep(0, nrow(df))
  df$labely <- seq_len(nrow(df))
  return(df)
}
hm.all.fli1 <- lapply(hm.all.fli, process_df)
str(hm.all.fli1)

# KEGG enrichment
km1 <- enricher(df2[df2$cluster == "TNFRSF1B_Tfh",]$gene, TERM2GENE=Kdb)
km2 <- enricher(df2[df2$cluster == "GFOD1_Tfh",]$gene, TERM2GENE=Kdb)
km3 <- enricher(df2[df2$cluster == "TOX2_Tfh",]$gene, TERM2GENE=Kdb)
km.all <- list(km1@result, km2@result, km3@result)
names(km.all) <- c("TNFRSF1B_Tfh","GFOD1_Tfh","TOX2_Tfh")
km.all.fli <- km.all
process_df.kegg <- function(km.res) {
  labels <- km.res[order(km.res$p.adjust, decreasing = TRUE), "ID"]
  km.res$ID <- factor(km.res$ID, levels = labels)
  km.res$ID <- substr(km.res$ID, 6, 100)
  km.res$ID <- gsub("_", " ", km.res$ID)
  return(km.res)
}
km.all.fli <- lapply(km.all.fli, process_df.kegg)
str(km.all.fli)
process_df <- function(km.res) {
  df <- km.res[1:10, ]
  df$LogP <- -log10(df$p.adjust)
  df$labelx <- rep(0, nrow(df))
  df$labely <- seq_len(nrow(df))
  return(df)
}
km.all.fli1 <- lapply(km.all.fli, process_df)
str(km.all.fli1)

# GO biological process enrichment
bm1 <- enricher(df2[df2$cluster == "TNFRSF1B_Tfh",]$gene, TERM2GENE=BPdb)
bm2 <- enricher(df2[df2$cluster == "GFOD1_Tfh",]$gene, TERM2GENE=BPdb)
bm3 <- enricher(df2[df2$cluster == "TOX2_Tfh",]$gene, TERM2GENE=BPdb)
bm.all <- list(bm1@result, bm2@result, bm3@result)
names(bm.all) <- c("TNFRSF1B_Tfh","GFOD1_Tfh","TOX2_Tfh")
bm.all.fli <- bm.all
process_df.bk <- function(bm.res) {
  labels <- bm.res[order(bm.res$p.adjust, decreasing = TRUE), "ID"]
  bm.res$ID <- factor(bm.res$ID, levels = labels)
  bm.res$ID <- substr(bm.res$ID, 6, 100)
  bm.res$ID <- gsub("_", " ", bm.res$ID)
  return(bm.res)
}
bm.all.fli <- lapply(bm.all.fli, process_df.bk)
str(bm.all.fli)
bm.all.fli1 <- lapply(bm.all.fli, process_df)
str(bm.all.fli1)
saveRDS(hm.all.fli, "sig/hm.all.tfh.RDS")
saveRDS(km.all.fli, "sig/km.all.tfh.RDS")
saveRDS(bm.all.fli, "sig/bm.all.tfh.RDS")

# Plot pathway results
library(ggplot2)
plot_pathway <- function(df, i) {
  if(!"LogP" %in% names(df) || !is.numeric(df$LogP)) {
    stop("The 'LogP' column must exist and be numeric.")
  }
  if(!"p.adjust" %in% names(df) || !is.numeric(df$p.adjust)) {
    stop("The 'p.adjust' column must exist and be numeric.")
  }
  df$LogP <- -log10(df$p.adjust)
  df$ID <- reorder(df$ID, df$LogP)
  p <- ggplot(data = df, aes(x = LogP, y = ID)) +
    geom_bar(stat = "identity", alpha = 0.5, fill = pal.tfh[i], width = 0.8) +
    geom_text(aes(label = ID, x = 0), size = 5, hjust = 0) +
    theme_classic() +
    theme(axis.text.y = element_blank(),
          axis.line.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.x = element_line(colour = 'black', linewidth = 1),
          axis.text.x = element_text(colour = 'black', size = 15),
          axis.ticks.x = element_line(colour = 'black', linewidth = 1),
          axis.title.x = element_text(colour = 'black', size = 16)) +
    xlab("-log10(p.adjust)") +
    theme(plot.title = element_text(size = 18, face = "bold")) +
    scale_x_continuous(expand = c(0, 0.5)) +
    theme(panel.border = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"))
  return(p)
}
pdf("sig/hallmark.tfh.pdf", width=7.5, height=5)
for (i in seq_along(hm.all.fli1)) {
  p <- plot_pathway(hm.all.fli1[[i]], i) + ggtitle(names(hm.all.fli1)[i])
  print(p)
}
dev.off()
pdf("sig/kegg.tfh.pdf", width=7.5, height=5)
for (i in seq_along(km.all.fli1)) {
  p <- plot_pathway(km.all.fli1[[i]], i) + ggtitle(names(km.all.fli1)[i])
  print(p)
}
dev.off()
pdf("sig/gobp.tfh.pdf", width=7.5, height=5)
for (i in seq_along(bm.all.fli1)) {
  p <- plot_pathway(bm.all.fli1[[i]], i) + ggtitle(names(bm.all.fli1)[i])
  print(p)
}
dev.off()

# 08. FM-Tfh versus GC-Tfh differential expression
library(Seurat)
library(EnhancedVolcano)
library(ggplot2)
objd.t.all <- readRDS("/mnt/data_nvme1/lzh/test/AH/atlas/recluster/CD4T/CD4T.all.RDS")
TFH <-subset(objd.t.all, celltype_t == "Tfh")
TFH$group <- ifelse(TFH$cellname_t == "TNFRSF1B_Tfh", "FM_Tfh", "GC_Tfh")
logFCfilter=0.50
adjPvalFilter=1e-5
diff_gene.tfh <- FindMarkers(TFH, min.pct = 0.25, logfc.threshold = 0.25, group.by = "group",
                            ident.1 = "FM_Tfh", ident.2 = "GC_Tfh", only.pos = FALSE, test.use="wilcox")
sig_gene.tfh <- subset(diff_gene.tfh, p_val_adj < adjPvalFilter & abs(avg_log2FC) > logFCfilter)
sig_gene.tfh <- sig_gene.tfh %>% arrange(desc(avg_log2FC), .by_group = TRUE)
saveRDS(sig_gene.tfh,"sig/diff_gene.tfh.RDS")
saveRDS(sig_gene.tfh,"sig/sig_gene.tfh.RDS")
write.csv(sig_gene.tfh,"sig/sig_gene.tfh.csv")

# Volcano plot
pdf("sig/diff_tfh.pdf",width = 6, height = 6.5)
EnhancedVolcano(diff_gene.tfh,lab = rownames(diff_gene.tfh),
                   x = 'avg_log2FC', y = 'p_val_adj',xlim = c(-2.5, 2.5),
                   pCutoff = adjPvalFilter, FCcutoff = logFCfilter,legendPosition = 'bottom',
                   pointSize = 1.5, labSize = 4, legendLabSize = 9,legendIconSize = 3.5, axisLabSize = 12,
                   title = "Tfh",subtitle = "",titleLabSize = 15, subtitleLabSize =12,
                   gridlines.major = FALSE,gridlines.minor = FALSE,
                   drawConnectors = TRUE,arrowheads = FALSE , max.overlaps = 16, shape = 16,
                   col=c("#D9D9D9", "#8DBA8F", "#4581B2", "#CB5C5B"), colAlpha = 1)
dev.off()

# 09. Monocle2 pseudotime analysis
library(Seurat)
library(igraph)
library(monocle)
library(tidyverse)
library(ggplot2)
rm(list = ls())
setwd("/mnt/data_nvme1/lzh/test/AH/atlas/recluster/CD4T/monocle")
set.seed(1234)
theme_global <- theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank()) + theme(aspect.ratio = 1)
objd.t.all <- readRDS("/mnt/data_nvme1/lzh/test/AH/atlas/recluster/CD4T/CD4T.all.RDS")
DefaultAssay(objd.t.all) <- "Spatial"
dim(objd.t.all)
str(objd.t.all@meta.data)
table(objd.t.all$cellname_t)
table(objd.t.all$celltype_t)
objd.t.fli <-  objd.t.all[,objd.t.all$celltype_t %in% c("Naive","Tfh")]
dim(objd.t.fli)
levels(objd.t.fli)
saveRDS(objd.t.fli, file = "/mnt/data_nvme1/lzh/test/AH/atlas/recluster/CD4T/objd.t.fli.RDS")

# Construct the CellDataSet
expr_matrix<-as(as.matrix(objd.t.fli@assays$Spatial@counts),'sparseMatrix')  # UMI count matrix
p_data<-objd.t.fli@meta.data
f_data<-data.frame(gene_short_name=row.names(objd.t.fli),row.names = row.names(objd.t.fli))
pd<-new("AnnotatedDataFrame",data=p_data)
fd<-new("AnnotatedDataFrame",data=f_data)
cds<-newCellDataSet(expr_matrix,phenoData = pd,featureData = fd,
                    lowerDetectionLimit = 0.5, expressionFamily = negbinomial.size())  # Negative binomial model
cds<-estimateSizeFactors(cds)
cds<-estimateDispersions(cds)
cds<-detectGenes(cds,min_expr = 0.1)

# Select ordering genes
objd.t.fli <- FindVariableFeatures(objd.t.fli, nfeatures = 1000)
express_genes <- VariableFeatures(objd.t.fli)
cds <- setOrderingFilter(cds,express_genes)
plot_ordering_genes(cds)

# Infer cell trajectories
cds<-reduceDimension(cds, reduction_method = "DDRTree", max_components = 2, num_dim = 6,
                     norm_method = 'log', residualModelFormulaStr = "~Patient")  # Adjust for patient effects
cds<-orderCells(cds)
library(ggplot2)
library(RColorBrewer)
library(paletteer)
library(ggpubr)
library(ggsci)
pData(cds)$CCR7 <- log2(exprs(cds)['CCR7',]+1)  # naive T
pData(cds)$LEF1 <- log2(exprs(cds)['LEF1',]+1)
pData(cds)$BCL6 <- log2(exprs(cds)['BCL6',]+1)  # Tfh
pData(cds)$CXCR5 <- log2(exprs(cds)['CXCR5',]+1)
pData(cds)$PDCD1 <- log2(exprs(cds)['PDCD1',]+1)
pData(cds)$CXCR3 <- log2(exprs(cds)['CXCR3',]+1)
pData(cds)$CCR6 <- log2(exprs(cds)['CCR6',]+1)
pData(cds)$KLF2 <- log2(exprs(cds)['KLF2',]+1)
pData(cds)$S1PR1 <- log2(exprs(cds)['S1PR1',]+1)
pData(cds)$S1PR2 <- log2(exprs(cds)['S1PR2',]+1)
pData(cds)$CCL19 <- log2(exprs(cds)['CCL19',]+1)
pData(cds)$CCL21 <- log2(exprs(cds)['CCL21',]+1)
pData(cds)$GPR183 <- log2(exprs(cds)['GPR183',]+1)
pData(cds)$SELL <- log2(exprs(cds)['SELL',]+1)
pData(cds)$IL7R <- log2(exprs(cds)['IL7R',]+1)
pData(cds)$BTLA <- log2(exprs(cds)['BTLA',]+1)
pData(cds)$SELPLG <- log2(exprs(cds)['SELPLG',]+1)
pData(cds)$GZMK <- log2(exprs(cds)['GZMK',]+1)
pData(cds)$GZMA <- log2(exprs(cds)['GZMA',]+1)
pData(cds)$TBX21 <- log2(exprs(cds)['TBX21',]+1)  # Th1
pData(cds)$GATA3 <- log2(exprs(cds)['GATA3',]+1)  # Th2
pData(cds)$RORC <- log2(exprs(cds)['RORC',]+1)  # Th17
pData(cds)$FOXP3 <- log2(exprs(cds)['FOXP3',]+1)
pData(cds)$ANK3 <- log2(exprs(cds)['ANK3',]+1)  # memory T

# Visualize candidate genes
pdf("gene.exp.vln.pdf",width = 6, height = 4)
for(i in c("CCR7","LEF1","BCL6","CXCR5","PDCD1","CXCR3","CCR6","KLF2","S1PR1","S1PR2","CCL19","CCL21","GPR183","SELL","IL7R","BTLA","SELPLG","GZMK","GZMA","TBX21","GATA3","RORC","FOXP3","ANK3")){
  p <- VlnPlot(objd.t.fli, features = i, group.by = "cellname_t", pt.size = 0.01, cols = pal1) + xlab("")
  print(p)
}
dev.off()
objd.t.fli$cellname_t <- factor(objd.t.fli$cellname_t,levels = c("FCRL3_Naive", "PDE3B_Naive", "TRABD2A_Naive", "TNFRSF1B_Tfh","GFOD1_Tfh","TOX2_Tfh"))
pdf("gene.exp.bubble.pdf", width = 8, height = 4)
DotPlot(objd.t.fli, features = c("CCR7","LEF1","BCL6","CXCR5","PDCD1","CXCR3","CCR6","KLF2","S1PR1","S1PR2","CCL19","CCL21","GPR183","SELL","IL7R","BTLA","SELPLG","GZMK","GZMA","TBX21","GATA3","RORC","FOXP3","ANK3"),
    dot.min = 0, dot.scale = 6, idents = NULL, group.by = "cellname_t",
        split.by = NULL, cluster.idents = FALSE, scale = TRUE, scale.by = "radius", scale.min = NA, scale.max = NA) +
    scale_color_distiller(palette = "RdYlBu") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_y_discrete(limits = levels(objd.t.fli$cellname_t)) + labs(x="", y="") +
  theme(panel.border = element_rect(color = "black", size = 0.5, fill = NA))
dev.off()
pdf("gene.exp.pdf",width = 5,height = 5)
plot_cell_trajectory(cds,color_by="CCR7", cell_size = 0.8) + scale_color_material("deep-orange")
plot_cell_trajectory(cds,color_by="LEF1", cell_size = 0.8) + scale_color_material("deep-orange")
plot_cell_trajectory(cds,color_by="BCL6", cell_size = 0.8) + scale_color_material("deep-orange")
plot_cell_trajectory(cds,color_by="CXCR5", cell_size = 0.8) + scale_color_material("deep-orange")
plot_cell_trajectory(cds,color_by="PDCD1", cell_size = 0.8) + scale_color_material("deep-orange")
plot_cell_trajectory(cds,color_by="KLF2", cell_size = 0.8) + scale_color_material("deep-orange")
plot_cell_trajectory(cds,color_by="GPR183", cell_size = 0.8) + scale_color_material("deep-orange")
plot_cell_trajectory(cds,color_by="SELPLG", cell_size = 0.8) + scale_color_material("deep-orange")
plot_cell_trajectory(cds,color_by="GZMK", cell_size = 0.8) + scale_color_material("deep-orange")
plot_cell_trajectory(cds,color_by="GZMA", cell_size = 0.8) + scale_color_material("deep-orange")
plot_cell_trajectory(cds,color_by="TBX21", cell_size = 0.8) + scale_color_material("deep-orange")
plot_cell_trajectory(cds,color_by="GATA3", cell_size = 0.8) + scale_color_material("deep-orange")
plot_cell_trajectory(cds,color_by="RORC", cell_size = 0.8) + scale_color_material("deep-orange")
plot_cell_trajectory(cds,color_by="FOXP3", cell_size = 0.8) + scale_color_material("deep-orange")
plot_cell_trajectory(cds,color_by="ANK3", cell_size = 0.8) + scale_color_material("deep-orange")
dev.off()

# Pseudotime trajectory plots
pal1 <-c("#8DD2C5","#8DBA8F","#00A087","#B3CDE3","#7FB2D5","#4581B2")
pdf("Pseudotime.pdf",width = 5, height = 5)
plot_cell_trajectory(cds,color_by = "Pseudotime", cell_size = 0.8)
plot_cell_trajectory(cds,color_by = "Pseudotime", cell_size = 0.8)+scale_color_gradientn(colors = rev(brewer.pal(10, "RdYlBu")))
plot_cell_trajectory(cds,color_by = "cellname_t", cell_size = 0.8)+scale_color_manual(values = pal1)
plot_cell_trajectory(cds,color_by = "TH_meta", cell_size = 0.8)+scale_color_manual(values = c("#9EC4BE","#FAC795"))
plot_cell_trajectory(cds,color_by = "celltype_t", cell_size = 0.8)+scale_color_manual(values = c("#8DBA8F", "#7FB2D5"))
plot_cell_trajectory(cds,color_by = "State", cell_size = 0.8)
dev.off()
pdf("Pseudotime_all_type.pdf",width = 12, height = 8)
plot_cell_trajectory(cds,color_by = "cellname_t", cell_size = 0.8) + facet_wrap(~cellname_t,nrow=2) + scale_color_manual(values = pal1)
dev.off()
pdf("Pseudotime_all_State.pdf",width = 35, height = 20)
plot_cell_trajectory(cds,color_by = "State", cell_size = 0.8) + facet_wrap(~State,nrow=2)
dev.off()

# Branching trajectory plots
pdf("Pseudotime_tree.pdf",width = 6,height = 6)
plot_complex_cell_trajectory(cds,x=1,y=2,color_by = "State")
dev.off()

# Cell density along pseudotime
df<-pData(cds)
pdf("Pseudotime_density.pdf",width = 8, height = 4)
ggplot(df, aes(Pseudotime,colour=cellname_t, fill=cellname_t))+
  scale_color_manual(values = pal1)+
  scale_fill_manual(values = pal1)+
  geom_density(bw=0.5,size=0.5,alpha=0.5)+theme_classic2()+
  theme(axis.text.x = element_text(colour="black",size = 8),
        axis.text.y = element_text(colour="black",size = 8),axis.title = element_text(size=10),
        legend.text = element_text(size=8), legend.title = element_text(size=10), legend.key.size = unit(10, "pt"))
dev.off()
df$cellname_t <- as.character(df$cellname_t)
df$cellname_t <- factor(df$cellname_t, levels = levels(objd.t.fli))
write.csv(pData(cds),"pseudotime.csv")
saveRDS(cds, file = "cds.RDS")
saveRDS(df, file = "pda.RDS")

# Branch-associated cell-type proportions
str(df)
table(df$State)
ddff <- as.data.frame(table(df$State, df$cellname_t))
colnames(ddff)[1:2] <- c("State", "cellname_t")
write.csv(ddff, "type_state.csv", row.names = F)
library(ggplot2)
library(dplyr)
ddff$State <- as.character(ddff$State)
cell_order <- c("FCRL3_Naive", "PDE3B_Naive", "TRABD2A_Naive", "TNFRSF1B_Tfh", "GFOD1_Tfh", "TOX2_Tfh")
pal1 <- c("FCRL3_Naive" = "#8DD2C5", "PDE3B_Naive" = "#8DBA8F", "TRABD2A_Naive" = "#00A087",
  "TNFRSF1B_Tfh" = "#B3CDE3", "GFOD1_Tfh" = "#7FB2D5", "TOX2_Tfh" = "#4581B2")
create_dat <- function(state_num) {
  data.frame(
    group = factor(cell_order, levels = cell_order),
    value = ddff[ddff$State == state_num, "Freq"]
  ) %>%
    mutate(
      prop = value / sum(value) * 100,
      ypos = cumsum(prop) - 0.5 * prop
    )
}
dat1 <- create_dat("1")
dat2 <- create_dat("6")
dat3 <- create_dat("4")
p1 <- ggplot(dat1, aes(x = "", y = prop, fill = group)) +
  geom_bar(stat = "identity", width = 0.25, color = "white") +
  coord_polar("y", start = 0) +
  theme_void() + theme(legend.position="none") +
  scale_fill_manual(values = pal1) +
  ggtitle("State1")
p2 <- ggplot(dat2, aes(x = "", y = prop, fill = group)) +
  geom_bar(stat = "identity", width = 0.25, color = "white") +
  coord_polar("y", start = 0) +
  theme_void() + theme(legend.position="none") +
  scale_fill_manual(values = pal1) +
  ggtitle("State6")
p3 <- ggplot(dat3, aes(x = "", y = prop, fill = group)) +
  geom_bar(stat = "identity", width = 0.25, color = "white") +
  coord_polar("y", start = 0) +
  theme_void() + theme(legend.position="none") +
  scale_fill_manual(values = pal1) +
  ggtitle("State4")
pdf("proportion.pdf", width = 9, height = 3)
gridExtra::grid.arrange(p1, p2, p3, ncol = 3)
dev.off()
library(ggridges)
library(RColorBrewer)
library(scales)
pdf("peak.pdf", width = 7, height = 4.5)
ggplot(df, aes(x=Pseudotime, y=cellname_t, fill = stat(x))) +
  geom_density_ridges_gradient(scale=1) +
  geom_vline(xintercept = c(4,8), linetype=2)+
  scale_fill_gradientn(name="Pseudotime",colors = colorRampPalette(rev(brewer.pal(10, "RdYlBu")))(99))+
  scale_y_discrete("")+
  theme_minimal()+
  theme(
    panel.grid = element_blank()
  )
dev.off()

# 10. Map pseudotime to UMAP and spatial positions
objd.t.fli <- AddMetaData(objd.t.fli, metadata = df$Pseudotime, col.name = "Pseudotime")
pdf('CD4T_fli_cellname_umap.pdf',useDingbats=F,width = 8, height = 6)
DimPlot(objd.t.fli, reduction = "umap", label = TRUE, repel = TRUE, cols = pal1,
    group.by = 'cellname_t', pt.size = 0.8)  + theme_global + ggtitle("Cluster")
dev.off()
pdf('CD4T_fli_Pseudotime_umap.pdf',useDingbats=F,width = 8, height = 6)
FeaturePlot(objd.t.fli, reduction = "umap", features = "Pseudotime", raster = FALSE, pt.size = 0.8) +
  theme_global + scale_color_gradientn(colors = rev(brewer.pal(10, "RdYlBu"))) + ggtitle("Pseudotime")
dev.off()
pal2 <- c("#7FB2D5","#8DBA8F","#00A087","#FFD569", "#984EA3","#F4CAE4")
pdf("CD4T_fli_spatial_byPatient_cellname.pdf", width = 8, height = 6)
for (p in levels(objd.t.fli$Patient)) {
  message("Plotting Patient: ", p)
  obj.p <- subset(objd.t.fli, subset = Patient == p)
  p2 <- DimPlot(obj.p, reduction = "spatial", group.by = "cellname_t", pt.size = 1, cols = pal2) + ggtitle(p) + theme_global
  print(p2)
}
dev.off()
pdf("CD4T_fli_spatial_byPatient_Pseudotime.pdf", width = 8, height = 6)
for (p in levels(objd.t.fli$Patient)) {
  message("Plotting Patient: ", p)
  obj.p <- subset(objd.t.fli, subset = Patient == p)
  p2 <- FeaturePlot(obj.p, reduction = "spatial", features = "Pseudotime", pt.size = 1) + ggtitle(p) + theme_global + scale_color_gradientn(colors = rev(brewer.pal(10, "RdYlBu")))
  print(p2)
}
dev.off()

# 11. Pseudotime and branch-dependent expression
Time_diff <- differentialGeneTest(cds[express_genes,], cores = 1, fullModelFormulaStr = "~sm.ns(Pseudotime)")
write.csv(Time_diff,"Time_diff_all.csv", row.names = F)
Time_genes <- Time_diff %>%pull(gene_short_name) %>% as.character()
plot_pseudotime_heatmap(cds[Time_genes,],num_clusters=4, show_rownames=T, return_heatmap=T)
expressed_genes = row.names(subset(fData(cds),num_cells_expressed>=10))
pseudotime_de <- differentialGeneTest(cds[expressed_genes,],  # Expressed genes
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
pseudotime_de <- pseudotime_de[order(pseudotime_de$qval), ]
write.csv(pseudotime_de, "pseudotime_de.csv", row.names = F)
saveRDS(pseudotime_de,"pseudotime_de.RDS")
pseudotime_de2 <- differentialGeneTest(cds[express_genes,],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
pseudotime_de2 <- pseudotime_de2[order(pseudotime_de2$qval), ]
write.csv(pseudotime_de2, "pseudotime_de.var.csv", row.names = F)
saveRDS(pseudotime_de2,"pseudotime_de.var.RDS")

# Pseudotime expression heatmaps
Time_genes <- pseudotime_de2 %>%pull(gene_short_name) %>% as.character()
gene.show <- head(Time_genes, 250)
pdf("Time_heatmap.pdf", width = 6, height = 10)
p1 <- plot_pseudotime_heatmap(cds[gene.show,], num_clusters=5,
                             show_rownames=T, return_heatmap=T,
                             hmcols = colorRampPalette(c("#7990C8", "white", "#F05662"))(100))
print(p1)
dev.off()
p1$tree_row
s.gene <- gene.show
pdf("gene.pse.pdf",width = 6, height = 4)
for (i in s.gene) {
  p <- plot_genes_in_pseudotime(cds[i,], color_by = "cellname_t") + scale_color_manual(values = pal1)
  print(p)
}
dev.off()

# Branch expression analysis
BEAM_res <- BEAM(cds[Time_genes,], branch_point = 2, progenitor_method = "duplicate", cores = 1)
BEAM_res <-BEAM_res[order(BEAM_res$qval),]
BEAM_res <-BEAM_res[,c("gene_short_name", "pval", "qval")]
head(BEAM_res)
write.csv(BEAM_res, "BEAM_res.csv", row.names = F)
saveRDS(BEAM_res, file = "BEAM_res.RDS")
plot_genes_branched_heatmap(cds[row.names(subset(BEAM_res, qval < 1e-4)),], branch_point = 2,
                            num_clusters = 3, cores = 1, use_gene_short_name= T, show_rownames = T)
BEAM_genes <- top_n(BEAM_res, n = 250, desc(qval)) %>% pull(gene_short_name) %>% as.character()
tmp1 = plot_genes_branched_heatmap(cds[BEAM_genes,],
                                 branch_point = 2, num_clusters = 3,
                                 cores = 1, branch_labels = c("Cell fate 1", "Cell fate 2"),
                                 hmcols = colorRampPalette(c("#4974a4", "white", "#b81316"))(66),
                                 branch_colors = c("lightgrey", "#8DBA8F", "#A473AE"),  # Branch colors
                                 use_gene_short_name = T,
                                 show_rownames = F,
                                 return_heatmap = T  # Return heatmap information
)
pdf("branched_heatmap.pdf", width = 6, height = 7)
tmp1$ph_res
dev.off()
# Export pseudotime gene groups.
gene_group=tmp1$annotation_row
gene_group$gene=rownames(gene_group)
write.csv(gene_group, "branched_gene_group.csv")

# Pathway enrichment of trajectory genes

# Trajectory heatmaps
library(jjAnno)
library(ClusterGVis)
library(monocle)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(ComplexHeatmap)
library(clusterProfiler)
setwd("/mnt/data_nvme1/lzh/test/AH/atlas/recluster/CD4T/monocle")
set.seed(1234)
cds <- readRDS("cds.RDS")
pseudotime_de <- readRDS("pseudotime_de.var.RDS")
gene.show <- pseudotime_de %>%pull(gene_short_name) %>% as.character() %>% head(250)
p2 <- plot_pseudotime_heatmap2(cds[gene.show,], num_clusters=5, cores = 1
                               )
str(p2)
gene2 <- read.csv("/mnt/data_nvme1/lzh/test/AH/atlas/recluster/CD4T/sig/top5T.csv")
gene2 <- gene2[gene2$cluster %in% c("FCRL3_Naive", "PDE3B_Naive", "TRABD2A_Naive", "TNFRSF1B_Tfh", "GFOD1_Tfh", "TOX2_Tfh"),] %>% pull(gene)
gene2 <- intersect(gene2, gene.show)
pdf("Time_heatmap.inte2.pdf", width = 6.5, height = 6)
visCluster(object = p2, plot.type = "line")
visCluster(object = p2, plot.type = "both"
           )
visCluster(object = p2 ,plot.type = "heatmap",
           markGenes = gene2
           )
dev.off()
library(org.Hs.eg.db)
enrich.k <- enrichCluster(object = p2, OrgDb = org.Hs.eg.db,
                        type = "KEGG", organism = "hsa",
                        pvalueCutoff = 0.5,
                        topn = 5, seed = 1234)
dim(enrich.k)
enrich.k$Description
enrich.b <- enrichCluster(object = p2, OrgDb = org.Hs.eg.db,
                          type = "BP", organism = "hsa",
                          pvalueCutoff = 0.5,
                          topn = 5, seed = 1234)
dim(enrich.b)
enrich.b$Description
pdf("Time_heatmap.inte3.pdf", width = 12, height = 8)
visCluster(object = p2, plot.type = "both",
           column_names_rot = 45, show_row_dend = F,
           markGenes = gene2, markGenes.side = "left",
           annoTerm.data = enrich.b,
           go.col = rep(jjAnno::useMyCol("calm",n = 5),each = 5),
           cluster.order = c(5:1),  line.side = "left")
dev.off()

# Branch-point heatmaps
library(jjAnno)
library(ClusterGVis)
library(monocle)
library(igraph)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(ComplexHeatmap)
library(clusterProfiler)
setwd("/mnt/data_nvme1/lzh/test/AH/atlas/recluster/CD4T/monocle")
set.seed(1234)
cds <- readRDS("cds.RDS")
BEAM_res <- readRDS("BEAM_res.RDS")
BEAM_genes <- top_n(BEAM_res, n = 250, desc(qval)) %>% pull(gene_short_name) %>% as.character()
p2 <- plot_genes_branched_heatmap2(cds[BEAM_genes,],
                                   branch_point = 2, num_clusters = 3,
                                   cores = 1,
                                   use_gene_short_name = T,
                                   show_rownames = T)
str(p2)
gene2 <- read.csv("/mnt/data_nvme1/lzh/test/AH/atlas/recluster/CD4T/sig/top10T.csv")
gene2 <- gene2[gene2$cluster %in% c("FCRL3_Naive", "PDE3B_Naive", "TRABD2A_Naive", "TNFRSF1B_Tfh", "GFOD1_Tfh", "TOX2_Tfh"),] %>% pull(gene)
gene2 <- intersect(gene2, BEAM_genes)
pdf("Time_heatmap_branched.inte2.pdf", width = 6.5, height = 6)
visCluster(object = p2, plot.type = "line")
visCluster(object = p2, plot.type = "both",pseudotime_col = c("#A473AE", "lightgrey", "#8DBA8F")
           )
visCluster(object = p2 ,plot.type = "heatmap",pseudotime_col = c("#A473AE", "lightgrey", "#8DBA8F"),
           markGenes = gene2
           )
dev.off()
library(org.Hs.eg.db)
enrich.k <- enrichCluster(object = p2, OrgDb = org.Hs.eg.db,
                        type = "KEGG", organism = "hsa",
                        pvalueCutoff = 0.5,
                        topn = 5, seed = 1234)
dim(enrich.k)
enrich.k$Description
enrich.b <- enrichCluster(object = p2, OrgDb = org.Hs.eg.db,
                          type = "BP", organism = "hsa",
                          pvalueCutoff = 0.5,
                          topn = 5, seed = 1234)
dim(enrich.b)
enrich.b$Description
pdf("Time_heatmap_branched.inte3.pdf", width = 12, height = 8)
visCluster(object = p2, plot.type = "both", pseudotime_col = c("#A473AE", "lightgrey", "#8DBA8F"),
           column_names_rot = 45, show_row_dend = F,
           markGenes = gene2, markGenes.side = "left",
           annoTerm.data = enrich.b,
           go.col = rep(jjAnno::useMyCol("calm",n = 3),each = 5),
           cluster.order = c(3:1),  line.side = "left")
dev.off()

# 12. Transcription factor analysis
set.seed(1234)
setwd("/mnt/data_nvme1/lzh/test/AH/atlas/recluster/CD4T/dorothea")
library(Seurat)
library(tidyverse)
library(VennDiagram)
library(grid)
library(futile.logger)
library(decoupleR)
library(OmnipathR)
library(dorothea)
library(viper)

# TFs among TNFRSF1B_Tfh differential genes
all_tfs <- read.table("/data/YQ/linzh/pub_data/ref/scenic/allTFs_hg38.txt")
all_tfs <- all_tfs$V1
sig_gene <- read.csv("/mnt/data_nvme1/lzh/test/AH/atlas/recluster/CD4T/sig/sig_gene.tfh.csv")
sig_gene_tfh <- sig_gene %>% filter(avg_log2FC > 0) %>% pull(X)
TF.intersect <- intersect(all_tfs, sig_gene_tfh)
pdf("tf_deg_venn.pdf", width = 6, height = 6)
p1=venn.diagram(x = list(all_tfs = all_tfs,
                         fm_tfh_degs = sig_gene_tfh),
                col = "transparent", fill=c("#FFE2CE", "#C4D8E9"),
                fontface = "bold",
                cat.col='black' ,label.col ='black', cat.default.pos = "outer",
                resolution = 300,filename = NULL
)
grid.draw(p1)
dev.off()

# TF regulatory network analysis
objd.t.fli <- readRDS("/mnt/data_nvme1/lzh/test/AH/atlas/recluster/CD4T/objd.t.fli.RDS")
objd.tfh <- subset(objd.t.fli, subset = celltype_t == "Tfh")
objd.tfh$cellname_t <- as.character(objd.tfh$cellname_t)
objd.tfh$cellname_t <- factor(objd.tfh$cellname_t, levels = c("TNFRSF1B_Tfh", "GFOD1_Tfh", "TOX2_Tfh"))
saveRDS(objd.tfh, "../objd.tfh.RDS")
dorothea_regulon_human <- get(data("dorothea_hs", package = "dorothea"))
# Use regulons with confidence levels A, B, and C.
regulon <- dorothea_regulon_human %>%
  dplyr::filter(confidence %in% c("A","B","C"))
curl::curl_fetch_memory("https://www.ensembl.org/info/about/species.html")
net <- decoupleR::get_collectri(split_complexes = T,organism='human' )
head(net)
regulon
save(list = c("net","regulon"), file = "dorothea_human.RData")
library(Matrix)
expr <- GetAssayData(objd.tfh, assay = "Spatial", slot = "data")
genes_use <- intersect(rownames(expr), unique(regulon2$target))
expr_use <- expr[genes_use, , drop = FALSE]
library(dplyr)
regulon2 <- regulon %>%
  mutate(
    mor = if ("mor" %in% names(.)) as.numeric(mor) else 1,
    likelihood = if ("likelihood" %in% names(.)) as.numeric(likelihood) else 1
  ) %>%
  rename(
    source = any_of(c("source", "regulator", "tf", "TF")),
    target = any_of(c("target", "gene", "target_gene"))
  ) %>%
  select(source, target, mor, likelihood) %>%
  filter(!is.na(source), !is.na(target))
dim(expr_use)
objd.tfh <- ScaleData(objd.tfh, assay="Spatial", features=genes_use, verbose=FALSE)
expr_use2 <- GetAssayData(objd.tfh, assay="Spatial", slot="scale.data")
viper_scores <- decoupleR::run_viper(
  mat = expr_use2,
  network = regulon2,
  minsize = 4,
  pleiotropy = FALSE,
  eset.filter = FALSE,
  verbose = FALSE
)
saveRDS(viper_scores, "viper_scores.RDS")

# Integrate regulon activities
viper_mat <- viper_scores %>%
  filter(statistic == "viper") %>%
  select(source, condition, score) %>%
  pivot_wider(names_from = condition, values_from = score) %>%
  tibble::column_to_rownames("source") %>%
  as.matrix()
cells <- colnames(objd.tfh)
viper_mat <- viper_mat[, intersect(colnames(viper_mat), cells), drop = FALSE]
viper_mat <- viper_mat[, cells, drop = FALSE]  # Match Seurat cell order
saveRDS(viper_mat, "viper_mat.RDS")
objd.tfh[["VIPER"]] <- CreateAssayObject(data = viper_mat)
DefaultAssay(objd.tfh) <- "VIPER"
DefaultAssay(objd.tfh)
dim(GetAssayData(objd.tfh, slot = "data"))
head(rownames(objd.tfh[["VIPER"]]))
objd.tfh@assays$VIPER[1:4,1:4]

# Visualize TF activity
library(dplyr)
library(tidyr)
library(tibble)
library(pheatmap)
top_n_per_group <- 20  # Top differential TFs per subtype
oversample <- 3
use_zscore <- TRUE  # Use z-score normalization
avg_mat <- summarized_viper_scores %>%
  dplyr::select(tf, cell_type, avg) %>%
  tidyr::pivot_wider(names_from = cell_type, values_from = avg) %>%
  tibble::column_to_rownames("tf") %>%
  as.matrix()
diff_all <- summarized_viper_scores %>%
  dplyr::select(tf, cell_type, avg) %>%
  dplyr::group_by(tf) %>%
  dplyr::mutate(
    other_mean = (sum(avg) - avg) / (dplyr::n() - 1),
    delta = avg - other_mean
  ) %>%
  dplyr::ungroup()
diff_top <- diff_all %>%
  dplyr::group_by(cell_type) %>%
  dplyr::slice_max(delta, n = top_n_per_group * oversample, with_ties = FALSE) %>%
  dplyr::ungroup()
tf_group <- diff_top %>%
  dplyr::group_by(tf) %>%
  dplyr::slice_max(delta, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(cell_type) %>%
  dplyr::slice_max(delta, n = top_n_per_group, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(cell_type, dplyr::desc(delta))
selected_tfs <- tf_group$tf
heat_mat <- avg_mat[selected_tfs, , drop = FALSE]
if (use_zscore) {
  heat_mat <- t(scale(t(heat_mat)))
}
row_anno <- data.frame(Subcluster = tf_group$cell_type)
rownames(row_anno) <- tf_group$tf
col_anno <- data.frame(Subcluster = colnames(heat_mat))
rownames(col_anno) <- colnames(heat_mat)
palette_length <- 100
my_color <- colorRampPalette(c("#4581B2",'white',"#CB5C5B"))(palette_length)
my_breaks <- c(
  seq(min(heat_mat, na.rm = TRUE), 0,
      length.out = ceiling(palette_length / 2) + 1),
  seq(max(heat_mat, na.rm = TRUE) / palette_length,
      max(heat_mat, na.rm = TRUE),
      length.out = floor(palette_length / 2))
)
pdf("VIPER_TF_activity_per_subcluster.pdf", width = 8, height = 10)
pheatmap(
  mat = heat_mat,
  color = my_color,
  breaks = my_breaks,
  cluster_cols = FALSE,
  fontsize = 10, fontsize_row = 10,
  cellheight = 10, cellwidth = 50,
  angle_col = 45,
  border_color = "white", border = 0.5,
  treeheight_col = 0
)
dev.off()

# 13. pySCENIC analysis

# Prepare expression matrices
rm(list = ls())
setwd("/mnt/data_nvme1/lzh/test/AH/atlas/recluster/CD4T/pyscenic")
set.seed(1234)
devtools::install_github("aertslab/SCopeLoomR")
BiocManager::install("RcisTarget")
devtools::install_github("aertslab/SCENIC")
suppressMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(biomaRt)
  library(SCopeLoomR)
  library(RcisTarget)
  library(GENIE3)
  library(AUCell)
  library(SCENIC)
  library(pheatmap)
  library(circlize)
  library(gridExtra)
})

# Load pySCENIC results
objd.tfh <- readRDS("/mnt/data_nvme1/lzh/test/AH/atlas/recluster/CD4T/objd.tfh.RDS")
counts <- as.matrix(objd.tfh@assays$RNA@layers$counts)
colnames(counts) <- colnames(objd.tfh)
rownames(counts) <- rownames(objd.tfh)
write.csv(counts,file = "/home/cwj/project/09_pyscenic/result/counts.csv")

# Plot regulon activity
library(Seurat)
library(dplyr)
library(limma)
library(Matrix)
library(magrittr)
library(stringr)
library(patchwork)
library(clustree)
library(paletteer)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(forcats)
set.seed(1234)
theme_global <- theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank()) + theme(aspect.ratio = 1)
filedir="/mnt/data_nvme1/lzh/test/AH/atlas/recluster/CD4T/pyscenic"
seuratObj <- readRDS("/mnt/data_nvme1/lzh/test/AH/atlas/recluster/CD4T/objd.tfh.RDS")
aucmtx = readRDS(file.path(filedir, "objd.tfh.auc_mtx.RDS"))
aucmtx = aucmtx[colnames(seuratObj), ]
genes <- intersect(colnames(aucmtx), rownames(seuratObj))
Clusters = unique(as.character(seuratObj$cellname_t))

# Bar plots
rssMat = read.csv(file.path(filedir, "objd.tfh.RSS.csv"))
rownames(rssMat) = rssMat$X
rssMat$X = NULL
pdf("tfh_rssMat_tfs_barplot.pdf", useDingbats = F, height = 5, width = 8)
for (i in colnames(rssMat)) {
  tmp = data.frame(RSS = rssMat[[i]],
                   TF = rownames(rssMat) ) %>% arrange(desc(RSS)) %>% head(25)
  tmp$TF = factor(tmp$TF, levels=tmp$TF)
  p = ggplot(data = tmp, aes(x = TF, y = RSS)) +
    geom_bar(stat="identity", fill='#f2c396', alpha=.8) +  xlab("") +
    theme_classic() + ggtitle(i) +
    theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))
  print(p)
}
dev.off()

# Scatter plots
plotlist <- lapply(Clusters, function(x) {
  library(dplyr)
  library(ggrepel)
  library(ggthemes)
  library(gridExtra)
  dat <- data.frame(
    rss = as.numeric(rssMat[, x]),
    tfs = rownames(rssMat)
  )
  dat1 <- dat %>% arrange(desc(rss))
  dat1$index <- 1:nrow(dat1)
  dat1$showtfs <- ifelse(dat1$index <= 25, "1", "0")
  dat2 <- dat1[1:25,]
  library(ggplot2)
  p <- ggplot(data = dat1[1:50,], aes(x = index, y = rss, label=tfs)) +
    geom_point(aes(color = showtfs)) + theme_classic() +
    scale_color_manual(values = c("grey", "red")) +
    geom_point(data = dat2, alpha = 1, size = 1.8, shape = 1, stroke = 1, color = "black") +
    geom_text_repel(data = dat2, show.legend = FALSE,
                    size = 6, box.padding = unit(0.4, "lines"), point.padding = unit(0.3, "lines")) +
    guides(color = guide_legend(title = NULL)) +
    xlab("tfs") + ylab("RSS") + ggtitle(x) +
    theme(legend.position = "none") +
    theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
    theme(aspect.ratio = 1)
  return(p)
})
pdf("tfh_scVI_rssMat_tfs.pdf", useDingbats = F, height = 12, width = 10)
for (i in 1:3) {
  print(plotlist[i])
}
dev.off()
tfs = readRDS(file.path(filedir, "objd.tfh_info_clusters.RDS"))
subtype_list <- sapply(Clusters, function(x) {
  df = tfs[tfs[[str_glue(x, "(sc_act_pvalue)")]] < 1e-30 &
             tfs[[str_glue(x, "(sc_exp_pvalue)")]] < 0.005 &
             tfs[[str_glue(x, "(sc_act_mean)")]] > 0.04 &
             tfs[[str_glue(x, "(rss)")]] > 0.1,]
  df$id
})

# Heatmaps
pdf("tfh_scVI_tf_significant_heatmap.pdf", width = 8, height = 10)
tmp = as.data.frame(tfs[, str_detect(colnames(tfs), "sc_act_mean")])
rownames(tmp) = tfs$id
colnames(tmp) = str_replace_all(colnames(tmp), "\\(sc_act_mean\\)", "")
pheatmap::pheatmap(tmp[unique(unlist(subtype_list)), ], scale = "row",
                   fontsize = 12, cellwidth = 36, cellheight = 12,
                   treeheight_row = 5, cluster_cols = FALSE,
                   clustering_method="ward.D2", border = NA,
                   color = colorRampPalette(c("#1272B4", "white", "#C5533F"))(100))
dev.off()
pdf("tfh_scVI_tf_significant_heatmap_TNFRSF1B_Tfh.pdf", width = 8, height = 10)
tmp = as.data.frame(tfs[, str_detect(colnames(tfs), "sc_act_mean")])
rownames(tmp) = tfs$id
colnames(tmp) = str_replace_all(colnames(tmp), "\\(sc_act_mean\\)", "")
pheatmap::pheatmap(tmp[unique(unlist(subtype_list$TNFRSF1B_Tfh)), ], scale = "row",
                   fontsize = 12, cellwidth = 20, cellheight = 15,
                   treeheight_row = 5, cluster_cols = FALSE,
                   clustering_method="ward.D2", border = NA,
                   color = colorRampPalette(c("#1272B4", "white", "#C5533F"))(100))
dev.off()

# Subtype-specific transcription factors
tfs2 <- readRDS("/mnt/data_nvme1/lzh/test/AH/atlas/recluster/CD4T/pyscenic/objd.tfh_info_clusters.RDS") %>% as.data.frame
rownames(tfs2) = tfs2$id
Clusters <- c("TNFRSF1B_Tfh","GFOD1_Tfh","TOX2_Tfh")
subtype_list2 <- sapply(Clusters, function(x) {
    df = tfs2[tfs2[[str_glue(x, "(sc_act_pvalue)")]] < 1e-30 &
        tfs2[[str_glue(x, "(sc_exp_pvalue)")]] < 0.005 &
        tfs2[[str_glue(x, "(sc_act_mean)")]] > 0.04 &
        tfs2[[str_glue(x, "(rss)")]] > 0.1,]
    rownames(df)
})
TNFRSF1B_TF <- subtype_list2$TNFRSF1B_Tfh
write.csv(TNFRSF1B_TF, "TNFRSF1B_TF.csv",row.names = FALSE)

# TF activity by subcluster
tf.act <- readRDS("/mnt/data_nvme1/lzh/test/AH/atlas/recluster/CD4T/pyscenic/objd.tfh.auc_mtx.RDS")
tf.act1 <- tf.act[,TNFRSF1B_TF]
objd.tfh <- readRDS("/mnt/data_nvme1/lzh/test/AH/atlas/recluster/CD4T/objd.tfh.RDS")
tf.act1 <- tf.act1[colnames(objd.tfh),]
colnames(tf.act1) <- paste0("Activity_", colnames(tf.act1))
objd.tfh <- AddMetaData(objd.tfh, tf.act1)
str(objd.tfh@meta.data)

# UMAP visualizations
pdf("activity.umap.pdf", width = 5, height = 5)
for (i in colnames(tf.act1)) {
  try({
    p = FeaturePlot(objd.tfh, features = i, pt.size = 0.8,
                    raster = F, reduction="umap") + theme_global +
                    scale_color_distiller(palette = "RdYlBu")
    print(p)
  },silent = T)
}
dev.off()
pdf("exp.umap.pdf", width = 5, height = 5)
for (i in TNFRSF1B_TF) {
  try({
    p = FeaturePlot(objd.tfh, features = i, pt.size = 0.8,
                    raster = F, reduction="umap") + theme_global +
                    scale_color_distiller(palette = "RdYlBu")
    print(p)
  },silent = T)
}
dev.off()

# Regulon distributions
pal <-c("#B3CDE3","#7FB2D5","#4581B2")
pdf('cluster.umap.pdf', useDingbats=F,width = 5, height = 5)
DimPlot(objd.tfh, reduction = "umap", label = T, repel = TRUE, cols = pal, raster = FALSE,
    group.by = 'cellname_t', pt.size = 0.8)  + theme_global + ggtitle("Tfh_Clustrer")
dev.off()

# Dot plots
pdf("activity.dot.pdf", width = 6, height = 5)
DotPlot(objd.tfh, assay = "Spatial", features = colnames(tf.act1),
    dot.min = 0, dot.scale = 6, idents = NULL, group.by = "cellname_t",
        split.by = NULL, cluster.idents = FALSE, scale = TRUE, scale.by = "radius", scale.min = NA, scale.max = NA) +
    scale_color_distiller(palette = "RdYlBu") +  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_y_discrete(limits = levels(objd.tfh$cellname_t)) + labs(x="", y="") +
  theme(panel.border = element_rect(color = "black", size = 0.5, fill = NA)) + coord_flip()
dev.off()
pdf("exp.dot.pdf", width = 5, height = 4.5)
DotPlot(objd.tfh, assay = "Spatial", features = TNFRSF1B_TF,
    dot.min = 0, dot.scale = 6, idents = NULL, group.by = "cellname_t",
        split.by = NULL, cluster.idents = FALSE, scale = TRUE, scale.by = "radius", scale.min = NA, scale.max = NA) +
    scale_color_distiller(palette = "RdYlBu") +  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_y_discrete(limits = levels(objd.tfh$cellname_t)) + labs(x="", y="") +
  theme(panel.border = element_rect(color = "black", size = 0.5, fill = NA)) + coord_flip()
dev.off()

# 14. Donor-level subcluster proportion comparison (AH vs ATH-A)
# Each donor is one independent replicate; absent subclusters count as 0%.
library(ggplot2)

sample_rds <- "/mnt/data_nvme1/lzh/test/AH/atlas/recluster/CD4T/CD4T.all.RDS"
sample_obj <- readRDS(sample_rds)
sample_meta <- sample_obj@meta.data[, c("Patient", "TH_meta", "cellname_t"), drop = FALSE]
colnames(sample_meta) <- c("Patient", "GroupCode", "Subtype")
sample_meta[] <- lapply(sample_meta, as.character)
sample_meta <- sample_meta[complete.cases(sample_meta), , drop = FALSE]
sample_meta <- subset(sample_meta, GroupCode %in% c("N", "Y"))
if (!nrow(sample_meta)) stop("No cells with valid Patient, TH_meta, and subcluster labels.")
if (!all(c("N", "Y") %in% unique(sample_meta$GroupCode))) {
  stop("Both N (AH) and Y (ATH-A) donors are required.")
}
donor_groups <- unique(sample_meta[, c("Patient", "GroupCode")])
if (anyDuplicated(donor_groups$Patient)) {
  stop("A Patient appears in both groups; check TH_meta before donor-level testing.")
}

subtype_levels <- levels(sample_obj$cellname_t)
if (is.null(subtype_levels)) subtype_levels <- unique(sample_meta$Subtype)
subtype_levels <- subtype_levels[subtype_levels %in% sample_meta$Subtype]
donor_levels <- sort(unique(sample_meta$Patient))
counts <- table(
  Patient = factor(sample_meta$Patient, levels = donor_levels),
  Subtype = factor(sample_meta$Subtype, levels = subtype_levels)
)
donor_tbl <- as.data.frame(counts, stringsAsFactors = FALSE)
names(donor_tbl)[3] <- "CellCount"
donor_tbl$TotalCells <- rowSums(counts)[match(donor_tbl$Patient, rownames(counts))]
donor_tbl$Proportion <- 100 * donor_tbl$CellCount / donor_tbl$TotalCells
donor_tbl$GroupCode <- donor_groups$GroupCode[
  match(donor_tbl$Patient, donor_groups$Patient)
]
donor_tbl$Group <- factor(donor_tbl$GroupCode, levels = c("N", "Y"),
                          labels = c("AH", "ATH-A"))
donor_tbl$Subtype <- factor(donor_tbl$Subtype, levels = subtype_levels)
donor_tbl <- donor_tbl[order(donor_tbl$Subtype, donor_tbl$Group,
                             donor_tbl$Patient), ]

group_summary <- aggregate(Proportion ~ Subtype + Group, data = donor_tbl,
  FUN = function(x) c(mean = mean(x), sem = sd(x) / sqrt(length(x))))
summary_tbl <- data.frame(
  Subtype = group_summary$Subtype,
  Group = group_summary$Group,
  Mean = group_summary$Proportion[, "mean"],
  SEM = group_summary$Proportion[, "sem"]
)

test_tbl <- do.call(rbind, lapply(subtype_levels, function(ct) {
  tmp <- donor_tbl[donor_tbl$Subtype == ct, ]
  x <- tmp$Proportion[tmp$Group == "AH"]
  y <- tmp$Proportion[tmp$Group == "ATH-A"]
  p_value <- if (length(x) >= 2 && length(y) >= 2) {
    tryCatch(t.test(x, y, paired = FALSE, var.equal = FALSE)$p.value,
             error = function(e) NA_real_)
  } else NA_real_
  data.frame(
    Subtype = ct, N_AH = length(x), N_ATHA = length(y),
    Mean_AH = mean(x), Mean_ATHA = mean(y), P_value = p_value
  )
}))
test_tbl$Q_BH <- p.adjust(test_tbl$P_value, method = "BH")
test_tbl$Label <- ifelse(
  is.na(test_tbl$P_value), "P = NA",
  paste0("P = ", sprintf("%.3g", test_tbl$P_value))
)
test_tbl$LabelY <- vapply(test_tbl$Subtype, function(ct) {
  max(donor_tbl$Proportion[donor_tbl$Subtype == ct]) * 1.15 + 0.5
}, numeric(1))
test_tbl$Subtype <- factor(test_tbl$Subtype, levels = subtype_levels)

sample_output <- file.path(dirname(sample_rds), "sig")
dir.create(sample_output, recursive = TRUE, showWarnings = FALSE)
write.csv(donor_tbl, file.path(sample_output,
          "CD4T_donor_level_subcluster_proportions.csv"), row.names = FALSE)
write.csv(test_tbl[, c("Subtype", "N_AH", "N_ATHA", "Mean_AH",
                       "Mean_ATHA", "P_value", "Q_BH")],
          file.path(sample_output, "CD4T_donor_level_subcluster_ttests.csv"),
          row.names = FALSE)

group_colors <- c("AH" = "#9EC4BE", "ATH-A" = "#FAC795")
plot_subtypes <- function(data_points, data_means, data_tests, columns = 4) {
  ggplot(data_means, aes(x = Group, y = Mean, fill = Group)) +
    geom_col(width = 0.6, color = "#52616B", linewidth = 0.25) +
    geom_errorbar(aes(ymin = pmax(0, Mean - SEM), ymax = Mean + SEM),
                  width = 0.15, na.rm = TRUE) +
    geom_point(data = data_points, aes(y = Proportion),
               position = position_jitter(width = 0.07, height = 0, seed = 1234),
               shape = 21, size = 2.4, stroke = 0.35, color = "#334155") +
    geom_text(data = data_tests,
              aes(x = 1.5, y = LabelY, label = Label),
              inherit.aes = FALSE, size = 3.2) +
    facet_wrap(~ Subtype, scales = "free_y", ncol = columns) +
    scale_fill_manual(values = group_colors, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
    labs(title = "CD4+ T subcluster proportions by donor",
         x = NULL, y = "Cells within lineage (%)", fill = NULL) +
    theme_classic(base_size = 11) +
    theme(legend.position = "top", strip.text = element_text(face = "bold"),
          axis.text.x = element_text(size = 9))
}

pdf(file.path(sample_output, "CD4T_donor_level_subcluster_comparison.pdf"),
    width = 3.5, height = 4, useDingbats = FALSE)
print(plot_subtypes(donor_tbl, summary_tbl, test_tbl))
for (ct in subtype_levels) {
  print(plot_subtypes(
    donor_tbl[donor_tbl$Subtype == ct, ],
    summary_tbl[summary_tbl$Subtype == ct, ],
    test_tbl[test_tbl$Subtype == ct, ],
    columns = 1
  ) + labs(title = paste("CD4+ T:", ct)))
}
dev.off()
