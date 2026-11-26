# Spatial transcriptomics: B-cell lineage subclustering and downstream analysis
.libPaths("/home/lzh/miniconda3/envs/singlecell/lib/R/library")
setwd("/mnt/data_nvme1/lzh/test/AH/atlas/recluster/B")
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
spatial.all <- readRDS("/data/YQ/linzh/test/AH/AH_all_spatial.RDS")
str(spatial.all@meta.data)
table(spatial.all$celltype)

# 01. Subset B-cell lineage populations
objd.b <- subset(spatial.all, celltype == "GC_B"| celltype == "Memory_B"| celltype == "Naive_B"| celltype == "Plasma_cell")
DefaultAssay(objd.b) <- "Spatial"

# 02. Normalize and recluster B-lineage cells
objd.b <- NormalizeData(objd.b, normalization.method = "LogNormalize", scale.factor = 10000) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>%
  ScaleData()
objd.b <- RunPCA(objd.b, npcs = 30)
pdf("pre/ElbowPlot.pdf", width = 8, height = 5)
ElbowPlot(objd.b, ndims = 30)
dev.off()
pc.select <- 1:20
objd.b@meta.data <- objd.b@meta.data %>% select(-starts_with("Spatial_snn_res"))

# Harmony integration
library(harmony)
pca_embed <- Embeddings(objd.b, reduction = "pca")[, pc.select, drop = FALSE]
h_embed <- harmony::HarmonyMatrix(data_mat = pca_embed,
  meta_data = objd.b@meta.data, vars_use = "Patient",  # Batch correction by patient
  do_pca = FALSE, verbose = TRUE)
objd.b[["harmony"]] <- CreateDimReducObject(embeddings = h_embed, key = "harmony_", assay = DefaultAssay(objd.b))
objd.b <- FindNeighbors(objd.b, reduction = "harmony", dims = pc.select)
objd.b <- FindClusters(objd.b, resolution = seq(0.1, 0.8, 0.1))
objd.b <- RunUMAP(objd.b, reduction = "harmony", dims = pc.select)
str(objd.b@meta.data)
table(objd.b$Spatial_snn_res.0.5)
table(objd.b$Spatial_snn_res.0.6)
table(objd.b$Spatial_snn_res.0.7)
table(objd.b$Spatial_snn_res.0.8)
# Confirm B-lineage identity, including plasma-cell markers.

# Check B-cell lineage identity
DefaultAssay(objd.b)<-"Spatial"
pdf("pre/identity.pdf", width = 6, height = 6)
for (i in c("MS4A1","CD79A","CD79B", "CD74", "CD19","MZB1","XBP1","IGHM","IGHD","IGKC","JCHAIN")) {
  try({
    p = FeaturePlot(objd.b, features = i, pt.size = 0.8, cols = c("grey", "red"),
                    raster = F, reduction="umap") + theme_global
    print(p)
  },silent = T)
}
dev.off()

# Set cluster order
library(dplyr)
objd.b$seurat_clusters <- objd.b$Spatial_snn_res.0.6
objd.b$seurat_clusters <- as.character(objd.b$seurat_clusters)
objd.b$seurat_clusters <- factor(objd.b$seurat_clusters, levels = c("0","1","2","3","4","5","6","7","8","9","10","11","12","13"))
objd.b$seurat_clusters <- as.character(objd.b$seurat_clusters)
objd.b$seurat_clusters <- paste("C", objd.b$seurat_clusters, sep ="")
objd.b$seurat_clusters <- factor(objd.b$seurat_clusters, levels = paste0("C", c(0:13)))
levels <- paste0("C", c(0:13))
objd.b$seurat_clusters <- factor(objd.b$seurat_clusters, levels = levels)
levels(objd.b$seurat_clusters)
table(objd.b$seurat_clusters)
pdf('pre/cluster_umap.pdf', useDingbats=F, width = 6, height = 6)
DimPlot(objd.b, reduction = "umap", label = TRUE, repel = TRUE, raster = FALSE,
    group.by = 'seurat_clusters', pt.size = 0.6)  + theme_global + ggtitle("Cluster")
dev.off()

# 03. B-lineage subtype annotation
markers<-c("MS4A1","TCL1A", "IGHD", "IGHM", "IL4R", "FCER2", "CCR7",  # Naive B
           "TRAF4","NME1","HMGB1","CD69", "JUN", "FOS","CD86", "NFKBIA",  # Activated B / pre-GC
           "MEF2B", "RGS13", "BCL6", "AICDA", "S1PR2", "SERPINA9",  # GC B
           "MKI67","TOP2A","HMGB2","TYMS","PCNA", "CXCR4","MYBL1",  # DZ GCB
           "CD83","MIR155HG","BATF","IRF4", "CXCR5","LMO2",  # LZ GCB
           "CD27","TNFRSF13B","FCRL4","AIM2", "TBX21","ITGAX","FCRL5","FCGR2B","CD72",  # Memory B
           "JCHAIN","MZB1","XBP1","HSPA5","SEC11C","FKBP11",  # Plasmablast
           "PRDM1","SDC1","TNFRSF17","DERL3","IGKC","IGHA1","IGHA2",  # Plasma cell
           "TRAC","CD3D","CD3E"
)
Idents(objd.b) <- "seurat_clusters"
pdf('pre/markerBubble2.pdf',useDingbats=F, width = 16, height = 8)
DotPlot(objd.b, features = markers, assay = 'Spatial') + scale_color_distiller(palette = "RdYlBu") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
dev.off()
pdf("pre/marker.umap.pdf", width = 6, height = 6)
for (i in markers) {
  try({
    p = FeaturePlot(objd.b, features = i, pt.size = 0.8, raster = F, reduction="umap") +
                    theme_global + scale_color_distiller(palette = "RdYlBu")
    print(p)
  },silent = T)
}
dev.off()
cluster2celltype_b <- c("C0"="Naive_B", "C1"="Memory_B",
                        "C2"="Memory_B", "C3"="Naive_B",
                        "C4"="GCB_DZ", "C5"="GCB_LZ",
                        "C6"="Memory_B", "C7"="Plasma_cell",
                        "C8"="Pre_GCB", "C9"="GCB_DZ",
                        "C10"="GCB_LZ", "C11"="Mix",
                        "C12"="Naive_B", "C13"="Naive_B")
objd.b[['celltype_b']]<-unname(cluster2celltype_b[objd.b@meta.data$seurat_clusters])
objd.b <- subset(objd.b, celltype_b != "Mix")
dim(objd.b)
table(objd.b$celltype_b)
objd.b$celltype_b <- as.character(objd.b$celltype_b)
objd.b$celltype_b <- factor(objd.b$celltype_b, levels = c("GCB_DZ","GCB_LZ","Memory_B", "Naive_B", "Plasma_cell", "Pre_GCB"))
pal1 <- c("#984EA3","#BC9DCC", "#8DBA8F","#7FB2D5", "#F5B47F","#BEDADA")
pdf('pre/celltype_umap.pdf',useDingbats=F,width = 6, height = 6)
DimPlot(objd.b, reduction = "umap", label = TRUE, repel = TRUE, cols = pal1, raster = FALSE,
    group.by = 'celltype_b', pt.size = 0.6)  + theme_global + ggtitle("Celltype")
dev.off()

# Visualize subtype markers
library(ggplot2)
library(grid)
markers.show<-list(
  GCB_DZ = c("BCL6","RGS13","MKI67","TOP2A"),
  GCB_LZ = c("LMO2","CXCR5"),
  Memory_B = c("TNFRSF13B","FCRL4","AIM2"),
  Naive_B = c("TCL1A","IGHD","IGHM"),
  Plasma_cell = c("MZB1","XBP1","PRDM1"),
  Pre_GCB = c("TRAF4","NME1","NFKBIA")
)
markers.show2 <- c("BCL6","RGS13","MKI67","TOP2A", "LMO2","CXCR5",
                   "TNFRSF13B","FCRL4","AIM2", "TCL1A","IGHD","IGHM",
                   "MZB1","XBP1","PRDM1","TRAF4","NME1","NFKBIA")
Idents(objd.b) <- "celltype_b"
pdf("pre/DotPlot.Ttype.pdf",width = 8, height = 4)
DotPlot(objd.b, features = markers.show2,  cols = c("grey", "red"),
    dot.min = 0, dot.scale = 6, idents = NULL, group.by = "celltype_b",
        split.by = NULL, cluster.idents = FALSE, scale = TRUE, scale.by = "radius", scale.min = NA, scale.max = NA)+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_y_discrete(limits = levels(objd.b$celltype_b)) +
  theme(panel.border = element_rect(color="black"),
        panel.spacing = unit(1, "mm"),
        strip.text = element_text(margin=margin(b=3, unit="mm")),
        strip.placement = 'outlet',
        axis.line = element_blank()) + labs(x="", y="")
DotPlot(objd.b, features = markers.show2,
    dot.min = 0, dot.scale = 6, idents = NULL, group.by = "celltype_b",
        split.by = NULL, cluster.idents = FALSE, scale = TRUE, scale.by = "radius", scale.min = NA, scale.max = NA) +
    scale_color_distiller(palette = "RdYlBu") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_y_discrete(limits = levels(objd.b$celltype_b)) + labs(x="", y="") +
  theme(panel.border = element_rect(color = "black", size = 0.5, fill = NA))
dev.off()

# 04. Differential markers by cluster
levels <- paste0("C", c(0:10,12:13))
objd.b$seurat_clusters <- factor(objd.b$seurat_clusters, levels = levels)
table(objd.b$seurat_clusters)
DefaultAssay(objd.b) <- "Spatial"
Idents(objd.b) <- "seurat_clusters"
logFCfilter=0.50
adjPvalFilter=1e-5
mrksB <- FindAllMarkers(objd.b, logfc.threshold = 0.25, min.pct = 0.2, only.pos = TRUE)
length(mrksB$gene)
sig.markers<-mrksB[(abs(as.numeric(mrksB$avg_log2FC))>logFCfilter & as.numeric(mrksB$p_val_adj)<adjPvalFilter),]
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
Idents(objd.b) <- "seurat_clusters"
pdf("sig/qc.pdf",width = 15, height = 5)
VlnPlot(objd.b, features=c("nFeature_Spatial","nCount_Spatial","percent.mt"),
        pt.size=0, ncol=3)
dev.off()

# Assign named subclusters
cluster2cellname <- c("C0"="IGHD_Naive_B", "C1"="PDE4D_Memory_B",
                      "C2"="PTPRJ_Memory_B", "C3"="STEAP1B_Naive_B",
                      "C4"="MKI67_GCB_DZ", "C5"="MAML3_GCB_LZ",
                      "C6"="FCRL4_Memory_B", "C7"="IGHG1_Plasma_cell",
                      "C8"="CD83_Pre_GCB", "C9"="DIAPH3_GCB_DZ",
                      "C10"="BACH2_GCB_LZ", "C12"="ACSM3_Naive_B",
                      "C13"="STAT4_Naive_B")
objd.b[['cellname_b']]<-unname(cluster2cellname[objd.b@meta.data$seurat_clusters])
objd.b$cellname_b <- factor(objd.b$cellname_b, levels = c("DIAPH3_GCB_DZ", "MKI67_GCB_DZ", "BACH2_GCB_LZ", "MAML3_GCB_LZ",
                            "FCRL4_Memory_B", "PDE4D_Memory_B", "PTPRJ_Memory_B", "ACSM3_Naive_B","IGHD_Naive_B",
                            "STAT4_Naive_B", "STEAP1B_Naive_B", "IGHG1_Plasma_cell","CD83_Pre_GCB"))
table(objd.b$cellname_b)

# Visualize annotated clusters
Idents(objd.b) <- "cellname_b"
pal <-c("#FD8052","#F5B47F","#FFF3CA","#FFE2CE","#8DD2C5","#8DBA8F","#00A087","#4581B2","#B3CDE3","#354e97","#7FB2D5","#F4CAE4","#984EA3")
pdf('sig/cellname_umap.pdf',useDingbats=F,width = 8, height = 6)
DimPlot(objd.b, reduction = "umap", label = TRUE, repel = TRUE, cols = pal, raster = FALSE,
    group.by = 'cellname_b', pt.size = 0.8)  + theme_global + ggtitle("Cluster")
dev.off()

# Spatial distribution by patient
pdf("sig/B_spatial_byPatient_cellname.pdf", width = 8, height = 6)
for (p in levels(objd.b$Patient)) {
  message("Plotting Patient: ", p)
  obj.p <- subset(objd.b, subset = Patient == p)
  p2 <- DimPlot(obj.p, reduction = "spatial", group.by = "cellname_b", pt.size = 1, cols = pal) + ggtitle(p) + theme_global
  print(p2)
}
dev.off()
pdf("sig/B_spatial_byPatient_cellname2.pdf", width = 6, height = 6)
for (p in levels(objd.b$Patient)) {
  message("Plotting Patient: ", p)
  obj.p <- subset(objd.b, subset = Patient == p)
  p2 <- DimPlot(obj.p, reduction = "spatial", group.by = "cellname_b", pt.size = 1, cols = pal) + ggtitle(p) + theme_global+ NoLegend()
  print(p2)
}
dev.off()

# 05. Differential markers by named subcluster
DefaultAssay(objd.b) <- "Spatial"
Idents(objd.b) <- "cellname_b"
logFCfilter=0.50
adjPvalFilter=1e-5
mrksB <- FindAllMarkers(objd.b, logfc.threshold = 0.25, min.pct = 0.2, only.pos = TRUE)
length(mrksB$gene)
sig.markers<-mrksB[(abs(as.numeric(mrksB$avg_log2FC))>logFCfilter & as.numeric(mrksB$p_val_adj)<adjPvalFilter),]
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
top50B<-sig.markers2 %>% group_by(cluster) %>% top_n(n=50, wt=avg_log2FC)
top30B<-sig.markers2 %>% group_by(cluster) %>% top_n(n=30, wt=avg_log2FC)
top20B<-sig.markers2 %>% group_by(cluster) %>% top_n(n=20, wt=avg_log2FC)
top10B<-sig.markers2 %>% group_by(cluster) %>% top_n(n=10, wt=avg_log2FC)
top5B<-sig.markers2 %>% group_by(cluster) %>% top_n(n=5, wt=avg_log2FC)
write.csv(top50B, "sig/top50B.csv", row.names = F)
write.csv(top30B, "sig/top30B.csv", row.names = F)
write.csv(top20B, "sig/top20B.csv", row.names = F)
write.csv(top10B, "sig/top10B.csv", row.names = F)
write.csv(top5B, "sig/top5B.csv", row.names = F)

# Marker heatmaps
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
objd.b <- ScaleData(objd.b, verbose = FALSE)
mapal <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256)
pdf(file="sig/top10HeatmapB.pdf", width=16, height=14)
DoHeatmap(objd.b, features = as.character(unique(top10B$gene)),
        slot = "scale.data", raster = F, size = 3,
        group.by = "cellname_b", group.colors = pal) +
        scale_fill_gradient2(low = '#0099CC',mid = 'white',high = '#CC0033',
                       name = 'Z-score')
dev.off()
pdf(file="sig/top5HeatmapB.pdf", width=12, height=10.5)
DoHeatmap(objd.b, features = as.character(unique(top5B$gene)),
        slot = "scale.data", raster = F, size = 3,
        group.by = "cellname_b", group.colors = pal) +
        scale_fill_gradient2(low = '#0099CC',mid = 'white',high = '#CC0033',
                       name = 'Z-score')
dev.off()

# Cluster-averaged expression
mapal <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256)
cluster.averages <- AverageExpression(objd.b,return.seurat = TRUE)
pdf(file="sig/top10HeatmapB_mean.pdf", width=10, height=15)
DoHeatmap(cluster.averages, features = top10B$gene, slot = "scale.data", group.by = "ident",
          group.colors = pal, size = 3, draw.lines = FALSE) + scale_fill_gradientn(colours = rev(mapal))
dev.off()
saveRDS(objd.b, file = '/mnt/data_nvme1/lzh/test/AH/atlas/recluster/B/B.all.RDS')

# 06. Cell-type composition by group
library(ggplot2)
table(objd.b$TH_meta, objd.b$cellname_b)
rt1 <- data.frame(count=c(1910, 5421, 1754, 5920, 4453, 9878, 8003, 833, 10151, 157, 6822, 2818, 2417),
                   category=levels(objd.b$cellname_b))
rt1$category <- factor(rt1$category, levels=levels(objd.b$cellname_b))
rt1$fraction <- rt1$count / sum(rt1$count)
rt1$ymax <- cumsum(rt1$fraction)
rt1$ymin <- c(0, head(rt1$ymax, n=-1))
rt1$labelPosition <- (rt1$ymax + rt1$ymin) / 2
rt1$label <- paste0(rt1$category, "\n value: ", round(rt1$fraction*100,2), "%")
rt2 <- data.frame(count=c(2576, 5461, 1180, 4887, 4624, 9907, 6681, 926, 11608, 1173, 7533, 2740, 2586),
                   category=levels(objd.b$cellname_b))
rt2$category <- factor(rt2$category, levels=levels(objd.b$cellname_b))
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
table(objd.b$TH_meta, objd.b$cellname_b)
num2 <- table(objd.b$TH_meta, objd.b$cellname_b)
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
  "DIAPH3_GCB_DZ", "MKI67_GCB_DZ", "BACH2_GCB_LZ", "MAML3_GCB_LZ",
  "FCRL4_Memory_B", "PDE4D_Memory_B", "PTPRJ_Memory_B", "ACSM3_Naive_B","IGHD_Naive_B",
  "STAT4_Naive_B", "STEAP1B_Naive_B", "IGHG1_Plasma_cell","CD83_Pre_GCB"
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

# 07. B-cell pathway enrichment
library(msigdbr)
library(clusterProfiler)
library(enrichplot)
msigdbr_collections()
Hdb <- msigdbr(species = "Homo sapiens", category = "H") %>% dplyr::select(gs_name, gene_symbol)
Kdb <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG_LEGACY") %>% dplyr::select(gs_name, gene_symbol)
Wdb <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:WIKIPATHWAYS") %>% dplyr::select(gs_name, gene_symbol)
BPdb <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP") %>% dplyr::select(gs_name, gene_symbol)
df <- readRDS("/mnt/data_nvme1/lzh/test/AH/atlas/recluster/B/sig/sig.markers2.name.RDS")
df <- read.csv("/mnt/data_nvme1/lzh/test/AH/atlas/recluster/B/sig/top50B.csv")
table(df$cluster)
df2 <- df %>% filter(cluster %in% c("DIAPH3_GCB_DZ","MKI67_GCB_DZ","BACH2_GCB_LZ","MAML3_GCB_LZ"))
pal.gcb <-c("#FD8052","#F5B47F","#FFF3CA","#FFE2CE")

# Hallmark enrichment
hm1 <- enricher(df2[df2$cluster == "DIAPH3_GCB_DZ",]$gene, TERM2GENE=Hdb)
hm2 <- enricher(df2[df2$cluster == "MKI67_GCB_DZ",]$gene, TERM2GENE=Hdb)
hm3 <- enricher(df2[df2$cluster == "BACH2_GCB_LZ",]$gene, TERM2GENE=Hdb)
hm4 <- enricher(df2[df2$cluster == "MAML3_GCB_LZ",]$gene, TERM2GENE=Hdb)
hm.all <- list(hm1@result, hm2@result, hm3@result, hm4@result)
names(hm.all) <- c("DIAPH3_GCB_DZ","MKI67_GCB_DZ","BACH2_GCB_LZ","MAML3_GCB_LZ")
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
km1 <- enricher(df2[df2$cluster == "DIAPH3_GCB_DZ",]$gene, TERM2GENE=Kdb)
km2 <- enricher(df2[df2$cluster == "MKI67_GCB_DZ",]$gene, TERM2GENE=Kdb)
km3 <- enricher(df2[df2$cluster == "BACH2_GCB_LZ",]$gene, TERM2GENE=Kdb)
km4 <- enricher(df2[df2$cluster == "MAML3_GCB_LZ",]$gene, TERM2GENE=Kdb)
km.all <- list(km1@result, km2@result, km3@result, km4@result)
names(km.all) <- c("DIAPH3_GCB_DZ","MKI67_GCB_DZ","BACH2_GCB_LZ","MAML3_GCB_LZ")
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
bm1 <- enricher(df2[df2$cluster == "DIAPH3_GCB_DZ",]$gene, TERM2GENE=BPdb)
bm2 <- enricher(df2[df2$cluster == "MKI67_GCB_DZ",]$gene, TERM2GENE=BPdb)
bm3 <- enricher(df2[df2$cluster == "BACH2_GCB_LZ",]$gene, TERM2GENE=BPdb)
bm4 <- enricher(df2[df2$cluster == "MAML3_GCB_LZ",]$gene, TERM2GENE=BPdb)
bm.all <- list(bm1@result, bm2@result, bm3@result, bm4@result)
names(bm.all) <- c("DIAPH3_GCB_DZ","MKI67_GCB_DZ","BACH2_GCB_LZ","MAML3_GCB_LZ")
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
    geom_bar(stat = "identity", alpha = 0.5, fill = pal.gcb[i], width = 0.8) +
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

# 08. Nonnegative matrix factorization (NMF)
chmod +x step1_prepare_cnmf_input.r
Rscript step1_prepare_cnmf_input.r
chmod +x step2_run_until_kselect.sh
bash step2_run_until_kselect.sh > step2_run_until_kselect.log 2>&1 &
tail -f step2_run_until_kselect.log
bash step3_rerun_consensus.sh > step3_rerun_consensus.log 2>&1 &

# NMF plots
.libPaths("/home/lzh/miniconda3/envs/singlecell/lib/R/library")
library(Seurat)
library(ggplot2)
library(pbapply)
library(RColorBrewer)
set.seed(1234)
setwd("/mnt/data_nvme1/lzh/test/AH/atlas/recluster/B/cnmf/result")
theme_global <- theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank()) + theme(aspect.ratio = 1)
filedir <- "/mnt/data_nvme1/lzh/test/AH/atlas/recluster/B"
outdir <- file.path(filedir, "cnmf")
name <- "AH_B"
k <- 5
dt <- "0.1"
objd.b <- readRDS("/mnt/data_nvme1/lzh/test/AH/atlas/recluster/B/B.all.RDS")
seuratObj <- objd.b
usages <- read.table("/mnt/data_nvme1/lzh/test/AH/atlas/recluster/B/cnmf/AH_B/AH_B.usages.k_5.dt_0_1.consensus.txt", row.names = 1, header = TRUE, check.names = FALSE)
usages_normalized <- as.data.frame(t(apply(usages, 1, function(x) x / sum(x))))
seuratObj[["cnmf"]] <- CreateAssayObject(data = t(as.matrix(usages_normalized)))
DefaultAssay(seuratObj) <- "cnmf"
colnames(usages_normalized)
FeaturePlot(seuratObj, features = "1", pt.size = 0.75) + theme_global
pdf(file.path(getwd(), paste0(name, "_cNMF_score.k_", k, ".dt_", dt, "_manual.pdf")), useDingbats = FALSE, width = 7, height = 6)
for (i in colnames(usages_normalized)) {
    p <- FeaturePlot(seuratObj, features = i, pt.size = 0.75) +
        theme_global + ggtitle(paste0("Module ", i))
    print(p)
}
dev.off()
gene_spectra_score <- t(read.table("/mnt/data_nvme1/lzh/test/AH/atlas/recluster/B/cnmf/AH_B/AH_B.gene_spectra_score.k_5.dt_0_1.txt", header = T, row.names = 1, sep = "\t", check.names = F))
colnames(gene_spectra_score)
top_genes <- lapply(colnames(gene_spectra_score), function(x) {
    tmp <- gene_spectra_score[, x]
    names(sort(tmp, decreasing = TRUE))[1:30]
})
names(top_genes) <- colnames(gene_spectra_score)
top_genes
for (i in colnames(usages_normalized)) {
    seuratObj[[paste0("Program", i)]] <- usages_normalized[, i]
}
pal <-c("#FD8052","#F5B47F","#FFF3CA","#FFE2CE","#8DD2C5","#8DBA8F","#00A087","#4581B2","#B3CDE3","#354e97","#7FB2D5","#F4CAE4","#984EA3")
pdf("Program_vlnplot.pdf", width=6, height=5)
for (i in paste0("Program", colnames(usages_normalized))) {
   p <- VlnPlot(seuratObj, features = i, group.by = "cellname_b", pt.size = 0, cols = pal) + xlab("") +
                theme(legend.position = "none", axis.text.x = element_text(angle = 60))
   print(p)
}
dev.off()
colnames(gene_spectra_score) <- paste0("Program", seq_len(ncol(gene_spectra_score)))
colnames(usages) <- paste0("Program", seq_len(ncol(usages)))
colnames(usages_normalized) <- paste0("Program", seq_len(ncol(usages_normalized)))
pdf("Program_featureplot.pdf", width=6, height=4.5)
for (i in colnames(usages_normalized)) {
   p <- FeaturePlot(seuratObj, reduction = "umap", features = i, raster = F, pt.size = 0.6) +
                    theme_global + scale_color_gradientn(colors = rev(brewer.pal(10, "RdYlBu")))
   print(p)
}
dev.off()
pdf("Program_featureplot2.pdf", width=6, height=4.5)
for (i in colnames(usages_normalized)) {
   p <- FeaturePlot(seuratObj, reduction = "umap", features = i, raster = T, pt.size = 1) +
                    theme_global + scale_color_gradientn(colors = rev(brewer.pal(10, "RdYlBu")))
   print(p)
}
dev.off()
pdf("Program_spatial_byPatient.pdf", width = 8, height = 6)
for (p in levels(seuratObj$Patient)) {
  message("Plotting Patient: ", p)
  obj.p <- subset(seuratObj, subset = Patient == p)
  for (i in colnames(usages_normalized)) {
    p2 <- FeaturePlot(obj.p, features = i, reduction = "spatial", pt.size = 1) +
                    theme_global + scale_color_gradientn(colors = rev(brewer.pal(10, "RdYlBu"))) + ggtitle(paste0(p, "_", i)) + theme_global
    print(p2)
  }
}
dev.off()

# Top genes by NMF module
dim(gene_spectra_score)
head(rownames(gene_spectra_score))
head(colnames(gene_spectra_score))
library(dplyr)
library(tidyverse)
top_n <- 15
top_genes15 <- lapply(colnames(gene_spectra_score), function(m) {
    x <- gene_spectra_score[, m]
    names(sort(x, decreasing = TRUE))[1:top_n]
})
names(top_genes15) <- colnames(gene_spectra_score)
top_genes15
gene_order <- unlist(top_genes15)
gene_group <- rep(names(top_genes15), each = top_n)
heatmat_list <- lapply(seq_along(gene_order), function(i) {
    g <- gene_order[i]
    as.numeric(gene_spectra_score[g, ])
})
heatmat1 <- do.call(rbind, heatmat_list)
rownames(heatmat1) <- make.unique(gene_order)
colnames(heatmat1) <- colnames(gene_spectra_score)
mat_topgene <- t(heatmat1)
anno_col <- data.frame(GeneGroup = factor(gene_group, levels = names(top_genes15)))
rownames(anno_col) <- colnames(mat_topgene)
library(pheatmap)
pdf("nmf_top15_heatmap.pdf", width = 20, height = 8)
pheatmap(mat_topgene,
    cluster_rows = FALSE, cluster_cols = FALSE,
    scale = "none",
    color = colorRampPalette(c("#4581B2", "#fdf5f1", "#CB5C5B"))(100),
    annotation_col = anno_col,
    gaps_col = seq(top_n, ncol(mat_topgene) - top_n, by = top_n),
    border_color = "white", angle_col = 90,
    cellwidth = 16, cellheight = 26,
    fontsize_row = 13, fontsize_col = 13)
dev.off()

# Correlations among NMF modules
cor_mat <- cor(usages_normalized, method = "pearson")
cor_mat
pdf("nmf_cor_heatmap.pdf", width = 10, height = 6)
pheatmap(cor_mat,
    cluster_rows = FALSE, cluster_cols = FALSE,
    color = colorRampPalette(c("#4581B2", "#fdf5f1", "#CB5C5B"))(100),
    border_color = "white", angle_col = 90,
    cellwidth = 30, cellheight = 30,
    fontsize_row = 20, fontsize_col = 20)
dev.off()

# NMF activity by cell subtype
common_cells <- intersect(rownames(usages_normalized), colnames(seuratObj))
length(common_cells)
usage_df <- usages_normalized[common_cells, , drop = FALSE]
meta_df <- seuratObj@meta.data[common_cells, , drop = FALSE]
cluster_usage_mean <- aggregate(usage_df, by = list(CellType = meta_df$cellname_b), FUN = mean)
rownames(cluster_usage_mean) <- cluster_usage_mean$CellType
cluster_usage_mean <- cluster_usage_mean[, -1, drop = FALSE]
pdf("nmf_activity_across_ct_heatmap.pdf", width = 7, height = 6)
pheatmap(cluster_usage_mean,
    scale = "none",
    cluster_rows = FALSE, cluster_cols = FALSE,
    color = colorRampPalette(rev(c("#A50026", "#D73027", "#F46D43", "#FDAE61", "#FEE090", "#FFFFBF", "#E0F3F8", "#ABD9E9", "#74ADD1", "#4575B4", "#313695")))(50),
    border_color = "white",angle_col = 90,
    fontsize_row = 15, fontsize_col = 15,
    cellheight = 20, cellwidth = 30
)
dev.off()

# 09. Pathway enrichment by NMF program
library(dplyr)
library(tibble)
library(msigdbr)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
set.seed(1234)
outdir2 <- "/mnt/data_nvme1/lzh/test/AH/atlas/recluster/B/cnmf/enrich"
dir.create(outdir2, recursive = TRUE, showWarnings = FALSE)
msigdbr_collections()
Kdb <- msigdbr(species = "Homo sapiens",
    category = "C2", subcategory = "CP:KEGG_LEGACY"
) %>% dplyr::select(gs_name, gene_symbol)
BPdb <- msigdbr(species = "Homo sapiens",
    category = "C5", subcategory = "GO:BP"
) %>% dplyr::select(gs_name, gene_symbol)
top_genes <- lapply(colnames(gene_spectra_score), function(x) {
    tmp <- gene_spectra_score[, x]
    names(sort(tmp, decreasing = TRUE))[1:30]
})
names(top_genes) <- colnames(gene_spectra_score)
top_genes
saveRDS(top_genes, file.path(outdir2, "program_top30_genes.RDS"))
top_genes_df <- do.call(rbind, lapply(names(top_genes), function(nm) {data.frame(Program = nm,
        Gene = top_genes[[nm]], stringsAsFactors = FALSE)
}))
write.table(top_genes_df,
    file = file.path(outdir2, "program_top30_genes.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE)

# GO and KEGG enrichment of NMF programs
km.res <- lapply(names(top_genes), function(nm) {
    enricher(top_genes[[nm]], TERM2GENE = Kdb)
})
names(km.res) <- names(top_genes)
km.all <- lapply(km.res, function(x) if (is.null(x)) NULL else as.data.frame(x))
names(km.all) <- names(top_genes)
process_df.kegg <- function(km.res) {
    if (is.null(km.res) || nrow(km.res) == 0) return(NULL)
    km.res <- km.res[order(km.res$p.adjust, decreasing = FALSE), , drop = FALSE]
    km.res$ID <- km.res$Description
    km.res$ID <- gsub("^KEGG_", "", km.res$ID)
    km.res$ID <- gsub("_", " ", km.res$ID)
    return(km.res)
}
process_df <- function(res) {
    if (is.null(res) || nrow(res) == 0) return(NULL)
    df <- res[1:min(10, nrow(res)), , drop = FALSE]
    df$LogP <- -log10(df$p.adjust)
    df$labelx <- rep(0, nrow(df))
    df$labely <- seq_len(nrow(df))
    return(df)
}
km.all.fli <- lapply(km.all, process_df.kegg)
km.all.fli1 <- lapply(km.all.fli, process_df)
saveRDS(km.all.fli, file.path(outdir2, "program_kegg_processed.RDS"))
saveRDS(km.all.fli1, file.path(outdir2, "program_kegg_plot_input.RDS"))
bm.res <- lapply(names(top_genes), function(nm) {
    enricher(top_genes[[nm]], TERM2GENE = BPdb)
})
names(bm.res) <- names(top_genes)
bm.all <- lapply(bm.res, function(x) if (is.null(x)) NULL else as.data.frame(x))
names(bm.all) <- names(top_genes)
process_df.bk <- function(bm.res) {
    if (is.null(bm.res) || nrow(bm.res) == 0) return(NULL)
    bm.res <- bm.res[order(bm.res$p.adjust, decreasing = FALSE), , drop = FALSE]
    bm.res$ID <- bm.res$Description
    bm.res$ID <- gsub("^GOBP_", "", bm.res$ID)
    bm.res$ID <- gsub("_", " ", bm.res$ID)
    return(bm.res)
}
bm.all.fli <- lapply(bm.all, process_df.bk)
bm.all.fli1 <- lapply(bm.all.fli, process_df)
saveRDS(bm.all.fli, file.path(outdir2, "program_gobp_processed.RDS"))
saveRDS(bm.all.fli1, file.path(outdir2, "program_gobp_plot_input.RDS"))

# Enrichment plotting function
pal.p <- c("#B3CDE3", "#CCEBC5", "#DECBE4", "#FED9A6", "#FDDAEC")
names(pal.p) <- names(top_genes)
pal.p
plot_pathway <- function(df, i, pal) {
    if (is.null(df) || nrow(df) == 0) {
        return(ggplot() + theme_void() + ggtitle("No enriched terms"))
    }
    df$p.adjust[df$p.adjust == 0] <- 1e-300
    df$LogP <- -log10(df$p.adjust)
    df$ID <- tools::toTitleCase(tolower(as.character(df$ID)))
    df$ID <- reorder(df$ID, df$LogP)
    xmax <- max(df$LogP, na.rm = TRUE) * 1.15
    x_offset <- xmax * 0.02
    p <- ggplot(data = df, aes(x = LogP, y = ID)) +
        geom_bar(stat = "identity", alpha = 0.5, fill = pal[i], width = 0.8) +
        geom_text(aes(label = ID, x = x_offset), size = 5, hjust = 0) +
        coord_cartesian(xlim = c(0, xmax), clip = "off") +
        theme_classic() +
        theme(
            axis.text.y = element_blank(),
            axis.line.y = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.line.x = element_line(colour = 'black', linewidth = 1),
            axis.text.x = element_text(colour = 'black', size = 15),
            axis.ticks.x = element_line(colour = 'black', linewidth = 1),
            axis.title.x = element_text(colour = 'black', size = 16),
            plot.title = element_text(size = 18, face = "bold"),
            panel.border = element_rect(fill = NA, colour = "black", linewidth = 1, linetype = "solid"),
            plot.margin = margin(5.5, 5.5, 5.5, 30)
        ) +
        xlab("-log10(p.adjust)") +
        scale_x_continuous(expand = c(0, 0))
    return(p)
}

# Plot program enrichment
pdf(file.path(outdir2, "program_KEGG_barplot.pdf"), width = 7, height = 5)
for (i in seq_along(km.all.fli1)) {
    p <- plot_pathway(km.all.fli1[[i]], i, pal = pal.p) + ggtitle(names(km.all.fli1)[i])
    print(p)
}
dev.off()
pdf(file.path(outdir2, "program_GO_BP_barplot.pdf"), width = 7, height = 5)
for (i in seq_along(bm.all.fli1)) {
    p <- plot_pathway(bm.all.fli1[[i]], i, pal = pal.p) + ggtitle(names(bm.all.fli1)[i])
    print(p)
}
dev.off()

# Export pathway enrichment
for (nm in names(bm.all)) {
    df <- bm.all[[nm]]
    if (is.null(df) || nrow(df) == 0) next
    write.csv(df, file.path(outdir2, paste0(nm, "_GOBP_raw.csv")), row.names = FALSE)
}
for (nm in names(km.all)) {
    df <- km.all[[nm]]
    if (is.null(df) || nrow(df) == 0) next
    write.csv(df, file.path(outdir2, paste0(nm, "_KEGG_raw.csv")), row.names = FALSE)
}

# 10. GSEA and pathway activity
setwd("/mnt/data_nvme1/lzh/test/AH/atlas/recluster/B/gsea")
```{r}
.libPaths("/home/lzh/miniconda3/envs/singlecell/lib/R/library")
rm(list = ls())
set.seed(1234)
library(Seurat)
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(readr)
library(tidyr)
library(limma)
library(forcats)
library(patchwork)
library(ggpubr)
library(scales)
library(presto)
library(GseaVis)
objd.b <- readRDS("/mnt/data_nvme1/lzh/test/AH/atlas/recluster/B/B.all.RDS")
table(objd.b$cellname_b)
```

# IGHG1 plasma cell GSEA
```{r}
cellName = "IGHG1_Plasma_cell"
dir.create("./IGHG1_Plasma_cell")
seuratObj = subset(objd.b, cells = colnames(objd.b)[objd.b$cellname_b == cellName])
{
  library(msigdbr)
  m_df <- msigdbr(species = "Homo sapiens")
  pathways <- m_df %>% filter(m_df$gs_subcollection %in% c("CP:KEGG_LEGACY", "GO:BP")) %>% select(gs_name, gene_symbol)
  colnames(pathways) <- c("ont", "term")
  genes = rownames(seuratObj)[rowSums(as.matrix(seuratObj@assays$Spatial@counts) > 0) > ncol(seuratObj) * 1e-3]
  seuratObj2 <- subset(seuratObj, features=genes)
  markers <- presto::wilcoxauc(seuratObj2, 'TH_meta', seurat_assay = "Spatial")  # Presto-based test
  genes <- markers %>% dplyr::filter(group == "Y") %>% arrange(desc(logFC)) %>%
    dplyr::select(feature, logFC)
  ranks <- tibble::deframe(genes)
  GSEA_res <- clusterProfiler::GSEA(ranks, TERM2GENE=pathways, pvalueCutoff = 1, maxGSSize = 1500,
                                    nPermSimple=10000, eps=0)
  write.csv(GSEA_res@result, "./IGHG1_Plasma_cell/Y_GSEA_res.csv")
  pdf("./IGHG1_Plasma_cell/Y_GSEA_res.pdf", useDingbats=FALSE, width = 5.5, height = 5.5)
  GSEA_res_sig <- GSEA_res@result %>% filter(abs(NES) > 1 & p.adjust < 0.05)
  for (i in rownames(GSEA_res_sig)) {
    p <- GseaVis::gseaNb(object = GSEA_res,geneSetID = i,addPval = T,
                        pCol = "steelblue", pvalX = 0.625, pvalY = 0.8, pHjust = 0,
                        nesDigit = 4, pDigit = 4,
                        subPlot = 2) + ggtitle(i)
    print(p)
  }
  dev.off()
}
```

# IGHG1 plasma cell module scores
```{r}
{
  GSEA_res_sig = read.csv("./IGHG1_Plasma_cell/Y_GSEA_res.csv") %>% filter(abs(NES) > 1 & p.adjust < 0.05)
  library(msigdbr)
  m_df <- msigdbr(species = "Homo sapiens")
  selected_pathways <- m_df %>% filter(gs_name %in% c(GSEA_res_sig$ID))
  pathway_list <- sapply(unique(selected_pathways$gs_name), function(x){
    tmp <- selected_pathways %>% filter(gs_name == x)
    return(as.character(tmp$gene_symbol))
  })
  seuratObj = AddModuleScore(seuratObj, pathway_list)
  ES = as.data.frame(seuratObj@meta.data[, str_detect(colnames(seuratObj@meta.data), "^Cluster")])
  colnames(ES) = names(pathway_list)
  cluster_cols <- c("N" = "#9EC4BE", "Y" = "#FAC795")
  my_comparisons <- list(c("N", "Y"))
  pdf("./IGHG1_Plasma_cell/pathway_AddModuleScore_violin.pdf", useDingbats = F, width = 5, height = 5)
  for (i in 1:ncol(ES)) {
    df = data.frame(Activity = as.numeric(ES[,i]), cluster = seuratObj$TH_meta)
    p <- ggplot(df, aes(x = cluster, y = Activity, fill = cluster)) +
    geom_violin(width = 0.8, trim = FALSE) +
    geom_boxplot(width = 0.3, outlier.shape = NA) +
    scale_fill_manual(values = cluster_cols) +
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") +
    ylab("Activity") + xlab("") + ggtitle(colnames(ES)[i]) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 0.8,legend.position = "none")
    print(p)
    }
  dev.off()
  pdf("./IGHG1_Plasma_cell/pathway_AddModuleScore_violin2.pdf", useDingbats = F)
  for (i in 1:ncol(ES)) {
    df = data.frame(Activity = as.numeric(ES[,i]), cluster = seuratObj$Patient)
    p <- ggplot(df, aes(x=cluster,y=Activity))+
      geom_violin(width=1, aes(fill = cluster)) + ylab("Activity") + xlab("Patient") +
      geom_boxplot(width = 0.3, outlier.shape = NA) + ggtitle(colnames(ES)[i]) +
      theme_classic() + theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank()) +
      theme(aspect.ratio = 0.8)
    print(p)
  }
  dev.off()
  pdf("./IGHG1_Plasma_cell/pathway_AddModuleScore_heatmap.pdf", useDingbats = F, width = 10, height = 15)
  df = ES
  df$cluster = seuratObj$Patient
  df2 = df %>% pivot_longer(-cluster, names_to = "pathway", values_to = "Activity") %>%
    group_by(cluster, pathway) %>% summarize(mean_ = mean(Activity)) %>%
    pivot_wider(names_from = "cluster", values_from = "mean_") %>% as.data.frame()
  rownames(df2) = df2$pathway
  df2$pathway = NULL
  pheatmap::pheatmap(df2, scale = "row", fontsize = 9, treeheight_row = 10, treeheight_col = 10, cluster_cols = F, border_color = "white",
                     cellwidth = 15, cellheight = 10, color = colorRampPalette(c("navy", "white", "firebrick3"))(100))
  dev.off()
}
```

# CD83 pre-GC B-cell GSEA
```{r}
cellName = "CD83_Pre_GCB"
dir.create(cellName)
seuratObj = subset(objd.b, cells = colnames(objd.b)[objd.b$cellname_b == cellName])
{
  library(msigdbr)
  m_df <- msigdbr(species = "Homo sapiens")
  pathways <- m_df %>% filter(m_df$gs_subcollection %in% c("CP:KEGG_LEGACY", "GO:BP")) %>% select(gs_name, gene_symbol)
  colnames(pathways) <- c("ont", "term")
  genes = rownames(seuratObj)[rowSums(as.matrix(seuratObj@assays$Spatial@counts) > 0) > ncol(seuratObj) * 1e-3]
  seuratObj2 <- subset(seuratObj, features=genes)
  markers <- presto::wilcoxauc(seuratObj2, 'TH_meta', seurat_assay = "Spatial")  # Presto-based test
  genes <- markers %>% dplyr::filter(group == "Y") %>% arrange(desc(logFC)) %>%
    dplyr::select(feature, logFC)
  ranks <- tibble::deframe(genes)
  GSEA_res <- clusterProfiler::GSEA(ranks, TERM2GENE=pathways, pvalueCutoff = 1, maxGSSize = 1500,
                                    nPermSimple=10000, eps=0)
  write.csv(GSEA_res@result, "./CD83_Pre_GCB/Y_GSEA_res.csv")
  pdf("./CD83_Pre_GCB/Y_GSEA_res.pdf", useDingbats=FALSE, width = 5.5, height = 5.5)
  GSEA_res_sig <- GSEA_res@result %>% filter(abs(NES) > 1 & pvalue < 0.01)
  for (i in rownames(GSEA_res_sig)) {
    p <- GseaVis::gseaNb(object = GSEA_res,geneSetID = i,addPval = T,
                        pCol = "steelblue", pvalX = 0.625, pvalY = 0.8, pHjust = 0,
                        nesDigit = 4, pDigit = 4,
                        subPlot = 2) + ggtitle(i)
    print(p)
  }
  dev.off()
}
```

# CD83 pre-GC B-cell module scores
```{r}
{
  GSEA_res_sig = read.csv("./CD83_Pre_GCB/Y_GSEA_res.csv") %>% filter(abs(NES) > 1 & pvalue < 0.01)
  library(msigdbr)
  m_df <- msigdbr(species = "Homo sapiens")
  selected_pathways <- m_df %>% filter(gs_name %in% c(GSEA_res_sig$ID))
  pathway_list <- sapply(unique(selected_pathways$gs_name), function(x){
    tmp <- selected_pathways %>% filter(gs_name == x)
    return(as.character(tmp$gene_symbol))
  })
  seuratObj = AddModuleScore(seuratObj, pathway_list)
  ES = as.data.frame(seuratObj@meta.data[, str_detect(colnames(seuratObj@meta.data), "^Cluster")])
  colnames(ES) = names(pathway_list)
  cluster_cols <- c("N" = "#9EC4BE", "Y" = "#FAC795")
  my_comparisons <- list(c("N", "Y"))
  pdf("./CD83_Pre_GCB/pathway_AddModuleScore_violin.pdf", useDingbats = F, width = 5, height = 5)
  for (i in 1:ncol(ES)) {
    df = data.frame(Activity = as.numeric(ES[,i]), cluster = seuratObj$TH_meta)
    p <- ggplot(df, aes(x = cluster, y = Activity, fill = cluster)) +
    geom_violin(width = 0.8, trim = FALSE) +
    geom_boxplot(width = 0.3, outlier.shape = NA) +
    scale_fill_manual(values = cluster_cols) +
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") +
    ylab("Activity") + xlab("") + ggtitle(colnames(ES)[i]) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 0.8,legend.position = "none")
    print(p)
    }
  dev.off()
  pdf("./CD83_Pre_GCB/pathway_AddModuleScore_violin2.pdf", useDingbats = F)
  for (i in 1:ncol(ES)) {
    df = data.frame(Activity = as.numeric(ES[,i]), cluster = seuratObj$Patient)
    p <- ggplot(df, aes(x=cluster,y=Activity))+
      geom_violin(width=1, aes(fill = cluster)) + ylab("Activity") + xlab("Patient") +
      geom_boxplot(width = 0.3, outlier.shape = NA) + ggtitle(colnames(ES)[i]) +
      theme_classic() + theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank()) +
      theme(aspect.ratio = 0.8)
    print(p)
  }
  dev.off()
  pdf("./CD83_Pre_GCB/pathway_AddModuleScore_heatmap.pdf", useDingbats = F, width = 10, height = 16)
  df = ES
  df$cluster = seuratObj$Patient
  df2 = df %>% pivot_longer(-cluster, names_to = "pathway", values_to = "Activity") %>%
    group_by(cluster, pathway) %>% summarize(mean_ = mean(Activity)) %>%
    pivot_wider(names_from = "cluster", values_from = "mean_") %>% as.data.frame()
  rownames(df2) = df2$pathway
  df2$pathway = NULL
  pheatmap::pheatmap(df2, scale = "row", fontsize = 9, treeheight_row = 10, treeheight_col = 10, cluster_cols = F, border_color = "white",
                     cellwidth = 15, cellheight = 10, color = colorRampPalette(c("navy", "white", "firebrick3"))(100))
  dev.off()
}
```

# 11. Differential expression of selected B-cell subtypes
setwd("/mnt/data_nvme1/lzh/test/AH/atlas/recluster/B/deg")
.libPaths("/home/lzh/miniconda3/envs/singlecell/lib/R/library")
rm(list = ls())
set.seed(1234)
library(Seurat)
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(readr)
library(tidyr)
library(limma)
library(forcats)
library(patchwork)
library(EnhancedVolcano)
objd.b <- readRDS("/mnt/data_nvme1/lzh/test/AH/atlas/recluster/B/B.all.RDS")
table(objd.b$cellname_b)
logFCfilter=log2(1.25)
adjPvalFilter=1e-5
group = "TH_meta"
cellName = "CD83_Pre_GCB"
dir.create(cellName)
seuratObj = subset(objd.b, cells = colnames(objd.b)[objd.b$cellname_b == cellName])

# Compare subtype gene expression
diff_gene <- FindMarkers(seuratObj, min.pct = 0.25, logfc.threshold = 0.25, group.by = group,
                            ident.1 = "Y", ident.2 = "N", only.pos = FALSE, test.use="wilcox")
diff_gene <- diff_gene %>% arrange(desc(avg_log2FC), .by_group = TRUE)
sig_gene <- subset(diff_gene, p_val_adj < adjPvalFilter & abs(avg_log2FC) > logFCfilter)
saveRDS(diff_gene,paste0(cellName,"/diff_gene.RDS"))
write.csv(diff_gene,paste0(cellName,"/diff_gene.csv"))
saveRDS(sig_gene,paste0(cellName,"/sig_gene.RDS"))
write.csv(sig_gene,paste0(cellName,"/sig_gene.csv"))

# Volcano plots
pdf(paste0(cellName,"/diff_gene.pdf"),width = 6, height = 6.5)
EnhancedVolcano(diff_gene,lab = rownames(diff_gene),
                   x = 'avg_log2FC', y = 'p_val_adj',xlim = c(-2.5, 2.5),
                   pCutoff = adjPvalFilter, FCcutoff = logFCfilter,legendPosition = 'bottom',
                   pointSize = 1.5, labSize = 4, legendLabSize = 9,legendIconSize = 3.5, axisLabSize = 12,
                   title = cellName,subtitle = "",titleLabSize = 15, subtitleLabSize =12,
                   gridlines.major = FALSE,gridlines.minor = FALSE,
                   drawConnectors = TRUE,arrowheads = FALSE , max.overlaps = 16, shape = 16,
                   col=c("#D9D9D9", "#8DBA8F", "#4581B2", "#CB5C5B"), colAlpha = 1)
dev.off()

# Expression differences across groups
IGHG1_Plasma_cell <- subset(objd.b, cellname_b == "IGHG1_Plasma_cell")
Idents(IGHG1_Plasma_cell) <-"TH_meta"
gene.show <- c("SDC1","MZB1","XBP1","PRDM1","JCHAIN","IGHM","IGHA1","HLA-DRA","CD74")
pdf('/mnt/data_nvme1/lzh/test/AH/atlas/recluster/B/deg/IGHG1_Plasma_cell/dotplot.pdf',useDingbats=F, width = 8, height = 4)
DotPlot(IGHG1_Plasma_cell, features = gene.show,
    dot.min = 0, dot.scale = 6, idents = NULL, group.by = "TH_meta",
        split.by = NULL, cluster.idents = FALSE, scale = TRUE, scale.by = "radius", scale.min = NA, scale.max = NA) +
    scale_color_distiller(palette = "RdYlBu") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_y_discrete(limits = levels(IGHG1_Plasma_cell$TH_meta)) + labs(x="", y="") +
  theme(panel.border = element_rect(color = "black", size = 0.5, fill = NA))
dev.off()

# 12. CytoTRACE trajectory analysis
library(Seurat)
library(dplyr)
library(limma)
library(Matrix)
library(magrittr)
library(stringr)
library(patchwork)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(forcats)
library(tidyverse)
library(reticulate)
library(CytoTRACE)
packageVersion("CytoTRACE")
set.seed(1234)
setwd("/mnt/data_nvme1/lzh/test/AH/atlas/recluster/B/trajectory/cytotrace")
BC<-readRDS("/mnt/data_nvme1/lzh/test/AH/atlas/recluster/B/B.all.RDS")
DefaultAssay(BC) <- "Spatial"
BC2 <- BC
BC2.counts<-as.matrix(BC2@assays$Spatial@counts)  # UMI count matrix
BC2.phe<-as.character(BC2$cellname_b)
names(BC2.phe) <- rownames(BC2@meta.data)
results <- CytoTRACE(mat = BC2.counts, ncores = 8, enableFast = T)
emb<-data.frame(BC2@reductions$umap@cell.embeddings)
plotCytoTRACE(results, phenotype = BC2.phe, emb = emb)
plotCytoGenes(results, numOfGenes = 20)  # Automatically save plots
saveRDS(results, "cyto.result.RDS")
saveRDS(BC2.phe, "BC2.phe.RDS")

# Visualize CytoTRACE scores
rm(list = ls())
results<-readRDS("cyto.result.RDS")
BC2.phe<-readRDS("BC2.phe.RDS")
cyto.res<-data.frame(cbind(BC2.phe, results$CytoTRACE))
colnames(cyto.res)[2]<-"CytoTRACE"
cyto.res$CytoTRACE<-as.numeric(cyto.res$CytoTRAC)
library(ggplot2)
library(ggpubr)
cyto.res$BC2.phe = factor(cyto.res$BC2.phe, levels = c(
  "DIAPH3_GCB_DZ", "MKI67_GCB_DZ", "BACH2_GCB_LZ", "MAML3_GCB_LZ",
  "FCRL4_Memory_B", "PDE4D_Memory_B", "PTPRJ_Memory_B", "ACSM3_Naive_B","IGHD_Naive_B",
  "STAT4_Naive_B", "STEAP1B_Naive_B", "IGHG1_Plasma_cell","CD83_Pre_GCB"
))
pal1 <- c("#FD8052","#F5B47F","#FFF3CA","#FFE2CE","#8DD2C5","#8DBA8F","#00A087",
          "#4581B2","#B3CDE3","#354e97","#7FB2D5","#F4CAE4","#984EA3")
pdf(file="boxpolt.cytotrace.pdf", width=6, height=4.5)
ggboxplot(cyto.res, x="BC2.phe", y="CytoTRACE", color="BC2.phe",
          palette = pal1, width = 0.6,
          xlab=" ", ylab="Predicted ordering by CytoTRACE",
          add = "none", bxp.errorbar=T) + theme(legend.position="none") +
  theme(axis.text.x = element_text(colour="black", angle = 90, size = 10, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(colour="black", size = 10),
        axis.title.y = element_text(colour="black", size = 12))
dev.off()
cyto.res$BC2.phe = factor(cyto.res$BC2.phe, levels = c(
  "MKI67_GCB_DZ", "MAML3_GCB_LZ", "BACH2_GCB_LZ", "DIAPH3_GCB_DZ",
  "FCRL4_Memory_B", "IGHG1_Plasma_cell", "PDE4D_Memory_B", "CD83_Pre_GCB","PTPRJ_Memory_B",
  "STEAP1B_Naive_B", "ACSM3_Naive_B", "STAT4_Naive_B","IGHD_Naive_B"
))
pal2 <- c("#F5B47F","#FFE2CE","#FFF3CA","#FD8052","#8DD2C5","#F4CAE4","#8DBA8F",
          "#984EA3","#00A087","#7FB2D5","#4581B2","#354e97","#B3CDE3")
pdf(file="boxpolt.order.cytotrace.pdf", width=6, height=4.5)
ggboxplot(cyto.res, x="BC2.phe", y="CytoTRACE", color="BC2.phe",
          palette = pal2, width = 0.6,
          xlab=" ", ylab="Predicted ordering by CytoTRACE",
          add = "none", bxp.errorbar=T) + theme(legend.position="none") +
  theme(axis.text.x = element_text(colour="black", angle = 90, size = 10, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(colour="black", size = 10),
        axis.title.y = element_text(colour="black", size = 12))
dev.off()



# 13. Donor-level subcluster proportion comparison (AH vs ATH-A)
# Each donor is one independent replicate; absent subclusters count as 0%.
library(ggplot2)

sample_rds <- "/mnt/data_nvme1/lzh/test/AH/atlas/recluster/B/B.all.RDS"
sample_obj <- readRDS(sample_rds)
sample_meta <- sample_obj@meta.data[, c("Patient", "TH_meta", "cellname_b"), drop = FALSE]
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

subtype_levels <- levels(sample_obj$cellname_b)
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
          "B_donor_level_subcluster_proportions.csv"), row.names = FALSE)
write.csv(test_tbl[, c("Subtype", "N_AH", "N_ATHA", "Mean_AH",
                       "Mean_ATHA", "P_value", "Q_BH")],
          file.path(sample_output, "B_donor_level_subcluster_ttests.csv"),
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
    labs(title = "B-cell lineage subcluster proportions by donor",
         x = NULL, y = "Cells within lineage (%)", fill = NULL) +
    theme_classic(base_size = 11) +
    theme(legend.position = "top", strip.text = element_text(face = "bold"),
          axis.text.x = element_text(size = 9))
}

pdf(file.path(sample_output, "B_donor_level_subcluster_comparison.pdf"),
    width = 3.5, height = 4, useDingbats = FALSE)
print(plot_subtypes(donor_tbl, summary_tbl, test_tbl))
for (ct in subtype_levels) {
  print(plot_subtypes(
    donor_tbl[donor_tbl$Subtype == ct, ],
    summary_tbl[summary_tbl$Subtype == ct, ],
    test_tbl[test_tbl$Subtype == ct, ],
    columns = 1
  ) + labs(title = paste("B-cell lineage:", ct)))
}
dev.off()
