# scRNA-seq preprocessing, integration, and cell annotation
.libPaths("/home/lzh/miniconda3/envs/singlecell/lib/R/library")
set.seed(1234)
setwd("/mnt/data_nvme1/lzh/test/AH/sc_atlas/pre")
suppressWarnings(library(Seurat))
suppressWarnings(library(ggplot2))
suppressWarnings(library(patchwork))
suppressWarnings(library(tidyverse))
suppressWarnings(library(dplyr))
suppressWarnings(library(stringr))
packageVersion("Seurat")
packageVersion("SeuratObject")
packageVersion("ggplot2")
packageVersion("patchwork")
packageVersion("Matrix")
theme_global <- theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank()) + theme(aspect.ratio = 1)

# 01. Load and merge scRNA-seq samples
base_dir <- "/mnt/sda1/lzh/data/inhouse/AH/sc/scRNA"
sample_dirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)
sample_dirs <- sample_dirs[
  file.exists(file.path(sample_dirs, "filtered_feature_bc_matrix/barcodes.tsv.gz")) &
  file.exists(file.path(sample_dirs, "filtered_feature_bc_matrix/features.tsv.gz")) &
  file.exists(file.path(sample_dirs, "filtered_feature_bc_matrix/matrix.mtx.gz"))
]
sample_names <- basename(sample_dirs)
print(sample_names)
length(sample_names)
seurat_list <- lapply(seq_along(sample_dirs), function(i) {
      sample_dir <- sample_dirs[i]
      sample_name <- sample_names[i]
      message("Reading sample: ", sample_name)
      expr <- Read10X(data.dir = file.path(sample_dir, "filtered_feature_bc_matrix"))
      obj <- CreateSeuratObject(counts = expr, project = sample_name, min.cells = 3, min.features = 200)
      obj$orig.ident <- sample_name
      obj$sample <- sample_name
      obj <- RenameCells(object = obj, new.names = paste0(sample_name, "_", colnames(obj)))
      return(obj)
})
names(seurat_list) <- sample_names
if (length(seurat_list) == 1) {
  merged_obj <- seurat_list[[1]]
} else {
  merged_obj <- merge(x = seurat_list[[1]], y = seurat_list[-1], project = "AH_scRNA")
}

print(merged_obj)
dim(merged_obj)
head(colnames(merged_obj))
table(merged_obj$orig.ident)
DefaultAssay(merged_obj) <- "RNA"
saveRDS(merged_obj, file = "/mnt/data_nvme1/lzh/test/AH/sc_atlas/pre/merged_obj.RDS")

# 02. Quality control and sample metadata
merged_obj <- readRDS("/mnt/data_nvme1/lzh/test/AH/sc_atlas/pre/merged_obj.RDS")
str(merged_obj@meta.data)
merged_obj$orig.ident <- factor(merged_obj$orig.ident, levels = c("AH97","AH99","AH100","AH101","AH102","AH103","TH99","TH100","TH102"))
Idents(merged_obj) <- "orig.ident"
levels(merged_obj$orig.ident)
merged_obj[["percent.mt"]] <- PercentageFeatureSet(merged_obj, pattern = "^MT-")
merged_obj[["percent.rb"]] <- PercentageFeatureSet(merged_obj, pattern = "^RPL|^RPS")

pdf("QC_bySample_pre.pdf", width = 10, height = 6)
for (i in c("nCount_RNA","nFeature_RNA","percent.mt","percent.rb")) {
  print(VlnPlot(merged_obj, features = i, group.by = "orig.ident", pt.size = 0.01))
}
dev.off()

pdf("QC_bySample_pre2.pdf", width = 10, height = 6)
for (i in c("nCount_RNA","nFeature_RNA","percent.mt","percent.rb")) {
  print(VlnPlot(merged_obj, features = i, group.by = "orig.ident", pt.size = 0))
}
dev.off()

info <- read.csv("/mnt/data_nvme1/lzh/test/AH/sc_atlas/phenotype.csv", stringsAsFactors = FALSE)
rownames(info) <- info$sampleid
setdiff(unique(merged_obj$orig.ident), info$sampleid)
setdiff(info$sampleid, unique(merged_obj$orig.ident))
cell_metadata <- info[merged_obj$orig.ident, c("Gender", "Group","Tissue","Patient","Age_grade"), drop = FALSE]
rownames(cell_metadata) <- colnames(merged_obj)
merged_obj <- AddMetaData(merged_obj, metadata = cell_metadata)
merged_obj$Gender <- factor(merged_obj$Gender, levels = c("Female", "Male"))
merged_obj$Group <- factor(merged_obj$Group, levels = c("AH", "ATH"))
merged_obj$Tissue <- factor(merged_obj$Tissue, levels = c("Adenoid", "Tonsil"))
merged_obj$Patient <- factor(merged_obj$Patient, levels = c("Patient97","Patient99", "Patient100","Patient101","Patient102","Patient103"))
merged_obj$Age_grade <- factor(merged_obj$Age_grade, levels = c("II","III"))
str(merged_obj@meta.data)
table(merged_obj$Group)
table(merged_obj$Tissue)
table(merged_obj$Age_grade)

pdf("FeatureScatter.pdf", width = 8, height = 6)
FeatureScatter(merged_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()
quantile(merged_obj$nFeature_RNA, probs = c(0.25, 0.5, 0.75, 0.9, 0.95), na.rm = TRUE)
quantile(merged_obj$nCount_RNA, probs = c(0.25, 0.5, 0.75, 0.9, 0.95), na.rm = TRUE)
quantile(merged_obj$percent.mt, probs = c(0.25, 0.5, 0.75, 0.9, 0.95), na.rm = TRUE)
quantile(merged_obj$percent.rb, probs = c(0.25, 0.5, 0.75, 0.9, 0.95), na.rm = TRUE)
objd <- subset(merged_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 &
               nCount_RNA > 200 & nCount_RNA < 20000 &
               percent.mt < 10 & percent.rb < 15)
dim(objd)
table(objd$orig.ident)
pal.sample <- c("#8da1c1","#F6B9C3","#DAECBD", "#FFD569","#D8BFDF", "#F68A8B","#de9db5","#85B293","#ab94bc")
pal.patient <- c("#8da1c1","#F6B9C3","#DAECBD", "#FFD569","#D8BFDF", "#F68A8B")

pdf("QC_bySample_post.pdf", width = 18, height = 5)
VlnPlot(objd, features = c("nCount_RNA","nFeature_RNA","percent.mt","percent.rb"),
        group.by = "orig.ident", pt.size = 0.001,  ncol = 4, cols = pal.sample)
dev.off()

pdf("QC_bySample_post2.pdf", width = 18, height = 5)
VlnPlot(objd, features = c("nCount_RNA","nFeature_RNA","percent.mt","percent.rb"),
        group.by = "orig.ident", pt.size = 0, ncol = 4, cols = pal.sample)
dev.off()

# 03. Normalization, Harmony integration, and clustering
objd <- NormalizeData(objd, normalization.method = "LogNormalize", scale.factor = 10000) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>%
  ScaleData()
objd <- RunPCA(objd, npcs = 50)

pdf("pca_ElbowPlot.pdf", width = 10, height = 6)
ElbowPlot(objd, ndims = 50)
dev.off()
pc.select <- 1:30

# Harmony integration by patient
library(harmony)
pca_embed <- Embeddings(objd, reduction = "pca")[, pc.select, drop = FALSE]
h_embed <- harmony::HarmonyMatrix(data_mat = pca_embed,
  meta_data = objd@meta.data, vars_use = "Patient",  # Harmony batch variable: Patient
  do_pca = FALSE, verbose = TRUE)
objd[["harmony"]] <- CreateDimReducObject(embeddings = h_embed, key = "harmony_", assay = DefaultAssay(objd))

# Graph-based clustering
objd <- FindNeighbors(objd, reduction = "harmony", dims = pc.select)
library(clustree)
objd<-FindClusters(objd,resolution=seq(0.4,1.4,by=0.2))

pdf("clustree.pdf", width = 15, height = 12)
clustree(objd)
dev.off()
objd <- RunUMAP(objd, reduction = "harmony", dims = pc.select)
saveRDS(objd, file = "/mnt/data_nvme1/lzh/test/AH/sc_atlas/objd.RDS")
str(objd@meta.data)
table(objd$RNA_snn_res.1)
table(objd$RNA_snn_res.1.2)
table(objd$RNA_snn_res.1.4)
str(objd@meta.data)
objd$seurat_clusters <- objd$RNA_snn_res.1.4
Idents(objd) <- "seurat_clusters"
pdf("cluster_UMAP.pdf", width = 8, height = 6)
DimPlot(objd, reduction = "umap", group.by = "seurat_clusters", pt.size = 0.8, raster = FALSE, label = TRUE) + theme_global
dev.off()

# 04. Cell-type annotation
DefaultAssay(objd) <- "RNA"
Idents(objd) <- "seurat_clusters"
table(objd$seurat_clusters)

# Cell-type marker panels
markers1 <- c(
  # Airway / Mucosal epithelium (subtypes)
  "KRT8","KRT18","KRT19",  # Pan-epithelial/Glandular epithelium
  "KRT5","KRT14","TP63",  # Basal epithelial
  "SCGB1A1","SCGB3A1","SCGB3A2",  # Club / Secretory
  "MUC5AC","SPDEF","AGR2",  # Goblet / mucous
  "FOXJ1","TPPP3","DNAH5",  # Ciliated
  "KRT13","KRT4","SFN",  # Hillock / squamous-like
  "SPIB","TNFRSF11A","CCL20"  # M cells (microfold)
)
markers2 <- c(
  # B cell lineage (subtypes)
  "MS4A1","IGHD","TCL1A",  # Naive B
  "CD27","TNFRSF13B","FCRL4",  # Memory / activated B (class-switch often shows IGHA/IGHG)
  "AICDA","RGS13","MEF2B",  # Germinal center B
  "MZB1","XBP1","JCHAIN"  # Plasma cells
)
markers3 <- c(
  # T / NK cells (subtypes)
  "CD3D","CD3E","TRAC",  # Pan-T
  "CD4","IL7R","TRBC2",  # CD4 T (naive/memory)
  "CD8A","CD8B","GZMK",  # CD8 T
  "FOXP3","IL2RA","CTLA4",  # Treg
  "CXCR5","PDCD1","ICOS", "CD40LG",  # Tfh
  "NKG7","GNLY","PRF1"  # NK
)
markers4 <- c(
  # Myeloid / antigen-presenting cells (subtypes)
  "LYZ","S100A8","S100A9",  # Monocytes (inflammatory)
  "C1QA","C1QB","APOE",  # Macrophages
  "FCER1A","CD1C","CLEC10A",  # cDC2
  "CLEC9A","XCR1","BATF3",  # cDC1
  "IL3RA","TCF4","IRF7",  # pDC
  "TPSAB1","TPSB2","KIT",  # Mast cells
  "FCGR3B","CSF3R","MPO"  # Neutrophils
)
markers5 <- c(
  # Stromal / vascular cells (subtypes)
  "COL1A1","COL3A1","DCN",  # Fibroblasts
  "ACTA2","TAGLN","MYH11",  # Smooth muscle / myofibroblast (MYH11 often cleaner than MYL9)
  "PECAM1","VWF","KDR",  # Endothelial
  "PROX1","PDPN","LYVE1",  # Lymphatic endothelial
  "RGS5","PDGFRB","CSPG4",  # Pericytes
  "MKI67","TOP2A","STMN1"  # Proliferative"
)
markers_all <- c(markers1, markers2, markers3, markers4, markers5)

pdf("AH_all_markerBubble.pdf", useDingbats=FALSE, width = 25, height = 12)
DotPlot(objd, features = markers_all, assay = "RNA", group.by = "seurat_clusters") + scale_color_distiller(palette = "RdYlBu") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
dev.off()

cluster2celltype_self <- c("0"="Naive_B", "1"="Memory_B",
                           "2"="Naive_B", "3"="GC_B",
                           "4"="Naive_B", "5"="Tfh",
                           "6"="Naive_B", "7"="Tfh",
                           "8"="GC_B", "9"="CD4_T",
                           "10"="Memory_B", "11"="Naive_B",
                           "12"="Plasma_cell", "13"="Naive_B",
                           "14"="Tfh", "15"="GC_B",
                           "16"="CD8_T", "17"="Treg",
                           "18"="CD4_T", "19"="Unknown",
                           "20"="Tfh", "21"="GC_B",
                           "22"="CD8_T", "23"="GC_B",
                           "24"="Naive_B", "25"="Macrophage",
                           "26"="Fibroblast", "27"="Glandular_epithelium",
                           "28"="NK", "29"="GC_B",
                           "30"="Unknown", "31"="Endothelial",
                           "32"="pDC", "33"="Unknown",
                           "34"="Unknown", "35"="Unknown")
objd[['celltype']]<-unname(cluster2celltype_self[objd@meta.data$seurat_clusters])
table(objd@meta.data$celltype)
dim(objd)
objd <- subset(objd, celltype != "Unknown")
dim(objd)
objd$celltype <- factor(objd$celltype, levels = c("CD4_T","CD8_T","Endothelial","Fibroblast",
                              "GC_B","Glandular_epithelium","Macrophage","Memory_B","Naive_B","NK","pDC","Plasma_cell","Tfh","Treg"))
Idents(objd) <- "celltype"
str(objd@meta.data)
table(objd$celltype)
table(objd$Group)
pal <- c("#7FB2D5","#B3CDE3","#BEDADA","#8DBA8F","#6A3D9A","#00A087","#F39B7FFF",
        "#984EA3","#BC9DCC","#57C3F3","#F5B47F","#D37295","#F4CAE4","#4581B2","#F2F2F2")
pal.sample <- c("#8da1c1","#F6B9C3","#DAECBD", "#FFD569","#D8BFDF", "#F68A8B","#de9db5","#85B293","#ab94bc")
pal.patient <- c("#8da1c1","#F6B9C3","#DAECBD", "#FFD569","#D8BFDF", "#F68A8B")

# 05. Sample and clinical metadata UMAPs
pdf("AH_all_INFO_UMAP.pdf", width = 8, height = 6)
DimPlot(objd, reduction = "umap", group.by = "orig.ident", cols = pal.sample, pt.size = 0.8, raster = FALSE, order = TRUE) + ggtitle("Sample") + theme_global
DimPlot(objd, reduction = "umap", group.by = "Patient", cols = pal.patient, pt.size = 0.8, raster = FALSE) + ggtitle("Patient") + theme_global
DimPlot(objd, reduction = "umap", group.by = "Group", cols = c("#9EC4BE","#FAC795"), pt.size = 0.8, raster = FALSE, order = TRUE) + ggtitle("Group") + theme_global
DimPlot(objd, reduction = "umap", group.by = "Tissue", cols = c("#F19EC2","#f39800"), pt.size = 0.8, raster = FALSE) + ggtitle("Tissue") + theme_global
DimPlot(objd, reduction = "umap", group.by = "Gender", cols = c("#E7C1D7","#AFB0D0"), pt.size = 0.8, raster = FALSE) + ggtitle("Gender") + theme_global
dev.off()

pdf("AH_all_INFO_UMAP2.pdf", width = 8, height = 6)
DimPlot(objd, reduction = "umap", group.by = "orig.ident", cols = pal.sample, pt.size = 0.8, raster = TRUE, order = TRUE) + ggtitle("Sample") + theme_global
DimPlot(objd, reduction = "umap", group.by = "Patient", cols = pal.patient, pt.size = 0.8, raster = TRUE) + ggtitle("Patient") + theme_global
DimPlot(objd, reduction = "umap", group.by = "Group", cols = c("#9EC4BE","#FAC795"), pt.size = 0.8, raster = TRUE) + ggtitle("Group") + theme_global
DimPlot(objd, reduction = "umap", group.by = "Tissue", cols = c("#F19EC2","#f39800"), pt.size = 0.8, raster = TRUE) + ggtitle("Tissue") + theme_global
DimPlot(objd, reduction = "umap", group.by = "Gender", cols = c("#E7C1D7","#AFB0D0"), pt.size = 0.8, raster = TRUE) + ggtitle("Gender") + theme_global
dev.off()

# 06. Post-annotation quality control
pdf("AH_all_QC_bySample.pdf", width = 18, height = 5)
VlnPlot(objd, features = c("nCount_RNA","nFeature_RNA","percent.mt","percent.rb"),
        group.by = "orig.ident", pt.size = 0, ncol = 4, cols = pal.sample)
dev.off()

# Cell-type marker visualization
library(grid)
Idents(objd) <- "celltype"
markers<-c("CD3E","CD4","IL7R",  # CD4_T
          "CD8A","CD8B","GZMK",  # CD8_T
          "PECAM1","VWF","KDR",  # Endothelial
          "COL1A1","COL3A1","DCN",  # Fibroblast
          "AICDA","RGS13","MEF2B",  # GC_B
          "KRT8","KRT18","KRT19",  # Glandular_epithelium
          "C1QA","C1QB","APOE",  # Macrophage
          "CD27","TNFRSF13B","FCRL4",  # Memory_B
          "MS4A1","IGHD","TCL1A",  # Naive_B
          "NKG7","GNLY","PRF1",  # NK
          "IL3RA","TCF4","IRF7",  # pDC
          "MZB1","XBP1","JCHAIN",  # Plasma_cell
          "CD40LG","PDCD1","ICOS",  # Tfh
          "FOXP3","IL2RA","CTLA4"  # Treg
)

pdf("DotPlot.celltype.pdf", width = 11, height = 5)
DotPlot(objd, features = markers,  cols = c("grey", "red"),
    dot.min = 0, dot.scale = 6, idents = NULL, group.by = "celltype",
        split.by = NULL, cluster.idents = FALSE, scale = TRUE, scale.by = "radius", scale.min = NA, scale.max = NA)+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_y_discrete(limits = levels(objd$celltype)) +
  theme(panel.border = element_rect(color="black"),
        panel.spacing = unit(1, "mm"),
        strip.text = element_text(margin=margin(b=3, unit="mm")),
        strip.placement = 'outlet',
        axis.line = element_blank()) + labs(x="", y="")
DotPlot(objd, features = markers,
    dot.min = 0, dot.scale = 6, idents = NULL, group.by = "celltype",
        split.by = NULL, cluster.idents = FALSE, scale = TRUE, scale.by = "radius", scale.min = NA, scale.max = NA) +
    scale_color_distiller(palette = "RdYlBu") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_y_discrete(limits = levels(objd$celltype)) + labs(x="", y="") +
  theme(panel.border = element_rect(color = "black", size = 0.5, fill = NA))
dev.off()

# Cell-type UMAPs
pdf("AH_all_Celltype_UMAP.pdf", width = 8, height = 6)
DimPlot(objd, reduction = 'umap', group.by = 'seurat_clusters', pt.size = 0.8, raster = FALSE,label = TRUE) + ggtitle("Cluster") + theme_global
DimPlot(objd, reduction = 'umap', group.by = 'celltype', pt.size = 0.8, cols = pal, raster = FALSE) + ggtitle("Celltype") + theme_global
dev.off()

pdf("AH_all_Celltype_UMAP2.pdf", width = 8, height = 6)
DimPlot(objd, reduction = 'umap', group.by = 'seurat_clusters', pt.size = 0.8, raster = TRUE,label = TRUE) + ggtitle("Cluster") + theme_global
DimPlot(objd, reduction = 'umap', group.by = 'celltype', pt.size = 0.8, cols = pal, raster = TRUE) + ggtitle("Celltype") + theme_global
dev.off()

# 07. Cell-cycle scoring and feature-level QC
str(objd@meta.data)
s.genes <- Seurat::cc.genes$s.genes
g2m.genes <- Seurat::cc.genes$g2m.genes
objd <- CellCycleScoring(objd,s.features = s.genes, g2m.features = g2m.genes,set.ident = T)  # Cell-cycle scoring

pdf("AH_all_nCount_UMAP.pdf", width = 8, height = 6)
FeaturePlot(objd, reduction = 'umap', features = 'nCount_RNA', pt.size = 0.8, raster = FALSE) + ggtitle("nCount") +
            theme_global + scale_color_distiller(palette = "RdYlBu")
FeaturePlot(objd, reduction = 'umap', features = 'nFeature_RNA', pt.size = 0.8, raster = FALSE) + ggtitle("nFeature") +
            theme_global + scale_color_distiller(palette = "RdYlBu")
dev.off()

Idents(objd) <- "celltype"
median(objd$nCount_RNA)
median(objd$nFeature_RNA)
table(objd$Group,objd$Tissue)
table(objd$orig.ident)
sample2group2 <- c("AH97"="AH", "AH99"="ATH-A",
                   "AH100"="ATH-A", "AH101"="AH",
                   "AH102"="ATH-A", "AH103"="AH",
                   "TH99"="ATH-T", "TH100"="ATH-T", "TH102"="ATH-T")
objd[['Group2']]<-unname(sample2group2[objd@meta.data$sample])
objd$Group2 <- factor(objd$Group2, levels = c("AH", "ATH-A", "ATH-T"))
table(objd$Group2)
table(objd$Group2,objd$celltype)
saveRDS(objd, "/mnt/data_nvme1/lzh/test/AH/sc_atlas/AH_all_rna.fli.RDS")

# 08. Sample-level PCA of mean RNA expression
# Aggregate expression per sample, select variable genes, and run PCA.
group_colors <- c("AH" = "#F7D8B7", "ATH-A" = "#E7C1D7", "ATH-T" = "#AFB0D0")
library(Seurat)
library(dplyr)
library(ggplot2)
library(matrixStats)
run_sample_pca <- function(seu, top_n = 3000, assay = "RNA", slot = "data",
                           group_colors = NULL) {
  expr <- GetAssayData(seu, assay = assay, slot = slot)
  sample_ids <- seu@meta.data$orig.ident
  sample_levels <- unique(as.character(sample_ids))
  avg_list <- lapply(sample_levels, function(s) {
    cells_use <- colnames(seu)[sample_ids == s]
    Matrix::rowMeans(expr[, cells_use, drop = FALSE])
  })
  names(avg_list) <- sample_levels
  avg_expr <- do.call(cbind, avg_list)
  colnames(avg_expr) <- sample_levels
  avg_expr <- as.matrix(avg_expr)
  gene_var <- matrixStats::rowVars(avg_expr)
  gene_var[is.na(gene_var)] <- 0
  top_n <- min(top_n, nrow(avg_expr))
  top_genes <- names(sort(gene_var, decreasing = TRUE))[1:top_n]
  avg_expr_top <- t(avg_expr[top_genes, , drop = FALSE])  # Rows are samples; columns are genes
  pca_res <- prcomp(avg_expr_top, center = TRUE, scale. = TRUE)
  pca_df <- as.data.frame(pca_res$x[, 1:2])
  pca_df$orig.ident <- rownames(pca_df)
  meta_df <- seu@meta.data %>% select(orig.ident, Group2) %>% distinct()
  meta_df$orig.ident <- as.character(meta_df$orig.ident)
  pca_df$orig.ident <- as.character(pca_df$orig.ident)
  pca_df <- left_join(pca_df, meta_df, by = "orig.ident")
  pca_df$Group2 <- factor(pca_df$Group2, levels = names(group_colors))
  percent_var <- round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2)), 1)
  p <- ggplot(pca_df, aes(x = PC1, y = PC2, fill = Group2, color = Group2)) +
    stat_ellipse(aes(fill = Group2), geom = "polygon", alpha = 0.25, color = NA, level = 0.95) +
    geom_point(shape = 21, size = 5, stroke = 1.1, alpha = 0.9) +
    ggrepel::geom_text_repel(aes(label = orig.ident), size = 4, color = "black",
                             box.padding = 0.4, point.padding = 0.3, max.overlaps = 100) +
    scale_fill_manual(values = group_colors, na.value = "grey80") +
    scale_color_manual(values = group_colors, na.value = "grey50") +
    labs(title = "Group", x = paste0("PC1 (", percent_var[1], "% variance)"),
         y = paste0("PC2 (", percent_var[2], "% variance)")) +
    theme_classic(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
          axis.title = element_text(face = "bold"), axis.text = element_text(color = "black"),
          legend.title = element_blank(), legend.text = element_text(size = 12))
  return(list(pca = pca_res, plot = p,df = pca_df, avg_expr = avg_expr))
}

res_group <- run_sample_pca(objd.all, group_colors = group_colors)
pdf("pre/sample_PCA_by_treatment.pdf", width = 6, height = 5)
print(res_group$plot)
dev.off()
