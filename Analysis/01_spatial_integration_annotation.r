# Pediatric adenoid spatial transcriptomics
# Integration, QC, clustering, cell-type annotation, and sample-level PCA
# Input/output paths are unchanged from the original analysis environment.

# 01. Environment and packages
.libPaths("/home/zihan/miniconda3/envs/singlecell/lib/R/library")

setwd("/data/YQ/linzh/test/AH/atlas")
set.seed(1234)

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

theme_global <- theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  theme(aspect.ratio = 1)

# 02. Spatial samples
samples <- list(
  list(
    sample = "AH5_AH7",
    matrix_dir = "/data/YQ/linzh/inhouse_data/ST/AH/20251230/25070088_AH5_AH7_expression/AH5_AH7_filtered_feature_bc_matrix",
    coord_file = "/data/YQ/linzh/inhouse_data/ST/AH/20251230/25070088_AH5_AH7_expression/AH5_AH7_filtered_feature_bc_matrix/cell_locations.tsv.gz",
    size_x = 55128, size_y = 19906,
    dapi_png = "/data/YQ/linzh/inhouse_data/ST/AH/20251230/25070088_AH5_AH7_expression/AH5_AH7_aligned_DAPI.png",
    he_png   = "/data/YQ/linzh/inhouse_data/ST/AH/20251230/25070088_AH5_AH7_expression/AH5_AH7_aligned_HE_TIMG.png"
  ),
  list(
    sample = "AH15_AH35",
    matrix_dir = "/data/YQ/linzh/inhouse_data/ST/AH/20251230/25070088_AH15_AH35_expression/AH15_AH35_filtered_feature_bc_matrix",
    coord_file = "/data/YQ/linzh/inhouse_data/ST/AH/20251230/25070088_AH15_AH35_expression/AH15_AH35_filtered_feature_bc_matrix/cell_locations.tsv.gz",
    size_x = 55128, size_y = 19906,
    dapi_png = "/data/YQ/linzh/inhouse_data/ST/AH/20251230/25070088_AH15_AH35_expression/AH15_AH35_aligned_DAPI.png",
    he_png   = "/data/YQ/linzh/inhouse_data/ST/AH/20251230/25070088_AH15_AH35_expression/AH15_AH35_aligned_HE_TIMG.png"
  ),
  list(
    sample = "AH17_AH8",
    matrix_dir = "/data/YQ/linzh/inhouse_data/ST/AH/20251230/25070088_AH17_AH8_expression/AH17_AH8_filtered_feature_bc_matrix",
    coord_file = "/data/YQ/linzh/inhouse_data/ST/AH/20251230/25070088_AH17_AH8_expression/AH17_AH8_filtered_feature_bc_matrix/cell_locations.tsv.gz",
    size_x = 55128, size_y = 19906,
    dapi_png = "/data/YQ/linzh/inhouse_data/ST/AH/20251230/25070088_AH17_AH8_expression/AH17_AH8_aligned_DAPI.png",
    he_png   = "/data/YQ/linzh/inhouse_data/ST/AH/20251230/25070088_AH17_AH8_expression/AH17_AH8_aligned_HE_TIMG.png"
  )
)

# 03. Read spatial expression, coordinates and tissue images
suppressWarnings(library(Seurat))
suppressWarnings(library(Matrix))
suppressWarnings(library(png))
suppressWarnings(library(grid))
suppressWarnings(library(base64enc))

read_seekgene_spatial <- function(one) {
  message("Reading: ", one$sample)

  # Count matrix
  mat <- Read10X(one$matrix_dir)
  obj <- CreateSeuratObject(counts = mat, project = "AH", assay = "Spatial")
  obj$Sample <- one$sample

  # Barcode-associated spatial coordinates
  spatial_df <- read.table(one$coord_file, row.names = 1, sep = "\t", header = TRUE)
  colnames(spatial_df) <- c("spatial_1", "spatial_2")
  spatial_matrix <- as.matrix(spatial_df)

  # Align coordinates to Seurat cell order.
  spatial_matrix_sorted <- spatial_matrix[match(rownames(obj@meta.data), rownames(spatial_matrix)), ]
  stopifnot(all(rownames(obj@meta.data) == rownames(spatial_matrix_sorted)))

  obj@reductions$spatial <- CreateDimReducObject(
    embeddings = spatial_matrix_sorted,
    key = "spatial_",
    assay = "Spatial"
  )

  # Store image data and slide dimensions.
  samplename <- one$sample
  obj@misc$info[[samplename]]$size_x <- as.integer(one$size_x)
  obj@misc$info[[samplename]]$size_y <- as.integer(one$size_y)

  # DAPI image
  obj@misc$info[[samplename]]$img <- base64enc::dataURI(file = one$dapi_png)
  img_gg <- png::readPNG(one$dapi_png)
  obj@misc$info[[samplename]]$img_gg <- grid::rasterGrob(
    img_gg, interpolate = FALSE,
    width = unit(1, "npc"), height = unit(1, "npc")
  )

  # H&E image
  obj@misc$info[[samplename]]$img_he <- base64enc::dataURI(file = one$he_png)
  img_he <- png::readPNG(one$he_png)
  obj@misc$info[[samplename]]$img_he_gg <- grid::rasterGrob(
    img_he, interpolate = FALSE,
    width = unit(1, "npc"), height = unit(1, "npc")
  )
  return(obj)
}

# 04. Merge spatial samples
objs <- lapply(samples, read_seekgene_spatial)
# Prefix cell barcodes with sample identifiers during merging.
spatial.all <- merge(x = objs[[1]], y = objs[-1],
  add.cell.ids = vapply(samples, `[[`, "", "sample"), project = "AH_all")

DefaultAssay(spatial.all) <- "Spatial"

table(spatial.all$Sample)
dim(spatial.all)

# 05. Quality control and sample metadata
spatial.all[["percent.mt"]] <- PercentageFeatureSet(spatial.all, pattern = "^MT-")
spatial.all[["percent.rb"]] <- PercentageFeatureSet(spatial.all, pattern = "^RPL|^RPS")

pdf("pre/AH_all_QC_bySample.pdf", width = 12, height = 6)
VlnPlot(
  spatial.all,
  features = c("nCount_Spatial", "nFeature_Spatial", "percent.mt", "percent.rb"),
  group.by = "Sample", pt.size = 0.01, ncol = 4
)
dev.off()

# Match cell metadata by barcode.
meta_df <- read.csv("/data/YQ/linzh/test/AH/sample_info.csv", row.names = 1, check.names = FALSE,stringsAsFactors = FALSE)

# Match metadata to Seurat cell order.
meta_df2 <- meta_df[colnames(spatial.all), , drop = FALSE]  # Unmatched cells receive NA metadata.
spatial.all <- AddMetaData(spatial.all, metadata = meta_df2)
spatial.all <- subset(spatial.all, subset = Patient != "NA")

dim(spatial.all)

# 06. Normalization, dimensionality reduction and clustering
spatial.all <- NormalizeData(spatial.all, normalization.method = "LogNormalize", scale.factor = 10000) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>%
  ScaleData()

spatial.all <- RunPCA(spatial.all, npcs = 30)

pdf("pre/AH_all_ElbowPlot.pdf", width = 8, height = 5)
ElbowPlot(spatial.all, ndims = 30)
dev.off()

pc.select <- 1:30

# Harmony integration by patient
library(harmony)
spatial.all$Patient <- factor(spatial.all$Patient, levels = c("AH5","AH7","AH8","AH15","AH17","AH35"))
pca_embed <- Embeddings(spatial.all, reduction = "pca")[, pc.select, drop = FALSE]

h_embed <- harmony::HarmonyMatrix(
  data_mat = pca_embed,
  meta_data = spatial.all@meta.data, vars_use = "Patient",  # Patient-level batch variable
  do_pca = FALSE, verbose = TRUE
)

spatial.all[["harmony"]] <- CreateDimReducObject(embeddings = h_embed, key = "harmony_", assay = DefaultAssay(spatial.all))

spatial.all <- FindNeighbors(spatial.all, reduction = "harmony", dims = pc.select)
spatial.all <- FindClusters(spatial.all, resolution = seq(0.4, 1.4, 0.2))
spatial.all <- RunUMAP(spatial.all, reduction = "harmony", dims = pc.select)

str(spatial.all@meta.data)

spatial.all$seurat_clusters <- spatial.all$Spatial_snn_res.1.4
spatial.all$seurat_clusters <- as.character(spatial.all$seurat_clusters)
spatial.all$seurat_clusters <- factor(spatial.all$seurat_clusters, levels = c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17",
                                      "18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36"))
levels(spatial.all$seurat_clusters)
table(spatial.all$seurat_clusters)
Idents(spatial.all) <- "seurat_clusters"

# Patient-level UMAP visualization
pdf("pre/AH_all_DimPlot_batcheffect.pdf", width = 8, height = 6)
DimPlot(spatial.all, reduction = "umap", group.by = "Patient", cols = pal.patient, pt.size = 0.8, raster = FALSE) + ggtitle("Before Harmony") + theme_global
dev.off()

# Spatial cell distributions
spatial.all$spatial_1 <- NA_real_
spatial.all$spatial_2 <- NA_real_

for (one in samples) {
  s <- one$sample
  message("Adding coords for: ", s)
  # Coordinate rows are indexed by barcode.
  df <- read.table(one$coord_file, row.names = 1, sep = "\t", header = TRUE)
  colnames(df) <- c("spatial_1", "spatial_2")
  # Retrieve cells from the current chip.
  cells <- colnames(spatial.all)[spatial.all$Sample == s]
  if (length(cells) == 0) next
  # Remove the merged sample prefix to recover each barcode.
  barcode <- sub("^.*_", "", cells)
  # Match coordinates to barcodes and assign to metadata.
  m <- as.matrix(df[match(barcode, rownames(df)), c("spatial_1","spatial_2"), drop = FALSE])
  rownames(m) <- cells
  spatial.all$spatial_1[cells] <- m[, 1]
  spatial.all$spatial_2[cells] <- m[, 2]
}

# Count missing spatial coordinates.
sum(is.na(spatial.all$spatial_1))
# Create the spatial reduction for plotting.
coords <- as.matrix(spatial.all@meta.data[, c("spatial_1", "spatial_2")])
rownames(coords) <- colnames(spatial.all)

spatial.all[["spatial"]] <- CreateDimReducObject(
  embeddings = coords,
  key = "spatial_",
  assay = "Spatial"
)

Reductions(spatial.all)

pdf("pre/AH_all_spatial_byPatient_clusters.pdf", width = 8, height = 6)
for (p in levels(spatial.all$Patient)) {
  message("Plotting Patient: ", p)
  obj.p <- subset(spatial.all, subset = Patient == p)
  p1 <- DimPlot(obj.p, reduction = "spatial", group.by = "seurat_clusters", pt.size = 1) + ggtitle(p) + theme_global
  print(p1)
}
dev.off()

# Add patient sex metadata.
patient2gender <- c("AH5"="Female", "AH7"="Female",
                    "AH8"="Male", "AH15"="Male",
                    "AH17"="Male", "AH35"="Female")

spatial.all[['gender']]<-unname(patient2gender[spatial.all@meta.data$Patient])
table(spatial.all$gender)

# 07. Marker-based cell-type annotation
DefaultAssay(spatial.all) <- "Spatial"
Idents(spatial.all) <- "seurat_clusters"
# Canonical marker panels
markers1 <- c(
  # Airway and mucosal epithelium
  "KRT8","KRT18","KRT19",  # Pan-epithelial/Glandular epithelium
  "KRT5","KRT14","TP63",  # Basal epithelial
  "SCGB1A1","SCGB3A1","SCGB3A2",  # Club / Secretory
  "MUC5AC","SPDEF","AGR2",  # Goblet / mucous
  "FOXJ1","TPPP3","DNAH5",  # Ciliated
  "KRT13","KRT4","SFN",  # Hillock / squamous-like
  "SPIB","TNFRSF11A","CCL20"  # M cells (microfold)
)
markers2 <- c(
  # B-cell lineage
  "MS4A1","IGHD","TCL1A",  # Naive B
  "CD27","TNFRSF13B","FCRL4",  # Memory / activated B (class-switch often shows IGHA/IGHG)
  "AICDA","RGS13","MEF2B",  # Germinal center B
  "MZB1","XBP1","JCHAIN"  # Plasma cells
)
markers3 <- c(
  # T- and NK-cell lineages
  "CD3D","CD3E","TRAC",  # Pan-T
  "CD4","IL7R","TRBC2",  # CD4 T (naive/memory)
  "CD8A","CD8B","GZMK",  # CD8 T
  "FOXP3","IL2RA","CTLA4",  # Treg
  "CXCR5","PDCD1","ICOS",
  "NKG7","GNLY","PRF1"
)
markers4 <- c(
  # Myeloid and antigen-presenting cells
  "LYZ","S100A8","S100A9",  # Monocytes (inflammatory)
  "C1QA","C1QB","APOE",  # Macrophages
  "FCER1A","CD1C","CLEC10A",  # cDC2
  "CLEC9A","XCR1","BATF3",  # cDC1
  "IL3RA","TCF4","IRF7",  # pDC
  "TPSAB1","TPSB2","KIT",  # Mast cells
  "FCGR3B","CSF3R","MPO"  # Neutrophils
)
markers5 <- c(
  # Stromal and vascular cells
  "COL1A1","COL3A1","DCN",  # Fibroblasts
  "ACTA2","TAGLN","MYH11",  # Smooth muscle / myofibroblast (MYH11 often cleaner than MYL9)
  "PECAM1","VWF","KDR",  # Endothelial
  "PROX1","PDPN","LYVE1",  # Lymphatic endothelial
  "RGS5","PDGFRB","CSPG4",  # Pericytes
  "MKI67","TOP2A","STMN1"  # Proliferative
)

markers_all <- c(markers1, markers2, markers3, markers4, markers5)

pdf("pre/AH_all_markerBubble.pdf", useDingbats=FALSE, width = 25, height = 12)
DotPlot(spatial.all, features = markers_all, assay = "Spatial", group.by = "seurat_clusters") + scale_color_distiller(palette = "RdYlBu") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
dev.off()

cluster2celltype_self <- c("0"="Naive_B", "1"="Naive_B",
                           "2"="CD4_T", "3"="CD4_T",
                           "4"="GC_B", "5"="Memory_B",
                           "6"="Memory_B", "7"="GC_B",
                           "8"="Tfh", "9"="Memory_B",
                           "10"="Memory_B", "11"="CD8_T",
                           "12"="Memory_B", "13"="Plasma_cell",
                           "14"="GC_B", "15"="Fibroblast",
                           "16"="Treg", "17"="Naive_B",
                           "18"="GC_B", "19"="Macrophage",
                           "20"="Macrophage", "21"="pDC",
                           "22"="Unknown", "23"="Glandular_epithelium",
                           "24"="Fibroblast", "25"="Fibroblast",
                           "26"="Tfh", "27"="CD4_T",
                           "28"="Endothelial", "29"="Ciliated_epithelium",
                           "30"="NK", "31"="Unknown",
                           "32"="CD8_T", "33"="cDC1",
                           "34"="GC_B", "35"="Naive_B",
                           "36"="Mast_cell")

spatial.all[['celltype']]<-unname(cluster2celltype_self[spatial.all@meta.data$seurat_clusters])

table(spatial.all@meta.data$celltype)

dim(spatial.all)

spatial.all$celltype <- factor(spatial.all$celltype, levels = c("CD4_T","CD8_T","cDC1","Ciliated_epithelium","Endothelial","Fibroblast",
                              "GC_B","Glandular_epithelium","Macrophage","Mast_cell","Memory_B","Naive_B","NK","pDC","Plasma_cell","Tfh","Treg","Unknown"))
Idents(spatial.all) <- "celltype"

str(spatial.all@meta.data)
table(spatial.all$celltype)
table(spatial.all$TH_meta)

pal <- c("#7FB2D5","#B3CDE3","#F7D8B7","#17c3b2","#BEDADA","#8DBA8F","#6A3D9A","#00A087","#F39B7FFF","#FFD569",
        "#984EA3","#BC9DCC","#57C3F3","#F5B47F","#D37295","#F4CAE4","#4581B2","#F2F2F2")

pal.patient <- c("#0072B599","#DC0000FF","#925E9F99","#DF8F4499","#4DBBD599","#008B4599")

# Metadata UMAPs
pdf("pre/AH_all_INFO_UMAP.pdf", width = 8, height = 6)
DimPlot(spatial.all, reduction = 'umap', group.by = 'gender', pt.size = 0.8, cols = c("#E7C1D7","#AFB0D0"), raster = FALSE) + ggtitle("Gender") + theme_global
DimPlot(spatial.all, reduction = 'umap', group.by = 'Patient', pt.size = 0.8, cols = pal.patient, raster = FALSE) + ggtitle("Patient") + theme_global
DimPlot(spatial.all, reduction = 'umap', group.by = 'TH_meta', pt.size = 0.8, cols = c("#9EC4BE","#FAC795"), raster = FALSE) + ggtitle("TH_meta") + theme_global
dev.off()

pdf("pre/AH_all_INFO_UMAP2.pdf", width = 8, height = 6)
DimPlot(spatial.all, reduction = 'umap', group.by = 'gender', pt.size = 0.8, cols = c("#E7C1D7","#AFB0D0"), raster = TRUE) + ggtitle("Gender") + theme_global
DimPlot(spatial.all, reduction = 'umap', group.by = 'Patient', pt.size = 0.8, cols = pal.patient, raster = TRUE) + ggtitle("Patient") + theme_global
DimPlot(spatial.all, reduction = 'umap', group.by = 'TH_meta', pt.size = 0.8, cols = c("#9EC4BE","#FAC795"), raster = TRUE) + ggtitle("TH_meta") + theme_global
dev.off()

# Patient-level QC
pdf("pre/AH_all_QC_byPatient.pdf", width = 18, height = 5)
VlnPlot(
  spatial.all,
  features = c("nCount_Spatial", "nFeature_Spatial", "percent.mt", "percent.rb"),
  group.by = "Patient", pt.size = 0, ncol = 4, cols = pal.patient
)
dev.off()

# Cell-type marker DotPlots
library(grid)
Idents(spatial.all) <- "celltype"

markers<-c("CD3E","CD4","IL7R",         #CD4_T
          "CD8A","CD8B","GZMK",         #CD8_T
          "CLEC9A","XCR1","BATF3",      #cDC1
          "FOXJ1","TPPP3","DNAH5",      #Ciliated_epithelium
          "PECAM1","VWF","KDR",         #Endothelial
          "COL1A1","COL3A1","DCN",      #Fibroblast
          "AICDA","RGS13","MEF2B",      #GC_B
          "KRT8","KRT18","KRT19",       #Glandular_epithelium
          "C1QA","C1QB","APOE",         #Macrophage
          "TPSAB1","TPSB2","KIT",       #Mast_cell
          "CD27","TNFRSF13B","FCRL4",   #Memory_B
          "MS4A1","IGHD","TCL1A",       #Naive_B
          "NKG7","GNLY","PRF1",         #NK
          "IL3RA","TCF4","IRF7",        #pDC
          "MZB1","XBP1","JCHAIN",       #Plasma_cell
          "CD40LG","PDCD1","ICOS",       #Tfh
          "FOXP3","IL2RA","CTLA4"       #Treg
)


pdf("/mnt/data_nvme1/lzh/test/AH/atlas/pre/AH_all_DotPlot.celltype.pdf", width = 13, height = 6)
DotPlot(
  spatial.all, features = markers, cols = c("grey", "red"),
  dot.min = 0, dot.scale = 6, idents = NULL, group.by = "celltype",
  split.by = NULL, cluster.idents = FALSE, scale = TRUE,
  scale.by = "radius", scale.min = NA, scale.max = NA
) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_y_discrete(limits = levels(spatial.all$celltype)) +
  theme(
    panel.border = element_rect(color = "black"),
    panel.spacing = unit(1, "mm"),
    strip.text = element_text(margin = margin(b = 3, unit = "mm")),
    strip.placement = 'outlet',
    axis.line = element_blank()
  ) +
  labs(x = "", y = "")

DotPlot(
  spatial.all, features = markers,
  dot.min = 0, dot.scale = 6, idents = NULL, group.by = "celltype",
  split.by = NULL, cluster.idents = FALSE, scale = TRUE,
  scale.by = "radius", scale.min = NA, scale.max = NA
) +
  scale_color_distiller(palette = "RdYlBu") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_y_discrete(limits = levels(spatial.all$celltype)) +
  labs(x = "", y = "") +
  theme(panel.border = element_rect(color = "black", size = 0.5, fill = NA))
dev.off()

# Cluster and cell-type UMAPs
pdf("pre/AH_all_Celltype_UMAP.pdf", width = 8, height = 6)
DimPlot(spatial.all, reduction = 'umap', group.by = 'seurat_clusters', pt.size = 0.8, raster = FALSE,label = TRUE) + ggtitle("Cluster") + theme_global
DimPlot(spatial.all, reduction = 'umap', group.by = 'celltype', pt.size = 0.8, cols = pal, raster = FALSE) + ggtitle("Celltype") + theme_global
dev.off()

pdf("pre/AH_all_Celltype_UMAP2.pdf", width = 8, height = 6)
DimPlot(spatial.all, reduction = 'umap', group.by = 'seurat_clusters', pt.size = 0.8, raster = TRUE,label = TRUE) + ggtitle("Cluster") + theme_global
DimPlot(spatial.all, reduction = 'umap', group.by = 'celltype', pt.size = 0.8, cols = pal, raster = TRUE) + ggtitle("Celltype") + theme_global
dev.off()

# Spatial cell distributions
pdf("pre/AH_all_spatial_byPatient_celltype.pdf", width = 8, height = 6)
for (p in levels(spatial.all$Patient)) {
  message("Plotting Patient: ", p)
  obj.p <- subset(spatial.all, subset = Patient == p)
  p2 <- DimPlot(obj.p, reduction = "spatial", group.by = "celltype", pt.size = 1, cols = pal) + ggtitle(p) + theme_global
  print(p2)
}
dev.off()

pdf("pre/AH_all_spatial_byPatient_celltype2.pdf", width = 6, height = 6)
for (p in levels(spatial.all$Patient)) {
  message("Plotting Patient: ", p)
  obj.p <- subset(spatial.all, subset = Patient == p)
  p2 <- DimPlot(obj.p, reduction = "spatial", group.by = "celltype", pt.size = 1, cols = pal) + ggtitle(p) + theme_global+ NoLegend()
  print(p2)
}
dev.off()

# 08. Cell-cycle scoring and Seurat object export
str(spatial.all@meta.data)
# Score cell-cycle phases.
s.genes <- Seurat::cc.genes$s.genes
g2m.genes <- Seurat::cc.genes$g2m.genes
spatial.all <- CellCycleScoring(spatial.all,s.features = s.genes, g2m.features = g2m.genes,set.ident = T)

pdf("pre/AH_all_nCount_Spatial_UMAP.pdf", width = 8, height = 6)
FeaturePlot(spatial.all, reduction = 'umap', features = 'nCount_Spatial', pt.size = 0.8, raster = FALSE) + ggtitle("nCount_Spatial") +
            theme_global + scale_color_distiller(palette = "RdYlBu")
FeaturePlot(spatial.all, reduction = 'umap', features = 'nFeature_Spatial', pt.size = 0.8, raster = FALSE) + ggtitle("nFeature_Spatial") +
            theme_global + scale_color_distiller(palette = "RdYlBu")
dev.off()

Idents(spatial.all) <- "celltype"

median(spatial.all$nCount_Spatial)
median(spatial.all$nFeature_Spatial)

saveRDS(spatial.all, "/mnt/data_nvme1/lzh/test/AH/AH_all_spatial.fli.RDS")

# 09. Donor-level PCA
setwd("/mnt/data_nvme1/lzh/test/AH/atlas")
objd.all <- readRDS("/mnt/data_nvme1/lzh/test/AH/AH_all_spatial.fli.RDS")
str(objd.all@meta.data)

# Aggregate expression by patient for group-level PCA.
objd.all$Group2 <- ifelse(objd.all$TH_meta == "N", "AH", "ATH-A")
group_colors <- c("AH" = "#F7D8B7", "ATH-A" = "#E7C1D7")

# Average expression within each patient and select variable genes.
library(Seurat)
library(dplyr)
library(ggplot2)
library(matrixStats)

run_sample_pca <- function(seu, top_n = 3000, assay = "Spatial", slot = "data",
                           group_colors = NULL) {
  expr <- GetAssayData(seu, assay = assay, slot = slot)
  sample_ids <- seu@meta.data$Patient

  # Calculate mean expression per patient.
  sample_levels <- unique(as.character(sample_ids))
  avg_list <- lapply(sample_levels, function(s) {
    cells_use <- colnames(seu)[sample_ids == s]
    Matrix::rowMeans(expr[, cells_use, drop = FALSE])
  })
  names(avg_list) <- sample_levels
  avg_expr <- do.call(cbind, avg_list)
  colnames(avg_expr) <- sample_levels

  avg_expr <- as.matrix(avg_expr)

  # Select genes with highest between-patient variance.
  gene_var <- matrixStats::rowVars(avg_expr)
  gene_var[is.na(gene_var)] <- 0
  top_n <- min(top_n, nrow(avg_expr))
  top_genes <- names(sort(gene_var, decreasing = TRUE))[1:top_n]
  avg_expr_top <- t(avg_expr[top_genes, , drop = FALSE])  # Rows: patients; columns: genes

  pca_res <- prcomp(avg_expr_top, center = TRUE, scale. = TRUE)
  pca_df <- as.data.frame(pca_res$x[, 1:2])
  pca_df$Patient <- rownames(pca_df)

  # Join patient-level group metadata.
  meta_df <- seu@meta.data %>% select(Patient, Group2) %>% distinct()
  meta_df$Patient <- as.character(meta_df$Patient)
  pca_df$Patient <- as.character(pca_df$Patient)
  pca_df <- left_join(pca_df, meta_df, by = "Patient")
  pca_df$Group2 <- factor(pca_df$Group2, levels = names(group_colors))

  # Variance explained by each principal component
  percent_var <- round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2)), 1)
  p <- ggplot(pca_df, aes(x = PC1, y = PC2, fill = Group2, color = Group2)) +
    stat_ellipse(aes(fill = Group2), geom = "polygon", alpha = 0.25, color = NA, level = 0.95) +
    geom_point(shape = 21, size = 5, stroke = 1.1, alpha = 0.9) +
    ggrepel::geom_text_repel(aes(label = Patient), size = 4, color = "black",
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

# Run donor-level PCA.
res_group <- run_sample_pca(objd.all, group_colors = group_colors)

# Export the PCA figure.
pdf("pre/sample_PCA_by_treatment.pdf", width = 6, height = 5)
print(res_group$plot)
dev.off()
