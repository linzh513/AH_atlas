# Paired adenoid–tonsil scRNA-seq differential analysis
.libPaths("/home/lzh/miniconda3/envs/singlecell/lib/R/library")
rm(list = ls())
setwd("/mnt/data_nvme1/lzh/test/AH/sc_atlas")
set.seed(1234)
library(Seurat)
library(dplyr)
library(tidyr)
library(limma)
library(Matrix)
library(magrittr)
library(stringr)
library(patchwork)
library(clustree)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(forcats)
theme_global <- theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank()) + theme(aspect.ratio = 1)
objd.all <- readRDS("/mnt/data_nvme1/lzh/test/AH/sc_atlas/AH_all_rna.fli.RDS")
str(objd.all@meta.data)
pal <- c("#7FB2D5","#B3CDE3","#BEDADA","#8DBA8F","#6A3D9A","#00A087","#F39B7FFF",
        "#984EA3","#BC9DCC","#57C3F3","#F5B47F","#D37295","#F4CAE4","#4581B2")
library(ggpie)
library(ggprism)
table(objd.all$celltype, objd.all$Group2)
table(objd.all$celltype, objd.all$Group)

# 01. Cell-type composition
Df_raw <- data.frame(Group = objd.all$celltype)
Df_count <- as.data.frame(table(objd.all$celltype))
colnames(Df_count) <- c("Group", "count")

pdf("diff/Celltype_pie.pdf", width = 10, height = 10)
ggpie(data = Df_raw,
            group_key = "Group",
            count_type = "full",
            label_info = c('group', 'ratio'),
            label_type = "horizon",
            label_split = NULL,
            label_size = 4,
            label_pos = "in",
            label_threshold = 5,
            fill_color = pal) +
      theme_void(base_size = 20) +
      theme(legend.position = "none")
dev.off()

pdf("diff/Celltype_pie2.pdf", width = 8, height = 8)
ggdonut(data = Df_raw,
              group_key = "Group",
              count_type = "full",
              fill_color = pal,
              label_info = c('group', 'ratio'),
              label_split = " ",
              label_color = "black", label_type = "horizon",
              label_pos = "out",
              label_gap = 0.2, label_threshold = 0.1, label_size = 4.75,
              border_color = "white", border_size = 0.5,
              r0 = 0.8,
              r1 = 2.75,
              donut.label = TRUE,
              donut.label.size = 5.25,
              donut.label.color = "darkblue") +
       theme(legend.position = "none")
dev.off()

# Cell-type composition by tissue group
table(objd.all$celltype, objd.all$Group2)
num2<-table(objd.all$celltype, objd.all$Group2)
write.csv(num2,"diff/num2.csv")
num2<-as.data.frame((read.csv("diff/num2.csv",header = T,row.names = 1)))
num2
pct <- num2 %>%
  mutate(across(everything(), ~ .x / sum(.x, na.rm = TRUE)))
pct2 <- as.data.frame(t(pct))
pct2$Group <- rownames(pct2)
dat2 <- pivot_longer(data = pct2, cols = c("CD4_T","CD8_T","Endothelial","Fibroblast","GC_B","Glandular_epithelium",
                              "Macrophage","Memory_B","Naive_B","NK","pDC","Plasma_cell","Tfh","Treg"),
                     names_to = "celltype", values_to = "Percentage")
dat2$Group <- factor(dat2$Group, levels = c("AH","ATH.A","ATH.T"))
dat2$celltype <- factor(dat2$celltype, levels=c("CD4_T","CD8_T","Endothelial","Fibroblast","GC_B","Glandular_epithelium",
                              "Macrophage","Memory_B","Naive_B","NK","pDC","Plasma_cell","Tfh","Treg"))
library(ggplot2)
dat2$label <- paste0(sprintf("%.2f", as.numeric(dat2$Percentage) * 100), "%")
dat2$percent <- 100*(dat2$Percentage)

pdf("diff/percent.celltype.pdf", width = 8, height = 5)
ggplot(dat2, aes(x = Group, y = percent, fill = celltype))+
  geom_col(width = 0.2)+
  geom_bar(stat="identity", position = position_stack(), width = 0.75)+
  scale_fill_manual(values = pal)+
  scale_y_continuous(labels = scales::percent_format(scale = 1))+
  labs(x="",y="Percentage",fill="")+
  geom_text(aes(label = label),vjust=3,size=0,color="black")+theme_classic()+
  theme(legend.position = "top",
        legend.text = element_text(size=12), legend.title = element_text(size = 14),
        axis.text.x = element_text(colour="black",size = 12,),
        axis.text.y = element_text(colour="black",size = 12),
        axis.title = element_text(size=14),
        legend.key.size = unit(10, "pt")) + coord_flip()
dev.off()

# Cell-type composition by individual sample
table(objd.all$celltype, objd.all$Patient)
table(objd.all$Patient, objd.all$Tissue)  # AH:97,101,103; ATH:99,100,102
objd.all$sample <- factor(objd.all$sample, levels = c("AH97","AH101","AH103","AH99","AH100","AH102","TH99","TH100","TH102"))
num3 <- table(objd.all$sample, objd.all$celltype)
write.csv(num3,"diff/num3.csv")
num3 <- as.data.frame.matrix(num3)
num3
pct <- num3 / rowSums(num3)
pct$sample <- rownames(pct)
dat2 <- pivot_longer(pct, cols = -sample,
  names_to = "celltype", values_to = "Percentage")
dat2$Percentage[is.na(dat2$Percentage)] <- 0
dat2$percent <- dat2$Percentage * 100
dat2$celltype <- factor(
  dat2$celltype,
  levels = c(
    "CD4_T","CD8_T","Endothelial","Fibroblast",
    "GC_B","Glandular_epithelium","Macrophage","Memory_B",
    "Naive_B","NK","pDC","Plasma_cell","Tfh","Treg"
  )
)
dat2$sample <- factor(dat2$sample, levels = c("AH97","AH101","AH103","AH99","AH100","AH102","TH99","TH100","TH102"))

pdf("diff/percent.celltype.patient.pdf", width = 8, height = 6)
ggplot(dat2, aes(x = sample, y = percent, fill = celltype)) +
  geom_col(width = 0.75) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  scale_fill_manual(values = pal) +
  labs(x = "", y = "Cell proportion", fill = "") +
  theme_classic() +
  theme(legend.position = "top",
    axis.text = element_text(color = "black", size = 12),
    axis.title = element_text(size = 14)
  )
dev.off()

# 02. Sample-level cell-type proportion comparisons
# Use samples as independent units in the original proportion analysis.
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
num <- table(objd.all$orig.ident, objd.all$celltype)
num <- as.data.frame.matrix(num)
pct <- num / rowSums(num)
pct <- as.data.frame(pct)
pct$orig.ident <- rownames(pct)
pct$group <- c("AH", "ATH-A", "ATH-A", "AH", "ATH-A", "AH", "ATH-T", "ATH-T", "ATH-T")
pct$group <- factor(pct$group, levels = c("AH", "ATH-A", "ATH-T"))
pct <- pct %>% arrange(group, orig.ident)
write.csv(pct, "diff/sample_cell_proportion.csv", row.names = FALSE)
meta_cols <- c("group", "orig.ident")
celltype_cols <- setdiff(colnames(pct), meta_cols)
df_long <- pct %>% pivot_longer(cols = all_of(celltype_cols), names_to = "celltype", values_to = "Proportion")
df_long$Proportion <- as.numeric(df_long$Proportion)
df_long$group <- factor(df_long$group, levels = c("AH", "ATH-A", "ATH-T"))
celltype_levels <- c("CD4_T","CD8_T","Endothelial","Fibroblast",
                      "GC_B","Glandular_epithelium","Macrophage","Memory_B",
                      "Naive_B","NK","pDC","Plasma_cell","Tfh","Treg")
group_colors <- c("AH" = "#F7D8B7", "ATH-A" = "#E7C1D7", "ATH-T" = "#AFB0D0")

# Convert P values to significance symbols.
p_to_symbol <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  return("ns")
}

# Plot sample-level means and donor points
plot_one_celltype <- function(dat, celltype_name, fill_colors = NULL) {
  dat <- dat %>% filter(!is.na(Proportion), !is.na(group))
  dat$Proportion <- as.numeric(dat$Proportion)
  dat$group <- droplevels(dat$group)
  stat_df <- dat %>% group_by(group) %>%
    summarise(mean_val = mean(Proportion, na.rm = TRUE),
              sd_val = sd(Proportion, na.rm = TRUE),
              n = n(),
              se_val = sd_val / sqrt(n),
              .groups = "drop")
  p <- ggplot(stat_df, aes(x = group, y = mean_val, fill = group)) +
    geom_col(width = 0.6, color = "black", alpha = 0.6) +
    geom_errorbar(aes(ymin = mean_val - se_val, ymax = mean_val + se_val), width = 0.18, linewidth = 0.6) +
    geom_point(data = dat, aes(x = group, y = Proportion, fill = group),
               shape = 21, color = "black", size = 3,
               position = position_jitter(width = 0.12, height = 0), inherit.aes = FALSE) +
    labs(title = celltype_name, x = "Group", y = "Cell proportion") +
    theme_classic(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
  if (!is.null(fill_colors)) p <- p + scale_fill_manual(values = fill_colors)
  group_levels <- levels(dat$group)
  group_levels <- group_levels[group_levels %in% unique(as.character(dat$group))]
  if (length(group_levels) >= 2) {
    comps <- combn(group_levels, 2, simplify = FALSE)
    y_max <- max(dat$Proportion, stat_df$mean_val + stat_df$se_val, na.rm = TRUE)
    y_step <- y_max * 0.12
    if (is.na(y_step) || y_step == 0) y_step <- 0.05
    anno_df <- lapply(seq_along(comps), function(i) {
      g1 <- comps[[i]][1]
      g2 <- comps[[i]][2]
      x1 <- which(group_levels == g1)
      x2 <- which(group_levels == g2)
      v1 <- dat$Proportion[dat$group == g1]
      v2 <- dat$Proportion[dat$group == g2]
      pval <- tryCatch(wilcox.test(v1, v2)$p.value, error = function(e) NA)
      data.frame(x1 = x1, x2 = x2, y = y_max + i * y_step, label = sprintf("%.2f", pval, 2))
    }) %>% bind_rows()
    if (nrow(anno_df) > 0) {
      p <- p +
        geom_segment(data = anno_df, aes(x = x1, xend = x2, y = y, yend = y), inherit.aes = FALSE) +
        geom_text(data = anno_df, aes(x = (x1 + x2) / 2, y = y + y_step * 0.08, label = label),
                  inherit.aes = FALSE, size = 5) +
        coord_cartesian(ylim = c(0, max(anno_df$y) + y_step * 0.3))
    }
  }
  return(p)
}

pdf("diff/celltype.pdf", width = 5.5, height = 6)
for (ct in celltype_levels) {
  sub_ct <- df_long %>% filter(celltype == ct)
  if (nrow(sub_ct) == 0) next
  p <- plot_one_celltype(sub_ct, celltype_name = ct, fill_colors = group_colors)
  print(p)
}
dev.off()

# 03. Differential gene expression
objd.all3 <- subset(objd.all, Group == "ATH")
library(EnhancedVolcano)
Idents(objd.all3) <- "celltype"
DefaultAssay(objd.all3) <- "RNA"
logFCfilter=0.50
adjPvalFilter=1e-5

# Marker genes across cell types
diff_gene.all <- FindAllMarkers(objd.all3, min.pct = 0.5, logfc.threshold = 0.25, group.by = "celltype",
                              only.pos = FALSE, test.use="wilcox")
sig_gene.all <- subset(diff_gene.all, p_val_adj < adjPvalFilter & abs(avg_log2FC) > logFCfilter)
sig_gene.all <- sig_gene.all %>%  group_by(cluster) %>% arrange(desc(avg_log2FC), .by_group = TRUE)
saveRDS(diff_gene.all,"diff2/diff_gene.all.ct.RDS")
saveRDS(sig_gene.all,"diff2/sig_gene.all.ct.RDS")
write.csv(sig_gene.all,"diff2/sig_gene.all.ct.csv", row.names = FALSE)

# Pooled tonsil-versus-adenoid comparison
diff_gene.all.ct <- FindMarkers(objd.all3, min.pct = 0.1, logfc.threshold = 0.25, group.by = "Tissue",
                              ident.1 = "Tonsil", ident.2 = "Adenoid", test.use="wilcox")
diff_gene.all.ct <- diff_gene.all.ct %>% arrange(desc(avg_log2FC), .by_group = TRUE)
dim(diff_gene.all.ct)
saveRDS(diff_gene.all.ct,"diff2/diff_gene.all.ct.RDS")
write.csv(diff_gene.all.ct,"diff2/diff_gene.all.ct.csv", row.names = TRUE)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(tidyverse)
diff_gene.all.ct <- readRDS("diff2/diff_gene.all.ct.RDS")
diff_gene.all.ct$names <- rownames(diff_gene.all.ct)
diff_gene.all.ct <- diff_gene.all.ct %>% mutate(Difference = pct.1-pct.2)
diff_gene.all.ct$group <- ifelse(diff_gene.all.ct$avg_log2FC>0,'UP','Down')
diff_gene.all.ct$group2 <- ifelse((abs(diff_gene.all.ct$avg_log2FC) < log2(1.25) | diff_gene.all.ct$p_val_adj > 1e-5),
                           'ns', diff_gene.all.ct$group)
top.up <- diff_gene.all.ct  %>% top_n(n = 15, wt = avg_log2FC) %>% row.names()
top.dn <- diff_gene.all.ct  %>% top_n(n = -15, wt = avg_log2FC) %>% row.names()
label <- c(top.up, top.dn)
diff_gene.all.ct$label <- ifelse(diff_gene.all.ct$names %in% label, diff_gene.all.ct$names,"")

pdf('diff2/degVolcano.all.ct.pdf',height=6, width=6.5)
ggplot(diff_gene.all.ct, aes(x=Difference, y=avg_log2FC))+
  geom_point(size = 2, aes(color = group2))+
  scale_color_manual(values = c("#1176bd","#7b7f80","#cb5a54"))+
  geom_label_repel(aes(label = label),fill = NA,
                   segment.size = 0.2, size = 2.75, max.overlaps = 60,
                   fontface="bold", label.size = NA) +
  geom_vline(xintercept = 0, linetype = 2)+
  geom_hline(yintercept = 0, linetype = 2)+
  theme_classic()+
  xlab("Percentage Difference")+
  ylab("Log-Fold Change")
dev.off()

# Cell-type-specific tonsil-versus-adenoid comparisons
celltypes <- unique(objd.all3$celltype)
celltypes
markers_list <- lapply(celltypes, function(ct) {
  celltype_data <- subset(objd.all3, subset = celltype == ct)
  markers <- FindMarkers(celltype_data, min.pct = 0.2, logfc.threshold = 0.25, group.by = "Tissue", ident.1 = "Tonsil", ident.2 = "Adenoid", test.use="wilcox")
  return(markers)
})
names(markers_list) <- celltypes
markers_df <- do.call(
  rbind,
  lapply(names(markers_list), function(ct) {
    df <- markers_list[[ct]]
    df$gene <- rownames(df)
    df$celltype <- ct
    rownames(df) <- NULL
    return(df)
  })
)
sig.markers_df <- subset(markers_df, p_val_adj < 1e-3 & abs(avg_log2FC) > 0.5)
sig.markers_df <- sig.markers_df %>%  group_by(celltype) %>% arrange(desc(avg_log2FC), .by_group = TRUE)
table(sig.markers_df$celltype)
write.csv(sig.markers_df, file = "diff2/sig.markers_df.csv", row.names = FALSE)
saveRDS(markers_list, file = "diff2/markers_list.RDS")
saveRDS(markers_df, file = "diff2/markers_df.RDS")

# 04. Multi-cell-type volcano plots
df <- read.csv("diff2/sig.markers_df.csv", header = TRUE)
df <- df %>% mutate(chr = as.character(celltype),
    type = if_else(avg_log2FC > 0, "UP_Highly", "Down_Highly"))
pal <- c("#7FB2D5","#B3CDE3","#8DBA8F","#6A3D9A","#00A087","#F39B7FFF",
        "#984EA3","#BC9DCC","#57C3F3","#F5B47F","#D37295","#F4CAE4","#4581B2")
df$celltype <- factor(df$celltype, levels = unique(df$celltype))
bg_df <- tibble(celltype = levels(df$celltype),
  xmin = seq_along(levels(df$celltype)) - 0.48,
  xmax = seq_along(levels(df$celltype)) + 0.48,
  ymin = -0.5, ymax = 0.5)
bg_vertical <- df %>%
  group_by(celltype) %>%
  summarise(ymin = min(avg_log2FC, na.rm = TRUE) - 0.1,
            ymax = max(avg_log2FC, na.rm = TRUE) + 0.1,
    .groups = "drop") %>%
  mutate(xmin = seq_along(celltype) - 0.48,
         xmax = seq_along(celltype) + 0.48)
label_df <- df %>%
  filter(p_val_adj < 1e-3, abs(avg_log2FC) > 1) %>%
  group_by(celltype) %>%
  group_modify(~ bind_rows(
    slice_max(.x, avg_log2FC, n = 5, with_ties = FALSE),
    slice_min(.x, avg_log2FC, n = 5, with_ties = FALSE)
  ) %>% distinct(gene, .keep_all = TRUE)) %>%  # Remove duplicate genes
  ungroup()
library(ggplot2)
library(ggprism)

pdf("diff2/volcano.pdf", width = 15, height = 8)
ggplot(df, aes(celltype, avg_log2FC, color = type)) +
  geom_rect(data = bg_vertical,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, fill = "grey95") +
  geom_rect(data = bg_df,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                fill = celltype),
            inherit.aes = FALSE) +
  geom_jitter(stroke = 0) +
  geom_text_repel(
    data = label_df,
    aes(celltype, avg_log2FC, label = gene),
    size = 2.5, color = "black", box.padding = 0.2) +
  geom_text(aes(celltype, 0, label = chr),
            size = 3, color = "white", show.legend = FALSE) +
  scale_fill_manual(values = pal, guide = "none") +
  scale_color_manual(values = c("#0073C2FF", "#EE0000FF")) +
  scale_y_continuous(
    limits = c(-4.5, 5.5),
    breaks = c(-4.5, -2.5, -0.5, 0, 1.5, 3.5, 5.5),
    guide = "prism_offset_minor") +
  labs(x = NULL, y = "average log2FC") +
  guides(color = guide_legend(
    override.aes = list(size = 5, shape = 19))) +
  theme_prism(base_line_size = 0.3) +
  theme(
    panel.background = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x  = element_blank(),
    axis.text.y  = element_text(size = 10),
    axis.title   = element_text(size = 11),
    legend.position = c(0.08, 0.9),
    legend.key = element_blank(),
    legend.background = element_blank(),
    legend.text = element_text(margin = margin(l = 0),size = 8)
)
dev.off()

# Differentially expressed gene counts
df <- read.csv("diff2/sig.markers_df.csv", header = TRUE)
df <- df %>% mutate(chr = as.character(celltype),
    type = if_else(avg_log2FC > 0, "UP_Highly", "Down_Highly"))
ddff <- df %>% filter(type == "Down_Highly") %>%
  group_by(celltype) %>% summarise(n = n_distinct(gene)) %>%
  arrange(desc(n), .by_group = TRUE)
library(ggplot2)
ddff$celltype <- factor(ddff$celltype, levels = rev(ddff$celltype))

pdf("diff2/ct_n_deg_dn.pdf", width = 6, height = 4)
ggplot(ddff, aes(x = n, y = celltype)) +
  geom_segment(aes(x = 0, xend = n, yend = celltype),
               color = "grey75", linewidth = 0.8) +
  geom_point(size = 5.5, color = "#4581B2") +  # "#CB5C5B", "#4581B2"
  geom_text(aes(label = n), hjust = 0.5, size = 3, color = "black") +
  labs(x = "Number of down-regulated genes", y = NULL) +
  expand_limits(x = max(ddff$n) * 1.1) +
  theme(axis.text = element_text(size = 10), axis.title   = element_text(size = 11))+
  theme_classic()
dev.off()

# 05. UpSet analysis of shared genes
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(UpSetR)
library(grid)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(grid)
sig.single <- read.csv("diff2/sig.markers_df.csv")
sig.df <- sig.single %>%
  filter(avg_log2FC > 0) %>%
  select(gene, celltype) %>%
  distinct()
gene_counts <- sig.df %>%
  group_by(gene) %>%
  summarise(n_celltypes = n(), .groups = "drop")
sig.df.filtered <- sig.df %>%
  inner_join(gene_counts, by = "gene") %>%
  filter(n_celltypes > 2) %>%
  select(gene, celltype)
sig <- sig.df.filtered %>%
  group_by(celltype) %>%
  summarise(gene = list(gene), .groups = "drop") %>%
  deframe()
sig <- lapply(sig, na.omit)
sig <- sig[sapply(sig, length) > 0]
pal <- pal
labels_all <- c("CD4_T","CD8_T","Fibroblast",
                "GC_B","Glandular_epithelium","Macrophage","Mast_cell","Memory_B",
                "Naive_B","NK","pDC","Plasma_cell","Tfh","Treg")
existing_celltypes <- names(sig)
idx <- match(existing_celltypes, labels_all)
pal <- pal[idx]
labels <- labels_all[idx]
m <- make_comb_mat(sig)
cs <- comb_size(m)

pdf("diff2/upset.sc2.pdf", width = 12, height = 6.5)
ht <- draw(
  UpSet(
    m,
    top_annotation = upset_top_annotation(
      m,
      ylim = c(0, max(cs) * 1.1),
      annotation_name_rot = 90,
      point.size = 1.25,
      line.size = 0.5,
      annotation_name_side = "left",
      height = unit(3, "cm")
    ),
    bg_col = rev(pal[seq_len(length(sig))]),
    set_order = labels
  )
)
od <- column_order(ht)
decorate_annotation("intersection_size", {
  grid.text(
    cs[od],
    x = seq_along(cs),
    y = unit(cs[od], "native") + unit(2, "pt"),
    default.units = "native",
    just = "left",
    gp = gpar(fontsize = 10),
    rot = 90
  )
})
dev.off()

binary_to_celltypes <- function(bin_string, labels) {
  bits <- strsplit(bin_string, "")[[1]]
  paste(labels[bits == "1"], collapse = ",")
}

all_comb <- comb_name(m)
all_size <- comb_size(m)
comb_ordered <- all_comb[od]
size_ordered <- all_size[od]
res_df <- bind_rows(lapply(seq_along(comb_ordered), function(i) {
  comb <- comb_ordered[i]
  genes <- extract_comb(m, comb)
  data.frame(
    plot_order = i,
    combination = comb,
    celltypes = binary_to_celltypes(comb, labels),
    gene_count = length(genes),
    genes = paste(sort(genes), collapse = ","),
    stringsAsFactors = FALSE
  )
}))
write.csv(res_df, "diff2/upset_intersection_genes_by_plot_order.csv", row.names = FALSE)
up.genes <- sig.single %>%
  filter(avg_log2FC > 0) %>%
  select(gene, celltype) %>%
  distinct()
gene.freq <- up.genes %>%
  count(gene, name = "n_celltypes") %>%
  arrange(desc(n_celltypes))
gene.freq <- gene.freq %>%
  left_join(
    sig.single %>%
      filter(avg_log2FC > 0) %>%
      group_by(gene) %>%
      summarise(mean_log2FC = mean(avg_log2FC), .groups = "drop"),
    by = "gene"
  ) %>%
  arrange(desc(n_celltypes), desc(mean_log2FC))
head(gene.freq, 15)

# Upregulated and downregulated gene counts
library(dplyr)
library(ggprism)
gene.ratio <- sig.single %>%
  select(gene, celltype, avg_log2FC) %>%
  distinct() %>%  # Avoid double counting
  group_by(celltype) %>%
  summarise(
    up_genes   = sum(avg_log2FC > 0),
    down_genes = sum(avg_log2FC < 0),
    ratio_up_down = ifelse(down_genes == 0,
                           NA,
                           up_genes / down_genes),
    .groups = "drop"
  ) %>%
  arrange(desc(ratio_up_down))
gene.ratio
dat4 <- gene.ratio %>%
  mutate(Group = celltype,
    high = pmax(up_genes, down_genes), low  = pmin(up_genes, down_genes),
    high_color = ifelse(up_genes > down_genes, "#CB5C5B", "#4581B2"),
    low_color  = ifelse(up_genes > down_genes, "#4581B2", "#CB5C5B")
  )

# Visualize upregulated and downregulated genes
pdf("diff2/up_down_gene_count.fc.pdf", width =5, height = 5.5)
ggplot(dat4) +
  geom_segment(
    aes(x = Group, xend = Group, y = high, yend = low),color = "lightgrey",size = 0.75) +
  geom_point(aes(x = Group, y = high),color = dat4$high_color,size = 2.5) +
  geom_point(aes(x = Group, y = low),color = dat4$low_color,size = 2.5) +
  theme_prism(palette = "pearl",base_size = 12, base_line_size = 0.6, axis_text_angle = 90) +
  theme(legend.position = "none") +
  xlab("") +
  ylab("Number of DE genes")
dev.off()

# 06. Gene Ontology enrichment analysis
library(dplyr)
library(ggplot2)
library(enrichR)
library(stringr)
input_file <- "/mnt/data_nvme1/lzh/test/AH/sc_atlas/diff2/sig_gene.all.ct.csv"
outdir <- "/mnt/data_nvme1/lzh/test/AH/sc_atlas/diff2/GO_BP_by_cluster"
csv_dir <- file.path(outdir, "csv")
rds_dir <- file.path(outdir, "rds")
dir.create(csv_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(rds_dir, showWarnings = FALSE, recursive = TRUE)
deg <- read.csv(input_file, check.names = FALSE)
db_go <- "GO_Biological_Process_2025"
dbs <- enrichR::listEnrichrDbs()
if(!db_go %in% dbs$libraryName){
  stop("Database not available: ", db_go)
}

run_enrichr_safe <- function(gene_vec, db_name){
  gene_vec <- unique(na.omit(gene_vec))
  if(length(gene_vec) == 0) return(data.frame())
  res <- tryCatch(
    enrichr(gene_vec, databases = db_name),
    error = function(e) NULL
  )
  if(is.null(res) || !db_name %in% names(res)) return(data.frame())
  as.data.frame(res[[db_name]])
}

process_enrichr_result <- function(res_df, top_n = 10, direction = "Up"){
  if(is.null(res_df) || nrow(res_df) == 0) return(data.frame())
  if("Adjusted.P.value" %in% colnames(res_df)){
    res_df <- res_df %>% arrange(Adjusted.P.value)
    res_df$plot_p <- res_df$Adjusted.P.value
  } else if("P.value" %in% colnames(res_df)){
    res_df <- res_df %>% arrange(P.value)
    res_df$plot_p <- res_df$P.value
  } else {
    return(data.frame())
  }
  res_df <- res_df[1:min(top_n, nrow(res_df)), , drop = FALSE]
  res_df$plot_p[res_df$plot_p == 0] <- 1e-300
  res_df$LogP <- -log10(res_df$plot_p)
  res_df$Direction <- direction
  res_df$PlotValue <- ifelse(direction == "Up", res_df$LogP, -res_df$LogP)
  res_df
}

plot_bidirectional_enrich <- function(plot_df, cluster_name){
  if(is.null(plot_df) || nrow(plot_df) == 0){
    return(
      ggplot() +
        annotate("text", x = 0, y = 0,
                 label = paste0(cluster_name, "\nNo enriched terms")) +
        theme_void()
    )
  }
  dt <- plot_df %>%
    mutate(
      Description = stringr::str_remove(Term, "\\s*\\(GO:\\d+\\)$"),
      Cluster = ifelse(Direction == "Up", "up", "down"),
      pvalue = plot_p
    ) %>%
    group_by(Cluster) %>%
    arrange(pvalue, .by_group = TRUE) %>%
    slice_head(n = 10) %>%
    ungroup()
  dt$pvalue[dt$pvalue == 0] <- 1e-300
  dt$pvalue_loc <- -log10(dt$pvalue)
  dt$pvalue_loc <- ifelse(dt$Cluster == "up", dt$pvalue_loc, -dt$pvalue_loc)
  dt <- bind_rows(
  dt %>%
    filter(Cluster == "up") %>%
    arrange(desc(-log10(pvalue))),
  dt %>%
      filter(Cluster == "down") %>%
      arrange(-log10(pvalue))
  )
  dt$Description_plot <- ifelse(
    dt$Cluster == "up",
    paste0(dt$Description, " "),
    dt$Description
  )
  dt$Description_plot <- factor(
    dt$Description_plot,
    levels = rev(dt$Description_plot)
  )
  dt$label_xloc <- ifelse(dt$Cluster == "up", -0.5, 0.5)
  dt$label_hjust <- ifelse(dt$Cluster == "up", 1, 0)
  xmax <- ceiling(max(abs(dt$pvalue_loc), na.rm = TRUE) * 1.25)
  if(xmax < 5) xmax <- 5
  ggplot(dt, aes(x = pvalue_loc, y = Description_plot, fill = Cluster)) +
    geom_col(width = 0.8, color = "white", linewidth = 0.2) +
    geom_text(
      aes(
        x = label_xloc,
        y = Description_plot,
        label = Description,
        hjust = label_hjust
      ),
      size = 4,
      color = "black"
    ) +
    scale_fill_manual(values = c(up = "#AFB0D0", down = "#E7C1D7")) +
    scale_x_continuous(
      limits = c(-xmax, xmax),
      breaks = pretty(c(-xmax, xmax), n = 5),
      labels = abs
    ) +
    labs(
      title = cluster_name,
      x = expression(-Log[10]~"(p.adjust)"),
      y = NULL
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      axis.text.x = element_text(size = 14, color = "black"),
      axis.title.x = element_text(size = 15, color = "black"),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank(),
      axis.title.y = element_blank(),
      legend.position = "none"
    )
}

# Enrich GO biological processes by cell type
clusters_use <- unique(deg$cluster)
all_raw <- list()
all_plot <- list()
pdf_file <- file.path(outdir, "GO_BP_UpDown_by_cluster.pdf")
pdf(pdf_file, width = 13, height = 7)
for(cl in clusters_use){
  message("Processing cluster: ", cl)
  deg_sub <- deg %>% filter(cluster == cl)
  up_genes <- deg_sub %>%
    filter(avg_log2FC > 0) %>%
    pull(gene) %>%
    unique()
  down_genes <- deg_sub %>%
    filter(avg_log2FC < 0) %>%
    pull(gene) %>%
    unique()
  up_raw <- run_enrichr_safe(up_genes, db_go)
  down_raw <- run_enrichr_safe(down_genes, db_go)
  up_plot <- process_enrichr_result(up_raw, top_n = 10, direction = "Up")
  down_plot <- process_enrichr_result(down_raw, top_n = 10, direction = "Down")
  plot_df <- bind_rows(up_plot, down_plot)
  safe_cl <- gsub("[/\\\\:*?\"<>| ]+", "_", cl)
  write.csv(up_raw,
            file.path(csv_dir, paste0(safe_cl, "_Up_GO_BP_raw.csv")),
            row.names = FALSE)
  write.csv(down_raw,
            file.path(csv_dir, paste0(safe_cl, "_Down_GO_BP_raw.csv")),
            row.names = FALSE)
  write.csv(plot_df,
            file.path(csv_dir, paste0(safe_cl, "_UpDown_GO_BP_top10_plot.csv")),
            row.names = FALSE)
  saveRDS(
    list(
      cluster = cl,
      up_genes = up_genes,
      down_genes = down_genes,
      up_raw = up_raw,
      down_raw = down_raw,
      up_plot = up_plot,
      down_plot = down_plot,
      plot_df = plot_df
    ),
    file = file.path(rds_dir, paste0(safe_cl, "_GO_BP_enrichr.RDS"))
  )
  all_raw[[cl]] <- list(
    up_raw = up_raw,
    down_raw = down_raw
  )
  all_plot[[cl]] <- plot_df
  print(plot_bidirectional_enrich(plot_df, cl))
}

dev.off()
saveRDS(all_raw, file.path(rds_dir, "ALL_clusters_GO_BP_raw.RDS"))
saveRDS(all_plot, file.path(rds_dir, "ALL_clusters_GO_BP_top10_plot.RDS"))
message("Done.")
message("PDF saved to: ", pdf_file)
message("CSV saved to: ", csv_dir)
message("RDS saved to: ", rds_dir)
