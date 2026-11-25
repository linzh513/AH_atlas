# Spatial transcriptomics
# Part I: Original AH_Atlas_diff.r workflow and outputs.
# Part II (at the end): Additional donor-level comparisons
# Original analytical operations and outputs are preserved in Part I.
# One initial object reference supports Part II without loading another dataset.

# 01 | Setup and spatial Seurat object
.libPaths("/home/zihan/miniconda3/envs/singlecell/lib/R/library")
setwd("/mnt/data_nvme1/lzh/test/AH/atlas")
dir.create("diff", showWarnings = FALSE, recursive = TRUE)
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
library(paletteer)
library(ggplot2)
library(ggpubr)
library(ggprism)
library(RColorBrewer)
library(forcats)

theme_global <- theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank()) + theme(aspect.ratio = 1)
spatial.all <- readRDS("/mnt/data_nvme1/lzh/test/AH/AH_all_spatial.fli.RDS")
spatial.review <- spatial.all
str(spatial.all@meta.data)

pal <- c("#7FB2D5","#B3CDE3","#F7D8B7","#17c3b2","#BEDADA","#8DBA8F","#6A3D9A","#00A087","#F39B7FFF","#FFD569",
        "#984EA3","#BC9DCC","#57C3F3","#F5B47F","#D37295","#F4CAE4","#4581B2")

library(ggpie)
table(spatial.all$celltype, spatial.all$TH_meta)


# 02 | Original cell-type composition analysis
Df_raw <- data.frame(Group = spatial.all$celltype)

Df_count <- as.data.frame(table(spatial.all$celltype))
colnames(Df_count) <- c("Group", "count")

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

table(spatial.all$celltype, spatial.all$TH_meta)
num2<-table(spatial.all$celltype, spatial.all$TH_meta)
write.csv(num2,"diff/num2.csv")
num2<-as.data.frame((read.csv("diff/num2.csv",header = T,row.names = 1)))
num2

pct <- num2 %>%
  mutate(across(everything(), ~ .x / sum(.x, na.rm = TRUE)))

pct2 <- as.data.frame(t(pct))

pct2$Group <- rownames(pct2)
dat2 <- pivot_longer(data = pct2, cols = c("CD4_T","CD8_T","cDC1","Ciliated_epithelium","Endothelial","Fibroblast",
                              "GC_B","Glandular_epithelium","Macrophage","Mast_cell","Memory_B","Naive_B","NK","pDC","Plasma_cell","Tfh","Treg"),
                     names_to = "celltype", values_to = "Percentage")

dat2$Group <- factor(dat2$Group, levels = c("N","Y"))
dat2$celltype <- factor(dat2$celltype, levels=c("CD4_T","CD8_T","cDC1","Ciliated_epithelium","Endothelial","Fibroblast",
                              "GC_B","Glandular_epithelium","Macrophage","Mast_cell","Memory_B","Naive_B","NK","pDC","Plasma_cell","Tfh","Treg"))

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

table(spatial.all$celltype, spatial.all$Patient)
table(spatial.all$Patient, spatial.all$TH_meta)
num3 <- table(spatial.all$Patient, spatial.all$celltype)
num3 <- as.data.frame.matrix(num3)
num3 <- num3[,-18]
num3

pct <- num3 / rowSums(num3)

pct$Patient <- rownames(pct)

dat2 <- pivot_longer(pct, cols = -Patient,
  names_to = "celltype", values_to = "Percentage")

dat2$Percentage[is.na(dat2$Percentage)] <- 0
dat2$percent <- dat2$Percentage * 100

dat2$celltype <- factor(
  dat2$celltype,
  levels = c(
    "CD4_T","CD8_T","cDC1","Ciliated_epithelium","Endothelial","Fibroblast",
    "GC_B","Glandular_epithelium","Macrophage","Mast_cell","Memory_B",
    "Naive_B","NK","pDC","Plasma_cell","Tfh","Treg"
  )
)

pdf("diff/percent.celltype.patient.pdf", width = 8, height = 6)
ggplot(dat2, aes(x = Patient, y = percent, fill = celltype)) +
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

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)


# 03 | Original cell-type proportion comparisons
num <- table(spatial.all$Patient, spatial.all$celltype)
num <- as.data.frame.matrix(num)

pct <- num / rowSums(num)
pct <- as.data.frame(pct)
pct$Patient <- rownames(pct)
pct$group <- c("ATH", "AH", "AH", "ATH", "ATH", "AH")

pct$group <- factor(pct$group, levels = c("AH", "ATH"))
pct <- pct %>% arrange(group, Patient)

write.csv(pct, "diff/sample_cell_proportion.csv", row.names = FALSE)

meta_cols <- c("group", "Patient")
celltype_cols <- setdiff(colnames(pct), meta_cols)
df_long <- pct %>% pivot_longer(cols = all_of(celltype_cols), names_to = "celltype", values_to = "Proportion")

df_long$Proportion <- as.numeric(df_long$Proportion)
df_long$group <- factor(df_long$group, levels = c("AH", "ATH"))
celltype_levels <- c("CD4_T","CD8_T","cDC1","Ciliated_epithelium","Endothelial","Fibroblast",
                     "GC_B","Glandular_epithelium", "Macrophage", "Mast_cell","Memory_B","Naive_B",
                     "NK","pDC","Plasma_cell", "Tfh","Treg")

group_colors <- c("AH" = "#F7D8B7", "ATH" = "#E7C1D7")

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

library(EnhancedVolcano)

Idents(spatial.all) <- "celltype"
DefaultAssay(spatial.all) <- "Spatial"


# 04 | Original differential-expression analyses
logFCfilter=0.50
adjPvalFilter=1e-5

diff_gene.all <- FindAllMarkers(spatial.all, min.pct = 0.5, logfc.threshold = 0.25, group.by = "celltype",
                              only.pos = FALSE, test.use="wilcox")
sig_gene.all <- subset(diff_gene.all, p_val_adj < adjPvalFilter & abs(avg_log2FC) > logFCfilter)
sig_gene.all <- sig_gene.all %>%  group_by(cluster) %>% arrange(desc(avg_log2FC), .by_group = TRUE)

saveRDS(diff_gene.all,"diff/diff_gene.all.ct.RDS")
saveRDS(sig_gene.all,"diff/sig_gene.all.ct.RDS")
write.csv(sig_gene.all,"diff/sig_gene.all.ct.csv", row.names = FALSE)

diff_gene.all.ct <- FindMarkers(spatial.all, min.pct = 0.1, logfc.threshold = 0.25, group.by = "TH_meta",
                              ident.1 = "Y", ident.2 = "N", test.use="wilcox")
diff_gene.all.ct <- diff_gene.all.ct %>% arrange(desc(avg_log2FC), .by_group = TRUE)
dim(diff_gene.all.ct)

saveRDS(diff_gene.all.ct,"diff/diff_gene.all.ct.RDS")
write.csv(diff_gene.all.ct,"diff/diff_gene.all.ct.csv", row.names = TRUE)

library(dplyr)
library(ggplot2)
library(ggrepel)
library(tidyverse)
diff_gene.all.ct <- readRDS("diff/diff_gene.all.ct.RDS")
diff_gene.all.ct$names <- rownames(diff_gene.all.ct)
diff_gene.all.ct <- diff_gene.all.ct %>% mutate(Difference = pct.1-pct.2)
diff_gene.all.ct$group <- ifelse(diff_gene.all.ct$avg_log2FC>0,'UP','Down')
diff_gene.all.ct$group2 <- ifelse((abs(diff_gene.all.ct$avg_log2FC) < log2(1.25) | diff_gene.all.ct$p_val_adj > 1e-5),
                           'ns', diff_gene.all.ct$group)

top.up <- diff_gene.all.ct  %>% top_n(n = 15, wt = avg_log2FC) %>% row.names()
top.dn <- diff_gene.all.ct  %>% top_n(n = -15, wt = avg_log2FC) %>% row.names()
label <- c(top.up, top.dn)
diff_gene.all.ct$label <- ifelse(diff_gene.all.ct$names %in% label, diff_gene.all.ct$names,"")

pdf('diff/degVolcano.all.ct.pdf',height=6, width=6.5)
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

celltypes <- unique(spatial.all$celltype)
celltypes

markers_list <- lapply(celltypes, function(ct) {
  celltype_data <- subset(spatial.all, subset = celltype == ct)
  markers <- FindMarkers(celltype_data, min.pct = 0.2, logfc.threshold = 0.25, group.by = "TH_meta", ident.1 = "Y", ident.2 = "N", test.use="wilcox")
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

write.csv(sig.markers_df, file = "diff/sig.markers_df.csv", row.names = FALSE)
saveRDS(markers_list, file = "diff/markers_list.RDS")
saveRDS(markers_df, file = "diff/markers_df.RDS")

df <- read.csv("diff/sig.markers_df.csv", header = TRUE)
df <- df[-which(df$celltype == "Unknown"),] %>% mutate(chr = as.character(celltype),
    type = if_else(avg_log2FC > 0, "UP_Highly", "Down_Highly"))

pal2 <- c("#7FB2D5","#B3CDE3","#F7D8B7","#17c3b2","#BEDADA","#8DBA8F","#6A3D9A","#00A087","#F39B7FFF","#FFD569",
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
  ) %>% distinct(gene, .keep_all = TRUE)) %>%
  ungroup()


# 05 | Differential-expression visualization
pdf("diff/volcano.pdf", width = 16, height = 8)
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
  scale_fill_manual(values = pal2, guide = "none") +
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

df <- read.csv("diff/sig.markers_df.csv", header = TRUE)
df <- df %>% mutate(chr = as.character(celltype),
    type = if_else(avg_log2FC > 0, "UP_Highly", "Down_Highly"))

ddff <- df %>% filter(type == "Down_Highly") %>%
  group_by(celltype) %>% summarise(n = n_distinct(gene)) %>%
  arrange(desc(n), .by_group = TRUE)

library(ggplot2)

ddff$celltype <- factor(ddff$celltype, levels = rev(ddff$celltype))

pdf("diff/ct_n_deg_dn.pdf", width = 6, height = 4)
ggplot(ddff, aes(x = n, y = celltype)) +
  geom_segment(aes(x = 0, xend = n, yend = celltype),
               color = "grey75", linewidth = 0.8) +
  geom_point(size = 5.5, color = "#4581B2") +
  geom_text(aes(label = n), hjust = 0.5, size = 3, color = "black") +
  labs(x = "Number of down-regulated genes", y = NULL) +
  expand_limits(x = max(ddff$n) * 1.1) +
  theme(axis.text = element_text(size = 10), axis.title   = element_text(size = 11))+
  theme_classic()
dev.off()

library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(UpSetR)
library(grid)

library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(grid)


# 06 | Gene intersections and differential-gene summaries
sig.single <- read.csv("diff/sig.markers_df.csv")

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

pal <- c(
  "#7FB2D5","#B3CDE3","#F7D8B7","#17c3b2","#BEDADA",
  "#8DBA8F","#6A3D9A","#00A087","#F39B7FFF","#FFD569",
  "#984EA3","#BC9DCC","#57C3F3","#F5B47F","#D37295",
  "#F4CAE4","#4581B2"
)

labels_all <- c("CD4_T","CD8_T","cDC1","Ciliated_epithelium","Endothelial","Fibroblast",
                "GC_B","Glandular_epithelium","Macrophage","Mast_cell","Memory_B",
                "Naive_B","NK","pDC","Plasma_cell","Tfh","Treg")

existing_celltypes <- names(sig)
idx <- match(existing_celltypes, labels_all)
pal <- pal[idx]
labels <- labels_all[idx]

m <- make_comb_mat(sig)
cs <- comb_size(m)

pdf("diff/upset.sc2.pdf", width = 12, height = 6.5)
ht <- draw(
  UpSet(m,
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

write.csv(res_df, "diff/upset_intersection_genes_by_plot_order.csv", row.names = FALSE)

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

library(dplyr)
library(ggprism)

gene.ratio <- sig.single %>%
  select(gene, celltype, avg_log2FC) %>%
  distinct() %>%
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

pdf("diff/up_down_gene_count.fc.pdf", width =5, height = 5.5)
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

library(msigdbr)
library(clusterProfiler)
library(enrichplot)

msigdbr_collections()


# 07 | Original pathway enrichment
Hdb <- msigdbr(species = "Homo sapiens", category = "H") %>% dplyr::select(gs_name, gene_symbol)
Kdb <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG_LEGACY") %>% dplyr::select(gs_name, gene_symbol)
Wdb <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:WIKIPATHWAYS") %>% dplyr::select(gs_name, gene_symbol)
BPdb <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP") %>% dplyr::select(gs_name, gene_symbol)

df2 <- df %>% filter(celltype %in% c("Macrophage","GC_B","Tfh","Plasma_cell")) %>% filter(type == "UP_Highly")

pal.fil <- c("#F39B7FFF","#6A3D9A","#F4CAE4","#D37295")

hm1 <- enricher(df2[df2$celltype == "Macrophage",]$gene, TERM2GENE=Hdb)
hm2 <- enricher(df2[df2$celltype == "GC_B",]$gene, TERM2GENE=Hdb)
hm3 <- enricher(df2[df2$celltype == "Tfh",]$gene, TERM2GENE=Hdb)
hm4 <- enricher(df2[df2$celltype == "Plasma_cell",]$gene, TERM2GENE=Hdb)

hm.all <- list(hm1@result, hm2@result, hm3@result, hm4@result)
names(hm.all) <- c("Macrophage","GC_B","Tfh","Plasma_cell")

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

km1 <- enricher(df2[df2$celltype == "Macrophage",]$gene, TERM2GENE=Kdb)
km2 <- enricher(df2[df2$celltype == "GC_B",]$gene, TERM2GENE=Kdb)
km3 <- enricher(df2[df2$celltype == "Tfh",]$gene, TERM2GENE=Kdb)
km4 <- enricher(df2[df2$celltype == "Plasma_cell",]$gene, TERM2GENE=Kdb)

km.all <- list(km1@result, km2@result, km3@result, km4@result)
names(km.all) <- c("Macrophage","GC_B","Tfh","Plasma_cell")

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

bm1 <- enricher(df2[df2$celltype == "Macrophage",]$gene, TERM2GENE=BPdb)
bm2 <- enricher(df2[df2$celltype == "GC_B",]$gene, TERM2GENE=BPdb)
bm3 <- enricher(df2[df2$celltype == "Tfh",]$gene, TERM2GENE=BPdb)
bm4 <- enricher(df2[df2$celltype == "Plasma_cell",]$gene, TERM2GENE=BPdb)

bm.all <- list(bm1@result, bm2@result, bm3@result, bm4@result)
names(bm.all) <- c("Macrophage","GC_B","Tfh","Plasma_cell")

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

saveRDS(hm.all.fli, "diff/hm.all.fli.RDS")
saveRDS(km.all.fli, "diff/km.all.fli.RDS")
saveRDS(bm.all.fli, "diff/bm.all.fli.RDS")

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
    geom_bar(stat = "identity", alpha = 0.5, fill = pal.fil[i], width = 0.8) +
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

pdf("diff/hallmark.all.pdf", width=7.5, height=5)
for (i in seq_along(hm.all.fli1)) {
  p <- plot_pathway(hm.all.fli1[[i]], i) + ggtitle(names(hm.all.fli1)[i])
  print(p)
}
dev.off()

pdf("diff/kegg.all.pdf", width=7.5, height=5)
for (i in seq_along(km.all.fli1)) {
  p <- plot_pathway(km.all.fli1[[i]], i) + ggtitle(names(km.all.fli1)[i])
  print(p)
}
dev.off()

pdf("diff/gobp.all.pdf", width=7.5, height=5)
for (i in seq_along(bm.all.fli1)) {
  p <- plot_pathway(bm.all.fli1[[i]], i) + ggtitle(names(bm.all.fli1)[i])
  print(p)
}
dev.off()

set.seed(1234)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(Augur)
packageVersion("Augur")

spatial.all <- readRDS("/data/YQ/linzh/test/AH/AH_all_spatial.RDS")
spatial.all <- subset(spatial.all, celltype != "Unknown")

spatial.all$celltype <- as.character(spatial.all$celltype)
spatial.all$TH_meta <- as.character(spatial.all$TH_meta)


# 08 | Augur cell-state prioritization
augur <- calculate_auc(spatial.all, cell_type_col = "celltype", label_col = "TH_meta", n_threads = 8)

saveRDS(augur$results, "/data/YQ/linzh/test/AH/atlas/diff/Augur_inte.RDS")

augur_res <- readRDS("/data/YQ/linzh/test/AH/atlas/diff/Augur_inte.RDS")
mat <- augur_res %>% filter(metric == "roc_auc")

mean_scores <- mat %>%
  group_by(cell_type) %>%
  summarise(Mean_Score = mean(estimate, na.rm = TRUE))

sorted_mean_scores <- mean_scores %>% arrange(desc(Mean_Score))

mat$cell_type <- factor(mat$cell_type, levels = sorted_mean_scores$cell_type)

sorted_mean_scores

custom_colors <- c(CD4_T = "#7FB2D5", CD8_T = "#B3CDE3", cDC1 = "#F7D8B7", Ciliated_epithelium = "#17c3b2", Endothelial = "#BEDADA", Fibroblast = "#8DBA8F",
                   GC_B = "#6A3D9A", Glandular_epithelium = "#00A087", Macrophage = "#F39B7FFF", Mast_cell = "#FFD569", Memory_B = "#984EA3", Naive_B = "#BC9DCC",
                   NK = "#57C3F3", pDC = "#F5B47F", Plasma_cell = "#D37295",  Tfh = "#F4CAE4", Treg = "#4581B2")

pdf("diff/Augur_AUC_rank_mean.pdf", useDingbats = F, height = 6, width = 9)
ggplot(mat, aes(x = cell_type, y = estimate)) +
    geom_violin(width = 1, aes(fill = cell_type)) +
    geom_boxplot(width = 0.3, outlier.shape = NA) +
    scale_fill_manual(values = custom_colors) +
    ylab("AUC") + xlab("") + theme_classic() + theme_bw() +
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),aspect.ratio = 0.5)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    legend.position = "top")
dev.off()

median_scores <- mat %>%
  group_by(cell_type) %>%
  summarise(Median_Score = median(estimate, na.rm = TRUE)) %>%
  arrange(desc(Median_Score))

mat$cell_type <- factor(mat$cell_type, levels = median_scores$cell_type)
median_scores

pdf("diff/Augur_AUC_rank_median.pdf", useDingbats = F, height = 6, width = 9)
ggplot(mat, aes(x = cell_type, y = estimate)) +
  geom_violin(width = 1, aes(fill = cell_type), trim = FALSE) +
  geom_boxplot(width = 0.3, outlier.shape = NA) +
  scale_fill_manual(values = custom_colors) +
  ylab("AUC") + xlab("") +
  theme_classic() + theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 0.5)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    legend.position = "top")
dev.off()

library(Seurat)
library(tidyverse)
library(msigdbr)
library(clusterProfiler)
library(AUCell)

Hypoxia <- clusterProfiler::read.gmt("/data/YQ/linzh/pub_data/ref/pathway/HALLMARK_HYPOXIA.v2026.1.Hs.gmt")
Inflammatory_response <- clusterProfiler::read.gmt("/data/YQ/linzh/pub_data/ref/pathway/HALLMARK_INFLAMMATORY_RESPONSE.v2026.1.Hs.gmt")
GeneSet <- list(Hypoxia = Hypoxia$gene, Inflammatory_response = Inflammatory_response$gene)

library(AUCell)


# 09 | AUCell pathway scoring and plots
run_AUCell <- function(seu, geneSets, assay = "Spatial", slot = "data",
                       aucMaxRank = 0.1, nCores = 1) {
  expr <- GetAssayData(seu, assay = assay, slot = slot)

  cells_ranking <- AUCell_buildRankings(expr, nCores = nCores, plotStats = FALSE)

  if (aucMaxRank < 1) {
    aucMaxRank <- nrow(cells_ranking) * aucMaxRank
  }

  cells_AUC <- AUCell_calcAUC(geneSets, cells_ranking,
                              aucMaxRank = aucMaxRank, nCores = nCores)
  auc <- t(getAUC(cells_AUC))
  auc <- as.data.frame(auc)
  return(auc)
}

auc_df1 <- run_AUCell(spatial.all, GeneSet, assay = "Spatial", slot = "data", aucMaxRank = 0.1)
colnames(auc_df1) <- paste0("AUC_", colnames(auc_df1))
spatial.all <- AddMetaData(spatial.all, metadata = auc_df1)

str(spatial.all@meta.data)

pal.cl <-c ("#9EC4BE","#FAC795")
ES = as.data.frame(spatial.all@meta.data[, str_detect(colnames(spatial.all@meta.data), "^AUC")])

pdf("diff/AUCell_AddModuleScore_violin.pdf", useDingbats = F, width = 4, height = 4)
for (i in 1:ncol(ES)) {
  df = data.frame(Activity = as.numeric(ES[,i]), cluster = spatial.all$TH_meta)
  p <- ggplot(df, aes(x=cluster, y=Activity))+
    geom_violin(width=0.8, aes(fill = cluster)) + ylab("Activity") + xlab("Groups") +
    geom_boxplot(width = 0.3, outlier.shape = NA) + ggtitle(colnames(ES)[i]) +
    theme_classic() + theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank()) +
    theme(aspect.ratio = 0.5)+scale_fill_manual(values = pal.cl)+
    theme(axis.text.x = element_text(colour="black",size = 10,angle = 30,hjust = 0.6,vjust = 0.75),
          axis.text.y = element_text(colour="black",size = 10),axis.title = element_text(size=10),
          legend.title = element_blank(), legend.key.size = unit(10, "pt"))+
    theme(plot.title = element_text(hjust = 0.5,size = 12, face = "bold"))
  print(p)
}
dev.off()

pdf("diff/AUCell_AddModuleScore_box.pdf", useDingbats = F,width = 3.5, height = 4)
for (i in 1:ncol(ES)) {
  df = data.frame(Activity = as.numeric(ES[,i]), cluster = spatial.all$TH_meta)
  p<-ggboxplot(df, x="cluster", y="Activity", width = 0.6,pt.size=0.5,
               color = "black",
               fill="cluster",
               palette = pal.cl,
               ylab="Activity",xlab = "" ,
               bxp.errorbar=T,
               bxp.errorbar.width=0.5,
               size=1,
               outlier.shape=NA,
               legend = "none")+ggtitle(colnames(ES)[i]) +
    stat_compare_means(comparisons = list(c("N","Y")), method = "wilcox.test")+
    theme(axis.text.x = element_text(colour="black",size = 10,angle = 30,hjust = 0.6,vjust = 0.75),
          axis.text.y = element_text(colour="black",size = 10),axis.title = element_text(size=10),
          legend.title = element_blank(), legend.key.size = unit(10, "pt"))+
    theme(plot.title = element_text(hjust = 0.5,size = 12, face = "bold"))
  print(p)
}
dev.off()

pdf("diff/AUCell_AddModuleScore_viplin2.pdf", useDingbats = F,width = 4, height = 4)
for (i in 1:ncol(ES)) {
  df = data.frame(Activity = as.numeric(ES[,i]), cluster = spatial.all$TH_meta)
  p<-ggviolin(df, x="cluster", y="Activity",fill = "cluster", width=0.8,
              ylab="Activity",xlab= "" ,legend = "none")+
    stat_compare_means(comparisons = list(c("N","Y")),
                       method = "wilcox.test")+
    geom_boxplot(width = 0.3, outlier.shape = NA) + ggtitle(colnames(ES)[i])+
    scale_fill_manual(values = pal.cl)
    theme_classic() + theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank()) +
    theme(axis.text.x = element_text(colour="black",size = 10,angle = 30, hjust = 0.6,vjust = 0.75),
          axis.text.y = element_text(colour="black",size = 10),axis.title = element_text(size=10),
          legend.title = element_blank(), legend.key.size = unit(10, "pt"))+
    theme(plot.title = element_text(hjust = 0.5,size = 12, face = "bold"))

  print(p)
}
dev.off()





# ==============================================================================
# 10 | Additional donor-level analyses
# Comparison: AH (TH_meta = N) vs ATH-A (TH_meta = Y).
# Biological replicates are patients, not individual cells.
# Two-sided unpaired Welch t-tests; Benjamini-Hochberg FDR shown alongside p.
# Figures retain the group colors from the original sample-level plots.

group_colors_sample <- c("AH" = "#F7D8B7", "ATH-A" = "#E7C1D7")
group_levels_sample <- c("AH", "ATH-A")
celltype_levels_sample <- c(
  "CD4_T", "CD8_T", "cDC1", "Ciliated_epithelium", "Endothelial",
  "Fibroblast", "GC_B", "Glandular_epithelium", "Macrophage",
  "Mast_cell", "Memory_B", "Naive_B", "NK", "pDC",
  "Plasma_cell", "Tfh", "Treg"
)

required_columns <- c("Patient", "TH_meta", "celltype")
if (!all(required_columns %in% colnames(spatial.review@meta.data))) {
  stop("The Seurat object must contain Patient, TH_meta, and celltype metadata.")
}

spatial_meta <- data.frame(
  barcode = rownames(spatial.review@meta.data),
  Patient = as.character(spatial.review$Patient),
  TH_meta = as.character(spatial.review$TH_meta),
  celltype = as.character(spatial.review$celltype),
  stringsAsFactors = FALSE
) %>%
  filter(
    !is.na(Patient), Patient != "", Patient != "NA",
    TH_meta %in% c("N", "Y"),
    !is.na(celltype), celltype != ""
  ) %>%
  mutate(Group = if_else(TH_meta == "N", "AH", "ATH-A"))

patient_group <- spatial_meta %>%
  distinct(Patient, Group)

if (anyDuplicated(patient_group$Patient)) {
  stop("Each patient must belong to exactly one AH/ATH-A group.")
}
if (!all(group_levels_sample %in% patient_group$Group)) {
  stop("Both AH and ATH-A patients must be present.")
}

print(patient_group %>% arrange(Group, Patient))
print(table(patient_group$Group))
if (any(table(patient_group$Group) < 2L)) {
  stop("At least two patients per group are required for a t-test.")
}

# Plot only the 17 annotated cell types, but use all QC-retained cells
# (including Unknown) in each patient's proportion denominator.
patient_totals <- spatial_meta %>%
  count(Patient, Group, name = "total_cells")

celltype_counts <- spatial_meta %>%
  filter(celltype %in% celltype_levels_sample) %>%
  count(Patient, Group, celltype, name = "n_cells")

proportion_donor <- tidyr::crossing(
  patient_group, celltype = celltype_levels_sample
) %>%
  left_join(celltype_counts, by = c("Patient", "Group", "celltype")) %>%
  left_join(patient_totals, by = c("Patient", "Group")) %>%
  mutate(
    n_cells = replace_na(n_cells, 0L),
    Proportion = 100 * n_cells / total_cells,
    Group = factor(Group, levels = group_levels_sample),
    celltype = factor(celltype, levels = celltype_levels_sample)
  )

write.csv(proportion_donor, "diff/sample_level_celltype_proportions.csv", row.names = FALSE)

# Mean log-normalized FKBP5 expression per patient and cell type.
# Zero-expression cells are included. No new input data are loaded.
if (!"FKBP5" %in% rownames(spatial.review[["Spatial"]])) {
  stop("FKBP5 is not present in the Spatial assay.")
}

fkbp5_matrix <- GetAssayData(spatial.review, assay = "Spatial", slot = "data")
spatial_meta$FKBP5 <- as.numeric(
  fkbp5_matrix["FKBP5", spatial_meta$barcode, drop = TRUE]
)
rm(fkbp5_matrix)

fkbp5_by_donor <- spatial_meta %>%
  filter(celltype %in% celltype_levels_sample) %>%
  group_by(Patient, Group, celltype) %>%
  summarise(
    mean_FKBP5 = mean(FKBP5, na.rm = TRUE),
    n_cells = n(),
    .groups = "drop"
  )

fkbp5_by_donor <- tidyr::crossing(
  patient_group, celltype = celltype_levels_sample
) %>%
  left_join(
    fkbp5_by_donor,
    by = c("Patient", "Group", "celltype")
  ) %>%
  mutate(
    Group = factor(Group, levels = group_levels_sample),
    celltype = factor(celltype, levels = celltype_levels_sample)
  )

write.csv(fkbp5_by_donor, "diff/sample_level_FKBP5_by_celltype.csv", row.names = FALSE)

# Compute one unpaired, two-sided Welch t-test per cell type.
# NA means that a comparison has insufficient/constant sample-level data.
format_p <- function(p) {
  if (length(p) != 1L || !is.finite(p)) return("NA")
  if (p < 1e-4) return(formatC(p, format = "e", digits = 2))
  sprintf("%.4f", p)
}

compare_by_donor <- function(dat, value_col) {
  result <- lapply(celltype_levels_sample, function(ct) {
    d <- dat[as.character(dat$celltype) == ct, , drop = FALSE]
    ah <- d[[value_col]][d$Group == "AH"]
    ath <- d[[value_col]][d$Group == "ATH-A"]
    ah <- ah[is.finite(ah)]
    ath <- ath[is.finite(ath)]

    p_val <- if (length(ah) >= 2L && length(ath) >= 2L) {
      tryCatch(
        stats::t.test(ah, ath, alternative = "two.sided",
                      paired = FALSE, var.equal = FALSE)$p.value,
        error = function(e) NA_real_
      )
    } else {
      NA_real_
    }

    data.frame(
      celltype = ct,
      n_AH = length(ah),
      n_ATH_A = length(ath),
      mean_AH = if (length(ah)) mean(ah) else NA_real_,
      mean_ATH_A = if (length(ath)) mean(ath) else NA_real_,
      p_value = p_val
    )
  })

  result <- bind_rows(result)
  result$BH_FDR <- p.adjust(result$p_value, method = "BH")
  result$celltype <- factor(
    result$celltype, levels = celltype_levels_sample
  )
  result
}

proportion_tests <- compare_by_donor(proportion_donor, "Proportion")
fkbp5_tests <- compare_by_donor(fkbp5_by_donor, "mean_FKBP5")

write.csv(proportion_tests, "diff/sample_level_celltype_proportion_ttests.csv", row.names = FALSE)
write.csv(fkbp5_tests, "diff/sample_level_FKBP5_ttests.csv", row.names = FALSE)

summary_by_donor <- function(dat, value_col) {
  dat %>%
    group_by(celltype, Group) %>%
    summarise(
      mean_value = if (all(is.na(.data[[value_col]]))) NA_real_
                   else mean(.data[[value_col]], na.rm = TRUE),
      n = sum(is.finite(.data[[value_col]])),
      se = if (sum(is.finite(.data[[value_col]])) > 1L) {
        sd(.data[[value_col]], na.rm = TRUE) /
          sqrt(sum(is.finite(.data[[value_col]])))
      } else {
        NA_real_
      },
      .groups = "drop"
    )
}


draw_overview <- function(dat, tests, value_col, y_label) {
  means <- summary_by_donor(dat, value_col)
  labels <- tests %>%
    mutate(p_label = paste0("p = ", vapply(p_value, format_p, character(1))))

  ggplot() +
    geom_col(data = means,
      aes(x = Group, y = mean_value, fill = Group),
      width = 0.6, alpha = 0.75, color = "#354B52",
      na.rm = TRUE) +
    geom_errorbar(data = means,
      aes(x = Group, ymin = pmax(0, mean_value - se),
        ymax = mean_value + se),
      width = 0.15, linewidth = 0.5, na.rm = TRUE) +
    geom_point(data = dat,
      aes(x = Group, y = .data[[value_col]], fill = Group),
      shape = 21, size = 2, color = "#354B52", stroke = 0.35,
      position = position_jitter(width = 0.08, height = 0, seed = 1234),
      na.rm = TRUE) +
    geom_text(data = labels,
      aes(x = 1.5, y = Inf, label = p_label),
      inherit.aes = FALSE, vjust = 1.1, size = 2.5) +
    facet_wrap(~celltype, ncol = 4, scales = "free_y") +
    scale_fill_manual(values = group_colors_sample, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.04, 0.24))) +
    labs(
      x = NULL, y = y_label,
      title = "AH vs ATH-A: donor-level comparisons",
      subtitle = "Bars: mean +/- SEM; points: individual patients; unpaired Welch t-test"
    ) +
    theme_classic(base_size = 11) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(face = "bold", size = 9),
      axis.text.x = element_text(size = 9),
      panel.spacing = unit(0.9, "lines"),
      plot.title = element_text(face = "bold")
    )
}

draw_one_celltype <- function(dat, tests, cell_type, value_col, y_label) {
  sub_dat <- dat[as.character(dat$celltype) == cell_type, , drop = FALSE]
  means <- summary_by_donor(sub_dat, value_col)
  test <- tests[as.character(tests$celltype) == cell_type, , drop = FALSE]

  vals <- c(sub_dat[[value_col]], means$mean_value + means$se)
  ymax <- suppressWarnings(max(vals[is.finite(vals)], 0))
  gap <- if (ymax > 0) 0.12 * ymax else 0.05
  y_bracket <- ymax + gap
  label_position <- y_bracket + 0.55 * gap
  y_limit <- label_position + 0.8 * gap

  p_label <- paste0("p = ", format_p(test$p_value), ", BH q = ", format_p(test$BH_FDR))

  ggplot() + geom_col(data = means,
      aes(x = Group, y = mean_value, fill = Group),
      width = 0.6, alpha = 0.75, color = "#354B52",
      na.rm = TRUE) +
    geom_errorbar(data = means,
      aes(x = Group, ymin = pmax(0, mean_value - se),
        ymax = mean_value + se),
      width = 0.15, linewidth = 0.5, na.rm = TRUE) +
    geom_point(data = sub_dat,
      aes(x = Group, y = .data[[value_col]], fill = Group),
      shape = 21, color = "#354B52", stroke = 0.35, size = 2,
      position = position_jitter(width = 0.08, height = 0, seed = 1234),
      na.rm = TRUE) +
    annotate("segment", x = 1, xend = 2, y = y_bracket, yend = y_bracket, linewidth = 0.55) +
    annotate("segment", x = c(1, 2), xend = c(1, 2),
      y = y_bracket - 0.1 * gap, yend = y_bracket,
      linewidth = 0.55) +
    annotate("text", x = 1.5, y = label_position, label = p_label, size = 4) +
    scale_fill_manual(values = group_colors_sample, drop = FALSE) +
    scale_y_continuous(limits = c(0, y_limit), expand = c(0, 0)) +
    labs(
      title = cell_type, x = NULL, y = y_label,
      subtitle = paste0(
        "AH: n = ", test$n_AH,
        " patients   |   ATH-A: n = ", test$n_ATH_A, " patients"
      )
    ) +
    theme_classic(base_size = 13) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
      plot.subtitle = element_text(hjust = 0.5, size = 12)
    )
}

# Each PDF has one overview page followed by one page per cell type.
pdf("diff/AH_vs_ATHA_celltype_proportion_samplelevel.pdf", width = 4, height = 4, onefile = TRUE, useDingbats = FALSE)
for (ct in celltype_levels_sample) {
  print(draw_one_celltype(
    proportion_donor, proportion_tests, ct,
    "Proportion", "Cell proportion (%)"
  ))
}
dev.off()

pdf("diff/AH_vs_ATHA_FKBP5_by_celltype_samplelevel.pdf", width = 4, height = 4, onefile = TRUE, useDingbats = FALSE)
for (ct in celltype_levels_sample) {
  print(draw_one_celltype(
    fkbp5_by_donor, fkbp5_tests, ct,
    "mean_FKBP5", "Mean log-normalized FKBP5 expression"
  ))
}
dev.off()



# The original exploratory cell-level DE and pathway analyses follow below.
pdf("diff/AH_vs_ATHA_FKBP5_proportion_samplelevel.pdf", width = 16, height = 18, onefile = TRUE, useDingbats = FALSE)
print(draw_overview(
  proportion_donor, proportion_tests, "Proportion", "Cell proportion (%)"
))

print(draw_overview(
  fkbp5_by_donor, fkbp5_tests,
  "mean_FKBP5", "Mean log-normalized FKBP5 expression"
))
dev.off()


