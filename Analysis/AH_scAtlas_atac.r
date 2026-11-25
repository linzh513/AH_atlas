# scATAC-seq preprocessing, common-peak integration, and RNA-to-ATAC label transfer
.libPaths("/home/lzh/miniconda3/envs/singlecell/lib/R/library")
set.seed(1234)
setwd("/mnt/data_nvme1/lzh/test/AH/sc_atlas_atac/pre")
rm(list = ls())
suppressPackageStartupMessages({
  library(Signac)
  library(Seurat)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(ensembldb)
  library(AnnotationHub)
  library(biovizBase)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tibble)
  library(future)
})
packageVersion("Seurat")
packageVersion("SeuratObject")
packageVersion("Signac")
theme_global <- theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(aspect.ratio = 1)

# 01. Input and output paths
base_dir <- "/mnt/sda1/lzh/data/inhouse/AH/sc/scATAC"
out_dir <- "/mnt/data_nvme1/lzh/test/AH/sc_atlas_atac/pre/batch_pre_commonpeaks"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "objects"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "plots"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "tables"), showWarnings = FALSE, recursive = TRUE)
plan("multisession", workers = 6)
options(future.globals.maxSize = 20 * 1024^3)
cat("Base directory:", base_dir, "\n")
cat("Output directory:", out_dir, "\n")

# 02. Genome annotation and blacklist
cat("\n--- Loading genome annotations ---\n")
ah <- AnnotationHub(localHub = FALSE)
query(ah, "EnsDb.Hsapiens.v98")
ensdb_v98 <- ah[["AH75011"]]
annotations <- GetGRangesFromEnsDb(ensdb = ensdb_v98)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "hg38"
blacklist_regions <- ah[["AH107305"]]  # hg38 blacklist
seqlevelsStyle(blacklist_regions) <- "UCSC"
blacklist_regions <- keepStandardChromosomes(blacklist_regions, pruning.mode = "coarse")
genome(blacklist_regions) <- "hg38"
cat("Annotation loaded.\n")
cat("Annotation seqlevels example:\n")
print(head(seqlevels(annotations)))
cat("Blacklist regions:", length(blacklist_regions), "\n")

# 03. Discover samples
sample_dirs <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)
sample_dirs <- sample_dirs[file.exists(sample_dirs)]
sample_ids <- basename(sample_dirs)
cat("Found samples:\n")
print(sample_ids)

# 04. Construct a common peak set
cat("\n--- Building common peak set across all samples ---\n")
safe_read_matrix <- function(mtx_dir) {
  counts <- Read10X(data.dir = mtx_dir)
  return(counts)
}

peak_names_to_granges <- function(peak_names) {
  peak_names <- peak_names[!is.na(peak_names)]
  gr <- tryCatch(
    {
      StringToGRanges(peak_names, sep = c(":", "-"))
    },
    error = function(e) {
      message("StringToGRanges with sep c(':', '-') failed. Trying c('-', '-') ...")
      StringToGRanges(peak_names, sep = c("-", "-"))
    }
  )
  return(gr)
}

get_sample_peaks <- function(sample_dir, sample_id) {
  mtx_dir <- file.path(sample_dir, "filtered_peaks_bc_matrix")
  if (!dir.exists(mtx_dir)) {
    stop("Missing filtered_peaks_bc_matrix in ", sample_dir)
  }
  counts <- safe_read_matrix(mtx_dir)
  peak_names <- rownames(counts)
  peaks_gr <- peak_names_to_granges(peak_names)
  seqlevelsStyle(peaks_gr) <- "UCSC"
  peaks_gr <- keepStandardChromosomes(peaks_gr, pruning.mode = "coarse")
  genome(peaks_gr) <- "hg38"
  cat(sample_id, "original peaks:", length(peaks_gr), "\n")
  return(peaks_gr)
}

peak_list <- lapply(seq_along(sample_dirs), function(i) {
  get_sample_peaks(sample_dirs[i], sample_ids[i])
})
names(peak_list) <- sample_ids
print(sapply(peak_list, class))
if (!all(sapply(peak_list, function(x) inherits(x, "GRanges")))) {
  stop("Some elements in peak_list are not GRanges.")
}

peak_grl <- GenomicRanges::GRangesList(peak_list)
all_peaks <- unlist(peak_grl, use.names = FALSE)
cat("Total peaks before reduce:", length(all_peaks), "\n")
combined_peaks <- GenomicRanges::reduce(all_peaks)
combined_peaks <- keepStandardChromosomes(combined_peaks, pruning.mode = "coarse")
seqlevelsStyle(combined_peaks) <- "UCSC"
genome(combined_peaks) <- "hg38"
peak_widths <- width(combined_peaks)
combined_peaks <- combined_peaks[peak_widths > 20 & peak_widths < 10000]
combined_peaks <- sortSeqlevels(combined_peaks)
combined_peaks <- sort(combined_peaks)
cat("Common peaks after reduce/filter:", length(combined_peaks), "\n")
saveRDS(combined_peaks, file = file.path(out_dir, "objects", "combined_peaks.RDS"))

# 05. Sample processing functions
read_one_sample <- function(sample_dir, sample_id, combined_peaks) {
  mtx_dir <- file.path(sample_dir, "filtered_peaks_bc_matrix")
  bed_dir <- file.path(sample_dir, "bed")
  frag_dir <- file.path(sample_dir, "fragments")
  cat("\nReading sample with common peaks:", sample_id, "\n")
  if (!dir.exists(mtx_dir)) {
    stop("Missing filtered_peaks_bc_matrix in ", sample_dir)
  }
  # Use the original peak matrix to identify filtered cell barcodes.
  raw_counts <- safe_read_matrix(mtx_dir)
  cells_use <- colnames(raw_counts)
  qc_file <- file.path(bed_dir, "per_barcode_metrics.csv")
  if (!file.exists(qc_file)) {
    stop("Missing per_barcode_metrics.csv in ", sample_dir)
  }
  qc_meta <- read.csv(qc_file, stringsAsFactors = FALSE)
  if (!"barcode" %in% colnames(qc_meta)) {
    colnames(qc_meta)[1] <- "barcode"
  }
  rownames(qc_meta) <- qc_meta$barcode
  # Keep barcodes shared by the count matrix and QC metadata.
  cells_use <- intersect(cells_use, rownames(qc_meta))
  qc_meta <- qc_meta[cells_use, , drop = FALSE]
  cat(sample_id, "cells used for FeatureMatrix:", length(cells_use), "\n")
  frag_file <- list.files(
    frag_dir,
    pattern = "fragments\\.tsv\\.gz$",
    full.names = TRUE
  )
  if (length(frag_file) == 0) {
    stop("Missing fragments.tsv.gz in ", sample_dir)
  }
  frag_file <- frag_file[1]
  if (!file.exists(paste0(frag_file, ".tbi"))) {
    stop("Missing fragments index .tbi for ", frag_file)
  }
  frag_obj <- CreateFragmentObject(
    path = frag_file,
    cells = cells_use
  )
  # Recount fragment overlaps against the common peak set.
  counts <- FeatureMatrix(
    fragments = frag_obj,
    features = combined_peaks,
    cells = cells_use
  )
  chrom_assay <- CreateChromatinAssay(
    counts = counts,
    fragments = frag_obj,
    annotation = annotations,
    genome = "hg38",
    min.cells = 10,
    min.features = 200
  )
  obj <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "ATAC",
    project = sample_id,
    meta.data = qc_meta
  )
  obj$sample_id <- sample_id
  obj$sample_dir <- sample_dir
  Annotation(obj) <- annotations
  return(obj)
}

add_qc_metrics <- function(obj) {
  DefaultAssay(obj) <- "ATAC"
  obj <- NucleosomeSignal(obj)
  obj <- TSSEnrichment(obj)
  obj$blacklist_ratio <- FractionCountsInRegion(
    object = obj,
    assay = "ATAC",
    regions = blacklist_regions
  )
  # Recalculate depth metrics from the common-peak matrix.
  if ("fragments_num" %in% colnames(obj@meta.data)) {
    obj$fragments_num <- as.numeric(obj$fragments_num)
    obj$pct_reads_in_peaks <- ifelse(
      obj$fragments_num > 0,
      obj$nCount_ATAC / obj$fragments_num * 100,
      NA
    )
  } else if ("fraction_frag_overlap_peak" %in% colnames(obj@meta.data)) {
    obj$fraction_frag_overlap_peak <- as.numeric(obj$fraction_frag_overlap_peak)
    obj$pct_reads_in_peaks <- obj$fraction_frag_overlap_peak * 100
  } else if (all(c("fragments_num_overlap_peaks", "fragments_num") %in% colnames(obj@meta.data))) {
    obj$fragments_num <- as.numeric(obj$fragments_num)
    obj$fragments_num_overlap_peaks <- as.numeric(obj$fragments_num_overlap_peaks)
    obj$pct_reads_in_peaks <- ifelse(
      obj$fragments_num > 0,
      obj$fragments_num_overlap_peaks / obj$fragments_num * 100,
      NA
    )
  } else if (all(c("events_overlap_peak", "atac_reads") %in% colnames(obj@meta.data))) {
    obj$events_overlap_peak <- as.numeric(obj$events_overlap_peak)
    obj$atac_reads <- as.numeric(obj$atac_reads)
    obj$pct_reads_in_peaks <- ifelse(
      obj$atac_reads > 0,
      obj$events_overlap_peak / obj$atac_reads * 100,
      NA
    )
  } else {
    obj$pct_reads_in_peaks <- NA
  }
  obj$nucleosome_group <- ifelse(obj$nucleosome_signal > 2, "NS > 2", "NS < 2")
  return(obj)
}

qc_plots <- function(obj, out_prefix) {
  pdf(paste0(out_prefix, "_qc_density.pdf"), width = 7, height = 6)
  print(
    DensityScatter(
      obj,
      x = "nCount_ATAC",
      y = "TSS.enrichment",
      log_x = TRUE,
      quantiles = TRUE
    )
  )
  dev.off()
  pdf(paste0(out_prefix, "_fragment_hist.pdf"), width = 8, height = 4)
  print(FragmentHistogram(object = obj, group.by = "nucleosome_group"))
  dev.off()
  plot_features <- c("nCount_ATAC", "TSS.enrichment", "nucleosome_signal", "pct_reads_in_peaks", "blacklist_ratio")
  plot_features <- plot_features[
    sapply(plot_features, function(f) {
      x <- obj[[f]][, 1]
      !all(is.na(x)) && length(unique(x[!is.na(x)])) > 1
    })
  ]
  pdf(paste0(out_prefix, "_quality_violin.pdf"), width = 15, height = 5)
  print(
    VlnPlot(
      obj,
      features = plot_features,
      pt.size = 0.1,
      ncol = length(plot_features)
    )
  )
  dev.off()
}

filter_sample <- function(obj) {
  obj_f <- subset(obj,
    subset =
      nCount_ATAC > 1000 &
      nCount_ATAC < 30000 &
      TSS.enrichment > 3 &
      nucleosome_signal < 2 &
      pct_reads_in_peaks > 25 &
      blacklist_ratio < 0.01
  )
  return(obj_f)
}

run_lsi <- function(obj) {
  DefaultAssay(obj) <- "ATAC"
  obj <- RunTFIDF(obj)
  obj <- FindTopFeatures(obj, min.cutoff = "q5")
  obj <- RunSVD(
    obj,
    assay = "ATAC",
    reduction.name = "lsi",
    reduction.key = "LSI_"
  )
  obj <- RunUMAP(obj, reduction = "lsi", dims = 2:30)
  obj <- FindNeighbors(obj, reduction = "lsi", dims = 2:30)
  obj <- FindClusters(obj, verbose = FALSE, algorithm = 3)
  return(obj)
}

# 06. Process individual samples
sample_objects <- list()
summary_list <- list()
for (sd in sample_dirs) {
  sid <- basename(sd)
  cat("\n==============================\n")
  cat("Processing sample:", sid, "\n")
  cat("==============================\n")
  obj <- read_one_sample(
    sample_dir = sd,
    sample_id = sid,
    combined_peaks = combined_peaks
  )
  saveRDS(obj, file = file.path(out_dir, "objects", paste0(sid, "_commonpeak_raw.RDS")))
  obj <- add_qc_metrics(obj)
  qc_plots(obj, file.path(out_dir, "plots", sid))
  obj_f <- filter_sample(obj)
  cat("Cells before filtering:", ncol(obj), "\n")
  cat("Cells after filtering:", ncol(obj_f), "\n")
  saveRDS(obj_f, file = file.path(out_dir, "objects", paste0(sid, "_commonpeak_filtered.RDS")))
  obj_f <- run_lsi(obj_f)
  pdf(file.path(out_dir, "plots", paste0(sid, "_commonpeak_atac_umap.pdf")), width = 7, height = 6)
  print(
    DimPlot(obj_f, label = TRUE, repel = TRUE) +
      ggtitle(paste0(sid, " ATAC UMAP common peaks")) +
      theme_global
  )
  dev.off()
  saveRDS(obj_f, file = file.path(out_dir, "objects", paste0(sid, "_commonpeak_lsi.RDS"))
  summary_tbl <- data.frame(
    sample_id = sid,
    cells_before = ncol(obj),
    cells_after = ncol(obj_f),
    peaks = nrow(obj_f),
    stringsAsFactors = FALSE
  )
  write.csv(
    summary_tbl,
    file = file.path(out_dir, "tables", paste0(sid, "_summary.csv")),
    row.names = FALSE
  )
  saveRDS(
    summary_tbl,
    file = file.path(out_dir, "tables", paste0(sid, "_summary.RDS"))
  )
  summary_list[[sid]] <- summary_tbl
  sample_objects[[sid]] <- obj_f
  cat("Finished sample:", sid, "\n")
}

# 07. Summarize sample-level QC
summary_all <- bind_rows(summary_list)
write.csv(summary_all, file = file.path(out_dir, "tables", "all_samples_summary.csv"), row.names = FALSE)
saveRDS(summary_all, file = file.path(out_dir, "tables", "all_samples_summary.RDS"))
cat("\nAll samples processed.\n")
print(summary_all)
saveRDS(sample_objects, file = file.path(out_dir, "objects", "all_samples_commonpeak_lsi_list.RDS"))
cat("\nBatch preprocessing with common peaks finished.\n")
cat("Outputs saved to:", out_dir, "\n")

# 08. Merge common-peak ATAC objects
cat("\n--- Merging all common-peak ATAC objects ---\n")
sample_objects <- sample_objects[names(sample_objects)]
atac_merged <- merge(
  x = sample_objects[[1]],
  y = sample_objects[2:length(sample_objects)],
  add.cell.ids = names(sample_objects)
)
DefaultAssay(atac_merged) <- "ATAC"
cat("Merged object cells:", ncol(atac_merged), "\n")
cat("Merged object peaks:", nrow(atac_merged), "\n")
saveRDS(atac_merged, file = file.path(out_dir, "objects", "all_samples_commonpeak_merged.RDS"))

# 09. Merged-object LSI, UMAP, and clustering
cat("\n--- Running merged LSI ---\n")
DefaultAssay(atac_merged) <- "ATAC"
atac_merged <- RunTFIDF(atac_merged)
atac_merged <- FindTopFeatures(atac_merged, min.cutoff = "q5")
atac_merged <- RunSVD(
  atac_merged,
  assay = "ATAC",
  reduction.name = "lsi",
  reduction.key = "LSI_"
)

pdf(file.path(out_dir, "plots", "all_samples_depthcor.pdf"), width = 7, height = 5)
print(DepthCor(atac_merged))
dev.off()
atac_merged <- RunUMAP(
  atac_merged,
  reduction = "lsi",
  dims = 2:30,
  reduction.name = "umap.raw",
  reduction.key = "UMAPRAW_"
)
atac_merged <- FindNeighbors(
  atac_merged,
  reduction = "lsi",
  dims = 2:30
)
atac_merged <- FindClusters(
  atac_merged,
  algorithm = 3,
  resolution = 0.5,
  verbose = FALSE
)

pdf(file.path(out_dir, "plots", "all_samples_raw_umap_by_sample.pdf"), width = 7, height = 6)
print(
  DimPlot(
    atac_merged,
    reduction = "umap.raw",
    group.by = "sample_id"
  ) +
    ggtitle("Merged ATAC raw LSI UMAP by sample") +
    theme_global
)
dev.off()

pdf(file.path(out_dir, "plots", "all_samples_raw_umap_by_cluster.pdf"), width = 7, height = 6)
print(
  DimPlot(
    atac_merged,
    reduction = "umap.raw",
    label = TRUE,
    repel = TRUE
  ) +
    ggtitle("Merged ATAC raw LSI UMAP by cluster") +
    theme_global
)
dev.off()

saveRDS(atac_merged, file = file.path(out_dir, "objects", "all_samples_commonpeak_merged_lsi.RDS"))

# 10. Prepare for rLSI integration
cat("\n--- Preparing objects for rLSI integration ---\n")

# Sequential processing avoids copying large objects across workers.
plan("sequential")
options(future.globals.maxSize = 100 * 1024^3)
if (!exists("sample_objects")) {
  sample_objects <- readRDS(
    file.path(out_dir, "objects", "all_samples_commonpeak_lsi_list.RDS")
  )
}

if (!exists("atac_merged")) {
  atac_merged <- readRDS(
    file.path(out_dir, "objects", "all_samples_commonpeak_merged_lsi.RDS")
  )
}

DefaultAssay(atac_merged) <- "ATAC"
cat("sample_objects:\n")
print(names(sample_objects))
print(sapply(sample_objects, ncol))
cat("atac_merged cells:", ncol(atac_merged), "\n")
cat("atac_merged peaks:", nrow(atac_merged), "\n")
cat("atac_merged reductions:\n")
print(Reductions(atac_merged))

# 11. Match cell barcodes across ATAC objects
cat("\n--- Renaming cells in sample_objects ---\n")
sample_ids <- names(sample_objects)
sample_objects <- lapply(sample_ids, function(sid) {
  obj <- sample_objects[[sid]]
  old_names <- colnames(obj)
  # Remove existing sample prefixes before applying consistent cell names.
  raw_names <- sub(paste0("^", sid, "_"), "", old_names)
  new_names <- paste0(sid, "_", raw_names)
  if (!identical(old_names, new_names)) {
    obj <- RenameCells(
      object = obj,
      new.names = new_names
    )
    cat(sid, ": renamed cells.\n")
  } else {
    cat(sid, ": cell names already correct.\n")
  }
  obj$sample_id <- sid
  DefaultAssay(obj) <- "ATAC"
  return(obj)
})
names(sample_objects) <- sample_ids
dup_n <- sum(duplicated(unlist(lapply(sample_objects, colnames))))
cat("Duplicated cell names after RenameCells:", dup_n, "\n")
if (dup_n > 0) {
  stop("Still found duplicated cell names after RenameCells.")
}

cat("Example cell names:\n")
print(head(colnames(sample_objects[[1]])))
saveRDS(
  sample_objects,
  file = file.path(out_dir, "objects", "all_samples_commonpeak_lsi_list_renamed.RDS")
)

# 12. Verify cell barcode consistency
cat("\n--- Checking cell name consistency ---\n")
cells_from_list <- unlist(lapply(sample_objects, colnames))
cells_from_merged <- colnames(atac_merged)
cat("Cells in sample_objects:", length(cells_from_list), "\n")
cat("Cells in atac_merged:", length(cells_from_merged), "\n")
same_cells <- setequal(cells_from_list, cells_from_merged)
cat("Same cell set between sample_objects and atac_merged:", same_cells, "\n")
if (!same_cells) {
  missing_in_merged <- setdiff(cells_from_list, cells_from_merged)
  missing_in_list <- setdiff(cells_from_merged, cells_from_list)
  cat("Cells in sample_objects but not in atac_merged:", length(missing_in_merged), "\n")
  cat("Cells in atac_merged but not in sample_objects:", length(missing_in_list), "\n")
  if (length(missing_in_merged) > 0) {
    print(head(missing_in_merged))
  }
  if (length(missing_in_list) > 0) {
    print(head(missing_in_list))
  }
  stop("Cell names do not match. Do not continue integration.")
}

cat("Cell names matched. Continue to rLSI integration.\n")

# 13. Find rLSI integration anchors
cat("\n--- Finding integration anchors with rLSI ---\n")
DefaultAssay(atac_merged) <- "ATAC"
if (!"lsi" %in% Reductions(atac_merged)) {
  stop("atac_merged does not contain lsi reduction.")
}

integration_anchors <- FindIntegrationAnchors(
  object.list = sample_objects,
  anchor.features = rownames(atac_merged),
  reduction = "rlsi",
  dims = 2:30
)
saveRDS(integration_anchors, file = file.path(out_dir, "objects", "integration_anchors_rlsi.RDS"))
cat("integration_anchors generated and saved.\n")

# 14. Integrate LSI embeddings
cat("\n--- Integrating LSI embeddings ---\n")
atac_integrated <- IntegrateEmbeddings(
  anchorset = integration_anchors,
  reductions = atac_merged[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:30
)
cat("Integrated object generated.\n")
cat("Reductions in atac_integrated:\n")
print(Reductions(atac_integrated))
saveRDS(atac_integrated, file = file.path(out_dir, "objects", "all_samples_commonpeak_integrated_before_umap.RDS"))

# 15. Integrated UMAP and clustering
cat("\n--- Running integrated UMAP / clustering ---\n")
DefaultAssay(atac_integrated) <- "ATAC"
atac_integrated <- RunUMAP(
  atac_integrated,
  reduction = "integrated_lsi",
  dims = 2:30,
  reduction.name = "umap.integrated",
  reduction.key = "UMAPI_"
)
atac_integrated <- FindNeighbors(
  atac_integrated,
  reduction = "integrated_lsi",
  dims = 2:30
)
atac_integrated <- FindClusters(
  atac_integrated,
  algorithm = 3,
  resolution = 0.5,
  verbose = FALSE
)

pdf(file.path(out_dir, "plots", "all_samples_integrated_umap_by_sample.pdf"), width = 7, height = 6)
print(
  DimPlot(
    atac_integrated,
    reduction = "umap.integrated",
    group.by = "sample_id"
  ) +
    ggtitle("Integrated ATAC UMAP by sample") +
    theme_global
)
dev.off()

pdf(file.path(out_dir, "plots", "all_samples_integrated_umap_by_cluster.pdf"), width = 7, height = 6)
print(
  DimPlot(
    atac_integrated,
    reduction = "umap.integrated",
    label = TRUE,
    repel = TRUE
  ) +
    ggtitle("Integrated ATAC UMAP by cluster") +
    theme_global
)
dev.off()

saveRDS(atac_integrated, file = file.path(out_dir, "objects", "all_samples_commonpeak_integrated.RDS"))
cat("\nAll common peak integration finished.\n")
cat("Final integrated object saved to:\n")
cat(file.path(out_dir, "objects", "all_samples_commonpeak_integrated.RDS"), "\n")

# 16. Export integrated ATAC summaries
cat("\n--- Saving integrated summary tables ---\n")
cluster_by_sample <- as.data.frame.matrix(table(atac_integrated$seurat_clusters, atac_integrated$sample_id))
write.csv(cluster_by_sample, file = file.path(out_dir, "tables", "integrated_cluster_by_sample.csv"))
sample_by_cluster <- as.data.frame.matrix(table(atac_integrated$sample_id, atac_integrated$seurat_clusters))
write.csv(sample_by_cluster, file = file.path(out_dir, "tables", "integrated_sample_by_cluster.csv"))
integrated_meta <- atac_integrated@meta.data
integrated_meta$cell <- rownames(integrated_meta)
write.csv(integrated_meta, file = file.path(out_dir, "tables", "integrated_metadata.csv"), row.names = FALSE)
cat("Integrated object cells:", ncol(atac_integrated), "\n")
cat("Integrated object peaks:", nrow(atac_integrated), "\n")
cat("Number of clusters:", length(unique(atac_integrated$seurat_clusters)), "\n")
cat("Tables saved.\n")

# 17. Configure RNA-to-ATAC label transfer
.libPaths("/home/lzh/miniconda3/envs/singlecell/lib/R/library")
set.seed(1234)
suppressPackageStartupMessages({
  library(Signac)
  library(Seurat)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(future)
})
atac_integrated_rds <- "/mnt/data_nvme1/lzh/test/AH/sc_atlas_atac/pre/batch_pre_commonpeaks/objects/all_samples_commonpeak_integrated.RDS"
rna_rds <- "/mnt/data_nvme1/lzh/test/AH/sc_atlas/AH_all_rna.fli.RDS"
rna_label_col <- "celltype"
lt_out_dir <- "/mnt/data_nvme1/lzh/test/AH/sc_atlas_atac/label_transfer"
dir.create(lt_out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(lt_out_dir, "objects"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(lt_out_dir, "plots"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(lt_out_dir, "tables"), showWarnings = FALSE, recursive = TRUE)
plan("sequential")
options(future.globals.maxSize = 100 * 1024^3)
theme_global <- theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(aspect.ratio = 1)
cat("ATAC object:", atac_integrated_rds, "\n")
cat("RNA object:", rna_rds, "\n")
cat("RNA label column:", rna_label_col, "\n")
cat("Label transfer output:", lt_out_dir, "\n")

# 18. Load RNA reference and ATAC query
cat("\n--- Loading RNA reference and ATAC query ---\n")
rna_ref <- readRDS(rna_rds)
atac_integrated <- readRDS(atac_integrated_rds)
rna_ref <- UpdateSeuratObject(rna_ref)
atac_integrated <- UpdateSeuratObject(atac_integrated)
DefaultAssay(rna_ref) <- "RNA"
DefaultAssay(atac_integrated) <- "ATAC"
cat("RNA cells:", ncol(rna_ref), "\n")
cat("RNA features:", nrow(rna_ref), "\n")
cat("ATAC cells:", ncol(atac_integrated), "\n")
cat("ATAC peaks:", nrow(atac_integrated), "\n")
cat("RNA assays:\n")
print(Assays(rna_ref))
cat("ATAC assays:\n")
print(Assays(atac_integrated))
cat("ATAC reductions:\n")
print(Reductions(atac_integrated))
if (!rna_label_col %in% colnames(rna_ref@meta.data)) {
  cat("Available RNA metadata columns:\n")
  print(colnames(rna_ref@meta.data))
  stop("rna_label_col not found in rna_ref@meta.data.")
}

cat("RNA label summary:\n")
print(table(rna_ref@meta.data[[rna_label_col]], useNA = "ifany"))

# 19. Prepare RNA reference
cat("\n--- Preparing RNA reference ---\n")
DefaultAssay(rna_ref) <- "RNA"

# Normalize RNA counts only if the normalized data slot is empty.
rna_data_dim <- dim(GetAssayData(rna_ref, assay = "RNA", slot = "data"))
rna_data_sum <- sum(GetAssayData(rna_ref, assay = "RNA", slot = "data") != 0)
cat("RNA data dim:", paste(rna_data_dim, collapse = " x "), "\n")
cat("Nonzero values in RNA data slot:", rna_data_sum, "\n")
if (rna_data_sum == 0) {
  cat("RNA data slot is empty. Running NormalizeData...\n")
  rna_ref <- NormalizeData(
    object = rna_ref,
    assay = "RNA",
    normalization.method = "LogNormalize",
    scale.factor = 10000
  )
}

# Identify variable features if insufficient features are available.
if (length(VariableFeatures(rna_ref)) < 1000) {
  cat("Finding variable features for RNA reference...\n")
  rna_ref <- FindVariableFeatures(
    object = rna_ref,
    assay = "RNA",
    nfeatures = 5000
  )
}

cat("Number of RNA variable features:", length(VariableFeatures(rna_ref)), "\n")
saveRDS(rna_ref, file = file.path(lt_out_dir, "objects", "rna_reference_prepared.RDS"))

# 20. Compute ATAC gene activity
cat("\n--- Computing ATAC gene activity ---\n")
DefaultAssay(atac_integrated) <- "ATAC"
if (!"integrated_lsi" %in% Reductions(atac_integrated)) {
  stop("integrated_lsi not found in atac_integrated. Please check ATAC integration result.")
}

# Use RNA variable genes to limit GeneActivity memory use.
rna_variable_genes <- VariableFeatures(rna_ref)
rna_variable_genes <- unique(rna_variable_genes)
cat("RNA variable genes:", length(rna_variable_genes), "\n")
gene_activities <- GeneActivity(object = atac_integrated, features = rna_variable_genes)
cat("Gene activity matrix dim:\n")
print(dim(gene_activities))
atac_integrated[["ACTIVITY"]] <- CreateAssayObject(counts = gene_activities)
DefaultAssay(atac_integrated) <- "ACTIVITY"

# Use median ATAC gene-activity depth for normalization.
activity_scale_factor <- median(atac_integrated$nCount_ACTIVITY)
cat("ACTIVITY scale factor:", activity_scale_factor, "\n")
atac_integrated <- NormalizeData(object = atac_integrated, assay = "ACTIVITY",
  normalization.method = "LogNormalize", scale.factor = activity_scale_factor)
activity_features <- intersect(rownames(atac_integrated[["ACTIVITY"]]), rownames(rna_ref[["RNA"]]))
cat("Shared genes between RNA and ATAC ACTIVITY:", length(activity_features), "\n")
if (length(activity_features) < 500) {
  stop("Too few shared genes between RNA and ATAC gene activity. Check gene naming, e.g. symbols vs Ensembl IDs.")
}

atac_integrated <- ScaleData(object = atac_integrated, assay = "ACTIVITY", features = activity_features)
saveRDS(atac_integrated, file = file.path(lt_out_dir, "objects", "atac_integrated_with_activity.RDS"))

# 21. Identify RNA-to-ATAC transfer anchors
cat("\n--- Finding transfer anchors RNA -> ATAC ---\n")
DefaultAssay(rna_ref) <- "RNA"
DefaultAssay(atac_integrated) <- "ACTIVITY"
transfer_features <- intersect( VariableFeatures(rna_ref), rownames(atac_integrated[["ACTIVITY"]]))
cat("Transfer features:", length(transfer_features), "\n")
if (length(transfer_features) < 500) {
  stop("Too few transfer features. Check gene names between RNA and ACTIVITY.")
}

transfer_anchors <- FindTransferAnchors(
  reference = rna_ref,
  query = atac_integrated,
  reference.assay = "RNA",
  query.assay = "ACTIVITY",
  features = transfer_features,
  reduction = "cca",
  dims = 1:30
)
saveRDS(transfer_anchors, file = file.path(lt_out_dir, "objects", "rna_to_atac_transfer_anchors.RDS"))
cat("Transfer anchors saved.\n")

# 22. Transfer RNA cell-type labels to ATAC
cat("\n--- Transferring RNA labels to ATAC ---\n")
predicted_labels <- TransferData(
  anchorset = transfer_anchors,
  refdata = rna_ref@meta.data[[rna_label_col]],
  weight.reduction = atac_integrated[["integrated_lsi"]],
  dims = 2:30
)
cat("Prediction result columns:\n")
print(colnames(predicted_labels))
atac_integrated <- AddMetaData(object = atac_integrated, metadata = predicted_labels)

# Keep a distinct column for the transferred RNA label.
atac_integrated$predicted_rna_celltype <- atac_integrated$predicted.id
atac_integrated$predicted_rna_score <- atac_integrated$prediction.score.max

# Apply the original prediction-confidence threshold.
atac_integrated$predicted_rna_celltype_filtered <- ifelse(atac_integrated$predicted_rna_score < 0.5, "Unknown", atac_integrated$predicted_rna_celltype)
DefaultAssay(atac_integrated) <- "ATAC"
saveRDS(atac_integrated, file = file.path(lt_out_dir, "objects", "atac_integrated_rna_labeltransfer.RDS"))
cat("Label transfer finished and object saved.\n")

# 23. Export label-transfer tables
cat("\n--- Saving label transfer tables ---\n")
lt_meta <- atac_integrated@meta.data
lt_meta$cell <- rownames(lt_meta)
write.csv(lt_meta, file = file.path(lt_out_dir, "tables", "atac_metadata_with_rna_labeltransfer.csv"), row.names = FALSE)
pred_table <- lt_meta[, c("cell", "sample_id", "seurat_clusters", "predicted_rna_celltype", "predicted_rna_score", "predicted_rna_celltype_filtered")]
write.csv(pred_table, file = file.path(lt_out_dir, "tables", "atac_rna_labeltransfer_predictions.csv"), row.names = FALSE)
celltype_by_sample <- as.data.frame.matrix(
  table(
    atac_integrated$sample_id,
    atac_integrated$predicted_rna_celltype_filtered
  )
)
write.csv(celltype_by_sample, file = file.path(lt_out_dir, "tables", "predicted_celltype_by_sample.csv"))
celltype_by_cluster <- as.data.frame.matrix(
  table(
    atac_integrated$seurat_clusters,
    atac_integrated$predicted_rna_celltype_filtered
  )
)
write.csv(celltype_by_cluster, file = file.path(lt_out_dir, "tables", "predicted_celltype_by_cluster.csv"))
cat("Tables saved to:", file.path(lt_out_dir, "tables"), "\n")

# 24. Export the subset without Unknown labels
cat("\n--- Creating no-Unknown label transfer object ---\n")
cat("Cells before removing Unknown:", ncol(atac_integrated), "\n")
atac_integrated_no_unknown <- subset(atac_integrated, subset = predicted_rna_celltype_filtered != "Unknown")
cat("Cells after removing Unknown:", ncol(atac_integrated_no_unknown), "\n")
saveRDS(atac_integrated_no_unknown, file = file.path(lt_out_dir, "objects", "atac_integrated_rna_labeltransfer_no_unknown.RDS"))
lt_meta_no_unknown <- atac_integrated_no_unknown@meta.data
lt_meta_no_unknown$cell <- rownames(lt_meta_no_unknown)
write.csv(lt_meta_no_unknown, file = file.path(lt_out_dir, "tables", "atac_metadata_with_rna_labeltransfer_no_unknown.csv"), row.names = FALSE)
pred_table_no_unknown <- lt_meta_no_unknown[, c(
  "cell",
  "sample_id",
  "seurat_clusters",
  "predicted_rna_celltype",
  "predicted_rna_score",
  "predicted_rna_celltype_filtered"
)]
write.csv(pred_table_no_unknown, file = file.path(lt_out_dir, "tables", "atac_rna_labeltransfer_predictions_no_unknown.csv"), row.names = FALSE)
celltype_by_sample_no_unknown <- as.data.frame.matrix(
  table(
    atac_integrated_no_unknown$sample_id,
    atac_integrated_no_unknown$predicted_rna_celltype_filtered
  )
)
write.csv(celltype_by_sample_no_unknown, file = file.path(lt_out_dir, "tables", "predicted_celltype_by_sample_no_unknown.csv"))
celltype_by_cluster_no_unknown <- as.data.frame.matrix(
  table(
    atac_integrated_no_unknown$seurat_clusters,
    atac_integrated_no_unknown$predicted_rna_celltype_filtered
  )
)
write.csv(celltype_by_cluster_no_unknown, file = file.path(lt_out_dir, "tables", "predicted_celltype_by_cluster_no_unknown.csv"))
cat("No-Unknown tables saved to:", file.path(lt_out_dir, "tables"), "\n")

# 25. Visualize label-transfer results
cat("\n--- Plotting label transfer results ---\n")
DefaultAssay(atac_integrated) <- "ATAC"
pal <- c("#7FB2D5","#B3CDE3","#BEDADA","#8DBA8F","#6A3D9A","#00A087","#F39B7FFF",
        "#984EA3","#BC9DCC","#57C3F3","#F5B47F","#D37295","#F4CAE4","#4581B2","#F2F2F2")

pdf(file.path(lt_out_dir, "plots", "atac_umap_predicted_rna_celltype.pdf"), width = 9, height = 7)
print(
  DimPlot(
    atac_integrated, reduction = "umap.integrated",
    group.by = "predicted_rna_celltype", cols = pal,
    label = TRUE, repel = TRUE) +
    ggtitle("scATAC labeled by scRNA reference") + theme_global
)
dev.off()

pdf(file.path(lt_out_dir, "plots", "atac_umap_predicted_rna_celltype_filtered.pdf"), width = 9, height = 7)
print(
  DimPlot(atac_integrated, reduction = "umap.integrated",
    group.by = "predicted_rna_celltype_filtered", cols = pal,
    label = TRUE, repel = TRUE) +
    ggtitle("scATAC labeled by scRNA reference, filtered") + theme_global
)
dev.off()

pdf(file.path(lt_out_dir, "plots", "atac_umap_prediction_score.pdf"), width = 7, height = 6)
print(
  FeaturePlot(atac_integrated,
    reduction = "umap.integrated",
    features = "predicted_rna_score") +
    ggtitle("RNA label transfer prediction score")
)
dev.off()

pdf(file.path(lt_out_dir, "plots", "prediction_score_by_predicted_celltype.pdf"), width = 12, height = 5)
print(
  VlnPlot(atac_integrated, features = "predicted_rna_score",
    group.by = "predicted_rna_celltype", cols = pal,
    pt.size = 0) + RotatedAxis() +
    ggtitle("Prediction score by predicted RNA cell type")
)
dev.off()

pal.sample <- c("#8da1c1","#F6B9C3","#DAECBD", "#FFD569","#D8BFDF", "#F68A8B","#de9db5","#85B293","#ab94bc")
atac_integrated$sample_id <- factor(atac_integrated$sample_id, levels = c("AH97","AH99","AH100","AH101","AH102","AH103","TH99","TH100","TH102"))
pdf(file.path(lt_out_dir, "plots", "atac_umap_by_sample_after_labeltransfer.pdf"), width = 8, height = 6)
print(
  DimPlot(atac_integrated,
    reduction = "umap.integrated",
    group.by = "sample_id", cols = pal.sample
  ) +
    ggtitle("Integrated scATAC UMAP by sample") + theme_global
)
dev.off()

pdf(file.path(lt_out_dir, "plots", "atac_umap_by_cluster_after_labeltransfer.pdf"), width = 8, height = 6)
print(
  DimPlot(atac_integrated, reduction = "umap.integrated",
    group.by = "seurat_clusters",  label = TRUE, repel = TRUE) +
    ggtitle("Integrated scATAC UMAP by cluster") + theme_global
)
dev.off()

cat("Plots saved to:", file.path(lt_out_dir, "plots"), "\n")

# 26. Visualize the subset without Unknown labels
cat("\n--- Plotting label transfer results without Unknown cells ---\n")
DefaultAssay(atac_integrated_no_unknown) <- "ATAC"
pdf(file.path(lt_out_dir, "plots", "atac_umap_predicted_rna_celltype_no_unknown.pdf"), width = 9, height = 7)
print(
  DimPlot(atac_integrated_no_unknown, reduction = "umap.integrated",
    group.by = "predicted_rna_celltype_filtered", cols = pal,
    label = TRUE, repel = TRUE) +
    ggtitle("scATAC labeled by scRNA reference, no Unknown") +
    theme_global
)
dev.off()

pdf(file.path(lt_out_dir, "plots", "atac_umap_prediction_score_no_unknown.pdf"), width = 7, height = 6)
print(
  FeaturePlot(
    atac_integrated_no_unknown,
    reduction = "umap.integrated",
    features = "predicted_rna_score") +
    ggtitle("RNA label transfer prediction score, no Unknown")
)
dev.off()

pdf(file.path(lt_out_dir, "plots", "prediction_score_by_predicted_celltype_no_unknown.pdf"), width = 12, height = 5)
print(
  VlnPlot(
    atac_integrated_no_unknown, features = "predicted_rna_score",
    group.by = "predicted_rna_celltype_filtered", cols = pal,
    pt.size = 0) + RotatedAxis() +
    ggtitle("Prediction score by predicted RNA cell type, no Unknown")
)
dev.off()

pdf(file.path(lt_out_dir, "plots", "atac_umap_by_sample_after_labeltransfer_no_unknown.pdf"), width = 8, height = 6)
print(
  DimPlot(
    atac_integrated_no_unknown, reduction = "umap.integrated",
    group.by = "sample_id", cols = pal.sample
  ) +
    ggtitle("Integrated scATAC UMAP by sample, no Unknown") +
    theme_global
)
dev.off()

pdf(file.path(lt_out_dir, "plots", "atac_umap_by_cluster_after_labeltransfer_no_unknown.pdf"), width = 8, height = 6)
print(
  DimPlot(
    atac_integrated_no_unknown,
    reduction = "umap.integrated",
    group.by = "seurat_clusters",
    label = TRUE, repel = TRUE) +
    ggtitle("Integrated scATAC UMAP by cluster, no Unknown") +
    theme_global
)
dev.off()

cat("No-Unknown plots saved to:", file.path(lt_out_dir, "plots"), "\n")

# 27. Report final label-transfer summaries
cat("\n--- Final label transfer summary ---\n")
cat("Final ATAC object:\n")
print(atac_integrated)
cat("\nPredicted cell type counts:\n")
print(table(atac_integrated$predicted_rna_celltype_filtered, useNA = "ifany"))
cat("\nPrediction score summary:\n")
print(summary(atac_integrated$predicted_rna_score))
cat("\nMain output object:\n")
cat(file.path(lt_out_dir, "objects", "atac_integrated_rna_labeltransfer.RDS"), "\n")
cat("\nMain output tables:\n")
cat(file.path(lt_out_dir, "tables", "atac_rna_labeltransfer_predictions.csv"), "\n")
cat(file.path(lt_out_dir, "tables", "predicted_celltype_by_cluster.csv"), "\n")
cat(file.path(lt_out_dir, "tables", "predicted_celltype_by_sample.csv"), "\n")
cat("\nMain output plots:\n")
cat(file.path(lt_out_dir, "plots", "atac_umap_predicted_rna_celltype.pdf"), "\n")
cat(file.path(lt_out_dir, "plots", "atac_umap_predicted_rna_celltype_filtered.pdf"), "\n")
cat(file.path(lt_out_dir, "plots", "atac_umap_prediction_score.pdf"), "\n")
cat(file.path(lt_out_dir, "plots", "prediction_score_by_predicted_celltype.pdf"), "\n")
cat(file.path(lt_out_dir, "plots", "atac_umap_by_sample_after_labeltransfer.pdf"), "\n")
cat(file.path(lt_out_dir, "plots", "atac_umap_by_cluster_after_labeltransfer.pdf"), "\n")

# 28. Report results after removing Unknown labels
if (exists("atac_integrated_no_unknown")) {
  cat("\n--- Final no-Unknown label transfer summary ---\n")
  cat("No-Unknown ATAC object:\n")
  print(atac_integrated_no_unknown)
  cat("\nPredicted cell type counts without Unknown:\n")
  print(table(atac_integrated_no_unknown$predicted_rna_celltype_filtered, useNA = "ifany"))
  cat("\nPrediction score summary without Unknown:\n")
  print(summary(atac_integrated_no_unknown$predicted_rna_score))
  cat("\nNo-Unknown output object:\n")
  cat(file.path(lt_out_dir, "objects", "atac_integrated_rna_labeltransfer_no_unknown.RDS"), "\n")
  cat("\nNo-Unknown output tables:\n")
  cat(file.path(lt_out_dir, "tables", "atac_rna_labeltransfer_predictions_no_unknown.csv"),  "\n")
  cat(file.path(lt_out_dir, "tables", "predicted_celltype_by_cluster_no_unknown.csv"), "\n")
  cat(file.path(lt_out_dir, "tables", "predicted_celltype_by_sample_no_unknown.csv"), "\n")
  cat("\nNo-Unknown output plots:\n")
  cat(file.path(lt_out_dir, "plots", "atac_umap_predicted_rna_celltype_no_unknown.pdf"), "\n")
  cat(file.path(lt_out_dir, "plots", "atac_umap_prediction_score_no_unknown.pdf"), "\n")
}
