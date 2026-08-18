#!/usr/bin/env Rscript

# Visium preprocessing, annotation, and UCell analyses
# ---------------------------------------------------
# Cleaned/reproducible version of the exploratory MergeDatasetsWorkflow.R used
# for the manuscript analyses.
#
# Analysis represented here:
#   - 7 Visium HNSCC tissue sections
#       GSE281978 primary tumors: GSM8633891, GSM8633893, GSM8633895
#       GSE181300 invasive-front sections: GSM5494475, GSM5494476,
#                                          GSM5494479, GSM5494480
#   - QC thresholds preserved from the original workflow
#   - SCTransform normalization, PCA, Harmony integration, clustering, UMAP
#   - cell-type labels preserved from the original cluster annotation
#   - nerve-associated spots: UCell score for S100B/SOX10/MPZ, cutoff >= 0.1
#   - literature-derived PNI and EMT UCell scores
#   - median split of PNI_UCell and Kruskal-Wallis comparison of EMT_UCell
#
# IMPORTANT PROVENANCE NOTE:
# The original exploratory workflow referenced a pre-existing `nerves_UCell`
# column but did not retain the line that created it. The UCell nerve-scoring
# step below reconstructs the method documented in the manuscript
# (S100B/SOX10/MPZ; cutoff 0.1). The PNI/EMT UCell code and gene sets are taken
# directly from the analysis snippet supplied by the analyst.
#
# Configure data locations with environment variables if desired:
#   export GSE281978_ROOT=/path/to/GSE281978
#   export GSE181300_ROOT=/path/to/GSE181300
# Default locations are project-relative under data/Visium/.

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(tidyverse)
  library(glmGamPoi)
  library(UCell)
  library(harmony)
})

set.seed(42)

OUTDIR <- "results/07_Visium_validation"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

GSE281978_ROOT <- Sys.getenv(
  "GSE281978_ROOT",
  unset = "data/Visium/GSE281978"
)
GSE181300_ROOT <- Sys.getenv(
  "GSE181300_ROOT",
  unset = "data/Visium/GSE181300"
)

# -----------------------------------------------------------------------------
# 1. Load the seven manuscript Visium sections
# -----------------------------------------------------------------------------
load_visium_sample <- function(root, gsm) {
  data_dir <- file.path(root, gsm)
  h5_file <- paste0(gsm, "_filtered_feature_bc_matrix.h5")
  obj <- Load10X_Spatial(data.dir = data_dir, filename = h5_file)
  obj$dataset <- gsm
  obj
}

object1 <- load_visium_sample(GSE281978_ROOT, "GSM8633891")
object2 <- load_visium_sample(GSE281978_ROOT, "GSM8633893")
object3 <- load_visium_sample(GSE281978_ROOT, "GSM8633895")
object4 <- load_visium_sample(GSE181300_ROOT, "GSM5494475")
object5 <- load_visium_sample(GSE181300_ROOT, "GSM5494476")
object6 <- load_visium_sample(GSE181300_ROOT, "GSM5494479")
object7 <- load_visium_sample(GSE181300_ROOT, "GSM5494480")

PNI <- merge(
  x = object1,
  y = list(object2, object3, object4, object5, object6, object7),
  add.cell.ids = c(
    "GSM8633891", "GSM8633893", "GSM8633895",
    "GSM5494475", "GSM5494476", "GSM5494479", "GSM5494480"
  ),
  project = "MergedProject"
)

rm(object1, object2, object3, object4, object5, object6, object7)

# Sample-to-image mapping produced by the merged Seurat object.
sample_image_map <- data.frame(
  dataset = c(
    "GSM8633891", "GSM8633893", "GSM8633895",
    "GSM5494475", "GSM5494476", "GSM5494479", "GSM5494480"
  ),
  image = c(
    "slice1", "slice1.2", "slice1.3",
    "slice1.4", "slice1.5", "slice1.6", "slice1.7"
  )
)
write.csv(
  sample_image_map,
  file.path(OUTDIR, "Visium_sample_image_map.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 2. Quality control and normalization
# -----------------------------------------------------------------------------
PNI[["percent.mt"]] <- PercentageFeatureSet(PNI, pattern = "^MT-")
PNI[["percent.ribo"]] <- PercentageFeatureSet(PNI, pattern = "^RP[SL]")

PNI <- subset(
  PNI,
  subset = nFeature_Spatial < 7500 &
           nFeature_Spatial > 200 &
           nCount_Spatial < 50000 &
           nCount_Spatial > 250 &
           percent.mt < 15 &
           percent.ribo < 40
)

# In the original workflow, the four GSE181300 sections had hires images only;
# their lowres scaling factors were therefore set equal to hires.
for (img in c("slice1.4", "slice1.5", "slice1.6", "slice1.7")) {
  PNI@images[[img]]@scale.factors$lowres <- PNI@images[[img]]@scale.factors$hires
}

PNI <- SCTransform(PNI, assay = "Spatial")

# -----------------------------------------------------------------------------
# 3. PCA, clustering, Harmony integration, and UMAP
# -----------------------------------------------------------------------------
PNI <- RunPCA(PNI, assay = "SCT", reduction.name = "pca.SCT")
PNI <- FindNeighbors(PNI, assay = "SCT", reduction = "pca.SCT", dims = 1:15)
PNI <- FindClusters(PNI, cluster.name = "seurat_cluster.SCT", resolution = 0.5)
PNI <- RunUMAP(
  PNI,
  reduction = "pca.SCT",
  reduction.name = "umap.SCT",
  return.model = TRUE,
  dims = 1:15
)

PNI <- RunHarmony(
  object = PNI,
  group.by.vars = "dataset",
  assay.use = "SCT",
  reduction = "pca.SCT",
  reduction.save = "harmony.SCT",
  theta = 8
)

PNI <- FindNeighbors(PNI, reduction = "harmony.SCT", dims = 1:15)
PNI <- FindClusters(
  PNI,
  resolution = 0.5,
  cluster.name = "seurat_cluster.harmony.SCT"
)
PNI <- RunUMAP(
  PNI,
  reduction = "harmony.SCT",
  reduction.name = "umap.harmony.SCT",
  dims = 1:15
)

Idents(PNI) <- "seurat_cluster.harmony.SCT"

# -----------------------------------------------------------------------------
# 4. Cell-type annotation used in the Visium analysis
# -----------------------------------------------------------------------------
PNI_new <- RenameIdents(
  PNI,
  "0"  = "B Cell",
  "1"  = "Epithelial",
  "2"  = "Epithelial",
  "3"  = "Unknown",
  "4"  = "T Cell",
  "5"  = "Fibroblast",
  "6"  = "Endothelial",
  "7"  = "Unknown",
  "8"  = "Epithelial",
  "9"  = "T Cell",
  "10" = "Fibroblast",
  "11" = "B Cell",
  "12" = "Epithelial",
  "13" = "Epithelial",
  "14" = "Unknown",
  "15" = "Unknown"
)

# -----------------------------------------------------------------------------
# 5. Nerve-associated spot annotation with UCell
# -----------------------------------------------------------------------------
# This reconstructs the documented nerve-associated score that generated the
# `nerves_UCell` metadata used later in the original exploratory workflow.
nerve_markers <- list(
  nerves = c("S100B", "SOX10", "MPZ")
)
PNI_new <- AddModuleScore_UCell(PNI_new, features = nerve_markers)

PNI_new$nerve_annotation <- ifelse(
  PNI_new$nerves_UCell >= 0.1,
  "Nerve",
  as.character(Idents(PNI_new))
)
Idents(PNI_new) <- "nerve_annotation"
PNI_new$cell_types <- as.character(Idents(PNI_new))

# -----------------------------------------------------------------------------
# 6. PNI-related and EMT UCell scores
# -----------------------------------------------------------------------------
# Gene sets supplied in the original analysis snippet.
markers <- list(
  PNI = c(
    "S100B", "SOX10",                 # nerve markers
    "NGF", "NTRK1",                 # neurotrophic factors
    "CXCR4", "CXCL12", "SEMA4D", "PLXNB1", "EFNA4", # axon guidance
    "MMP2", "MMP9",                 # ECM-degrading molecules
    "NCAM1", "L1CAM"                # adhesion factors
  ),
  EMT = c("VIM", "CDH2", "FN1", "SNAI1", "TWIST1")
)

PNI_new <- AddModuleScore_UCell(PNI_new, features = markers)

# Median dichotomization used for the spot-level PNI-vs-EMT comparison.
pni_median <- median(PNI_new$PNI_UCell, na.rm = TRUE)
PNI_new$PNI_status <- ifelse(
  PNI_new$PNI_UCell < pni_median,
  "LOW",
  "HIGH"
)
PNI_new$PNI_status <- factor(PNI_new$PNI_status, levels = c("LOW", "HIGH"))

kw <- kruskal.test(EMT_UCell ~ PNI_status, data = PNI_new@meta.data)

kw_out <- data.frame(
  statistic = unname(kw$statistic),
  df = unname(kw$parameter),
  p_value = kw$p.value,
  median_PNI_UCell_cutoff = pni_median,
  total_spots = ncol(PNI_new),
  nerve_associated_spots = sum(PNI_new$nerve_annotation == "Nerve", na.rm = TRUE)
)
write.csv(
  kw_out,
  file.path(OUTDIR, "PNI_high_low_EMT_kruskal_wallis.csv"),
  row.names = FALSE
)

write.csv(
  as.data.frame(table(PNI_new$dataset, PNI_new$nerve_annotation)),
  file.path(OUTDIR, "nerve_annotation_counts_by_dataset.csv"),
  row.names = FALSE
)

write.csv(
  PNI_new@meta.data,
  file.path(OUTDIR, "Visium_spot_metadata_with_UCell_scores.csv")
)

# -----------------------------------------------------------------------------
# 7. Manuscript-supporting visualizations (Figure 5A-D components)
# -----------------------------------------------------------------------------
p_umap <- DimPlot(
  PNI_new,
  reduction = "umap.harmony.SCT",
  group.by = "nerve_annotation",
  label = FALSE
) + ggtitle("Visium spot annotations")

ggsave(
  file.path(OUTDIR, "Visium_UMAP_cell_annotations.pdf"),
  p_umap,
  width = 7,
  height = 5
)

image_names <- sample_image_map$image

p_spatial <- SpatialDimPlot(
  PNI_new,
  images = image_names,
  group.by = "nerve_annotation",
  label = FALSE,
  combine = TRUE,
  ncol = 4
)
ggsave(
  file.path(OUTDIR, "Visium_spatial_cell_annotations.pdf"),
  p_spatial,
  width = 14,
  height = 8
)

p_pni <- SpatialFeaturePlot(
  PNI_new,
  images = image_names,
  features = "PNI_UCell",
  combine = TRUE,
  ncol = 4
)
ggsave(
  file.path(OUTDIR, "Visium_spatial_PNI_UCell.pdf"),
  p_pni,
  width = 14,
  height = 8
)

p_emt <- ggplot(
  PNI_new@meta.data,
  aes(x = PNI_status, y = EMT_UCell)
) +
  geom_boxplot(staplewidth = 0.2) +
  theme_minimal() +
  xlab("PNI gene-set status") +
  ylab("EMT UCell score") +
  ggtitle("EMT score by PNI gene-set status")

ggsave(
  file.path(OUTDIR, "Visium_PNI_high_low_EMT_boxplot.pdf"),
  p_emt,
  width = 5.5,
  height = 4.5
)

# Save the analysis-ready object for the SpatialCellChat script.
saveRDS(
  PNI_new,
  file.path(OUTDIR, "PNI_new_visium_annotated.rds")
)

cat("\nVisium preprocessing/UCell analysis complete.\n")
cat("Total spots:", ncol(PNI_new), "\n")
cat(
  "Nerve-associated spots:",
  sum(PNI_new$nerve_annotation == "Nerve", na.rm = TRUE),
  "\n"
)
cat("Kruskal-Wallis p-value:", kw$p.value, "\n")
cat("Saved outputs to:", OUTDIR, "\n")
