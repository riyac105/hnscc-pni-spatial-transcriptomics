#!/usr/bin/env Rscript

# Visium SpatialCellChat analysis
# -------------------------------
# Cleaned/reproducible version of the exploratory spatialcellchat.R used for
# the manuscript analyses.
#
# Confirmed sample/object mapping from the analyst:
#   original `chat`  = GSM8633893, image slice1.2
#   original `chat1` = GSM5494476, image slice1.5
#   original `chat3` = GSM5494475, image slice1.4
#
# The original chat objects used numeric group selections in some exploratory
# plotting code; the analyst confirmed that 3 represented epithelial/tumor and
# 5 represented nerve cells. This cleaned script uses the stable group names
# "Epithelial" and "Nerve" for downstream filtering/visualization so the same
# populations are selected regardless of internal factor ordering.
#
# Spatial parameters are preserved exactly from the analysis workflow that
# generated the manuscript results:
#   spot.size = 65
#   spot_diameter_fullres = 225.54629266892906
#   tol = 100
#   interaction.range = 220
#   contact.range = 100
#   scale.distance = 0.2
#   nboot = 100
#
# This script intentionally documents the historical analysis rather than
# changing those parameters after the fact.

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(SpatialCellChat)
  library(Matrix)
  library(patchwork)
})

set.seed(42)

INDIR <- "results/07_Visium_validation"
OUTDIR <- file.path(INDIR, "SpatialCellChat")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

PNI_RDS <- file.path(INDIR, "PNI_new_visium_annotated.rds")
if (!file.exists(PNI_RDS)) {
  stop(
    "Missing ", PNI_RDS,
    ". Run 11_visium_preprocessing_ucell.R first, or point PNI_RDS to the " ,
    "analysis-ready PNI_new object."
  )
}

PNI_new <- readRDS(PNI_RDS)

# -----------------------------------------------------------------------------
# 1. Explicitly reconstruct the three sample subsets used by SpatialCellChat
# -----------------------------------------------------------------------------
sample_chat <- subset(PNI_new, subset = dataset == "GSM8633893")
sample_chat1 <- subset(PNI_new, subset = dataset == "GSM5494476")
sample_chat3 <- subset(PNI_new, subset = dataset == "GSM5494475")

sample_map <- data.frame(
  historical_object = c("chat", "chat1", "chat3"),
  dataset = c("GSM8633893", "GSM5494476", "GSM5494475"),
  image = c("slice1.2", "slice1.5", "slice1.4")
)
write.csv(
  sample_map,
  file.path(OUTDIR, "SpatialCellChat_sample_map.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 2. Spatial geometry and model parameters preserved from original workflow
# -----------------------------------------------------------------------------
spot.size <- 65
spot_diameter_fullres <- 225.54629266892906
conversion.factor <- spot.size / spot_diameter_fullres
spatial.factors <- list(ratio = conversion.factor, tol = 100)

INTERACTION_RANGE <- 220
CONTACT_RANGE <- 100
SCALE_DISTANCE <- 0.2
NBOOT <- 100

# Coordinate transform preserved exactly from the original script.
prepare_coordinates <- function(seurat_obj, image_name) {
  spatial.locs <- GetTissueCoordinates(
    seurat_obj,
    image = image_name,
    scale = NULL,
    cols = c("row", "col")
  )

  spatial.locs[, 2] <- max(spatial.locs[, 2]) - spatial.locs[, 2]
  temp_coordinates <- spatial.locs
  spatial.locs[, 1] <- temp_coordinates[, 2]
  spatial.locs[, 2] <- temp_coordinates[, 1]
  spatial.locs[, 1] <- max(spatial.locs[, 1]) - spatial.locs[, 1]
  colnames(spatial.locs) <- c("x", "y")
  spatial.locs[, 1:2, drop = FALSE]
}

# -----------------------------------------------------------------------------
# 3. Run SpatialCellChat for one tissue section
# -----------------------------------------------------------------------------
run_spatialcellchat <- function(seurat_obj, image_name, sample_id) {
  spatial.locs <- prepare_coordinates(seurat_obj, image_name)

  chat_obj <- createSpatialCellChat(
    object = seurat_obj,
    group.by = "nerve_annotation",
    assay = "SCT",
    datatype = "spatial",
    coordinates = spatial.locs,
    spatial.factors = spatial.factors
  )

  # Human signaling database; retain Secreted, ECM-Receptor, and Cell-Cell
  # Contact classes, as in the original workflow.
  db_use <- subsetDB(
    CellChatDB.human,
    search = c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact"),
    non_protein = FALSE
  )
  chat_obj@DB <- db_use

  chat_obj <- subsetData(chat_obj)
  chat_obj <- preProcessing(chat_obj)
  chat_obj <- identifyOverExpressedGenes(
    chat_obj,
    selection.method = "meringue",
    do.grid = FALSE
  )
  chat_obj <- identifyOverExpressedInteractions(
    chat_obj,
    variable.both = FALSE
  )

  chat_obj <- computeCommunProb(
    chat_obj,
    distance.use = TRUE,
    scale.distance = SCALE_DISTANCE,
    contact.dependent = TRUE,
    interaction.range = INTERACTION_RANGE,
    contact.range = CONTACT_RANGE
  )

  chat_obj <- filterProbability(chat_obj)
  chat_obj <- filterCommunication(
    chat_obj,
    min.cells = NULL,
    min.links = 10,
    min.cells.sr = 10
  )

  chat_obj <- computeAvgCommunProb(
    chat_obj,
    avg.type = "avg",
    nboot = NBOOT,
    do.permutation = TRUE
  )

  chat_obj <- computeCommunProbPathway(chat_obj)
  chat_obj <- aggregateNet(chat_obj)

  # Network centrality calculations retained from the original workflow.
  chat_obj <- netAnalysis_computeCentrality(
    chat_obj,
    slot.name = "net",
    do.group = FALSE,
    degree.only = TRUE
  )
  chat_obj <- netAnalysis_computeCentrality(
    chat_obj,
    slot.name = "net",
    do.group = TRUE,
    degree.only = TRUE
  )
  chat_obj <- netAnalysis_computeCentrality(
    chat_obj,
    slot.name = "netP",
    do.group = FALSE,
    degree.only = TRUE
  )
  chat_obj <- netAnalysis_computeCentrality(
    chat_obj,
    slot.name = "netP",
    do.group = TRUE,
    degree.only = TRUE
  )

  saveRDS(
    chat_obj,
    file.path(OUTDIR, paste0("SpatialCellChat_", sample_id, ".rds"))
  )

  chat_obj
}

chat <- run_spatialcellchat(
  sample_chat,
  image_name = "slice1.2",
  sample_id = "GSM8633893"
)
chat1 <- run_spatialcellchat(
  sample_chat1,
  image_name = "slice1.5",
  sample_id = "GSM5494476"
)
chat3 <- run_spatialcellchat(
  sample_chat3,
  image_name = "slice1.4",
  sample_id = "GSM5494475"
)

# -----------------------------------------------------------------------------
# 4. Epithelial <-> nerve interactions used for manuscript visualization
# -----------------------------------------------------------------------------
# Original selected pathway sets for the three chat objects.
pathways_chat <- c("WNT", "PTN", "EPHA", "NECTIN", "CEACAM")
pathways_chat1 <- c("MK", "CEACAM", "LAMININ", "SEMA3", "WNT")
pathways_chat3 <- c("WNT", "SEMA7", "SEMA4", "SEMA3", "EPHA", "EPHB")

make_epi_nerve_bubble <- function(chat_obj, pathways, sample_id) {
  pairLR <- extractEnrichedLR(
    chat_obj,
    signaling = pathways,
    geneLR.return = FALSE
  )

  # Stable named populations replace the historical numeric selections
  # sources.use/targets.use = c(3, 5), where 3=epithelial/tumor and 5=nerve.
  p <- netVisual_bubble(
    chat_obj,
    sources.use = c("Epithelial", "Nerve"),
    targets.use = c("Epithelial", "Nerve"),
    remove.isolate = TRUE,
    pairLR.use = pairLR
  )

  pdf(
    file.path(OUTDIR, paste0("bubble_epithelial_nerve_", sample_id, ".pdf")),
    width = 8,
    height = 7
  )
  print(p)
  dev.off()

  # Save the corresponding interaction table for transparent reporting.
  epi_to_nerve <- subsetCommunication(
    chat_obj,
    sources.use = "Epithelial",
    targets.use = "Nerve",
    thresh = 0.05
  )
  nerve_to_epi <- subsetCommunication(
    chat_obj,
    sources.use = "Nerve",
    targets.use = "Epithelial",
    thresh = 0.05
  )

  write.csv(
    epi_to_nerve,
    file.path(OUTDIR, paste0("interactions_Epithelial_to_Nerve_", sample_id, ".csv")),
    row.names = FALSE
  )
  write.csv(
    nerve_to_epi,
    file.path(OUTDIR, paste0("interactions_Nerve_to_Epithelial_", sample_id, ".csv")),
    row.names = FALSE
  )

  invisible(p)
}

p_chat <- make_epi_nerve_bubble(chat, pathways_chat, "GSM8633893")
p_chat1 <- make_epi_nerve_bubble(chat1, pathways_chat1, "GSM5494476")
p_chat3 <- make_epi_nerve_bubble(chat3, pathways_chat3, "GSM5494475")

# Combined three-sample bubble figure, corresponding to the type of visualization
# described for manuscript Figure 5E.
p_combined <- wrap_plots(
  p_chat + ggtitle("GSM8633893"),
  p_chat1 + ggtitle("GSM5494476"),
  p_chat3 + ggtitle("GSM5494475"),
  ncol = 1
)

ggsave(
  file.path(OUTDIR, "SpatialCellChat_epithelial_nerve_three_samples.pdf"),
  p_combined,
  width = 9,
  height = 18
)

# Save parameter table so the exact distance/permutation settings are explicit.
parameter_table <- data.frame(
  parameter = c(
    "spot.size",
    "spot_diameter_fullres",
    "conversion.factor",
    "tol",
    "interaction.range",
    "contact.range",
    "scale.distance",
    "nboot",
    "min.links",
    "min.cells.sr"
  ),
  value = c(
    spot.size,
    spot_diameter_fullres,
    conversion.factor,
    spatial.factors$tol,
    INTERACTION_RANGE,
    CONTACT_RANGE,
    SCALE_DISTANCE,
    NBOOT,
    10,
    10
  )
)
write.csv(
  parameter_table,
  file.path(OUTDIR, "SpatialCellChat_parameters.csv"),
  row.names = FALSE
)

cat("\nSpatialCellChat analysis complete.\n")
cat("Sample mapping:\n")
print(sample_map)
cat("\nParameters preserved from original workflow:\n")
print(parameter_table)
cat("\nSaved outputs to:", OUTDIR, "\n")
